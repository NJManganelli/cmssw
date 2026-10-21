// Per-track track-quality heads (conifer BDTs) as an extension table on a track table: N named heads, each
// with one single-output model per window, read from JSON files named in the configuration (joint-head models
// swap in by configuration). The features are computed from the nano tables' STORED values (mantissa reduction
// included), so they equal the ones the models were trained on, bit for bit in float32
// (/Volumes/WDMac/smartpixels_scratch/heads_cmssw/SPEC.md, the tagger's heads):
//
//   collection "OT"           features 1-13 from the OT track table (otTables), one window "all"
//   collection "SmartPixels"  features 1-38: 1-13 from the OT track the refit started from (row aligned), 14-31
//                             from the SmartPixels track table (spTables), 32-38 aggregated over the track's
//                             accepted rows of the refit hit table (hitTable); windows "prompt" (|d0| <= dPromptMax,
//                             d0 the SmartPixels track's) and "displaced"
//
// A feature that is missing or not finite is -999 (the models' missing value); a refit column <= -998 too.
// Output, per head: the margin (init_predict + the trees' leaf values), or 1 / (1 + exp(-margin)) for a head
// whose output is "probability" (clean).
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/NanoAOD/interface/FlatTable.h"

#include "ap_fixed.h"
#include "conifer.h"
#include "nlohmann/json.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace {
  constexpr float kMissing = -999.f;
  constexpr unsigned kOTFeatures = 13, kSPFeatures = 38;

  // the columns of one logical table, spread over its main table and its extensions
  class TableView {
  public:
    void add(const nanoaod::FlatTable& t) {
      if (!tables_.empty() && t.nRows() != tables_.front()->nRows())
        throw cms::Exception("L1SmartPixelsTrackHeads")
            << t.name() << ": " << t.nRows() << " rows in one product, " << tables_.front()->nRows() << " in another";
      tables_.push_back(&t);
    }
    size_t rows() const { return tables_.empty() ? 0 : tables_.front()->nRows(); }
    // a column as doubles (the stored value, converted exactly)
    std::vector<double> column(const std::string& name) const {
      for (const auto* t : tables_) {
        const int c = t->columnIndex(name);
        if (c < 0)
          continue;
        std::vector<double> v(t->nRows());
        using CT = nanoaod::FlatTable::ColumnType;
        const auto fill = [&](auto data) { std::copy(data.begin(), data.end(), v.begin()); };
        switch (t->columnType(c)) {
          case CT::UInt8:
            fill(t->columnData<uint8_t>(c));
            break;
          case CT::Int16:
            fill(t->columnData<int16_t>(c));
            break;
          case CT::UInt16:
            fill(t->columnData<uint16_t>(c));
            break;
          case CT::Int32:
            fill(t->columnData<int32_t>(c));
            break;
          case CT::UInt32:
            fill(t->columnData<uint32_t>(c));
            break;
          case CT::Int64:
            fill(t->columnData<int64_t>(c));
            break;
          case CT::UInt64:
            fill(t->columnData<uint64_t>(c));
            break;
          case CT::Bool:
            fill(t->columnData<bool>(c));
            break;
          case CT::Float:
            fill(t->columnData<float>(c));
            break;
          case CT::Double:
            fill(t->columnData<double>(c));
            break;
        }
        return v;
      }
      throw cms::Exception("L1SmartPixelsTrackHeads")
          << "column " << name << " is in none of the "
          << (tables_.empty() ? std::string("(no)") : tables_.front()->name()) << " table products";
    }

  private:
    std::vector<const nanoaod::FlatTable*> tables_;
  };

  float featureValue(double x) { return std::isfinite(x) ? static_cast<float>(x) : kMissing; }
  float refitColumn(double x) { return (std::isfinite(x) && x > -998.) ? static_cast<float>(x) : kMissing; }
  double inv(double x) { return 1.0 / std::max(x, 1e-3); }

  // zero bits strictly between the lowest and the highest set bit of the lowest 7 bits (trkquality.nlaymiss_interior)
  int nmissInterior(int64_t hitPattern) {
    if (hitPattern <= 0)
      return 0;
    constexpr int width = 7;
    int first = 0, last = width - 1;
    while (first < width && !((hitPattern >> first) & 1))
      ++first;
    while (last >= 0 && !((hitPattern >> last) & 1))
      --last;
    if (first == width)  // no bit set among the lowest 7: the reference's argmax gives the full range
      return width;
    int n = 0;
    for (int b = first; b <= last; ++b)
      n += !((hitPattern >> b) & 1);
    return n;
  }

  std::string resolve(const std::string& path) {
    if (std::filesystem::exists(path))
      return path;
    return edm::FileInPath(path).fullPath();
  }
}  // namespace

class L1SmartPixelsTrackHeadsTableProducer : public edm::global::EDProducer<> {
public:
  using BDT = conifer::BDT<float, float>;
  struct Head {
    std::string name;
    bool probability = false;
    std::map<std::string, std::unique_ptr<BDT>> models;  // window -> model
  };

  explicit L1SmartPixelsTrackHeadsTableProducer(const edm::ParameterSet& cfg)
      : collection_(cfg.getParameter<std::string>("collection")),
        tableName_(cfg.getParameter<std::string>("tableName")),
        prefix_(cfg.getParameter<std::string>("columnPrefix")),
        dPromptMax_(cfg.getParameter<double>("dPromptMax")),
        writeFeatures_(cfg.getParameter<bool>("writeFeatures")) {
    if (collection_ != "OT" && collection_ != "SmartPixels")
      throw cms::Exception("Configuration") << "collection must be OT or SmartPixels, not " << collection_;
    const bool sp = collection_ == "SmartPixels";
    for (const auto& tag : cfg.getParameter<std::vector<edm::InputTag>>("otTables"))
      otTokens_.push_back(consumes<nanoaod::FlatTable>(tag));
    if (sp) {
      for (const auto& tag : cfg.getParameter<std::vector<edm::InputTag>>("spTables"))
        spTokens_.push_back(consumes<nanoaod::FlatTable>(tag));
      hitToken_ = consumes<nanoaod::FlatTable>(cfg.getParameter<edm::InputTag>("hitTable"));
    }
    if (otTokens_.empty() || (sp && spTokens_.empty()))
      throw cms::Exception("Configuration") << "otTables (and, for SmartPixels, spTables) must not be empty";
    const std::vector<std::string> windows = sp ? std::vector<std::string>{"prompt", "displaced"}
                                                : std::vector<std::string>{"all"};
    const unsigned nFeatures = sp ? kSPFeatures : kOTFeatures;
    for (const auto& h : cfg.getParameter<std::vector<edm::ParameterSet>>("heads")) {
      Head head;
      head.name = h.getParameter<std::string>("name");
      const auto output = h.getParameter<std::string>("output");
      if (output != "margin" && output != "probability")
        throw cms::Exception("Configuration") << "head " << head.name << ": output must be margin or probability";
      head.probability = output == "probability";
      for (const auto& m : h.getParameter<std::vector<edm::ParameterSet>>("models")) {
        const auto window = m.getParameter<std::string>("window");
        if (std::find(windows.begin(), windows.end(), window) == windows.end())
          throw cms::Exception("Configuration") << "head " << head.name << ": window " << window << " is not one of "
                                                << collection_ << "'s";
        const auto file = resolve(m.getParameter<std::string>("file"));
        // conifer keeps the model's shape private: read it from the JSON (n_classes 2 = one binary output)
        std::ifstream in(file);
        const auto j = nlohmann::json::parse(in);
        const auto nf = j.at("n_features").get<unsigned>(), nc = j.at("n_classes").get<unsigned>();
        if (nf != nFeatures || (nc != 1 && nc != 2))
          throw cms::Exception("Configuration")
              << "head " << head.name << " / " << window << ": " << file << " has " << nf << " features and "
              << nc << " classes, " << collection_ << " needs " << nFeatures << " features and one output";
        head.models[window] = std::make_unique<BDT>(file);
      }
      for (const auto& w : windows)
        if (!head.models.count(w))
          throw cms::Exception("Configuration") << "head " << head.name << ": no model for window " << w;
      heads_.push_back(std::move(head));
    }
    produces<nanoaod::FlatTable>();
  }

  void produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup&) const override {
    TableView ot;
    for (const auto& t : otTokens_)
      ot.add(iEvent.get(t));
    const size_t n = ot.rows();
    std::vector<std::vector<float>> X(n);

    // 1-13: the OT track
    {
      const auto pt = ot.column("pt"), tanL = ot.column("tanL"), eta = ot.column("eta"), nStubs = ot.column("nStubs"),
                 hp = ot.column("hitPattern"), cxy = ot.column("chi2XYRed"), cz = ot.column("chi2ZRed"),
                 cb = ot.column("chi2BendRed"), mva = ot.column("trkMVA1"), crr = ot.column("cov_rInv_rInv"),
                 cpp = ot.column("cov_phi_phi"), ctt = ot.column("cov_tanL_tanL"), czz = ot.column("cov_z0_z0");
      for (size_t i = 0; i < n; ++i)
        X[i] = {featureValue(inv(pt[i])),
                featureValue(tanL[i]),
                featureValue(std::abs(eta[i])),
                featureValue(nStubs[i]),
                featureValue(nmissInterior(static_cast<int64_t>(hp[i]))),
                featureValue(cxy[i]),
                featureValue(cz[i]),
                featureValue(cb[i]),
                featureValue(mva[i]),
                featureValue(std::sqrt(std::abs(crr[i]))),
                featureValue(std::sqrt(std::abs(cpp[i]))),
                featureValue(std::sqrt(std::abs(ctt[i]))),
                featureValue(std::sqrt(std::abs(czz[i])))};
    }

    std::vector<std::string> window(n, "all");
    if (collection_ == "SmartPixels") {
      TableView sp;
      for (const auto& t : spTokens_)
        sp.add(iEvent.get(t));
      if (sp.rows() != n)
        throw cms::Exception("L1SmartPixelsTrackHeads")
            << "the SmartPixels track table has " << sp.rows() << " rows, the OT track table " << n;
      const auto bpt = ot.column("pt"), btanL = ot.column("tanL"), bz0 = ot.column("z0");
      const auto rpt = sp.column("pt"), rtanL = sp.column("tanL"), rz0 = sp.column("z0"), rd0 = sp.column("d0");
      static const char* const kRefit[] = {"spixRefitPerformed", "spixAnyWindowTruncated", "spixNCrossings",
                                           "spixNAcceptedHits", "spixLayerHitMask", "spixMaxWindowMult",
                                           "spixChi2IncXTot", "spixChi2IncYTot", "spixChi2IncAlphaTot",
                                           "spixChi2ITAtOT", "spixChi2ITAtRefit", "spixShiftChi2",
                                           "spixLogDetRatio"};
      std::vector<std::vector<double>> refit;
      for (const auto* c : kRefit)
        refit.push_back(sp.column(c));
      for (size_t i = 0; i < n; ++i) {
        auto& x = X[i];
        x.push_back(featureValue(inv(rpt[i])));
        x.push_back(featureValue(rtanL[i]));
        x.push_back(featureValue(std::abs(rz0[i] - bz0[i])));
        x.push_back(featureValue(std::abs(rtanL[i] - btanL[i])));
        x.push_back(featureValue(inv(rpt[i]) - inv(bpt[i])));
        for (const auto& c : refit)
          x.push_back(refitColumn(c[i]));
        window[i] = std::abs(rd0[i]) <= dPromptMax_ ? "prompt" : "displaced";
      }

      // 32-38: the track's accepted refit hits
      TableView hits;
      hits.add(iEvent.get(hitToken_));
      const auto trk = hits.column("trackIdx"), acc = hits.column("hitAccepted"), hasA = hits.column("hasAlpha"),
                 c2x = hits.column("chi2IncX"), c2y = hits.column("chi2IncY"), c2a = hits.column("chi2IncAlpha"),
                 wm = hits.column("windowMult"), pa = hits.column("pullAlpha"), px = hits.column("pullX"),
                 py = hits.column("pullY");
      std::vector<double> nacc(n, 0.), nangle(n, 0.), sumAngle(n, 0.), sumPos(n, 0.), sumWm(n, 0.), max1(n, 0.),
          max2(n, 0.), maxPos(n, 0.);
      std::vector<int> nPull(n, 0);
      for (size_t h = 0; h < trk.size(); ++h) {
        if (!(acc[h] > 0))
          continue;
        const auto t = static_cast<long>(trk[h]);
        if (t < 0 || static_cast<size_t>(t) >= n)
          throw cms::Exception("L1SmartPixelsTrackHeads") << "hit trackIdx " << t << " outside " << n << " tracks";
        nacc[t] += 1.;
        nangle[t] += hasA[h];
        sumAngle[t] += c2a[h] * hasA[h];
        sumPos[t] += c2x[h] + c2y[h];
        sumWm[t] += wm[h];
        const double a = hasA[h] > 0 ? std::abs(pa[h]) : 0.;
        // the largest and the largest after removing ONE occurrence of the maximum
        if (nPull[t] == 0 || a > max1[t]) {
          max2[t] = nPull[t] == 0 ? 0. : max1[t];
          max1[t] = a;
        } else if (nPull[t] == 1 || a > max2[t])
          max2[t] = a;
        ++nPull[t];
        maxPos[t] = std::max(maxPos[t], std::max(std::abs(px[h]), std::abs(py[h])));
      }
      for (size_t i = 0; i < n; ++i) {
        auto& x = X[i];
        x.push_back(featureValue(nangle[i]));
        x.push_back(featureValue(sumAngle[i] / std::max(nangle[i], 1.)));
        x.push_back(featureValue(sumPos[i] / std::max(nacc[i], 1.)));
        x.push_back(featureValue(sumWm[i] / std::max(nacc[i], 1.)));
        x.push_back(featureValue(max1[i]));
        x.push_back(featureValue(nPull[i] < 2 ? 0. : max2[i]));
        x.push_back(featureValue(maxPos[i]));
      }
    }

    auto table = std::make_unique<nanoaod::FlatTable>(n, tableName_, false, true);
    for (const auto& head : heads_) {
      std::vector<float> out(n);
      for (size_t i = 0; i < n; ++i) {
        const float m = head.models.at(window[i])->decision_function(X[i]).at(0);
        out[i] = head.probability ? 1.f / (1.f + std::exp(-m)) : m;
      }
      table->addColumn<float>(prefix_ + head.name,
                              out,
                              collection_ + " track-quality head " + head.name +
                                  (head.probability ? ": P = 1/(1+exp(-margin))" : ": margin (log sigma_hat)"));
    }
    if (writeFeatures_ && n > 0)
      for (size_t f = 0; f < X.front().size(); ++f) {
        std::vector<float> v(n);
        for (size_t i = 0; i < n; ++i)
          v[i] = X[i][f];
        table->addColumn<float>(prefix_ + "feature" + std::to_string(f), v, "TEST ONLY: model input " +
                                                                                std::to_string(f) + " (SPEC.md order)");
      }
    iEvent.put(std::move(table));
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription model;
    model.add<std::string>("window")->setComment("all (OT) | prompt | displaced (SmartPixels)");
    model.add<std::string>("file")->setComment("conifer JSON; a path, or relative to the CMSSW search path");
    edm::ParameterSetDescription head;
    head.add<std::string>("name");
    head.add<std::string>("output", "margin")->setComment("margin | probability (1/(1+exp(-margin)))");
    head.addVPSet("models", model, {});
    edm::ParameterSetDescription desc;
    desc.add<std::string>("collection", "SmartPixels")->setComment("OT (features 1-13) | SmartPixels (1-38)");
    desc.add<std::string>("tableName", "L1TSmartPixelsExtTrack")->setComment("the extension table's name");
    desc.add<std::string>("columnPrefix", "head_");
    desc.add<std::vector<edm::InputTag>>("otTables", {})
        ->setComment("the FlatTable products of the OT track table (main and extensions)");
    desc.add<std::vector<edm::InputTag>>("spTables", {})
        ->setComment("the FlatTable products of the SmartPixels track table (SmartPixels only)");
    desc.add<edm::InputTag>("hitTable", edm::InputTag(""))->setComment("the refit hit table (SmartPixels only)");
    desc.add<double>("dPromptMax", 0.02)->setComment("|d0| [cm] up to which a SmartPixels track is prompt");
    desc.addVPSet("heads", head, {});
    desc.add<bool>("writeFeatures", false)->setComment("TEST ONLY: also write the model inputs (parity checks)");
    descriptions.addWithDefaultLabel(desc);
  }

private:
  const std::string collection_, tableName_, prefix_;
  const double dPromptMax_;
  const bool writeFeatures_;
  std::vector<edm::EDGetTokenT<nanoaod::FlatTable>> otTokens_, spTokens_;
  edm::EDGetTokenT<nanoaod::FlatTable> hitToken_;
  std::vector<Head> heads_;
};

DEFINE_FWK_MODULE(L1SmartPixelsTrackHeadsTableProducer);
