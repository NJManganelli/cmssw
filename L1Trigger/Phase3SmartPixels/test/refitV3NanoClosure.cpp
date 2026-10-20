// Development closure of the refit v3 engine against the replay on IDENTICAL inputs: the engine is fed the
// nano's OT tracks and SmartPixels clusters (the values the replay reads), the replay's geometry JSON and
// its chunk-median kpt, and compared track by track with the replay's reference output (spec/reference).
//   refitV3NanoClosure <nano.root> <reference.jsonl> <geometry.json> <manifest.json> <mask AAAA|AAAI|...>
//                      [nEvents=10] [knownBad=""]
// Exit code 0 when every compared quantity agrees within the development tolerances.
#include "L1Trigger/Phase3SmartPixels/interface/SmartPixelsRefitV3.h"

#include <nlohmann/json.hpp>

#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <tuple>

using namespace smartpixels::v3;

namespace {
  Geometry loadGeometry(const std::string& path) {
    std::ifstream in(path);
    const auto j = nlohmann::json::parse(in);
    Geometry g;
    for (const auto& m : j.at("modules")) {
      const std::string sub = m.at("subdet").get<std::string>();
      int S = 0;
      if (sub == "PXB")
        S = m.at("layer").get<int>();
      else if (sub == "PXF")
        S = 10 + m.at("disk").get<int>();
      if (!((S >= 1 && S <= 4) || (S >= 11 && S <= 14)))
        continue;
      Module mod;
      mod.detId = m.at("detId").get<uint32_t>();
      mod.surface = S;
      for (int k = 0; k < 3; ++k) {
        mod.o[k] = m.at("origin").at(k).get<double>();
        mod.ex[k] = m.at("ex").at(k).get<double>();
        mod.ey[k] = m.at("ey").at(k).get<double>();
        mod.ez[k] = m.at("ez").at(k).get<double>();
      }
      mod.hw = m.at("halfWidth").get<double>();
      mod.hl = m.at("halfLength").get<double>();
      g.add(mod);
    }
    g.finalize();
    return g;
  }

  std::array<bool, 23> activeOf(const std::string& mask) {
    std::array<bool, 23> a{};
    for (int k = 0; k < 4; ++k)
      if (mask.at(k) == 'A' || mask.at(k) == '1')
        a[1 + k] = a[11 + k] = true;
    return a;
  }

  struct Tally {
    std::map<std::string, long> n;
    std::map<std::string, std::vector<std::string>> ex;
    void bad(const std::string& what, const std::string& msg) {
      ++n[what];
      if (ex[what].size() < 5)
        ex[what].push_back(msg);
    }
  };
}  // namespace

int main(int argc, char** argv) {
  if (argc < 6) {
    std::cerr << "usage: refitV3NanoClosure nano.root reference.jsonl geometry.json manifest.json mask [nEvents] [knownBad]\n";
    return 2;
  }
  const std::string nano = argv[1], refPath = argv[2], geoPath = argv[3], manifest = argv[4], mask = argv[5];
  const int nEvents = argc > 6 ? std::stoi(argv[6]) : 10;
  const std::string knownBad = argc > 7 ? argv[7] : "";

  const Geometry geo = loadGeometry(geoPath);
  Payloads pay;
  pay.load(manifest);
  Config cfg;
  for (auto& c : cfg.use)
    for (auto& a : c)
      a = {{true, true}};
  cfg.knownBad = knownBad;
  const Engine eng(cfg, pay, geo);
  const auto active = activeOf(mask);

  // the reference, keyed by (run, lumi, event, track)
  std::map<std::tuple<long, long, long, int>, nlohmann::json> ref;
  {
    std::ifstream in(refPath);
    std::string line;
    while (std::getline(in, line)) {
      auto j = nlohmann::json::parse(line);
      ref[{j.at("run").get<long>(), j.at("lumi").get<long>(), j.at("event").get<long>(), j.at("track").get<int>()}] =
          std::move(j);
    }
  }

  TFile f(nano.c_str());
  TTreeReader r("Events", &f);
  TTreeReaderValue<UInt_t> run(r, "run"), lumi(r, "luminosityBlock");
  TTreeReaderValue<ULong64_t> evt(r, "event");
  const char* hp[5] = {"rInv", "phi", "tanL", "z0", "d0"};
  std::vector<std::unique_ptr<TTreeReaderArray<float>>> tp;
  for (const char* h : hp)
    tp.emplace_back(new TTreeReaderArray<float>(r, (std::string("L1TTrack_") + h).c_str()));
  TTreeReaderArray<float> tpt(r, "L1TTrack_pt");
  std::vector<std::unique_ptr<TTreeReaderArray<float>>> tc;  // upper triangle, row-major
  for (int i = 0; i < 5; ++i)
    for (int j = i; j < 5; ++j)
      tc.emplace_back(new TTreeReaderArray<float>(r, (std::string("L1TTrack_cov_") + hp[i] + "_" + hp[j]).c_str()));
  TTreeReaderArray<UInt_t> cdet(r, "L1TSmartPixelsCluster_detId");
  TTreeReaderArray<float> cx(r, "L1TSmartPixelsCluster_localX"), cy(r, "L1TSmartPixelsCluster_localY"),
      cca(r, "L1TSmartPixelsCluster_localCotAlpha"), ccb(r, "L1TSmartPixelsCluster_localCotBeta"),
      csx(r, "L1TSmartPixelsCluster_sigX"), csy(r, "L1TSmartPixelsCluster_sigY"), csa(r, "L1TSmartPixelsCluster_sigAlpha"),
      csb(r, "L1TSmartPixelsCluster_sigBeta");
  TTreeReaderArray<UChar_t> cha(r, "L1TSmartPixelsCluster_hasAlpha"), chb(r, "L1TSmartPixelsCluster_hasBeta"),
      csizx(r, "L1TSmartPixelsCluster_sizeX"), csizy(r, "L1TSmartPixelsCluster_sizeY");

  // kpt: the chunk median of pT |rInv| over the chunk's OT tracks (ORP.pt_constant)
  std::vector<double> ks;
  for (int e = 0; e < nEvents && r.Next(); ++e)
    for (size_t t = 0; t < tpt.GetSize(); ++t) {
      const double k = static_cast<double>(tpt[t]) * std::abs(static_cast<double>((*tp[0])[t]));
      if (std::isfinite(k) && k > 0)
        ks.push_back(k);
    }
  std::sort(ks.begin(), ks.end());
  const size_t nk = ks.size();
  const double kpt = nk % 2 ? ks[nk / 2] : 0.5 * (ks[nk / 2 - 1] + ks[nk / 2]);
  std::cout << "kpt " << std::setprecision(12) << kpt << " over " << nk << " tracks\n";

  Tally T;
  long nTracks = 0, nSteps = 0;
  double maxDa = 0., maxDC = 0., maxDchi2 = 0., maxDseed = 0.;
  r.Restart();
  for (int e = 0; e < nEvents && r.Next(); ++e) {
    EventClusters cls(geo.modules().size());
    for (size_t q = 0; q < cdet.GetSize(); ++q) {
      const int mi = geo.index(cdet[q]);
      if (mi < 0)
        continue;
      Cluster c;
      c.x = cx[q];
      c.y = cy[q];
      c.cotA = cca[q];
      c.cotB = ccb[q];
      c.sigX = csx[q];
      c.sigY = csy[q];
      c.sigA = csa[q];
      c.sigB = csb[q];
      c.hasA = cha[q] > 0 && std::abs(c.cotA) <= 12.0;
      c.hasB = chb[q] > 0 && std::abs(c.cotB) <= 12.0;
      c.sizeX = csizx[q];
      c.sizeY = csizy[q];
      c.index = static_cast<int>(q);
      cls[mi].push_back(c);
    }
    for (size_t t = 0; t < tpt.GetSize(); ++t) {
      Vec5 a;
      Mat5 C{};
      for (int i = 0; i < 5; ++i)
        a[i] = (*tp[i])[t];
      int n = 0;
      for (int i = 0; i < 5; ++i)
        for (int j = i; j < 5; ++j) {
          C[i][j] = C[j][i] = (*tc[n])[t];
          ++n;
        }
      const auto key = std::make_tuple(static_cast<long>(*run), static_cast<long>(*lumi), static_cast<long>(*evt),
                                       static_cast<int>(t));
      const auto it = ref.find(key);
      if (it == ref.end()) {
        T.bad("missing reference", "track " + std::to_string(t));
        continue;
      }
      const auto& R = it->second;
      ++nTracks;
      const TrackResult res = eng.run(a, C, kpt, active, cls);
      const std::string id = "ev " + std::to_string(*evt) + " trk " + std::to_string(t);
      if (res.refit != R.at("refit").get<bool>())
        T.bad("refit flag", id);
      if (res.nacc != R.at("nacc").get<int>())
        T.bad("nacc", id + " got " + std::to_string(res.nacc) + " want " + std::to_string(R.at("nacc").get<int>()));
      if (res.discGoing != R.at("disc_going").get<bool>())
        T.bad("disc_going", id);
      if (res.covOK) {
        for (int i = 0; i < 5; ++i)
          for (int j = 0; j < 5; ++j) {
            const double w = R.at("seed_C").at(i).at(j).get<double>();
            const double d = std::abs(res.seedC[i][j] - w) /
                             std::sqrt(std::abs(R.at("seed_C").at(i).at(i).get<double>() * R.at("seed_C").at(j).at(j).get<double>()));
            maxDseed = std::max(maxDseed, d);
            if (d > 1e-9)
              T.bad("seed C", id);
          }
      }
      const auto& steps = R.at("steps");
      bool firstDiff = true;
      for (size_t j = 0; j < steps.size(); ++j) {
        const auto& s = steps[j];
        const bool vis = s.at("visited").get<bool>();
        const StepRecord* m = j < res.steps.size() ? &res.steps[j] : nullptr;
        const bool mvis = m && m->visited;
        if (vis != mvis) {
          T.bad("visited", id + " step " + std::to_string(j));
          continue;
        }
        if (!vis)
          continue;
        ++nSteps;
        auto diff = [&](const std::string& what, const std::string& msg) {
          if (firstDiff)
            T.bad(what + " (first per track)", id + " step " + std::to_string(j) + " " + msg);
          else
            T.bad(what + " (later)", id + " step " + std::to_string(j) + " " + msg);
          firstDiff = false;
        };
        if (m->surface != s.at("surface").get<int>())
          diff("surface", "");
        if (m->ncand != s.at("ncand").get<int>())
          diff("ncand", std::to_string(m->ncand) + " vs " + std::to_string(s.at("ncand").get<int>()));
        if (m->accepted != s.at("accepted").get<bool>() || (m->accepted && m->cluster != s.at("cluster").get<int>()))
          diff("cluster", std::to_string(m->cluster) + " vs " + std::to_string(s.at("cluster").get<int>()));
        if (m->kchoice != s.at("kchoice").get<int>())
          diff("kchoice", std::to_string(m->kchoice) + " vs " + std::to_string(s.at("kchoice").get<int>()));
        if (m->accepted)
          for (int k = 0; k < 4; ++k) {
            const double w = s.at("chi2").at(k).get<double>();
            const double d = std::abs(m->chi2[k] - w) / std::max(1.0, std::abs(w));
            maxDchi2 = std::max(maxDchi2, d);
            if (d > 1e-6)
              diff("chi2", std::to_string(k) + " " + std::to_string(m->chi2[k]) + " vs " + std::to_string(w));
          }
      }
      for (size_t j = steps.size(); j < res.steps.size(); ++j)
        if (res.steps[j].visited)
          T.bad("visited", id + " extra step " + std::to_string(j));
      const char* fn[5] = {"rInv", "phi0", "tanL", "z0", "d0"};
      const auto& FC = R.at("final_C");
      for (int i = 0; i < 5; ++i) {
        const double w = R.at("final").at(fn[i]).get<double>();
        const double sg = std::sqrt(std::abs(FC.at(i).at(i).get<double>()));
        const double d = std::abs(res.a[i] - w) / (sg > 0 ? sg : 1.0);
        maxDa = std::max(maxDa, d);
        if (d > 1e-6)
          T.bad("final a", id + " " + fn[i] + " " + std::to_string(d) + " sigma");
        for (int k = 0; k < 5; ++k) {
          const double wc = FC.at(i).at(k).get<double>();
          const double nrm = std::sqrt(std::abs(FC.at(i).at(i).get<double>() * FC.at(k).at(k).get<double>()));
          const double dc = std::abs(res.C[i][k] - wc) / (nrm > 0 ? nrm : 1.0);
          maxDC = std::max(maxDC, dc);
          if (dc > 1e-6)
            T.bad("final C", id + " (" + std::to_string(i) + "," + std::to_string(k) + ") " + std::to_string(dc));
        }
      }
    }
  }
  std::cout << "tracks " << nTracks << ", visited steps " << nSteps << "\n"
            << "max |da|/sigma " << maxDa << ", max |dC|/sqrt(CiiCjj) " << maxDC << ", max rel dchi2 " << maxDchi2
            << ", max seed dC " << maxDseed << "\n";
  for (const auto& [k, v] : T.n) {
    std::cout << "MISMATCH " << k << ": " << v << "\n";
    for (const auto& x : T.ex[k])
      std::cout << "    " << x << "\n";
  }
  std::cout << (T.n.empty() ? "CLOSURE OK\n" : "CLOSURE FAILED\n");
  return T.n.empty() ? 0 : 1;
}
