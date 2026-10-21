// Extension table on a SmartPixels track table: l1TrackIdx, the row of the OT
// track (in otTrackTableName, whose source collection is `otTracks`) each
// SmartPixels track was built from. Read from the per-track OT-track links the
// L1SmartPixelsTrackProducer puts next to its tracks (same InputTag, type
// std::vector<TTTrackRef>), so the link is explicit rather than inferred from
// the producer's current 1:1 ordering.
//
//   >= 0  row in otTrackTableName
//     -1  the track was not built from an OT track (e.g. a future joint IT+OT seeding producer)
//     -2  the otTracks product is unavailable in this event
//
// A link into any other collection than otTracks is a configuration error
// (the OT track table would not describe those tracks) and throws.
//
// With otTruthTable (the OT track table's truth extension, L1TrackTruthTableProducer), the OT track's
// matched-TP truth is copied onto the SmartPixels track ("otlink": the SmartPixels track takes its OT
// track's TP), column by column and bit for bit: ot<Column> (otTp_pt, otTpIdx, otGenuine, ...). A track
// without an OT link gets the OT table's unmatched defaults.
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"

#include <cctype>
#include <string>
#include <type_traits>
#include <vector>

namespace {
  // the OT truth columns copied onto the SmartPixels track, with the OT producer's unmatched default
  struct TruthColumn {
    const char* name;
    nanoaod::FlatTable::ColumnType type;
    float fDefault;
    int32_t iDefault;
  };
  using CT = nanoaod::FlatTable::ColumnType;
  constexpr TruthColumn kTruthColumns[] = {
      {"tpIdx", CT::Int32, 0.f, -1},         {"tpPdgId", CT::Int32, 0.f, 0},
      {"genuine", CT::Bool, 0.f, 0},         {"looselyGenuine", CT::Bool, 0.f, 0},
      {"combinatoric", CT::Bool, 0.f, 0},    {"unknown", CT::Bool, 0.f, 0},
      {"tpFromHardInteraction", CT::Bool, 0.f, 0},
      {"tp_pt", CT::Float, -1.f, 0},         {"tp_eta", CT::Float, 0.f, 0},
      {"tp_phi0", CT::Float, -999.f, 0},     {"tp_tanL", CT::Float, -999.f, 0},
      {"tp_charge", CT::Float, -999.f, 0},   {"tp_rInv", CT::Float, -999.f, 0},
      {"tp_d0", CT::Float, -999.f, 0},       {"tp_z0", CT::Float, -999.f, 0},
  };

  std::string otName(const char* c) {
    std::string n(c);
    n[0] = static_cast<char>(std::toupper(n[0]));
    return "ot" + n;
  }
}  // namespace

class L1SmartPixelsOTTrackLinkTableProducer : public edm::global::EDProducer<> {
public:
  using L1TrackCollection = std::vector<TTTrack<Ref_Phase2TrackerDigi_>>;

  explicit L1SmartPixelsOTTrackLinkTableProducer(const edm::ParameterSet& cfg)
      : tracksToken_(consumes<L1TrackCollection>(cfg.getParameter<edm::InputTag>("src"))),
        linksToken_(consumes<std::vector<TTTrackRef>>(cfg.getParameter<edm::InputTag>("src"))),
        otTracksToken_(consumes<L1TrackCollection>(cfg.getParameter<edm::InputTag>("otTracks"))),
        trackTableName_(cfg.getParameter<std::string>("trackTableName")),
        otTrackTableName_(cfg.getParameter<std::string>("otTrackTableName")) {
    const auto truthTag = cfg.getParameter<edm::InputTag>("otTruthTable");
    if (!truthTag.label().empty())
      otTruthToken_ = consumes<nanoaod::FlatTable>(truthTag);
    produces<nanoaod::FlatTable>();
  }

  void produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup&) const override {
    const auto& tracks = iEvent.get(tracksToken_);
    const auto& links = iEvent.get(linksToken_);
    if (links.size() != tracks.size())
      throw cms::Exception("SmartPixelsOTTrackLinkSize")
          << trackTableName_ << ": " << links.size() << " OT-track links for " << tracks.size() << " tracks";
    const auto otTracks = iEvent.getHandle(otTracksToken_);

    std::vector<int32_t> idx(links.size(), otTracks.isValid() ? -1 : -2);
    if (otTracks.isValid()) {
      for (size_t i = 0; i < links.size(); ++i) {
        const auto& ref = links[i];
        if (ref.isNull())
          continue;
        if (ref.id() != otTracks.id())
          throw cms::Exception("SmartPixelsOTTrackLinkMismatch")
              << trackTableName_ << " row " << i << " is built from product " << ref.id() << ", not from the "
              << otTrackTableName_ << " source collection (" << otTracks.id() << ")";
        idx[i] = ref.key();
      }
    }

    auto table = std::make_unique<nanoaod::FlatTable>(links.size(), trackTableName_, false, true);
    table->addColumn<int32_t>("l1TrackIdx",
                              idx,
                              "row of the " + otTrackTableName_ + " track this SmartPixels track was built from (-1 if "
                                  "not built from an " + otTrackTableName_ +
                                  " track, e.g. joint IT+OT seeding; -2 if the " + otTrackTableName_ +
                                  " product is unavailable)");
    if (!otTruthToken_.isUninitialized())
      addOTTruth(*table, iEvent.get(otTruthToken_), idx);
    iEvent.put(std::move(table));
  }

  // The OT truth columns at each SmartPixels track's OT row, copied without conversion (the stored values,
  // mantissa reduction included).
  void addOTTruth(nanoaod::FlatTable& table, const nanoaod::FlatTable& ot, const std::vector<int32_t>& idx) const {
    if (ot.name() != otTrackTableName_)
      throw cms::Exception("SmartPixelsOTTruthTable")
          << "otTruthTable describes " << ot.name() << ", not " << otTrackTableName_;
    for (const auto& c : kTruthColumns) {
      const int col = ot.columnIndex(c.name);
      if (col < 0 || ot.columnType(col) != c.type)
        throw cms::Exception("SmartPixelsOTTruthTable")
            << otTrackTableName_ << " truth: column " << c.name << " missing or of another type";
      const std::string doc = "otlink truth: " + otTrackTableName_ + "_" + c.name + " of the OT track at l1TrackIdx (" +
                              ot.columnDoc(col) + ")";
      const auto pick = [&](auto data, auto dflt) {
        std::vector<std::remove_cv_t<std::remove_reference_t<decltype(data[0])>>> v(idx.size(), dflt);
        for (size_t i = 0; i < idx.size(); ++i)
          if (idx[i] >= 0) {
            if (static_cast<unsigned>(idx[i]) >= ot.nRows())
              throw cms::Exception("SmartPixelsOTTruthTable")
                  << trackTableName_ << " row " << i << " links OT row " << idx[i] << " of " << ot.nRows();
            v[i] = data[idx[i]];
          }
        return v;
      };
      if (c.type == CT::Float)
        table.addColumn<float>(otName(c.name), pick(ot.columnData<float>(col), c.fDefault), doc);
      else if (c.type == CT::Int32)
        table.addColumn<int32_t>(otName(c.name), pick(ot.columnData<int32_t>(col), c.iDefault), doc);
      else {
        const auto v = pick(ot.columnData<bool>(col), uint8_t(c.iDefault));
        table.addColumn<bool>(otName(c.name), std::vector<bool>(v.begin(), v.end()), doc);
      }
    }
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("src", edm::InputTag("l1tSmartPixelsTrackProducer", "Level1TTTracks"))
        ->setComment("SmartPixels tracks; the OT-track links are read from the same tag");
    desc.add<edm::InputTag>("otTracks", edm::InputTag("l1tTTTracksFromTrackletEmulation", "Level1TTTracks"))
        ->setComment("source collection of the otTrackTableName table");
    desc.add<std::string>("trackTableName", "L1TSmartPixelsTrack");
    desc.add<std::string>("otTrackTableName", "L1TTrack");
    desc.add<edm::InputTag>("otTruthTable", edm::InputTag(""))
        ->setComment("the otTrackTableName truth extension (L1TrackTruthTableProducer); empty = no otlink truth columns");
    descriptions.addWithDefaultLabel(desc);
  }

private:
  const edm::EDGetTokenT<L1TrackCollection> tracksToken_;
  const edm::EDGetTokenT<std::vector<TTTrackRef>> linksToken_;
  const edm::EDGetTokenT<L1TrackCollection> otTracksToken_;
  const std::string trackTableName_;
  const std::string otTrackTableName_;
  edm::EDGetTokenT<nanoaod::FlatTable> otTruthToken_;
};

DEFINE_FWK_MODULE(L1SmartPixelsOTTrackLinkTableProducer);
