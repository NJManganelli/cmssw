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

#include <string>
#include <vector>

class L1SmartPixelsOTTrackLinkTableProducer : public edm::global::EDProducer<> {
public:
  using L1TrackCollection = std::vector<TTTrack<Ref_Phase2TrackerDigi_>>;

  explicit L1SmartPixelsOTTrackLinkTableProducer(const edm::ParameterSet& cfg)
      : tracksToken_(consumes<L1TrackCollection>(cfg.getParameter<edm::InputTag>("src"))),
        linksToken_(consumes<std::vector<TTTrackRef>>(cfg.getParameter<edm::InputTag>("src"))),
        otTracksToken_(consumes<L1TrackCollection>(cfg.getParameter<edm::InputTag>("otTracks"))),
        trackTableName_(cfg.getParameter<std::string>("trackTableName")),
        otTrackTableName_(cfg.getParameter<std::string>("otTrackTableName")) {
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
    iEvent.put(std::move(table));
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("src", edm::InputTag("l1tSmartPixelsTrackProducer", "Level1TTTracks"))
        ->setComment("SmartPixels tracks; the OT-track links are read from the same tag");
    desc.add<edm::InputTag>("otTracks", edm::InputTag("l1tTTTracksFromTrackletEmulation", "Level1TTTracks"))
        ->setComment("source collection of the otTrackTableName table");
    desc.add<std::string>("trackTableName", "L1TSmartPixelsTrack");
    desc.add<std::string>("otTrackTableName", "L1TTrack");
    descriptions.addWithDefaultLabel(desc);
  }

private:
  const edm::EDGetTokenT<L1TrackCollection> tracksToken_;
  const edm::EDGetTokenT<std::vector<TTTrackRef>> linksToken_;
  const edm::EDGetTokenT<L1TrackCollection> otTracksToken_;
  const std::string trackTableName_;
  const std::string otTrackTableName_;
};

DEFINE_FWK_MODULE(L1SmartPixelsOTTrackLinkTableProducer);
