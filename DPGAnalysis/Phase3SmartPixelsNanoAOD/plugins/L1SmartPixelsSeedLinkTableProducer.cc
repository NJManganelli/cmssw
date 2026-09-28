// Extension table on a SmartPixels track table: l1TrackIdx, the row of the OT
// track (in seedTableName, whose source collection is `seedTracks`) each
// SmartPixels track was seeded from. Read from the per-track seed links the
// L1SmartPixelsTrackProducer puts next to its tracks (same InputTag, type
// std::vector<TTTrackRef>), so the link is explicit rather than inferred from
// the producer's current 1:1 ordering.
//
//   >= 0  row in seedTableName
//     -1  the track has no OT seed (e.g. a future joint IT+OT seeding producer)
//     -2  the seedTracks product is unavailable in this event
//
// A link into any other collection than seedTracks is a configuration error
// (the seed table would not describe the seeds) and throws.
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

class L1SmartPixelsSeedLinkTableProducer : public edm::global::EDProducer<> {
public:
  using L1TrackCollection = std::vector<TTTrack<Ref_Phase2TrackerDigi_>>;

  explicit L1SmartPixelsSeedLinkTableProducer(const edm::ParameterSet& cfg)
      : tracksToken_(consumes<L1TrackCollection>(cfg.getParameter<edm::InputTag>("src"))),
        linksToken_(consumes<std::vector<TTTrackRef>>(cfg.getParameter<edm::InputTag>("src"))),
        seedTracksToken_(consumes<L1TrackCollection>(cfg.getParameter<edm::InputTag>("seedTracks"))),
        trackTableName_(cfg.getParameter<std::string>("trackTableName")),
        seedTableName_(cfg.getParameter<std::string>("seedTableName")) {
    produces<nanoaod::FlatTable>();
  }

  void produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup&) const override {
    const auto& tracks = iEvent.get(tracksToken_);
    const auto& links = iEvent.get(linksToken_);
    if (links.size() != tracks.size())
      throw cms::Exception("SmartPixelsSeedLinkSize")
          << trackTableName_ << ": " << links.size() << " seed links for " << tracks.size() << " tracks";
    const auto seeds = iEvent.getHandle(seedTracksToken_);

    std::vector<int32_t> idx(links.size(), seeds.isValid() ? -1 : -2);
    if (seeds.isValid()) {
      for (size_t i = 0; i < links.size(); ++i) {
        const auto& ref = links[i];
        if (ref.isNull())
          continue;
        if (ref.id() != seeds.id())
          throw cms::Exception("SmartPixelsSeedLinkMismatch")
              << trackTableName_ << " row " << i << " is seeded from product " << ref.id() << ", not from the "
              << seedTableName_ << " source collection (" << seeds.id() << ")";
        idx[i] = ref.key();
      }
    }

    auto table = std::make_unique<nanoaod::FlatTable>(links.size(), trackTableName_, false, true);
    table->addColumn<int32_t>("l1TrackIdx",
                              idx,
                              "row of the " + seedTableName_ + " track this SmartPixels track was seeded from (-1 if "
                                  "not seeded from an " + seedTableName_ +
                                  " track, e.g. joint IT+OT seeding; -2 if the " + seedTableName_ +
                                  " product is unavailable)");
    iEvent.put(std::move(table));
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("src", edm::InputTag("l1tSmartPixelsTrackProducer", "Level1TTTracks"))
        ->setComment("SmartPixels tracks; the seed links are read from the same tag");
    desc.add<edm::InputTag>("seedTracks", edm::InputTag("l1tTTTracksFromTrackletEmulation", "Level1TTTracks"))
        ->setComment("source collection of the seedTableName table");
    desc.add<std::string>("trackTableName", "L1TSmartPixelsTrack");
    desc.add<std::string>("seedTableName", "L1TTrack");
    descriptions.addWithDefaultLabel(desc);
  }

private:
  const edm::EDGetTokenT<L1TrackCollection> tracksToken_;
  const edm::EDGetTokenT<std::vector<TTTrackRef>> linksToken_;
  const edm::EDGetTokenT<L1TrackCollection> seedTracksToken_;
  const std::string trackTableName_;
  const std::string seedTableName_;
};

DEFINE_FWK_MODULE(L1SmartPixelsSeedLinkTableProducer);
