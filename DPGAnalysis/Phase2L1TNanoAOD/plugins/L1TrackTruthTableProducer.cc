// Emits an extension table on a L1 track table with:
//   - nStubs (always available; needed by the track-quality GBDT)
//   - MC truth columns from the TTTrackAssociationMap when present:
//     genuine / looselyGenuine / combinatoric / unknown flags and matched
//     TrackingParticle pt/eta/pdgId (defaults when the associator did not
//     run, e.g. data or central workflows without the truth sequence).
//
// Index/row alignment relies on the track table using the same source
// collection with an empty cut and native ordering.
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "SimDataFormats/Associations/interface/TTTrackAssociationMap.h"

class L1TrackTruthTableProducer : public edm::stream::EDProducer<> {
public:
  using L1Track = TTTrack<Ref_Phase2TrackerDigi_>;
  using L1TrackCollection = std::vector<L1Track>;
  using AssocMap = TTTrackAssociationMap<Ref_Phase2TrackerDigi_>;

  explicit L1TrackTruthTableProducer(const edm::ParameterSet& cfg)
      : tracksToken_(consumes<L1TrackCollection>(cfg.getParameter<edm::InputTag>("tracks"))),
        assocToken_(consumes<AssocMap>(cfg.getParameter<edm::InputTag>("truthAssociation"))),
        trackTableName_(cfg.getParameter<std::string>("trackTableName")) {
    produces<nanoaod::FlatTable>();
  }

  void produce(edm::Event& iEvent, const edm::EventSetup&) override {
    auto tracksHandle = iEvent.getHandle(tracksToken_);
    auto assocHandle = iEvent.getHandle(assocToken_);

    const size_t nTracks = tracksHandle.isValid() ? tracksHandle->size() : 0;
    std::vector<uint8_t> nStubs(nTracks, 0);
    std::vector<bool> genuine(nTracks, false), loose(nTracks, false), combinatoric(nTracks, false),
        unknown(nTracks, false);
    std::vector<float> tpPt(nTracks, -1.f), tpEta(nTracks, 0.f);
    std::vector<int32_t> tpPdgId(nTracks, 0);
    std::vector<bool> tpFromHardInteraction(nTracks, false);

    for (size_t i = 0; i < nTracks; ++i) {
      const auto& track = (*tracksHandle)[i];
      nStubs[i] = track.getStubRefs().size();
      if (assocHandle.isValid()) {
        edm::Ptr<L1Track> ptr(tracksHandle, i);
        genuine[i] = assocHandle->isGenuine(ptr);
        loose[i] = assocHandle->isLooselyGenuine(ptr);
        combinatoric[i] = assocHandle->isCombinatoric(ptr);
        unknown[i] = assocHandle->isUnknown(ptr);
        const auto& tp = assocHandle->findTrackingParticlePtr(ptr);
        if (tp.isNonnull()) {
          tpPt[i] = tp->pt();
          tpEta[i] = tp->eta();
          tpPdgId[i] = tp->pdgId();
          tpFromHardInteraction[i] = (tp->eventId().event() == 0) && (tp->eventId().bunchCrossing() == 0);
        }
      }
    }

    auto table = std::make_unique<nanoaod::FlatTable>(nTracks, trackTableName_, false, true);
    table->addColumn<uint8_t>("nStubs", nStubs, "number of stubs on the track");
    table->addColumn<bool>("genuine", genuine, "all stubs from the same TrackingParticle (TTTrackAssociationMap)");
    table->addColumn<bool>("looselyGenuine", loose, "at most one stub unmatched (TTTrackAssociationMap)");
    table->addColumn<bool>("combinatoric", combinatoric, "combinatoric fake (TTTrackAssociationMap)");
    table->addColumn<bool>("unknown", unknown, "no truth association (TTTrackAssociationMap)");
    table->addColumn<float>("tpPt", tpPt, "matched TrackingParticle pt (-1 if none)");
    table->addColumn<float>("tpEta", tpEta, "matched TrackingParticle eta");
    table->addColumn<int32_t>("tpPdgId", tpPdgId, "matched TrackingParticle pdgId");
    table->addColumn<bool>(
        "tpFromHardInteraction", tpFromHardInteraction, "matched TrackingParticle from the hard interaction");
    iEvent.put(std::move(table));
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("tracks", edm::InputTag("l1tTTTracksFromTrackletEmulation", "Level1TTTracks"));
    desc.add<edm::InputTag>("truthAssociation", edm::InputTag("TTTrackAssociatorFromPixelDigis", "Level1TTTracks"));
    desc.add<std::string>("trackTableName", "L1TTrack");
    descriptions.addWithDefaultLabel(desc);
  }

private:
  const edm::EDGetTokenT<L1TrackCollection> tracksToken_;
  const edm::EDGetTokenT<AssocMap> assocToken_;
  const std::string trackTableName_;
};

DEFINE_FWK_MODULE(L1TrackTruthTableProducer);
