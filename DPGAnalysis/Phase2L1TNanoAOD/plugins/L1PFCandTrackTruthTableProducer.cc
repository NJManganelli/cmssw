// Emits an extension table on a L1 PF candidate table with track MC truth
// obtained by following the candidate's PFCandidate->PFTrack->TTTrack ref
// chain into the TTTrackAssociationMap. This makes track genuine/fake truth
// available in candidate-only tiers (L1PFNano) that do not carry the track
// tables themselves. Columns default when the candidate has no track or the
// associator did not run.
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "DataFormats/L1TParticleFlow/interface/PFCandidate.h"
#include "DataFormats/L1TParticleFlow/interface/PFTrack.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "SimDataFormats/Associations/interface/TTTrackAssociationMap.h"

class L1PFCandTrackTruthTableProducer : public edm::stream::EDProducer<> {
public:
  using L1Track = l1t::PFTrack::L1TTTrackType;
  using L1TrackCollection = std::vector<L1Track>;
  using AssocMap = TTTrackAssociationMap<Ref_Phase2TrackerDigi_>;

  explicit L1PFCandTrackTruthTableProducer(const edm::ParameterSet& cfg)
      : candsToken_(consumes<l1t::PFCandidateCollection>(cfg.getParameter<edm::InputTag>("cands"))),
        tracksToken_(consumes<L1TrackCollection>(cfg.getParameter<edm::InputTag>("tracks"))),
        assocToken_(consumes<AssocMap>(cfg.getParameter<edm::InputTag>("truthAssociation"))),
        candTableName_(cfg.getParameter<std::string>("candTableName")) {
    produces<nanoaod::FlatTable>();
  }

  void produce(edm::Event& iEvent, const edm::EventSetup&) override {
    auto candsHandle = iEvent.getHandle(candsToken_);
    auto tracksHandle = iEvent.getHandle(tracksToken_);
    auto assocHandle = iEvent.getHandle(assocToken_);

    const size_t nCands = candsHandle.isValid() ? candsHandle->size() : 0;
    std::vector<bool> genuine(nCands, false), loose(nCands, false), combinatoric(nCands, false),
        unknown(nCands, false);
    std::vector<float> tpPt(nCands, -1.f);
    std::vector<bool> tpFromHardInteraction(nCands, false);

    // Sentinel convention (see DPGAnalysis/Phase3SmartPixelsNanoAOD/README.md):
    //   -2 = PRODUCT-level failure: a collection needed to resolve the
    //        PFCandidate->PFTrack->TTTrack->truth chain is unavailable in this event
    //        (the PFTrack collection is dropped so its refs dangle, or the TTTrack /
    //        TTTrackAssociationMap products are absent), so truth could not run for
    //        ANY candidate. This replaces the former hard-deref of unstored PFTrack
    //        refs (which threw ProductNotFound and forced the table to be dropped);
    //        the table now self-guards and can always be scheduled.
    //   -1 = ELEMENT-level no-match: the products ARE available but this specific
    //        candidate has no resolvable track ref (e.g. a neutral candidate).
    //    0 = the candidate's underlying track was resolved; the truth columns carry
    //        the association result for it.
    // Product-level probe for the PFTrack collection: isNonnull() only means the Ref
    // carries a ProductID+key; isAvailable() detects whether the referenced PFTrack
    // product is actually present. A single dangling ref => the product is missing.
    bool pfTrackProductAvailable = true;
    if (candsHandle.isValid()) {
      for (size_t i = 0; i < nCands; ++i) {
        const auto& trkRef = (*candsHandle)[i].pfTrack();
        if (trkRef.isNonnull() && !trkRef.isAvailable()) {
          pfTrackProductAvailable = false;
          break;
        }
      }
    }
    const bool truthResolvable =
        candsHandle.isValid() && tracksHandle.isValid() && assocHandle.isValid() && pfTrackProductAvailable;

    std::vector<int32_t> truthStatus(nCands, truthResolvable ? -1 : -2);

    if (truthResolvable) {
      for (size_t i = 0; i < nCands; ++i) {
        const auto& trkRef = (*candsHandle)[i].pfTrack();
        if (trkRef.isNull() || !trkRef.isAvailable())
          continue;
        const auto& ttRef = trkRef->track();
        if (ttRef.isNull() || ttRef.id() != tracksHandle.id())
          continue;
        edm::Ptr<L1Track> ptr(tracksHandle, ttRef.key());
        truthStatus[i] = 0;
        genuine[i] = assocHandle->isGenuine(ptr);
        loose[i] = assocHandle->isLooselyGenuine(ptr);
        combinatoric[i] = assocHandle->isCombinatoric(ptr);
        unknown[i] = assocHandle->isUnknown(ptr);
        const auto& tp = assocHandle->findTrackingParticlePtr(ptr);
        if (tp.isNonnull()) {
          tpPt[i] = tp->pt();
          tpFromHardInteraction[i] = (tp->eventId().event() == 0) && (tp->eventId().bunchCrossing() == 0);
        }
      }
    }

    auto table = std::make_unique<nanoaod::FlatTable>(nCands, candTableName_, false, true);
    table->addColumn<int32_t>(
        "trkTruthStatus", truthStatus,
        "candidate->track truth resolution: 0 = resolved (truth columns valid); -1 = candidate has no "
        "resolvable track ref; -2 = the PFTrack/TTTrack/association product is unavailable this event");
    table->addColumn<bool>("trkGenuine", genuine, "underlying track is genuine (TTTrackAssociationMap)");
    table->addColumn<bool>("trkLooselyGenuine", loose, "underlying track is loosely genuine");
    table->addColumn<bool>("trkCombinatoric", combinatoric, "underlying track is a combinatoric fake");
    table->addColumn<bool>("trkUnknown", unknown, "underlying track has no truth association");
    table->addColumn<float>("trkTpPt", tpPt, "pt of the TrackingParticle matched to the underlying track (-1 if none)");
    table->addColumn<bool>("trkTpFromHardInteraction",
                           tpFromHardInteraction,
                           "matched TrackingParticle is from the hard interaction");
    iEvent.put(std::move(table));
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("cands", edm::InputTag("l1tLayer2Deregionizer", "Puppi"));
    desc.add<edm::InputTag>("tracks", edm::InputTag("l1tTTTracksFromTrackletEmulation", "Level1TTTracks"));
    desc.add<edm::InputTag>("truthAssociation", edm::InputTag("TTTrackAssociatorFromPixelDigis", "Level1TTTracks"));
    desc.add<std::string>("candTableName", "L1PuppiCand");
    descriptions.addWithDefaultLabel(desc);
  }

private:
  const edm::EDGetTokenT<l1t::PFCandidateCollection> candsToken_;
  const edm::EDGetTokenT<L1TrackCollection> tracksToken_;
  const edm::EDGetTokenT<AssocMap> assocToken_;
  const std::string candTableName_;
};

DEFINE_FWK_MODULE(L1PFCandTrackTruthTableProducer);
