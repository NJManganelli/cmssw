// Produces jet <-> constituent association tables for L1 PF jets, in the
// style of the PFNano JetPFCands tables (flat-ntuple analogue of the
// podio/edm4hep one-to-many relation):
//   - a link table with one row per (jet, constituent) pair, carrying the
//     index of the jet in the jet table, the index of the candidate in the
//     candidate table, the constituent slot (position in the stored
//     constituent list, which is the order consumed by the L1TSC4NGJetID
//     tagger) and whether the slot is within the tagger constituent window
//   - optionally, an extension table on the candidate table with a backref
//     to the first jet clustering each candidate and the index of the
//     underlying TTTrack in the L1 track table
//
// Index validity relies on the referenced tables using the same source
// collections with an empty cut and native ordering.
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/L1TParticleFlow/interface/PFCandidate.h"
#include "DataFormats/L1TParticleFlow/interface/PFJet.h"
#include "DataFormats/L1TParticleFlow/interface/PFTrack.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"

class L1JetCandLinkTableProducer : public edm::stream::EDProducer<> {
public:
  explicit L1JetCandLinkTableProducer(const edm::ParameterSet& cfg)
      : jetsToken_(consumes<edm::View<l1t::PFJet>>(cfg.getParameter<edm::InputTag>("jets"))),
        candsToken_(consumes<l1t::PFCandidateCollection>(cfg.getParameter<edm::InputTag>("cands"))),
        tracksToken_(consumes<std::vector<l1t::PFTrack::L1TTTrackType>>(cfg.getParameter<edm::InputTag>("tracks"))),
        linkTableName_(cfg.getParameter<std::string>("linkTableName")),
        candTableName_(cfg.getParameter<std::string>("candTableName")),
        jetTableName_(cfg.getParameter<std::string>("jetTableName")),
        trackTableName_(cfg.getParameter<std::string>("trackTableName")),
        maxConstituentsTagger_(cfg.getParameter<uint32_t>("maxConstituentsTagger")),
        writeCandExtension_(cfg.getParameter<bool>("writeCandExtension")) {
    produces<nanoaod::FlatTable>("link");
    if (writeCandExtension_) {
      produces<nanoaod::FlatTable>("cands");
    }
  }

  void produce(edm::Event& iEvent, const edm::EventSetup&) override {
    const auto& jets = iEvent.get(jetsToken_);
    auto candsHandle = iEvent.getHandle(candsToken_);
    auto tracksHandle = iEvent.getHandle(tracksToken_);

    std::vector<int32_t> linkJetIdx, linkCandIdx;
    std::vector<int16_t> linkSlot;
    std::vector<bool> linkInTagger;

    const size_t nCands = candsHandle.isValid() ? candsHandle->size() : 0;
    std::vector<int32_t> candJetIdx(nCands, -1);

    // Sentinel convention (see DPGAnalysis/Phase3SmartPixelsNanoAOD/README.md):
    //   -2 = PRODUCT-level failure (the L1PuppiCand collection is not available in
    //        this event, so no jet-constituent link can resolve for ANY row);
    //   -1 = ELEMENT-level no-match (the collection IS available but this specific
    //        constituent Ptr does not resolve into it).
    const bool candsProductAvailable = candsHandle.isValid();

    for (size_t ijet = 0; ijet < jets.size(); ++ijet) {
      const auto& constituents = jets[ijet].constituents();
      for (size_t islot = 0; islot < constituents.size(); ++islot) {
        const auto& ptr = constituents[islot];
        int32_t candIdx = candsProductAvailable ? -1 : -2;
        if (candsProductAvailable && ptr.isNonnull() && ptr.id() == candsHandle.id()) {
          candIdx = ptr.key();
          if (candJetIdx[ptr.key()] < 0) {
            candJetIdx[ptr.key()] = ijet;
          }
        }
        linkJetIdx.push_back(ijet);
        linkCandIdx.push_back(candIdx);
        linkSlot.push_back(islot);
        linkInTagger.push_back(islot < maxConstituentsTagger_);
      }
    }

    auto linkTable = std::make_unique<nanoaod::FlatTable>(linkJetIdx.size(), linkTableName_, false, false);
    linkTable->addColumn<int32_t>("jetIdx", linkJetIdx, "index of the jet in the " + jetTableName_ + " table");
    linkTable->addColumn<int32_t>(
        "candIdx", linkCandIdx,
        "index of the constituent in the " + candTableName_ +
            " table (-1 if this constituent does not resolve; -2 if the " + candTableName_ +
            " collection is unavailable in this event)");
    linkTable->addColumn<int16_t>(
        "slot", linkSlot, "position in the stored constituent list (order consumed by the jet tagger)");
    linkTable->addColumn<bool>(
        "inTagger", linkInTagger, "constituent is within the tagger window of the first N constituents");
    linkTable->setDoc("Association between " + jetTableName_ + " jets and their " + candTableName_ + " constituents");
    iEvent.put(std::move(linkTable), "link");

    if (writeCandExtension_) {
      // Sentinel convention (see DPGAnalysis/Phase3SmartPixelsNanoAOD/README.md):
      //   -2 = PRODUCT-level failure: the PFTrack (or the TTTrack) collection needed
      //        to resolve ANY candidate->track link is unavailable in this event, so
      //        linking could not run for anyone. isNonnull() only checks that a Ref
      //        carries a ProductID+key; the referenced PFTrack collection may not be
      //        persisted (e.g. PU files that keep l1tLayer1:PuppiRegional but drop the
      //        PFTrack collection it points to). isAvailable() detects that.
      //   -1 = ELEMENT-level no-match: the products ARE available but this specific
      //        candidate has no resolvable track ref (e.g. a neutral candidate).
      // Product-level probe: if any nonnull pfTrack ref is unavailable, the whole
      // PFTrack product is missing. Also treat an unavailable TTTrack collection as a
      // product-level failure (nothing can index into it).
      bool pfTrackProductAvailable = true;
      if (candsHandle.isValid()) {
        for (size_t icand = 0; icand < nCands; ++icand) {
          const auto& trkRef = (*candsHandle)[icand].pfTrack();
          if (trkRef.isNonnull() && !trkRef.isAvailable()) {
            pfTrackProductAvailable = false;
            break;
          }
        }
      }
      const bool trackLinkResolvable = candsHandle.isValid() && tracksHandle.isValid() && pfTrackProductAvailable;

      std::vector<int32_t> candTrackIdx(nCands, trackLinkResolvable ? -1 : -2);
      if (trackLinkResolvable) {
        for (size_t icand = 0; icand < nCands; ++icand) {
          const auto& trkRef = (*candsHandle)[icand].pfTrack();
          if (trkRef.isNonnull() && trkRef.isAvailable()) {
            const auto& ttRef = trkRef->track();
            if (ttRef.isNonnull() && ttRef.id() == tracksHandle.id()) {
              candTrackIdx[icand] = ttRef.key();
            }
          }
        }
      }
      auto candTable = std::make_unique<nanoaod::FlatTable>(nCands, candTableName_, false, true);
      candTable->addColumn<int32_t>(
          "jetIdx", candJetIdx, "index of the first " + jetTableName_ + " jet clustering this candidate (-1 if none)");
      candTable->addColumn<int32_t>(
          "l1TrackIdx", candTrackIdx,
          "index of the underlying track in the " + trackTableName_ +
              " table (-1 if this candidate has no resolvable track; -2 if the track/PFTrack product is unavailable)");
      iEvent.put(std::move(candTable), "cands");
    }
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("jets", edm::InputTag("l1tSC4NGJetProducer", "l1tSC4NGJets"));
    desc.add<edm::InputTag>("cands", edm::InputTag("l1tLayer2Deregionizer", "Puppi"));
    desc.add<edm::InputTag>("tracks", edm::InputTag("l1tTTTracksFromTrackletEmulation", "Level1TTTracks"));
    desc.add<std::string>("linkTableName", "L1SC4NGJetCands");
    desc.add<std::string>("candTableName", "L1PuppiCand");
    desc.add<std::string>("jetTableName", "L1puppiJetSC4NG");
    desc.add<std::string>("trackTableName", "L1TTrack");
    desc.add<uint32_t>("maxConstituentsTagger", 16);
    desc.add<bool>("writeCandExtension", true);
    descriptions.addWithDefaultLabel(desc);
  }

private:
  const edm::EDGetTokenT<edm::View<l1t::PFJet>> jetsToken_;
  const edm::EDGetTokenT<l1t::PFCandidateCollection> candsToken_;
  const edm::EDGetTokenT<std::vector<l1t::PFTrack::L1TTTrackType>> tracksToken_;
  const std::string linkTableName_;
  const std::string candTableName_;
  const std::string jetTableName_;
  const std::string trackTableName_;
  const uint32_t maxConstituentsTagger_;
  const bool writeCandExtension_;
};

DEFINE_FWK_MODULE(L1JetCandLinkTableProducer);
