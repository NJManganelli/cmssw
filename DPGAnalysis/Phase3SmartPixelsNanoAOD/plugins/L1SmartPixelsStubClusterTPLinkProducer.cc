// The TrackingParticles the OT stub table links to through its clusters: the owner
// (smartpixels::dominantTrackingParticle) of both clusters of every stub, once each -- the TPs
// behind L1TOTStub_tpIdxInner / tpIdxOuter, which also cover the stub clusters entering the
// SmartPixels track tpIdx. Input of the "linked" L1TTP selection
// (L1TrackingParticleTableProducer linkedTrackingParticles).

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/L1TrackTrigger/interface/TTStub.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "L1Trigger/Phase3SmartPixels/interface/SmartPixelsOTClusterTruth.h"
#include "SimDataFormats/Associations/interface/TTClusterAssociationMap.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"

#include <memory>
#include <set>

class L1SmartPixelsStubClusterTPLinkProducer : public edm::global::EDProducer<> {
public:
  explicit L1SmartPixelsStubClusterTPLinkProducer(const edm::ParameterSet& cfg)
      : stubsToken_(consumes<TTStubDetSetVec>(cfg.getParameter<edm::InputTag>("stubs"))),
        clusterTruthToken_(consumes<TTClusterAssociationMap<Ref_Phase2TrackerDigi_>>(
            cfg.getParameter<edm::InputTag>("clusterTruth"))),
        tpToken_(consumes<TrackingParticleCollection>(cfg.getParameter<edm::InputTag>("trackingParticles"))) {
    produces<TrackingParticleRefVector>();
  }

  void produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup&) const override {
    const auto& clusterTruth = iEvent.get(clusterTruthToken_);
    // Refs from the TP handle: a Ptr read back from a file carries no product getter.
    const auto tps = iEvent.getHandle(tpToken_);
    auto out = std::make_unique<TrackingParticleRefVector>();
    std::set<size_t> seen;
    for (const auto& det : iEvent.get(stubsToken_))
      for (const auto& stub : det)
        for (unsigned int member = 0; member < 2; ++member) {
          const auto tp = smartpixels::dominantTrackingParticle(clusterTruth, stub.clusterRef(member));
          if (tp.isNull() || !seen.insert(tp.key()).second)
            continue;
          if (tp.id() != tps.id())
            throw cms::Exception("L1SmartPixelsStubClusterTPLinkProducer")
                << "the cluster association links TrackingParticle product " << tp.id() << ", not "
                << tps.id() << " (trackingParticles)";
          out->push_back(TrackingParticleRef(tps, tp.key()));
        }
    iEvent.put(std::move(out));
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("stubs", edm::InputTag("TTStubsFromPhase2TrackerDigis", "StubAccepted"))
        ->setComment("the stubs of the OT stub table");
    desc.add<edm::InputTag>("clusterTruth", edm::InputTag("TTClusterAssociatorFromPixelDigis", "ClusterAccepted"))
        ->setComment("the OT stub table's cluster association");
    desc.add<edm::InputTag>("trackingParticles", edm::InputTag("mix", "MergedTrackTruth"))
        ->setComment("the TrackingParticles the cluster association links (the L1TTP table's input)");
    descriptions.addWithDefaultLabel(desc);
  }

private:
  const edm::EDGetTokenT<TTStubDetSetVec> stubsToken_;
  const edm::EDGetTokenT<TTClusterAssociationMap<Ref_Phase2TrackerDigi_>> clusterTruthToken_;
  const edm::EDGetTokenT<TrackingParticleCollection> tpToken_;
};

DEFINE_FWK_MODULE(L1SmartPixelsStubClusterTPLinkProducer);
