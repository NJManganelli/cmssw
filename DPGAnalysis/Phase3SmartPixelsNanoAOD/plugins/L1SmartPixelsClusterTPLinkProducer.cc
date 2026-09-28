// The TrackingParticles the IT cluster table links to: the dominant charge contributor of
// every SmartPixels rec hit (the TP behind L1TSmartPixelsCluster_tpIdx), once each. Input
// of the "linked" L1TTP selection (L1TrackingParticleTableProducer linkedTrackingParticles).

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "DataFormats/Phase3SmartPixels/interface/SmartPixelsRecHitTruth.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"

#include <memory>
#include <set>

class L1SmartPixelsClusterTPLinkProducer : public edm::global::EDProducer<> {
public:
  explicit L1SmartPixelsClusterTPLinkProducer(const edm::ParameterSet& cfg)
      : truthToken_(consumes<SmartPixelsRecHitTruthCollection>(cfg.getParameter<edm::InputTag>("smartPixelsRecHits"))) {
    produces<TrackingParticleRefVector>();
  }

  void produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup&) const override {
    auto out = std::make_unique<TrackingParticleRefVector>();
    std::set<size_t> seen;
    for (const auto& det : iEvent.get(truthToken_))
      for (const auto& truth : det)
        if (truth.hasTp() && seen.insert(truth.dominantTp().key()).second)
          out->push_back(truth.dominantTp());
    iEvent.put(std::move(out));
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("smartPixelsRecHits", edm::InputTag("spixSmartPixelsRecHits"))
        ->setComment("SmartPixelsRecHitTruth of the rec hits the cluster table is built from");
    descriptions.addWithDefaultLabel(desc);
  }

private:
  const edm::EDGetTokenT<SmartPixelsRecHitTruthCollection> truthToken_;
};

DEFINE_FWK_MODULE(L1SmartPixelsClusterTPLinkProducer);
