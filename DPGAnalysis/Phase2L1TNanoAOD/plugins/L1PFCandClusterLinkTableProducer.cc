// Emits an extension column on a L1 PF candidate table with the index of
// the calo cluster the candidate's caloPtr points to, guarded by ProductID
// match against the configured cluster collection (candidates whose caloPtr
// targets a different collection, or have none, get -1). The endcap
// correlator layer-1 consumes the HGCal multiclusters directly, so for
// endcap candidates this links to the HGCal multicluster table.
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/L1TParticleFlow/interface/PFCandidate.h"
#include "DataFormats/L1Trigger/interface/L1Candidate.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"

class L1PFCandClusterLinkTableProducer : public edm::stream::EDProducer<> {
public:
  explicit L1PFCandClusterLinkTableProducer(const edm::ParameterSet& cfg)
      : candsToken_(consumes<l1t::PFCandidateCollection>(cfg.getParameter<edm::InputTag>("cands"))),
        clustersToken_(consumes<edm::View<l1t::L1Candidate>>(cfg.getParameter<edm::InputTag>("clusters"))),
        candTableName_(cfg.getParameter<std::string>("candTableName")),
        columnName_(cfg.getParameter<std::string>("columnName")),
        clusterTableName_(cfg.getParameter<std::string>("clusterTableName")) {
    produces<nanoaod::FlatTable>();
  }

  void produce(edm::Event& iEvent, const edm::EventSetup&) override {
    auto candsHandle = iEvent.getHandle(candsToken_);
    auto clustersHandle = iEvent.getHandle(clustersToken_);

    const size_t nCands = candsHandle.isValid() ? candsHandle->size() : 0;
    std::vector<int32_t> clusterIdx(nCands, -1);
    if (candsHandle.isValid() && clustersHandle.isValid()) {
      for (size_t icand = 0; icand < nCands; ++icand) {
        const auto& caloPtr = (*candsHandle)[icand].caloPtr();
        if (caloPtr.isNonnull() && caloPtr.id() == clustersHandle.id()) {
          clusterIdx[icand] = caloPtr.key();
        }
      }
    }

    auto candTable = std::make_unique<nanoaod::FlatTable>(nCands, candTableName_, false, true);
    candTable->addColumn<int32_t>(
        columnName_, clusterIdx, "index of the calo cluster in the " + clusterTableName_ + " table (-1 if none)");
    iEvent.put(std::move(candTable));
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("cands", edm::InputTag("l1tLayer2Deregionizer", "Puppi"));
    desc.add<edm::InputTag>("clusters",
                            edm::InputTag("l1tHGCalBackEndLayer2Producer", "HGCalBackendLayer2Processor3DClustering"));
    desc.add<std::string>("candTableName", "L1PuppiCand");
    desc.add<std::string>("columnName", "hgcClusterIdx");
    desc.add<std::string>("clusterTableName", "L1HGCCluster");
    descriptions.addWithDefaultLabel(desc);
  }

private:
  const edm::EDGetTokenT<l1t::PFCandidateCollection> candsToken_;
  const edm::EDGetTokenT<edm::View<l1t::L1Candidate>> clustersToken_;
  const std::string candTableName_;
  const std::string columnName_;
  const std::string clusterTableName_;
};

DEFINE_FWK_MODULE(L1PFCandClusterLinkTableProducer);
