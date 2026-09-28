// SmartPixels track -> IT cluster LINK table: one row per (track, cluster) pair,
//   trackIdx   -> row of the track in trackTableName
//   clusterIdx -> row of the cluster in clusterTableName (L1TSmartPixelsCluster)
// The clusters of a track are the IT hits its variant used: for digiRefit the hits
// accepted into the Kalman refit, read from the refit sidecar. A variant that uses no IT
// hits (passthrough) has no sidecar: leave `sidecar` empty and the table has no rows, so
// the schema does not depend on the variant.

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

#include "L1Trigger/Phase3SmartPixels/interface/SmartPixelsRefitSidecar.h"

#include <cstdint>
#include <string>
#include <vector>

class L1SmartPixelsTrackClusterLinkTableProducer : public edm::global::EDProducer<> {
public:
  using TTTrackCollection = std::vector<TTTrack<Ref_Phase2TrackerDigi_>>;

  explicit L1SmartPixelsTrackClusterLinkTableProducer(const edm::ParameterSet& cfg)
      : tracksToken_(consumes<TTTrackCollection>(cfg.getParameter<edm::InputTag>("tracks"))),
        linkTableName_(cfg.getParameter<std::string>("linkTableName")),
        trackTableName_(cfg.getParameter<std::string>("trackTableName")),
        clusterTableName_(cfg.getParameter<std::string>("clusterTableName")) {
    const auto& sidecar = cfg.getParameter<edm::InputTag>("sidecar");
    if (!sidecar.label().empty())
      sidecarToken_ = consumes<smartpixels::SmartPixelsRefitSidecar>(sidecar);
    produces<nanoaod::FlatTable>();
  }

  void produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup&) const override {
    const auto& tracks = iEvent.get(tracksToken_);
    std::vector<int32_t> trackIdx, clusterIdx;
    if (!sidecarToken_.isUninitialized()) {
      const auto& sidecar = iEvent.get(sidecarToken_);
      if (sidecar.hitInfo.size() != tracks.size())
        throw cms::Exception("SmartPixelsSyncBroken")
            << linkTableName_ << ": sidecar has " << sidecar.hitInfo.size() << " track records for "
            << tracks.size() << " tracks";
      for (size_t it = 0; it < tracks.size(); ++it) {
        for (const auto& hit : sidecar.hitInfo[it]) {
          if (!(hit.flags & smartpixels::hitflag::kHitAccepted))
            continue;
          if (hit.selClusterIdx < 0)
            throw cms::Exception("SmartPixelsSyncBroken")
                << linkTableName_ << ": track " << it << " has an accepted hit without a cluster row";
          trackIdx.push_back(static_cast<int32_t>(it));
          clusterIdx.push_back(hit.selClusterIdx);
        }
      }
    }

    auto table = std::make_unique<nanoaod::FlatTable>(trackIdx.size(), linkTableName_, false, false);
    table->addColumn<int32_t>("trackIdx", trackIdx, "row of the track in " + trackTableName_);
    table->addColumn<int32_t>("clusterIdx", clusterIdx, "row of the cluster in " + clusterTableName_);
    table->setDoc("IT clusters used by the " + trackTableName_ +
                  " tracks (digiRefit: the hits accepted into the refit; one row per track-cluster pair)");
    iEvent.put(std::move(table));
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("tracks", edm::InputTag("l1tSmartPixelsTrackProducer", "Level1TTTracks"));
    desc.add<edm::InputTag>("sidecar", edm::InputTag())
        ->setComment("refit sidecar of the tracks; empty for a variant that uses no IT hits (no rows)");
    desc.add<std::string>("linkTableName", "L1TSmartPixelsTrackCluster");
    desc.add<std::string>("trackTableName", "L1TSmartPixelsTrack");
    desc.add<std::string>("clusterTableName", "L1TSmartPixelsCluster");
    descriptions.addWithDefaultLabel(desc);
  }

private:
  const edm::EDGetTokenT<TTTrackCollection> tracksToken_;
  edm::EDGetTokenT<smartpixels::SmartPixelsRefitSidecar> sidecarToken_;
  const std::string linkTableName_;
  const std::string trackTableName_;
  const std::string clusterTableName_;
};

DEFINE_FWK_MODULE(L1SmartPixelsTrackClusterLinkTableProducer);
