// Nano table of EVERY IT pixel cluster, deliberately UNTRUNCATED.
//
// WHY UNTRUNCATED. The refit's match stage faces ~20-40 clusters on the crossed
// module at PU200 (measured: 40.6 / 32.7 / 31.6 / 19.7 per occupied module on
// L1-L4, p95 up to 72). The static window admits ~2 of those and
// maxHitsPerWindow truncates at 8. Every combinatorics question -- how wide a
// window can be afforded, what a covariance-derived window would admit, how much
// sequential refitting reduces the candidate pool -- is a question about the
// clusters the current window THREW AWAY. A truncated or window-filtered table
// therefore cannot answer any of them, and a truth-pT-filtered table is worse
// than useless: it removes precisely the soft clusters that do the confusing,
// which would make the measured combinatorics optimistic by ~45x.
//
// It is consequently a large table (~26.5k rows/event at PU200) and belongs only
// to the L1PFTrkNanoSmartPixClusters tier. It is never added to a physics tier.
// Measured cost is 0.29 MB/event -- see the payload breakdown below.
//
// WHAT A "MODULE" IS (D121 geometry, measured by SmartPixelsClusterCensusAnalyzer).
// One module == one PixelGeomDetUnit == one DetId == the unit this table's detId
// column identifies and the unit the refit candidate loop scans:
//
//   layer  modules  rows x cols  pixels/module  ROCs(x,y)  rows,cols/ROC  pitch [um]
//     L1     216     672 x 216      145 152        1 x 1     672 x 216     25 x 100
//     L2     216     672 x 434      291 648        1 x 2     672 x 217     25 x 100
//     L3     180    1354 x 434      587 636        2 x 2     677 x 217     25 x 100
//     L4     252    1354 x 434      587 636        2 x 2     677 x 217     25 x 100
//                                                        (864 TBPX modules total)
//
// There is NO separate "sensor" DetId in CMSSW: the sensor and its readout-chip
// array are a single detUnit addressed as one pixel matrix. The nearest thing to a
// sensor subdivision is the topology's ROC tiling above (rocsX x rocsY) -- an L1
// module is one tile, an L3/L4 module is a 2x2 array of tiles ~677x217 pixels each.
// Treat that tiling as CMSSW's pixel-addressing granularity, NOT as a verified 1:1
// map to physical RD53 chips; nothing here establishes the latter.
//
// Physical extent follows from pitch: an L1 module is 672*25um x 216*100um =
// 16.8 x 21.6 mm; an L3/L4 module is 33.9 x 43.4 mm.
//
// PER-CLUSTER PAYLOAD, measured on a real file (not estimated):
//
//   column             stored bits/cluster
//   localX                  17.8
//   localY                  17.7
//   charge                  12.6
//   truthPt                 12.2   TRUTH-ONLY
//   detId                    7.5
//   sizeX                    4.4
//   sizeY                    4.1
//   sigY                     3.9
//   sigX                     3.8
//   truthChargeFrac          2.2   TRUTH-ONLY
//   layer                    1.2
//   truthLinked              1.1   TRUTH-ONLY
//   TOTAL                   88.6 bits = 11.1 B/cluster stored (37.0 B raw)
//
// Floats are written with 10-bit mantissa precision, which is why sigX/sigY (nearly
// constant per module) cost under 4 bits while localX/localY (genuinely uniform
// across the module) cost ~18. At the measured PU200 occupancy of 26 479
// clusters/event that is 0.29 MB/event stored -- 3.4x smaller than a naive
// 37 B/cluster estimate, which is why it was measured rather than assumed.
//
// truthPt is the pT of the PARENT of the cluster's DOMINANT charge contributor,
// assigned by the SAME charge-share logic L1SmartPixelsTrackProducer uses
// (per-(eventId, SimTrackId) ADC sum over the cluster's pixels via
// PixelDigi::pixelToChannel). TRUTH-ONLY: it exists to bound what an ideal
// pT-discriminating on-sensor readout could buy, and must never be used as a
// selection in anything claiming to be deployable.

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "DataFormats/Phase3SmartPixels/interface/SmartPixelsRecHit.h"
#include "DataFormats/Phase3SmartPixels/interface/SmartPixelsRecHitTruth.h"
#include "Geometry/CommonTopologies/interface/PixelGeomDetUnit.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "L1Trigger/Phase3SmartPixels/interface/SmartPixelsParentMap.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/TrackerDigiSimLink/interface/PixelDigiSimLink.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"

#include <cmath>
#include <cstdint>
#include <map>
#include <string>
#include <vector>

class L1SmartPixelsClusterTableProducer : public edm::stream::EDProducer<> {
public:
  explicit L1SmartPixelsClusterTableProducer(const edm::ParameterSet& cfg)
      : recHitToken_(consumes<SmartPixelsRecHitCollection>(cfg.getParameter<edm::InputTag>("smartPixelsRecHits"))),
        truthToken_(consumes<SmartPixelsRecHitTruthCollection>(cfg.getParameter<edm::InputTag>("smartPixelsRecHits"))),
        topoToken_(esConsumes()),
        tableName_(cfg.getParameter<std::string>("tableName")),
        maxLayer_(cfg.getParameter<unsigned>("maxLayer")),
        doTruth_(cfg.getParameter<bool>("doTruth")) {
    produces<nanoaod::FlatTable>();
  }

  void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override {
    const auto& topo = iSetup.getData(topoToken_);
    const auto& recHits = iEvent.get(recHitToken_);
    const auto& truthColl = iEvent.get(truthToken_);

    std::vector<uint8_t> layer, sizeX, sizeY, truthLinked;
    std::vector<uint16_t> size;
    std::vector<uint32_t> detId;
    std::vector<float> localX, localY, sigX, sigY, charge;
    std::vector<float> recoCotAlpha, recoCotBeta, sigAlpha, sigBeta;
    std::vector<uint8_t> hasAlpha, hasBeta;
    std::vector<float> truthPt, truthChargeFrac, truthCotAlpha, truthCotBeta;
    std::vector<int32_t> truthTpIdx;

    for (const auto& dsv : recHits) {
      const DetId did(dsv.detId());
      const unsigned lay = topo.pxbLayer(did);
      if (lay < 1 || lay > maxLayer_)
        continue;
      const auto tsv = truthColl.find(dsv.detId());
      const bool haveTruth = doTruth_ && tsv != truthColl.end() && tsv->size() == dsv.size();
      if (doTruth_ && !haveTruth)
        throw cms::Exception("SmartPixelsRecHitTruthMisaligned")
            << "cluster table: rechit/truth disagree on det " << dsv.detId();

      for (size_t j = 0; j < dsv.size(); ++j) {
        const auto& rh = dsv[j];
        layer.push_back(static_cast<uint8_t>(lay));
        detId.push_back(did.rawId());
        localX.push_back(rh.localPosition().x());
        localY.push_back(rh.localPosition().y());
        sigX.push_back(std::sqrt(std::max(0.f, static_cast<float>(rh.localPositionError().xx()))));
        sigY.push_back(std::sqrt(std::max(0.f, static_cast<float>(rh.localPositionError().yy()))));
        sizeX.push_back(static_cast<uint8_t>(std::min<unsigned>(rh.sizeX(), 255)));
        sizeY.push_back(static_cast<uint8_t>(std::min<unsigned>(rh.sizeY(), 255)));
        size.push_back(rh.size());
        charge.push_back(rh.charge());
        recoCotAlpha.push_back(rh.cotAlpha());
        recoCotBeta.push_back(rh.cotBeta());
        sigAlpha.push_back(rh.sigAlpha());
        sigBeta.push_back(rh.sigBeta());
        hasAlpha.push_back(rh.hasAlpha() ? 1 : 0);
        hasBeta.push_back(rh.hasBeta() ? 1 : 0);

        if (haveTruth) {
          const auto& tr = (*tsv)[j];
          truthLinked.push_back(tr.hasTp() ? 1 : 0);
          truthTpIdx.push_back(tr.hasTp() ? static_cast<int32_t>(tr.dominantTp().key()) : -1);
          truthPt.push_back(tr.hasTp() ? static_cast<float>(tr.dominantTp()->pt()) : -999.f);
          truthChargeFrac.push_back(tr.chargeFrac());
          truthCotAlpha.push_back(tr.trueCotAlpha());
          truthCotBeta.push_back(tr.trueCotBeta());
        }
      }
    }

    auto tab = std::make_unique<nanoaod::FlatTable>(layer.size(), tableName_, false, false);
    tab->addColumn<uint8_t>("layer", layer, "TBPX layer 1..4");
    tab->addColumn<uint32_t>("detId", detId, "module rawId; join key to the refit hit table detId");
    tab->addColumn<float>("localX", localX, "cluster position, module-local x [cm]", 10);
    tab->addColumn<float>("localY", localY, "cluster position, module-local y [cm]", 10);
    tab->addColumn<float>("sigX", sigX, "CPE position uncertainty, local x [cm]", 10);
    tab->addColumn<float>("sigY", sigY, "CPE position uncertainty, local y [cm]", 10);
    tab->addColumn<uint8_t>("sizeX", sizeX, "cluster bounding-box extent in pixels, local x");
    tab->addColumn<uint8_t>("sizeY", sizeY, "cluster bounding-box extent in pixels, local y");
    tab->addColumn<uint16_t>("size", size, "FIRED-PIXEL COUNT (not the bounding box): charge/size is "
                                           "the real charge density");
    tab->addColumn<float>("charge", charge, "cluster charge [ADC]", 10);
    tab->addColumn<float>("recoCotAlpha", recoCotAlpha,
                          "SENSOR angle estimate, module-local cotAlpha (-999 if none)", 12);
    tab->addColumn<float>("recoCotBeta", recoCotBeta, "SENSOR angle estimate, cotBeta", 12);
    tab->addColumn<float>("sigAlpha", sigAlpha, "angle-estimator resolution, alpha", 10);
    tab->addColumn<float>("sigBeta", sigBeta, "angle-estimator resolution, beta", 10);
    tab->addColumn<uint8_t>("hasAlpha", hasAlpha,
                            "sensor reports an alpha at all (payload validity gate / grazing clamp)");
    tab->addColumn<uint8_t>("hasBeta", hasBeta, "sensor reports a beta at all");
    if (doTruth_) {
      tab->addColumn<uint8_t>("truthLinked", truthLinked, "TRUTH-ONLY: cluster has a dominant TP");
      tab->addColumn<int32_t>("truthTpIdx", truthTpIdx,
                              "TRUTH-ONLY: TrackingParticle index of the dominant charge contributor, "
                              "or -1. Join key against spixMatchedTpIdx");
      tab->addColumn<float>("truthPt", truthPt, "TRUTH-ONLY: pT [GeV] of the dominant TP", 10);
      tab->addColumn<float>("truthCotAlpha", truthCotAlpha,
                            "TRUTH-ONLY: TRUE incidence cotAlpha at this module (helix-propagated), "
                            "i.e. what the sensor is trying to measure", 12);
      tab->addColumn<float>("truthCotBeta", truthCotBeta, "TRUTH-ONLY: true incidence cotBeta", 12);
      tab->addColumn<float>("truthChargeFrac", truthChargeFrac,
                            "TRUTH-ONLY: dominant contributor share of the cluster charge", 10);
    }
    iEvent.put(std::move(tab));
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("smartPixelsRecHits", edm::InputTag("spixSmartPixelsRecHits"))
        ->setComment("SmartPixelsRecHit + Truth (same label): the single source of the "
                     "sensor angle, shared with the refit so selClusterIdx is exact");
    desc.add<std::string>("tableName", "L1TSmartPixelsCluster");
    desc.add<unsigned>("maxLayer", 4)->setComment("highest TBPX layer kept (SmartPixels scope is 1..4)");
    desc.add<bool>("doTruth", true)->setComment("attach TRUTH-ONLY truthPt/truthChargeFrac/truthLinked");
    descriptions.addWithDefaultLabel(desc);
  }

private:
  const edm::EDGetTokenT<SmartPixelsRecHitCollection> recHitToken_;
  const edm::EDGetTokenT<SmartPixelsRecHitTruthCollection> truthToken_;
  const edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> topoToken_;
  const std::string tableName_;
  const unsigned maxLayer_;
  const bool doTruth_;
};

DEFINE_FWK_MODULE(L1SmartPixelsClusterTableProducer);
