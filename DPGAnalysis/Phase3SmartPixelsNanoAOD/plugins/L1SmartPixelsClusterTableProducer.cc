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
// PER-CLUSTER PAYLOAD. Measured 88.6 stored bits/cluster (11.1 B) on the original
// 12-column table; the table has since grown the global-frame geometry block and
// raised localX/localY from 10 to 16 mantissa bits, so it is now larger and the
// figure needs re-measuring before being quoted. What has NOT changed is why the
// precision is what it is: globalPhi and the direction columns are stored at 16
// mantissa bits (~4.8e-5 rad) because reconstructing phi from lower-precision
// Cartesian columns would give ~1 mrad, several times the sensor resolution, and
// would make any angular result an artefact of storage rather than of physics.
//
// NAMING. Frame prefix only for sensor quantities (localX/localY, localCotAlpha/
// Beta, globalR/Phi/Z, globalClusterPhi/CotTheta); "reco" is redundant because a
// cluster carries ONLY reco quantities and truth is reachable solely through the
// tp link. Truth carries the tp prefix, matching the convention already used on
// L1TTrack (tpPt, tpEta, tpPdgId). cotAlpha/cotBeta are module-frame BY
// DEFINITION (PixelAV) and so are named local*, not given a global twin.
//
// tpPt is the pT of the PARENT of the cluster's DOMINANT charge contributor,
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
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "DataFormats/Phase3SmartPixels/interface/SmartPixelsFrames.h"
#include "DataFormats/Phase3SmartPixels/interface/SmartPixelsRecHit.h"
#include "DataFormats/Phase3SmartPixels/interface/SmartPixelsRecHitTruth.h"
#include "Geometry/CommonTopologies/interface/PixelGeomDetUnit.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonTopologies/interface/PixelGeomDetUnit.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
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
        geomToken_(esConsumes()),
        tableName_(cfg.getParameter<std::string>("tableName")),
        maxLayer_(cfg.getParameter<unsigned>("maxLayer")),
        doTruth_(cfg.getParameter<bool>("doTruth")) {
    produces<nanoaod::FlatTable>();
  }

  void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override {
    const auto& topo = iSetup.getData(topoToken_);
    const auto& geom = iSetup.getData(geomToken_);
    const auto& recHits = iEvent.get(recHitToken_);
    const auto& truthColl = iEvent.get(truthToken_);
    std::vector<uint8_t> layer, sizeX, sizeY;
    std::vector<uint16_t> size;
    std::vector<uint32_t> detId;
    std::vector<float> localX, localY, sigX, sigY, charge;
    std::vector<float> globalR, globalPhi, globalZ;
    std::vector<float> globalClusterPhi, globalClusterCotTheta, tpGlobalClusterPhi, tpGlobalClusterCotTheta;
    std::vector<float> sigGlobalClusterPhi, sigGlobalClusterCotTheta;
    unsigned closureFail = 0;
    std::vector<float> localCotAlpha, localCotBeta, sigAlpha, sigBeta;
    std::vector<uint8_t> hasAlpha, hasBeta;
    std::vector<float> tpPt, tpChargeFrac, tpLocalCotAlpha, tpLocalCotBeta;
    std::vector<int32_t> tpIdx;

    for (const auto& dsv : recHits) {
      const DetId did(dsv.detId());
      const unsigned lay = topo.pxbLayer(did);
      if (lay < 1 || lay > maxLayer_)
        continue;
      const auto* pdu = dynamic_cast<const PixelGeomDetUnit*>(geom.idToDet(did));
      const auto tsv = truthColl.find(dsv.detId());
      const bool haveTruth = doTruth_ && tsv != truthColl.end() && tsv->size() == dsv.size();
      if (doTruth_ && !haveTruth)
        throw cms::Exception("SmartPixelsRecHitTruthMisaligned")
            << "cluster table: rechit/truth disagree on det " << dsv.detId();

      for (size_t j = 0; j < dsv.size(); ++j) {
        const auto& rh = dsv[j];
        layer.push_back(static_cast<uint8_t>(lay));
        detId.push_back(did.rawId());
        // GLOBAL position and direction, computed HERE because this is where the
        // geometry is. nano carries none, so without these every analysis would
        // re-derive them from detId with its own copy of the tracker geometry --
        // duplicated work and a silent divergence risk on TILTED modules (TBPX
        // tilt reaches 16.5 deg, so the rotation is per-module, not per-layer).
        //
        // CYLINDRICAL, not Cartesian, and deliberately: phi stored directly at 16
        // mantissa bits gives ~1e-4 rad, whereas phi RECONSTRUCTED from two
        // Cartesian columns inherits their relative precision -- at the 10 bits
        // used for localX/localY that is ~1 mrad, several times the sensor
        // resolution, which would make any angular study an artefact of storage.
        if (pdu != nullptr) {
          const auto gp = pdu->toGlobal(rh.localPosition());
          globalR.push_back(std::hypot(gp.x(), gp.y()));
          globalPhi.push_back(std::atan2(gp.y(), gp.x()));
          globalZ.push_back(gp.z());
          const auto gr = smartpixels::toGlobalDirection(*pdu, rh.cotAlpha(), rh.cotBeta(),
                                                         rh.sigAlpha(), rh.sigBeta());
          globalClusterPhi.push_back(gr.valid ? gr.dirPhi : -999.f);
          globalClusterCotTheta.push_back(gr.valid ? gr.dirCotTheta : -999.f);
          sigGlobalClusterPhi.push_back(gr.valid ? gr.sigDirPhi : -999.f);
          sigGlobalClusterCotTheta.push_back(gr.valid ? gr.sigDirCotTheta : -999.f);
          // CLOSURE: rotate the global direction back and require the module-frame
          // angles to reappear. A mis-applied rotation on a tilted module would
          // otherwise be a large, silent error.
          if (rh.hasAlpha() && rh.hasBeta() &&
              !smartpixels::closesBackToModule(*pdu, gr, rh.cotAlpha(), rh.cotBeta()))
            ++closureFail;
        } else {
          globalR.push_back(-999.f); globalPhi.push_back(-999.f); globalZ.push_back(-999.f);
          globalClusterPhi.push_back(-999.f); globalClusterCotTheta.push_back(-999.f);
          sigGlobalClusterPhi.push_back(-999.f); sigGlobalClusterCotTheta.push_back(-999.f);
        }
        localX.push_back(rh.localPosition().x());
        localY.push_back(rh.localPosition().y());
        sigX.push_back(std::sqrt(std::max(0.f, static_cast<float>(rh.localPositionError().xx()))));
        sigY.push_back(std::sqrt(std::max(0.f, static_cast<float>(rh.localPositionError().yy()))));
        sizeX.push_back(static_cast<uint8_t>(std::min<unsigned>(rh.sizeX(), 255)));
        sizeY.push_back(static_cast<uint8_t>(std::min<unsigned>(rh.sizeY(), 255)));
        size.push_back(rh.size());
        charge.push_back(rh.charge());
        localCotAlpha.push_back(rh.cotAlpha());
        localCotBeta.push_back(rh.cotBeta());
        sigAlpha.push_back(rh.sigAlpha());
        sigBeta.push_back(rh.sigBeta());
        hasAlpha.push_back(rh.hasAlpha() ? 1 : 0);
        hasBeta.push_back(rh.hasBeta() ? 1 : 0);

        if (haveTruth) {
          const auto& tr = (*tsv)[j];
          tpIdx.push_back(tr.hasTp() ? static_cast<int32_t>(tr.dominantTp().key()) : -1);
          tpPt.push_back(tr.hasTp() ? static_cast<float>(tr.dominantTp()->pt()) : -999.f);
          tpChargeFrac.push_back(tr.chargeFrac());
          tpLocalCotAlpha.push_back(tr.trueCotAlpha());
          tpLocalCotBeta.push_back(tr.trueCotBeta());
          if (pdu != nullptr && tr.trueCotAlpha() > -900.f) {
            const auto gt = smartpixels::toGlobalDirection(*pdu, tr.trueCotAlpha(), tr.trueCotBeta());
            tpGlobalClusterPhi.push_back(gt.valid ? gt.dirPhi : -999.f);
            tpGlobalClusterCotTheta.push_back(gt.valid ? gt.dirCotTheta : -999.f);
          } else {
            tpGlobalClusterPhi.push_back(-999.f); tpGlobalClusterCotTheta.push_back(-999.f);
          }
        }
      }
    }

    auto tab = std::make_unique<nanoaod::FlatTable>(layer.size(), tableName_, false, false);
    tab->addColumn<uint8_t>("layer", layer, "TBPX layer 1..4");
    tab->addColumn<uint32_t>("detId", detId, "module rawId; join key to the refit hit table detId");
    tab->addColumn<float>("localX", localX, "cluster position, module-local x [cm]", 16);
    tab->addColumn<float>("localY", localY, "cluster position, module-local y [cm]", 16);
    tab->addColumn<float>("globalR", globalR,
                          "POSITION. Cylindrical radius of WHERE THE CLUSTER IS [cm]", 16);
    tab->addColumn<float>("globalPhi", globalPhi,
                          "POSITION. Azimuth of WHERE THE CLUSTER IS, in CMS global coordinates "
                          "[rad]. This is a location on the detector surface and has NOTHING to do "
                          "with any direction or angle estimate: it is atan2(y,x) of the cluster "
                          "centroid. Pairs with globalR and globalZ to give the full position. "
                          "Do NOT confuse with globalClusterPhi, which is a DIRECTION. Stored "
                          "directly at 16 mantissa bits (~5e-5 rad) because reconstructing it from "
                          "Cartesian columns would inherit their precision (~1 mrad)", 16);
    tab->addColumn<float>("globalZ", globalZ,
                          "POSITION. CMS global z of WHERE THE CLUSTER IS [cm]", 16);
    tab->addColumn<float>("globalClusterPhi", globalClusterPhi,
                          "DIRECTION. Azimuth of the SMART-PIXEL ML ANGLE ESTIMATE for this "
                          "cluster, rotated into CMS global coordinates [rad]. This is the "
                          "estimated direction of the PARTICLE that made the cluster -- the ML "
                          "regressor output, carrying the PixelAV response -- NOT a position. It "
                          "is the global counterpart of localCotAlpha/localCotBeta, which are "
                          "module-frame by definition (PixelAV) and have no meaningful global "
                          "variant. Do NOT confuse with globalPhi, which is WHERE the cluster is. "
                          "Uses the per-module surface rotation (TBPX tilt reaches 16.5 deg, so it "
                          "is not a per-layer constant). Uncertainty: sigGlobalClusterPhi", 16);
    tab->addColumn<float>("globalClusterCotTheta", globalClusterCotTheta,
                          "DIRECTION. cot(theta) = pz/pt of the SMART-PIXEL ML ANGLE ESTIMATE, CMS "
                          "global frame. Chosen over eta because the r-z Hough wants "
                          "z = z0 + r*cotTheta directly. Uncertainty: sigGlobalClusterCotTheta", 16);
    tab->addColumn<float>("sigGlobalClusterPhi", sigGlobalClusterPhi,
                          "Uncertainty on globalClusterPhi [rad], propagated from the module-frame "
                          "sigAlpha/sigBeta through the SAME per-module rotation by a numerical "
                          "Jacobian, added in quadrature. Alpha and beta are treated as "
                          "independent, matching how the refit applies them as two independent "
                          "scalar Kalman updates", 12);
    tab->addColumn<float>("sigGlobalClusterCotTheta", sigGlobalClusterCotTheta,
                          "Uncertainty on globalClusterCotTheta, same propagation", 12);
    tab->addColumn<float>("sigX", sigX, "CPE position uncertainty, local x [cm]", 10);
    tab->addColumn<float>("sigY", sigY, "CPE position uncertainty, local y [cm]", 10);
    tab->addColumn<uint8_t>("sizeX", sizeX, "cluster bounding-box extent in pixels, local x");
    tab->addColumn<uint8_t>("sizeY", sizeY, "cluster bounding-box extent in pixels, local y");
    tab->addColumn<uint16_t>("size", size, "FIRED-PIXEL COUNT (not the bounding box): charge/size is "
                                           "the real charge density");
    tab->addColumn<float>("charge", charge, "cluster charge [ADC]", 10);
    tab->addColumn<float>("localCotAlpha", localCotAlpha,
                          "SENSOR angle estimate, module-local cotAlpha (-999 if none)", 12);
    tab->addColumn<float>("localCotBeta", localCotBeta, "SENSOR angle estimate, cotBeta", 12);
    tab->addColumn<float>("sigAlpha", sigAlpha, "angle-estimator resolution, alpha", 10);
    tab->addColumn<float>("sigBeta", sigBeta, "angle-estimator resolution, beta", 10);
    tab->addColumn<uint8_t>("hasAlpha", hasAlpha,
                            "sensor reports an alpha at all (payload validity gate / grazing clamp)");
    tab->addColumn<uint8_t>("hasBeta", hasBeta, "sensor reports a beta at all");
    if (doTruth_) {
      tab->addColumn<int32_t>("tpIdx", tpIdx,
                              "TRUTH-ONLY: TrackingParticle index of the dominant charge contributor, "
                              "or -1. Join key against spixMatchedTpIdx");
      tab->addColumn<float>("tpPt", tpPt, "TRUTH-ONLY: pT [GeV] of the dominant TP", 10);
      tab->addColumn<float>("tpLocalCotAlpha", tpLocalCotAlpha,
                            "TRUTH-ONLY: TRUE incidence cotAlpha at this module (helix-propagated), "
                            "i.e. what the sensor is trying to measure", 12);
      tab->addColumn<float>("tpLocalCotBeta", tpLocalCotBeta, "TRUTH-ONLY: true incidence cotBeta", 12);
      tab->addColumn<float>("tpGlobalClusterPhi", tpGlobalClusterPhi,
                            "TRUTH-ONLY: true direction, global phi [rad]", 16);
      tab->addColumn<float>("tpGlobalClusterCotTheta", tpGlobalClusterCotTheta,
                            "TRUTH-ONLY: true direction, global cot(theta)", 16);
      tab->addColumn<float>("tpChargeFrac", tpChargeFrac,
                            "TRUTH-ONLY: dominant contributor share of the cluster charge", 10);
    }
    if (closureFail)
        throw cms::Exception("SmartPixelsFrameClosureFailed")
            << closureFail << " clusters failed the module<->global direction closure test "
            << "(rotate the stored global direction back, require cotAlpha/cotBeta to reappear). "
            << "The per-module rotation is being mis-applied; on a 16-degree-tilted TBPX module "
            << "that is a large silent error, so this refuses to emit the table.";
    iEvent.put(std::move(tab));
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("smartPixelsRecHits", edm::InputTag("spixSmartPixelsRecHits"))
        ->setComment("SmartPixelsRecHit + Truth (same label): the single source of the "
                     "sensor angle, shared with the refit so selClusterIdx is exact");
    desc.add<std::string>("tableName", "L1TSmartPixelsCluster");
    desc.add<unsigned>("maxLayer", 4)->setComment("highest TBPX layer kept (SmartPixels scope is 1..4)");
    desc.add<bool>("doTruth", true)
        ->setComment("attach the TRUTH-ONLY tp* block (tpIdx, tpPt, tpChargeFrac, "
                     "tpLocalCot*, tpGlobalCluster*)");
    descriptions.addWithDefaultLabel(desc);
  }

private:
  const edm::EDGetTokenT<SmartPixelsRecHitCollection> recHitToken_;
  const edm::EDGetTokenT<SmartPixelsRecHitTruthCollection> truthToken_;
  const edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> topoToken_;
  const edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> geomToken_;
  const std::string tableName_;
  const unsigned maxLayer_;
  const bool doTruth_;
};

DEFINE_FWK_MODULE(L1SmartPixelsClusterTableProducer);
