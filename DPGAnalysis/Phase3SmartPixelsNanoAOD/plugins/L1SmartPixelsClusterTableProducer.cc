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
// It is consequently a large table (~74k rows/event at PU200 with ID1-ID4, measured
// 3.1 MB/event) and belongs only to the @SPClustersStubs nano flavour. It is never
// added to a physics tier. See the payload breakdown below.
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
// ONE TABLE FOR BARREL AND DISCS. Rows are every TBPX cluster (layers 1..maxLayer)
// FIRST, then every disc cluster in discScope (default sp = ID1-ID4, TFPX disks 1-4),
// in the order of smartpixels::clusterTableOrder, which the refit's selClusterIdx
// uses too. The barrel prefix is the historical barrel-only table row for row.
// surface (1-4 = IL1-IL4, 10 + disk on the discs) is the one column to select on;
// (detId, idxInDet) is a stable per-event key. The table emits EVERY rec hit and
// refuses to run if the rec hits were made with a different disc scope or layer
// range, since the refit's row indices assume table == rec hits.
//
// DISC COLUMN SEMANTICS. On a disc module local x is along r-phi and localY runs
// along the radius, pointing INWARD, so sigY is sigma_r (not sigma_z as on TBPX).
// localCotAlpha/localCotBeta stay SIGNED in the module frame; the sign of the disc
// cotBeta flips with the module facing (the two blade faces point opposite ways
// along z), so compare directions through globalClusterCotTheta, not raw cotBeta.
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

#include "FWCore/Common/interface/Provenance.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DPGAnalysis/Phase2L1TNanoAOD/interface/L1TPRowLink.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "DataFormats/Phase3SmartPixels/interface/SmartPixelsFrames.h"
#include "DataFormats/Phase3SmartPixels/interface/SmartPixelsRecHit.h"
#include "DataFormats/Phase3SmartPixels/interface/SmartPixelsRecHitTruth.h"
#include "Geometry/CommonTopologies/interface/PixelGeomDetUnit.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "L1Trigger/Phase3SmartPixels/interface/SmartPixelsParentMap.h"
#include "L1Trigger/Phase3SmartPixels/interface/SmartPixelsSurface.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/TrackerDigiSimLink/interface/PixelDigiSimLink.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"

#include <cmath>
#include <cstdint>
#include <map>
#include <string>
#include <vector>

namespace {
  constexpr double kClosureMaxAbsCot = 12.;  // the refit's measAngleMaxAbs
}

class L1SmartPixelsClusterTableProducer : public edm::stream::EDProducer<> {
public:
  explicit L1SmartPixelsClusterTableProducer(const edm::ParameterSet& cfg)
      : recHitToken_(consumes<SmartPixelsRecHitCollection>(cfg.getParameter<edm::InputTag>("smartPixelsRecHits"))),
        truthToken_(consumes<SmartPixelsRecHitTruthCollection>(cfg.getParameter<edm::InputTag>("smartPixelsRecHits"))),
        topoToken_(esConsumes()),
        geomToken_(esConsumes()),
        tableName_(cfg.getParameter<std::string>("tableName")),
        discScope_(cfg.getParameter<std::string>("discScope")),
        maxLayer_(cfg.getParameter<unsigned>("maxLayer")),
        maxDisk_(smartpixels::discScopeMaxDisk(discScope_)),
        doTruth_(cfg.getParameter<bool>("doTruth")),
        closurePrecision_(cfg.getParameter<bool>("closurePrecision")) {
    if (doTruth_)
      tpRowMapToken_ = consumes<edm::ValueMap<int>>(cfg.getParameter<edm::InputTag>("tpRowMap"));
    produces<nanoaod::FlatTable>();
  }

  void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override {
    const auto& topo = iSetup.getData(topoToken_);
    const auto& geom = iSetup.getData(geomToken_);
    const auto recHitHandle = iEvent.getHandle(recHitToken_);
    const auto& recHits = *recHitHandle;
    const auto& truthColl = iEvent.get(truthToken_);
    edm::Handle<edm::ValueMap<int>> tpRowMapHandle;
    if (doTruth_)
      iEvent.getByToken(tpRowMapToken_, tpRowMapHandle);
    const edm::ValueMap<int>* tpRowMap = tpRowMapHandle.isValid() ? tpRowMapHandle.product() : nullptr;
    if (!scopeChecked_) {
      const auto& rhCfg = edm::parameterSet(recHitHandle.provenance()->stable(), iEvent.processHistory());
      if (rhCfg.getParameter<std::string>("discScope") != discScope_ ||
          rhCfg.getParameter<unsigned>("maxLayer") != maxLayer_)
        throw cms::Exception("Configuration")
            << "cluster table discScope/maxLayer (" << discScope_ << ", " << maxLayer_
            << ") differ from the rec hits' (" << rhCfg.getParameter<std::string>("discScope") << ", "
            << rhCfg.getParameter<unsigned>("maxLayer") << "): selClusterIdx would point at the wrong rows";
      scopeChecked_ = true;
    }
    const auto order = smartpixels::clusterTableOrder(recHits, topo, maxLayer_, maxDisk_);
    if (order.size() != recHits.size())
      throw cms::Exception("SmartPixelsClusterTableScope")
          << recHits.size() - order.size() << " rec-hit modules lie outside the table scope; the table must "
          << "emit every rec hit for selClusterIdx to be a row index";
    std::vector<uint8_t> subdet, layer, disk, side, ring, surface, sizeX, sizeY;
    std::vector<uint16_t> size, idxInDet;
    std::vector<uint32_t> detId;
    std::vector<float> localX, localY, sigX, sigY, charge;
    std::vector<float> globalR, globalPhi, globalZ;
    std::vector<float> globalClusterPhi, globalClusterCotTheta, tpGlobalClusterPhi, tpGlobalClusterCotTheta;
    std::vector<float> sigGlobalClusterPhi, sigGlobalClusterCotTheta;
    unsigned closureFail = 0;
    std::vector<float> localCotAlpha, localCotBeta, sigAlpha, sigBeta;
    std::vector<uint8_t> hasAlpha, hasBeta, satAlpha, satBeta;
    std::vector<int16_t> tpEmuClass;
    std::vector<float> tpPt, tpVx, tpVy, tpVz, tpEta, tpPhi, tpChargeFrac, tpLocalCotAlpha, tpLocalCotBeta;
    std::vector<int32_t> tpIdx;

    for (const auto& [dsv, surf] : order) {
      const DetId did(dsv.detId());
      const bool disc = smartpixels::isDiscSurface(surf);
      const auto* pdu = dynamic_cast<const PixelGeomDetUnit*>(geom.idToDet(did));
      const auto tsv = truthColl.find(dsv.detId());
      const bool haveTruth = doTruth_ && tsv != truthColl.end() && tsv->size() == dsv.size();
      if (doTruth_ && !haveTruth)
        throw cms::Exception("SmartPixelsRecHitTruthMisaligned")
            << "cluster table: rechit/truth disagree on det " << dsv.detId();

      for (size_t j = 0; j < dsv.size(); ++j) {
        const auto& rh = dsv[j];
        subdet.push_back(static_cast<uint8_t>(did.subdetId()));
        layer.push_back(static_cast<uint8_t>(disc ? 0 : surf));
        disk.push_back(static_cast<uint8_t>(disc ? topo.pxfDisk(did) : 0));
        side.push_back(static_cast<uint8_t>(disc ? topo.pxfSide(did) : 0));
        ring.push_back(static_cast<uint8_t>(disc ? topo.pxfBlade(did) : 0));
        surface.push_back(static_cast<uint8_t>(surf));
        detId.push_back(did.rawId());
        idxInDet.push_back(static_cast<uint16_t>(j));
        // GLOBAL position and direction, computed HERE because this is where the
        // geometry is. nano carries none, so without these every analysis would
        // re-derive them from detId with its own copy of the tracker geometry --
        // duplicated work and a silent divergence risk on TILTED modules (TBPX
        // tilt reaches 16.5 deg, so the rotation is per-module, not per-layer).
        //
        // CYLINDRICAL, not Cartesian. globalR and globalPhi are stored at full float: at 16
        // mantissa bits (canon_v3) the phi step reached 4.9 um at IL4, up to 0.53x the CPE
        // sigX the census KF weights hits with, and the r step reached 0.14x sigY on discs,
        // where r is the measured coordinate.
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
          // otherwise be a large, silent error. The absolute tolerance is meaningless
          // at grazing |cot| (the emu source reports up to ~5730), so the test keeps
          // the range it had when hasAlpha/hasBeta still carried the |cot| <= 12 cap.
          if (rh.hasAlpha() && rh.hasBeta() && std::abs(rh.cotAlpha()) <= kClosureMaxAbsCot &&
              std::abs(rh.cotBeta()) <= kClosureMaxAbsCot &&
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
        satAlpha.push_back(rh.satAlpha() ? 1 : 0);
        satBeta.push_back(rh.satBeta() ? 1 : 0);

        if (haveTruth) {
          const auto& tr = (*tsv)[j];
          tpIdx.push_back(l1tnano::tpIdx(tpRowMap, tr.dominantTp()));
          tpPt.push_back(tr.hasTp() ? static_cast<float>(tr.dominantTp()->pt()) : -999.f);
          tpVx.push_back(tr.hasTp() ? static_cast<float>(tr.dominantTp()->vx()) : -999.f);
          tpVy.push_back(tr.hasTp() ? static_cast<float>(tr.dominantTp()->vy()) : -999.f);
          tpVz.push_back(tr.hasTp() ? static_cast<float>(tr.dominantTp()->vz()) : -999.f);
          tpEta.push_back(tr.hasTp() ? static_cast<float>(tr.dominantTp()->eta()) : -999.f);
          tpPhi.push_back(tr.hasTp() ? static_cast<float>(tr.dominantTp()->phi()) : -999.f);
          tpChargeFrac.push_back(tr.chargeFrac());
          tpLocalCotAlpha.push_back(tr.trueCotAlpha());
          tpLocalCotBeta.push_back(tr.trueCotBeta());
          tpEmuClass.push_back(tr.emuClass());
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

    auto tab = std::make_unique<nanoaod::FlatTable>(detId.size(), tableName_, false, false);
    tab->addColumn<uint8_t>("subdet", subdet, "1 = TBPX barrel, 2 = pixel disc");
    tab->addColumn<uint8_t>("layer", layer, "TBPX layer 1..4; 0 on a disc");
    tab->addColumn<uint8_t>("disk", disk, "pxfDisk (1-8 TFPX, 9-12 TEPX); 0 on the barrel");
    tab->addColumn<uint8_t>("side", side, "pxfSide: 1 = -z, 2 = +z; 0 on the barrel");
    tab->addColumn<uint8_t>("ring", ring, "ring (pxfBlade in the Phase-2 numbering); 0 on the barrel");
    tab->addColumn<uint8_t>("surface", surface,
                            "IT surface: 1-4 = IL1-IL4 (TBPX layer), 10 + disk on the discs (ID1-ID4 = 11-14). "
                            "Rows are all barrel clusters, then all disc clusters");
    tab->addColumn<uint32_t>("detId", detId, "module rawId; join key to the refit hit table detId");
    tab->addColumn<uint16_t>("idxInDet", idxInDet,
                             "position of the cluster within its module; (detId, idxInDet) is a stable key");
    tab->addColumn<float>("localX", localX, "cluster position, module-local x [cm]", bits(16));
    tab->addColumn<float>("localY", localY,
                          "cluster position, module-local y [cm]; along z on TBPX, along the radius "
                          "(pointing inward) on a disc", bits(16));
    tab->addColumn<float>("globalR", globalR,
                          "POSITION. Cylindrical radius of WHERE THE CLUSTER IS [cm], full float");
    tab->addColumn<float>("globalPhi", globalPhi,
                          "POSITION. Azimuth of WHERE THE CLUSTER IS, in CMS global coordinates "
                          "[rad]. This is a location on the detector surface and has NOTHING to do "
                          "with any direction or angle estimate: it is atan2(y,x) of the cluster "
                          "centroid. Pairs with globalR and globalZ to give the full position. "
                          "Do NOT confuse with globalClusterPhi, which is a DIRECTION. Full float");
    tab->addColumn<float>("globalZ", globalZ,
                          "POSITION. CMS global z of WHERE THE CLUSTER IS [cm]", bits(16));
    tab->addColumn<float>("globalClusterPhi", globalClusterPhi,
                          "DIRECTION. Azimuth of the SMART-PIXEL ML ANGLE ESTIMATE for this "
                          "cluster, rotated into CMS global coordinates [rad]. This is the "
                          "estimated direction of the PARTICLE that made the cluster -- the ML "
                          "regressor output, carrying the PixelAV response -- NOT a position. It "
                          "is the global counterpart of localCotAlpha/localCotBeta, which are "
                          "module-frame by definition (PixelAV) and have no meaningful global "
                          "variant. Do NOT confuse with globalPhi, which is WHERE the cluster is. "
                          "Uses the per-module surface rotation (TBPX tilt reaches 16.5 deg, so it "
                          "is not a per-layer constant). Uncertainty: sigGlobalClusterPhi", bits(16));
    tab->addColumn<float>("globalClusterCotTheta", globalClusterCotTheta,
                          "DIRECTION. cot(theta) = pz/pt of the SMART-PIXEL ML ANGLE ESTIMATE, CMS "
                          "global frame. Chosen over eta because the r-z Hough wants "
                          "z = z0 + r*cotTheta directly; measured median RMS of z - r*cotTheta "
                          "along one track is 0.20 cm, against 3.47 cm for using z alone. SENSE is "
                          "supplied by physics, not by the module frame: cotAlpha/cotBeta are ratios "
                          "and so are blind to direction reversal, and the sign is fixed by requiring "
                          "the direction to point AWAY from the beam axis on TBPX and AWAY from z = 0 "
                          "on the discs. TBPX tracks that curl back inward violate that and keep an "
                          "irrecoverable sign ambiguity -- 15% of tracks below 0.5 GeV, falling to 2% "
                          "above 5 GeV. "
                          "Uncertainty: sigGlobalClusterCotTheta", bits(16));
    tab->addColumn<float>("sigGlobalClusterPhi", sigGlobalClusterPhi,
                          "Uncertainty on globalClusterPhi [rad], propagated from the module-frame "
                          "sigAlpha/sigBeta through the SAME per-module rotation by a numerical "
                          "Jacobian, added in quadrature. Alpha and beta are treated as "
                          "independent, matching how the refit applies them as two independent "
                          "scalar Kalman updates", bits(12));
    tab->addColumn<float>("sigGlobalClusterCotTheta", sigGlobalClusterCotTheta,
                          "Uncertainty on globalClusterCotTheta, same propagation", bits(12));
    tab->addColumn<float>("sigX", sigX, "CPE position uncertainty, local x [cm]", bits(10));
    tab->addColumn<float>("sigY", sigY, "CPE position uncertainty, local y [cm] (sigma_r on a disc)", bits(10));
    tab->addColumn<uint8_t>("sizeX", sizeX, "cluster bounding-box extent in pixels, local x");
    tab->addColumn<uint8_t>("sizeY", sizeY, "cluster bounding-box extent in pixels, local y");
    tab->addColumn<uint16_t>("size", size, "FIRED-PIXEL COUNT (not the bounding box): charge/size is "
                                           "the real charge density");
    tab->addColumn<float>("charge", charge, "cluster charge [electrons] (SiPixelCluster::charge)", bits(10));
    tab->addColumn<float>("localCotAlpha", localCotAlpha,
                          "SENSOR angle estimate, module-local cotAlpha (-999 if none)", bits(12));
    tab->addColumn<float>("localCotBeta", localCotBeta,
                          "SENSOR angle estimate, cotBeta (signed; on a disc the sign flips with the module "
                          "facing, use globalClusterCotTheta to compare directions)", bits(12));
    tab->addColumn<float>("sigAlpha", sigAlpha, "angle-estimator resolution of localCotAlpha (a sigma of cotAlpha, dimensionless)", bits(10));
    tab->addColumn<float>("sigBeta", sigBeta, "angle-estimator resolution of localCotBeta (a sigma of cotBeta, dimensionless)", bits(10));
    tab->addColumn<uint8_t>("hasAlpha", hasAlpha,
                            "sensor reports an alpha at all, whatever its value: angleSource emu ALWAYS 1 (a "
                            "flag that depended on the value would tell class C from class T); payload 0 outside "
                            "the PixelAV validity gate and on discs. Algorithms apply their own |cot| cap (the "
                            "refit: measAngleMaxAbs = 12)");
    tab->addColumn<uint8_t>("hasBeta", hasBeta, "sensor reports a beta; as hasAlpha");
    tab->addColumn<uint8_t>("satAlpha", satAlpha,
                            "the emulated alpha (angleSource emu) hit its clamp at +-(90 - 0.01) deg: the "
                            "reported value is the bound (|cot| ~ 5730; hasAlpha stays 1); always 0 for the "
                            "payload source");
    tab->addColumn<uint8_t>("satBeta", satBeta, "as satAlpha, for beta");
    if (doTruth_) {
      tab->addColumn<int32_t>("tpIdx", tpIdx,
                              std::string("TRUTH-ONLY, dominant charge contributor: ") + l1tnano::tpIdxDoc());
      tab->addColumn<float>("tpPt", tpPt, "TRUTH-ONLY: pT [GeV] of the dominant TP", bits(10));
      // PRODUCTION STATE of the dominant TP, at reduced precision (vertex 16 bits, eta/phi 12
      // bits): enough to tell a displaced TP from one that merely scattered, NOT enough to
      // recompute a helix projection at IT or OT precision. For projections use the full-precision
      // L1TTP row (L1TTP[tpIdx]). tpLocalCot{Alpha,Beta} are helix-propagated in the producer and
      // ignore multiple scattering. NOTE these are at the PRODUCTION VERTEX, not at the module --
      // do not use them as an incidence angle.
      tab->addColumn<float>("tpVx", tpVx, "TRUTH-ONLY: dominant TP production vertex x [cm]", bits(16));
      tab->addColumn<float>("tpVy", tpVy, "TRUTH-ONLY: dominant TP production vertex y [cm]", bits(16));
      tab->addColumn<float>("tpVz", tpVz, "TRUTH-ONLY: dominant TP production vertex z [cm]", bits(16));
      tab->addColumn<float>("tpEta", tpEta,
                            "TRUTH-ONLY: dominant TP eta AT PRODUCTION, 12 bits; for projections use L1TTP[tpIdx]", bits(12));
      tab->addColumn<float>("tpPhi", tpPhi,
                            "TRUTH-ONLY: dominant TP phi AT PRODUCTION, 12 bits; for projections use L1TTP[tpIdx]", bits(12));
      tab->addColumn<float>("tpLocalCotAlpha", tpLocalCotAlpha,
                            "TRUTH-ONLY: TRUE incidence cotAlpha at this module (helix-propagated), "
                            "i.e. what the sensor is trying to measure", bits(12));
      tab->addColumn<float>("tpLocalCotBeta", tpLocalCotBeta, "TRUTH-ONLY: true incidence cotBeta", bits(12));
      tab->addColumn<float>("tpGlobalClusterPhi", tpGlobalClusterPhi,
                            "TRUTH-ONLY: true direction, global phi [rad]", bits(16));
      tab->addColumn<float>("tpGlobalClusterCotTheta", tpGlobalClusterCotTheta,
                            "TRUTH-ONLY: true direction, global cot(theta)", bits(16));
      tab->addColumn<float>("tpChargeFrac", tpChargeFrac,
                            "TRUTH-ONLY: dominant contributor share of the cluster charge", bits(10));
      tab->addColumn<int16_t>("tpEmuClass", tpEmuClass,
                             "TRUTH-ONLY, VALIDATION ONLY (a truth tell by construction; never a feature): "
                             "SmartPixels-output emulator class, 1 = T (angle emulated from the dominant TP's "
                             "truth), 0 = C (CNN + drawn residual), -1 = not the emu angle source");
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
    desc.add<std::string>("discScope", "sp")
        ->setComment("disc rows: none | sp (ID1-ID4) | tfpx | all. MUST equal the rec hits' discScope (checked)");
    desc.add<unsigned>("maxLayer", 4)
        ->setComment("highest TBPX layer kept (SmartPixels scope is 1..4); MUST equal the rec hits' (checked)");
    desc.add<bool>("closurePrecision", false)
        ->setComment("TEST ONLY (refit v3 closure): every float column at full precision");
    desc.add<bool>("doTruth", true)
        ->setComment("attach the TRUTH-ONLY tp* block (tpIdx, tpPt, tpChargeFrac, "
                     "tpLocalCot*, tpGlobalCluster*)");
    desc.add<edm::InputTag>("tpRowMap", edm::InputTag("l1tTrackingParticleTable"))
        ->setComment("L1TTP key -> row map from L1TrackingParticleTableProducer; tpIdx = -2 when absent");
    descriptions.addWithDefaultLabel(desc);
  }

private:
  const edm::EDGetTokenT<SmartPixelsRecHitCollection> recHitToken_;
  const edm::EDGetTokenT<SmartPixelsRecHitTruthCollection> truthToken_;
  const edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> topoToken_;
  const edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> geomToken_;
  const std::string tableName_;
  const std::string discScope_;
  const unsigned maxLayer_;
  const unsigned maxDisk_;
  const bool doTruth_;
  const bool closurePrecision_;
  int bits(int n) const { return closurePrecision_ ? -1 : n; }
  edm::EDGetTokenT<edm::ValueMap<int>> tpRowMapToken_;
  bool scopeChecked_ = false;
};

DEFINE_FWK_MODULE(L1SmartPixelsClusterTableProducer);
