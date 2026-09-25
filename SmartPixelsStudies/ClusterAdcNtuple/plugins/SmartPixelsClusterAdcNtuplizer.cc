// SmartPixelsClusterAdcNtuplizer: one TTree row per TBPX (L1-L4) pixel cluster, carrying
// the SPARSE pixel map plus everything needed to interpret it and to train a
// per-cluster angle/position regressor on TrackingParticle truth.
//
// SCOPE. A study ntuple, deliberately NOT a nano table. It reads the SAME
// spixSmartPixelsRecHits (SmartPixelsRecHit + SmartPixelsRecHitTruth) that feed the
// L1TSmartPixelsCluster nano table, iterated in the same order, so row k of an event
// here is row k of that table (clIdx), and (event, detId, localX, localY) agree.
//
// ROWS: every cluster on TBPX layers 1..maxLayer -- TP-linked, simlinked-without-TP
// and simlink-free, any pT. No truncation of the pixel list.
//
// TRUTH.
//  * tpLocalCotAlpha/Beta: copied from SmartPixelsRecHitTruth, i.e. the producer's
//    helix-propagated direction evaluated at the RECO hit position (identical to the
//    nano column of the same name).
//  * hx*: the dominant TP's helix intersected with the module mid-plane (local z = 0):
//    true local x/y and the direction there. Pure helix from the TP production state
//    in the uniform field bz(0,0,0) -- the producer's model -- so no multiple
//    scattering and no energy loss.
//  * sh*: the highest-|p| PSimHit of the cluster's dominant SimTrack on this module, when one
//    exists (g4SimHits = SIGNAL crossing only, ~2% of PU200 clusters): entry/exit
//    derived direction and position interpolated to local z = 0.
//
// LINK CLASS (linkClass): 0 = no pixel carries a PixelDigiSimLink; 1 = simlinked, but
// no contributing SimTrack belongs to a TrackingParticle (e.g. sub-threshold particles);
// 2 = TP-linked (tpIdx >= 0). Classification logic is the producer's (winner-takes-pixel
// on the link fraction, charge summed per TP index) and is cross-checked against the
// producer's own TP choice (mismatch count printed at endJob; must be 0).

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CondFormats/DataRecord/interface/SiPixelLorentzAngleRcd.h"
#include "CondFormats/DataRecord/interface/SiPixelLorentzAngleSimRcd.h"
#include "CondFormats/SiPixelObjects/interface/SiPixelLorentzAngle.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Phase3SmartPixels/interface/SmartPixelsRecHit.h"
#include "DataFormats/Phase3SmartPixels/interface/SmartPixelsRecHitTruth.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "Geometry/CommonTopologies/interface/PixelGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/TrackerDigiSimLink/interface/PixelDigiSimLink.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"

#include "TTree.h"

#include <cmath>
#include <cstdint>
#include <map>
#include <unordered_map>
#include <vector>

namespace {
  // Pure helix of a TrackingParticle from its production state in a uniform field
  // along z. Parameter s = TRANSVERSE arc length from the production vertex. Same
  // centre/sign convention as SmartPixelsRecHitProducer.
  struct Helix {
    bool curved = false;
    double R = 0, sgn = 0, cx = 0, cy = 0, th0 = 0;
    double vx = 0, vy = 0, vz = 0, ux = 0, uy = 0, cotTh = 0;

    Helix(const TrackingParticle& tp, double bz) {
      const double pt = std::hypot(tp.px(), tp.py());
      vx = tp.vx();
      vy = tp.vy();
      vz = tp.vz();
      ux = pt > 0 ? tp.px() / pt : 1.;
      uy = pt > 0 ? tp.py() / pt : 0.;
      cotTh = pt > 0 ? tp.pz() / pt : 0.;
      if (tp.charge() != 0 && pt > 1e-6 && std::abs(bz) > 1e-6) {
        curved = true;
        R = pt / (0.29979246 * std::abs(bz)) * 100.0;  // cm
        sgn = (tp.charge() * bz > 0.) ? +1.0 : -1.0;
        cx = vx + sgn * R * uy;
        cy = vy - sgn * R * ux;
        th0 = std::atan2(vy - cy, vx - cx);
      }
    }
    GlobalPoint pos(double s) const {
      if (!curved)
        return GlobalPoint(vx + s * ux, vy + s * uy, vz + s * cotTh);
      const double th = th0 - sgn * s / R;
      return GlobalPoint(cx + R * std::cos(th), cy + R * std::sin(th), vz + s * cotTh);
    }
    GlobalVector dir(double s) const {  // d(pos)/ds; transverse part has unit length
      if (!curved)
        return GlobalVector(ux, uy, cotTh);
      const double th = th0 - sgn * s / R;
      return GlobalVector(sgn * std::sin(th), -sgn * std::cos(th), cotTh);
    }
  };
}  // namespace

class SmartPixelsClusterAdcNtuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit SmartPixelsClusterAdcNtuplizer(const edm::ParameterSet&);
  static void fillDescriptions(edm::ConfigurationDescriptions&);
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

private:
  void book();

  const edm::EDGetTokenT<SmartPixelsRecHitCollection> hitToken_;
  const edm::EDGetTokenT<SmartPixelsRecHitTruthCollection> truthToken_;
  const edm::EDGetTokenT<edm::DetSetVector<PixelDigiSimLink>> simLinkToken_;
  const edm::EDGetTokenT<edm::DetSetVector<PixelDigi>> digiToken_;
  const edm::EDGetTokenT<std::vector<TrackingParticle>> tpToken_;
  std::vector<edm::EDGetTokenT<std::vector<PSimHit>>> simHitTokens_;
  const edm::EDGetTokenT<edm::SimTrackContainer> simTrackToken_;
  const edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> geomToken_;
  const edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> topoToken_;
  const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> fieldToken_;
  edm::ESGetToken<SiPixelLorentzAngle, SiPixelLorentzAngleRcd> laToken_;
  edm::ESGetToken<SiPixelLorentzAngle, SiPixelLorentzAngleSimRcd> laSimToken_;
  const unsigned maxLayer_;
  const bool doLorentz_;
  const double helixStep_;

  unsigned long long nClusters_ = 0, nTpMismatch_ = 0, nHelixFail_ = 0;

  TTree* tree_ = nullptr;
  // --- identity
  unsigned run_, lumi_;
  unsigned long long event_;
  int clIdx_, idxInDet_;
  unsigned detId_;
  int layer_, ladder_, module_;
  // --- reco
  float localX_, localY_, sigX_, sigY_, gX_, gY_, gZ_, gR_, gPhi_;
  int sizeX_, sizeY_, size_;
  float charge_;
  int minRow_, minCol_;
  float originLocalX_, originLocalY_;  // local coords of the (minRow, minCol) pixel CORNER
  std::vector<short> pxRow_, pxCol_, pxAdc_;
  std::vector<int> pxQ_;
  // --- module
  float pitchX_, pitchY_, thickness_;
  int nRows_, nCols_;
  float modX_, modY_, modZ_;
  float axXx_, axXy_, axXz_, axYx_, axYy_, axYz_, axZx_, axZy_, axZz_;
  float bModX_, bModY_, bModZ_, laPerT_, driftX_, driftY_, driftZ_;
  float laSimPerT_, tanLAx_, tanLAy_;  // digitizer: SimRcd LA, module-centre B
  float bLocX_, bLocY_, bLocZ_, bMag_;  // B at the CLUSTER position, local frame
  int zOutward_, is3D_;
  // --- payload (synthesised) angle, for the baseline
  float pCotA_, pCotB_, pSigA_, pSigB_;
  int pHasA_, pHasB_;
  // --- link classification
  int linkClass_, nTp_, nSimTrk_;
  float qFracLinked_, qFracTp_;
  int domSimTrackId_, domSimEvent_, domSimBx_, domSimPdg_, domSimHasTp_;
  float domSimFrac_, domSimPt_, domSimEta_;
  // --- dominant TP
  int tpIdx_, tpPdgId_, tpCharge_, tpBx_, tpEvent_, tpMerged_;
  float tpPt_, tpEta_, tpPhi_, tpVx_, tpVy_, tpVz_, tpPx_, tpPy_, tpPz_, tpChargeFrac_, tpSecondFrac_;
  float tpLocalCotAlpha_, tpLocalCotBeta_;
  // --- helix at module plane
  int hxOk_, hxNRoots_;
  float hxLocalX_, hxLocalY_, hxCotAlpha_, hxCotBeta_, hxS_, hxDist_, hxTurns_;
  // --- PSimHit
  int shN_, shPdg_, shProcess_;
  float shLocalX_, shLocalY_, shCotAlpha_, shCotBeta_, shPabs_, shTof_, shEloss_, shPathZ_;
};

SmartPixelsClusterAdcNtuplizer::SmartPixelsClusterAdcNtuplizer(const edm::ParameterSet& cfg)
    : hitToken_(consumes<SmartPixelsRecHitCollection>(cfg.getParameter<edm::InputTag>("smartPixelsRecHits"))),
      truthToken_(consumes<SmartPixelsRecHitTruthCollection>(cfg.getParameter<edm::InputTag>("smartPixelsRecHits"))),
      simLinkToken_(consumes<edm::DetSetVector<PixelDigiSimLink>>(cfg.getParameter<edm::InputTag>("pixelDigiSimLink"))),
      digiToken_(consumes<edm::DetSetVector<PixelDigi>>(cfg.getParameter<edm::InputTag>("pixelDigis"))),
      tpToken_(consumes<std::vector<TrackingParticle>>(cfg.getParameter<edm::InputTag>("trackingParticles"))),
      simTrackToken_(consumes<edm::SimTrackContainer>(cfg.getParameter<edm::InputTag>("simTracks"))),
      geomToken_(esConsumes()),
      topoToken_(esConsumes()),
      fieldToken_(esConsumes()),
      maxLayer_(cfg.getParameter<unsigned>("maxLayer")),
      doLorentz_(cfg.getParameter<bool>("doLorentz")),
      helixStep_(cfg.getParameter<double>("helixStep")) {
  for (const auto& t : cfg.getParameter<std::vector<edm::InputTag>>("simHits"))
    simHitTokens_.push_back(consumes<std::vector<PSimHit>>(t));
  if (doLorentz_)
  {
    laToken_ = esConsumes(edm::ESInputTag("", cfg.getParameter<std::string>("lorentzAngleLabel")));
    laSimToken_ = esConsumes();
  }
  usesResource("TFileService");
  book();
}

void SmartPixelsClusterAdcNtuplizer::book() {
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("clusters", "one row per TBPX pixel cluster");
  auto I = [&](const char* n, int* v) { tree_->Branch(n, v, (std::string(n) + "/I").c_str()); };
  auto F = [&](const char* n, float* v) { tree_->Branch(n, v, (std::string(n) + "/F").c_str()); };
  tree_->Branch("run", &run_, "run/i");
  tree_->Branch("lumi", &lumi_, "lumi/i");
  tree_->Branch("event", &event_, "event/l");
  I("clIdx", &clIdx_);
  I("idxInDet", &idxInDet_);
  tree_->Branch("detId", &detId_, "detId/i");
  I("layer", &layer_);
  I("ladder", &ladder_);
  I("module", &module_);
  F("localX", &localX_);
  F("localY", &localY_);
  F("sigX", &sigX_);
  F("sigY", &sigY_);
  F("globalX", &gX_);
  F("globalY", &gY_);
  F("globalZ", &gZ_);
  F("globalR", &gR_);
  F("globalPhi", &gPhi_);
  I("sizeX", &sizeX_);
  I("sizeY", &sizeY_);
  I("size", &size_);
  F("charge", &charge_);
  I("minRow", &minRow_);
  I("minCol", &minCol_);
  F("originLocalX", &originLocalX_);
  F("originLocalY", &originLocalY_);
  tree_->Branch("pxRow", &pxRow_);
  tree_->Branch("pxCol", &pxCol_);
  tree_->Branch("pxQ", &pxQ_);
  tree_->Branch("pxAdc", &pxAdc_);
  F("pitchX", &pitchX_);
  F("pitchY", &pitchY_);
  F("thickness", &thickness_);
  I("nRows", &nRows_);
  I("nCols", &nCols_);
  F("modX", &modX_);
  F("modY", &modY_);
  F("modZ", &modZ_);
  F("axXx", &axXx_);
  F("axXy", &axXy_);
  F("axXz", &axXz_);
  F("axYx", &axYx_);
  F("axYy", &axYy_);
  F("axYz", &axYz_);
  F("axZx", &axZx_);
  F("axZy", &axZy_);
  F("axZz", &axZz_);
  I("zOutward", &zOutward_);
  I("is3D", &is3D_);
  F("bModLocalX", &bModX_);
  F("bModLocalY", &bModY_);
  F("bModLocalZ", &bModZ_);
  F("bLocalX", &bLocX_);
  F("bLocalY", &bLocY_);
  F("bLocalZ", &bLocZ_);
  F("bMag", &bMag_);
  F("laPerTesla", &laPerT_);
  F("laSimPerTesla", &laSimPerT_);
  F("tanLAx", &tanLAx_);
  F("tanLAy", &tanLAy_);
  F("driftX", &driftX_);
  F("driftY", &driftY_);
  F("driftZ", &driftZ_);
  F("pCotAlpha", &pCotA_);
  F("pCotBeta", &pCotB_);
  F("pSigAlpha", &pSigA_);
  F("pSigBeta", &pSigB_);
  I("pHasAlpha", &pHasA_);
  I("pHasBeta", &pHasB_);
  I("linkClass", &linkClass_);
  I("nTp", &nTp_);
  I("nSimTrk", &nSimTrk_);
  F("qFracLinked", &qFracLinked_);
  F("qFracTp", &qFracTp_);
  I("domSimTrackId", &domSimTrackId_);
  I("domSimEvent", &domSimEvent_);
  I("domSimBx", &domSimBx_);
  I("domSimPdg", &domSimPdg_);
  I("domSimHasTp", &domSimHasTp_);
  F("domSimFrac", &domSimFrac_);
  F("domSimPt", &domSimPt_);
  F("domSimEta", &domSimEta_);
  I("tpIdx", &tpIdx_);
  I("tpPdgId", &tpPdgId_);
  I("tpCharge", &tpCharge_);
  I("tpBx", &tpBx_);
  I("tpEvent", &tpEvent_);
  I("tpMerged", &tpMerged_);
  F("tpPt", &tpPt_);
  F("tpEta", &tpEta_);
  F("tpPhi", &tpPhi_);
  F("tpVx", &tpVx_);
  F("tpVy", &tpVy_);
  F("tpVz", &tpVz_);
  F("tpPx", &tpPx_);
  F("tpPy", &tpPy_);
  F("tpPz", &tpPz_);
  F("tpChargeFrac", &tpChargeFrac_);
  F("tpSecondFrac", &tpSecondFrac_);
  F("tpLocalCotAlpha", &tpLocalCotAlpha_);
  F("tpLocalCotBeta", &tpLocalCotBeta_);
  I("hxOk", &hxOk_);
  I("hxNRoots", &hxNRoots_);
  F("hxLocalX", &hxLocalX_);
  F("hxLocalY", &hxLocalY_);
  F("hxCotAlpha", &hxCotAlpha_);
  F("hxCotBeta", &hxCotBeta_);
  F("hxS", &hxS_);
  F("hxDist", &hxDist_);
  F("hxTurns", &hxTurns_);
  I("shN", &shN_);
  I("shPdg", &shPdg_);
  I("shProcess", &shProcess_);
  F("shLocalX", &shLocalX_);
  F("shLocalY", &shLocalY_);
  F("shCotAlpha", &shCotAlpha_);
  F("shCotBeta", &shCotBeta_);
  F("shPabs", &shPabs_);
  F("shTof", &shTof_);
  F("shEloss", &shEloss_);
  F("shPathZ", &shPathZ_);
}

void SmartPixelsClusterAdcNtuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  const auto& topo = iSetup.getData(topoToken_);
  const auto& geom = iSetup.getData(geomToken_);
  const auto& field = iSetup.getData(fieldToken_);
  const SiPixelLorentzAngle* la = doLorentz_ ? &iSetup.getData(laToken_) : nullptr;
  const SiPixelLorentzAngle* laSim = doLorentz_ ? &iSetup.getData(laSimToken_) : nullptr;
  const auto& hits = iEvent.get(hitToken_);
  const auto& truths = iEvent.get(truthToken_);
  const auto& simLinks = iEvent.get(simLinkToken_);
  const auto& digis = iEvent.get(digiToken_);
  const auto& tps = iEvent.get(tpToken_);
  const auto& simTracks = iEvent.get(simTrackToken_);
  const double bz = field.inTesla(GlobalPoint(0, 0, 0)).z();

  std::map<std::pair<unsigned, unsigned>, int> tpIndex;  // (eventId raw, g4 trackId) -> TP key
  for (size_t i = 0; i < tps.size(); ++i)
    for (const auto& g4 : tps[i].g4Tracks())
      tpIndex.emplace(std::make_pair(tps[i].eventId().rawId(), g4.trackId()), static_cast<int>(i));
  std::unordered_map<unsigned, const SimTrack*> sigSimTrack;  // signal-crossing SimTracks by id
  for (const auto& st : simTracks)
    if (st.eventId().event() == 0 && st.eventId().bunchCrossing() == 0)
      sigSimTrack[st.trackId()] = &st;
  std::unordered_map<unsigned, std::vector<const PSimHit*>> simHitsByDet;
  for (const auto& tok : simHitTokens_) {
    edm::Handle<std::vector<PSimHit>> h;
    iEvent.getByToken(tok, h);
    if (h.isValid())
      for (const auto& sh : *h)
        simHitsByDet[sh.detUnitId()].push_back(&sh);
  }

  run_ = iEvent.id().run();
  lumi_ = iEvent.id().luminosityBlock();
  event_ = iEvent.id().event();
  int clIdx = 0;

  for (const auto& dsv : hits) {
    const DetId did(dsv.detId());
    if (did.subdetId() != PixelSubdetector::PixelBarrel)
      continue;
    const unsigned lay = topo.pxbLayer(did);
    if (lay < 1 || lay > maxLayer_)
      continue;
    const auto* pdu = dynamic_cast<const PixelGeomDetUnit*>(geom.idToDet(did));
    if (pdu == nullptr)
      continue;
    const auto tsv = truths.find(dsv.detId());
    if (tsv == truths.end() || tsv->size() != dsv.size())
      throw cms::Exception("ClusterAdcNtuple") << "rechit/truth misaligned on det " << dsv.detId();
    const auto truthSet = *tsv;  // DetSet view into the truth collection

    // ---- per-module constants
    const auto& ptopo = pdu->specificTopology();
    const auto& surf = pdu->surface();
    const auto pitch = ptopo.pitch();
    detId_ = did.rawId();
    layer_ = lay;
    ladder_ = topo.pxbLadder(did);
    module_ = topo.pxbModule(did);
    pitchX_ = pitch.first;
    pitchY_ = pitch.second;
    thickness_ = surf.bounds().thickness();
    nRows_ = ptopo.nrows();
    nCols_ = ptopo.ncolumns();
    const GlobalPoint mc = surf.position();
    modX_ = mc.x();
    modY_ = mc.y();
    modZ_ = mc.z();
    const GlobalVector ex = surf.toGlobal(LocalVector(1, 0, 0)), ey = surf.toGlobal(LocalVector(0, 1, 0)),
                       ez = surf.toGlobal(LocalVector(0, 0, 1));
    axXx_ = ex.x(); axXy_ = ex.y(); axXz_ = ex.z();
    axYx_ = ey.x(); axYy_ = ey.y(); axYz_ = ey.z();
    axZx_ = ez.x(); axZy_ = ez.y(); axZz_ = ez.z();
    zOutward_ = (ez.x() * mc.x() + ez.y() * mc.y()) > 0 ? 1 : -1;
    // 3D sensor => Pixel3DDigitizerAlgorithm (usePseudoPixel3DAlgo = False by default)
    is3D_ = geom.getDetectorType(did) == TrackerGeometry::ModuleType::Ph2PXB3D ? 1 : 0;
    // Field at the MODULE CENTRE: the point both the Phase-2 digitizer
    // (Phase2TrackerDigitizer::accumulate) and the CPE (PixelCPEBase) evaluate it at.
    const LocalVector bl = surf.toLocal(field.inTesla(mc));
    bModX_ = bl.x(); bModY_ = bl.y(); bModZ_ = bl.z();
    // Digitizer drift: Phase2TrackerDigitizerAlgorithm::driftDirection with the
    // SiPixelLorentzAngleSimRcd value (LorentzAngle_DB = True for the IT), Alpha2Order on;
    // tan components = dir/scale (TanLorenzAngleY used because Alpha2Order is on).
    laSimPerT_ = laSim ? laSim->getLorentzAngle(did.rawId()) : -999.f;
    if (laSim) {
      const float a2 = laSimPerT_ * laSimPerT_;
      const float scale = 1.f + a2 * bl.z() * bl.z();
      tanLAx_ = -(laSimPerT_ * bl.y() + a2 * bl.z() * bl.x()) / scale;
      tanLAy_ = (laSimPerT_ * bl.x() - a2 * bl.z() * bl.y()) / scale;
    } else {
      tanLAx_ = tanLAy_ = -999.f;
    }
    // Drift direction exactly as PixelCPEBase::driftDirection (alpha2Order on).
    laPerT_ = la ? la->getLorentzAngle(did.rawId()) : -999.f;
    if (la) {
      const float a2 = laPerT_ * laPerT_;
      driftX_ = -(laPerT_ * bl.y() + a2 * bl.z() * bl.x());
      driftY_ = (laPerT_ * bl.x() - a2 * bl.z() * bl.y());
      driftZ_ = -(1.f + a2 * bl.z() * bl.z());
    } else {
      driftX_ = driftY_ = driftZ_ = -999.f;
    }

    std::unordered_map<unsigned, const PixelDigiSimLink*> linkByChannel;
    const auto dsl = simLinks.find(did);
    if (dsl != simLinks.end())
      for (const auto& lk : *dsl) {
        auto it = linkByChannel.find(lk.channel());
        if (it == linkByChannel.end() || it->second->fraction() < lk.fraction())
          linkByChannel[lk.channel()] = &lk;
      }
    std::unordered_map<unsigned, int> adcByChannel;
    const auto ddi = digis.find(did);
    if (ddi != digis.end())
      for (const auto& d : *ddi)
        adcByChannel[d.channel()] = d.adc();
    const auto shIt = simHitsByDet.find(did.rawId());

    for (size_t j = 0; j < dsv.size(); ++j, ++clIdx) {
      const auto& rh = dsv[j];
      const auto& tr = truthSet[j];
      const SiPixelCluster& cl = *rh.cluster();
      ++nClusters_;
      clIdx_ = clIdx;
      idxInDet_ = tr.idxInDet();

      localX_ = rh.localPosition().x();
      localY_ = rh.localPosition().y();
      sigX_ = std::sqrt(std::max(0.f, float(rh.localPositionError().xx())));
      sigY_ = std::sqrt(std::max(0.f, float(rh.localPositionError().yy())));
      const GlobalPoint gp = surf.toGlobal(rh.localPosition());
      gX_ = gp.x(); gY_ = gp.y(); gZ_ = gp.z();
      const GlobalVector bg = field.inTesla(gp);
      const LocalVector bc = surf.toLocal(bg);
      bLocX_ = bc.x(); bLocY_ = bc.y(); bLocZ_ = bc.z();
      bMag_ = bg.mag();
      gR_ = std::hypot(gp.x(), gp.y());
      gPhi_ = std::atan2(gp.y(), gp.x());
      sizeX_ = cl.sizeX();
      sizeY_ = cl.sizeY();
      size_ = cl.size();
      charge_ = cl.charge();
      minRow_ = cl.minPixelRow();
      minCol_ = cl.minPixelCol();
      const LocalPoint org = ptopo.localPosition(MeasurementPoint(minRow_, minCol_));
      originLocalX_ = org.x();
      originLocalY_ = org.y();

      pCotA_ = rh.cotAlpha(); pCotB_ = rh.cotBeta();
      pSigA_ = rh.sigAlpha(); pSigB_ = rh.sigBeta();
      pHasA_ = rh.hasAlpha(); pHasB_ = rh.hasBeta();

      // ---- pixels + truth aggregation (producer logic, re-derived)
      pxRow_.clear(); pxCol_.clear(); pxQ_.clear(); pxAdc_.clear();
      std::map<int, double> qByTp;
      std::map<std::pair<unsigned, unsigned>, double> qBySim;
      double qTot = 0, qLinked = 0, qTpSum = 0;
      for (const auto& px : cl.pixels()) {
        pxRow_.push_back(static_cast<short>(px.x - minRow_));
        pxCol_.push_back(static_cast<short>(px.y - minCol_));
        pxQ_.push_back(px.adc);
        const unsigned ch = PixelDigi::pixelToChannel(px.x, px.y);
        const auto ai = adcByChannel.find(ch);
        pxAdc_.push_back(ai == adcByChannel.end() ? -1 : static_cast<short>(ai->second));
        qTot += px.adc;
        const auto lit = linkByChannel.find(ch);
        if (lit == linkByChannel.end())
          continue;
        qLinked += px.adc;
        const auto key = std::make_pair(lit->second->eventId().rawId(), lit->second->SimTrackId());
        qBySim[key] += px.adc;
        const auto tit = tpIndex.find(key);
        if (tit != tpIndex.end()) {
          qByTp[tit->second] += px.adc;
          qTpSum += px.adc;
        }
      }
      int myDom = -1;
      double qDom = 0, qSecond = 0;
      for (const auto& kv : qByTp) {
        if (kv.second > qDom) {
          qSecond = qDom;
          qDom = kv.second;
          myDom = kv.first;
        } else if (kv.second > qSecond)
          qSecond = kv.second;
      }
      nTp_ = qByTp.size();
      nSimTrk_ = qBySim.size();
      qFracLinked_ = qTot > 0 ? qLinked / qTot : 0.f;
      qFracTp_ = qTot > 0 ? qTpSum / qTot : 0.f;
      linkClass_ = qBySim.empty() ? 0 : (qByTp.empty() ? 1 : 2);

      // dominant SimTrack (any TP or none)
      std::pair<unsigned, unsigned> domSim{0, 0};
      double qDomSim = 0;
      for (const auto& kv : qBySim)
        if (kv.second > qDomSim) {
          qDomSim = kv.second;
          domSim = kv.first;
        }
      domSimTrackId_ = domSimEvent_ = domSimBx_ = -1;
      domSimPdg_ = 0;
      domSimHasTp_ = 0;
      domSimFrac_ = domSimPt_ = domSimEta_ = -999.f;
      EncodedEventId domSimEid;
      if (qDomSim > 0) {
        domSimEid = EncodedEventId(domSim.first);
        domSimTrackId_ = domSim.second;
        domSimEvent_ = domSimEid.event();
        domSimBx_ = domSimEid.bunchCrossing();
        domSimFrac_ = qDomSim / qTot;
        domSimHasTp_ = tpIndex.count(domSim) ? 1 : 0;
        if (domSimEvent_ == 0 && domSimBx_ == 0) {
          const auto sit = sigSimTrack.find(domSim.second);
          if (sit != sigSimTrack.end()) {
            domSimPdg_ = sit->second->type();
            domSimPt_ = sit->second->momentum().pt();
            domSimEta_ = sit->second->momentum().eta();
          }
        }
      }

      // ---- dominant TP: the producer's choice is authoritative (it defines tpIdx)
      tpIdx_ = tr.hasTp() ? static_cast<int>(tr.dominantTp().key()) : -1;
      if (tpIdx_ != myDom)
        ++nTpMismatch_;
      tpChargeFrac_ = tr.chargeFrac();
      tpSecondFrac_ = qTot > 0 ? qSecond / qTot : 0.f;
      tpMerged_ = tr.merged();
      tpLocalCotAlpha_ = tr.trueCotAlpha();
      tpLocalCotBeta_ = tr.trueCotBeta();
      hxOk_ = 0;
      hxNRoots_ = 0;
      hxLocalX_ = hxLocalY_ = hxCotAlpha_ = hxCotBeta_ = hxS_ = hxDist_ = hxTurns_ = -999.f;
      if (tpIdx_ >= 0) {
        const TrackingParticle& tp = tps[tpIdx_];
        tpPdgId_ = tp.pdgId();
        tpCharge_ = tp.charge();
        tpBx_ = tp.eventId().bunchCrossing();
        tpEvent_ = tp.eventId().event();
        tpPt_ = tp.pt(); tpEta_ = tp.eta(); tpPhi_ = tp.phi();
        tpVx_ = tp.vx(); tpVy_ = tp.vy(); tpVz_ = tp.vz();
        tpPx_ = tp.px(); tpPy_ = tp.py(); tpPz_ = tp.pz();

        // ---- helix x module mid-plane: scan + bisection, keep the root nearest the hit
        const Helix hx(tp, bz);
        const GlobalVector n = surf.normalVector();
        auto f = [&](double s) {
          const GlobalPoint p = hx.pos(s);
          return n.x() * (p.x() - mc.x()) + n.y() * (p.y() - mc.y()) + n.z() * (p.z() - mc.z());
        };
        const double sMax = hx.curved ? std::min(4. * M_PI * hx.R, 300.) : 300.;
        double bestD = 1e9, bestS = -1;
        double s0 = 0, f0 = f(0);
        while (s0 < sMax) {
          const double s1 = std::min(s0 + helixStep_, sMax);
          const double f1 = f(s1);
          if ((f0 <= 0 && f1 > 0) || (f0 >= 0 && f1 < 0)) {
            double a = s0, b = s1, fa = f0;
            for (int it = 0; it < 50; ++it) {
              const double m = 0.5 * (a + b), fm = f(m);
              if ((fa <= 0 && fm > 0) || (fa >= 0 && fm < 0))
                b = m;
              else {
                a = m;
                fa = fm;
              }
            }
            const double sr = 0.5 * (a + b);
            ++hxNRoots_;
            const double d = (hx.pos(sr) - gp).mag();
            if (d < bestD) {
              bestD = d;
              bestS = sr;
            }
          }
          s0 = s1;
          f0 = f1;
        }
        if (bestS >= 0) {
          hxOk_ = 1;
          const LocalPoint lp = surf.toLocal(hx.pos(bestS));
          const LocalVector ld = surf.toLocal(hx.dir(bestS));
          const double lz = std::abs(ld.z()) > 1e-9 ? ld.z() : 1e-9;
          hxLocalX_ = lp.x();
          hxLocalY_ = lp.y();
          hxCotAlpha_ = ld.x() / lz;
          hxCotBeta_ = ld.y() / lz;
          hxS_ = bestS;
          hxDist_ = bestD;
          hxTurns_ = hx.curved ? bestS / (2. * M_PI * hx.R) : 0.f;
        } else
          ++nHelixFail_;
      } else {
        tpPdgId_ = tpCharge_ = tpBx_ = tpEvent_ = 0;
        tpPt_ = tpEta_ = tpPhi_ = tpVx_ = tpVy_ = tpVz_ = tpPx_ = tpPy_ = tpPz_ = -999.f;
      }

      // ---- PSimHit of the dominant SimTrack on this module (signal only in practice).
      // Secondaries that Geant4 does not save (delta rays) inherit their ancestor's
      // trackId, so one det can hold several hits of "the same" track; the one with the
      // highest momentum is the primary crossing, not a soft delta segment.
      shN_ = 0;
      shPdg_ = shProcess_ = 0;
      shLocalX_ = shLocalY_ = shCotAlpha_ = shCotBeta_ = shPabs_ = shTof_ = shEloss_ = shPathZ_ = -999.f;
      if (qDomSim > 0 && shIt != simHitsByDet.end()) {
        const PSimHit* best = nullptr;
        double bd = 1e9;
        for (const PSimHit* sh : shIt->second) {
          if (sh->trackId() != domSim.second || sh->eventId().rawId() != domSim.first)
            continue;
          ++shN_;
          const double d = -sh->pabs();
          if (d < bd) {
            bd = d;
            best = sh;
          }
        }
        if (best) {
          const LocalPoint en = best->entryPoint(), exi = best->exitPoint();
          const double dz = exi.z() - en.z();
          const double t = std::abs(dz) > 1e-7 ? -en.z() / dz : 0.5;
          shLocalX_ = en.x() + t * (exi.x() - en.x());
          shLocalY_ = en.y() + t * (exi.y() - en.y());
          shCotAlpha_ = std::abs(dz) > 1e-7 ? (exi.x() - en.x()) / dz : -999.f;
          shCotBeta_ = std::abs(dz) > 1e-7 ? (exi.y() - en.y()) / dz : -999.f;
          shPabs_ = best->pabs();
          shTof_ = best->tof();
          shEloss_ = best->energyLoss();
          shPdg_ = best->particleType();
          shProcess_ = best->processType();
          shPathZ_ = dz;
        }
      }
      tree_->Fill();
    }
  }
}

void SmartPixelsClusterAdcNtuplizer::endJob() {
  edm::LogPrint("ClusterAdcNtuple") << "clusters=" << nClusters_ << " tpIdx mismatches vs producer=" << nTpMismatch_
                                    << " helix-plane failures=" << nHelixFail_;
}

void SmartPixelsClusterAdcNtuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("smartPixelsRecHits", edm::InputTag("spixSmartPixelsRecHits"));
  desc.add<edm::InputTag>("pixelDigiSimLink", edm::InputTag("simSiPixelDigis", "Pixel"));
  desc.add<edm::InputTag>("pixelDigis", edm::InputTag("simSiPixelDigis", "Pixel"));
  desc.add<edm::InputTag>("trackingParticles", edm::InputTag("mix", "MergedTrackTruth"));
  desc.add<edm::InputTag>("simTracks", edm::InputTag("g4SimHits"));
  desc.add<std::vector<edm::InputTag>>("simHits",
                                       {edm::InputTag("g4SimHits", "TrackerHitsPixelBarrelLowTof"),
                                        edm::InputTag("g4SimHits", "TrackerHitsPixelBarrelHighTof")});
  desc.add<unsigned>("maxLayer", 4);
  desc.add<bool>("doLorentz", true);
  desc.add<std::string>("lorentzAngleLabel", "");
  desc.add<double>("helixStep", 0.2)->setComment("transverse scan step [cm] for the helix-plane root search");
  descriptions.addWithDefaultLabel(desc);
}

DEFINE_FWK_MODULE(SmartPixelsClusterAdcNtuplizer);
