// -*- C++ -*-
//
// Package:    L1Trigger/Phase3SmartPixels
// Class:      SmartPixelsPayloadAnalyzer
//
/**\class SmartPixelsPayloadAnalyzer SmartPixelsPayloadAnalyzer.cc L1Trigger/Phase3SmartPixels/plugins/SmartPixelsPayloadAnalyzer.cc

 Description: Phase 1 of the Tier-2 (digiRefit) refit plan. Derives the raw
 inputs for the "smarthit_true" (Stack A) and "smarthit_fake" (Stack B)
 correctionlib payloads from a RelVal, by measuring — per L1Track crossing of an
 Inner-Tracker BARREL pixel layer — the in-window pixel-digi population and its
 per-digi truth class (via PixelDigiSimLink). Angle information is recorded in
 the module-LOCAL frame (PixelAV cotAlpha/cotBeta convention) together with the
 local B-field, so downstream payloads are correct on tilted/endcap modules.

 This is a PURE file-reading analyzer: it remakes nothing. It consumes the
 already-persisted TTTracks, TTTrackAssociationMap, TrackingParticles, the
 simSiPixelDigis:Pixel PixelDigi + PixelDigiSimLink collections, and the
 g4SimHits SimTracks/SimVertices, all present in RelValTTbar_D121_noPU.root.

 Propagation is intentionally LIGHT: analytic helix -> cylinder intersection to
 find the crossing phi/z, module lookup from TrackerGeometry, then transform to
 module-local coordinates. No reco propagators are pulled in (BuildFile stays
 lean). Material effects (multiple scattering, energy loss) are IGNORED for v0 —
 documented, and adequate for barrel window sizing at this stage.

 Local-frame conventions (module GeomDet frame):
   x_local = measurement (r-phi) direction; y_local = along module length
   (~z in barrel); z_local = outward sensor normal.
   PixelAV angle convention: cotAlpha = p_x_local / p_z_local,
                             cotBeta  = p_y_local / p_z_local
   (matches the CMS pixel-template cotalpha/cotbeta definition). The track's
   local direction uses the global helix momentum at the crossing; each linked
   digi's parent angle uses that SimTrack's ORIGIN momentum (origin-momentum
   approximation, adequate for v0 — noted; a per-crossing propagation of the
   parent is the upgrade path).

 Original Author:  Nick Manganelli (Phase-1 payload analyzer)
*/

// Framework
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// Data formats
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/GeometryVector/interface/LocalVector.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/MeasurementPoint.h"

#include "SimDataFormats/TrackerDigiSimLink/interface/PixelDigiSimLink.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "SimDataFormats/Associations/interface/TTTrackAssociationMap.h"

// Geometry + field
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/CommonTopologies/interface/PixelGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "Geometry/CommonTopologies/interface/GeomDet.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

// ROOT
#include "TTree.h"

#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <vector>
#include <string>

namespace {
  // Maximum number of in-window digis stored per crossing in the flat arrays.
  constexpr int kMaxWinDigis = 128;
  // PixelDigiSimLink truth classes.
  constexpr int kClassSameTP = 0;  // digi linked to the matched TP's g4 track(s)
  constexpr int kClassOtherTP = 1;  // digi linked to a DIFFERENT SimTrack
  constexpr int kClassNoise = 2;  // no usable link (noise / unlinked)
}  // namespace

class SmartPixelsPayloadAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit SmartPixelsPayloadAnalyzer(const edm::ParameterSet&);
  ~SmartPixelsPayloadAnalyzer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void analyze(const edm::Event&, const edm::EventSetup&) override;

  // ---- configuration ----
  const edm::EDGetTokenT<std::vector<TTTrack<Ref_Phase2TrackerDigi_>>> trackToken_;
  const edm::EDGetTokenT<TTTrackAssociationMap<Ref_Phase2TrackerDigi_>> assocToken_;
  const edm::EDGetTokenT<std::vector<TrackingParticle>> tpToken_;
  const edm::EDGetTokenT<edm::DetSetVector<PixelDigi>> digiToken_;
  const edm::EDGetTokenT<edm::DetSetVector<PixelDigiSimLink>> simlinkToken_;
  const edm::EDGetTokenT<edm::SimTrackContainer> simTrackToken_;
  const edm::EDGetTokenT<edm::SimVertexContainer> simVertexToken_;

  const edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> geomToken_;
  const edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> topoToken_;
  const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> fieldToken_;

  const double windowRPhi_;  // r-phi window half-width [cm] (local x)
  const double windowZ_;     // z window half-width [cm] (local y)
  const int nLayers_;        // number of TBPX layers to consider (1..4)
  const bool debug_;

  // ---- output tree (one row per track-layer crossing) ----
  TTree* tree_ = nullptr;

  // track-level
  int b_event_ = 0;
  float b_trk_pt_ = 0, b_trk_eta_ = 0, b_trk_phi_ = 0, b_trk_chi2_ = 0;
  int b_trk_charge_ = 0, b_trk_nstub_ = 0;
  int b_trk_tpMatched_ = 0;
  float b_tp_pt_ = 0, b_tp_eta_ = 0;

  // crossing-level (module-local)
  int b_layer_ = 0;
  unsigned int b_detid_ = 0;
  float b_exp_localx_ = 0, b_exp_localy_ = 0;      // expected crossing, local frame [cm]
  float b_trk_cotAlpha_ = 0, b_trk_cotBeta_ = 0;   // PixelAV convention, local frame
  float b_bx_ = 0, b_by_ = 0, b_bz_ = 0;           // local B-field [T]

  // per in-window digi arrays
  int b_nWinDigi_ = 0;
  int b_digi_class_[kMaxWinDigis];
  float b_digi_localx_[kMaxWinDigis], b_digi_localy_[kMaxWinDigis];
  float b_digi_resx_[kMaxWinDigis], b_digi_resy_[kMaxWinDigis];
  int b_digi_adc_[kMaxWinDigis];
  float b_digi_parCotAlpha_[kMaxWinDigis], b_digi_parCotBeta_[kMaxWinDigis];  // -999 if none
};

SmartPixelsPayloadAnalyzer::SmartPixelsPayloadAnalyzer(const edm::ParameterSet& iConfig)
    : trackToken_(consumes<std::vector<TTTrack<Ref_Phase2TrackerDigi_>>>(
          iConfig.getParameter<edm::InputTag>("l1TracksInputTag"))),
      assocToken_(consumes<TTTrackAssociationMap<Ref_Phase2TrackerDigi_>>(
          iConfig.getParameter<edm::InputTag>("mcTruthTrackInputTag"))),
      tpToken_(consumes<std::vector<TrackingParticle>>(
          iConfig.getParameter<edm::InputTag>("trackingParticleInputTag"))),
      digiToken_(consumes<edm::DetSetVector<PixelDigi>>(
          iConfig.getParameter<edm::InputTag>("pixelDigiInputTag"))),
      simlinkToken_(consumes<edm::DetSetVector<PixelDigiSimLink>>(
          iConfig.getParameter<edm::InputTag>("pixelDigiSimLinkInputTag"))),
      simTrackToken_(consumes<edm::SimTrackContainer>(
          iConfig.getParameter<edm::InputTag>("simTrackInputTag"))),
      simVertexToken_(consumes<edm::SimVertexContainer>(
          iConfig.getParameter<edm::InputTag>("simVertexInputTag"))),
      geomToken_(esConsumes<TrackerGeometry, TrackerDigiGeometryRecord>()),
      topoToken_(esConsumes<TrackerTopology, TrackerTopologyRcd>()),
      fieldToken_(esConsumes<MagneticField, IdealMagneticFieldRecord>()),
      windowRPhi_(iConfig.getParameter<double>("windowRPhi")),
      windowZ_(iConfig.getParameter<double>("windowZ")),
      nLayers_(iConfig.getParameter<int>("nLayers")),
      debug_(iConfig.getParameter<bool>("debug")) {
  usesResource(TFileService::kSharedResource);

  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("crossings", "SmartPixels payload: one row per L1Track x TBPX-layer crossing");

  tree_->Branch("event", &b_event_, "event/I");
  tree_->Branch("trk_pt", &b_trk_pt_, "trk_pt/F");
  tree_->Branch("trk_eta", &b_trk_eta_, "trk_eta/F");
  tree_->Branch("trk_phi", &b_trk_phi_, "trk_phi/F");
  tree_->Branch("trk_charge", &b_trk_charge_, "trk_charge/I");
  tree_->Branch("trk_chi2", &b_trk_chi2_, "trk_chi2/F");
  tree_->Branch("trk_nstub", &b_trk_nstub_, "trk_nstub/I");
  tree_->Branch("trk_tpMatched", &b_trk_tpMatched_, "trk_tpMatched/I");
  tree_->Branch("tp_pt", &b_tp_pt_, "tp_pt/F");
  tree_->Branch("tp_eta", &b_tp_eta_, "tp_eta/F");

  tree_->Branch("layer", &b_layer_, "layer/I");
  tree_->Branch("detid", &b_detid_, "detid/i");
  tree_->Branch("exp_localx", &b_exp_localx_, "exp_localx/F");
  tree_->Branch("exp_localy", &b_exp_localy_, "exp_localy/F");
  tree_->Branch("trk_cotAlpha", &b_trk_cotAlpha_, "trk_cotAlpha/F");
  tree_->Branch("trk_cotBeta", &b_trk_cotBeta_, "trk_cotBeta/F");
  tree_->Branch("b_localx", &b_bx_, "b_localx/F");
  tree_->Branch("b_localy", &b_by_, "b_localy/F");
  tree_->Branch("b_localz", &b_bz_, "b_localz/F");

  tree_->Branch("nWinDigi", &b_nWinDigi_, "nWinDigi/I");
  tree_->Branch("digi_class", b_digi_class_, "digi_class[nWinDigi]/I");
  tree_->Branch("digi_localx", b_digi_localx_, "digi_localx[nWinDigi]/F");
  tree_->Branch("digi_localy", b_digi_localy_, "digi_localy[nWinDigi]/F");
  tree_->Branch("digi_resx", b_digi_resx_, "digi_resx[nWinDigi]/F");
  tree_->Branch("digi_resy", b_digi_resy_, "digi_resy[nWinDigi]/F");
  tree_->Branch("digi_adc", b_digi_adc_, "digi_adc[nWinDigi]/I");
  tree_->Branch("digi_parCotAlpha", b_digi_parCotAlpha_, "digi_parCotAlpha[nWinDigi]/F");
  tree_->Branch("digi_parCotBeta", b_digi_parCotBeta_, "digi_parCotBeta[nWinDigi]/F");
}

void SmartPixelsPayloadAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  const TrackerGeometry& geom = iSetup.getData(geomToken_);
  const TrackerTopology& topo = iSetup.getData(topoToken_);
  const MagneticField& field = iSetup.getData(fieldToken_);

  edm::Handle<std::vector<TTTrack<Ref_Phase2TrackerDigi_>>> tracks;
  iEvent.getByToken(trackToken_, tracks);
  edm::Handle<TTTrackAssociationMap<Ref_Phase2TrackerDigi_>> assoc;
  iEvent.getByToken(assocToken_, assoc);
  edm::Handle<edm::DetSetVector<PixelDigi>> digis;
  iEvent.getByToken(digiToken_, digis);
  edm::Handle<edm::DetSetVector<PixelDigiSimLink>> simlinks;
  iEvent.getByToken(simlinkToken_, simlinks);
  edm::Handle<edm::SimTrackContainer> simTracks;
  iEvent.getByToken(simTrackToken_, simTracks);

  b_event_ = iEvent.id().event();

  // Fast SimTrack lookup: trackId -> ORIGIN momentum direction (global).
  // The SimTrack momentum is at production; we use it as the origin-momentum
  // approximation for the parent angle at the layer (v0; documented).
  std::map<unsigned int, math::XYZTLorentzVectorD> simMomById;
  for (const auto& st : *simTracks) {
    simMomById[st.trackId()] = st.momentum();
  }

  // The B-field magnitude (global, at origin) sets the helix curvature scale.
  // rInv from the TTTrack already encodes 1/R in cm^-1; we use it directly.
  const double bz0 = field.inTesla(GlobalPoint(0, 0, 0)).z();

  // Enumerate TBPX (Inner Tracker barrel pixel) modules once per event is
  // wasteful; but the geometry is const across the job. Build a per-layer list
  // of (mean radius, detids) lazily on the first event only.
  static std::vector<std::vector<const GeomDetUnit*>> layerModules;  // [layer-1] -> units
  static std::vector<double> layerRadius;                            // [layer-1] mean r
  if (layerModules.empty()) {
    layerModules.resize(nLayers_);
    layerRadius.assign(nLayers_, 0.);
    std::vector<int> nPerLayer(nLayers_, 0);
    int nUnits = 0, nPixBarrel = 0;
    for (const auto* det : geom.detUnits()) {
      ++nUnits;
      const DetId id = det->geographicalId();
      if (id.det() != DetId::Tracker)
        continue;
      if (id.subdetId() != static_cast<int>(PixelSubdetector::PixelBarrel))
        continue;
      ++nPixBarrel;
      const unsigned int lay = topo.pxbLayer(id);
      if (lay < 1 || static_cast<int>(lay) > nLayers_)
        continue;
      layerModules[lay - 1].push_back(det);
      layerRadius[lay - 1] += det->position().perp();
      nPerLayer[lay - 1] += 1;
    }
    std::cout << "[SmartPixelsPayloadAnalyzer] geom.detUnits()=" << nUnits << " pixelBarrel=" << nPixBarrel
              << std::endl;
    for (int l = 0; l < nLayers_; ++l) {
      if (nPerLayer[l] > 0)
        layerRadius[l] /= nPerLayer[l];
      std::cout << "[SmartPixelsPayloadAnalyzer] TBPX layer " << (l + 1) << ": " << nPerLayer[l]
                << " modules, mean r = " << layerRadius[l] << " cm" << std::endl;
    }
  }

  int nCrossingsThisEvent = 0;
  std::cout << "[SmartPixelsPayloadAnalyzer] event " << b_event_ << " nTracks=" << tracks->size() << std::endl;

  for (size_t itrk = 0; itrk < tracks->size(); ++itrk) {
    const auto& trk = (*tracks)[itrk];
    edm::Ptr<TTTrack<Ref_Phase2TrackerDigi_>> trkPtr(tracks, itrk);

    b_trk_pt_ = trk.momentum().perp();
    b_trk_eta_ = trk.momentum().eta();
    b_trk_phi_ = trk.momentum().phi();
    b_trk_charge_ = (trk.rInv() > 0) ? 1 : -1;
    b_trk_chi2_ = trk.chi2();
    b_trk_nstub_ = static_cast<int>(trk.getStubRefs().size());

    // Truth match (same association the producer uses).
    std::set<unsigned int> matchedSimIds;
    EncodedEventId matchedEvtId;
    b_trk_tpMatched_ = 0;
    b_tp_pt_ = -999;
    b_tp_eta_ = -999;
    edm::Ptr<TrackingParticle> tp = assoc->findTrackingParticlePtr(trkPtr);
    if (tp.isNonnull()) {
      b_trk_tpMatched_ = 1;
      b_tp_pt_ = tp->p4().pt();
      b_tp_eta_ = tp->p4().eta();
      matchedEvtId = tp->eventId();
      for (const auto& g4 : tp->g4Tracks())
        matchedSimIds.insert(g4.trackId());
    }

    // Helix parameters (r-phi circle) from the TTTrack.
    const double rInv = trk.rInv();            // signed 1/R [cm^-1]
    const double phi0 = trk.momentum().phi();  // momentum phi at POCA
    const double x0 = trk.POCA().x();
    const double y0 = trk.POCA().y();
    const double z0 = trk.z0();
    const double tanL = trk.tanL();
    const double cosPhi0 = std::cos(phi0);
    const double sinPhi0 = std::sin(phi0);
    (void)bz0;  // curvature taken directly from rInv; bz0 kept for provenance

    for (int l = 0; l < nLayers_; ++l) {
      const double Rl = layerRadius[l];
      if (Rl <= 0.)
        continue;

      // Analytic helix -> cylinder(R=Rl) intersection in the transverse plane.
      // Arc-length parametrization of a circle of curvature rInv starting at
      // (x0,y0) with direction (cosPhi0,sinPhi0). Solve for the turning angle
      // psi at which the transverse radius equals Rl. Small-|rInv| safe.
      double sTransverse;  // transverse path length to the crossing
      if (std::abs(rInv) < 1e-6) {
        // straight line: |(x0,y0) + s*dir| = Rl
        const double b = x0 * cosPhi0 + y0 * sinPhi0;
        const double c = x0 * x0 + y0 * y0 - Rl * Rl;
        const double disc = b * b - c;
        if (disc < 0)
          continue;
        sTransverse = -b + std::sqrt(disc);
        if (sTransverse < 0)
          continue;
      } else {
        // Circle center is offset by 1/rInv perpendicular to the momentum.
        const double R = 1.0 / rInv;  // signed
        const double cx = x0 - R * sinPhi0;
        const double cy = y0 + R * cosPhi0;
        const double dcenter = std::hypot(cx, cy);
        // crossing exists only if the circle reaches radius Rl
        const double absR = std::abs(R);
        if (Rl > dcenter + absR || Rl < std::abs(dcenter - absR))
          continue;
        // Law of cosines for the turning angle from POCA to the crossing.
        // NB: must use |R| — a signed R flips cosArg for negative charge and
        // sends the crossing to the far side of the circle (psi -> pi - psi).
        const double cosArg = (absR * absR + dcenter * dcenter - Rl * Rl) / (2.0 * absR * dcenter);
        if (cosArg < -1.0 || cosArg > 1.0)
          continue;
        // psi = turning angle; take the forward (smallest positive) solution.
        const double phiCenterToStart = std::atan2(y0 - cy, x0 - cx);
        const double delta = std::acos(std::max(-1.0, std::min(1.0, cosArg)));
        // pick the turning direction consistent with the charge sign
        const double dpsi = (rInv > 0) ? delta : -delta;
        (void)phiCenterToStart;
        sTransverse = std::abs(dpsi) * absR;
      }

      // Global crossing point via helix stepping in transverse arc length.
      double gx, gy;
      if (std::abs(rInv) < 1e-6) {
        gx = x0 + sTransverse * cosPhi0;
        gy = y0 + sTransverse * sinPhi0;
      } else {
        const double psi = rInv * sTransverse;  // turning angle
        gx = x0 + (std::sin(phi0 + psi) - sinPhi0) / rInv;
        gy = y0 - (std::cos(phi0 + psi) - cosPhi0) / rInv;
      }
      const double gz = z0 + tanL * sTransverse;
      const GlobalPoint gp(gx, gy, gz);

      // Global momentum direction at the crossing (rotate phi0 by psi).
      const double psi = rInv * sTransverse;
      const double pt = b_trk_pt_;
      const GlobalVector gmom(pt * std::cos(phi0 + psi), pt * std::sin(phi0 + psi), pt * tanL);

      // Find the module in this layer whose surface the crossing lands on.
      // Modules are tilted/staggered, so the mean-radius crossing point has a
      // nonzero local-z on each module; we (1) pick the geometrically closest
      // module center, then (2) require the crossing to be within that module's
      // TRANSVERSE (x,y) bounds only (ignore thickness/local-z), which is the
      // correct containment test for a thin sensor. This tolerates the O(mm)
      // radial spread of a tilted layer without a full plane re-propagation.
      const GeomDetUnit* best = nullptr;
      double bestDist2 = 1e18;
      LocalPoint bestLocal;
      for (const auto* det : layerModules[l]) {
        const LocalPoint lp = det->toLocal(gp);
        if (!det->surface().bounds().inside(LocalPoint(lp.x(), lp.y(), 0.)))
          continue;
        const auto dc = det->position() - gp;
        const double d2 = dc.mag2();
        if (d2 < bestDist2) {
          bestDist2 = d2;
          best = det;
          bestLocal = lp;
        }
      }
      if (!best)
        continue;

      const PixelGeomDetUnit* pixDet = dynamic_cast<const PixelGeomDetUnit*>(best);
      if (!pixDet)
        continue;
      const PixelTopology& pixTopo = pixDet->specificTopology();
      const DetId detId = best->geographicalId();

      // Track local direction -> PixelAV cotAlpha/cotBeta.
      const LocalVector lmom = best->toLocal(gmom);
      const double lpz = (std::abs(lmom.z()) > 1e-9) ? lmom.z() : 1e-9;
      b_trk_cotAlpha_ = lmom.x() / lpz;
      b_trk_cotBeta_ = lmom.y() / lpz;

      // Local B-field at the module position.
      const GlobalVector gB = field.inTesla(best->position());
      const LocalVector lB = best->toLocal(gB);
      b_bx_ = lB.x();
      b_by_ = lB.y();
      b_bz_ = lB.z();

      b_layer_ = l + 1;
      b_detid_ = detId.rawId();
      b_exp_localx_ = bestLocal.x();
      b_exp_localy_ = bestLocal.y();

      // ---- collect in-window digis on this module ----
      b_nWinDigi_ = 0;
      const auto digiSet = digis->find(detId.rawId());
      if (digiSet != digis->end()) {
        // Build a channel->simlink map for this module for fast classification.
        std::map<unsigned int, const PixelDigiSimLink*> linkByChannel;
        const auto linkSet = simlinks->find(detId.rawId());
        if (linkSet != simlinks->end()) {
          for (const auto& lk : *linkSet) {
            // keep the highest-fraction link per channel
            auto it = linkByChannel.find(lk.channel());
            if (it == linkByChannel.end() || it->second->fraction() < lk.fraction())
              linkByChannel[lk.channel()] = &lk;
          }
        }

        for (const auto& digi : *digiSet) {
          const LocalPoint dlp =
              pixTopo.localPosition(MeasurementPoint(digi.row() + 0.5, digi.column() + 0.5));
          const double dx = dlp.x() - bestLocal.x();
          const double dy = dlp.y() - bestLocal.y();
          if (std::abs(dx) > windowRPhi_ || std::abs(dy) > windowZ_)
            continue;
          if (b_nWinDigi_ >= kMaxWinDigis)
            break;

          int cls = kClassNoise;
          float parCotA = -999.f, parCotB = -999.f;
          auto lit = linkByChannel.find(digi.channel());
          if (lit != linkByChannel.end()) {
            const PixelDigiSimLink* lk = lit->second;
            const bool sameTP =
                b_trk_tpMatched_ && (lk->eventId() == matchedEvtId) && matchedSimIds.count(lk->SimTrackId());
            cls = sameTP ? kClassSameTP : kClassOtherTP;
            // parent local angle from the linked SimTrack origin momentum
            auto mit = simMomById.find(lk->SimTrackId());
            if (mit != simMomById.end()) {
              const auto& pm = mit->second;
              const GlobalVector pgv(pm.px(), pm.py(), pm.pz());
              const LocalVector plv = best->toLocal(pgv);
              const double ppz = (std::abs(plv.z()) > 1e-9) ? plv.z() : 1e-9;
              parCotA = plv.x() / ppz;
              parCotB = plv.y() / ppz;
            }
          }

          const int idx = b_nWinDigi_;
          b_digi_class_[idx] = cls;
          b_digi_localx_[idx] = dlp.x();
          b_digi_localy_[idx] = dlp.y();
          b_digi_resx_[idx] = dx;
          b_digi_resy_[idx] = dy;
          b_digi_adc_[idx] = digi.adc();
          b_digi_parCotAlpha_[idx] = parCotA;
          b_digi_parCotBeta_[idx] = parCotB;
          ++b_nWinDigi_;
        }
      }

      if (debug_) {
        edm::LogInfo("SmartPixelsPayloadAnalyzer")
            << "trk " << itrk << " layer " << (l + 1) << " detid " << detId.rawId() << " matched "
            << b_trk_tpMatched_ << " nWinDigi " << b_nWinDigi_ << " cotAlpha " << b_trk_cotAlpha_ << " cotBeta "
            << b_trk_cotBeta_;
      }

      tree_->Fill();
      ++nCrossingsThisEvent;
    }  // layers
  }  // tracks
  std::cout << "[SmartPixelsPayloadAnalyzer] event " << b_event_ << " crossings filled=" << nCrossingsThisEvent
            << std::endl;
}

void SmartPixelsPayloadAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("l1TracksInputTag", edm::InputTag("l1tTTTracksFromTrackletEmulation", "Level1TTTracks"));
  desc.add<edm::InputTag>("mcTruthTrackInputTag", edm::InputTag("TTTrackAssociatorFromPixelDigis", "Level1TTTracks"));
  desc.add<edm::InputTag>("trackingParticleInputTag", edm::InputTag("mix", "MergedTrackTruth"));
  desc.add<edm::InputTag>("pixelDigiInputTag", edm::InputTag("simSiPixelDigis", "Pixel"));
  desc.add<edm::InputTag>("pixelDigiSimLinkInputTag", edm::InputTag("simSiPixelDigis", "Pixel"));
  desc.add<edm::InputTag>("simTrackInputTag", edm::InputTag("g4SimHits"));
  desc.add<edm::InputTag>("simVertexInputTag", edm::InputTag("g4SimHits"));
  // window half-widths default to DIGIREFIT_DEFAULTS (windowRPhi=0.05, windowZ=0.10 cm)
  desc.add<double>("windowRPhi", 0.05);
  desc.add<double>("windowZ", 0.10);
  desc.add<int>("nLayers", 4);
  desc.add<bool>("debug", false);
  descriptions.add("smartPixelsPayloadAnalyzer", desc);
}

DEFINE_FWK_MODULE(SmartPixelsPayloadAnalyzer);
