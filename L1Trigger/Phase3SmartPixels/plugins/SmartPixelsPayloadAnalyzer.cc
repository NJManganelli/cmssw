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

// Shared helix propagation + TBPX module lookup (one impl, no drift with producer)
#include "L1Trigger/Phase3SmartPixels/interface/SmartPixelsHelixProjector.h"
#include "L1Trigger/Phase3SmartPixels/interface/SmartPixelsParentMap.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "FWCore/Utilities/interface/Exception.h"

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
  constexpr int kMaxWinDigis = 256;  // wide L3/L4 windows at PU can exceed 128
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
  const edm::EDGetTokenT<std::vector<PSimHit>> simHitToken_;  // true crossings (diagnostic)

  const edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> geomToken_;
  const edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> topoToken_;
  const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> fieldToken_;

  // Per-layer window half-widths [cm]. The rphi (local x) spread of the
  // extrapolation GROWS outward: the beamline-constrained 4-par fit is pinned
  // at the origin and at the OT stubs, and the MS-compensation bulge peaks
  // between them, right on L3/L4 (PSimHit-verified: q68 0.04/0.16/0.49/0.90 cm
  // at 2-5 GeV). z (local y) shrinks outward (no vertex constraint in z).
  const std::vector<double> windowRPhi_;  // local x half-widths, [layer-1]
  const std::vector<double> windowZ_;     // local y half-widths, [layer-1]
  const int nLayers_;        // number of TBPX layers to consider (1..4)
  const bool debug_;

  // Shared propagation/module-lookup (built lazily on the first event).
  smartpixels::HelixProjector projector_;

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
  int b_nSimHit_ = 0;                                 // matched-TP PSimHits on the assigned module
  float b_sim_localx_ = -999.f, b_sim_localy_ = -999.f;  // nearest such PSimHit (true crossing)
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
      simHitToken_(consumes<std::vector<PSimHit>>(
          iConfig.getParameter<edm::InputTag>("simHitInputTag"))),
      geomToken_(esConsumes<TrackerGeometry, TrackerDigiGeometryRecord>()),
      topoToken_(esConsumes<TrackerTopology, TrackerTopologyRcd>()),
      fieldToken_(esConsumes<MagneticField, IdealMagneticFieldRecord>()),
      windowRPhi_(iConfig.getParameter<std::vector<double>>("windowRPhi")),
      windowZ_(iConfig.getParameter<std::vector<double>>("windowZ")),
      nLayers_(iConfig.getParameter<int>("nLayers")),
      debug_(iConfig.getParameter<bool>("debug")) {
  usesResource(TFileService::kSharedResource);

  if (windowRPhi_.size() != 4 || windowZ_.size() != 4)
    throw cms::Exception("Configuration")
        << "windowRPhi/windowZ must each have exactly 4 per-layer entries (TBPX L1-L4), got "
        << windowRPhi_.size() << "/" << windowZ_.size() << ".";

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

  tree_->Branch("nSimHit", &b_nSimHit_, "nSimHit/I");
  tree_->Branch("sim_localx", &b_sim_localx_, "sim_localx/F");
  tree_->Branch("sim_localy", &b_sim_localy_, "sim_localy/F");
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
  edm::Handle<std::vector<PSimHit>> simHits;
  iEvent.getByToken(simHitToken_, simHits);
  std::map<unsigned int, std::vector<const PSimHit*>> simHitsByDet;
  for (const auto& sh : *simHits)
    simHitsByDet[sh.detUnitId()].push_back(&sh);

  b_event_ = iEvent.id().event();

  // Parent lookup keyed (eventId, trackId): TrackingParticles' embedded
  // g4Tracks cover signal AND pileup parents (the g4SimHits SimTrack container
  // is signal-only); the signal container remains a fallback for pruned TPs.
  edm::Handle<std::vector<TrackingParticle>> tps;
  iEvent.getByToken(tpToken_, tps);
  const smartpixels::ParentMomentumMap parentMom = smartpixels::buildParentMomentumMap(*tps, simTracks.product());

  // The B-field magnitude (global, at origin) sets the helix curvature scale.
  // rInv from the TTTrack already encodes 1/R in cm^-1; we use it directly.
  const double bz0 = field.inTesla(GlobalPoint(0, 0, 0)).z();

  // Build the shared per-layer TBPX module cache (idempotent across events).
  projector_.build(geom, topo, nLayers_);

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
    (void)bz0;  // curvature taken directly from rInv; bz0 kept for provenance

    smartpixels::HelixParams hp;
    hp.rInv = rInv;
    hp.phi0 = phi0;
    hp.x0 = x0;
    hp.y0 = y0;
    hp.z0 = z0;
    hp.tanL = tanL;
    hp.pt = b_trk_pt_;

    for (int l = 0; l < nLayers_; ++l) {
      const smartpixels::Crossing cx = projector_.crossLayer(hp, l + 1, field);
      if (!cx.valid)
        continue;

      const PixelGeomDetUnit* pixDet = cx.det;
      const PixelTopology& pixTopo = pixDet->specificTopology();
      const DetId detId(cx.detId);
      const LocalPoint bestLocal = cx.local;
      const GeomDet* best = pixDet;

      b_trk_cotAlpha_ = cx.cotAlpha;
      b_trk_cotBeta_ = cx.cotBeta;
      b_bx_ = cx.bLocalX;
      b_by_ = cx.bLocalY;
      b_bz_ = cx.bLocalZ;

      b_layer_ = l + 1;
      b_detid_ = cx.detId;
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
          if (std::abs(dx) > windowRPhi_[l] || std::abs(dy) > windowZ_[l])
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
            // parent local angle from the linked parent's origin momentum
            auto mit = parentMom.find({lk->eventId().rawId(), lk->SimTrackId()});
            if (mit != parentMom.end()) {
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

      // Diagnostic: the matched TP's TRUE crossing(s) on the assigned module
      // (signal-only PSimHits; nearest to the extrapolation if several).
      b_nSimHit_ = 0;
      b_sim_localx_ = -999.f;
      b_sim_localy_ = -999.f;
      const auto shIt = simHitsByDet.find(cx.detId);
      if (b_trk_tpMatched_ && shIt != simHitsByDet.end()) {
        double bestD2 = 1e18;
        for (const PSimHit* sh : shIt->second) {
          if (!matchedSimIds.count(sh->trackId()))
            continue;
          ++b_nSimHit_;
          const auto slp = sh->localPosition();
          const double dxs = slp.x() - bestLocal.x();
          const double dys = slp.y() - bestLocal.y();
          const double d2 = dxs * dxs + dys * dys;
          if (d2 < bestD2) {
            bestD2 = d2;
            b_sim_localx_ = slp.x();
            b_sim_localy_ = slp.y();
          }
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
  desc.add<edm::InputTag>("simHitInputTag", edm::InputTag("g4SimHits", "TrackerHitsPixelBarrelLowTof"));
  // Measurement windows must capture the FULL residual distribution
  // (truncation belongs to the refit). Per-layer, ~3x the PSimHit-measured
  // q68 at 2-5 GeV: rphi bulges outward (beamline-constrained fit + MS),
  // z shrinks outward.
  desc.add<std::vector<double>>("windowRPhi", std::vector<double>{0.15, 0.5, 1.5, 2.7});
  desc.add<std::vector<double>>("windowZ", std::vector<double>{0.7, 0.6, 0.5, 0.4});
  desc.add<int>("nLayers", 4);
  desc.add<bool>("debug", false);
  descriptions.add("smartPixelsPayloadAnalyzer", desc);
}

DEFINE_FWK_MODULE(SmartPixelsPayloadAnalyzer);
