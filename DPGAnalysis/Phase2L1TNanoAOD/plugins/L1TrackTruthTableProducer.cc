// Emits an extension table on a L1 track table with:
//   - nStubs (always available; needed by the track-quality GBDT)
//   - MC truth columns from the TTTrackAssociationMap when present:
//     genuine / looselyGenuine / combinatoric / unknown flags and matched
//     TrackingParticle pt/eta/pdgId (defaults when the associator did not
//     run, e.g. data or central workflows without the truth sequence).
//   - the matched TrackingParticle's FULL helix parameters propagated to its
//     POCA to the beamline (tp_pt/tp_eta/tp_phi/tp_tanL/tp_charge/tp_d0/tp_z0)
//     plus its production vertex (tp_vx/tp_vy/tp_vz). These enable
//     resolution-vs-truth for DISPLACED / b-decay tracks whose TRUE d0 is
//     nonzero and is NOT captured by GenVtx (the PV). Sentinel -999 when the
//     track is not genuinely matched to a TP.
//
// TP d0/z0 use the exact curvature-aware POCA propagation from
// L1Trigger/TrackFindingTracklet/test/L1TrackNtupleMaker.cc (which itself
// propagates back to the IP rather than trusting TrackingParticle::d0()/z0()).
// The sign convention matches TTTrack::d0(): TTTrack's ctor sets
// thePOCA_(d0*sin(phi0), -d0*cos(phi0), z0), i.e. d0 = x0*sin(phi)-y0*cos(phi),
// so a positive tp_d0 is directly comparable to L1TTrack_d0 for the resolution.
//
// Index/row alignment relies on the track table using the same source
// collection with an empty cut and native ordering.
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "SimDataFormats/Associations/interface/TTTrackAssociationMap.h"

#include "CLHEP/Units/PhysicalConstants.h"

#include <cmath>

class L1TrackTruthTableProducer : public edm::stream::EDProducer<> {
public:
  using L1Track = TTTrack<Ref_Phase2TrackerDigi_>;
  using L1TrackCollection = std::vector<L1Track>;
  using AssocMap = TTTrackAssociationMap<Ref_Phase2TrackerDigi_>;

  explicit L1TrackTruthTableProducer(const edm::ParameterSet& cfg)
      : tracksToken_(consumes<L1TrackCollection>(cfg.getParameter<edm::InputTag>("tracks"))),
        assocToken_(consumes<AssocMap>(cfg.getParameter<edm::InputTag>("truthAssociation"))),
        bFieldToken_(esConsumes<MagneticField, IdealMagneticFieldRecord>()),
        trackTableName_(cfg.getParameter<std::string>("trackTableName")) {
    produces<nanoaod::FlatTable>();
  }

  void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override {
    auto tracksHandle = iEvent.getHandle(tracksToken_);
    auto assocHandle = iEvent.getHandle(assocToken_);

    // Conversion factor curvature-radius <-> pT (same as L1TrackNtupleMaker).
    const auto& bField = iSetup.getData(bFieldToken_);
    const float b_field = bField.inTesla(GlobalPoint(0.f, 0.f, 0.f)).z();
    const float cLight = CLHEP::c_light / 1.0E5;
    const float convertRtoPt = cLight * b_field;

    const size_t nTracks = tracksHandle.isValid() ? tracksHandle->size() : 0;
    std::vector<uint8_t> nStubs(nTracks, 0);
    std::vector<bool> genuine(nTracks, false), loose(nTracks, false), combinatoric(nTracks, false),
        unknown(nTracks, false);
    std::vector<float> tpPt(nTracks, -1.f), tpEta(nTracks, 0.f);
    std::vector<int32_t> tpPdgId(nTracks, 0);
    std::vector<bool> tpFromHardInteraction(nTracks, false);
    // Matched-TP full helix params @ POCA + production vertex (-999 sentinel for
    // tracks without a genuine TP match; consumers gate on genuine / tp_pt>0).
    constexpr float kSentinel = -999.f;
    std::vector<float> tpPhi(nTracks, kSentinel), tpTanL(nTracks, kSentinel), tpCharge(nTracks, kSentinel);
    std::vector<float> tpD0(nTracks, kSentinel), tpZ0(nTracks, kSentinel);
    std::vector<float> tpVx(nTracks, kSentinel), tpVy(nTracks, kSentinel), tpVz(nTracks, kSentinel);

    for (size_t i = 0; i < nTracks; ++i) {
      const auto& track = (*tracksHandle)[i];
      nStubs[i] = track.getStubRefs().size();
      if (assocHandle.isValid()) {
        edm::Ptr<L1Track> ptr(tracksHandle, i);
        genuine[i] = assocHandle->isGenuine(ptr);
        loose[i] = assocHandle->isLooselyGenuine(ptr);
        combinatoric[i] = assocHandle->isCombinatoric(ptr);
        unknown[i] = assocHandle->isUnknown(ptr);
        const auto& tp = assocHandle->findTrackingParticlePtr(ptr);
        if (tp.isNonnull()) {
          const float pt = tp->pt();
          const float phi = tp->phi();
          const float charge = static_cast<float>(tp->charge());
          const float vx = tp->vx();
          const float vy = tp->vy();
          const float vz = tp->vz();
          tpPt[i] = pt;
          tpEta[i] = tp->eta();
          tpPdgId[i] = tp->pdgId();
          tpFromHardInteraction[i] = (tp->eventId().event() == 0) && (tp->eventId().bunchCrossing() == 0);
          tpPhi[i] = phi;
          tpTanL[i] = std::sinh(tp->eta());  // tanL = pz/pt = sinh(eta)
          tpCharge[i] = charge;
          tpVx[i] = vx;
          tpVy[i] = vy;
          tpVz[i] = vz;
          // Curvature-aware POCA propagation to the beamline (L1TrackNtupleMaker).
          // Guard against pt==0 / charge==0 (would give r2_inv==0); such TPs keep
          // the -999 sentinel for d0/z0 only (params above are still valid).
          if (pt > 0.f && charge != 0.f) {
            const float delx = -vx;
            const float dely = -vy;
            const float r2_inv = charge * convertRtoPt / pt / 2.0f;
            const float x0p = delx - (1.f / (2.f * r2_inv) * std::sin(phi));
            const float y0p = dely + (1.f / (2.f * r2_inv) * std::cos(phi));
            const float rp = std::sqrt(x0p * x0p + y0p * y0p);
            tpD0[i] = charge * rp - (1.f / (2.f * r2_inv));
            const float delphi = reco::deltaPhi(phi, std::atan2(-r2_inv * x0p, r2_inv * y0p));
            tpZ0[i] = vz + std::sinh(tp->eta()) * delphi / (2.0f * r2_inv);
          }
        }
      }
    }

    auto table = std::make_unique<nanoaod::FlatTable>(nTracks, trackTableName_, false, true);
    table->addColumn<uint8_t>("nStubs", nStubs, "number of stubs on the track");
    table->addColumn<bool>("genuine", genuine, "all stubs from the same TrackingParticle (TTTrackAssociationMap)");
    table->addColumn<bool>("looselyGenuine", loose, "at most one stub unmatched (TTTrackAssociationMap)");
    table->addColumn<bool>("combinatoric", combinatoric, "combinatoric fake (TTTrackAssociationMap)");
    table->addColumn<bool>("unknown", unknown, "no truth association (TTTrackAssociationMap)");
    table->addColumn<float>("tpPt", tpPt, "matched TrackingParticle pt (-1 if none)");
    table->addColumn<float>("tpEta", tpEta, "matched TrackingParticle eta");
    table->addColumn<int32_t>("tpPdgId", tpPdgId, "matched TrackingParticle pdgId");
    table->addColumn<bool>(
        "tpFromHardInteraction", tpFromHardInteraction, "matched TrackingParticle from the hard interaction");
    // Matched-TP full helix params for resolution-vs-truth (displaced/b-tracks).
    table->addColumn<float>("tp_pt", tpPt, "matched TP pt (=tpPt; -999 if unmatched via tp_d0 gate)");
    table->addColumn<float>("tp_eta", tpEta, "matched TP eta (=tpEta)");
    table->addColumn<float>("tp_phi", tpPhi, "matched TP phi at production (rad); -999 if unmatched");
    table->addColumn<float>("tp_tanL", tpTanL, "matched TP tanLambda = sinh(eta); -999 if unmatched");
    table->addColumn<float>("tp_charge", tpCharge, "matched TP charge; -999 if unmatched");
    table->addColumn<float>(
        "tp_d0",
        tpD0,
        "matched TP transverse impact parameter (cm) at POCA to the beamline, curvature-aware; "
        "sign matches L1TTrack_d0 (d0=x0*sin(phi)-y0*cos(phi)); nonzero for displaced/b-tracks; -999 if unmatched");
    table->addColumn<float>(
        "tp_z0", tpZ0, "matched TP longitudinal impact parameter z0 (cm) at POCA to the beamline; -999 if unmatched");
    table->addColumn<float>("tp_vx", tpVx, "matched TP production vertex x (cm); -999 if unmatched");
    table->addColumn<float>("tp_vy", tpVy, "matched TP production vertex y (cm); -999 if unmatched");
    table->addColumn<float>("tp_vz", tpVz, "matched TP production vertex z (cm); -999 if unmatched");
    iEvent.put(std::move(table));
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("tracks", edm::InputTag("l1tTTTracksFromTrackletEmulation", "Level1TTTracks"));
    desc.add<edm::InputTag>("truthAssociation", edm::InputTag("TTTrackAssociatorFromPixelDigis", "Level1TTTracks"));
    desc.add<std::string>("trackTableName", "L1TTrack");
    descriptions.addWithDefaultLabel(desc);
  }

private:
  const edm::EDGetTokenT<L1TrackCollection> tracksToken_;
  const edm::EDGetTokenT<AssocMap> assocToken_;
  const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> bFieldToken_;
  const std::string trackTableName_;
};

DEFINE_FWK_MODULE(L1TrackTruthTableProducer);
