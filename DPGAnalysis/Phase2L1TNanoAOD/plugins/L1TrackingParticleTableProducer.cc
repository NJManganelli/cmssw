// Flat table of TrackingParticles: the efficiency denominator for L1 track studies.
//
// One row per TrackingParticle passing (charge != 0, pt >= minPt). `idx` is the
// TrackingParticle collection key, the same index as L1TOTStub_tpIdx and
// L1TSmartPixelsCluster_tpIdx, so TP identity joins exactly across tables and
// across jobs reading the same input events.
//
// phi0/d0/z0 are at the POCA to the beamline with the same curvature-aware propagation
// and sign convention as L1TrackTruthTableProducer's tp_phi0/tp_d0/tp_z0, so a TP row and
// a matched track's tp_* columns agree. `phi` is at the production vertex; build a helix
// from phi0, never from phi.
//
// GENERATOR LINK (doGen = true; needs the nano gen chain, i.e. finalGenParticles):
//   genPartIdx >= 0   the TP's own generator particle is GenPart row genPartIdx
//   genPartIdx == -1  the TP has no generator particle (pileup, or produced in Geant4)
//   genPartIdx == -2  it has one, but neither it nor any ancestor survives the pruning
//   genPartIdx <= -3  its own particle was pruned; the nearest KEPT ancestor (first-mother
//                     chain) is GenPart row (-3 - genPartIdx)
// The GenPart row index is the key in the pruned finalGenParticles collection, which holds
// because genParticleTable applies no cut. genFromB / genFromC walk the FULL (unpruned)
// generator record, so pruning never hides an ancestor: true when any ancestor (not the
// particle itself) is a b / c hadron. A b-hadron -> c-hadron -> pion chain sets both.
// Geant4 secondaries (genPartIdx == -1) carry no ancestry.
//
// SIZE: at PU200 with minPt = 1 GeV this is ~15k rows/event and roughly half of
// an L1 track nano, hence opt-in (addPh2L1TrackingParticles), not part of the
// default truth task.
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"

#include "CLHEP/Units/PhysicalConstants.h"

#include <cmath>
#include <cstdlib>

class L1TrackingParticleTableProducer : public edm::stream::EDProducer<> {
public:
  using GenMap = edm::Association<reco::GenParticleCollection>;

  explicit L1TrackingParticleTableProducer(const edm::ParameterSet& cfg)
      : tpToken_(consumes<TrackingParticleCollection>(cfg.getParameter<edm::InputTag>("trackingParticles"))),
        bFieldToken_(esConsumes<MagneticField, IdealMagneticFieldRecord>()),
        name_(cfg.getParameter<std::string>("name")),
        minPt_(cfg.getParameter<double>("minPt")),
        doGen_(cfg.getParameter<bool>("doGen")) {
    if (doGen_)
      genMapToken_ = consumes<GenMap>(cfg.getParameter<edm::InputTag>("genParticlePruningMap"));
    produces<nanoaod::FlatTable>();
  }

  void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override {
    const auto& tps = iEvent.get(tpToken_);
    const GenMap* genMap = doGen_ ? &iEvent.get(genMapToken_) : nullptr;
    const float b_field = iSetup.getData(bFieldToken_).inTesla(GlobalPoint(0.f, 0.f, 0.f)).z();
    const float convertRtoPt = CLHEP::c_light / 1.0E5 * b_field;

    std::vector<int32_t> idx, pdgId, bx, evt;
    std::vector<float> pt, eta, phi, charge, vx, vy, vz, phi0, d0, z0;
    std::vector<int32_t> genPartIdx;
    std::vector<bool> genFromB, genFromC;
    for (size_t i = 0; i < tps.size(); ++i) {
      const TrackingParticle& tp = tps[i];
      if (tp.charge() == 0 || tp.pt() < minPt_)
        continue;
      const float tpt = tp.pt(), tphi = tp.phi(), q = static_cast<float>(tp.charge());
      const float tvx = tp.vx(), tvy = tp.vy(), tvz = tp.vz();
      // POCA propagation, identical to L1TrackTruthTableProducer.
      const float r2_inv = q * convertRtoPt / tpt / 2.0f;
      const float x0p = -tvx - (1.f / (2.f * r2_inv) * std::sin(tphi));
      const float y0p = -tvy + (1.f / (2.f * r2_inv) * std::cos(tphi));
      const float rp = std::sqrt(x0p * x0p + y0p * y0p);
      const float tphi0 = std::atan2(-r2_inv * x0p, r2_inv * y0p);
      const float delphi = reco::deltaPhi(tphi, tphi0);
      idx.push_back(static_cast<int32_t>(i));
      pdgId.push_back(tp.pdgId());
      bx.push_back(tp.eventId().bunchCrossing());
      evt.push_back(tp.eventId().event());
      pt.push_back(tpt);
      eta.push_back(tp.eta());
      phi.push_back(tphi);
      charge.push_back(q);
      vx.push_back(tvx);
      vy.push_back(tvy);
      vz.push_back(tvz);
      phi0.push_back(tphi0);
      d0.push_back(q * rp - (1.f / (2.f * r2_inv)));
      z0.push_back(tvz + std::sinh(tp.eta()) * delphi / (2.0f * r2_inv));
      if (genMap) {
        int32_t link = -1;
        bool fromB = false, fromC = false;
        if (!tp.genParticles().empty()) {
          const reco::GenParticleRef g = tp.genParticles()[0];
          if (!genMap->contains(g.id()))
            throw cms::Exception("Configuration")
                << "L1TrackingParticleTableProducer: the TrackingParticle generator links point into a "
                   "GenParticle product the pruning map was not built from; set finalGenParticles.src to "
                   "the collection the TrackingParticles reference (useGenParticlesFromFile).";
          const reco::GenParticleRef kept = (*genMap)[g];
          if (kept.isNonnull()) {
            link = static_cast<int32_t>(kept.key());
          } else {
            link = -2;
            reco::GenParticleRef a = g;
            for (int depth = 0; depth < kMaxDepth && a->numberOfMothers() > 0; ++depth) {
              a = a->motherRef(0);
              const reco::GenParticleRef k = (*genMap)[a];
              if (k.isNonnull()) {
                link = -3 - static_cast<int32_t>(k.key());
                break;
              }
            }
          }
          reco::GenParticleRef a = g;
          for (int depth = 0; depth < kMaxDepth && a->numberOfMothers() > 0; ++depth) {
            a = a->motherRef(0);
            fromB |= isHadronOf(a->pdgId(), 5);
            fromC |= isHadronOf(a->pdgId(), 4);
          }
        }
        genPartIdx.push_back(link);
        genFromB.push_back(fromB);
        genFromC.push_back(fromC);
      }
    }
    auto t = std::make_unique<nanoaod::FlatTable>(idx.size(), name_, false, false);
    t->addColumn<int32_t>("idx", idx, "TrackingParticle collection index (= L1TOTStub_tpIdx)");
    t->addColumn<int32_t>("pdgId", pdgId, "pdgId");
    t->addColumn<int32_t>("bx", bx, "eventId().bunchCrossing()");
    t->addColumn<int32_t>("evt", evt, "eventId().event() (0 with bx 0 = hard interaction)");
    t->addColumn<float>("pt", pt, "pt (GeV)");
    t->addColumn<float>("eta", eta, "eta");
    t->addColumn<float>("phi", phi, "phi at production (rad)");
    t->addColumn<float>("charge", charge, "charge");
    t->addColumn<float>("vx", vx, "production vertex x (cm)");
    t->addColumn<float>("vy", vy, "production vertex y (cm)");
    t->addColumn<float>("vz", vz, "production vertex z (cm)");
    t->addColumn<float>("phi0", phi0, "phi (rad) at the POCA to the beamline, the helix phi0 that goes with d0/z0");
    t->addColumn<float>("d0", d0, "d0 (cm) at POCA to the beamline, L1TTrack_d0 sign convention");
    t->addColumn<float>("z0", z0, "z0 (cm) at POCA to the beamline");
    if (doGen_) {
      t->addColumn<int32_t>("genPartIdx",
                            genPartIdx,
                            "GenPart row of this TP's generator particle (>= 0); -1 no generator particle "
                            "(pileup or Geant4); -2 generator particle with no kept ancestor; <= -3 the particle "
                            "was pruned and its nearest kept ancestor is GenPart row (-3 - genPartIdx)");
      t->addColumn<bool>("genFromB", genFromB, "a b hadron is among the generator ancestors (full record)");
      t->addColumn<bool>("genFromC", genFromC, "a c hadron is among the generator ancestors (full record)");
    }
    iEvent.put(std::move(t));
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("trackingParticles", edm::InputTag("mix", "MergedTrackTruth"));
    desc.add<std::string>("name", "L1TTP");
    desc.add<double>("minPt", 1.0)->setComment("keep TPs with pt >= minPt [GeV]");
    desc.add<bool>("doGen", false)->setComment("add genPartIdx / genFromB / genFromC (needs the gen chain)");
    desc.add<edm::InputTag>("genParticlePruningMap", edm::InputTag("finalGenParticles"))
        ->setComment("GenParticlePruner behind the GenPart table; its full->pruned Association is used");
    descriptions.addWithDefaultLabel(desc);
  }

private:
  static constexpr int kMaxDepth = 256;  // guard against malformed mother chains

  // b (q = 5) or c (q = 4) hadron: meson |id| = ..0 q x x, baryon |id| = ..q x x x.
  // Quarks (|id| < 10) and diquarks (tens digit 0) are not hadrons; quarkonia count.
  static bool isHadronOf(int pdgId, int q) {
    const int a = std::abs(pdgId) % 10000;
    if (a < 100 || (a / 10) % 10 == 0)
      return false;
    return (a / 1000) % 10 == q || ((a / 1000) % 10 == 0 && (a / 100) % 10 == q);
  }

  const edm::EDGetTokenT<TrackingParticleCollection> tpToken_;
  const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> bFieldToken_;
  const std::string name_;
  const double minPt_;
  const bool doGen_;
  edm::EDGetTokenT<GenMap> genMapToken_;
};

DEFINE_FWK_MODULE(L1TrackingParticleTableProducer);
