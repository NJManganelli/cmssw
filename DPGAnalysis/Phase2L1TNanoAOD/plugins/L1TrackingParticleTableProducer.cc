// Flat table of TrackingParticles: the efficiency denominator for L1 track studies.
//
// One row per TrackingParticle passing (charge != 0, pt >= minPt). `idx` is the
// TrackingParticle collection key, the same index as L1TOTStub_tpIdx and
// L1TSmartPixelsCluster_tpIdx, so TP identity joins exactly across tables and
// across jobs reading the same input events.
//
// d0/z0 are at the POCA to the beamline with the same curvature-aware propagation
// and sign convention as L1TrackTruthTableProducer's tp_d0/tp_z0, so a TP row and
// a matched track's tp_* columns agree.
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

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"

#include "CLHEP/Units/PhysicalConstants.h"

#include <cmath>

class L1TrackingParticleTableProducer : public edm::stream::EDProducer<> {
public:
  explicit L1TrackingParticleTableProducer(const edm::ParameterSet& cfg)
      : tpToken_(consumes<TrackingParticleCollection>(cfg.getParameter<edm::InputTag>("trackingParticles"))),
        bFieldToken_(esConsumes<MagneticField, IdealMagneticFieldRecord>()),
        name_(cfg.getParameter<std::string>("name")),
        minPt_(cfg.getParameter<double>("minPt")) {
    produces<nanoaod::FlatTable>();
  }

  void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override {
    const auto& tps = iEvent.get(tpToken_);
    const float b_field = iSetup.getData(bFieldToken_).inTesla(GlobalPoint(0.f, 0.f, 0.f)).z();
    const float convertRtoPt = CLHEP::c_light / 1.0E5 * b_field;

    std::vector<int32_t> idx, pdgId, bx, evt;
    std::vector<float> pt, eta, phi, charge, vx, vy, vz, d0, z0;
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
      const float delphi = reco::deltaPhi(tphi, std::atan2(-r2_inv * x0p, r2_inv * y0p));
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
      d0.push_back(q * rp - (1.f / (2.f * r2_inv)));
      z0.push_back(tvz + std::sinh(tp.eta()) * delphi / (2.0f * r2_inv));
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
    t->addColumn<float>("d0", d0, "d0 (cm) at POCA to the beamline, L1TTrack_d0 sign convention");
    t->addColumn<float>("z0", z0, "z0 (cm) at POCA to the beamline");
    iEvent.put(std::move(t));
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("trackingParticles", edm::InputTag("mix", "MergedTrackTruth"));
    desc.add<std::string>("name", "L1TTP");
    desc.add<double>("minPt", 1.0)->setComment("keep TPs with pt >= minPt [GeV]");
    descriptions.addWithDefaultLabel(desc);
  }

private:
  const edm::EDGetTokenT<TrackingParticleCollection> tpToken_;
  const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> bFieldToken_;
  const std::string name_;
  const double minPt_;
};

DEFINE_FWK_MODULE(L1TrackingParticleTableProducer);
