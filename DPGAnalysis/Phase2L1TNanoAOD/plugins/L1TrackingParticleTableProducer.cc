// Flat table of TrackingParticles: the efficiency denominator for L1 track studies.
//
// SELECTION (`selection`), one row per TrackingParticle that passes:
//   charged1GeV          charge != 0 and 3-momentum |p| > minP (not pT, so forward particles
//                        that reach the discs are kept)
//   linked               TPs linked, through the truth associations, to an object the nano
//                        writes: the TPs of the stub / track association maps and of the
//                        TrackingParticleRefVectors listed in the linked* parameters
//                        (findTrackingParticlePtr of every associated stub / track, i.e. the
//                        TP a tpIdx of that object points to)
//   charged1GeVOrLinked  the union of the two
//   all                  every TrackingParticle in the input, neutrals included
// The table doc states the selection a file was made with.
//
// TP LINKS. Alongside the table this puts an edm::ValueMap<int> over the TrackingParticle
// collection: the L1TTP row of each TP, or -1 when the selection pruned it. It is the ONLY
// place the selection lives; every table writing a tpIdx link (L1TTrack, L1TExtTrack,
// L1TOTStub, L1TSmartPixelsCluster, L1TSmartPixelsTrack) consumes it and encodes
//   tpIdx >= 0   row in L1TTP
//   tpIdx == -1  not truth-matched
//   tpIdx == -2  TrackingParticle product (or this map) unavailable
//   tpIdx <= -3  matched TP not selected into L1TTP (never with `all`); -3 - tpIdx is its
//                index in the unwritten TrackingParticle collection
// through DPGAnalysis/Phase2L1TNanoAOD/interface/L1TPRowLink.h. Being a ValueMap, a lookup
// with a Ref into a different TrackingParticle product throws instead of mis-indexing.
//
// phi0/d0/z0 are at the POCA to the beamline with the same curvature-aware propagation
// and sign convention as L1TrackTruthTableProducer's tp_phi0/tp_d0/tp_z0, so a TP row and
// a matched track's tp_* columns agree; a neutral TP is a straight line. `phi` is at the
// production vertex; build a helix from phi0, never from phi.
//
// GENERATOR LINK (doGen = true; needs the pruning chain from the GenParticles the TPs
// reference to finalGenParticles, genParticlePruningMaps, all run in the job):
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
// SIZE: at PU200 a 1 GeV pT cut gave ~15k rows/event, roughly half of an L1 track
// nano (|p| > 1 GeV keeps more), hence opt-in (addPh2L1TrackingParticles), not part
// of the default truth task.
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "SimDataFormats/Associations/interface/TTStubAssociationMap.h"
#include "SimDataFormats/Associations/interface/TTTrackAssociationMap.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"

#include "CLHEP/Units/PhysicalConstants.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

class L1TrackingParticleTableProducer : public edm::stream::EDProducer<> {
public:
  using GenMap = edm::Association<reco::GenParticleCollection>;
  using StubMap = TTStubAssociationMap<Ref_Phase2TrackerDigi_>;
  using TrackMap = TTTrackAssociationMap<Ref_Phase2TrackerDigi_>;

  explicit L1TrackingParticleTableProducer(const edm::ParameterSet& cfg)
      : tpToken_(consumes<TrackingParticleCollection>(cfg.getParameter<edm::InputTag>("trackingParticles"))),
        bFieldToken_(esConsumes<MagneticField, IdealMagneticFieldRecord>()),
        name_(cfg.getParameter<std::string>("name")),
        selection_(cfg.getParameter<std::string>("selection")),
        minP_(cfg.getParameter<double>("minP")),
        doGen_(cfg.getParameter<bool>("doGen")) {
    if (selection_ != "charged1GeV" && selection_ != "linked" && selection_ != "charged1GeVOrLinked" &&
        selection_ != "all")
      throw cms::Exception("Configuration")
          << "L1TrackingParticleTableProducer: selection '" << selection_
          << "' is not one of charged1GeV, linked, charged1GeVOrLinked, all";
    const bool linked = selection_ == "linked" || selection_ == "charged1GeVOrLinked";
    if (linked) {
      for (const auto& tag : cfg.getParameter<std::vector<edm::InputTag>>("linkedStubAssociations"))
        stubMapTokens_.push_back(consumes<StubMap>(tag));
      for (const auto& tag : cfg.getParameter<std::vector<edm::InputTag>>("linkedTrackAssociations"))
        trackMapTokens_.push_back(consumes<TrackMap>(tag));
      for (const auto& tag : cfg.getParameter<std::vector<edm::InputTag>>("linkedTrackingParticles"))
        tpRefTokens_.push_back(consumes<TrackingParticleRefVector>(tag));
      if (stubMapTokens_.empty() && trackMapTokens_.empty() && tpRefTokens_.empty())
        throw cms::Exception("Configuration")
            << "L1TrackingParticleTableProducer: selection '" << selection_ << "' with no linked objects";
    }
    if (doGen_)
      for (const auto& tag : cfg.getParameter<std::vector<edm::InputTag>>("genParticlePruningMaps"))
        genMapTokens_.push_back(consumes<GenMap>(tag));
    produces<nanoaod::FlatTable>();
    produces<edm::ValueMap<int>>();
  }

  void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override {
    const auto tpHandle = iEvent.getHandle(tpToken_);
    const auto& tps = *tpHandle;
    std::vector<const GenMap*> genMaps;
    for (const auto& token : genMapTokens_)
      genMaps.push_back(&iEvent.get(token));
    const float b_field = iSetup.getData(bFieldToken_).inTesla(GlobalPoint(0.f, 0.f, 0.f)).z();
    const float convertRtoPt = CLHEP::c_light / 1.0E5 * b_field;
    const std::vector<bool> linked = linkedKeys(iEvent, tpHandle);

    std::vector<int> keyToRow(tps.size(), -1);
    std::vector<int32_t> pdgId, bx, evt;
    std::vector<float> pt, eta, phi, charge, vx, vy, vz, phi0, d0, z0;
    std::vector<int32_t> genPartIdx;
    std::vector<bool> genFromB, genFromC;
    for (size_t i = 0; i < tps.size(); ++i) {
      const TrackingParticle& tp = tps[i];
      if (!selected(tp, linked[i]))
        continue;
      const float tpt = tp.pt(), tphi = tp.phi(), q = static_cast<float>(tp.charge());
      const float tvx = tp.vx(), tvy = tp.vy(), tvz = tp.vz();
      float tphi0 = tphi, td0, tz0;
      if (q != 0.f && tpt > 0.f) {
        // POCA propagation, identical to L1TrackTruthTableProducer.
        const float r2_inv = q * convertRtoPt / tpt / 2.0f;
        const float x0p = -tvx - (1.f / (2.f * r2_inv) * std::sin(tphi));
        const float y0p = -tvy + (1.f / (2.f * r2_inv) * std::cos(tphi));
        const float rp = std::sqrt(x0p * x0p + y0p * y0p);
        tphi0 = std::atan2(-r2_inv * x0p, r2_inv * y0p);
        td0 = q * rp - (1.f / (2.f * r2_inv));
        tz0 = tvz + std::sinh(tp.eta()) * reco::deltaPhi(tphi, tphi0) / (2.0f * r2_inv);
      } else {
        // straight line through the production vertex, same sign convention
        td0 = tvx * std::sin(tphi) - tvy * std::cos(tphi);
        tz0 = tvz - (tvx * std::cos(tphi) + tvy * std::sin(tphi)) * std::sinh(tp.eta());
      }
      keyToRow[i] = static_cast<int>(pdgId.size());
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
      d0.push_back(td0);
      z0.push_back(tz0);
      if (doGen_) {
        int32_t link = -1;
        bool fromB = false, fromC = false;
        if (!tp.genParticles().empty()) {
          const reco::GenParticleRef g = tp.genParticles()[0];
          const reco::GenParticleRef kept = pruned(genMaps, g);
          if (kept.isNonnull()) {
            link = static_cast<int32_t>(kept.key());
          } else {
            link = -2;
            reco::GenParticleRef a = g;
            for (int depth = 0; depth < kMaxDepth && a->numberOfMothers() > 0; ++depth) {
              a = a->motherRef(0);
              const reco::GenParticleRef k = pruned(genMaps, a);
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
    auto t = std::make_unique<nanoaod::FlatTable>(pdgId.size(), name_, false, false);
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
    t->setDoc("TrackingParticles: " + selectionDoc() +
              ". Other tables link here through tpIdx (row), with -1 unmatched, -2 truth unavailable, <= -3 "
              "matched but not selected into this table");
    iEvent.put(std::move(t));
    auto rowMap = std::make_unique<edm::ValueMap<int>>();
    edm::ValueMap<int>::Filler filler(*rowMap);
    filler.insert(tpHandle, keyToRow.begin(), keyToRow.end());
    filler.fill();
    iEvent.put(std::move(rowMap));
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("trackingParticles", edm::InputTag("mix", "MergedTrackTruth"));
    desc.add<std::string>("name", "L1TTP");
    desc.add<std::string>("selection", "charged1GeV")
        ->setComment("charged1GeV (charged, |p| > minP) | linked | charged1GeVOrLinked | all (neutrals included)");
    desc.add<double>("minP", 1.0)->setComment("charged1GeV: keep charged TPs with 3-momentum |p| > minP [GeV]");
    desc.add<std::vector<edm::InputTag>>("linkedStubAssociations", {})
        ->setComment("linked: TTStubAssociationMaps of the stubs the nano writes");
    desc.add<std::vector<edm::InputTag>>("linkedTrackAssociations", {})
        ->setComment("linked: TTTrackAssociationMaps of the tracks the nano writes");
    desc.add<std::vector<edm::InputTag>>("linkedTrackingParticles", {})
        ->setComment("linked: TrackingParticleRefVectors of other objects the nano writes (e.g. IT clusters)");
    desc.add<bool>("doGen", false)->setComment("add genPartIdx / genFromB / genFromC (needs the gen chain)");
    desc.add<std::vector<edm::InputTag>>("genParticlePruningMaps",
                                         {edm::InputTag("prunedGenParticlesWithStatusOne"),
                                          edm::InputTag("prunedGenParticles"),
                                          edm::InputTag("finalGenParticles")})
        ->setComment("GenParticlePruner chain, in order, from the GenParticles the TrackingParticles reference to the GenPart table");
    descriptions.addWithDefaultLabel(desc);
  }

private:
  static constexpr int kMaxDepth = 256;  // guard against malformed mother chains

  // The GenPart row of a full-record generator particle: follow the pruning chain from the
  // collection the TrackingParticles reference to finalGenParticles. Null if pruned.
  static reco::GenParticleRef pruned(const std::vector<const GenMap*>& maps, reco::GenParticleRef g) {
    for (const GenMap* map : maps) {
      if (!map->contains(g.id()))
        throw cms::Exception("Configuration")
            << "L1TrackingParticleTableProducer: genParticlePruningMaps is not a chain from the GenParticle "
               "product the TrackingParticles reference ("
            << g.id() << ") to the GenPart collection";
      g = (*map)[g];
      if (g.isNull())
        break;
    }
    return g;
  }

  bool selected(const TrackingParticle& tp, bool linked) const {
    if (selection_ == "all")
      return true;
    const bool charged1GeV = tp.charge() != 0 && tp.p() > minP_;
    if (selection_ == "charged1GeV")
      return charged1GeV;
    if (selection_ == "linked")
      return linked;
    return charged1GeV || linked;
  }

  std::string selectionDoc() const {
    char p[32];
    std::snprintf(p, sizeof(p), "%g", minP_);
    const std::string charged = std::string("charged with |p| > ") + p + " GeV";
    const std::string linked = "linked by the truth associations to a stub, track or cluster this nano writes";
    if (selection_ == "all")
      return "every TrackingParticle of the input, neutrals included";
    if (selection_ == "charged1GeV")
      return charged;
    if (selection_ == "linked")
      return linked;
    return charged + ", or " + linked;
  }

  // TP key -> linked to a written object. A link into another TrackingParticle product
  // would index the wrong particle, so it throws.
  std::vector<bool> linkedKeys(const edm::Event& iEvent, const edm::Handle<TrackingParticleCollection>& tps) const {
    std::vector<bool> linked(tps->size(), false);
    auto mark = [&](const edm::ProductID& id, size_t key, const char* what) {
      if (id != tps.id())
        throw cms::Exception("L1TrackingParticleTableProducer")
            << what << " links a TrackingParticle of product " << id << ", not of the table's " << tps.id();
      linked[key] = true;
    };
    for (const auto& token : stubMapTokens_)
      for (const auto& [stub, tp] : iEvent.get(token).getTTStubToTrackingParticleMap())
        if (tp.isNonnull())
          mark(tp.id(), tp.key(), "a stub association");
    for (const auto& token : trackMapTokens_)
      for (const auto& [track, tp] : iEvent.get(token).getTTTrackToTrackingParticleMap())
        if (tp.isNonnull())
          mark(tp.id(), tp.key(), "a track association");
    for (const auto& token : tpRefTokens_)
      for (const auto& tp : iEvent.get(token))
        mark(tp.id(), tp.key(), "a TrackingParticleRefVector");
    return linked;
  }

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
  const std::string selection_;
  const double minP_;
  const bool doGen_;
  std::vector<edm::EDGetTokenT<StubMap>> stubMapTokens_;
  std::vector<edm::EDGetTokenT<TrackMap>> trackMapTokens_;
  std::vector<edm::EDGetTokenT<TrackingParticleRefVector>> tpRefTokens_;
  std::vector<edm::EDGetTokenT<GenMap>> genMapTokens_;
};

DEFINE_FWK_MODULE(L1TrackingParticleTableProducer);
