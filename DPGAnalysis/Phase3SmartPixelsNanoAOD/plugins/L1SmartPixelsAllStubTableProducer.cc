// Nano producer for the COMPLETE Outer-Tracker stub collection -- every stub the
// L1 track finder was offered, not only the ones that ended up on a track.
//
// WHY THIS EXISTS, AND WHY THE EXISTING STUB TABLE COULD NOT ANSWER IT.
// L1SmartPixelsStubPosTableProducer walks track->getStubRefs(), so it is an
// ON-TRACK table by construction: measured 1004 rows/event at PU200, which is
// ~5.3 stubs on each of ~190 found tracks. That is the right table for
// visualising or refitting a track that already exists. It is the WRONG table
// for any combinatorics or cost question, because those are questions about the
// stubs the finder had to consider and reject -- exactly the rows getStubRefs()
// cannot contain. Comparing an Inner-Tracker seeding design against the OT
// tracklet baseline needs the OT INPUT multiplicity, and no table in this
// package carried it.
//
// It cannot be recovered offline either: the persisted collection is an
// edmNew::DetSetVector, which uproot refuses ("memberwise serialization of
// edmNew::dstvdetails::DetSetVectorTrans::Item"), so there is no path from the
// DIGI-RAW tier to a Python study without a producer like this one.
//
// WHICH STUBS. TTStubsFromPhase2TrackerDigis:StubAccepted -- the same label the
// tracklet chain consumes (see customizeSmartPixels_cff.py, "the PU RelVal
// persists the entire stub tier WITH pileup"). Consuming the same collection the
// finder used is deliberate and is the same discipline the cluster table follows:
// a table built from a second, differently-configured producer would silently
// stop describing the objects the algorithm actually saw.
//
// Note these are POST front-end-gate stubs. The Phase-2 pT modules only form a
// stub when the local bend is consistent with pT above roughly 2 GeV, so this
// collection is already the reduced one -- which is the correct baseline, since
// it is what the track finder receives.
//
// Output: one standalone FlatTable, one row per stub, NO trackIdx (there is no
// owning track for most of them -- that is the point). Columns are otherwise
// deliberately identical to the on-track table so the two are directly
// comparable row-for-row.
//
// SIZE. ~15k rows/event at PU200 is expected, against ~26.5k for the IT cluster
// table, so this roughly doubles the Clusters-tier payload. That is accepted
// on this tier specifically: a tier built to study IT clusters may as well carry
// the nearest OT equivalent, so that studies stop splicing two data tiers
// together to get one comparison.

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/L1TrackTrigger/interface/TTStub.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "SimDataFormats/Associations/interface/TTStubAssociationMap.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"

#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "Geometry/CommonTopologies/interface/GeomDet.h"
#include "Geometry/CommonTopologies/interface/Topology.h"
// Phase-2 OT reuses the strip-tracker subdetId numbering: barrel == TOB (5),
// endcap == TID (4). Same choice, and same reasoning, as the on-track producer.
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"

#include <cmath>
#include <cstdint>
#include <vector>

class L1SmartPixelsAllStubTableProducer : public edm::stream::EDProducer<> {
public:
  explicit L1SmartPixelsAllStubTableProducer(const edm::ParameterSet& cfg)
      : stubsToken_(consumes<TTStubDetSetVec>(cfg.getParameter<edm::InputTag>("stubs"))),
        geomToken_(esConsumes<TrackerGeometry, TrackerDigiGeometryRecord>()),
        topoToken_(esConsumes<TrackerTopology, TrackerTopologyRcd>()),
        tableName_(cfg.getParameter<std::string>("tableName")),
        barrelOnly_(cfg.getParameter<bool>("barrelOnly")),
        doTruth_(cfg.getParameter<bool>("doTruth")) {
    if (doTruth_)
      truthToken_ = consumes<TTStubAssociationMap<Ref_Phase2TrackerDigi_>>(
          cfg.getParameter<edm::InputTag>("stubTruth"));
    produces<nanoaod::FlatTable>();
  }

  void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override {
    // Handle, not iEvent.get: forming a TTStubRef via edmNew::makeRefTo is the
    // only way to query the association map, and that needs the handle.
    edm::Handle<TTStubDetSetVec> stubHandle;
    iEvent.getByToken(stubsToken_, stubHandle);
    const auto& stubs = *stubHandle;
    const TrackerGeometry& geom = iSetup.getData(geomToken_);
    const TrackerTopology& topo = iSetup.getData(topoToken_);

    edm::Handle<TTStubAssociationMap<Ref_Phase2TrackerDigi_>> truth;
    if (doTruth_)
      iEvent.getByToken(truthToken_, truth);

    // TRUTH-ONLY, and deliberately mirroring the IT cluster table's tp* block so
    // an IT-vs-OT comparison can use the same quantities on both sides: tpIdx is
    // the TrackingParticle collection index (.key()), exactly as the cluster
    // table derives it from dominantTp(). The three quality flags have no IT
    // analogue and come free from the association map.
    std::vector<int32_t> tpIdx;
    std::vector<float> tpPt;
    std::vector<bool> tpGenuine, tpCombinatoric, tpUnknown;

    std::vector<uint8_t> layer;
    std::vector<bool> isBarrel;
    std::vector<uint32_t> detId;
    std::vector<float> x, y, z, r, phi, bend;

    for (auto ds = stubs.begin(); ds != stubs.end(); ++ds) {
      for (auto st = ds->begin(); st != ds->end(); ++st) {
        // Global position from the INNER cluster, the L1TrackNtupleMaker pattern
        // and identical to the on-track table so the two agree stub-for-stub.
        const DetId did = geom.idToDet(st->clusterRef(0)->getDetId())->geographicalId();
        const MeasurementPoint coords = st->clusterRef(0)->findAverageLocalCoordinatesCentered();
        const GeomDet* det = geom.idToDet(did);
        const GlobalPoint pos = det->surface().toGlobal(det->topology().localPosition(coords));

        const bool barrel = (did.subdetId() == StripSubdetector::TOB);
        if (barrelOnly_ && !barrel)
          continue;
        const int lay = static_cast<int>(topo.layer(did));

        detId.push_back(did.rawId());
        isBarrel.push_back(barrel);
        layer.push_back(static_cast<uint8_t>(barrel ? lay : 10 + lay));
        x.push_back(pos.x());
        y.push_back(pos.y());
        z.push_back(pos.z());
        r.push_back(std::hypot(pos.x(), pos.y()));
        phi.push_back(std::atan2(pos.y(), pos.x()));
        bend.push_back(static_cast<float>(st->bendFE()));

        if (doTruth_) {
          const TTStubRef ref = edmNew::makeRefTo(stubHandle, st);
          const edm::Ptr<TrackingParticle> tp = truth->findTrackingParticlePtr(ref);
          tpIdx.push_back(tp.isNull() ? -1 : static_cast<int32_t>(tp.key()));
          tpPt.push_back(tp.isNull() ? -999.f : static_cast<float>(tp->pt()));
          tpGenuine.push_back(truth->isGenuine(ref));
          tpCombinatoric.push_back(truth->isCombinatoric(ref));
          tpUnknown.push_back(truth->isUnknown(ref));
        }
      }
    }

    auto tab = std::make_unique<nanoaod::FlatTable>(detId.size(), tableName_, false, false);
    tab->addColumn<uint8_t>("layer", layer,
                            "tracker layer/disk id: barrel L1..L6 = 1..6, endcap disk d = 10+d");
    tab->addColumn<bool>("isBarrel", isBarrel, "stub in the OT barrel (TOB) vs endcap");
    tab->addColumn<uint32_t>("detId", detId, "module rawId carrying the stub");
    tab->addColumn<float>("x", x, "stub global x [cm]");
    tab->addColumn<float>("y", y, "stub global y [cm]");
    tab->addColumn<float>("z", z, "stub global z [cm]");
    tab->addColumn<float>("r", r, "stub global cylindrical r = hypot(x,y) [cm]");
    tab->addColumn<float>("phi", phi, "stub global phi = atan2(y,x) [rad]", /*mantissaBits=*/12);
    tab->addColumn<float>("bend", bend,
                          "stub FE bend (full-strip units); the local r-phi angle, i.e. the OT's "
                          "analogue of the SmartPixels cotAlpha",
                          /*mantissaBits=*/12);
    if (doTruth_) {
      tab->addColumn<int32_t>("tpIdx", tpIdx,
                              "TRUTH-ONLY: TrackingParticle collection index of the associated TP, "
                              "or -1. Same meaning as L1TSmartPixelsCluster_tpIdx, so the two "
                              "tables join on it");
      tab->addColumn<float>("tpPt", tpPt, "TRUTH-ONLY: pT [GeV] of the associated TP", 10);
      tab->addColumn<bool>("tpGenuine", tpGenuine,
                           "TRUTH-ONLY: both clusters from the same TP");
      tab->addColumn<bool>("tpCombinatoric", tpCombinatoric,
                           "TRUTH-ONLY: clusters from DIFFERENT TPs -- a genuinely fake stub");
      tab->addColumn<bool>("tpUnknown", tpUnknown,
                           "TRUTH-ONLY: no TP association (noise, or out-of-time pileup)");
    }
    tab->setDoc(
        "EVERY Outer-Tracker stub offered to the L1 track finder (no trackIdx: most of these "
        "are on no track, which is the point). Post front-end pT gate, so already the reduced "
        "collection the finder receives. Use with L1TSmartPixelsCluster for IT-vs-OT "
        "combinatorics; the on-track subset is L1TTrackStub / L1TExtTrackStub.");
    iEvent.put(std::move(tab));
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("stubs", edm::InputTag("TTStubsFromPhase2TrackerDigis", "StubAccepted"))
        ->setComment("the SAME stub collection the tracklet chain consumes; do not point this at a "
                     "second stub producer or the table stops describing what the finder saw");
    desc.add<std::string>("tableName", "L1TOTStub");
    desc.add<bool>("doTruth", true)
        ->setComment("attach the TRUTH-ONLY tp* block. Needed for any efficiency or fake-rate "
                     "statement; without it the OT side can only be compared on combinatorics "
                     "while the IT side has truth, which is an asymmetry that biases conclusions");
    desc.add<edm::InputTag>("stubTruth",
                            edm::InputTag("TTStubAssociatorFromPixelDigis", "StubAccepted"))
        ->setComment("NO process name, so it resolves whether the associator ran in-job or its "
                     "output was persisted -- the same convention customizeSmartPixels_cff uses");
    desc.add<bool>("barrelOnly", false)
        ->setComment("keep only TOB stubs. Default false: the endcap rows are what make the "
                     "|eta| dependence of any OT-side cost visible");
    descriptions.addWithDefaultLabel(desc);
  }

private:
  const edm::EDGetTokenT<TTStubDetSetVec> stubsToken_;
  const edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> geomToken_;
  const edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> topoToken_;
  edm::EDGetTokenT<TTStubAssociationMap<Ref_Phase2TrackerDigi_>> truthToken_;
  const std::string tableName_;
  const bool barrelOnly_;
  const bool doTruth_;
};

DEFINE_FWK_MODULE(L1SmartPixelsAllStubTableProducer);
