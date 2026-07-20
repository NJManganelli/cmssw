// Nano producer that persists, per reference track, the REAL global positions of
// its Outer-Tracker (OT) stubs so the refit-replay visualization can anchor on the
// TRUE OT constraints (not a schematic) and refit research has the actual stub
// geometry available offline.
//
// It consumes a TTTrack collection (a SmartPixels REFERENCE collection: the prompt
// or extended tracklet tracks whose helix + stub refs seed the refit) and the
// TrackerGeometry / TrackerTopology from the EventSetup. For each track it walks
// getStubRefs(); for each TTStub it decodes the stub's global position from the
// inner cluster + detId via TrackerGeometry (the L1TrackNtupleMaker pattern,
// L1Trigger/TrackFindingTracklet/test/L1TrackNtupleMaker.cc ~1312-1346) and its
// tracker layer/disk id via TrackerTopology.
//
// Output: a single per-stub LINK table (one row per stub across ALL tracks, the
// L1SmartPixelsRefitTableProducer / L1SC4NGJetCands link-table pattern):
//   trackIdx  -> index of the owning track in the reference track table
//   x, y, z   -> stub global position [cm]
//   r, phi    -> cylindrical convenience (r = hypot(x,y), phi = atan2(y,x))
//   layer     -> tracker layer/disk id: barrel L1..L6 = 1..6, endcap disk d -> 10+d
//   isBarrel  -> 1 barrel, 0 endcap
//   bend      -> stub FE bend (full-strip units), cheap TTStub accessor
//
// KEY DESIGN CHOICE (spec RefitSidecarSpec.md §1 output-sync invariant; task Part B.2):
// the stub table is keyed to the REFERENCE L1TTrack / L1TExtTrack collection, NOT
// the refit variant collections. The digiRefit producer reuses the same stub refs on
// its output (setStubRefs), so the stubs are identical either way -- but the SEED's
// stubs are the physics constraint being visualized, and the reference table is the
// untouched anchor that trackIdx resolves against in both coexist and coopt. One
// table per reference collection (prompt + extended), suffixed, so the viz can show
// stubs for whichever seed it renders.
//
// The stub table is a standalone FlatTable (not an extension): trackIdx links back
// to the reference track table exactly like the refit "hit" table links to its
// variant track table.

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/L1TrackTrigger/interface/TTStub.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"

#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "Geometry/CommonTopologies/interface/GeomDet.h"
#include "Geometry/CommonTopologies/interface/Topology.h"
// The Phase-2 Outer Tracker reuses the strip-tracker subdetId numbering:
// OT barrel == StripSubdetector::TOB (5), OT endcap == StripSubdetector::TID (4).
// (Phase2Tracker::Subdetector::Barrel, used in L1TrackNtupleMaker/StubKiller, is
// exactly this TOB value; using StripSubdetector avoids a numbering-pkg include.)
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"

#include <cmath>
#include <cstdint>
#include <vector>

class L1SmartPixelsStubPosTableProducer : public edm::stream::EDProducer<> {
public:
  using L1Track = TTTrack<Ref_Phase2TrackerDigi_>;
  using TTTrackCollection = std::vector<L1Track>;

  explicit L1SmartPixelsStubPosTableProducer(const edm::ParameterSet& cfg)
      : tracksToken_(consumes<TTTrackCollection>(cfg.getParameter<edm::InputTag>("tracks"))),
        geomToken_(esConsumes<TrackerGeometry, TrackerDigiGeometryRecord>()),
        topoToken_(esConsumes<TrackerTopology, TrackerTopologyRcd>()),
        stubTableName_(cfg.getParameter<std::string>("stubTableName")),
        trackTableName_(cfg.getParameter<std::string>("trackTableName")) {
    produces<nanoaod::FlatTable>("stub");
  }

  void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override {
    const auto& tracks = iEvent.get(tracksToken_);
    const TrackerGeometry& geom = iSetup.getData(geomToken_);
    const TrackerTopology& topo = iSetup.getData(topoToken_);

    // one row per stub across all tracks (the LINK-table pattern)
    std::vector<int32_t> trackIdx;
    std::vector<uint8_t> layer;
    std::vector<bool> isBarrel;
    std::vector<uint32_t> detId;
    std::vector<float> x, y, z, r, phi, bend;

    for (size_t it = 0; it < tracks.size(); ++it) {
      const auto stubRefs = tracks[it].getStubRefs();
      for (const auto& stubRef : stubRefs) {
        // detId + global position from the inner cluster (L1TrackNtupleMaker pattern)
        const DetId detIdStub = geom.idToDet(stubRef->clusterRef(0)->getDetId())->geographicalId();
        const MeasurementPoint coords = stubRef->clusterRef(0)->findAverageLocalCoordinatesCentered();
        const GeomDet* theGeomDet = geom.idToDet(detIdStub);
        const GlobalPoint pos = theGeomDet->surface().toGlobal(theGeomDet->topology().localPosition(coords));

        const bool barrel = (detIdStub.subdetId() == StripSubdetector::TOB);
        const int lay = static_cast<int>(topo.layer(detIdStub));

        trackIdx.push_back(static_cast<int32_t>(it));
        detId.push_back(detIdStub.rawId());
        isBarrel.push_back(barrel);
        // barrel L1..L6 -> 1..6 ; endcap disk d -> 10+d (L1TrackNtupleMaker layerdisk)
        layer.push_back(static_cast<uint8_t>(barrel ? lay : 10 + lay));
        x.push_back(pos.x());
        y.push_back(pos.y());
        z.push_back(pos.z());
        r.push_back(std::hypot(pos.x(), pos.y()));
        phi.push_back(std::atan2(pos.y(), pos.x()));
        bend.push_back(static_cast<float>(stubRef->bendFE()));
      }
    }

    const size_t nStubRows = trackIdx.size();
    auto stubTable = std::make_unique<nanoaod::FlatTable>(nStubRows, stubTableName_, false, false);
    stubTable->addColumn<int32_t>(
        "trackIdx", trackIdx, "index of the owning track in the " + trackTableName_ + " table");
    stubTable->addColumn<uint8_t>(
        "layer", layer, "tracker layer/disk id: barrel L1..L6 = 1..6, endcap disk d = 10+d");
    stubTable->addColumn<bool>("isBarrel", isBarrel, "stub in the OT barrel (TOB) vs endcap");
    stubTable->addColumn<uint32_t>("detId", detId, "module rawId carrying the stub");
    stubTable->addColumn<float>("x", x, "stub global x [cm]");
    stubTable->addColumn<float>("y", y, "stub global y [cm]");
    stubTable->addColumn<float>("z", z, "stub global z [cm]");
    stubTable->addColumn<float>("r", r, "stub global cylindrical r = hypot(x,y) [cm]");
    stubTable->addColumn<float>("phi", phi, "stub global phi = atan2(y,x) [rad]", /*mantissaBits=*/12);
    stubTable->addColumn<float>("bend", bend, "stub FE bend (full-strip units)", /*mantissaBits=*/12);
    stubTable->setDoc("Real global positions of the OT stubs on the " + trackTableName_ +
                      " reference tracks (one row per stub; trackIdx links to that track table)");
    iEvent.put(std::move(stubTable), "stub");
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("tracks", edm::InputTag("l1tSmartPixelsTrackProducer", "Level1TTTracks"));
    desc.add<std::string>("stubTableName", "L1TSmartPixelsStub");
    desc.add<std::string>("trackTableName", "L1TSmartPixelsTrack");
    descriptions.addWithDefaultLabel(desc);
  }

private:
  const edm::EDGetTokenT<TTTrackCollection> tracksToken_;
  const edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> geomToken_;
  const edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> topoToken_;
  const std::string stubTableName_;
  const std::string trackTableName_;
};

DEFINE_FWK_MODULE(L1SmartPixelsStubPosTableProducer);
