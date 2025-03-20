// -*- C++ -*-
//
// Package:    L1Trigger/TrackFindingTracklet
// Class:      L1SmartPixelsTrackProducer
//
/**\class L1SmartPixelsTrackProducer L1SmartPixelsTrackProducer.cc L1Trigger/TrackFindingTracklet/plugins/L1SmartPixelsTrackProducer.cc

 Description: Produces a collection of TTTrackRefs (TTTrack_TrackWord info included) from OTTF L1tracks and SmartPixels parameterization/hits

 Implementation:
     Inputs:
         std::vector<TTTrack> - Each floating point TTTrack inside this collection inherits from
                                a bit-accurate TTTrack_TrackWord, used for emulation purposes.
          TODO: Add more inputs
     Outputs:
         std::vector<TTTrackRef> - A collection of TTTracks
*/
//
// Original Author:  Nick Manganelli
//         Created:  Thu, 24 Feb 2025 12:01:32 GMT

//#include "FWCore/PluginManager/interface/ModuleDef.h"
//#include "FWCore/Framework/interface/MakerMacros.h"
// #include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

///////////////////////
// DATA FORMATS HEADERS
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ref.h"

#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/L1TrackTrigger/interface/TTCluster.h"
#include "DataFormats/L1TrackTrigger/interface/TTStub.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTClusterAssociationMap.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTStubAssociationMap.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTTrackAssociationMap.h"
#include "Geometry/Records/interface/StackedTrackerGeometryRecord.h"

// #include "DataFormats/JetReco/interface/GenJetCollection.h"
// #include "DataFormats/JetReco/interface/GenJet.h"

////////////////////////////
// DETECTOR GEOMETRY HEADERS
// #include "MagneticField/Engine/interface/MagneticField.h"
// #include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
// #include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
// #include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
// #include "Geometry/TrackerGeometryBuilder/interface/RectangularPixelTopology.h"
// #include "Geometry/CommonDetUnit/interface/GeomDetType.h"
// #include "Geometry/CommonDetUnit/interface/GeomDet.h"

// #include "Geometry/CommonTopologies/interface/PixelGeomDetUnit.h"
// #include "Geometry/CommonTopologies/interface/PixelGeomDetType.h"
// #include "Geometry/TrackerGeometryBuilder/interface/PixelTopologyBuilder.h"
// #include "Geometry/Records/interface/StackedTrackerGeometryRecord.h"

////////////////
// PHYSICS TOOLS
#include "L1Trigger/TrackFindingTracklet/interface/HitPatternHelper.h"
// #include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CLHEP/Units/PhysicalConstants.h"

///////////////
// ROOT HEADERS
// #include <TROOT.h>
// #include <TCanvas.h>
// #include <TTree.h>
// #include <TFile.h>
// #include <TF1.h>
// #include <TH2F.h>
// #include <TH1F.h>

// Xilinx HLS includes
#include <ap_fixed.h>
#include <ap_int.h>

//////////////
// STD HEADERS
#include <memory>
#include <string>
#include <iostream>

//////////////
// NAMESPACES
using namespace std;
using namespace edm;

//////////////////////////////
//                          //
//     CLASS DEFINITION     //
//                          //
//////////////////////////////

// class L1TrackNtupleMaker : public one::EDAnalyzer<one::WatchRuns, one::SharedResources> {
class L1SmartPixelsTrackProducer : public global::EDProducer<> {
public:
  // Constructor/destructor
  explicit L1SmartPixelsTrackProducer(const edm::ParameterSet& iConfig);
  ~L1SmartPixelsTrackProducer() override;

  // Mandatory methods
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  // void beginJob() override;
  // void endJob() override;
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const;

  // void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) override;
  // void beginRun(const Run& iEvent, const EventSetup& iSetup) override {}
  // void endRun(const Run& iEvent, const EventSetup& iSetup) override {}

protected:
private:
  // // ----------constants, enums and typedefs ---------
  // // Relevant constants for the converted track word
  // enum TrackBitWidths {
  //   kPtSize = TTTrack_TrackWord::TrackBitWidths::kRinvSize - 1,  // Width of pt
  //   kPtMagSize = 9,                                              // Width of pt magnitude (unsigned)
  //   kEtaSize = TTTrack_TrackWord::TrackBitWidths::kTanlSize,     // Width of eta
  //   kEtaMagSize = 3,                                             // Width of eta magnitude (signed)
  // };
  static constexpr float MagConstant =
      CLHEP::c_light / 1.0E3;  //constant is 0.299792458; who knew c_light was in mm/ns?

  typedef TTTrack<Ref_Phase2TrackerDigi_> L1Track;
  typedef std::vector<L1Track> TTTrackCollection;
  // typedef edm::Handle<TTTrackCollection> TTTrackCollectionHandle;
  // typedef edm::Ref<TTTrackCollection> TTTrackRef;
  // typedef edm::RefVector<TTTrackCollection> TTTrackRefCollection;
  // typedef std::unique_ptr<TTTrackRefCollection> TTTrackRefCollectionUPtr;
  
  //-----------------------------------------------------------------------------------------------
  // Containers of parameters passed by python configuration file
  edm::ParameterSet config;

  int MyProcess;       // 11/13/211 for single electrons/muons/pions, 6/15 for pions from ttbar/taus, 1 for inclusive
  bool DebugMode;      // lots of debug printout statements
  int L1Tk_nPar;       // use 4 or 5 parameter track fit?
  int L1Tk_minNStub;     // require L1 tracks to have >= minNStub (this is mostly for tracklet purposes)
  std::string outputCollectionName_; // name of the output collection
  std::string smartPixelsEmulatorMode_; // mode of the emulator, e.g. passthrough, passthroughFloat, passthroughHW, trackingParticleTruth, ...

  edm::InputTag L1TrackInputTag;       // L1 track collection
  edm::InputTag MCTruthTrackInputTag;  // MC truth collection
  edm::InputTag MCTruthClusterInputTag;
  edm::InputTag L1StubInputTag;
  edm::InputTag MCTruthStubInputTag;
  edm::InputTag TrackingParticleInputTag;

  edm::EDGetTokenT<edmNew::DetSetVector<TTCluster<Ref_Phase2TrackerDigi_> > > ttClusterToken_;
  edm::EDGetTokenT<edmNew::DetSetVector<TTStub<Ref_Phase2TrackerDigi_> > > ttStubToken_;
  edm::EDGetTokenT<TTClusterAssociationMap<Ref_Phase2TrackerDigi_> > ttClusterMCTruthToken_;
  edm::EDGetTokenT<TTStubAssociationMap<Ref_Phase2TrackerDigi_> > ttStubMCTruthToken_;

  edm::EDGetTokenT<TTTrackCollection> ttTrackToken_;
  edm::EDGetTokenT<TTTrackAssociationMap<Ref_Phase2TrackerDigi_> > ttTrackMCTruthToken_;

  edm::EDGetTokenT<std::vector<TrackingParticle> > TrackingParticleToken_;
  edm::EDGetTokenT<std::vector<TrackingVertex> > TrackingVertexToken_;

  edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> getTokenTrackerGeom_;
  edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> getTokenTrackerTopo_;
  edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> getTokenBField_;
  edm::ESGetToken<hph::Setup, hph::SetupRcd> getTokenHPHSetup_;
};

//////////////////////////////////
//                              //
//     CLASS IMPLEMENTATION     //
//                              //
//////////////////////////////////

//////////////
// CONSTRUCTOR
L1SmartPixelsTrackProducer::L1SmartPixelsTrackProducer(edm::ParameterSet const& iConfig) : config(iConfig) {
  MyProcess = iConfig.getParameter<int>("MyProcess");
  DebugMode = iConfig.getParameter<bool>("DebugMode");
  L1Tk_nPar = iConfig.getParameter<int>("L1Tk_nPar");
  L1TrackInputTag = iConfig.getParameter<edm::InputTag>("L1TrackInputTag");
  MCTruthTrackInputTag = iConfig.getParameter<edm::InputTag>("MCTruthTrackInputTag");
  L1Tk_minNStub = iConfig.getParameter<int>("L1Tk_minNStub");
  outputCollectionName_ = iConfig.getParameter<std::string>("outputCollectionName");
  smartPixelsEmulatorMode_ = iConfig.getParameter<std::string>("smartPixelsEmulatorMode");

  L1StubInputTag = iConfig.getParameter<edm::InputTag>("L1StubInputTag");
  MCTruthClusterInputTag = iConfig.getParameter<edm::InputTag>("MCTruthClusterInputTag");
  MCTruthStubInputTag = iConfig.getParameter<edm::InputTag>("MCTruthStubInputTag");
  TrackingParticleInputTag = iConfig.getParameter<edm::InputTag>("TrackingParticleInputTag");

  ttTrackToken_ = consumes<TTTrackCollection>(L1TrackInputTag);
  ttTrackMCTruthToken_ = consumes<TTTrackAssociationMap<Ref_Phase2TrackerDigi_> >(MCTruthTrackInputTag);
  ttStubToken_ = consumes<edmNew::DetSetVector<TTStub<Ref_Phase2TrackerDigi_> > >(L1StubInputTag);
  ttClusterMCTruthToken_ = consumes<TTClusterAssociationMap<Ref_Phase2TrackerDigi_> >(MCTruthClusterInputTag);
  ttStubMCTruthToken_ = consumes<TTStubAssociationMap<Ref_Phase2TrackerDigi_> >(MCTruthStubInputTag);

  TrackingParticleToken_ = consumes<std::vector<TrackingParticle> >(TrackingParticleInputTag);

  getTokenTrackerGeom_ = esConsumes<TrackerGeometry, TrackerDigiGeometryRecord>();
  getTokenTrackerTopo_ = esConsumes<TrackerTopology, TrackerTopologyRcd>();
  getTokenBField_ = esConsumes<MagneticField, IdealMagneticFieldRecord>();
  getTokenHPHSetup_ = esConsumes<hph::Setup, hph::SetupRcd>();

  produces<TTTrackCollection>(outputCollectionName_);
}

// DESTRUCTOR
L1SmartPixelsTrackProducer::~L1SmartPixelsTrackProducer() {}

// PRODUCE
void L1SmartPixelsTrackProducer::produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const {

  if (!(MyProcess == 13 || MyProcess == 11 || MyProcess == 211 || MyProcess == 6 || MyProcess == 15 ||
        MyProcess == 1)) {
    edm::LogVerbatim("SmartPixelsTrackProducer") << "The specified MyProcess is invalid! Exiting...";
    return;
  }

  if (!(L1Tk_nPar == 4 || L1Tk_nPar == 5)) {
    edm::LogVerbatim("SmartPixelsTrackProducer") << "Invalid number of track parameters, specified L1Tk_nPar == " << L1Tk_nPar
                                 << " but only 4/5 are valid options! Exiting...";
    return;
  }

  // -----------------------------------------------------------------------------------------------
  // retrieve various containers
  // -----------------------------------------------------------------------------------------------

  // L1 tracks
  edm::Handle<TTTrackCollection> TTTrackHandle;
  iEvent.getByToken(ttTrackToken_, TTTrackHandle);

  // L1 stubs
  edm::Handle<edmNew::DetSetVector<TTStub<Ref_Phase2TrackerDigi_> > > TTStubHandle;
  iEvent.getByToken(ttStubToken_, TTStubHandle);

  // MC truth association maps
  edm::Handle<TTClusterAssociationMap<Ref_Phase2TrackerDigi_> > MCTruthTTClusterHandle;
  iEvent.getByToken(ttClusterMCTruthToken_, MCTruthTTClusterHandle);
  edm::Handle<TTStubAssociationMap<Ref_Phase2TrackerDigi_> > MCTruthTTStubHandle;
  iEvent.getByToken(ttStubMCTruthToken_, MCTruthTTStubHandle);
  edm::Handle<TTTrackAssociationMap<Ref_Phase2TrackerDigi_> > MCTruthTTTrackHandle;
  iEvent.getByToken(ttTrackMCTruthToken_, MCTruthTTTrackHandle);

  // tracking particles
  //edm::Handle<std::vector<TrackingParticle> > TrackingParticleHandle;
  //edm::Handle<std::vector<TrackingVertex> > TrackingVertexHandle;
  //iEvent.getByToken(TrackingParticleToken_, TrackingParticleHandle);
  //iEvent.getByToken(TrackingVertexToken_, TrackingVertexHandle);

  // -----------------------------------------------------------------------------------------------
  // more for TTStubs
  edm::ESHandle<TrackerGeometry> tGeomHandle = iSetup.getHandle(getTokenTrackerGeom_);

  edm::ESHandle<TrackerTopology> tTopoHandle = iSetup.getHandle(getTokenTrackerTopo_);

  edm::ESHandle<MagneticField> bFieldHandle = iSetup.getHandle(getTokenBField_);

  edm::ESHandle<hph::Setup> hphHandle = iSetup.getHandle(getTokenHPHSetup_);

  const TrackerTopology* const tTopo = tTopoHandle.product();
  const TrackerGeometry* const theTrackerGeom = tGeomHandle.product();
  const hph::Setup* hphSetup = hphHandle.product();

  // ----------------------------------------------------------------------------------------------
  // loop over L1 stubs
  // ----------------------------------------------------------------------------------------------

  for (auto gd = theTrackerGeom->dets().begin(); gd != theTrackerGeom->dets().end(); gd++) {
    DetId detid = (*gd)->geographicalId();
    if (detid.subdetId() != StripSubdetector::TOB && detid.subdetId() != StripSubdetector::TID)
      continue;
    if (!tTopo->isLower(detid))
      continue;                              // loop on the stacks: choose the lower arbitrarily
    DetId stackDetid = tTopo->stack(detid);  // Stub module detid

    if (TTStubHandle->find(stackDetid) == TTStubHandle->end())
      continue;

    // Get the DetSets of the Clusters
    edmNew::DetSet<TTStub<Ref_Phase2TrackerDigi_> > stubs = (*TTStubHandle)[stackDetid];
    const GeomDetUnit* det0 = theTrackerGeom->idToDetUnit(detid);
    const auto* theGeomDet = dynamic_cast<const PixelGeomDetUnit*>(det0);
    const PixelTopology* topol = dynamic_cast<const PixelTopology*>(&(theGeomDet->specificTopology()));

    // loop over stubs
    for (auto stubIter = stubs.begin(); stubIter != stubs.end(); ++stubIter) {
      edm::Ref<edmNew::DetSetVector<TTStub<Ref_Phase2TrackerDigi_> >, TTStub<Ref_Phase2TrackerDigi_> > tempStubPtr =
          edmNew::makeRefTo(TTStubHandle, stubIter);

      int isBarrel = 0;
      int layer = -999999;
      if (detid.subdetId() == StripSubdetector::TOB) {
        isBarrel = 1;
        layer = static_cast<int>(tTopo->layer(detid));
      } else if (detid.subdetId() == StripSubdetector::TID) {
        isBarrel = 0;
        layer = static_cast<int>(tTopo->layer(detid));
      } else {
        edm::LogVerbatim("SmartPixelsTrackProducer") << "WARNING -- neither TOB or TID stub, shouldn't happen...";
        layer = -1;
      }

      if (DebugMode) {
        edm::LogVerbatim("SmartPixelsTrackProducer") << "\n Stubs: layer = " << layer << "\n";
      }

      int isPSmodule = 0;
      if (topol->nrows() == 960)
        isPSmodule = 1;

      const unsigned int tobSide = tTopo->tobSide(detid);  // nonBarrel = 0, tiltedMinus = 1, tiltedPlus = 2, flat = 3
      int isTiltedBarrel = 0;
      if (isBarrel == 1 && (tobSide == 1 || tobSide == 2))
        isTiltedBarrel = 1;

      if (DebugMode) {
        edm::LogVerbatim("SmartPixelsTrackProducer") << "\n Stubs: isPSmodule = " << isPSmodule
                                      << " isTiltedBarrel = " << isTiltedBarrel << "\n";
      }

      MeasurementPoint coords = tempStubPtr->clusterRef(0)->findAverageLocalCoordinatesCentered();
      LocalPoint clustlp = topol->localPosition(coords);
      GlobalPoint posStub = theGeomDet->surface().toGlobal(clustlp);

      double tmp_stub_x = posStub.x();
      double tmp_stub_y = posStub.y();
      double tmp_stub_z = posStub.z();

      float trigDisplace = tempStubPtr->rawBend();
      float trigOffset = tempStubPtr->bendOffset();
      float trigPos = tempStubPtr->innerClusterPosition();
      float trigBend = tempStubPtr->bendFE();

      // matched to tracking particle?
      edm::Ptr<TrackingParticle> my_tp = MCTruthTTStubHandle->findTrackingParticlePtr(tempStubPtr);

      int myTP_pdgid = -999;
      float myTP_pt = -999;
      float myTP_eta = -999;
      float myTP_phi = -999;

      if (my_tp.isNull() == false) {
        int tmp_eventid = my_tp->eventId().event();

        if (tmp_eventid > 0)
          continue;  // this means stub from pileup track

        myTP_pdgid = my_tp->pdgId();
        myTP_pt = my_tp->p4().pt();
        myTP_eta = my_tp->p4().eta();
        myTP_phi = my_tp->p4().phi();
      }

      int tmp_stub_genuine = 0;
      if (MCTruthTTStubHandle->isGenuine(tempStubPtr))
        tmp_stub_genuine = 1;

      if (DebugMode)
      edm::LogVerbatim("SmartPixelsTrackProducer")
        << "myTP (pdgId, pt, eta, phi): " << myTP_pdgid << " " << myTP_pt << " " << myTP_eta << " " << myTP_phi
        << " isGenuine: " << tmp_stub_genuine;
    }
  }

  
  // ----------------------------------------------------------------------------------------------
  // loop over L1 tracks
  // ----------------------------------------------------------------------------------------------

  if (DebugMode) {
    edm::LogVerbatim("SmartPixelsTrackProducer") << "\n Loop over L1 tracks!";
    edm::LogVerbatim("SmartPixelsTrackProducer") << "\n Looking at " << L1Tk_nPar << "-parameter tracks!";
  }
  //output collection
  auto outputTracks = std::make_unique<TTTrackCollection>();
  int this_l1track = 0;
  TTTrackCollection::const_iterator iterL1Track;
  for (iterL1Track = TTTrackHandle->begin(); iterL1Track != TTTrackHandle->end(); iterL1Track++) {
    edm::Ptr<L1Track > l1track_ptr(TTTrackHandle, this_l1track);
    this_l1track++;

    float tmp_trk_pt = iterL1Track->momentum().perp();
    float tmp_trk_eta = iterL1Track->momentum().eta();
    float tmp_trk_phi = iterL1Track->momentum().phi();
    float tmp_trk_z0 = iterL1Track->z0();  //cm
    float tmp_trk_tanL = iterL1Track->tanL();
    bool usingNewKF = hphSetup->useNewKF();
    if (usingNewKF) {
      // Skip crazy tracks to avoid crash (as NewKF applies no cuts to kill them).
      constexpr float crazy_z0_cut = 30.;  // Cut to kill any crazy tracks found by New KF (which applies no cuts)
      if (fabs(tmp_trk_z0) > crazy_z0_cut)
        continue;
    }

    int tmp_trk_hitpattern = 0;
    tmp_trk_hitpattern = (int)iterL1Track->hitPattern();
    hph::HitPatternHelper hph(hphSetup, tmp_trk_hitpattern, tmp_trk_tanL, tmp_trk_z0);
    std::vector<int> hitpattern_expanded_binary = hph.binary();
    int tmp_trk_lhits_hitpattern = 0;
    int tmp_trk_dhits_hitpattern = 0;
    for (int i = 0; i < (int)hitpattern_expanded_binary.size(); i++) {
      if (hitpattern_expanded_binary[i]) {
        if (i < 6) {
          tmp_trk_lhits_hitpattern += pow(10, i);
        } else {
          tmp_trk_dhits_hitpattern += pow(10, i - 6);
        }
      }
    }
    int tmp_trk_nPSstub_hitpattern = hph.numPS();
    int tmp_trk_n2Sstub_hitpattern = hph.num2S();
    int tmp_trk_nLostPSstub_hitpattern = hph.numMissingPS();
    int tmp_trk_nLost2Sstub_hitpattern = hph.numMissing2S();
    int tmp_trk_nLoststub_V1_hitpattern = hph.numMissingInterior1();
    int tmp_trk_nLoststub_V2_hitpattern = hph.numMissingInterior2();

    float tmp_trk_d0 = -999;
    if (L1Tk_nPar == 5) {
      float tmp_trk_x0 = iterL1Track->POCA().x();
      float tmp_trk_y0 = iterL1Track->POCA().y();
      tmp_trk_d0 = tmp_trk_x0 * sin(tmp_trk_phi) - tmp_trk_y0 * cos(tmp_trk_phi);
    }

    float tmp_trk_chi2 = iterL1Track->chi2();
    float tmp_trk_chi2rphi = iterL1Track->chi2XY();
    float tmp_trk_chi2rz = iterL1Track->chi2Z();
    float tmp_trk_bendchi2 = iterL1Track->stubPtConsistency();
    float tmp_trk_MVA1 = iterL1Track->trkMVA1();

    std::vector<edm::Ref<edmNew::DetSetVector<TTStub<Ref_Phase2TrackerDigi_> >, TTStub<Ref_Phase2TrackerDigi_> > >
        stubRefs = iterL1Track->getStubRefs();
    int tmp_trk_nstub = (int)stubRefs.size();
    int ndof = 2 * tmp_trk_nstub - L1Tk_nPar;
    int ndofrphi = tmp_trk_nstub - L1Tk_nPar + 2;
    int ndofrz = tmp_trk_nstub - 2;
    float tmp_trk_chi2_dof = (float)tmp_trk_chi2 / ndof;
    float tmp_trk_chi2rphi_dof = (float)tmp_trk_chi2rphi / ndofrphi;
    float tmp_trk_chi2rz_dof = (float)tmp_trk_chi2rz / ndofrz;

    int tmp_trk_seed = (int)iterL1Track->trackSeedType();

    unsigned int tmp_trk_phiSector = iterL1Track->phiSector();
    int tmp_trk_etaSector = hph.etaSector();

    // ----------------------------------------------------------------------------------------------
    // loop over stubs on tracks

    //float tmp_trk_bend_chi2 = 0;
    int tmp_trk_dhits = 0;
    int tmp_trk_lhits = 0;

    if (true) {
      // loop over stubs
      for (int is = 0; is < tmp_trk_nstub; is++) {
        //detID of stub
        DetId detIdStub = theTrackerGeom->idToDet((stubRefs.at(is)->clusterRef(0))->getDetId())->geographicalId();

        MeasurementPoint coords = stubRefs.at(is)->clusterRef(0)->findAverageLocalCoordinatesCentered();
        const GeomDet* theGeomDet = theTrackerGeom->idToDet(detIdStub);
        Global3DPoint posStub = theGeomDet->surface().toGlobal(theGeomDet->topology().localPosition(coords));

        double x = posStub.x();
        double y = posStub.y();
        double z = posStub.z();

        int layer = -999999;
        if (detIdStub.subdetId() == StripSubdetector::TOB) {
          layer = static_cast<int>(tTopo->layer(detIdStub));
          if (DebugMode)
            edm::LogVerbatim("SmartPixelsTrackProducer")
                << "   stub in layer " << layer << " at position x y z = " << x << " " << y << " " << z;
          tmp_trk_lhits += pow(10, layer - 1);
        } else if (detIdStub.subdetId() == StripSubdetector::TID) {
          layer = static_cast<int>(tTopo->layer(detIdStub));
          if (DebugMode)
            edm::LogVerbatim("SmartPixelsTrackProducer")
                << "   stub in disk " << layer << " at position x y z = " << x << " " << y << " " << z;
          tmp_trk_dhits += pow(10, layer - 1);
        }

      }  //end loop over stubs
    }
    // ----------------------------------------------------------------------------------------------

    int tmp_trk_genuine = 0;
    int tmp_trk_loose = 0;
    int tmp_trk_unknown = 0;
    int tmp_trk_combinatoric = 0;
    if (MCTruthTTTrackHandle->isLooselyGenuine(l1track_ptr))
      tmp_trk_loose = 1;
    if (MCTruthTTTrackHandle->isGenuine(l1track_ptr))
      tmp_trk_genuine = 1;
    if (MCTruthTTTrackHandle->isUnknown(l1track_ptr))
      tmp_trk_unknown = 1;
    if (MCTruthTTTrackHandle->isCombinatoric(l1track_ptr))
      tmp_trk_combinatoric = 1;

    if (DebugMode) {
      edm::LogVerbatim("SmartPixelsTrackProducer") << "L1 track,"
                                    << " pt: " << tmp_trk_pt << " eta: " << tmp_trk_eta << " phi: " << tmp_trk_phi
                                    << " z0: " << tmp_trk_z0 << " d0: " << tmp_trk_d0 << " chi2: " << tmp_trk_chi2
                                    << " chi2rphi: " << tmp_trk_chi2rphi << " chi2rz: " << tmp_trk_chi2rz
                                    << " nstub: " << tmp_trk_nstub;
      if (tmp_trk_genuine)
        edm::LogVerbatim("SmartPixelsTrackProducer") << "    (is genuine)";
      if (tmp_trk_loose)
        edm::LogVerbatim("SmartPixelsTrackProducer") << "    (is loose)";
      if (tmp_trk_unknown)
        edm::LogVerbatim("SmartPixelsTrackProducer") << "    (is unknown)";
      if (tmp_trk_combinatoric)
        edm::LogVerbatim("SmartPixelsTrackProducer") << "    (is combinatoric)";
    }

    // ----------------------------------------------------------------------------------------------
    // for studying the fake rate
    // ----------------------------------------------------------------------------------------------

    edm::Ptr<TrackingParticle> my_tp = MCTruthTTTrackHandle->findTrackingParticlePtr(l1track_ptr);

    int myFake = 0;

    int tmp_matchtp_pdgid = -999;
    float tmp_matchtp_pt = -999;
    float tmp_matchtp_eta = -999;
    float tmp_matchtp_phi = -999;
    float tmp_matchtp_z0 = -999;
    float tmp_matchtp_dxy = -999;
    float tmp_matchtp_d0 = -999;

    if (my_tp.isNull()) {
      myFake = 0;
      if (DebugMode) {
        edm::LogVerbatim("SmartPixelsTrackProducer") << "TP not matched to track: myFake = " << myFake << " pdgId = " << tmp_matchtp_pdgid;
      }
    }
    else {
      int tmp_eventid = my_tp->eventId().event();

      if (tmp_eventid > 0)
        myFake = 2;
      else
        myFake = 1;

      tmp_matchtp_pdgid = my_tp->pdgId();
      tmp_matchtp_pt = my_tp->pt();
      tmp_matchtp_eta = my_tp->eta();
      tmp_matchtp_phi = my_tp->phi();

      float tmp_matchtp_vz = my_tp->vz();
      float tmp_matchtp_vx = my_tp->vx();
      float tmp_matchtp_vy = my_tp->vy();
      tmp_matchtp_dxy = sqrt(tmp_matchtp_vx * tmp_matchtp_vx + tmp_matchtp_vy * tmp_matchtp_vy);

      // ----------------------------------------------------------------------------------------------
      // get d0/z0 propagated back to the IP

      float tmp_matchtp_t = 1.0 / tan(2.0 * atan(exp(-tmp_matchtp_eta)));

      float delx = -tmp_matchtp_vx;
      float dely = -tmp_matchtp_vy;

      float b_field = bFieldHandle.product()->inTesla(GlobalPoint(0, 0, 0)).z();
      float c_converted = CLHEP::c_light / 1.0E5;
      float r2_inv = my_tp->charge() * c_converted * b_field / tmp_matchtp_pt / 2.0;

      float tmp_matchtp_x0p = delx - (1. / (2. * r2_inv) * sin(tmp_matchtp_phi));
      float tmp_matchtp_y0p = dely + (1. / (2. * r2_inv) * cos(tmp_matchtp_phi));
      float tmp_matchtp_rp = sqrt(tmp_matchtp_x0p * tmp_matchtp_x0p + tmp_matchtp_y0p * tmp_matchtp_y0p);
      tmp_matchtp_d0 = my_tp->charge() * tmp_matchtp_rp - (1. / (2. * r2_inv));

      static double pi = M_PI;
      float delphi = tmp_matchtp_phi - atan2(-r2_inv * tmp_matchtp_x0p, r2_inv * tmp_matchtp_y0p);
      if (delphi < -pi)
        delphi += 2.0 * pi;
      if (delphi > pi)
        delphi -= 2.0 * pi;
      tmp_matchtp_z0 = tmp_matchtp_vz + tmp_matchtp_t * delphi / (2.0 * r2_inv);
      // ----------------------------------------------------------------------------------------------

      if (DebugMode) {
        edm::LogVerbatim("SmartPixelsTrackProducer") << "TP matched to track has pt = " << my_tp->p4().pt()
                                      << " eta = " << my_tp->momentum().eta() << " phi = " << my_tp->momentum().phi()
                                      << " z0 = " << my_tp->vertex().z() << " pdgid = " << my_tp->pdgId()
                                      << " dxy = " << tmp_matchtp_dxy;
      }
    }


    // ----------------------------------------------------------------------------------------------
    // create new L1Tracks through one of the valid modes
    // ----------------------------------------------------------------------------------------------
    // Need all of https://github.com/cms-sw/cmssw/blob/b043e2f617238b96701cf4ffa10531133d1f3f7c/L1Trigger/TrackFindingTracklet/plugins/L1FPGATrackProducer.cc#L702-L764
    // to be able to create a new L1Track object from the TTTrack object in a completely functional way
    if(smartPixelsEmulatorMode_ == "passthrough") {
      //add the track to the output collection
      outputTracks->push_back(*iterL1Track);
    }
    else if(smartPixelsEmulatorMode_ == "passthroughFloat") {
      //emulate the track
      //create a L1Track object from the TTTrack by using the float value constructor of the L1Track class itself
      L1Track track = L1Track(iterL1Track->rInv(),
                              iterL1Track->phi(),
                              iterL1Track->tanL(),
                              iterL1Track->z0(),
                              iterL1Track->d0(),
                              iterL1Track->chi2XY(), //or chi2XYRed()
                              iterL1Track->chi2Z(), //or chi2ZRed()
                              iterL1Track->trkMVA1(),
                              iterL1Track->trkMVA2(),
                              iterL1Track->trkMVA3(),
                              iterL1Track->hitPattern(),
                              iterL1Track->nFitPars(), //or L1Tk_nPar,
                              3.8
      );
      track.setPhiSector(iterL1Track->phiSector());
      track.setTrackSeedType(iterL1Track->trackSeedType());
      track.setStubRefs(iterL1Track->getStubRefs());
      // pt consistency
      /*track.setStubPtConsistency(
        StubPtConsistency::getConsistency(track, theTrackerGeom, tTopo, settings_.bfield(), settings_.nHelixPar()));*/ //TOTEST w/ settings
      track.setTrackWordBits();

      // redo track quality?
      /*
      if (trackQuality_) {
        trackQualityModel_->setL1TrackQuality(aTrack);
      }
      // set track word again to set MVA variable from TTTrack into track word
      track.setTrackWordBits();
      */

      //add the track to the output collection
      outputTracks->push_back(track);
    }
    else if(smartPixelsEmulatorMode_ == "passthroughHW") {
      //emulate the track
      //create a L1Track object from the TTTrack by using the hw value constructor of the TrackWord parent class and copying the trackWord private member
      //this type of track might be vulnerable to having the method trackRef.setTrackWordBits() being called on it, which is 
      // supposed to set the trackword from float constructor values; this misses some important member variables though and 
      // will fail, so e.g. GTTInputProducer needs to have the flag "setTrackWordBits" set to false
      //From GTTFileReader.cc
      L1Track track = L1Track(iterL1Track->rInv(),
                              iterL1Track->phi(),
                              iterL1Track->tanL(),
                              iterL1Track->z0(),
                              iterL1Track->d0(),
                              iterL1Track->chi2XY(), //or chi2XYRed()
                              iterL1Track->chi2Z(), //or chi2ZRed()
                              iterL1Track->trkMVA1(),
                              iterL1Track->trkMVA2(),
                              iterL1Track->trkMVA3(),
                              iterL1Track->hitPattern(),
                              iterL1Track->nFitPars(), //or L1Tk_nPar,
                              3.8
      );
      track.setPhiSector(iterL1Track->phiSector());
      track.setTrackSeedType(iterL1Track->trackSeedType());
      track.setStubRefs(iterL1Track->getStubRefs());
      // pt consistency
      /*track.setStubPtConsistency(
        StubPtConsistency::getConsistency(track, theTrackerGeom, tTopo, settings_.bfield(), settings_.nHelixPar()));*/ //TOTEST w/ settings
      track.setTrackWordBits();

      // different from passthroughFloat we directly reset the trackWord_ member variable
      track.trackWord_ = static_cast<TTTrack_TrackWord::tkword_bs_t>(iterL1Track->getTrackWord());

      // redo track quality?
      /*
      if (trackQuality_) {
        trackQualityModel_->setL1TrackQuality(aTrack);
      }
      // set track word again to set MVA variable from TTTrack into track word
      // Instead of setTrackWordBits(), directly update the track word...
      */

      //add the track to the output collection
      outputTracks->push_back(track);
    }
    else if(smartPixelsEmulatorMode_ == "trackingParticleTruth") {
      // replace the relevant track parameters with the tracking particle truth
      if (my_tp.isNull()) {
        //if no tracking particle is matched to the track, make no changes, i.e. passthroughFloat behavior
        L1Track track = L1Track(iterL1Track->rInv(),
                                iterL1Track->phi(),
                                iterL1Track->tanL(),
                                iterL1Track->z0(),
                                iterL1Track->d0(),
                                iterL1Track->chi2XY(), //or chi2XYRed()
                                iterL1Track->chi2Z(), //or chi2ZRed()
                                iterL1Track->trkMVA1(),
                                iterL1Track->trkMVA2(),
                                iterL1Track->trkMVA3(),
                                iterL1Track->hitPattern(),
                                iterL1Track->nFitPars(), //or L1Tk_nPar,
                                3.8
        );
        track.setPhiSector(iterL1Track->phiSector());
        track.setTrackSeedType(iterL1Track->trackSeedType());
        track.setStubRefs(iterL1Track->getStubRefs());
        // pt consistency
        /*track.setStubPtConsistency(
          StubPtConsistency::getConsistency(track, theTrackerGeom, tTopo, settings_.bfield(), settings_.nHelixPar()));*/ //TOTEST w/ settings
        track.setTrackWordBits();
        //add the track to the output collection
        outputTracks->push_back(track);
      }
      else {
        // double thePT = std::abs(MagConstant / theRInv_ * aBField / 100.0);  // Rinv is in cm-1
        auto tmp_matcht_rInv = my_tp->charge() * MagConstant * 3.8 / (tmp_matchtp_pt * 100.0);
        auto tmp_matchtp_tanL = my_tp->p4().pz() / tmp_matchtp_pt;
        L1Track track = L1Track(tmp_matcht_rInv,
                                tmp_matchtp_phi,
                                tmp_matchtp_tanL, //TODO: replace
                                tmp_matchtp_z0,
                                tmp_matchtp_d0,
                                iterL1Track->chi2XY(), // Keep from original track?
                                iterL1Track->chi2Z(), // Keep from original track?
                                iterL1Track->trkMVA1(), // Keep from original track?
                                iterL1Track->trkMVA2(), // Keep from original track?
                                iterL1Track->trkMVA3(), // Keep from original track?
                                iterL1Track->hitPattern(), // Keep from original track?
                                iterL1Track->nFitPars(), // Keep from original track?
                                3.8
        );
        track.setPhiSector(iterL1Track->phiSector());
        track.setTrackSeedType(iterL1Track->trackSeedType());
        track.setStubRefs(iterL1Track->getStubRefs());
        // pt consistency
        /*track.setStubPtConsistency(
          StubPtConsistency::getConsistency(track, theTrackerGeom, tTopo, settings_.bfield(), settings_.nHelixPar()));*/ //TOTEST w/ settings
        track.setTrackWordBits();
        //add the track to the output collection
        outputTracks->push_back(track);
        std::cout << "track input pt: " << iterL1Track->momentum().perp() << " rInv: " << iterL1Track->rInv() << std::endl
                  << "tp input pt: " << tmp_matchtp_pt << " rInv: " << tmp_matcht_rInv << std::endl
                  << " track output pt: " << track.momentum().perp() << " rInv: " << track.rInv() << std::endl;
      }
    }
    else if(smartPixelsEmulatorMode_ == "toyDetectorParameterized") {
      // correct the track parameters towards the tracking particle truth level based on the parameterized toy detector
      continue;
    }
  }  //end track loop
  std::cout << "this_l1track counter = " << this_l1track << std::endl;
  std::cout << "outputTracks->size() = " << outputTracks->size() << std::endl;
  iEvent.put(std::move(outputTracks), outputCollectionName_);
}  // end of produce()

// FILLDESCRIPTIONS
void L1SmartPixelsTrackProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<int>("MyProcess", 1)->setComment("Process ID");
  desc.add<bool>("DebugMode", false)->setComment("Printout lots of debug statements");
  desc.add<int>("L1Tk_nPar", 4)->setComment("Use 4 or 5-parameter L1 tracking?");
  desc.add<int>("L1Tk_minNStub", 4)->setComment("L1 tracks with >= 4 stubs");
  desc.add<edm::InputTag>("L1TrackInputTag", edm::InputTag("l1tTTTracksFromTrackletEmulation", "Level1TTTracks"))->setComment("TTTrack input");
  desc.add<edm::InputTag>("MCTruthTrackInputTag", edm::InputTag("TTTrackAssociatorFromPixelDigis",  "Level1TTTracks"))->setComment("MCTruth input");
  desc.add<edm::InputTag>("L1StubInputTag", edm::InputTag("TTStubsFromPhase2TrackerDigis", "StubAccepted"));
  desc.add<edm::InputTag>("MCTruthClusterInputTag", edm::InputTag("TTClusterAssociatorFromPixelDigis", "ClusterAccepted"));
  desc.add<edm::InputTag>("MCTruthStubInputTag", edm::InputTag("TTStubAssociatorFromPixelDigis", "StubAccepted"));
  desc.add<edm::InputTag>("TrackingParticleInputTag", edm::InputTag("mix", "MergedTrackTruth"));
  desc.add<std::string>("outputCollectionName", "Level1TTTracksEmulation");
  desc.add<std::string>("smartPixelsEmulatorMode", "passthrough")->setComment("passthrough, passthroughFloat, passthroughHW, trackingParticleTruth");
  descriptions.addWithDefaultLabel(desc);
}
///////////////////////////
// DEFINE THIS AS A PLUG-IN
DEFINE_FWK_MODULE(L1SmartPixelsTrackProducer);
