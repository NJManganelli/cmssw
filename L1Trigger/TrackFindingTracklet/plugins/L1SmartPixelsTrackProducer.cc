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
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;
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

  // typedef TTTrack<Ref_Phase2TrackerDigi_> L1Track;
  // typedef std::vector<L1Track> TTTrackCollection;
  // typedef edm::Handle<TTTrackCollection> TTTrackCollectionHandle;
  // typedef edm::Ref<TTTrackCollection> TTTrackRef;
  // typedef edm::RefVector<TTTrackCollection> TTTrackRefCollection;
  // typedef std::unique_ptr<TTTrackRefCollection> TTTrackRefCollectionUPtr;
  
  //-----------------------------------------------------------------------------------------------
  // Containers of parameters passed by python configuration file
  edm::ParameterSet config;

  int MyProcess;       // 11/13/211 for single electrons/muons/pions, 6/15 for pions from ttbar/taus, 1 for inclusive
  bool DebugMode;      // lots of debug printout statements
  bool SaveStubs;      // option to save also stubs in the ntuples (makes them large...)
  int L1Tk_nPar;       // use 4 or 5 parameter track fit?
  int TP_minNStub;  // require TPs to have >= minNStub (defining efficiency denominator) (==0 means to only require >= 1 cluster)
  int TP_minNStubLayer;  // require TPs to have stubs in >= minNStubLayer layers/disks (defining efficiency denominator)
  double TP_minPt;       // save TPs with pt > minPt
  double TP_maxEta;      // save TPs with |eta| < maxEta
  double TP_maxZ0;       // save TPs with |z0| < maxZ0
  int L1Tk_minNStub;     // require L1 tracks to have >= minNStub (this is mostly for tracklet purposes)

  //bool TrackingInJets;  // do tracking in jets?

  edm::InputTag L1TrackInputTag;       // L1 track collection
  edm::InputTag MCTruthTrackInputTag;  // MC truth collection
  edm::InputTag MCTruthClusterInputTag;
  edm::InputTag L1StubInputTag;
  edm::InputTag MCTruthStubInputTag;
  edm::InputTag TrackingParticleInputTag;
  //edm::InputTag TrackingVertexInputTag;
  //edm::InputTag GenJetInputTag;

  edm::EDGetTokenT<edmNew::DetSetVector<TTCluster<Ref_Phase2TrackerDigi_> > > ttClusterToken_;
  edm::EDGetTokenT<edmNew::DetSetVector<TTStub<Ref_Phase2TrackerDigi_> > > ttStubToken_;
  edm::EDGetTokenT<TTClusterAssociationMap<Ref_Phase2TrackerDigi_> > ttClusterMCTruthToken_;
  edm::EDGetTokenT<TTStubAssociationMap<Ref_Phase2TrackerDigi_> > ttStubMCTruthToken_;

  edm::EDGetTokenT<std::vector<TTTrack<Ref_Phase2TrackerDigi_> > > ttTrackToken_;
  edm::EDGetTokenT<TTTrackAssociationMap<Ref_Phase2TrackerDigi_> > ttTrackMCTruthToken_;

  edm::EDGetTokenT<std::vector<TrackingParticle> > TrackingParticleToken_;
  edm::EDGetTokenT<std::vector<TrackingVertex> > TrackingVertexToken_;

  //edm::EDGetTokenT<std::vector<reco::GenJet> > GenJetToken_;

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
  SaveStubs = iConfig.getParameter<bool>("SaveStubs");
  L1Tk_nPar = iConfig.getParameter<int>("L1Tk_nPar");
  TP_minNStub = iConfig.getParameter<int>("TP_minNStub");
  TP_minNStubLayer = iConfig.getParameter<int>("TP_minNStubLayer");
  TP_minPt = iConfig.getParameter<double>("TP_minPt");
  TP_maxEta = iConfig.getParameter<double>("TP_maxEta");
  TP_maxZ0 = iConfig.getParameter<double>("TP_maxZ0");
  L1TrackInputTag = iConfig.getParameter<edm::InputTag>("L1TrackInputTag");
  MCTruthTrackInputTag = iConfig.getParameter<edm::InputTag>("MCTruthTrackInputTag");
  L1Tk_minNStub = iConfig.getParameter<int>("L1Tk_minNStub");

  //TrackingInJets = iConfig.getParameter<bool>("TrackingInJets");

  L1StubInputTag = iConfig.getParameter<edm::InputTag>("L1StubInputTag");
  MCTruthClusterInputTag = iConfig.getParameter<edm::InputTag>("MCTruthClusterInputTag");
  MCTruthStubInputTag = iConfig.getParameter<edm::InputTag>("MCTruthStubInputTag");
  TrackingParticleInputTag = iConfig.getParameter<edm::InputTag>("TrackingParticleInputTag");
  //TrackingVertexInputTag = iConfig.getParameter<edm::InputTag>("TrackingVertexInputTag");
  //GenJetInputTag = iConfig.getParameter<edm::InputTag>("GenJetInputTag");

  ttTrackToken_ = consumes<std::vector<TTTrack<Ref_Phase2TrackerDigi_> > >(L1TrackInputTag);
  ttTrackMCTruthToken_ = consumes<TTTrackAssociationMap<Ref_Phase2TrackerDigi_> >(MCTruthTrackInputTag);
  ttStubToken_ = consumes<edmNew::DetSetVector<TTStub<Ref_Phase2TrackerDigi_> > >(L1StubInputTag);
  ttClusterMCTruthToken_ = consumes<TTClusterAssociationMap<Ref_Phase2TrackerDigi_> >(MCTruthClusterInputTag);
  ttStubMCTruthToken_ = consumes<TTStubAssociationMap<Ref_Phase2TrackerDigi_> >(MCTruthStubInputTag);

  TrackingParticleToken_ = consumes<std::vector<TrackingParticle> >(TrackingParticleInputTag);
  //TrackingVertexToken_ = consumes<std::vector<TrackingVertex> >(TrackingVertexInputTag);
  //GenJetToken_ = consumes<std::vector<reco::GenJet> >(GenJetInputTag);

  getTokenTrackerGeom_ = esConsumes<TrackerGeometry, TrackerDigiGeometryRecord>();
  getTokenTrackerTopo_ = esConsumes<TrackerTopology, TrackerTopologyRcd>();
  getTokenBField_ = esConsumes<MagneticField, IdealMagneticFieldRecord>();
  getTokenHPHSetup_ = esConsumes<hph::Setup, hph::SetupRcd>();
}

// DESTRUCTOR
L1SmartPixelsTrackProducer::~L1SmartPixelsTrackProducer() {}

// PRODUCE
void L1SmartPixelsTrackProducer::produce(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  if (!(MyProcess == 13 || MyProcess == 11 || MyProcess == 211 || MyProcess == 6 || MyProcess == 15 ||
        MyProcess == 1)) {
    edm::LogVerbatim("Tracklet") << "The specified MyProcess is invalid! Exiting...";
    return;
  }

  if (!(L1Tk_nPar == 4 || L1Tk_nPar == 5)) {
    edm::LogVerbatim("Tracklet") << "Invalid number of track parameters, specified L1Tk_nPar == " << L1Tk_nPar
                                 << " but only 4/5 are valid options! Exiting...";
    return;
  }

  // -----------------------------------------------------------------------------------------------
  // retrieve various containers
  // -----------------------------------------------------------------------------------------------

  // L1 tracks
  edm::Handle<std::vector<TTTrack<Ref_Phase2TrackerDigi_> > > TTTrackHandle;
  iEvent.getByToken(ttTrackToken_, TTTrackHandle);

  // L1 stubs
  edm::Handle<edmNew::DetSetVector<TTStub<Ref_Phase2TrackerDigi_> > > TTStubHandle;
  if (SaveStubs)
    iEvent.getByToken(ttStubToken_, TTStubHandle);

  // MC truth association maps
  edm::Handle<TTClusterAssociationMap<Ref_Phase2TrackerDigi_> > MCTruthTTClusterHandle;
  iEvent.getByToken(ttClusterMCTruthToken_, MCTruthTTClusterHandle);
  edm::Handle<TTStubAssociationMap<Ref_Phase2TrackerDigi_> > MCTruthTTStubHandle;
  iEvent.getByToken(ttStubMCTruthToken_, MCTruthTTStubHandle);
  edm::Handle<TTTrackAssociationMap<Ref_Phase2TrackerDigi_> > MCTruthTTTrackHandle;
  iEvent.getByToken(ttTrackMCTruthToken_, MCTruthTTTrackHandle);

  // tracking particles
  edm::Handle<std::vector<TrackingParticle> > TrackingParticleHandle;
  edm::Handle<std::vector<TrackingVertex> > TrackingVertexHandle;
  iEvent.getByToken(TrackingParticleToken_, TrackingParticleHandle);
  iEvent.getByToken(TrackingVertexToken_, TrackingVertexHandle);

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

  if (SaveStubs) {
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
          edm::LogVerbatim("Tracklet") << "WARNING -- neither TOB or TID stub, shouldn't happen...";
          layer = -1;
        }

        int isPSmodule = 0;
        if (topol->nrows() == 960)
          isPSmodule = 1;

        const unsigned int tobSide = tTopo->tobSide(detid);  // nonBarrel = 0, tiltedMinus = 1, tiltedPlus = 2, flat = 3
        int isTiltedBarrel = 0;
        if (isBarrel == 1 && (tobSide == 1 || tobSide == 2))
          isTiltedBarrel = 1;

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

        m_allstub_matchTP_pdgid->push_back(myTP_pdgid);
        m_allstub_matchTP_pt->push_back(myTP_pt);
        m_allstub_matchTP_eta->push_back(myTP_eta);
        m_allstub_matchTP_phi->push_back(myTP_phi);

        int tmp_stub_genuine = 0;
        if (MCTruthTTStubHandle->isGenuine(tempStubPtr))
          tmp_stub_genuine = 1;

        m_allstub_genuine->push_back(tmp_stub_genuine);
      }
    }
  }

  
  // ----------------------------------------------------------------------------------------------
  // loop over L1 tracks
  // ----------------------------------------------------------------------------------------------

  if (DebugMode) {
    edm::LogVerbatim("Tracklet") << "\n Loop over L1 tracks!";
    edm::LogVerbatim("Tracklet") << "\n Looking at " << L1Tk_nPar << "-parameter tracks!";
  }

  int this_l1track = 0;
  std::vector<TTTrack<Ref_Phase2TrackerDigi_> >::const_iterator iterL1Track;
  for (iterL1Track = TTTrackHandle->begin(); iterL1Track != TTTrackHandle->end(); iterL1Track++) {
    edm::Ptr<TTTrack<Ref_Phase2TrackerDigi_> > l1track_ptr(TTTrackHandle, this_l1track);
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

    int tmp_trk_seed = 0;
    tmp_trk_seed = (int)iterL1Track->trackSeedType();

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
            edm::LogVerbatim("Tracklet")
                << "   stub in layer " << layer << " at position x y z = " << x << " " << y << " " << z;
          tmp_trk_lhits += pow(10, layer - 1);
        } else if (detIdStub.subdetId() == StripSubdetector::TID) {
          layer = static_cast<int>(tTopo->layer(detIdStub));
          if (DebugMode)
            edm::LogVerbatim("Tracklet")
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
      edm::LogVerbatim("Tracklet") << "L1 track,"
                                    << " pt: " << tmp_trk_pt << " eta: " << tmp_trk_eta << " phi: " << tmp_trk_phi
                                    << " z0: " << tmp_trk_z0 << " chi2: " << tmp_trk_chi2
                                    << " chi2rphi: " << tmp_trk_chi2rphi << " chi2rz: " << tmp_trk_chi2rz
                                    << " nstub: " << tmp_trk_nstub;
      if (tmp_trk_genuine)
        edm::LogVerbatim("Tracklet") << "    (is genuine)";
      if (tmp_trk_unknown)
        edm::LogVerbatim("Tracklet") << "    (is unknown)";
      if (tmp_trk_combinatoric)
        edm::LogVerbatim("Tracklet") << "    (is combinatoric)";
    }

    m_trk_pt->push_back(tmp_trk_pt);
    m_trk_eta->push_back(tmp_trk_eta);
    m_trk_phi->push_back(tmp_trk_phi);
    m_trk_z0->push_back(tmp_trk_z0);
    if (L1Tk_nPar == 5)
      m_trk_d0->push_back(tmp_trk_d0);
    else
      m_trk_d0->push_back(999.);
    m_trk_chi2->push_back(tmp_trk_chi2);
    m_trk_chi2_dof->push_back(tmp_trk_chi2_dof);
    m_trk_chi2rphi->push_back(tmp_trk_chi2rphi);
    m_trk_chi2rphi_dof->push_back(tmp_trk_chi2rphi_dof);
    m_trk_chi2rz->push_back(tmp_trk_chi2rz);
    m_trk_chi2rz_dof->push_back(tmp_trk_chi2rz_dof);
    m_trk_bendchi2->push_back(tmp_trk_bendchi2);
    m_trk_MVA1->push_back(tmp_trk_MVA1);
    m_trk_nstub->push_back(tmp_trk_nstub);
    m_trk_dhits->push_back(tmp_trk_dhits);
    m_trk_lhits->push_back(tmp_trk_lhits);
    m_trk_seed->push_back(tmp_trk_seed);
    m_trk_hitpattern->push_back(tmp_trk_hitpattern);
    m_trk_lhits_hitpattern->push_back(tmp_trk_lhits_hitpattern);
    m_trk_dhits_hitpattern->push_back(tmp_trk_dhits_hitpattern);
    m_trk_nPSstub_hitpattern->push_back(tmp_trk_nPSstub_hitpattern);
    m_trk_n2Sstub_hitpattern->push_back(tmp_trk_n2Sstub_hitpattern);
    m_trk_nLostPSstub_hitpattern->push_back(tmp_trk_nLostPSstub_hitpattern);
    m_trk_nLost2Sstub_hitpattern->push_back(tmp_trk_nLost2Sstub_hitpattern);
    m_trk_nLoststub_V1_hitpattern->push_back(tmp_trk_nLoststub_V1_hitpattern);
    m_trk_nLoststub_V2_hitpattern->push_back(tmp_trk_nLoststub_V2_hitpattern);
    m_trk_phiSector->push_back(tmp_trk_phiSector);
    m_trk_etaSector->push_back(tmp_trk_etaSector);
    m_trk_genuine->push_back(tmp_trk_genuine);
    m_trk_loose->push_back(tmp_trk_loose);
    m_trk_unknown->push_back(tmp_trk_unknown);
    m_trk_combinatoric->push_back(tmp_trk_combinatoric);

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

    if (my_tp.isNull())
      myFake = 0;
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
        edm::LogVerbatim("Tracklet") << "TP matched to track has pt = " << my_tp->p4().pt()
                                      << " eta = " << my_tp->momentum().eta() << " phi = " << my_tp->momentum().phi()
                                      << " z0 = " << my_tp->vertex().z() << " pdgid = " << my_tp->pdgId()
                                      << " dxy = " << tmp_matchtp_dxy;
      }
    }

    m_trk_fake->push_back(myFake);

    m_trk_matchtp_pdgid->push_back(tmp_matchtp_pdgid);
    m_trk_matchtp_pt->push_back(tmp_matchtp_pt);
    m_trk_matchtp_eta->push_back(tmp_matchtp_eta);
    m_trk_matchtp_phi->push_back(tmp_matchtp_phi);
    m_trk_matchtp_z0->push_back(tmp_matchtp_z0);
    m_trk_matchtp_dxy->push_back(tmp_matchtp_dxy);
    m_trk_matchtp_d0->push_back(tmp_matchtp_d0);

  }  //end track loop


  // ----------------------------------------------------------------------------------------------
  // loop over tracking particles
  // ----------------------------------------------------------------------------------------------

  if (DebugMode)
    edm::LogVerbatim("Tracklet") << "\n Loop over tracking particles!";

  int this_tp = 0;
  std::vector<TrackingParticle>::const_iterator iterTP;
  for (iterTP = TrackingParticleHandle->begin(); iterTP != TrackingParticleHandle->end(); ++iterTP) {
    edm::Ptr<TrackingParticle> tp_ptr(TrackingParticleHandle, this_tp);
    this_tp++;

    int tmp_eventid = iterTP->eventId().event();
    if (MyProcess != 1 && tmp_eventid > 0)
      continue;  //only care about tracking particles from the primary interaction (except for MyProcess==1, i.e. looking at all TPs)

    float tmp_tp_pt = iterTP->pt();
    float tmp_tp_charge = iterTP->charge();
    float tmp_tp_eta = iterTP->eta();

    if (tmp_tp_pt < TP_minPt)  // Save CPU by applying these cuts here.
      continue;
    if (tmp_tp_charge == 0.)
      continue;
    if (std::abs(tmp_tp_eta) > TP_maxEta)
      continue;

    float tmp_tp_phi = iterTP->phi();
    float tmp_tp_vz = iterTP->vz();
    float tmp_tp_vx = iterTP->vx();
    float tmp_tp_vy = iterTP->vy();
    int tmp_tp_pdgid = iterTP->pdgId();
    float tmp_tp_z0_prod = tmp_tp_vz;
    float tmp_tp_d0_prod = tmp_tp_vx * sin(tmp_tp_phi) - tmp_tp_vy * cos(tmp_tp_phi);

    // ----------------------------------------------------------------------------------------------
    // get d0/z0 propagated back to the IP

    float tmp_tp_t = 1.0 / tan(2.0 * atan(exp(-tmp_tp_eta)));

    float delx = -tmp_tp_vx;
    float dely = -tmp_tp_vy;

    float b_field = bFieldHandle.product()->inTesla(GlobalPoint(0, 0, 0)).z();
    float c_converted = CLHEP::c_light / 1.0E5;
    float r2_inv = tmp_tp_charge * c_converted * b_field / tmp_tp_pt / 2.0;

    float tmp_tp_x0p = delx - (1. / (2. * r2_inv) * sin(tmp_tp_phi));
    float tmp_tp_y0p = dely + (1. / (2. * r2_inv) * cos(tmp_tp_phi));
    float tmp_tp_rp = sqrt(tmp_tp_x0p * tmp_tp_x0p + tmp_tp_y0p * tmp_tp_y0p);
    float tmp_tp_d0 = tmp_tp_charge * tmp_tp_rp - (1. / (2. * r2_inv));

    static double pi = M_PI;
    float delphi = tmp_tp_phi - atan2(-r2_inv * tmp_tp_x0p, r2_inv * tmp_tp_y0p);
    if (delphi < -pi)
      delphi += 2.0 * pi;
    if (delphi > pi)
      delphi -= 2.0 * pi;
    float tmp_tp_z0 = tmp_tp_vz + tmp_tp_t * delphi / (2.0 * r2_inv);
    // ----------------------------------------------------------------------------------------------

    if (MyProcess == 13 && abs(tmp_tp_pdgid) != 13)
      continue;
    if (MyProcess == 11 && abs(tmp_tp_pdgid) != 11)
      continue;
    if ((MyProcess == 6 || MyProcess == 15 || MyProcess == 211) && abs(tmp_tp_pdgid) != 211)
      continue;

    if (std::abs(tmp_tp_z0) > TP_maxZ0)
      continue;

    // for pions in ttbar, only consider TPs coming from near the IP!
    float dxy = sqrt(tmp_tp_vx * tmp_tp_vx + tmp_tp_vy * tmp_tp_vy);
    float tmp_tp_dxy = dxy;
    if (MyProcess == 6 && (dxy > 1.0))
      continue;

    if (DebugMode)
      edm::LogVerbatim("Tracklet") << "Tracking particle, pt: " << tmp_tp_pt << " eta: " << tmp_tp_eta
                                   << " phi: " << tmp_tp_phi << " z0: " << tmp_tp_z0 << " d0: " << tmp_tp_d0
                                   << " z_prod: " << tmp_tp_z0_prod << " d_prod: " << tmp_tp_d0_prod
                                   << " pdgid: " << tmp_tp_pdgid << " eventID: " << iterTP->eventId().event()
                                   << " ttclusters " << MCTruthTTClusterHandle->findTTClusterRefs(tp_ptr).size()
                                   << " ttstubs " << MCTruthTTStubHandle->findTTStubRefs(tp_ptr).size() << " tttracks "
                                   << MCTruthTTTrackHandle->findTTTrackPtrs(tp_ptr).size();

    // ----------------------------------------------------------------------------------------------
    // only consider TPs associated with >= 1 cluster, or >= X stubs, or have stubs in >= X layers (configurable options)

    if (MCTruthTTClusterHandle->findTTClusterRefs(tp_ptr).empty()) {
      if (DebugMode)
        edm::LogVerbatim("Tracklet") << "No matching TTClusters for TP, continuing...";
      continue;
    }

    std::vector<edm::Ref<edmNew::DetSetVector<TTStub<Ref_Phase2TrackerDigi_> >, TTStub<Ref_Phase2TrackerDigi_> > >
        theStubRefs = MCTruthTTStubHandle->findTTStubRefs(tp_ptr);
    int nStubTP = (int)theStubRefs.size();

    // how many layers/disks have stubs?
    int hasStubInLayer[11] = {0};
    for (auto& theStubRef : theStubRefs) {
      DetId detid(theStubRef->getDetId());

      int layer = -1;
      if (detid.subdetId() == StripSubdetector::TOB) {
        layer = static_cast<int>(tTopo->layer(detid)) - 1;  //fill in array as entries 0-5
      } else if (detid.subdetId() == StripSubdetector::TID) {
        layer = static_cast<int>(tTopo->layer(detid)) + 5;  //fill in array as entries 6-10
      }

      //bool isPS = (theTrackerGeom->getDetectorType(detid)==TrackerGeometry::ModuleType::Ph2PSP);

      //treat genuine stubs separately (==2 is genuine, ==1 is not)
      if (MCTruthTTStubHandle->findTrackingParticlePtr(theStubRef).isNull() && hasStubInLayer[layer] < 2)
        hasStubInLayer[layer] = 1;
      else
        hasStubInLayer[layer] = 2;
    }

    int nStubLayerTP = 0;
    int nStubLayerTP_g = 0;
    for (int isum : hasStubInLayer) {
      if (isum >= 1)
        nStubLayerTP += 1;
      if (isum == 2)
        nStubLayerTP_g += 1;
    }

    if (DebugMode)
      edm::LogVerbatim("Tracklet") << "TP is associated with " << nStubTP << " stubs, and has stubs in " << nStubLayerTP
                                   << " different layers/disks, and has GENUINE stubs in " << nStubLayerTP_g
                                   << " layers ";

    if (TP_minNStub > 0) {
      if (DebugMode)
        edm::LogVerbatim("Tracklet") << "Only consider TPs with >= " << TP_minNStub << " stubs";
      if (nStubTP < TP_minNStub) {
        if (DebugMode)
          edm::LogVerbatim("Tracklet") << "TP fails minimum nbr stubs requirement! Continuing...";
        continue;
      }
    }
    if (TP_minNStubLayer > 0) {
      if (DebugMode)
        edm::LogVerbatim("Tracklet") << "Only consider TPs with stubs in >= " << TP_minNStubLayer << " layers/disks";
      if (nStubLayerTP < TP_minNStubLayer) {
        if (DebugMode)
          edm::LogVerbatim("Tracklet") << "TP fails stubs in minimum nbr of layers/disks requirement! Continuing...";
        continue;
      }
    }

    // ----------------------------------------------------------------------------------------------
    // look for L1 tracks matched to the tracking particle

    std::vector<edm::Ptr<TTTrack<Ref_Phase2TrackerDigi_> > > matchedTracks =
        MCTruthTTTrackHandle->findTTTrackPtrs(tp_ptr);

    int nMatch = 0;
    int i_track = -1;
    float i_chi2dof = 99999;

    if (!matchedTracks.empty()) {
      if (DebugMode && (matchedTracks.size() > 1))
        edm::LogVerbatim("Tracklet") << "TrackingParticle has more than one matched L1 track!";

      // ----------------------------------------------------------------------------------------------
      // loop over matched L1 tracks
      // here, "match" means tracks that can be associated to a TrackingParticle with at least one hit of at least one of its clusters
      // https://twiki.cern.ch/twiki/bin/viewauth/CMS/SLHCTrackerTriggerSWTools#MC_truth_for_TTTrack

      for (int it = 0; it < (int)matchedTracks.size(); it++) {
        bool tmp_trk_genuine = false;
        bool tmp_trk_loosegenuine = false;
        if (MCTruthTTTrackHandle->isGenuine(matchedTracks.at(it)))
          tmp_trk_genuine = true;
        if (MCTruthTTTrackHandle->isLooselyGenuine(matchedTracks.at(it)))
          tmp_trk_loosegenuine = true;
        if (!tmp_trk_loosegenuine)
          continue;

        if (DebugMode) {
          if (MCTruthTTTrackHandle->findTrackingParticlePtr(matchedTracks.at(it)).isNull()) {
            edm::LogVerbatim("Tracklet") << "track matched to TP is NOT uniquely matched to a TP";
          } else {
            edm::Ptr<TrackingParticle> my_tp = MCTruthTTTrackHandle->findTrackingParticlePtr(matchedTracks.at(it));
            edm::LogVerbatim("Tracklet") << "TP matched to track matched to TP ... tp pt = " << my_tp->p4().pt()
                                         << " eta = " << my_tp->momentum().eta() << " phi = " << my_tp->momentum().phi()
                                         << " z0 = " << my_tp->vertex().z();
          }
          edm::LogVerbatim("Tracklet") << "   ... matched L1 track has pt = " << matchedTracks.at(it)->momentum().perp()
                                       << " eta = " << matchedTracks.at(it)->momentum().eta()
                                       << " phi = " << matchedTracks.at(it)->momentum().phi()
                                       << " chi2 = " << matchedTracks.at(it)->chi2()
                                       << " consistency = " << matchedTracks.at(it)->stubPtConsistency()
                                       << " z0 = " << matchedTracks.at(it)->z0()
                                       << " nstub = " << matchedTracks.at(it)->getStubRefs().size();
          if (tmp_trk_genuine)
            edm::LogVerbatim("Tracklet") << "    (genuine!) ";
          if (tmp_trk_loosegenuine)
            edm::LogVerbatim("Tracklet") << "    (loose genuine!) ";
        }

        // ----------------------------------------------------------------------------------------------
        // further require L1 track to be (loosely) genuine, that there is only one TP matched to the track
        // + have >= L1Tk_minNStub stubs for it to be a valid match (only relevant is your track collection
        // e.g. stores 3-stub tracks but at plot level you require >= 4 stubs (--> tracklet case)

        std::vector<edm::Ref<edmNew::DetSetVector<TTStub<Ref_Phase2TrackerDigi_> >, TTStub<Ref_Phase2TrackerDigi_> > >
            stubRefs = matchedTracks.at(it)->getStubRefs();
        int tmp_trk_nstub = stubRefs.size();

        if (tmp_trk_nstub < L1Tk_minNStub)
          continue;

        /*
        // PS stubs
        int tmp_trk_nPSstub = 0;
        for (int is=0; is<tmp_trk_nstub; is++) {
          DetId detIdStub = theTrackerGeom->idToDet( (stubRefs.at(is)->clusterRef(0))->getDetId() )->geographicalId();
          DetId stackDetid = tTopo->stack(detIdStub);
          bool isPS = (theTrackerGeom->getDetectorType(stackDetid)==TrackerGeometry::ModuleType::Ph2PSP);
          if (isPS) tmp_trk_nPSstub++;
        }
        */

        float dmatch_pt = 999;
        float dmatch_eta = 999;
        float dmatch_phi = 999;
        int match_id = 999;

        edm::Ptr<TrackingParticle> my_tp = MCTruthTTTrackHandle->findTrackingParticlePtr(matchedTracks.at(it));
        dmatch_pt = std::abs(my_tp->p4().pt() - tmp_tp_pt);
        dmatch_eta = std::abs(my_tp->p4().eta() - tmp_tp_eta);
        dmatch_phi = std::abs(my_tp->p4().phi() - tmp_tp_phi);
        match_id = my_tp->pdgId();

        float tmp_trk_chi2dof = (matchedTracks.at(it)->chi2()) / (2 * tmp_trk_nstub - L1Tk_nPar);

        // ensure that track is uniquely matched to the TP we are looking at!
        if (dmatch_pt < 0.1 && dmatch_eta < 0.1 && dmatch_phi < 0.1 && tmp_tp_pdgid == match_id && tmp_trk_genuine) {
          nMatch++;
          if (i_track < 0 || tmp_trk_chi2dof < i_chi2dof) {
            i_track = it;
            i_chi2dof = tmp_trk_chi2dof;
          }
        }

      }  // end loop over matched L1 tracks

    }  // end has at least 1 matched L1 track
    // ----------------------------------------------------------------------------------------------

    float tmp_matchtrk_pt = -999;
    float tmp_matchtrk_eta = -999;
    float tmp_matchtrk_phi = -999;
    float tmp_matchtrk_z0 = -999;
    float tmp_matchtrk_d0 = -999;
    float tmp_matchtrk_chi2 = -999;
    float tmp_matchtrk_chi2_dof = -999;
    float tmp_matchtrk_chi2rphi = -999;
    float tmp_matchtrk_chi2rphi_dof = -999;
    float tmp_matchtrk_chi2rz = -999;
    float tmp_matchtrk_chi2rz_dof = -999;
    float tmp_matchtrk_bendchi2 = -999;
    float tmp_matchtrk_MVA1 = -999;
    int tmp_matchtrk_nstub = -999;
    int tmp_matchtrk_dhits = -999;
    int tmp_matchtrk_lhits = -999;
    int tmp_matchtrk_seed = -999;
    int tmp_matchtrk_hitpattern = -999;

    if (nMatch > 1 && DebugMode)
      edm::LogVerbatim("Tracklet") << "WARNING *** 2 or more matches to genuine L1 tracks ***";

    if (nMatch > 0) {
      tmp_matchtrk_pt = matchedTracks.at(i_track)->momentum().perp();
      tmp_matchtrk_eta = matchedTracks.at(i_track)->momentum().eta();
      tmp_matchtrk_phi = matchedTracks.at(i_track)->momentum().phi();
      tmp_matchtrk_z0 = matchedTracks.at(i_track)->z0();

      if (L1Tk_nPar == 5) {
        float tmp_matchtrk_x0 = matchedTracks.at(i_track)->POCA().x();
        float tmp_matchtrk_y0 = matchedTracks.at(i_track)->POCA().y();
        tmp_matchtrk_d0 = tmp_matchtrk_x0 * sin(tmp_matchtrk_phi) - tmp_matchtrk_y0 * cos(tmp_matchtrk_phi);
      }

      tmp_matchtrk_chi2 = matchedTracks.at(i_track)->chi2();
      tmp_matchtrk_chi2rphi = matchedTracks.at(i_track)->chi2XY();
      tmp_matchtrk_chi2rz = matchedTracks.at(i_track)->chi2Z();
      tmp_matchtrk_bendchi2 = matchedTracks.at(i_track)->stubPtConsistency();
      tmp_matchtrk_MVA1 = matchedTracks.at(i_track)->trkMVA1();
      tmp_matchtrk_nstub = (int)matchedTracks.at(i_track)->getStubRefs().size();
      tmp_matchtrk_seed = (int)matchedTracks.at(i_track)->trackSeedType();
      tmp_matchtrk_hitpattern = (int)matchedTracks.at(i_track)->hitPattern();

      int ndof = 2 * tmp_matchtrk_nstub - L1Tk_nPar;
      int ndofrphi = tmp_matchtrk_nstub - L1Tk_nPar + 2;
      int ndofrz = tmp_matchtrk_nstub - 2;
      tmp_matchtrk_chi2_dof = (float)tmp_matchtrk_chi2 / ndof;
      tmp_matchtrk_chi2rphi_dof = (float)tmp_matchtrk_chi2rphi / ndofrphi;
      tmp_matchtrk_chi2rz_dof = (float)tmp_matchtrk_chi2rz / ndofrz;

      // ------------------------------------------------------------------------------------------

      //float tmp_matchtrk_bend_chi2 = 0;

      tmp_matchtrk_dhits = 0;
      tmp_matchtrk_lhits = 0;

      std::vector<edm::Ref<edmNew::DetSetVector<TTStub<Ref_Phase2TrackerDigi_> >, TTStub<Ref_Phase2TrackerDigi_> > >
          stubRefs = matchedTracks.at(i_track)->getStubRefs();
      int tmp_nstub = stubRefs.size();

      for (int is = 0; is < tmp_nstub; is++) {
        DetId detIdStub = theTrackerGeom->idToDet((stubRefs.at(is)->clusterRef(0))->getDetId())->geographicalId();
        /*
        MeasurementPoint coords = stubRefs.at(is)->clusterRef(0)->findAverageLocalCoordinatesCentered();
        const GeomDet* theGeomDet = theTrackerGeom->idToDet(detIdStub);
        Global3DPoint posStub = theGeomDet->surface().toGlobal( theGeomDet->topology().localPosition(coords) );
        */

        int layer = -999999;
        if (detIdStub.subdetId() == StripSubdetector::TOB) {
          layer = static_cast<int>(tTopo->layer(detIdStub));
          tmp_matchtrk_lhits += pow(10, layer - 1);
        } else if (detIdStub.subdetId() == StripSubdetector::TID) {
          layer = static_cast<int>(tTopo->layer(detIdStub));
          tmp_matchtrk_dhits += pow(10, layer - 1);
        }

        // ------------------------------------------------------------------------------------------
      }
    }

    m_tp_pt->push_back(tmp_tp_pt);
    m_tp_eta->push_back(tmp_tp_eta);
    m_tp_phi->push_back(tmp_tp_phi);
    m_tp_dxy->push_back(tmp_tp_dxy);
    m_tp_z0->push_back(tmp_tp_z0);
    m_tp_d0->push_back(tmp_tp_d0);
    m_tp_z0_prod->push_back(tmp_tp_z0_prod);
    m_tp_d0_prod->push_back(tmp_tp_d0_prod);
    m_tp_pdgid->push_back(tmp_tp_pdgid);
    m_tp_nmatch->push_back(nMatch);
    m_tp_nstub->push_back(nStubTP);
    m_tp_eventid->push_back(tmp_eventid);
    m_tp_charge->push_back(tmp_tp_charge);

    m_matchtrk_pt->push_back(tmp_matchtrk_pt);
    m_matchtrk_eta->push_back(tmp_matchtrk_eta);
    m_matchtrk_phi->push_back(tmp_matchtrk_phi);
    m_matchtrk_z0->push_back(tmp_matchtrk_z0);
    m_matchtrk_d0->push_back(tmp_matchtrk_d0);
    m_matchtrk_chi2->push_back(tmp_matchtrk_chi2);
    m_matchtrk_chi2rphi->push_back(tmp_matchtrk_chi2rphi);
    m_matchtrk_chi2rz->push_back(tmp_matchtrk_chi2rz);
    m_matchtrk_bendchi2->push_back(tmp_matchtrk_bendchi2);
    m_matchtrk_MVA1->push_back(tmp_matchtrk_MVA1);
    m_matchtrk_nstub->push_back(tmp_matchtrk_nstub);
    m_matchtrk_dhits->push_back(tmp_matchtrk_dhits);
    m_matchtrk_lhits->push_back(tmp_matchtrk_lhits);
    m_matchtrk_seed->push_back(tmp_matchtrk_seed);
    m_matchtrk_hitpattern->push_back(tmp_matchtrk_hitpattern);
    m_matchtrk_chi2_dof->push_back(tmp_matchtrk_chi2_dof);
    m_matchtrk_chi2rphi_dof->push_back(tmp_matchtrk_chi2rphi_dof);
    m_matchtrk_chi2rz_dof->push_back(tmp_matchtrk_chi2rz_dof);

  }  //end loop tracking particles

//   if (TrackingInJets) {
//     for (int ij = 0; ij < (int)v_jets.size(); ij++) {
//       if (ij < NJETS) {
//         m_jet_eta->push_back((v_jets.at(ij)).eta());
//         m_jet_phi->push_back((v_jets.at(ij)).phi());
//         m_jet_pt->push_back((v_jets.at(ij)).pt());
//         m_jet_tp_sumpt->push_back(jets_tp_sumpt[ij]);
//         m_jet_trk_sumpt->push_back(jets_trk_sumpt[ij]);
//         m_jet_matchtrk_sumpt->push_back(jets_matchtrk_sumpt[ij]);
//       }
//     }
//   }
// }  // end of produce()

// FILLDESCRIPTIONS
void L1SmartPixelsTrackProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //L1SmartPixelsTrackProducer
  //TODO: modify from L1TrackSelectionProducer.cc plugin
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("l1TracksInputTag", edm::InputTag("TTTracksFromTrackletEmulation", "Level1TTTracks"));
  desc.add<std::string>("outputCollectionName", "Level1TTTracksSelected");
  {
    edm::ParameterSetDescription descCutSet;
    descCutSet.add<double>("ptMin", 2.0)->setComment("pt must be greater than this value, [GeV]");
    descCutSet.add<double>("absEtaMax", 2.4)->setComment("absolute value of eta must be less than this value");
    descCutSet.add<double>("absZ0Max", 15.0)->setComment("z0 must be less than this value, [cm]");
    descCutSet.add<int>("nStubsMin", 4)->setComment("number of stubs must be greater than or equal to this value");
    descCutSet.add<int>("nPSStubsMin", 0)
        ->setComment("number of stubs in the PS Modules must be greater than or equal to this value");

    descCutSet.add<double>("promptMVAMin", -1.0)->setComment("MVA must be greater than this value");
    descCutSet.add<double>("reducedBendChi2Max", 2.25)->setComment("bend chi2 must be less than this value");
    descCutSet.add<double>("reducedChi2RZMax", 5.0)->setComment("chi2rz/dof must be less than this value");
    descCutSet.add<double>("reducedChi2RPhiMax", 20.0)->setComment("chi2rphi/dof must be less than this value");
    descCutSet.add<double>("reducedChi2RZMaxNstub4", 999.9)
        ->setComment("chi2rz/dof must be less than this value in nstub==4");
    descCutSet.add<double>("reducedChi2RZMaxNstub5", 999.9)
        ->setComment("chi2rz/dof must be less than this value in nstub>4");
    descCutSet.add<double>("reducedChi2RPhiMaxNstub4", 999.9)
        ->setComment("chi2rphi/dof must be less than this value in nstub==4");
    descCutSet.add<double>("reducedChi2RPhiMaxNstub5", 999.9)
        ->setComment("chi2rphi/dof must be less than this value in nstub>4");
    descCutSet.add<double>("reducedBendChi2MaxNstub4", 999.9)
        ->setComment("bend chi2 must be less than this value in nstub==4");
    descCutSet.add<double>("reducedBendChi2MaxNstub5", 999.9)
        ->setComment("bend chi2 must be less than this value in nstub>4");

    desc.add<edm::ParameterSetDescription>("cutSet", descCutSet);
  }
  desc.add<bool>("processSimulatedTracks", true)
      ->setComment("return selected tracks after cutting on the floating point values");
  desc.add<bool>("processEmulatedTracks", true)
      ->setComment("return selected tracks after cutting on the bitwise emulated values");
  desc.add<int>("debug", 0)->setComment("Verbosity levels: 0, 1, 2, 3");
  descriptions.addWithDefaultLabel(desc);
}
///////////////////////////
// DEFINE THIS AS A PLUG-IN
DEFINE_FWK_MODULE(L1SmartPixelsTrackProducer);
