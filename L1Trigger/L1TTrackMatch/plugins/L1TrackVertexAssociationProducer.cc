// -*- C++ -*-
//
// Package:    L1Trigger/L1TTrackMatch
// Class:      L1TrackVertexAssociationProducer
//
/**\class L1TrackVertexAssociationProducer L1TrackVertexAssociationProducer.cc L1Trigger/L1TTrackMatch/plugins/L1TrackVertexAssociationProducer.cc

 Description: Selects a set of L1Tracks based on a set of predefined criteria.

 Implementation:
     Inputs:
         std::vector<TTTrack> - Each floating point TTTrack inside this collection inherits from
                                a bit-accurate TTTrack_TrackWord, used for emulation purposes.
     Outputs:
         std::vector<TTTrack> - A collection of TTTracks selected from cuts on the TTTrack properties
         std::vector<TTTrack> - A collection of TTTracks selected from cuts on the TTTrack_TrackWord properties
*/
//
// Original Author:  Alexx Perloff
//         Created:  Thu, 16 Dec 2021 19:02:50 GMT
// Derivative Author: Nick Manganelli
//         Created: Wed, 1 Feb 2023 14:23:56 GMT
//
//

// system include files
#include <algorithm>
#include <memory>
#include <string>
#include <vector>

// Xilinx HLS includes
#include <ap_fixed.h>
#include <ap_int.h>

// user include files
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefVector.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/L1Trigger/interface/Vertex.h"
#include "DataFormats/L1Trigger/interface/VertexWord.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "CommonTools/Utils/interface/AndSelector.h"
#include "CommonTools/Utils/interface/EtaRangeSelector.h"
#include "CommonTools/Utils/interface/MinSelector.h"
#include "CommonTools/Utils/interface/MinFunctionSelector.h"
#include "CommonTools/Utils/interface/MinNumberSelector.h"
#include "CommonTools/Utils/interface/PtMinSelector.h"
#include "CommonTools/Utils/interface/Selection.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"

//
// class declaration
//

class L1TrackVertexAssociationProducer : public edm::global::EDProducer<> {
public:
  explicit L1TrackVertexAssociationProducer(const edm::ParameterSet&);
  ~L1TrackVertexAssociationProducer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  // ----------constants, enums and typedefs ---------
  // Relevant constants for the converted track word
  enum TrackBitWidths {
    kPtSize = TTTrack_TrackWord::TrackBitWidths::kRinvSize - 1,  // Width of pt
    kPtMagSize = 9,                                              // Width of pt magnitude (unsigned)
    kEtaSize = TTTrack_TrackWord::TrackBitWidths::kTanlSize,     // Width of eta
    kEtaMagSize = 3,                                             // Width of eta magnitude (signed)
  };

  typedef TTTrack<Ref_Phase2TrackerDigi_> L1Track;
  typedef std::vector<L1Track> TTTrackCollection;
  typedef edm::Handle<TTTrackCollection> TTTrackCollectionHandle;
  typedef edm::Ref<TTTrackCollection> TTTrackRef;
  typedef edm::RefVector<TTTrackCollection> TTTrackRefCollection;
  typedef std::unique_ptr<TTTrackRefCollection> TTTrackRefCollectionUPtr;

  // ----------member functions ----------------------
  void printDebugInfo(const TTTrackCollectionHandle& l1SelectedTracksHandle, 
		      const TTTrackCollectionHandle& l1SelectedTracksEmulationHandle, 
                      const TTTrackRefCollectionUPtr& vTTTrackAssociatedOutput,
                      const TTTrackRefCollectionUPtr& vTTTrackAssociatedEmulationOutput) const;
  void printTrackInfo(edm::LogInfo& log, const L1Track& track, bool printEmulation = false) const;
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  // ----------selectors -----------------------------
  // Based on recommendations from https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideGenericSelectors
  struct TTTrackWordAbsZ0MaxSelector {
    TTTrackWordAbsZ0MaxSelector(double absZ0Max) : absZ0Max_(absZ0Max) {}
    TTTrackWordAbsZ0MaxSelector(const edm::ParameterSet& cfg)
        : absZ0Max_(cfg.template getParameter<double>("absZ0Max")) {}
    bool operator()(const L1Track& t) const {
      double floatZ0 = t.undigitizeSignedValue(
          t.getZ0Bits(), TTTrack_TrackWord::TrackBitWidths::kZ0Size, TTTrack_TrackWord::stepZ0, 0.0);
      return std::abs(floatZ0) <= absZ0Max_;
    }

  private:
    double absZ0Max_;
  };

  struct TTTrackDeltaZMaxSelector {
    TTTrackDeltaZMaxSelector(const std::vector<double>& deltaZMaxEtaBounds, const std::vector<double>& deltaZMax)
        : deltaZMaxEtaBounds_(deltaZMaxEtaBounds), deltaZMax_(deltaZMax) {}
    TTTrackDeltaZMaxSelector(const edm::ParameterSet& cfg)
        : deltaZMaxEtaBounds_(cfg.template getParameter<double>("deltaZMaxEtaBounds")),
          deltaZMax_(cfg.template getParameter<double>("deltaZMax")) {}
    bool operator()(const L1Track& t, const l1t::Vertex& v) const {
      size_t etaIndex =
          std::upper_bound(deltaZMaxEtaBounds_.begin(), deltaZMaxEtaBounds_.end(), std::abs(t.momentum().eta())) -
          deltaZMaxEtaBounds_.begin() - 1;
      if (etaIndex > deltaZMax_.size() - 1)
        etaIndex = deltaZMax_.size() - 1;
      return std::abs(v.z0() - t.z0()) <= deltaZMax_[etaIndex];
    }

  private:
    std::vector<double> deltaZMaxEtaBounds_;
    std::vector<double> deltaZMax_;
  };
  struct TTTrackWordDeltaZMaxSelector {
    TTTrackWordDeltaZMaxSelector(const std::vector<double>& deltaZMaxEtaBounds, const std::vector<double>& deltaZMax)
        : deltaZMaxEtaBounds_(deltaZMaxEtaBounds), deltaZMax_(deltaZMax) {}
    TTTrackWordDeltaZMaxSelector(const edm::ParameterSet& cfg)
        : deltaZMaxEtaBounds_(cfg.template getParameter<double>("deltaZMaxEtaBounds")),
          deltaZMax_(cfg.template getParameter<double>("deltaZMax")) {}
    bool operator()(const L1Track& t, const l1t::VertexWord& v) const {
      TTTrack_TrackWord::tanl_t etaEmulationBits = t.getTanlWord();
      ap_fixed<TrackBitWidths::kEtaSize, TrackBitWidths::kEtaMagSize> etaEmulation;
      etaEmulation.V = etaEmulationBits.range();
      size_t etaIndex =
          std::upper_bound(deltaZMaxEtaBounds_.begin(), deltaZMaxEtaBounds_.end(), std::abs(etaEmulation.to_double())) -
          deltaZMaxEtaBounds_.begin() - 1;
      if (etaIndex > deltaZMax_.size() - 1)
        etaIndex = deltaZMax_.size() - 1;
      l1t::VertexWord::vtxz0_t fixedTkZ0 = t.undigitizeSignedValue(
          t.getZ0Bits(), TTTrack_TrackWord::TrackBitWidths::kZ0Size, TTTrack_TrackWord::stepZ0, 0.0);

      ap_uint<TrackBitWidths::kPtSize> ptEmulationBits = t.getTrackWord()(
          TTTrack_TrackWord::TrackBitLocations::kRinvMSB - 1, TTTrack_TrackWord::TrackBitLocations::kRinvLSB);
      ap_ufixed<TrackBitWidths::kPtSize, TrackBitWidths::kPtMagSize> ptEmulation;
      ptEmulation.V = ptEmulationBits.range();
      return std::abs(v.z0() - fixedTkZ0.to_double()) <= deltaZMax_[etaIndex];
    }

  private:
    std::vector<double> deltaZMaxEtaBounds_;
    std::vector<double> deltaZMax_;
  };

  // ----------member data ---------------------------
  edm::EDGetTokenT<TTTrackCollection> l1SelectedTracksToken_;
  edm::EDGetTokenT<TTTrackCollection> l1SelectedTracksEmulationToken_;
  edm::EDGetTokenT<l1t::VertexCollection> l1VerticesToken_;
  edm::EDGetTokenT<l1t::VertexWordCollection> l1VerticesEmulationToken_;
  const std::string outputCollectionName_;
  const edm::ParameterSet cutSet_;
  std::vector<double> deltaZMaxEtaBounds_, deltaZMax_;
  const double useDisplacedTracksDeltaZOverride_;
  bool processSimulatedTracks_, processEmulatedTracks_, doDeltaZCutSim_, doDeltaZCutEmu_;
  int debug_;
};

//
// constructors and destructor
//
L1TrackVertexAssociationProducer::L1TrackVertexAssociationProducer(const edm::ParameterSet& iConfig)
    : outputCollectionName_(iConfig.getParameter<std::string>("outputCollectionName")),
      cutSet_(iConfig.getParameter<edm::ParameterSet>("cutSet")),
      deltaZMaxEtaBounds_(cutSet_.getParameter<std::vector<double>>("deltaZMaxEtaBounds")),
      deltaZMax_(cutSet_.getParameter<std::vector<double>>("deltaZMax")),

      useDisplacedTracksDeltaZOverride_(iConfig.getParameter<double>("useDisplacedTracksDeltaZOverride")),
      processSimulatedTracks_(iConfig.getParameter<bool>("processSimulatedTracks")),
      processEmulatedTracks_(iConfig.getParameter<bool>("processEmulatedTracks")),
      debug_(iConfig.getParameter<int>("debug")) {
  // Confirm the the configuration makes sense
  if (!processSimulatedTracks_ && !processEmulatedTracks_) {
    throw cms::Exception("You must process at least one of the track collections (simulated or emulated).");
  }

  if (deltaZMax_.size() != deltaZMaxEtaBounds_.size() - 1) {
    throw cms::Exception("The number of deltaZ cuts does not match the number of eta bins!");
  }

  if (useDisplacedTracksDeltaZOverride_ >= 0) {
    deltaZMax_ = std::vector<double>(deltaZMax_.size(), useDisplacedTracksDeltaZOverride_);
  }

  // Get additional input tags and define the EDM output based on the previous configuration parameters
  doDeltaZCutSim_ = false;
  doDeltaZCutEmu_ = false;
  if (processSimulatedTracks_) {
    if (iConfig.exists("l1VerticesInputTag")) {
      l1SelectedTracksToken_ = consumes<TTTrackCollection>(iConfig.getParameter<edm::InputTag>("l1SelectedTracksInputTag"));
      l1VerticesToken_ = consumes<l1t::VertexCollection>(iConfig.getParameter<edm::InputTag>("l1VerticesInputTag"));
      doDeltaZCutSim_ = true;
      produces<TTTrackRefCollection>(outputCollectionName_);
    }
  }
  if (processEmulatedTracks_) {
    if (iConfig.exists("l1VerticesEmulationInputTag")) {
      l1SelectedTracksEmulationToken_ = consumes<TTTrackCollection>(iConfig.getParameter<edm::InputTag>("l1SelectedTracksEmulationInputTag"));
      l1VerticesEmulationToken_ =
	consumes<l1t::VertexWordCollection>(iConfig.getParameter<edm::InputTag>("l1VerticesEmulationInputTag"));
      doDeltaZCutEmu_ = true;
      produces<TTTrackRefCollection>(outputCollectionName_ + "Emulation");
    }
  }
}

L1TrackVertexAssociationProducer::~L1TrackVertexAssociationProducer() {}

//
// member functions
//

void L1TrackVertexAssociationProducer::printDebugInfo(const TTTrackCollectionHandle& l1SelectedTracksHandle, 
						      const TTTrackCollectionHandle& l1SelectedTracksEmulationHandle, 
						      const TTTrackRefCollectionUPtr& vTTTrackAssociatedOutput,
						      const TTTrackRefCollectionUPtr& vTTTrackAssociatedEmulationOutput
						      ) const {
  edm::LogInfo log("L1TrackVertexAssociationProducer");

  if (processSimulatedTracks_) {
    log << "\t---\n\tNumber of tracks in this selection = " << l1SelectedTracksHandle->size() << "\n\n";
    log << "The vertex associated track collection (pt, eta, phi, nstub, bendchi2, chi2rz, chi2rphi, z0) values are ... \n";
    for (const auto& track : *vTTTrackAssociatedOutput) {
      printTrackInfo(log, *track, debug_ >= 4);
    }
    log << "\t---\n\tNumber of tracks in this association = " << vTTTrackAssociatedOutput->size() << "\n\n";
  }

  if (processEmulatedTracks_) {
    log << "\t---\n\tNumber of emulated tracks in this selection = " << l1SelectedTracksEmulationHandle->size() << "\n\n";
    log << "The emulation vertex associated track collection (pt, eta, phi, nstub, bendchi2, chi2rz, chi2rphi, z0) values are "
           "... \n";
    for (const auto& track : *vTTTrackAssociatedEmulationOutput) {
      printTrackInfo(log, *track, debug_ >= 4);
    }
    log << "\t---\n\tNumber of emulated tracks in this association = " << vTTTrackAssociatedEmulationOutput->size() << "\n\n";
  }

  if (processSimulatedTracks_ && processEmulatedTracks_) {
    TTTrackRefCollection inSimButNotEmu;
    TTTrackRefCollection inEmuButNotSim;
    std::set_difference(vTTTrackAssociatedOutput->begin(),
                        vTTTrackAssociatedOutput->end(),
                        vTTTrackAssociatedEmulationOutput->begin(),
                        vTTTrackAssociatedEmulationOutput->end(),
                        std::back_inserter(inSimButNotEmu));
    std::set_difference(vTTTrackAssociatedEmulationOutput->begin(),
                        vTTTrackAssociatedEmulationOutput->end(),
                        vTTTrackAssociatedOutput->begin(),
                        vTTTrackAssociatedOutput->end(),
                        std::back_inserter(inEmuButNotSim));
    log << "The set of vertex associated tracks selected via cuts on the simulated values which are not in the set of tracks selected "
           "by cutting on the emulated values ... \n";
    for (const auto& track : inSimButNotEmu) {
      printTrackInfo(log, *track, debug_ >= 3);
    }
    log << "\t---\n\tNumber of tracks in this selection = " << inSimButNotEmu.size() << "\n\n"
        << "The set of vertex associated tracks selected via cuts on the emulated values which are not in the set of tracks selected "
           "by cutting on the simulated values ... \n";
    for (const auto& track : inEmuButNotSim) {
      printTrackInfo(log, *track, debug_ >= 3);
    }
    log << "\t---\n\tNumber of tracks in this selection = " << inEmuButNotSim.size() << "\n\n";
  }
}

void L1TrackVertexAssociationProducer::printTrackInfo(edm::LogInfo& log, const L1Track& track, bool printEmulation) const {
  log << "\t(" << track.momentum().perp() << ", " << track.momentum().eta() << ", " << track.momentum().phi() << ", "
      << track.getStubRefs().size() << ", " << track.stubPtConsistency() << ", " << track.chi2ZRed() << ", "
      << track.chi2XYRed() << ", " << track.z0() << ")\n";

  if (printEmulation) {
    ap_uint<TrackBitWidths::kPtSize> ptEmulationBits = track.getTrackWord()(
        TTTrack_TrackWord::TrackBitLocations::kRinvMSB - 1, TTTrack_TrackWord::TrackBitLocations::kRinvLSB);
    ap_ufixed<TrackBitWidths::kPtSize, TrackBitWidths::kPtMagSize> ptEmulation;
    ptEmulation.V = ptEmulationBits.range();
    TTTrack_TrackWord::tanl_t etaEmulationBits = track.getTanlWord();
    ap_fixed<TrackBitWidths::kEtaSize, TrackBitWidths::kEtaMagSize> etaEmulation;
    etaEmulation.V = etaEmulationBits.range();
    double floatTkZ0 = track.undigitizeSignedValue(
        track.getZ0Bits(), TTTrack_TrackWord::TrackBitWidths::kZ0Size, TTTrack_TrackWord::stepZ0, 0.0);
    double floatTkPhi = track.undigitizeSignedValue(
        track.getPhiBits(), TTTrack_TrackWord::TrackBitWidths::kPhiSize, TTTrack_TrackWord::stepPhi0, 0.0);
    log << "\t\t(" << ptEmulation.to_double() << ", " << etaEmulation.to_double() << ", " << floatTkPhi << ", "
        << track.getNStubs() << ", " << track.getBendChi2() << ", " << track.getChi2RZ() << ", " << track.getChi2RPhi()
        << ", " << floatTkZ0 << ")\n";
  }
}

// ------------ method called to produce the data  ------------
void L1TrackVertexAssociationProducer::produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const {
  // auto vTTTrackOutput = std::make_unique<TTTrackRefCollection>(); //FIXME:REMOVE
  auto vTTTrackAssociatedOutput = std::make_unique<TTTrackRefCollection>();
  // auto vTTTrackEmulationOutput = std::make_unique<TTTrackRefCollection>(); //FIXME:REMOVE
  auto vTTTrackAssociatedEmulationOutput = std::make_unique<TTTrackRefCollection>();

  TTTrackCollectionHandle l1SelectedTracksHandle;
  TTTrackCollectionHandle l1SelectedTracksEmulationHandle;
  edm::Handle<l1t::VertexCollection> l1VerticesHandle;
  edm::Handle<l1t::VertexWordCollection> l1VerticesEmulationHandle;

  l1t::Vertex leadingVertex;
  l1t::VertexWord leadingEmulationVertex;

  // TTTrackPtMinEtaMaxZ0MaxNStubsMinSelector kinSel(ptMin_, absEtaMax_, absZ0Max_, nStubsMin_); //FIXME:REMOVE
  // TTTrackWordPtMinEtaMaxZ0MaxNStubsMinSelector kinSelEmu(ptMin_, absEtaMax_, absZ0Max_, nStubsMin_); //FIXME:REMOVE
  // TTTrackBendChi2Chi2RZChi2RPhiMaxSelector chi2Sel(bendChi2Max_, reducedChi2RZMax_, reducedChi2RPhiMax_); //FIXME:REMOVE
  // TTTrackWordBendChi2Chi2RZChi2RPhiMaxSelector chi2SelEmu(bendChi2Max_, reducedChi2RZMax_, reducedChi2RPhiMax_); //FIXME:REMOVE
  TTTrackDeltaZMaxSelector deltaZSel(deltaZMaxEtaBounds_, deltaZMax_);
  TTTrackWordDeltaZMaxSelector deltaZSelEmu(deltaZMaxEtaBounds_, deltaZMax_);

  size_t nOutputApproximate = 0;
  if (processSimulatedTracks_) {
    iEvent.getByToken(l1SelectedTracksToken_, l1SelectedTracksHandle);
    nOutputApproximate = l1SelectedTracksHandle->size();
    if (doDeltaZCutSim_) {
      iEvent.getByToken(l1VerticesToken_, l1VerticesHandle);
      leadingVertex = l1VerticesHandle->at(0);
      if (debug_ >= 2) {
        edm::LogInfo("L1TrackVertexAssociationProducer") << "leading vertex z0 = " << leadingVertex.z0();
      }
    }
    // vTTTrackOutput->reserve(nOutputApproximate); //FIXME:REMOVE
    vTTTrackAssociatedOutput->reserve(nOutputApproximate);
  

    for (size_t i = 0; i < nOutputApproximate; i++) {
      const auto& track = l1SelectedTracksHandle->at(i);
      // Select tracks based on the floating point TTTrack
      if (processSimulatedTracks_) { //  && kinSel(track) && nPSStubsSel(track) && chi2Sel(track)) { //FIXME:REMOVE
        // vTTTrackOutput->push_back(TTTrackRef(l1TracksHandle, i)); //FIXME:REMOVE
        // if (doDeltaZCutSim_ && deltaZSel(track, leadingVertex)) { //FIXME:REMOVE
        if (deltaZSel(track, leadingVertex)) {
          vTTTrackAssociatedOutput->push_back(TTTrackRef(l1SelectedTracksHandle, i));
        }
      }
    }
    // Put the outputs into the event
    iEvent.put(std::move(vTTTrackAssociatedOutput), outputCollectionName_);
  }

  size_t nOutputEmulationApproximate = 0;
  if (processEmulatedTracks_) {
    iEvent.getByToken(l1SelectedTracksEmulationToken_, l1SelectedTracksEmulationHandle);
    nOutputEmulationApproximate = l1SelectedTracksEmulationHandle->size();
    if (doDeltaZCutEmu_) {
      iEvent.getByToken(l1VerticesEmulationToken_, l1VerticesEmulationHandle);
      leadingEmulationVertex = l1VerticesEmulationHandle->at(0);
      if (debug_ >= 2) {
        edm::LogInfo("L1TrackVertexAssociationProducer") << "leading emulation vertex z0 = " << leadingEmulationVertex.z0();
      }
    }
    // vTTTrackEmulationOutput->reserve(nOutputApproximate); //FIXME:REMOVE
    vTTTrackAssociatedEmulationOutput->reserve(nOutputEmulationApproximate);
    for (size_t i = 0; i < nOutputEmulationApproximate; i++) {
      const auto& track = l1SelectedTracksEmulationHandle->at(i);
      // Select tracks based on the bitwise accurate TTTrack_TrackWord
      if (processEmulatedTracks_) { // && kinSelEmu(track) && chi2SelEmu(track)) { //FIXME:REMOVE
	// vTTTrackEmulationOutput->push_back(TTTrackRef(l1TracksHandle, i)); //FIXME:REMOVE
	// if (doDeltaZCutEmu_ && deltaZSelEmu(track, leadingEmulationVertex)) { //FIXME:REMOVE
	if (deltaZSelEmu(track, leadingEmulationVertex)) {
	  vTTTrackAssociatedEmulationOutput->push_back(TTTrackRef(l1SelectedTracksEmulationHandle, i));
	}
      }
    }
    // Put the outputs into the event
    iEvent.put(std::move(vTTTrackAssociatedEmulationOutput), outputCollectionName_ + "Emulation");
  }

  if (debug_ >= 2 && processSimulatedTracks_ && processEmulatedTracks_) {
    printDebugInfo(l1SelectedTracksHandle,
		   l1SelectedTracksEmulationHandle,
                   vTTTrackAssociatedOutput,
                   vTTTrackAssociatedEmulationOutput);
  }
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void L1TrackVertexAssociationProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //L1TrackVertexAssociationProducer
  edm::ParameterSetDescription desc;
  desc.addOptional<edm::InputTag>("l1SelectedTracksInputTag", edm::InputTag("L1TrackSelectionProducer", "Level1TTTracksSelected"));
  desc.addOptional<edm::InputTag>("l1SelectedTracksEmulationInputTag", 
				  edm::InputTag("L1TrackSelectionProducer", "Level1TTTracksSelectedEmulation"));
  desc.addOptional<edm::InputTag>("l1VerticesInputTag", edm::InputTag("L1VertexFinder", "l1vertices"));
  desc.addOptional<edm::InputTag>("l1VerticesEmulationInputTag",
                                  edm::InputTag("L1VertexFinderEmulator", "l1verticesEmulation"));
  desc.add<std::string>("outputCollectionName", "Level1TTTracksSelectedAssociated");
  {
    edm::ParameterSetDescription descCutSet;
    descCutSet.add<std::vector<double>>("deltaZMaxEtaBounds", {0.0, 0.7, 1.0, 1.2, 1.6, 2.0, 2.4})
        ->setComment("these values define the bin boundaries in |eta|");
    descCutSet.add<std::vector<double>>("deltaZMax", {0.37, 0.50, 0.60, 0.75, 1.00, 1.60})
        ->setComment(
            "delta z must be less than these values, there will be one less value here than in deltaZMaxEtaBounds, "
            "[cm]");
    desc.add<edm::ParameterSetDescription>("cutSet", descCutSet);
  }
  desc.add<double>("useDisplacedTracksDeltaZOverride", -1.0)
      ->setComment("override the deltaZ cut value for displaced tracks");
  desc.add<bool>("processSimulatedTracks", true)
      ->setComment("return selected tracks after cutting on the floating point values");
  desc.add<bool>("processEmulatedTracks", true)
      ->setComment("return selected tracks after cutting on the bitwise emulated values");
  desc.add<int>("debug", 0)->setComment("Verbosity levels: 0, 1, 2, 3");
  descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1TrackVertexAssociationProducer);
