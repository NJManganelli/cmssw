// -*- C++ -*-
//
// Package:    L1Trigger/Phase3SmartPixels
// Class:      L1SmartPixelsTrackProducer
//
/**\class L1SmartPixelsTrackProducer L1SmartPixelsTrackProducer.cc L1Trigger/Phase3SmartPixels/plugins/L1SmartPixelsTrackProducer.cc

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

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/Exception.h"
// digiRefit: shared projector + pixel digi/simlink formats. The angle-synthesis
// throw is FACTORIZED OUT of compiled code into the payload (HashPRNG fused shift
// compounds, doc/PixelAVAngleResponseSpec.md §3): synthesized angles are pure
// correctionlib evaluates, deterministic per input tuple, no engine involved.
// A local CLHEP engine remains ONLY for the optional Stack-B noise-angle
// inverse-CDF draws (smarthitFakeSet quantile inputs) -- that payload's schema is
// out of scope here; the engine is seeded from hash(label,run,lumi,event) so it
// stays event-order-independent and split-job invariant, but its draw sequence is
// code-path-ordered (see doc/Phase2Acceptance.md §1, doc/RefitSidecarSpec.md).
// The RandomNumberGeneratorService stream engine is deliberately NOT used.
#include "CLHEP/Random/MixMaxRng.h"
#include "CLHEP/Random/RandFlat.h"
#include "L1Trigger/Phase3SmartPixels/interface/SmartPixelsHelixProjector.h"
#include "L1Trigger/Phase3SmartPixels/interface/SmartPixelsParentMap.h"
#include "L1Trigger/Phase3SmartPixels/interface/SmartPixelsRefitSidecar.h"
#include "L1Trigger/Phase3SmartPixels/interface/SmartPixelsTransmittedSubset.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "Geometry/CommonTopologies/interface/PixelGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/MeasurementPoint.h"
#include "SimDataFormats/TrackerDigiSimLink/interface/PixelDigiSimLink.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"

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
#include "SimDataFormats/Associations/interface/TTClusterAssociationMap.h"
#include "SimDataFormats/Associations/interface/TTStubAssociationMap.h"
#include "SimDataFormats/Associations/interface/TTTrackAssociationMap.h"
#include "Geometry/Records/interface/StackedTrackerGeometryRecord.h"

// #include "DataFormats/JetReco/interface/GenJetCollection.h"
// #include "DataFormats/JetReco/interface/GenJet.h"

////////////////////////////
// DETECTOR GEOMETRY HEADERS
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
// #include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
// #include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
// #include "Geometry/TrackerGeometryBuilder/interface/RectangularPixelTopology.h"
// #include "Geometry/CommonTopologies/interface/GeomDetType.h"
// #include "Geometry/CommonTopologies/interface/GeomDet.h"

// #include "Geometry/CommonTopologies/interface/PixelGeomDetUnit.h"
// #include "Geometry/CommonTopologies/interface/PixelGeomDetType.h"
// #include "Geometry/TrackerGeometryBuilder/interface/PixelTopologyBuilder.h"
// #include "Geometry/Records/interface/StackedTrackerGeometryRecord.h"

////////////////
// PHYSICS TOOLS
#include "L1Trigger/TrackFindingTracklet/interface/HitPatternHelper.h"
// Branch C (covMatrix omnibus): hph::Setup is keyed on trackerDTC::SetupRcd
// (EVENTSETUP_DATA_DEFAULT_RECORD in HitPatternHelper.h). No HitPatternHelperRcd.h
// in the #51503/L1TK-dev-15_1 track layer.
#include "L1Trigger/TrackerDTC/interface/SetupRcd.h"
// #include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include "L1Trigger/TrackTrigger/interface/StubPtConsistency.h"

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

// Correctionlib headers
#include "correction.h"

//////////////
// STD HEADERS
#include <map>
#include <set>
#include <memory>
#include <string>
#include <iostream>
#include <cstdint>
#include <cmath>
#include <fstream>
#include <limits>

// Refit-quality BDT (spec §6a): conifer GBDT evaluated in-flight on the
// REFIT_BDT_FEATURES v0 feature vector. ap_fixed.h pulls the HLS ap types that
// conifer.h references; we instantiate the float/float emulation backend.
#include "ap_fixed.h"
#include "conifer.h"
#include "nlohmann/json.hpp"

//////////////
// NAMESPACES
using namespace std;
using namespace edm;

namespace {
  // Deterministic per-event seed for the digiRefit local RNG. FNV-1a 64-bit hash
  // over the module label chars, then run/lumi/event folded in the SAME 1a step
  // (xor-then-multiply) so the result depends on (label, run, lumi, event) only.
  // This makes the angle-synthesis draw sequence identical for a given physics
  // event regardless of file splitting, skipped events, or stream scheduling
  // (bitwise split-job reproducibility for training productions). Reduced to a
  // nonzero 31-bit CLHEP seed for CLHEP::MixMaxRng (good short-seed behavior).
  inline long digiRefitSeed(const std::string& label,
                            unsigned int run,
                            unsigned int lumi,
                            unsigned long long event) {
    constexpr uint64_t kFnvOffset = 1469598103934665603ULL;
    constexpr uint64_t kFnvPrime = 1099511628211ULL;
    uint64_t h = kFnvOffset;
    for (unsigned char c : label) {
      h ^= static_cast<uint64_t>(c);
      h *= kFnvPrime;
    }
    for (int shift = 0; shift < 64; shift += 8) {
      h ^= static_cast<uint64_t>((static_cast<uint64_t>(run) >> shift) & 0xFFu);
      h *= kFnvPrime;
    }
    for (int shift = 0; shift < 64; shift += 8) {
      h ^= static_cast<uint64_t>((static_cast<uint64_t>(lumi) >> shift) & 0xFFu);
      h *= kFnvPrime;
    }
    for (int shift = 0; shift < 64; shift += 8) {
      h ^= static_cast<uint64_t>((event >> shift) & 0xFFu);
      h *= kFnvPrime;
    }
    // Fold the high half into the low half, then mask to a nonzero 31-bit seed.
    uint32_t s = static_cast<uint32_t>((h ^ (h >> 32)) & 0x7FFFFFFFu);
    if (s == 0u)
      s = 1u;
    return static_cast<long>(s);
  }

  // REFIT_BDT_FEATURES v1 (spec §6a v1): decode the INPUT track's classic-7
  // TrackQuality hw features EXACTLY as ngtagger-train/train/trkquality.py does
  // (so the producer's in-flight v1 vector is bit-identical to the training
  // vector). twos_complement mirrors trkquality.twos_complement; nlaymissInterior
  // mirrors trkquality.nlaymiss_interior (width 7).
  inline double twosComplement(unsigned int bits, int width) {
    const long v = static_cast<long>(bits);
    const long half = (1L << (width - 1));
    return static_cast<double>(v >= half ? v - (1L << width) : v);
  }
  inline double nlaymissInterior(unsigned int hitpattern, int width = 7) {
    if (hitpattern == 0u)
      return 0.0;
    int first = -1, last = -1;
    for (int b = 0; b < width; ++b)
      if ((hitpattern >> b) & 1u) {
        if (first < 0)
          first = b;
        last = b;
      }
    int miss = 0;
    for (int b = first; b <= last; ++b)
      if (!((hitpattern >> b) & 1u))
        ++miss;
    return static_cast<double>(miss);
  }
}  // namespace

//////////////////////////////
//                          //
//     CLASS DEFINITION     //
//                          //
//////////////////////////////

class L1SmartPixelsTrackProducer : public edm::stream::EDProducer<edm::RunSummaryCache< std::map< std::string, std::map< int, double > > > > {
public:
  // Constructor/destructor
  explicit L1SmartPixelsTrackProducer(const edm::ParameterSet& iConfig);
  ~L1SmartPixelsTrackProducer() override;

  // Mandatory methods
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  void beginStream(edm::StreamID) override;
  void endStream() override;
  void produce(edm::Event&, const edm::EventSetup&) override;

  static std::unique_ptr< std::map< std::string, std::map< int, double > > > initializeGlobalCache(edm::ParameterSet const& iPSet) {
    return std::unique_ptr< std::map< std::string, std::map< int, double > > >( new std::map< std::string, std::map< int, double > >);
  }

  static std::shared_ptr< std::map< std::string, std::map< int, double > > > globalBeginRunSummary(edm::Run const&,
												   edm::EventSetup const&,
												   RunContext const*) {
    return std::make_unique< std::map< std::string, std::map< int, double > > >(
      std::map<std::string, std::map<int, double>>{
	{"npars_4_entries", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
	{"npars_5_entries", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
	{"resolution_entries", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
	{"float_pt_diff_sum_orig", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
	{"float_pt_diff_sum2_orig", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
	{"float_pt_diff_sum_new", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
	{"float_pt_diff_sum2_new", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
	{"float_eta_diff_sum_orig", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
	{"float_eta_diff_sum2_orig", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
	{"float_eta_diff_sum_new", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
	{"float_eta_diff_sum2_new", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
	{"float_phi_diff_sum_orig", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
	{"float_phi_diff_sum2_orig", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
	{"float_phi_diff_sum_new", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
	{"float_phi_diff_sum2_new", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
	{"float_z0_diff_sum_orig", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
	{"float_z0_diff_sum2_orig", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
	{"float_z0_diff_sum_new", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
	{"float_z0_diff_sum2_new", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
	{"float_d0_diff_sum_orig", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
	{"float_d0_diff_sum2_orig", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
	{"float_d0_diff_sum_new", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
	{"float_d0_diff_sum2_new", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
	{"hw_pt_diff_sum_orig", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
	{"hw_pt_diff_sum2_orig", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
	{"hw_pt_diff_sum_new", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
	{"hw_pt_diff_sum2_new", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
	{"hw_eta_diff_sum_orig", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
	{"hw_eta_diff_sum2_orig", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
	{"hw_eta_diff_sum_new", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
	{"hw_eta_diff_sum2_new", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
	{"hw_phi_diff_sum_orig", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
	{"hw_phi_diff_sum2_orig", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
	{"hw_phi_diff_sum_new", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
	{"hw_phi_diff_sum2_new", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
	{"hw_z0_diff_sum_orig", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
	{"hw_z0_diff_sum2_orig", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
	{"hw_z0_diff_sum_new", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
	{"hw_z0_diff_sum2_new", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
	{"hw_d0_diff_sum_orig", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
	{"hw_d0_diff_sum2_orig", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
	{"hw_d0_diff_sum_new", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
	{"hw_d0_diff_sum2_new", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
      });
  }

  // void endRunSummary(edm::Run const& iRun,
  void endRunSummary(edm::Run const&,
		     edm::EventSetup const&,
		     std::map< std::string, std::map < int, double > >* iSummary) const override {
    //add the Stream's partial information to the full information
    for (const auto& outer_map : track_summary_) {
      for (const auto& inner_map : outer_map.second) {
	iSummary->at(outer_map.first).at(inner_map.first) += inner_map.second; // add the resolution_map's summary statistic to the global summary map
	track_summary_[outer_map.first][inner_map.first] = 0; // reset the resolution_map at the end of the run
      }
    }
    std::cout << std::endl;
  }
  static void globalEndRunSummary(edm::Run const&,
				  edm::EventSetup const&,
				  RunContext const* iContext,
				  std::map< std::string, std::map < int, double > >* iSummary) {
    std::cout << "===L1SmartPixelsTrackProducer Track Resolution Global Summary Start===" << std::endl;
    unsigned int skip_counter = 0;
    unsigned int print_counter = 0;
    std::string std_tag = "sum2";
    std::map<int, std::string> match_class_to_name = {
      {0, "//UNKNOWN[0]"},
      {1, "//UNKNOWN[1]"},
      {2, "//UNKNOWN[2]"},
      {4, "isLooselyGenuine"},
      {8, "isGenuine"}, // but is not isLooselyGenuine, may not occur due to tight/loose being inclusive
      {12, "isGenuine"},
      {16, "//UNKNOWN[16]"}
    };
    for (const auto& outer_map : *iSummary) {
      bool do_std = (outer_map.first.find(std_tag) != std::string::npos);
      for (const auto& inner_map : outer_map.second) {
	if (iSummary->at("resolution_entries").at(inner_map.first) < 0.1){
	  skip_counter++;
	  continue;
	}
	print_counter++;
	double output = 0;
	if( do_std )
	  output = std::sqrt(iSummary->at(outer_map.first).at(inner_map.first) / iSummary->at("resolution_entries").at(inner_map.first));
	else
	  output = iSummary->at(outer_map.first).at(inner_map.first) / iSummary->at("resolution_entries").at(inner_map.first);    
	std::cout << "\ttrack_summary[" << outer_map.first << "][" << match_class_to_name[inner_map.first] << "] = ";
	if( do_std )
	  std::cout << "std(";
	else
	  std::cout << "mean(";	
	std::cout << inner_map.second << ", nsamples = " << iSummary->at("resolution_entries").at(inner_map.first) << ") = " << output << std::endl;
      }
    }
    if (print_counter == 0)
      std::cout << "Skipped " << skip_counter << " resolution summary statistics because their corresponding resolution_entries were less than 1" << std::endl;
    std::cout << "===L1SmartPixelsTrackProducer Track Resolution Global Summary End===" << std::endl;
  }

  // void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) override;
  // void beginRun(const Run& iEvent, const EventSetup& iSetup) override {}
  // void endRun(const Run& iEvent, const EventSetup& iSetup) override {}

protected:
private:
  // // ----------constants, enums and typedefs ---------
  // // Relevant constants for the converted track word
  static constexpr float MagConstant =
      CLHEP::c_light / 1.0E3;  //constant is 0.299792458; who knew c_light was in mm/ns?

  typedef TTTrack<Ref_Phase2TrackerDigi_> L1Track;
  typedef std::vector<L1Track> TTTrackCollection;
  
  //-----------------------------------------------------------------------------------------------
  // Containers of parameters passed by python configuration file
  edm::ParameterSet config;

  int MyProcess;       // 11/13/211 for single electrons/muons/pions, 6/15 for pions from ttbar/taus, 1 for inclusive
  bool DebugMode;      // lots of debug printout statements
  bool DebugModeDetailed; // extremely detailed information printouts
  int L1Tk_minNStub;     // require L1 tracks to have >= minNStub (this is mostly for tracklet purposes)
  std::string outputCollectionName_; // name of the output collection
  std::string smartPixelsEmulatorMode_; // mode of the emulator, e.g. passthrough, passthroughFloat, passthroughHW, trackingParticleTruth, correctionlibRegression, correctionlibTPToySmear...
  std::string smartPixelsActiveLayers_;
  std::string smartPixelsCorrectionSet_;

  // -----------------------------------------------------------------------------------------------
  // digiRefit (Tier 2) configuration + machinery. Only used when
  // smartPixelsEmulatorMode_ == "digiRefit"; see mem:smartpixels-tier2-refit-plan
  // and doc/PixelAVAngleResponseSpec.md / doc/Phase2Acceptance.md.
  // Per-layer search half-widths [cm] ([layer-1], TBPX L1-L4). The rphi spread
  // of the beamline-constrained extrapolation GROWS outward (MS-compensation
  // bulge between origin anchor and OT stubs; PSimHit-verified), z shrinks.
  std::vector<double> digiRefitWindowRPhi_;  // module-local x
  std::vector<double> digiRefitWindowZ_;     // module-local y
  int digiRefitMinHits_ = 1;            // min attached IT hits to emit a refit track
  std::string digiRefitUseAngles_ = "alpha";  // "none" | "alpha" | "alphaBeta"
  int digiRefitMaxHitsPerWindow_ = 8;   // combinatorics truncation
  int digiRefitMaxKFUpdates_ = 4;       // KF update / layer cap
  std::string digiRefitGainMode_ = "full";     // "full" | "lut" (RESERVED -> throws)
  // KF numerical guards (spec §6b, FPGA-fidelity handles). Defaults chosen from
  // the PU chi2-tail investigation: physical |H| tops out ~200 and the physical
  // (incl. wrong-hit) chi2-increment ceiling is ~1.9e6, so 1e4/2e6 remove only the
  // non-physical grazing-crossing tail and leave the physical bulk untouched.
  double digiRefitJacobianMaxAbs_ = 1.0e4;  // |H[k][j]| above this (or non-finite) zeroes the column
  double digiRefitChi2UpdateGate_ = 2.0e6;  // scalar update with r^2/S above this is skipped entirely
  // Grazing-angle clamps (spec v0.4 §6b). measAngleMaxAbs is the LOAD-BEARING
  // clamp: a synthesized measured cotAlpha/cotBeta beyond it invalidates THAT
  // ANGLE (the hit keeps its position). predAngleMaxAbs is secondary hygiene,
  // passed to the projector to reject non-physical predicted crossings.
  double digiRefitMeasAngleMaxAbs_ = 12.0;  // |synthesized measured cot| above this clears has{Alpha,Beta}
  double digiRefitPredAngleMaxAbs_ = 12.0;  // |predicted crossing cot| above this invalidates the crossing
  int digiRefitSeedNPar_ = 5;           // 4 (prompt) | 5 (extended + covariance seed)
  // CMSSW_17 backport (Branch B): TTTrack carries NO helix covariance (Ian's
  // covMatrix PR cms-sw/cmssw#51503 is Branch C, not applied here). Default is
  // Branch C: "trackCov" (seed from TTTrack::helixCovMat(), default) | "parametrized".
  std::string digiRefitSeedCovMode_ = "trackCov";  // covMatrix present (cms-sw#51503 omnibus)
  std::vector<double> digiRefitParamSigmas_;       // parametrized-mode seed sigmas: (rInv[cm^-1], phi0, tanL, z0[cm], d0[cm])
  std::string digiRefitPixelavAngleSet_ = "";  // PixelAV angle sigma/bias/valid payload path
  std::string digiRefitSmarthitTrueSet_ = "";  // Stack A (RESERVED): NOT consumed by Tier-2; warns if set
  std::string digiRefitSmarthitFakeSet_ = "";  // Stack B inclusive noise payload path
  std::string digiRefitBdtModel_ = "";  // refit-quality BDT (conifer JSON); empty = keep original trkMVA1

  // Refit-quality BDT (spec §6a). Loaded at construction iff digiRefitBdtModel_
  // is non-empty; evaluated per refit track. The assembly is selected by the
  // loaded model's n_features: 17 = REFIT_BDT_FEATURES v0; 24 = v1 (v0 + the
  // classic-7 TrackQuality hw features of the INPUT track, spec §6a v1); any
  // other count throws at construction.
  static constexpr unsigned kRefitBdtNFeaturesV0 = 17;
  static constexpr unsigned kRefitBdtNFeaturesV1 = 24;
  unsigned digiRefitBdtNFeatures_ = 0;  // set at load: 17 or 24
  std::unique_ptr<conifer::BDT<float, float>> digiRefitBdt_;

  // Correctionlib refs for the PixelAV angle response (loaded iff digiRefit).
  // The synthesis throw lives in the payload (spec §3 fused shift compounds):
  //   cotX_meas = cotX_true + spx_angle_X_shift(layer,cotA,cotB,bLocalY, 1.0)
  //   accept iff spx_angle_valid_flat(in) < spx_angle_valid_prob(in)
  // -> bias + sigma*N(0,1) via HashPRNG, deterministic per input tuple, NO RNG here.
  correction::Correction::Ref corrAlphaSigma_, corrBetaSigma_;
  correction::Correction::Ref corrValidProb_, corrValidFlat_;
  correction::CompoundCorrection::Ref corrAlphaShift_, corrBetaShift_;
  // Stack B inclusive noise-angle distribution (used for no-link digis).
  correction::Correction::Ref corrNoiseCotAlpha_, corrNoiseCotBeta_;

  // Shared helix propagation + TBPX module lookup (built lazily; one impl w/ analyzer).
  mutable smartpixels::HelixProjector projector_;

  // Extra tokens for digiRefit (pixel digis, simlinks, sim tracks).
  edm::EDGetTokenT<edm::DetSetVector<PixelDigi>> pixelDigiToken_;
  edm::EDGetTokenT<edm::DetSetVector<PixelDigiSimLink>> pixelSimLinkToken_;
  edm::EDGetTokenT<edm::SimTrackContainer> simTrackToken_;

  // Loud-failure guard: count TP-matched tracks across the stream; if a
  // truth-required mode sees zero over a reasonable window, throw (no silent
  // passthrough). Reset/accumulated per stream; checked in endStream().
  mutable unsigned long long digiRefitTracksSeen_ = 0;
  mutable unsigned long long digiRefitTruthMatches_ = 0;
  mutable unsigned long long digiRefitZeroCovTracks_ = 0;  // trackCov seeds with an all-zero helixCovMat
  // KF numerical-guard activation counters (spec §6b), reported at endStream.
  mutable unsigned long long digiRefitGatedJacCols_ = 0;   // Jacobian columns zeroed by jacobianMaxAbs
  mutable unsigned long long digiRefitGatedUpdates_ = 0;   // scalar updates skipped by chi2UpdateGate
  mutable unsigned long long digiRefitClampedMeasAngles_ = 0;  // synthesized measured angles cleared by measAngleMaxAbs (spec v0.4 §6b)
  mutable unsigned long long digiRefitClampedPredCross_ = 0;   // predicted crossings rejected by predAngleMaxAbs (spec v0.4 §6b)

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
  edm::ESGetToken<hph::Setup, trackerDTC::SetupRcd> getTokenHPHSetup_;  // Branch C: trackerDTC::SetupRcd

  // correction::Correction::Ref cMap_pt, cMap_phi, cMap_d0;
  correction::CompoundCorrection::Ref cMap_pt, cMap_eta, cMap_phi, cMap_z0, cMap_d0;

  // diagnostics
  mutable std::set<int> tp_track_match_set_ = {};
  mutable std::map< std::string, std::map< int, double > > track_summary_ = {
    {"npars_4_entries", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
    {"npars_5_entries", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
    {"resolution_entries", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
    {"float_pt_diff_sum_orig", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
    {"float_pt_diff_sum2_orig", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
    {"float_pt_diff_sum_new", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
    {"float_pt_diff_sum2_new", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
    {"float_eta_diff_sum_orig", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
    {"float_eta_diff_sum2_orig", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
    {"float_eta_diff_sum_new", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
    {"float_eta_diff_sum2_new", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
    {"float_phi_diff_sum_orig", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
    {"float_phi_diff_sum2_orig", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
    {"float_phi_diff_sum_new", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
    {"float_phi_diff_sum2_new", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
    {"float_z0_diff_sum_orig", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
    {"float_z0_diff_sum2_orig", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
    {"float_z0_diff_sum_new", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
    {"float_z0_diff_sum2_new", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
    {"float_d0_diff_sum_orig", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
    {"float_d0_diff_sum2_orig", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
    {"float_d0_diff_sum_new", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
    {"float_d0_diff_sum2_new", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
    {"hw_pt_diff_sum_orig", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
    {"hw_pt_diff_sum2_orig", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
    {"hw_pt_diff_sum_new", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
    {"hw_pt_diff_sum2_new", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
    {"hw_z0_diff_sum_orig", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
    {"hw_z0_diff_sum2_orig", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
    {"hw_z0_diff_sum_new", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
    {"hw_z0_diff_sum2_new", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
    {"hw_d0_diff_sum_orig", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
    {"hw_d0_diff_sum2_orig", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
    {"hw_d0_diff_sum_new", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
    {"hw_d0_diff_sum2_new", {{0, 0}, {1, 0}, {2, 0}, {4, 0}, {8, 0}, {12, 0}, {16, 0}}},
  };
  
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
  DebugModeDetailed = iConfig.getParameter<bool>("DebugModeDetailed");
  L1TrackInputTag = iConfig.getParameter<edm::InputTag>("L1TrackInputTag");
  MCTruthTrackInputTag = iConfig.getParameter<edm::InputTag>("MCTruthTrackInputTag");
  L1Tk_minNStub = iConfig.getParameter<int>("L1Tk_minNStub");
  outputCollectionName_ = iConfig.getParameter<std::string>("outputCollectionName");
  smartPixelsEmulatorMode_ = iConfig.getParameter<std::string>("smartPixelsEmulatorMode");
  smartPixelsActiveLayers_ = iConfig.getParameter<std::string>("smartPixelsActiveLayers");
  smartPixelsCorrectionSet_ = iConfig.getParameter<edm::FileInPath>("smartPixelsCorrectionSet").fullPath();

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
  getTokenHPHSetup_ = esConsumes<hph::Setup, trackerDTC::SetupRcd>();  // Branch C

  // -----------------------------------------------------------------------------------------------
  // digiRefit (Tier 2) configuration + payload loading. Kept entirely off the
  // other modes' code paths.
  if (smartPixelsEmulatorMode_ == "digiRefit") {
    digiRefitWindowRPhi_ = iConfig.getParameter<std::vector<double>>("digiRefitWindowRPhi");
    digiRefitWindowZ_ = iConfig.getParameter<std::vector<double>>("digiRefitWindowZ");
    if (digiRefitWindowRPhi_.size() != 4 || digiRefitWindowZ_.size() != 4)
      throw cms::Exception("Configuration")
          << "digiRefitWindowRPhi/digiRefitWindowZ must each have exactly 4 per-layer entries, got "
          << digiRefitWindowRPhi_.size() << "/" << digiRefitWindowZ_.size() << ".";
    digiRefitMinHits_ = iConfig.getParameter<int>("digiRefitMinHits");
    digiRefitUseAngles_ = iConfig.getParameter<std::string>("digiRefitUseAngles");
    digiRefitMaxHitsPerWindow_ = iConfig.getParameter<int>("digiRefitMaxHitsPerWindow");
    digiRefitMaxKFUpdates_ = iConfig.getParameter<int>("digiRefitMaxKFUpdates");
    digiRefitGainMode_ = iConfig.getParameter<std::string>("digiRefitGainMode");
    digiRefitJacobianMaxAbs_ = iConfig.getParameter<double>("digiRefitJacobianMaxAbs");
    digiRefitChi2UpdateGate_ = iConfig.getParameter<double>("digiRefitChi2UpdateGate");
    digiRefitMeasAngleMaxAbs_ = iConfig.getParameter<double>("digiRefitMeasAngleMaxAbs");
    digiRefitPredAngleMaxAbs_ = iConfig.getParameter<double>("digiRefitPredAngleMaxAbs");
    digiRefitSeedNPar_ = iConfig.getParameter<int>("digiRefitSeedNPar");
    digiRefitSeedCovMode_ = iConfig.getParameter<std::string>("digiRefitSeedCovMode");
    digiRefitParamSigmas_ = iConfig.getParameter<std::vector<double>>("digiRefitParamSigmas");
    {
      const std::string av = iConfig.getParameter<std::string>("digiRefitPixelavAngleSet");
      digiRefitPixelavAngleSet_ = av.empty() ? std::string() : edm::FileInPath(av).fullPath();
    }
    {
      const std::string fk = iConfig.getParameter<std::string>("digiRefitSmarthitFakeSet");
      digiRefitSmarthitFakeSet_ = fk.empty() ? std::string() : edm::FileInPath(fk).fullPath();
    }
    // Stack A "smarthit_true" (RESERVED): validated but deliberately NOT consumed
    // by Tier-2. Tier-2 takes position from the real digis and angle from the
    // PixelAV response; Stack A only CHARACTERIZES true hits. The key is kept for
    // a future SmartPixels ASIC on-chip readout-inefficiency model (smarthit_true_eff),
    // which the digitizer cannot express. Resolve the path (so a bad path still fails
    // loudly), warn LOUDLY if set, and load nothing.
    {
      const std::string st = iConfig.getParameter<std::string>("digiRefitSmarthitTrueSet");
      digiRefitSmarthitTrueSet_ = st.empty() ? std::string() : edm::FileInPath(st).fullPath();
      if (!digiRefitSmarthitTrueSet_.empty())
        edm::LogWarning("SmartPixelsStackAUnused")
            << "digiRefit: smarthitTrueSet was provided (" << digiRefitSmarthitTrueSet_
            << ") but is NOT consumed by Tier-2 and is being IGNORED. Tier-2 takes hit position "
               "from the real pixel digis and hit angle from the PixelAV angle-response payload; "
               "Stack A ('smarthit_true') only characterizes true hits. This config key is RESERVED "
               "for a future SmartPixels ASIC on-chip readout-inefficiency model (smarthit_true_eff), "
               "which is not yet wired. No Stack A payload will be loaded.";
    }
    {
      const std::string bm = iConfig.getParameter<std::string>("digiRefitBdtModel");
      digiRefitBdtModel_ = bm.empty() ? std::string() : edm::FileInPath(bm).fullPath();
    }

    if (digiRefitGainMode_ == "lut")
      throw cms::Exception("NotImplemented")
          << "digiRefit gainMode='lut' (table-driven Kalman gains) is a reserved placeholder; "
             "use gainMode='full'.";
    if (digiRefitSeedCovMode_ != "trackCov" && digiRefitSeedCovMode_ != "parametrized")
      throw cms::Exception("Configuration")
          << "digiRefit seedCovMode='" << digiRefitSeedCovMode_
          << "' invalid; use 'parametrized' (default) or 'trackCov'.";
    if (digiRefitSeedCovMode_ == "parametrized" && digiRefitParamSigmas_.size() != 5)
      throw cms::Exception("Configuration")
          << "digiRefitParamSigmas must have exactly 5 entries (rInv, phi0, tanL, z0, d0), got "
          << digiRefitParamSigmas_.size() << ".";
    if (digiRefitSeedNPar_ != 4 && digiRefitSeedNPar_ != 5)
      throw cms::Exception("Configuration") << "digiRefit seedNPar must be 4 or 5.";
    if (digiRefitUseAngles_ != "none" && digiRefitUseAngles_ != "alpha" && digiRefitUseAngles_ != "alphaBeta")
      throw cms::Exception("Configuration") << "digiRefit useAngles='" << digiRefitUseAngles_ << "' invalid.";
    if (!(digiRefitJacobianMaxAbs_ > 0.))
      throw cms::Exception("Configuration")
          << "digiRefitJacobianMaxAbs must be a positive double, got " << digiRefitJacobianMaxAbs_ << ".";
    if (!(digiRefitChi2UpdateGate_ > 0.))
      throw cms::Exception("Configuration")
          << "digiRefitChi2UpdateGate must be a positive double, got " << digiRefitChi2UpdateGate_ << ".";
    if (!(digiRefitMeasAngleMaxAbs_ > 0.))
      throw cms::Exception("Configuration")
          << "digiRefitMeasAngleMaxAbs must be a positive double, got " << digiRefitMeasAngleMaxAbs_ << ".";
    if (!(digiRefitPredAngleMaxAbs_ > 0.))
      throw cms::Exception("Configuration")
          << "digiRefitPredAngleMaxAbs must be a positive double, got " << digiRefitPredAngleMaxAbs_ << ".";

    // Refit-quality BDT (spec §6a): load the conifer model at construction and
    // validate its feature count against REFIT_BDT_FEATURES v0 (17). conifer's
    // own decision_function() only assert()s the count (compiled out in opt
    // builds), so we read n_features from the JSON and throw LOUDLY on mismatch.
    if (!digiRefitBdtModel_.empty()) {
      std::ifstream ifs(digiRefitBdtModel_);
      if (!ifs.good())
        throw cms::Exception("Configuration")
            << "digiRefitBdtModel could not be opened: " << digiRefitBdtModel_;
      unsigned nFeat = 0;
      try {
        const nlohmann::json j = nlohmann::json::parse(ifs);
        nFeat = j.at("n_features").get<unsigned>();
      } catch (const std::exception& e) {
        throw cms::Exception("Configuration")
            << "digiRefitBdtModel " << digiRefitBdtModel_
            << " is not a valid conifer JSON (n_features unreadable): " << e.what();
      }
      if (nFeat != kRefitBdtNFeaturesV0 && nFeat != kRefitBdtNFeaturesV1)
        throw cms::Exception("Configuration")
            << "digiRefitBdtModel " << digiRefitBdtModel_ << " has n_features=" << nFeat
            << " but REFIT_BDT_FEATURES requires exactly " << kRefitBdtNFeaturesV0 << " (v0) or "
            << kRefitBdtNFeaturesV1 << " (v1).";
      digiRefitBdtNFeatures_ = nFeat;
      digiRefitBdt_ = std::make_unique<conifer::BDT<float, float>>(digiRefitBdtModel_);
    }

    // Angle-synthesis randomness is PAYLOAD-SIDE (HashPRNG fused shift compounds,
    // spec §3): deterministic per input tuple, replay-identical, no engine. A local
    // engine (seeded per produce() from hash(module label, run, lumi, event), see
    // digiRefitSeed()) is created ONLY when the optional smarthitFakeSet noise-angle
    // model is configured -- its inverse-CDF quantile draws still need a uniform
    // stream. The RandomNumberGeneratorService is deliberately NOT used
    // (its per-stream engine makes outputs order-dependent).

    pixelDigiToken_ =
        consumes<edm::DetSetVector<PixelDigi>>(iConfig.getParameter<edm::InputTag>("pixelDigiInputTag"));
    pixelSimLinkToken_ =
        consumes<edm::DetSetVector<PixelDigiSimLink>>(iConfig.getParameter<edm::InputTag>("pixelDigiSimLinkInputTag"));
    simTrackToken_ = consumes<edm::SimTrackContainer>(iConfig.getParameter<edm::InputTag>("simTrackInputTag"));

    // PixelAV angle-response payload (required, non-empty). See
    // doc/PixelAVAngleResponseSpec.md — the load-bearing correction/compound names.
    if (digiRefitPixelavAngleSet_.empty())
      throw cms::Exception("Configuration")
          << "digiRefit requires a non-empty pixelavAngleSet (PixelAV angle-response payload).";
    auto aset = correction::CorrectionSet::from_file(digiRefitPixelavAngleSet_);
    corrAlphaSigma_ = aset->at("spx_angle_alpha_sigma");
    corrBetaSigma_ = aset->at("spx_angle_beta_sigma");
    corrValidProb_ = aset->at("spx_angle_valid_prob");
    try {
      corrValidFlat_ = aset->at("spx_angle_valid_flat");
      corrAlphaShift_ = aset->compound().at("spx_angle_alpha_shift");
      corrBetaShift_ = aset->compound().at("spx_angle_beta_shift");
    } catch (const std::out_of_range&) {
      throw cms::Exception("Configuration")
          << "pixelavAngleSet '" << digiRefitPixelavAngleSet_
          << "' predates the HashPRNG synthesis-throw factorization (missing "
          << "spx_angle_valid_flat / spx_angle_{alpha,beta}_shift). Regenerate it with "
          << "ngtagger-train/eval_spixel_angles/extract_pixelav_angle_payload.py "
          << "(additions are purely additive; plain corrections stay bit-identical).";
    }

    // Optional Stack B inclusive noise-angle distribution. When absent, no-link
    // (noise) digis contribute position only (angles disabled for that hit).
    if (!digiRefitSmarthitFakeSet_.empty()) {
      auto fset = correction::CorrectionSet::from_file(digiRefitSmarthitFakeSet_);
      try {
        corrNoiseCotAlpha_ = fset->at("smarthit_noise_cotAlpha");
      } catch (const std::out_of_range&) {
      }
      try {
        corrNoiseCotBeta_ = fset->at("smarthit_noise_cotBeta");
      } catch (const std::out_of_range&) {
      }
    }
  }

  produces<TTTrackCollection>(outputCollectionName_);

  // digiRefit sidecar: one product per refit track collection (spec §2), same
  // instance label as the tracks so downstream adapters resolve them together
  // and rely on the 1:1 row-sync invariant (spec §1).
  if (smartPixelsEmulatorMode_ == "digiRefit")
    produces<smartpixels::SmartPixelsRefitSidecar>(outputCollectionName_);

  // -----------------------------------------------------------------------------------------------
  // Correctionlib loading
  if( smartPixelsEmulatorMode_ == "correctionlibRegression") {
    auto cSet = correction::CorrectionSet::from_file(smartPixelsCorrectionSet_);
    std::string cmap_pt_name = "pt_relative_smear_compound_" + smartPixelsActiveLayers_;
    std::string cmap_eta_name = "eta_relative_smear_compound_" + smartPixelsActiveLayers_;
    std::string cmap_phi_name = "phi_relative_smear_compound_" + smartPixelsActiveLayers_;
    std::string cmap_z0_name = "z0_relative_smear_compound_" + smartPixelsActiveLayers_;
    std::string cmap_d0_name = "d0_relative_smear_compound_" + smartPixelsActiveLayers_;
    cMap_pt = cSet->compound().at(cmap_pt_name);
    cMap_eta = cSet->compound().at(cmap_eta_name);
    cMap_phi = cSet->compound().at(cmap_phi_name);
    cMap_z0 = cSet->compound().at(cmap_z0_name);
    cMap_d0 = cSet->compound().at(cmap_d0_name);
  }
  else if( smartPixelsEmulatorMode_ == "correctionlibTPToySmear") {
    auto cSet = correction::CorrectionSet::from_file(smartPixelsCorrectionSet_);
    std::string cmap_pt_name = "pt_smear_compound_" + smartPixelsActiveLayers_;
    std::string cmap_eta_name = "eta_smear_compound_" + smartPixelsActiveLayers_;
    std::string cmap_phi_name = "phi_smear_compound_" + smartPixelsActiveLayers_;
    std::string cmap_z0_name = "z0_smear_compound_" + smartPixelsActiveLayers_;
    std::string cmap_d0_name = "d0_smear_compound_" + smartPixelsActiveLayers_;
    cMap_pt = cSet->compound().at(cmap_pt_name);
    cMap_eta = cSet->compound().at(cmap_eta_name);
    cMap_phi = cSet->compound().at(cmap_phi_name);
    cMap_z0 = cSet->compound().at(cmap_z0_name);
    cMap_d0 = cSet->compound().at(cmap_d0_name);
  }
}

// DESTRUCTOR
L1SmartPixelsTrackProducer::~L1SmartPixelsTrackProducer() {
}
void L1SmartPixelsTrackProducer::beginStream(edm::StreamID) {
}
void L1SmartPixelsTrackProducer::endStream() {
  // Loud-failure guard: truth-required digiRefit must never silently degrade to
  // per-track passthrough (stale/mismatched maps or wrong input posture would
  // otherwise silently break every downstream study).
  if (smartPixelsEmulatorMode_ == "digiRefit" && digiRefitTracksSeen_ > 0) {
    // KF numerical-guard activation summary (spec §6b). These count the
    // non-physical grazing-crossing pathologies removed by the two guards.
    edm::LogVerbatim("SmartPixelsTrackProducer")
        << "digiRefit KF guards: jacobianMaxAbs zeroed " << digiRefitGatedJacCols_
        << " Jacobian column(s); chi2UpdateGate skipped " << digiRefitGatedUpdates_
        << " scalar update(s); measAngleMaxAbs cleared " << digiRefitClampedMeasAngles_
        << " synthesized angle(s); predAngleMaxAbs rejected " << digiRefitClampedPredCross_
        << " predicted crossing(s) over " << digiRefitTracksSeen_ << " tracks.";
    if (digiRefitTruthMatches_ == 0)
      throw cms::Exception("SmartPixelsTruthMissing")
          << "digiRefit processed " << digiRefitTracksSeen_
          << " tracks with ZERO truth (TP) matches. The TTTrackAssociationMap almost certainly "
             "does not correspond to the input track collection (stale map + remade tracks, or "
             "wrong truthSource posture). Run posture B (DIGI + L1TrackTrigger in-job) so tracks, "
             "digis, simlinks and maps are self-consistent.";
    if (digiRefitSeedCovMode_ == "trackCov" && digiRefitZeroCovTracks_ == digiRefitTracksSeen_)
      throw cms::Exception("SmartPixelsSeedCovMissing")
          << "digiRefit seedCovMode='trackCov' but ALL " << digiRefitTracksSeen_
          << " seed tracks carried an all-zero helixCovMat. These look like schema-evolved "
             "old-layout file tracks - run posture B so the in-job fit fills the covariance, "
             "or use seedCovMode='parametrized'.";
  }
}

// PRODUCE
// void L1SmartPixelsTrackProducer::produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const {
void L1SmartPixelsTrackProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  if (!(MyProcess == 13 || MyProcess == 11 || MyProcess == 211 || MyProcess == 6 || MyProcess == 15 ||
        MyProcess == 1)) {
    edm::LogVerbatim("SmartPixelsTrackProducer") << "The specified MyProcess is invalid! Exiting...";
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
  float b_field = bFieldHandle.product()->inTesla(GlobalPoint(0, 0, 0)).z();

  edm::ESHandle<hph::Setup> hphHandle = iSetup.getHandle(getTokenHPHSetup_);

  const TrackerTopology* const tTopo = tTopoHandle.product();
  const TrackerGeometry* const theTrackerGeom = tGeomHandle.product();
  const hph::Setup* hphSetup = hphHandle.product();

  // -----------------------------------------------------------------------------------------------
  // digiRefit (Tier 2) per-event inputs: real pixel digis + simlinks, SimTrack
  // origin momenta (parent-angle synthesis), shared projector, local RNG engine,
  // and the per-collection refit sidecar (spec §2).
  edm::Handle<edm::DetSetVector<PixelDigi>> drDigis;
  edm::Handle<edm::DetSetVector<PixelDigiSimLink>> drSimlinks;
  smartpixels::ParentMomentumMap drParentMom;
  std::unique_ptr<CLHEP::MixMaxRng> drEngineOwned;
  CLHEP::HepRandomEngine* drEngine = nullptr;
  std::array<bool, 4> drActiveLayer{{false, false, false, false}};
  auto sidecar = std::make_unique<smartpixels::SmartPixelsRefitSidecar>();
  if (smartPixelsEmulatorMode_ == "digiRefit") {
    iEvent.getByToken(pixelDigiToken_, drDigis);
    iEvent.getByToken(pixelSimLinkToken_, drSimlinks);
    edm::Handle<edm::SimTrackContainer> drSimTracks;
    iEvent.getByToken(simTrackToken_, drSimTracks);
    edm::Handle<std::vector<TrackingParticle>> drTPs;
    iEvent.getByToken(TrackingParticleToken_, drTPs);
    // (eventId, trackId)-keyed parent momenta: TPs cover signal + pileup
    // parents; the signal-only SimTrack container backstops pruned TPs.
    drParentMom = smartpixels::buildParentMomentumMap(*drTPs, drSimTracks.product());
    for (size_t i = 0; i < drActiveLayer.size() && i < smartPixelsActiveLayers_.size(); ++i)
      drActiveLayer[i] = (smartPixelsActiveLayers_[i] == '1');
    if (std::none_of(drActiveLayer.begin(), drActiveLayer.end(), [](bool b) { return b; }))
      throw cms::Exception("Configuration")
          << "digiRefit requires at least one active layer (smartPixelsActiveLayers='"
          << smartPixelsActiveLayers_ << "').";
    projector_.build(*theTrackerGeom, *tTopo, 4);
    // Angle synthesis is engine-free (payload-side HashPRNG, spec §3). A local
    // engine is needed ONLY for the optional Stack-B noise-angle model's
    // inverse-CDF quantile draws; seeded deterministically from (module label,
    // run, lumi, event): event-order-independent and split-job invariant (see
    // digiRefitSeed()). The EDM module label (distinct for prompt vs extended)
    // decorrelates the two producers' draw streams.
    if (corrNoiseCotAlpha_ && corrNoiseCotBeta_) {
      const long drSeed = digiRefitSeed(moduleDescription().moduleLabel(),
                                        iEvent.id().run(),
                                        iEvent.id().luminosityBlock(),
                                        static_cast<unsigned long long>(iEvent.id().event()));
      drEngineOwned = std::make_unique<CLHEP::MixMaxRng>(drSeed);
      drEngine = drEngineOwned.get();
    }
  }

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

      if (DebugModeDetailed)
	edm::LogVerbatim("SmartPixelsTrackProducer")
	  << "stub info " << stubIter << " tmp_stub (x,y,z) = (" << tmp_stub_x << ", " << tmp_stub_y << ", " << tmp_stub_z << ")\t"
	  << "trigDisplace= " << trigDisplace << " trigOffset= " << trigOffset  << " trigPos= " << trigPos << " trigBend= " << trigBend << std::endl;
      

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
  }
  unsigned int matched_tracks = 0;
  unsigned int unmatched_tracks = 0;
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
    if (iterL1Track->nFitPars() == 5) {
      float tmp_trk_x0 = iterL1Track->POCA().x();
      float tmp_trk_y0 = iterL1Track->POCA().y();
      tmp_trk_d0 = tmp_trk_x0 * sin(tmp_trk_phi) - tmp_trk_y0 * cos(tmp_trk_phi);
    }

    float tmp_trk_chi2 = iterL1Track->chi2();
    float tmp_trk_chi2rphi = iterL1Track->chi2XY();
    float tmp_trk_chi2rz = iterL1Track->chi2Z();
    float tmp_trk_bendchi2 = iterL1Track->chi2BendRed();
    float tmp_trk_MVA1 = iterL1Track->trkMVA1();

    std::vector<edm::Ref<edmNew::DetSetVector<TTStub<Ref_Phase2TrackerDigi_> >, TTStub<Ref_Phase2TrackerDigi_> > >
        stubRefs = iterL1Track->getStubRefs();
    int tmp_trk_nstub = (int)stubRefs.size();
    int ndof = 2 * tmp_trk_nstub - iterL1Track->nFitPars();
    int ndofrphi = tmp_trk_nstub - iterL1Track->nFitPars() + 2;
    int ndofrz = tmp_trk_nstub - 2;
    float tmp_trk_chi2_dof = (float)tmp_trk_chi2 / ndof;
    float tmp_trk_chi2rphi_dof = (float)tmp_trk_chi2rphi / ndofrphi;
    float tmp_trk_chi2rz_dof = (float)tmp_trk_chi2rz / ndofrz;

    int tmp_trk_seed = (int)iterL1Track->trackSeedType();

    unsigned int tmp_trk_phiSector = iterL1Track->phiSector();
    int tmp_trk_etaSector = hph.etaSector();

    if (DebugModeDetailed)
	edm::LogVerbatim("SmartPixelsTrackProducer")
	  << "tmp_trk_nPSstub_hitpattern\ttmp_trk_n2Sstub_hitpattern\ttmp_trk_nLostPSstub_hitpattern\ttmp_trk_nLost2Sstub_hitpattern"
	  << tmp_trk_nPSstub_hitpattern << "\t" << tmp_trk_n2Sstub_hitpattern << "\t" << tmp_trk_nLostPSstub_hitpattern << "\t" << tmp_trk_nLost2Sstub_hitpattern << "\t"
	  << std::endl
	  << "\ttmp_trk_nLoststub_V1_hitpattern\ttmp_trk_nLoststub_V2_hitpattern\ttmp_trk_bendchi2\ttmp_trk_MVA1"
	  << "\t" << tmp_trk_nLoststub_V1_hitpattern << "\t" << tmp_trk_nLoststub_V2_hitpattern << "\t" << tmp_trk_bendchi2 << "\t" << tmp_trk_MVA1
	  << std::endl
	  << "\ttmp_trk_chi2_dof\ttmp_trk_chi2rphi_dof\ttmp_trk_chi2rz_dof\ttmp_trk_seed\ttmp_trk_phiSector\ttmp_trk_etaSector"
	  << "\t" << tmp_trk_chi2_dof << "\t" << tmp_trk_chi2rphi_dof << "\t" << tmp_trk_chi2rz_dof << "\t" << tmp_trk_seed
	  << "\t" << tmp_trk_phiSector << "\t" << tmp_trk_etaSector
	  << std::endl;

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
    int tp_track_match = 0;
    if (MCTruthTTTrackHandle->isUnknown(l1track_ptr)) {
      // Appears that EVERY track is this class, no isLooselyGenuine or isGenuine or isCombinatoric yet seen in SmartPixelsTrackProducer. Is this a bug or the current state of this function everywhere?
      tmp_trk_unknown = 1;
      tp_track_match += (1 << 0);
    }
    // Reserve tp_track_match += (1 << 1) for isFake if this is distinguishable
    if (MCTruthTTTrackHandle->isLooselyGenuine(l1track_ptr)) {
      tmp_trk_loose = 1;
      tp_track_match += (1 << 2);
    }
    if (MCTruthTTTrackHandle->isGenuine(l1track_ptr)) {
      tmp_trk_genuine = 1;
      tp_track_match += (1 << 3);
    }
    if (MCTruthTTTrackHandle->isCombinatoric(l1track_ptr)) {
      tmp_trk_combinatoric = 1;
      tp_track_match += (1 << 4);
    }
    tp_track_match_set_.insert(tp_track_match);

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
      unmatched_tracks++;
      myFake = 0;
      if (DebugMode) {
        edm::LogVerbatim("SmartPixelsTrackProducer") << "TP not matched to track: myFake = " << myFake << " pdgId = " << tmp_matchtp_pdgid;
      }
    }
    else {
      matched_tracks++;
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
                              iterL1Track->nFitPars(),
                              b_field
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
                              iterL1Track->nFitPars(),
                              b_field
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
      // NOTE: L1 Tracks should be valid only above 2GeV; the spec permits forming tracks below this, but around 1.9GeV there's conversion failures picked up 
      // by GTTInputProducer, so we should avoid altering any such tracks and leave them as-is. Physics wise there should be no impact, only tracks above 2GeV should be used
      if (my_tp.isNull() || (!my_tp.isNull() && my_tp->pt() <= 1.95)) {
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
                                iterL1Track->nFitPars(),
                                b_field
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
        auto tmp_matchtp_rInv = my_tp->charge() * MagConstant * b_field / (tmp_matchtp_pt * 100.0);
        auto tmp_matchtp_tanL = my_tp->p4().pz() / tmp_matchtp_pt;
        L1Track track = L1Track(tmp_matchtp_rInv,
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
                                b_field
        );
        track.setPhiSector(iterL1Track->phiSector());
        track.setEtaSector(iterL1Track->etaSector());
        track.setTrackSeedType(iterL1Track->trackSeedType());
        track.setStubRefs(iterL1Track->getStubRefs());
        // pt consistency
        track.setChi2BendRed(
          StubPtConsistency::getConsistency(track, theTrackerGeom, tTopo, b_field, iterL1Track->nFitPars()));
        track.setTrackWordBits();
        //auto oldChi2Bend = iterL1Track->stubPtConsistency();
        //auto newChi2Bend = track.stubPtConsistency();
        //std::cout << "change bendChi2 % " << ((newChi2Bend - oldChi2Bend) / oldChi2Bend) << std::endl;
        /* bchi2 = pd.read_csv("change_bendchi2.txt")
        bchi2 = bchi2.tp_bendch2_minus_oldchi2_dividedby_oldchi2
        np.min(bchi2), np.max(bchi2), np.mean(bchi2), np.std(bchi2)
        (np.float64(-0.945288), np.float64(7.85453), np.float64(0.01088359448426543), np.float64(0.2033839357486676)) */
        //add the track to the output collection
        outputTracks->push_back(track);
        /*
        if (tmp_matchtp_pt > 1.95 && abs(iterL1Track->momentum().perp() - tmp_matchtp_pt)/tmp_matchtp_pt > 0.5)
          std::cout << "track input pt: " << iterL1Track->momentum().perp() << " rInv: " << iterL1Track->rInv() << std::endl
                    << "tp input pt: " << tmp_matchtp_pt << " rInv: " << tmp_matcht_rInv << std::endl
                    << " track output pt: " << track.momentum().perp() << " rInv: " << track.rInv() << std::endl;
        */
      }
    }
    else if(smartPixelsEmulatorMode_ == "correctionlibRegression") {
      if (my_tp.isNull() || (!my_tp.isNull() && my_tp->pt() <= 1.95)) {
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
                                iterL1Track->nFitPars(),
                                b_field
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
	// correct the track parameters towards the tracking particle truth level based on the parameterized toy detector
        // double thePT = std::abs(MagConstant / theRInv_ * aBField / 100.0);  // Rinv is in cm-1

        // auto tmp_matchtp_rInv = my_tp->charge() * MagConstant * b_field / (tmp_matchtp_pt * 100.0);
        double iterL1Track_pt = my_tp->charge() * MagConstant * b_field / (iterL1Track->rInv() * 100.0);
        // auto tmp_matchtp_tanL = my_tp->p4().pz() / tmp_matchtp_pt;
        double iterL1Track_eta = std::asinh(iterL1Track->tanL());;

	std::vector<std::variant<int,double,std::string>> inputs = {};
	// inputs.emplace_back(fast_relevant_syst_for_shape_corr.at(hadronFlavour.at(i)));
	// inputs.emplace_back(std::clamp(abs(static_cast<double>(eta.at(i))), 0., eta_max-0.001)); //copy nanoAOD-tools clamping, but give it the proper eta max... eps=1.e-3
	// inputs.emplace_back(std::clamp(static_cast<double>(pt.at(i)), 20.0001, 9999.9999)); //need to parse the corrector to get the true pt bounds...eps=1.e-4
	// inputs.emplace_back(disc.at(i));
	// auto val = cMap->evaluate(inputs);
	inputs.emplace_back(tmp_matchtp_pt);
	inputs.emplace_back(tmp_matchtp_eta);
	inputs.emplace_back(tmp_matchtp_phi);
	inputs.emplace_back(tmp_matchtp_z0);
	inputs.emplace_back(tmp_matchtp_d0);
	inputs.emplace_back(tp_track_match); //encodes the matching flags as bits, 0:isUnknown,1:isLooselyGenuine,2:isGenuine,3:isCombinatoric
	// inputs.emplace_back(iterL1Track->pt());
	inputs.emplace_back(iterL1Track_pt);
	// inputs.emplace_back(iterL1Track->eta());
	inputs.emplace_back(iterL1Track_eta);
	inputs.emplace_back(iterL1Track->phi());
	inputs.emplace_back(iterL1Track->z0());
	inputs.emplace_back(iterL1Track->d0());


	// match the original tp phi
	// auto clib_phi = tmp_matchtp_phi;
	auto clib_phi = tmp_matchtp_phi + cMap_phi->evaluate(inputs);

	// match the original tanL under the assumption that the pt_new / pt_old scaling applies to the pz_new / pz_old, ergo cancelling in tanL calculation
        // double iterL1Track_eta = std::asinh(iterL1Track->tanL()); --> tanL = sinh(eta)
	auto clib_eta = my_tp->p4().eta() + cMap_eta->evaluate(inputs);
	auto clib_tanL = std::sinh(clib_eta);

	// auto tmp_pt = tmp_matchtp_pt; //cMap_pt->evaluate(inputs);
	auto clib_pt = tmp_matchtp_pt + cMap_pt->evaluate(inputs);
	auto clib_rInv = my_tp->charge() * MagConstant * b_field / (clib_pt * 100.0);

	//FIXME: update to dedicated z0 estimate when not using swizzled d0_smear * z0_track_tp_difference
	auto clib_z0 = tmp_matchtp_z0 + cMap_z0->evaluate(inputs);

	// auto tmp_d0 = tmp_matchtp_d0; //cMap_d0->evaluate(inputs);
	auto clib_d0 = tmp_matchtp_d0 + cMap_d0->evaluate(inputs);
	if (iterL1Track->nFitPars() == 4)
	  clib_d0 = 0; // Reset to 0 if we're looking at a 4-parameter track


        L1Track track = L1Track(clib_rInv,
                                clib_phi,
                                clib_tanL,
                                clib_z0,
                                clib_d0,
                                iterL1Track->chi2XY(), // Keep from original track?
                                iterL1Track->chi2Z(), // Keep from original track?
                                iterL1Track->trkMVA1(), // Keep from original track?
                                iterL1Track->trkMVA2(), // Keep from original track?
                                iterL1Track->trkMVA3(), // Keep from original track?
                                iterL1Track->hitPattern(), // Keep from original track?
                                iterL1Track->nFitPars(), // Keep from original track for sure
                                b_field
        );
        track.setPhiSector(iterL1Track->phiSector());
        track.setEtaSector(iterL1Track->etaSector());
        track.setTrackSeedType(iterL1Track->trackSeedType());
        track.setStubRefs(iterL1Track->getStubRefs());
        // pt consistency
        track.setChi2BendRed(
			     StubPtConsistency::getConsistency(track, theTrackerGeom, tTopo, b_field, iterL1Track->nFitPars()));
        track.setTrackWordBits();
        //auto oldChi2Bend = iterL1Track->stubPtConsistency();
        //auto newChi2Bend = track.stubPtConsistency();
        //std::cout << "change bendChi2 % " << ((newChi2Bend - oldChi2Bend) / oldChi2Bend) << std::endl;
        /* bchi2 = pd.read_csv("change_bendchi2.txt")
        bchi2 = bchi2.tp_bendch2_minus_oldchi2_dividedby_oldchi2
        np.min(bchi2), np.max(bchi2), np.mean(bchi2), np.std(bchi2)
        (np.float64(-0.945288), np.float64(7.85453), np.float64(0.01088359448426543), np.float64(0.2033839357486676)) */
        //add the track to the output collection
        outputTracks->push_back(track);
        /*
        if (tmp_matchtp_pt > 1.95 && abs(iterL1Track->momentum().perp() - tmp_matchtp_pt)/tmp_matchtp_pt > 0.5)
          std::cout << "track input pt: " << iterL1Track->momentum().perp() << " rInv: " << iterL1Track->rInv() << std::endl
                    << "tp input pt: " << tmp_matchtp_pt << " rInv: " << tmp_matcht_rInv << std::endl
                    << " track output pt: " << track.momentum().perp() << " rInv: " << track.rInv() << std::endl;
        */

	//DEBUG prints
	if( DebugMode && this_l1track < 5 ) {
	  edm::LogVerbatim("SmartPixelsTrackProducer") << "track " << this_l1track << " inputs:"
		    << "(pt_tp, eta_tp, phi_tp, z0_tp, d0_tp, tp_track_match, pt_track, eta_track, phi_track, z0_track, d0_track)";
	  unsigned int visit_counter = 0;
	  for (auto&& inp : inputs) {
	    if (visit_counter == 0 || visit_counter == 5 || visit_counter == 6)
	      edm::LogVerbatim("SmartPixelsTrackProducer") << std::endl;
	    std::visit([](auto&& arg){ edm::LogVerbatim("SmartPixelsTrackProducer") << " " << arg;}, inp);
	    visit_counter++;
	  }
	  edm::LogVerbatim("SmartPixelsTrackProducer") << std::endl;
	  edm::LogVerbatim("SmartPixelsTrackProducer") << "outputs:"
		    // << " clib_pt=" << cMap_pt->evaluate(inputs) << " clib_rInv=" <<  my_tp->charge() * MagConstant * b_field / (cMap_pt->evaluate(inputs) * 100.0)
		    // << " clib_phi=" << cMap_phi->evaluate(inputs)
		    // << " clib_d0" << cMap_d0->evaluate(inputs)
		    << " clib_pt= " << clib_pt << " clib_rInv= " <<  clib_rInv
		    << " clib_phi= " << clib_phi
		    << " clib_d0= " << clib_d0
		    << " clib_z0= " << clib_z0
		    << std::endl;
	  edm::LogVerbatim("SmartPixelsTrackProducer") << "input track " << this_l1track << " L1Track parameters:" << std::endl;
	  
	  edm::LogVerbatim("SmartPixelsTrackProducer") << "\trInv original :: new\t"
		    << (iterL1Track->rInv()) << "    -    "
		    << (track.rInv()) << std::endl;
	  edm::LogVerbatim("SmartPixelsTrackProducer") << "\tpt original :: new\t"
		    << (my_tp->charge() * MagConstant * b_field / (iterL1Track->rInv() * 100.0)) << "    -    "
		    << (my_tp->charge() * MagConstant * b_field / (track.rInv() * 100.0)) << std::endl;
	  edm::LogVerbatim("SmartPixelsTrackProducer") << "\ttanL original :: new\t"
		    << (iterL1Track->tanL()) << "    -    "
		    << (track.tanL()) << std::endl;
	  edm::LogVerbatim("SmartPixelsTrackProducer") << "\teta original :: new\t"
		    << (std::asinh(iterL1Track->tanL())) << "    -    "
		    << (std::asinh(track.tanL())) << std::endl;
	  edm::LogVerbatim("SmartPixelsTrackProducer") << "\tphi original :: new\t"
		    << (iterL1Track->phi()) << "    -    "
		    << (track.phi()) << std::endl;
	  edm::LogVerbatim("SmartPixelsTrackProducer") << "\tz0 original :: new\t"
		    << (iterL1Track->z0()) << "    -    "
		    << (track.z0()) << std::endl;
	  edm::LogVerbatim("SmartPixelsTrackProducer") << "\td0 original :: new\t"
		    << (iterL1Track->d0()) << "    -    "
		    << (track.d0()) << std::endl;
	} // end debug logging
      }
    } // end correctionlibRegression path
    else if(smartPixelsEmulatorMode_ == "digiRefit") {
      // ------------------------------------------------------------------------------------------
      // Tier-2 interim refit: project the OT track into BPIX, window-collect the
      // REAL pixel digis, classify each via PixelDigiSimLink, synthesize angle
      // measurements (parent-truth incidence ⊕ PixelAV response payload), then
      // Kalman-update the OT helix seeded from the track and its covariance.
      // Exactly one output track per input track (refit or passthrough fallback).
      ++digiRefitTracksSeen_;
      if (!my_tp.isNull())
        ++digiRefitTruthMatches_;

      const bool useAlpha = (digiRefitUseAngles_ != "none");
      const bool useBeta = (digiRefitUseAngles_ == "alphaBeta");

      // Truth identity of THIS track, for same-TP vs other-TP digi classification.
      // (Class does not change the measurement model here — other-TP digis carry
      // their OWN parent's angle, which is the correct fake phenomenology.)
      std::set<unsigned int> drMatchedSimIds;
      EncodedEventId drMatchedEvtId;
      if (!my_tp.isNull()) {
        drMatchedEvtId = my_tp->eventId();
        for (const auto& g4 : my_tp->g4Tracks())
          drMatchedSimIds.insert(g4.trackId());
      }

      // ---- sidecar accumulators (spec §2): one hitInfo entry per valid layer
      // crossing attempted (whether or not a hit is accepted), one trackInfo
      // entry per track. Filled below and appended to the sidecar under the 1:1
      // invariant. layerHitMask (bit0=L1..bit3=L4) feeds the compact word.
      std::vector<smartpixels::SmartPixelsRefitHitInfo> drHitInfo;
      smartpixels::SmartPixelsRefitTrackInfo drTrackInfo;
      uint8_t drLayerHitMask = 0;
      unsigned drMaxWindowMult = 0;
      bool drAnyWindowTruncated = false;
      // Refit-BDT sum-of-squared-pull accumulators (spec §6a features 5-8),
      // summed over all applied scalar updates across the track's crossings.
      double drSumPullX2 = 0., drSumPullY2 = 0., drSumPullAlpha2 = 0., drSumPullBeta2 = 0.;

      // ---- seed state (INVR, PHI0, TANL, Z0, D0) and covariance ----
      ROOT::Math::SVector<double, 5> a;
      a[0] = iterL1Track->rInv();
      a[1] = iterL1Track->momentum().phi();
      a[2] = iterL1Track->tanL();
      a[3] = iterL1Track->z0();
      a[4] = (digiRefitSeedNPar_ == 5 && iterL1Track->nFitPars() == 5) ? tmp_trk_d0 : 0.0;
      // Seed snapshot for the refit-BDT parameter deltas (spec §6a features 11-15).
      const ROOT::Math::SVector<double, 5> aSeed = a;

      ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>> C;
      bool seedCovOK = true;
      if (digiRefitSeedCovMode_ == "trackCov") {
        // Branch C: seed the KF from the TTTrack helix covariance (added by
        // cms-sw#51503 / L1TK-dev-15_1 omnibus). Re-emulated tracklet tracks
        // carry a real 5x5 helixCovMat (filled in L1FPGATrackProducer via the
        // 18-arg TTTrack ctor). Order matches enum Hpar{INVR,PHI0,TANL,Z0,D0}.
        const auto& tc = iterL1Track->helixCovMat();
        double diagSum = 0.;
        for (int i = 0; i < 5; ++i) {
          for (int j = 0; j <= i; ++j)
            C(i, j) = tc(i, j);
          diagSum += std::abs(tc(i, i));
        }
        if (!(diagSum > 0.)) {
          // Schema-evolved old-layout file tracks carry a default (all-zero)
          // covariance; refitting from a singular seed would be garbage. Fall
          // back to passthrough for this track and count it -- endStream throws
          // if EVERY track looked like this (wrong input posture).
          seedCovOK = false;
          ++digiRefitZeroCovTracks_;
        }
      } else {  // "parametrized" fallback/ablation seed
        for (int i = 0; i < 5; ++i)
          C(i, i) = digiRefitParamSigmas_[i] * digiRefitParamSigmas_[i];
      }
      // 4-par seeds carry no d0 information: weak prior so the IT hits determine d0.
      if (digiRefitSeedNPar_ == 4 || !(C(4, 4) > 0.))
        C(4, 4) = std::max(C(4, 4), 0.25);  // (0.5 cm)^2

      const MagneticField& drField = *bFieldHandle.product();
      const auto makeHelix = [&](const ROOT::Math::SVector<double, 5>& s) {
        smartpixels::HelixParams hp;
        hp.rInv = s[0];
        hp.phi0 = s[1];
        hp.tanL = s[2];
        hp.z0 = s[3];
        // POCA position from (d0, phi0); satisfies d0 = x0*sin(phi0) - y0*cos(phi0).
        hp.x0 = s[4] * std::sin(s[1]);
        hp.y0 = -s[4] * std::cos(s[1]);
        hp.pt = std::abs(MagConstant * b_field / (s[0] * 100.0));
        return hp;
      };

      int nAcceptedHits = 0;
      int nUpdates = 0;

      for (int layer = 1; layer <= projector_.nLayers(); ++layer) {
        if (!seedCovOK || nUpdates >= digiRefitMaxKFUpdates_)
          break;
        if (!drActiveLayer[layer - 1])
          continue;

        const smartpixels::Crossing cx =
            projector_.crossLayer(makeHelix(a), layer, drField, digiRefitPredAngleMaxAbs_);
        if (!cx.valid) {
          // Predicted-angle clamp (spec v0.4 §6b hygiene) or genuine off-acceptance.
          // Count the clamp activations for the endStream report (attribute a
          // rejected crossing to the clamp when the unclamped projection WOULD have
          // been valid but exceeded the predicted-angle bound).
          const smartpixels::Crossing cxUnclamped = projector_.crossLayer(makeHelix(a), layer, drField, 0.);
          if (cxUnclamped.valid &&
              (std::abs(cxUnclamped.cotAlpha) > digiRefitPredAngleMaxAbs_ ||
               std::abs(cxUnclamped.cotBeta) > digiRefitPredAngleMaxAbs_))
            ++digiRefitClampedPredCross_;
          continue;
        }
        const PixelGeomDetUnit* pixDet = cx.det;
        const PixelTopology& pixTopo = pixDet->specificTopology();
        const std::pair<float, float> pitch = pixTopo.pitch();
        const double sigX = pitch.first / std::sqrt(12.0);
        const double sigY = pitch.second / std::sqrt(12.0);

        // A valid crossing was attempted: seed its sidecar record now (spec §2 -
        // one hitInfo entry per valid crossing, whether or not a hit is accepted).
        smartpixels::SmartPixelsRefitHitInfo hi;
        hi.layer = static_cast<uint8_t>(layer);
        hi.detId = cx.detId;

        // ---- window-collect + classify digis (readout order, FPGA truncation) ----
        struct HitCand {
          double x = 0., y = 0.;
          double cotA = -999., cotB = -999.;
          double sigA = 0., sigB = 0.;
          bool hasA = false, hasB = false;  // per-angle validity (a real payload may be alpha-only)
          double sel = 0.;
          int8_t cls = 2;                   // simlink class: 0 sameTP, 1 otherTP, 2 noise
          double parCotA = -999., parCotB = -999.;  // parent local angles (truth-only)
        };
        std::vector<HitCand> cands;
        bool windowTruncated = false;
        const auto digiSet = drDigis->find(cx.detId);
        if (digiSet != drDigis->end()) {
          std::map<unsigned int, const PixelDigiSimLink*> linkByChannel;
          const auto linkSet = drSimlinks->find(cx.detId);
          if (linkSet != drSimlinks->end()) {
            for (const auto& lk : *linkSet) {
              auto it = linkByChannel.find(lk.channel());
              if (it == linkByChannel.end() || it->second->fraction() < lk.fraction())
                linkByChannel[lk.channel()] = &lk;
            }
          }
          for (const auto& digi : *digiSet) {
            if (static_cast<int>(cands.size()) >= digiRefitMaxHitsPerWindow_) {
              windowTruncated = true;  // combinatorics truncation in readout order (hardware-like)
              break;
            }
            const LocalPoint dlp =
                pixTopo.localPosition(MeasurementPoint(digi.row() + 0.5, digi.column() + 0.5));
            if (std::abs(dlp.x() - cx.local.x()) > digiRefitWindowRPhi_[layer - 1] ||
                std::abs(dlp.y() - cx.local.y()) > digiRefitWindowZ_[layer - 1])
              continue;

            HitCand cand;
            cand.x = dlp.x();
            cand.y = dlp.y();

            const auto lit = linkByChannel.find(digi.channel());
            if (lit != linkByChannel.end()) {
              // Linked digi (same-TP or other-TP): angle from ITS OWN parent's
              // origin-momentum incidence (payload-analyzer convention), smeared
              // with the PixelAV response; validity coin from valid_prob.
              // Class (truth-only): same TP as THIS track vs a different TP.
              cand.cls = (lit->second->eventId() == drMatchedEvtId &&
                          drMatchedSimIds.count(lit->second->SimTrackId()))
                             ? int8_t(0)
                             : int8_t(1);
              const auto mit = drParentMom.find({lit->second->eventId().rawId(), lit->second->SimTrackId()});
              if (mit != drParentMom.end()) {
                const auto& pm = mit->second;
                const LocalVector plv = pixDet->toLocal(GlobalVector(pm.px(), pm.py(), pm.pz()));
                const double ppz = (std::abs(plv.z()) > 1e-9) ? plv.z() : 1e-9;
                const double trueCotA = plv.x() / ppz;
                const double trueCotB = plv.y() / ppz;
                cand.parCotA = trueCotA;
                cand.parCotB = trueCotB;
                // Synthesis contract (spec §3, fused): validity gate + one shift
                // evaluate per angle. Everything stochastic is payload-side HashPRNG
                // -- deterministic per (layer, cotA, cotB, bLocalY) input tuple.
                const std::vector<std::variant<int, double, std::string>> pin = {
                    layer, trueCotA, trueCotB, static_cast<double>(cx.bLocalY)};
                if (corrValidFlat_->evaluate(pin) < corrValidProb_->evaluate(pin)) {
                  cand.sigA = corrAlphaSigma_->evaluate(pin);
                  cand.sigB = corrBetaSigma_->evaluate(pin);
                  const std::vector<std::variant<int, double, std::string>> pinAcc = {
                      layer, trueCotA, trueCotB, static_cast<double>(cx.bLocalY), 1.0};
                  cand.cotA = trueCotA + corrAlphaShift_->evaluate(pinAcc);
                  cand.cotB = trueCotB + corrBetaShift_->evaluate(pinAcc);
                  cand.hasA = (cand.sigA > 0.);
                  cand.hasB = (cand.sigB > 0.);
                }
              }
            } else if (corrNoiseCotAlpha_ && corrNoiseCotBeta_) {
              cand.cls = 2;  // no simlink -> noise
              // No-link (noise) digi with a Stack-B inclusive angle model:
              // inverse-CDF draw (layer, uniform quantile). Without the model,
              // noise digis contribute position only.
              const std::vector<std::variant<int, double, std::string>> na = {
                  layer, CLHEP::RandFlat::shoot(drEngine)};
              const std::vector<std::variant<int, double, std::string>> nb = {
                  layer, CLHEP::RandFlat::shoot(drEngine)};
              cand.cotA = corrNoiseCotAlpha_->evaluate(na);
              cand.cotB = corrNoiseCotBeta_->evaluate(nb);
              const std::vector<std::variant<int, double, std::string>> pin = {
                  layer, cand.cotA, cand.cotB, static_cast<double>(cx.bLocalY)};
              cand.sigA = corrAlphaSigma_->evaluate(pin);
              cand.sigB = corrBetaSigma_->evaluate(pin);
              cand.hasA = (cand.sigA > 0.);
              cand.hasB = (cand.sigB > 0.);
            }

            // Measured-angle grazing clamp (spec v0.4 §6b, LOAD-BEARING): a
            // synthesized measured cotAlpha/cotBeta beyond the physical bound is a
            // near-grazing-parent synthesis breakdown (p_z at its 1e-9 floor ->
            // |cot| up to ~2274). Invalidate THAT ANGLE ONLY (clear hasA/hasB); the
            // hit keeps its position measurement. This removes the entire chi2
            // pathology at source (77/77 gated updates were measurement-driven);
            // chi2UpdateGate remains as the numerical backstop.
            if (cand.hasA && std::abs(cand.cotA) > digiRefitMeasAngleMaxAbs_) {
              cand.hasA = false;
              ++digiRefitClampedMeasAngles_;
            }
            if (cand.hasB && std::abs(cand.cotB) > digiRefitMeasAngleMaxAbs_) {
              cand.hasB = false;
              ++digiRefitClampedMeasAngles_;
            }

            // Selection chi2 against the prediction (position always; angles when enabled+valid).
            double sel = std::pow((cand.x - cx.local.x()) / sigX, 2) +
                         std::pow((cand.y - cx.local.y()) / sigY, 2);
            if (useAlpha && cand.hasA)
              sel += std::pow((cand.cotA - cx.cotAlpha) / cand.sigA, 2);
            if (useBeta && cand.hasB)
              sel += std::pow((cand.cotB - cx.cotBeta) / cand.sigB, 2);
            cand.sel = sel;
            cands.push_back(cand);
          }
        }

        // Crossing-level occupancy (post-truncation) recorded on the hit info.
        hi.windowMult = static_cast<uint16_t>(cands.size());
        if (windowTruncated) {
          hi.flags |= smartpixels::hitflag::kWindowTruncated;
          drAnyWindowTruncated = true;
        }
        drMaxWindowMult = std::max<unsigned>(drMaxWindowMult, cands.size());
        ++drTrackInfo.nCrossings;

        if (cands.empty()) {
          // Valid crossing, no hit: sidecar record kept (selected-hit fields stay
          // at their -999.f / class -1 sentinels), no KF update.
          drHitInfo.push_back(hi);
          continue;
        }
        const auto bestIt = std::min_element(
            cands.begin(), cands.end(), [](const HitCand& p, const HitCand& q) { return p.sel < q.sel; });
        const HitCand& best = *bestIt;
        // selChi2Margin (spec v0.4 §2): runner-up minus best selection chi2 (>=0),
        // a hardware-plausible measure of how unambiguous the hit choice was.
        // Sentinel -999.f when the window held fewer than 2 candidates (no runner-up).
        double drSelChi2Margin = -999.;
        if (cands.size() >= 2) {
          double runnerUp = std::numeric_limits<double>::infinity();
          for (auto it = cands.begin(); it != cands.end(); ++it)
            if (it != bestIt && it->sel < runnerUp)
              runnerUp = it->sel;
          drSelChi2Margin = runnerUp - best.sel;  // >= 0 by construction
        }

        // ---- linearized measurement model h(a) = (localx, localy, cotAlpha, cotBeta) ----
        // Numerical Jacobian against the SAME shared projector (consistency by
        // construction). One-sided differences; if a perturbation falls off the
        // nominal module, try the other side; else the column stays zero
        // (parameter locally unobservable).
        const double h0[4] = {cx.local.x(), cx.local.y(), cx.cotAlpha, cx.cotBeta};
        double H[4][5] = {{0.}};
        constexpr std::array<double, 5> kEps{{1e-6, 1e-5, 1e-5, 1e-3, 1e-3}};
        for (int j = 0; j < 5; ++j) {
          for (const double sgn : {+1., -1.}) {
            ROOT::Math::SVector<double, 5> ap = a;
            ap[j] += sgn * kEps[j];
            const smartpixels::Crossing cp =
                projector_.crossLayer(makeHelix(ap), layer, drField, digiRefitPredAngleMaxAbs_);
            if (!cp.valid || cp.detId != cx.detId)
              continue;
            const double inv = 1.0 / (sgn * kEps[j]);
            H[0][j] = (cp.local.x() - h0[0]) * inv;
            H[1][j] = (cp.local.y() - h0[1]) * inv;
            H[2][j] = (cp.cotAlpha - h0[2]) * inv;
            H[3][j] = (cp.cotBeta - h0[3]) * inv;
            break;
          }
        }
        // jacobianMaxAbs (spec §6b): any Jacobian entry above the bound (or
        // non-finite) makes that PARAMETER column locally unobservable -- zero it
        // exactly like the same-module guard above. Physical |H| stays <~200 on the
        // PU study; near-grazing crossings blow a column up (observed to ~900) or
        // yield NaN/Inf. Column-zeroing keeps the update well-posed for the other
        // parameters. FPGA-fidelity handle: finite-precision hardware bounds |H|.
        for (int j = 0; j < 5; ++j) {
          bool bad = false;
          for (int k = 0; k < 4; ++k)
            if (!std::isfinite(H[k][j]) || std::abs(H[k][j]) > digiRefitJacobianMaxAbs_) {
              bad = true;
              break;
            }
          if (bad) {
            ++digiRefitGatedJacCols_;
            for (int k = 0; k < 4; ++k)
              H[k][j] = 0.;
          }
        }

        // ---- sequential scalar Kalman updates (diagonal R; FPGA-friendly form) ----
        // Per-scalar pull (r/sqrt(S), computed BEFORE the state update) and chi2
        // increment (r^2/S), captured for the sidecar. k: 0=x 1=y 2=alpha 3=beta.
        const ROOT::Math::SVector<double, 5> aLin = a;
        std::array<double, 4> pull{{-999., -999., -999., -999.}};
        std::array<double, 4> chi2inc{{0., 0., 0., 0.}};
        std::array<bool, 4> applied{{false, false, false, false}};
        const auto scalarUpdate = [&](int k, double m, double sig) {
          ROOT::Math::SVector<double, 5> Hrow;
          for (int j = 0; j < 5; ++j)
            Hrow[j] = H[k][j];
          const ROOT::Math::SVector<double, 5> v = C * Hrow;
          const double S = ROOT::Math::Dot(Hrow, v) + sig * sig;
          if (!(S > 0.))
            return;
          const double r = m - h0[k] - ROOT::Math::Dot(Hrow, a - aLin);
          const double chi2candidate = r * r / S;
          // chi2UpdateGate (spec §6b): a scalar update whose r^2/S exceeds the gate
          // is a numerical pathology (near-grazing crossing -> non-physical predicted
          // angle and/or blown-up Jacobian column, innovation up to ~1e3 cot-units).
          // Skip it ENTIRELY: no state/covariance change, no chi2 contribution, pull
          // stays sentinel. This is a numerical-pathology gate, NOT an outlier/quality
          // gate -- physical wrong-hit contamination (r^2/S <= ~2e6 on the PU study) is
          // the BDT's signal and is left untouched. The hit still counts as accepted.
          if (chi2candidate > digiRefitChi2UpdateGate_) {
            ++digiRefitGatedUpdates_;
            return;
          }
          pull[k] = r / std::sqrt(S);
          chi2inc[k] = chi2candidate;
          applied[k] = true;
          const ROOT::Math::SVector<double, 5> K = v / S;
          a += K * r;
          for (int i = 0; i < 5; ++i)
            for (int j = 0; j <= i; ++j)
              C(i, j) -= K[i] * v[j];  // K = v/S -> K_i v_j symmetric
        };
        scalarUpdate(0, best.x, sigX);
        scalarUpdate(1, best.y, sigY);
        if (useAlpha && best.hasA)
          scalarUpdate(2, best.cotA, best.sigA);
        if (useBeta && best.hasB)
          scalarUpdate(3, best.cotB, best.sigB);

        // ---- fill the accepted-hit sidecar record (spec §2) ----
        hi.flags |= smartpixels::hitflag::kHitAccepted;
        if (best.hasA)
          hi.flags |= smartpixels::hitflag::kHasAlpha;
        if (best.hasB)
          hi.flags |= smartpixels::hitflag::kHasBeta;
        hi.resX = static_cast<float>(best.x - cx.local.x());
        hi.resY = static_cast<float>(best.y - cx.local.y());
        hi.cotAlphaMeas = best.hasA ? static_cast<float>(best.cotA) : -999.f;
        hi.cotBetaMeas = best.hasB ? static_cast<float>(best.cotB) : -999.f;
        hi.sigAlpha = best.hasA ? static_cast<float>(best.sigA) : -999.f;
        hi.sigBeta = best.hasB ? static_cast<float>(best.sigB) : -999.f;
        hi.pullX = applied[0] ? static_cast<float>(pull[0]) : -999.f;
        hi.pullY = applied[1] ? static_cast<float>(pull[1]) : -999.f;
        hi.pullAlpha = applied[2] ? static_cast<float>(pull[2]) : -999.f;
        hi.pullBeta = applied[3] ? static_cast<float>(pull[3]) : -999.f;
        // Refit-BDT sum-of-squared-pull features (spec §6a features 5-8).
        if (applied[0]) drSumPullX2 += pull[0] * pull[0];
        if (applied[1]) drSumPullY2 += pull[1] * pull[1];
        if (applied[2]) drSumPullAlpha2 += pull[2] * pull[2];
        if (applied[3]) drSumPullBeta2 += pull[3] * pull[3];
        // rphi terms: x + alpha; rz terms: y + beta (spec §2).
        const double incRPhi = (applied[0] ? chi2inc[0] : 0.) + (applied[2] ? chi2inc[2] : 0.);
        const double incRZ = (applied[1] ? chi2inc[1] : 0.) + (applied[3] ? chi2inc[3] : 0.);
        hi.chi2IncRPhi = static_cast<float>(incRPhi);
        hi.chi2IncRZ = static_cast<float>(incRZ);
        hi.selChi2Margin = (drSelChi2Margin > -900.) ? static_cast<float>(drSelChi2Margin) : -999.f;
        hi.selHitClass = best.cls;
        hi.parCotAlpha = static_cast<float>(best.parCotA);
        hi.parCotBeta = static_cast<float>(best.parCotB);
        drHitInfo.push_back(hi);

        // Compact-word layer bitmask + per-track chi2 totals (spec §2/§3).
        drLayerHitMask |= static_cast<uint8_t>(1u << (layer - 1));
        drTrackInfo.chi2IncRPhiTot =
            (drTrackInfo.chi2IncRPhiTot < -900.f ? 0.f : drTrackInfo.chi2IncRPhiTot) + static_cast<float>(incRPhi);
        drTrackInfo.chi2IncRZTot =
            (drTrackInfo.chi2IncRZTot < -900.f ? 0.f : drTrackInfo.chi2IncRZTot) + static_cast<float>(incRZ);

        ++nAcceptedHits;
        ++nUpdates;
      }  // layer loop

      if (!seedCovOK || nAcceptedHits < digiRefitMinHits_) {
        // 1:1 passthrough fallback (exact copy, original covariance included).
        // Sidecar: status bit0 unset, EMPTY hitInfo, zeros elsewhere (spec §2).
        outputTracks->push_back(*iterL1Track);
        smartpixels::SmartPixelsRefitTrackInfo pt;  // all zero / status bit0 unset
        if (seedCovOK)
          pt.status |= smartpixels::trackstatus::kSeedCovOK;
        if (digiRefitSeedCovMode_ == "parametrized")
          pt.status |= smartpixels::trackstatus::kParametrizedSeed;
        pt.chi2IncRPhiTot = 0.f;
        pt.chi2IncRZTot = 0.f;
        sidecar->trackInfo.push_back(pt);
        sidecar->hitInfo.emplace_back();  // empty per spec
      } else {
        // CMSSW_17 covMatrix variant (Branch C): the refit updates the helix
        // params (INVR/PHI0/TANL/Z0/D0), packs them into the track word below, AND
        // persists the KF-updated 5x5 covariance on the output TTTrack via the
        // 18-arg ctor (enabled by the covMatrix omnibus, cms-sw#51503). The seed
        // covariance is read from the input track's helixCovMat() in trackCov mode.

        // Refit-quality BDT score (spec §6a). When a model is loaded, assemble
        // REFIT_BDT_FEATURES v0 (17 ordered features, all in-flight/full-fidelity)
        // and evaluate; the score REPLACES the trkMVA1 ctor slot on this refit
        // track. Empty model -> keep the input track's trkMVA1 (current behavior).
        // Passthrough tracks are never scored (handled in the fallback branch).
        double refitTrkMVA1 = iterL1Track->trkMVA1();
        if (digiRefitBdt_) {
          const double chi2RPhiTot = (drTrackInfo.chi2IncRPhiTot < -900.f) ? 0. : drTrackInfo.chi2IncRPhiTot;
          const double chi2RZTot = (drTrackInfo.chi2IncRZTot < -900.f) ? 0. : drTrackInfo.chi2IncRZTot;
          std::vector<float> feats = {
              static_cast<float>(drTrackInfo.nCrossings),               // 0
              static_cast<float>(nAcceptedHits),                        // 1
              static_cast<float>(drLayerHitMask),                       // 2
              static_cast<float>(drMaxWindowMult),                      // 3
              static_cast<float>(drAnyWindowTruncated ? 1 : 0),         // 4
              static_cast<float>(drSumPullX2),                          // 5
              static_cast<float>(drSumPullY2),                          // 6
              static_cast<float>(drSumPullAlpha2),                      // 7
              static_cast<float>(drSumPullBeta2),                       // 8
              static_cast<float>(chi2RPhiTot),                          // 9
              static_cast<float>(chi2RZTot),                            // 10
              static_cast<float>(a[0] - aSeed[0]),                      // 11 dRinv
              static_cast<float>(a[1] - aSeed[1]),                      // 12 dPhi
              static_cast<float>(a[2] - aSeed[2]),                      // 13 dTanl
              static_cast<float>(a[3] - aSeed[3]),                      // 14 dZ0
              static_cast<float>(a[4] - aSeed[4]),                      // 15 dD0
              static_cast<float>(iterL1Track->trkMVA1())};              // 16 seedTrkMVA1
          // v1 (spec §6a v1): append the INPUT track's classic-7 TrackQuality hw
          // features (indices 17-23), decoded bit-exactly as trkquality.py does
          // on the raw nano hw columns (getXxxBits() -> two's-complement / raw
          // bin index; nStubs = getStubRefs().size(); nlaymiss over hitPattern()).
          if (digiRefitBdtNFeatures_ == kRefitBdtNFeaturesV1) {
            feats.push_back(static_cast<float>(twosComplement(iterL1Track->getTanlBits(), 16)));   // 17 tanl
            feats.push_back(static_cast<float>(twosComplement(iterL1Track->getZ0Bits(), 12)));     // 18 z0_scaled
            feats.push_back(static_cast<float>(iterL1Track->getBendChi2Bits()));                   // 19 bendchi2_bin
            feats.push_back(static_cast<float>(iterL1Track->getStubRefs().size()));                // 20 nstub
            feats.push_back(static_cast<float>(nlaymissInterior(iterL1Track->hitPattern(), 7)));   // 21 nlaymiss_interior
            feats.push_back(static_cast<float>(iterL1Track->getChi2RPhiBits()));                   // 22 chi2rphi_bin
            feats.push_back(static_cast<float>(iterL1Track->getChi2RZBits()));                     // 23 chi2rz_bin
          }
          refitTrkMVA1 = digiRefitBdt_->decision_function(feats).at(0);
        }

        // Branch C: persist the KF-updated 5x5 helix covariance as the 18th
        // TTTrack ctor arg (cms-sw#51503 covMat overload). C is the running KF
        // covariance after all accepted scalar updates. Order = enum Hpar
        // {INVR,PHI0,TANL,Z0,D0}; float copy into L1Track::CovMat (ErrorF<5>).
        L1Track::CovMat outCov;
        for (int i = 0; i < 5; ++i)
          for (int j = 0; j <= i; ++j)
            outCov(i, j) = static_cast<float>(C(i, j));
        // 18-arg covMat ctor: phiSector/etaSector/chi2BendRed/trackSeedType are
        // ctor args here (no post-hoc setters needed for those). refitTrkMVA1 is
        // the refit-BDT score (Branch-B feature) in the trkMVA1 slot.
        L1Track track = L1Track(a[0],
                                a[1],
                                a[2],
                                a[3],
                                a[4],
                                iterL1Track->chi2XY(),
                                iterL1Track->chi2Z(),
                                refitTrkMVA1,
                                iterL1Track->trkMVA2(),
                                iterL1Track->trkMVA3(),
                                iterL1Track->hitPattern(),
                                iterL1Track->nFitPars(),
                                b_field,
                                iterL1Track->phiSector(),
                                iterL1Track->etaSector(),
                                iterL1Track->chi2BendRed(),
                                iterL1Track->trackSeedType(),
                                outCov);
        track.setStubRefs(iterL1Track->getStubRefs());
        track.setTrackWordBits();
        outputTracks->push_back(track);

        // Finalize the per-track sidecar record (spec §2).
        drTrackInfo.status |= smartpixels::trackstatus::kRefitPerformed;
        drTrackInfo.status |= smartpixels::trackstatus::kSeedCovOK;
        if (digiRefitSeedCovMode_ == "parametrized")
          drTrackInfo.status |= smartpixels::trackstatus::kParametrizedSeed;
        if (drAnyWindowTruncated)
          drTrackInfo.status |= smartpixels::trackstatus::kAnyWindowTruncated;
        drTrackInfo.nAcceptedHits = static_cast<uint8_t>(nAcceptedHits);
        drTrackInfo.nKFUpdates = static_cast<uint8_t>(nUpdates);
        drTrackInfo.layerHitMask = drLayerHitMask;  // popcount == nAcceptedHits (spec v0.1)
        drTrackInfo.maxWindowMult = static_cast<uint16_t>(drMaxWindowMult);
        if (drTrackInfo.chi2IncRPhiTot < -900.f)
          drTrackInfo.chi2IncRPhiTot = 0.f;
        if (drTrackInfo.chi2IncRZTot < -900.f)
          drTrackInfo.chi2IncRZTot = 0.f;
        sidecar->trackInfo.push_back(drTrackInfo);
        sidecar->hitInfo.push_back(std::move(drHitInfo));
      }
    } // end digiRefit path

    // Summary info for track changes
    if (my_tp.isNull() == false) {
      if (my_tp->pt() > 1.95) {
	auto o_trk = *iterL1Track;
	auto n_trk = outputTracks->back();
	if( o_trk.nFitPars() == 4 )
	  track_summary_["npars_4_entries"][tp_track_match]++;
	if( o_trk.nFitPars() == 5 )
	  track_summary_["npars_5_entries"][tp_track_match]++;
	track_summary_["resolution_entries"][tp_track_match]++;
	
	track_summary_["float_pt_diff_sum_orig"][tp_track_match] += (o_trk.momentum().perp() - my_tp->pt());
	track_summary_["float_pt_diff_sum2_orig"][tp_track_match] += std::pow((o_trk.momentum().perp() - my_tp->pt()), 2);
	track_summary_["float_pt_diff_sum_new"][tp_track_match] += (n_trk.momentum().perp() - my_tp->pt());
	track_summary_["float_pt_diff_sum2_new"][tp_track_match] += std::pow((n_trk.momentum().perp() - my_tp->pt()), 2);
	
	track_summary_["float_eta_diff_sum_orig"][tp_track_match] += (o_trk.eta() - my_tp->eta());
	track_summary_["float_eta_diff_sum2_orig"][tp_track_match] += std::pow((o_trk.eta() - my_tp->eta()), 2);
	track_summary_["float_eta_diff_sum_new"][tp_track_match] += (n_trk.eta() - my_tp->eta());
	track_summary_["float_eta_diff_sum2_new"][tp_track_match] += std::pow((n_trk.eta() - my_tp->eta()), 2);

	track_summary_["float_phi_diff_sum_orig"][tp_track_match] += (o_trk.phi() - my_tp->phi());
	track_summary_["float_phi_diff_sum2_orig"][tp_track_match] += std::pow((o_trk.phi() - my_tp->phi()), 2);
	track_summary_["float_phi_diff_sum_new"][tp_track_match] += (n_trk.phi() - my_tp->phi());
	track_summary_["float_phi_diff_sum2_new"][tp_track_match] += std::pow((n_trk.phi() - my_tp->phi()), 2);

	track_summary_["float_z0_diff_sum_orig"][tp_track_match] += (o_trk.z0() - my_tp->z0());
	track_summary_["float_z0_diff_sum2_orig"][tp_track_match] += std::pow((o_trk.z0() - my_tp->z0()), 2);
	track_summary_["float_z0_diff_sum_new"][tp_track_match] += (n_trk.z0() - my_tp->z0());
	track_summary_["float_z0_diff_sum2_new"][tp_track_match] += std::pow((n_trk.z0() - my_tp->z0()), 2);	
	
	track_summary_["float_d0_diff_sum_orig"][tp_track_match] += (o_trk.d0() - my_tp->d0());
	track_summary_["float_d0_diff_sum2_orig"][tp_track_match] += std::pow((o_trk.d0() - my_tp->d0()), 2);
	track_summary_["float_d0_diff_sum_new"][tp_track_match] += (n_trk.d0() - my_tp->d0());
	track_summary_["float_d0_diff_sum2_new"][tp_track_match] += std::pow((n_trk.d0() - my_tp->d0()), 2);

	auto o_trk_pt_diff = (my_tp->charge() * MagConstant * b_field / (o_trk.getRinv() * 100.0) - my_tp->pt());
	auto n_trk_pt_diff = (my_tp->charge() * MagConstant * b_field / (n_trk.getRinv() * 100.0) - my_tp->pt());
	track_summary_["hw_pt_diff_sum_orig"][tp_track_match] += o_trk_pt_diff;
	track_summary_["hw_pt_diff_sum2_orig"][tp_track_match] += std::pow(o_trk_pt_diff, 2);
	track_summary_["hw_pt_diff_sum_new"][tp_track_match] += n_trk_pt_diff;
	track_summary_["hw_pt_diff_sum2_new"][tp_track_match] += std::pow(n_trk_pt_diff, 2);

	auto o_trk_eta_diff = std::asinh(o_trk.getTanl()) - tmp_matchtp_eta;
	auto n_trk_eta_diff = std::asinh(n_trk.getTanl()) - tmp_matchtp_eta;
	track_summary_["hw_eta_diff_sum_orig"][tp_track_match] += o_trk_eta_diff;
	track_summary_["hw_eta_diff_sum2_orig"][tp_track_match] += std::pow(o_trk_eta_diff, 2);
	track_summary_["hw_eta_diff_sum_new"][tp_track_match] += n_trk_eta_diff;
	track_summary_["hw_eta_diff_sum2_new"][tp_track_match] += std::pow(n_trk_eta_diff, 2);

	auto o_trk_phi_diff = o_trk.getPhi() - tmp_matchtp_phi;
	auto n_trk_phi_diff = n_trk.getPhi() - tmp_matchtp_phi;
	track_summary_["hw_phi_diff_sum_orig"][tp_track_match] += o_trk_phi_diff;
	track_summary_["hw_phi_diff_sum2_orig"][tp_track_match] += std::pow(o_trk_phi_diff, 2);
	track_summary_["hw_phi_diff_sum_new"][tp_track_match] += n_trk_phi_diff;
	track_summary_["hw_phi_diff_sum2_new"][tp_track_match] += std::pow(n_trk_phi_diff, 2);

	auto o_trk_z0_diff = o_trk.getZ0() - tmp_matchtp_z0;
	auto n_trk_z0_diff = n_trk.getZ0() - tmp_matchtp_z0;
	track_summary_["hw_z0_diff_sum_orig"][tp_track_match] += o_trk_z0_diff;
	track_summary_["hw_z0_diff_sum2_orig"][tp_track_match] += std::pow(o_trk_z0_diff, 2);
	track_summary_["hw_z0_diff_sum_new"][tp_track_match] += n_trk_z0_diff;
	track_summary_["hw_z0_diff_sum2_new"][tp_track_match] += std::pow(n_trk_z0_diff, 2);	
	
	auto o_trk_d0_diff = o_trk.getD0() - tmp_matchtp_d0;
	auto n_trk_d0_diff = n_trk.getD0() - tmp_matchtp_d0;
	track_summary_["hw_d0_diff_sum_orig"][tp_track_match] += o_trk_d0_diff;
	track_summary_["hw_d0_diff_sum2_orig"][tp_track_match] += std::pow(o_trk_d0_diff, 2);
	track_summary_["hw_d0_diff_sum_new"][tp_track_match] += n_trk_d0_diff;
	track_summary_["hw_d0_diff_sum2_new"][tp_track_match] += std::pow(n_trk_d0_diff, 2);
      }
      else {
	// std::cout << "Skipping matched track resolution for pt below threshold: pt=" << my_tp->pt() << " eta=" << my_tp->eta() << std::endl;
      }
    } // end if (
  }  //end track loop
  //DEBUG
  if (DebugMode && (long unsigned int)this_l1track != outputTracks->size()){
    edm::LogVerbatim("SmartPixelsTrackProducer") << "(L1SmartPixelsTrackProducer.cc) Mismatch!!\n\tthis_l1track counter = " << this_l1track << std::endl
						 << "\n\toutputTracks->size() = " << outputTracks->size() << std::endl;
  
  //DEBUG
    edm::LogVerbatim("SmartPixelsTrackProducer") << "(L1SmartPixelsTrackProducer.cc) Event unmatched_tracks: " << unmatched_tracks << " matched_tracks: " << matched_tracks << std::endl;
  }
  // 1:1 output-sync invariant (spec §1): exactly one output track per input
  // track, and (in digiRefit) sidecar sizes equal to the track count. This
  // contract is what keeps every externally-synced collection (PF/Puppi track
  // refs, DispVertex indices, the truth table) valid row-wise. Fail LOUDLY.
  const size_t nInput = TTTrackHandle->size();
  if (outputTracks->size() != nInput)
    throw cms::Exception("SmartPixelsSyncBroken")
        << "output track count " << outputTracks->size() << " != input track count " << nInput
        << " (1:1 output-sync invariant broken).";
  if (smartPixelsEmulatorMode_ == "digiRefit") {
    if (sidecar->trackInfo.size() != nInput || sidecar->hitInfo.size() != nInput)
      throw cms::Exception("SmartPixelsSyncBroken")
          << "digiRefit sidecar sizes (trackInfo=" << sidecar->trackInfo.size()
          << ", hitInfo=" << sidecar->hitInfo.size() << ") != input track count " << nInput
          << " (1:1 output-sync invariant broken).";
    iEvent.put(std::move(sidecar), outputCollectionName_);
  }

  iEvent.put(std::move(outputTracks), outputCollectionName_);
}  // end of produce()

// FILLDESCRIPTIONS
void L1SmartPixelsTrackProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<int>("MyProcess", 1)->setComment("Process ID");
  desc.add<bool>("DebugMode", false)->setComment("Printout lots of debug statements");
  desc.add<bool>("DebugModeDetailed", false)->setComment("Printout extremely detailed information");
  desc.add<int>("L1Tk_minNStub", 4)->setComment("L1 tracks with >= 4 stubs");
  desc.add<edm::InputTag>("L1TrackInputTag", edm::InputTag("l1tTTTracksFromTrackletEmulation", "Level1TTTracks"))->setComment("TTTrack input");
  desc.add<edm::InputTag>("MCTruthTrackInputTag", edm::InputTag("TTTrackAssociatorFromPixelDigis",  "Level1TTTracks"))->setComment("MCTruth input");
  desc.add<edm::InputTag>("L1StubInputTag", edm::InputTag("TTStubsFromPhase2TrackerDigis", "StubAccepted"));
  desc.add<edm::InputTag>("MCTruthClusterInputTag", edm::InputTag("TTClusterAssociatorFromPixelDigis", "ClusterAccepted"));
  desc.add<edm::InputTag>("MCTruthStubInputTag", edm::InputTag("TTStubAssociatorFromPixelDigis", "StubAccepted"));
  desc.add<edm::InputTag>("TrackingParticleInputTag", edm::InputTag("mix", "MergedTrackTruth"));
  desc.add<std::string>("outputCollectionName", "Level1TTTracksEmulation");
  desc.add<std::string>("smartPixelsEmulatorMode", "passthrough")->setComment("passthrough, passthroughFloat, passthroughHW, trackingParticleTruth, correctionlibRegression, correctionlibTPToySmear");
  desc.add<std::string>("smartPixelsActiveLayers", "0000")->setComment("Active layers of SmartPixels, '0000' is none active, '1111' is all active, '0110' has active 2nd and 3rd layers, and '1000' is only innermost layer active");
  desc.add<edm::FileInPath>("smartPixelsCorrectionSet")->setComment("Name of json file containing the correctionlib set used for smartPixelsEmulatorMode (correctionlibRegression | correctionlibTPToySmear )");
  // --- digiRefit (Tier 2) ---
  desc.add<std::vector<double>>("digiRefitWindowRPhi", std::vector<double>{0.05, 0.17, 0.5, 0.9})
      ->setComment("per-layer r-phi search-window half-widths [cm] (TBPX L1-L4; ~q68 of the 2-5 GeV extrapolation spread)");
  desc.add<std::vector<double>>("digiRefitWindowZ", std::vector<double>{0.45, 0.35, 0.25, 0.2})
      ->setComment("per-layer z search-window half-widths [cm] (TBPX L1-L4)");
  desc.add<int>("digiRefitMinHits", 1)->setComment("min attached IT hits to emit a refit track (else passthrough copy)");
  desc.add<std::string>("digiRefitUseAngles", "alpha")->setComment("which synthesized angles enter the refit: none | alpha | alphaBeta");
  desc.add<int>("digiRefitMaxHitsPerWindow", 8)->setComment("combinatorics truncation: max in-window digis kept per layer (readout order)");
  desc.add<int>("digiRefitMaxKFUpdates", 4)->setComment("max Kalman layer updates applied");
  desc.add<std::string>("digiRefitGainMode", "full")->setComment("full (exact gains) | lut (RESERVED table-driven placeholder)");
  desc.add<double>("digiRefitJacobianMaxAbs", 1.0e4)
      ->setComment("FPGA-fidelity KF guard: any |H[k][j]| above this bound (or non-finite) zeroes that "
                   "Jacobian column (parameter locally unobservable). Physical |H|<~200 on the PU study; "
                   "1e4 leaves the physical+fake bulk untouched and only bounds grazing-crossing divergence.");
  desc.add<double>("digiRefitChi2UpdateGate", 2.0e6)
      ->setComment("FPGA-fidelity KF guard: a scalar update with r^2/S above this gate is skipped entirely "
                   "(no state/cov change, no chi2 contribution, pull sentinel; the hit still counts as accepted). "
                   "Numerical-pathology gate ONLY - the physical wrong-hit-contamination ceiling is ~1.9e6, so "
                   "2e6 removes only the non-physical tail (up to ~1e10) and does NOT gate genuine wrong hits. "
                   "Backstop to the v0.4 grazing-angle clamps.");
  desc.add<double>("digiRefitMeasAngleMaxAbs", 12.0)
      ->setComment("LOAD-BEARING grazing clamp (spec v0.4 §6b): a synthesized measured |cotAlpha|/|cotBeta| "
                   "above this bound invalidates THAT ANGLE (hasAlpha/hasBeta cleared; the hit keeps its "
                   "position). Removes the measurement-driven chi2 pathology at source (77/77 gated updates "
                   "on the reference PU sample); default 12 preserves the O(1-6) wrong-hit signal.");
  desc.add<double>("digiRefitPredAngleMaxAbs", 12.0)
      ->setComment("secondary hygiene grazing clamp (spec v0.4 §6b): a predicted crossing |cotAlpha|/|cotBeta| "
                   "above this bound invalidates the crossing at the projector (no window, no sidecar record). "
                   "Rejects the ~18 non-physical predicted crossings (up to |cotAlpha| 51.8); NOT the gate driver.");
  desc.add<int>("digiRefitSeedNPar", 5)->setComment("seed-track parametrization for the KF: 4 | 5");
  desc.add<std::string>("digiRefitSeedCovMode", "trackCov")
      ->setComment("seed covariance: trackCov (TTTrack helixCovMat, default; needs cms-sw#51503 covMatrix omnibus) | parametrized");
  desc.add<std::vector<double>>("digiRefitParamSigmas", std::vector<double>{1e-4, 1e-3, 2e-3, 0.06, 0.05})
      ->setComment("parametrized-mode seed sigmas (rInv[cm^-1], phi0, tanL, z0[cm], d0[cm])");
  desc.add<std::string>("digiRefitPixelavAngleSet", "")->setComment("PixelAV angle-response correctionlib payload path (REQUIRED for digiRefit)");
  desc.add<std::string>("digiRefitSmarthitFakeSet", "")->setComment("optional Stack B smarthit_fake payload (inclusive noise-angle model)");
  desc.add<std::string>("digiRefitSmarthitTrueSet", "")
      ->setComment("RESERVED: Stack A smarthit_true payload path. NOT consumed by Tier-2 (position from real "
                   "digis, angle from PixelAV; Stack A only characterizes true hits). Kept for a future "
                   "SmartPixels ASIC on-chip readout-inefficiency model (smarthit_true_eff). If non-empty, "
                   "the producer emits a LogWarning (category SmartPixelsStackAUnused) and loads nothing.");
  desc.add<std::string>("digiRefitBdtModel", "")
      ->setComment("optional refit-quality BDT: conifer JSON over REFIT_BDT_FEATURES (spec §6a). The assembly is "
                   "selected by the model's n_features: 17 = v0, 24 = v1 (v0 + the classic-7 TrackQuality hw "
                   "features of the input track). Empty = keep the input trkMVA1. When set, the score replaces the "
                   "refit track's trkMVA1 ctor slot (+ track-word MVA bits); n_features must be 17 or 24 or it throws.");
  desc.add<edm::InputTag>("pixelDigiInputTag", edm::InputTag("simSiPixelDigis", "Pixel"));
  desc.add<edm::InputTag>("pixelDigiSimLinkInputTag", edm::InputTag("simSiPixelDigis", "Pixel"));
  desc.add<edm::InputTag>("simTrackInputTag", edm::InputTag("g4SimHits"));
  descriptions.addWithDefaultLabel(desc);
}
///////////////////////////
// DEFINE THIS AS A PLUG-IN
DEFINE_FWK_MODULE(L1SmartPixelsTrackProducer);
