// Original simulation Author:  Rishi Patel
//
// Rewritting/improvements:  George Karathanasis,
//                           georgios.karathanasis@cern.ch, CU Boulder
//
//         Created:  Wed, 01 Aug 2018 14:01:41 GMT
//         Latest update: Nov 2022 (by GK)
//
// Track jets are clustered in a two-layer process, first by clustering in phi,
// then by clustering in eta. The code proceeds as following: putting all tracks
// in a grid of eta vs phi space, and then cluster them. Finally we merge the cl
// usters when needed. The code is improved to use the same module between emula
// tion and simulation was also improved, with bug fixes and being faster.
// Introduction to object (p10-13):
// https://indico.cern.ch/event/791517/contributions/3341650/attachments/1818736/2973771/TrackBasedAlgos_L1TMadrid_MacDonald.pdf
// New and improved version: https://indico.cern.ch/event/1203796/contributions/5073056/attachments/2519806/4333006/trackjet_emu.pdf

// system include files
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/L1TCorrelator/interface/TkJet.h"
#include "DataFormats/L1TCorrelator/interface/TkJetFwd.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "DataFormats/L1Trigger/interface/VertexWord.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonTopologies/interface/PixelGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/PixelGeomDetType.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "DataFormats/L1Trigger/interface/Vertex.h"

//own headers
#include "L1TrackJetEmulatorProducer.h"

//general
#include <TMath.h>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>

using namespace std;
using namespace edm;
using namespace l1t;
using namespace l1ttrackjet;

class L1TrackJetProducer : public stream::EDProducer<> {
public:
  explicit L1TrackJetProducer(const ParameterSet &);
  ~L1TrackJetProducer() override;
  typedef TTTrack<Ref_Phase2TrackerDigi_> L1TTTrackType;
  typedef vector<L1TTTrackType> L1TTTrackCollectionType;
  static void fillDescriptions(ConfigurationDescriptions &descriptions);

private:
  void beginStream(StreamID) override;
  void produce(Event &, const EventSetup &) override;
  void endStream() override;

  // ----------member data ---------------------------

  edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> tTopoToken_;
  const EDGetTokenT<vector<TTTrack<Ref_Phase2TrackerDigi_>>> trackToken_;
  const EDGetTokenT<l1t::VertexWordCollection> PVtxToken_;
  vector<Ptr<L1TTTrackType>> L1TrkPtrs_;
  vector<int> tdtrk_;
  float trkZMax_;
  float trkPtMax_;
  float trkPtMin_;
  float trkEtaMax_;
  float nStubs4PromptChi2_;
  float nStubs5PromptChi2_;
  float nStubs4PromptBend_;
  float nStubs5PromptBend_;
  int trkNPSStubMin_;
  int lowpTJetMinTrackMultiplicity_;
  float lowpTJetThreshold_;
  int highpTJetMinTrackMultiplicity_;
  float highpTJetThreshold_;
  int zBins_;
  int etaBins_;
  int phiBins_;
  double minTrkJetpT_;
  float zStep_;
  float etaStep_;
  float phiStep_;
  bool displaced_;
  float d0CutNStubs4_;
  float d0CutNStubs5_;
  float nStubs4DisplacedChi2_;
  float nStubs5DisplacedChi2_;
  float nStubs4DisplacedBend_;
  float nStubs5DisplacedBend_;
  int nDisplacedTracks_;
  float dzPVTrk_;
};

L1TrackJetProducer::L1TrackJetProducer(const ParameterSet &iConfig)
    : tTopoToken_(esConsumes<TrackerTopology, TrackerTopologyRcd>(edm::ESInputTag("", ""))),
      trackToken_(consumes<vector<TTTrack<Ref_Phase2TrackerDigi_>>>(iConfig.getParameter<InputTag>("L1TrackInputTag"))),
      PVtxToken_(consumes<l1t::VertexWordCollection>(iConfig.getParameter<InputTag>("VertexInputTag"))),
      trkZMax_(iConfig.getParameter<double>("trk_zMax")),
      trkPtMax_(iConfig.getParameter<double>("trk_ptMax")),
      trkPtMin_(iConfig.getParameter<double>("trk_ptMin")),
      trkEtaMax_(iConfig.getParameter<double>("trk_etaMax")),
      nStubs4PromptChi2_(iConfig.getParameter<double>("nStubs4PromptChi2")),
      nStubs5PromptChi2_(iConfig.getParameter<double>("nStubs5PromptChi2")),
      nStubs4PromptBend_(iConfig.getParameter<double>("nStubs4PromptBend")),
      nStubs5PromptBend_(iConfig.getParameter<double>("nStubs5PromptBend")),
      trkNPSStubMin_(iConfig.getParameter<int>("trk_nPSStubMin")),
      lowpTJetMinTrackMultiplicity_(iConfig.getParameter<int>("lowpTJetMinTrackMultiplicity")),
      lowpTJetThreshold_(iConfig.getParameter<double>("lowpTJetThreshold")),
      highpTJetMinTrackMultiplicity_(iConfig.getParameter<int>("highpTJetMinTrackMultiplicity")),
      highpTJetThreshold_(iConfig.getParameter<double>("highpTJetThreshold")),
      zBins_(iConfig.getParameter<int>("zBins")),
      etaBins_(iConfig.getParameter<int>("etaBins")),
      phiBins_(iConfig.getParameter<int>("phiBins")),
      minTrkJetpT_(iConfig.getParameter<double>("minTrkJetpT")),
      displaced_(iConfig.getParameter<bool>("displaced")),
      d0CutNStubs4_(iConfig.getParameter<double>("d0_cutNStubs4")),
      d0CutNStubs5_(iConfig.getParameter<double>("d0_cutNStubs5")),
      nStubs4DisplacedChi2_(iConfig.getParameter<double>("nStubs4DisplacedChi2")),
      nStubs5DisplacedChi2_(iConfig.getParameter<double>("nStubs5DisplacedChi2")),
      nStubs4DisplacedBend_(iConfig.getParameter<double>("nStubs4DisplacedBend")),
      nStubs5DisplacedBend_(iConfig.getParameter<double>("nStubs5DisplacedBend")),
      nDisplacedTracks_(iConfig.getParameter<int>("nDisplacedTracks")),
      dzPVTrk_(iConfig.getParameter<double>("MaxDzTrackPV")) {
  zStep_ = 2.0 * trkZMax_ / (zBins_ + 1);  // added +1 in denom
  etaStep_ = 2.0 * trkEtaMax_ / etaBins_;  //etaStep is the width of an etabin
  phiStep_ = 2 * M_PI / phiBins_;          ////phiStep is the width of a phibin

  if (displaced_)
    produces<TkJetCollection>("L1TrackJetsExtended");
  else
    produces<TkJetCollection>("L1TrackJets");
}

L1TrackJetProducer::~L1TrackJetProducer() {}

void L1TrackJetProducer::produce(Event &iEvent, const EventSetup &iSetup) {
  unique_ptr<TkJetCollection> L1L1TrackJetProducer(new TkJetCollection);

  // Read inputs
  const TrackerTopology &tTopo = iSetup.getData(tTopoToken_);

  edm::Handle<vector<TTTrack<Ref_Phase2TrackerDigi_>>> TTTrackHandle;
  iEvent.getByToken(trackToken_, TTTrackHandle);

  edm::Handle<l1t::VertexWordCollection> PVtx;
  iEvent.getByToken(PVtxToken_, PVtx);
  float PVz = (PVtx->at(0)).z0();

  L1TrkPtrs_.clear();
  tdtrk_.clear();

  // track selection
  for (unsigned int this_l1track = 0; this_l1track < TTTrackHandle->size(); this_l1track++) {
    edm::Ptr<L1TTTrackType> trkPtr(TTTrackHandle, this_l1track);
    float trk_pt = trkPtr->momentum().perp();
    int trk_nstubs = (int)trkPtr->getStubRefs().size();
    float trk_chi2dof = trkPtr->chi2Red();
    float trk_d0 = trkPtr->d0();
    float trk_bendchi2 = trkPtr->stubPtConsistency();

    int trk_nPS = 0;
    for (int istub = 0; istub < trk_nstubs; istub++) {  // loop over the stubs
      DetId detId(trkPtr->getStubRefs().at(istub)->getDetId());
      if (detId.det() == DetId::Detector::Tracker) {
        if ((detId.subdetId() == StripSubdetector::TOB && tTopo.tobLayer(detId) <= 3) ||
            (detId.subdetId() == StripSubdetector::TID && tTopo.tidRing(detId) <= 9))
          trk_nPS++;
      }
    }
    // PS track cut
    if (trk_nPS < trkNPSStubMin_)
      continue;
    if (!TrackQualitySelection(trk_nstubs,
                               trk_chi2dof,
                               trk_bendchi2,
                               nStubs4PromptBend_,
                               nStubs5PromptBend_,
                               nStubs4PromptChi2_,
                               nStubs5PromptChi2_,
                               nStubs4DisplacedBend_,
                               nStubs5DisplacedBend_,
                               nStubs4DisplacedChi2_,
                               nStubs5DisplacedChi2_,
                               displaced_))
      continue;
    if (std::abs(PVz - trkPtr->z0()) > dzPVTrk_ && dzPVTrk_ > 0)
      continue;
    if (std::abs(trkPtr->z0()) > trkZMax_)
      continue;
    if (std::abs(trkPtr->momentum().eta()) > trkEtaMax_)
      continue;
    if (trk_pt < trkPtMin_)
      continue;
    L1TrkPtrs_.push_back(trkPtr);

    if ((std::abs(trk_d0) > d0CutNStubs5_ && trk_nstubs >= 5 && d0CutNStubs5_ >= 0) ||
        (trk_nstubs == 4 && std::abs(trk_d0) > d0CutNStubs4_ && d0CutNStubs4_ >= 0))
      tdtrk_.push_back(1);  //displaced track
    else
      tdtrk_.push_back(0);  // not displaced track
  }

  // if no tracks pass selection return empty containers
  if (L1TrkPtrs_.empty()) {
    if (displaced_)
      iEvent.put(std::move(L1L1TrackJetProducer), "L1TrackJetsExtended");
    else
      iEvent.put(std::move(L1L1TrackJetProducer), "L1TrackJets");
    return;
  }

  MaxZBin mzb;
  mzb.znum = -1;
  mzb.nclust = -1;
  mzb.ht = -1;

  // create 2D grid of eta phi bins
  EtaPhiBin epbins_default[phiBins_][etaBins_];
  float phi = -1.0 * M_PI;
  float eta = -1.0 * trkEtaMax_;
  for (int i = 0; i < phiBins_; ++i) {
    eta = -1.0 * trkEtaMax_;
    for (int j = 0; j < etaBins_; ++j) {
      epbins_default[i][j].phi = (phi + (phi + phiStep_)) / 2.0;
      epbins_default[i][j].eta = (eta + (eta + etaStep_)) / 2.0;
      eta += etaStep_;
    }  // for each etabin
    phi += phiStep_;
  }  // for each phibin (finished creating bins)

  // create z-bins (might be useful for displaced if we run w/o dz cut)
  std::vector<float> zmins, zmaxs;
  for (int zbin = 0; zbin < zBins_; zbin++) {
    zmins.push_back(-1.0 * trkZMax_ + zStep_ * zbin);
    zmaxs.push_back(-1.0 * trkZMax_ + zStep_ * zbin + zStep_ * 2);
  }

  // create vectors of clusters
  std::vector<std::vector<EtaPhiBin>> L1clusters;
  L1clusters.reserve(phiBins_);
  std::vector<EtaPhiBin> L2clusters;

  // default (empty) grid
  for (int i = 0; i < phiBins_; ++i) {
    for (int j = 0; j < etaBins_; ++j) {
      epbins_default[i][j].pTtot = 0;
      epbins_default[i][j].used = false;
      epbins_default[i][j].ntracks = 0;
      epbins_default[i][j].nxtracks = 0;
      epbins_default[i][j].trackidx.clear();
    }
  }

  for (unsigned int zbin = 0; zbin < zmins.size(); ++zbin) {
    // initialize grid
    float zmin = zmins[zbin];
    float zmax = zmaxs[zbin];
    EtaPhiBin epbins[phiBins_][etaBins_];
    std::copy(&epbins_default[0][0], &epbins_default[0][0] + phiBins_ * etaBins_, &epbins[0][0]);

    //clear cluster containers
    L1clusters.clear();
    L2clusters.clear();

    // fill grid with tracks
    for (unsigned int k = 0; k < L1TrkPtrs_.size(); ++k) {
      float trkpt = L1TrkPtrs_[k]->momentum().perp();
      float trketa = L1TrkPtrs_[k]->momentum().eta();
      float trkphi = L1TrkPtrs_[k]->momentum().phi();
      float trkZ = L1TrkPtrs_[k]->z0();
      if (zmax < trkZ)
        continue;
      if (zbin == 0) {
        if (zmin > trkZ)
          continue;
      } else {
        if (zmin >= trkZ)
          continue;
      }

      for (int i = 0; i < phiBins_; ++i) {
        for (int j = 0; j < etaBins_; ++j) {
          float eta_min = epbins[i][j].eta - etaStep_ / 2.0;  //eta min
          float eta_max = epbins[i][j].eta + etaStep_ / 2.0;  //eta max
          float phi_min = epbins[i][j].phi - phiStep_ / 2.0;  //phi min
          float phi_max = epbins[i][j].phi + phiStep_ / 2.0;  //phi max

          if ((trketa < eta_min) || (trketa > eta_max) || (trkphi < phi_min) || (trkphi > phi_max))
            continue;
          if (trkpt < trkPtMax_)
            epbins[i][j].pTtot += trkpt;
          else
            epbins[i][j].pTtot += trkPtMax_;
          epbins[i][j].nxtracks += tdtrk_[k];
          epbins[i][j].trackidx.push_back(k);
          ++epbins[i][j].ntracks;
        }  // for each phibin
      }    // for each phibin
    }      //end loop over tracks

    // cluster tracks in eta (first layer) using grid
    for (int phibin = 0; phibin < phiBins_; ++phibin) {
      L1clusters.push_back(L1_clustering<EtaPhiBin, float, float, float>(epbins[phibin], etaBins_, etaStep_));
    }

    // second layer clustering in phi for using eta clusters
    L2clusters = L2_clustering<EtaPhiBin, float, float, float>(L1clusters, phiBins_, phiStep_, etaStep_);

    // sum all cluster pTs in this zbin to find max
    float sum_pt = 0;
    for (unsigned int k = 0; k < L2clusters.size(); ++k) {
      if (L2clusters[k].pTtot > lowpTJetThreshold_ && L2clusters[k].ntracks < lowpTJetMinTrackMultiplicity_)
        continue;
      if (L2clusters[k].pTtot > highpTJetThreshold_ && L2clusters[k].ntracks < highpTJetMinTrackMultiplicity_)
        continue;
      if (L2clusters[k].pTtot > minTrkJetpT_)
        sum_pt += L2clusters[k].pTtot;
    }

    if (sum_pt < mzb.ht)
      continue;

    // if ht is larger than previous max, this is the new vertex zbin
    mzb.ht = sum_pt;
    mzb.znum = zbin;
    mzb.clusters = L2clusters;
    mzb.nclust = L2clusters.size();
    mzb.zbincenter = (zmin + zmax) / 2.0;
  }  //zbin loop

  // output
  vector<Ptr<L1TTTrackType>> L1TrackAssocJet;
  for (unsigned int j = 0; j < mzb.clusters.size(); ++j) {
    float jetEta = mzb.clusters[j].eta;
    float jetPhi = mzb.clusters[j].phi;
    float jetPt = mzb.clusters[j].pTtot;
    float jetPx = jetPt * cos(jetPhi);
    float jetPy = jetPt * sin(jetPhi);
    float jetPz = jetPt * sinh(jetEta);
    float jetP = jetPt * cosh(jetEta);
    int totalDisptrk = mzb.clusters[j].nxtracks;
    bool isDispJet = false;
    if (totalDisptrk > nDisplacedTracks_ || totalDisptrk == nDisplacedTracks_)
      isDispJet = true;

    math::XYZTLorentzVector jetP4(jetPx, jetPy, jetPz, jetP);
    L1TrackAssocJet.clear();
    for (unsigned int itrk = 0; itrk < mzb.clusters[j].trackidx.size(); itrk++)
      L1TrackAssocJet.push_back(L1TrkPtrs_[mzb.clusters[j].trackidx[itrk]]);

    TkJet trkJet(jetP4, L1TrackAssocJet, mzb.zbincenter, mzb.clusters[j].ntracks, 0, totalDisptrk, 0, isDispJet);

    L1L1TrackJetProducer->push_back(trkJet);
  }

  std::sort(
      L1L1TrackJetProducer->begin(), L1L1TrackJetProducer->end(), [](auto &a, auto &b) { return a.pt() > b.pt(); });
  if (displaced_)
    iEvent.put(std::move(L1L1TrackJetProducer), "L1TrackJetsExtended");
  else
    iEvent.put(std::move(L1L1TrackJetProducer), "L1TrackJets");
}

void L1TrackJetProducer::beginStream(StreamID) {}

void L1TrackJetProducer::endStream() {}

void L1TrackJetProducer::fillDescriptions(ConfigurationDescriptions &descriptions) {
  ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1TrackJetProducer);
