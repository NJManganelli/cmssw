// -*- C++ -*-
//
// SmartPixelsClusterCensusAnalyzer
//
// How many IT pixel clusters exist per TBPX layer, and how many survive a
// truth-pT requirement on the cluster's DOMINANT charge contributor?
//
// This is the input to two decisions that cannot be made without it:
//   1. whether an untruncated per-cluster nano table is affordable at PU200, and
//      at what pT threshold (a threshold is a physics choice with a size
//      consequence, so it has to be measured before it is picked); and
//   2. what the refit's match stage is actually up against. windowMult in the
//      sidecar is POST-truncation (maxHitsPerWindow) and post-window, so it
//      cannot answer "how many clusters exist to be confused with".
//
// Truth assignment is DELIBERATELY IDENTICAL to L1SmartPixelsTrackProducer's:
// charge share per (eventId, SimTrackId) over the cluster's pixels via
// PixelDigi::pixelToChannel, dominant contributor wins. If the two ever diverge
// this census stops describing what the refit sees, so the logic is duplicated
// rather than approximated.
//
// The pT is the dominant contributor's PARENT momentum from
// buildParentMomentumMap (TPs cover signal + pileup; the signal-only SimTrack
// container backstops pruned TPs) -- the same map the producer uses for angles.

#include <algorithm>
#include <array>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <map>
#include <string>
#include <vector>

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "Geometry/CommonTopologies/interface/PixelGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "L1Trigger/Phase3SmartPixels/interface/SmartPixelsParentMap.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/TrackerDigiSimLink/interface/PixelDigiSimLink.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"

class SmartPixelsClusterCensusAnalyzer : public edm::one::EDAnalyzer<> {
public:
  explicit SmartPixelsClusterCensusAnalyzer(const edm::ParameterSet&);
  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  const edm::EDGetTokenT<SiPixelRecHitCollection> recHitToken_;
  const edm::EDGetTokenT<edm::DetSetVector<PixelDigiSimLink>> simLinkToken_;
  const edm::EDGetTokenT<std::vector<TrackingParticle>> tpToken_;
  const edm::EDGetTokenT<edm::SimTrackContainer> simTrackToken_;
  const edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> topoToken_;
  const edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> geomToken_;
  const std::vector<double> ptThresholds_;

  static constexpr int kNLayers = 4;
  long nEvents_ = 0;
  std::array<long, kNLayers> nTotal_{};     // all clusters
  std::array<long, kNLayers> nNoLink_{};    // no simlink on any pixel -> noise
  std::array<long, kNLayers> nNoParent_{};  // linked but parent momentum unknown
  // [layer][threshold] surviving counts
  std::vector<std::array<long, kNLayers>> nPass_;
  // charge-share bookkeeping, to show how often the "dominant" contributor is marginal
  std::array<long, kNLayers> nShared_{};    // dominant contributor holds < 90% of charge
  // Per-MODULE occupancy. This, not the per-layer total, is what the refit's
  // candidate loop faces: it only ever examines clusters on the crossed module.
  std::array<long, kNLayers> nModulesOcc_{};        // modules with >=1 cluster
  std::array<long, kNLayers> nClusOnOcc_{};         // clusters on those modules
  std::array<long, kNLayers> nModulesGe4_{};        // modules with >=4 clusters
  std::array<std::vector<int>, kNLayers> perModule_;  // for quantiles
  // Geometry census (filled once): what a "module" physically is.
  bool geomDone_ = false;
  std::array<long, kNLayers> geomModules_{};
  std::array<int, kNLayers> geomRows_{}, geomCols_{};
  std::array<int, kNLayers> geomRocsX_{}, geomRocsY_{}, geomRowsPerRoc_{}, geomColsPerRoc_{};
  std::array<float, kNLayers> geomPitchX_{}, geomPitchY_{};
};

SmartPixelsClusterCensusAnalyzer::SmartPixelsClusterCensusAnalyzer(const edm::ParameterSet& iConfig)
    : recHitToken_(consumes<SiPixelRecHitCollection>(iConfig.getParameter<edm::InputTag>("pixelRecHitInputTag"))),
      simLinkToken_(
          consumes<edm::DetSetVector<PixelDigiSimLink>>(iConfig.getParameter<edm::InputTag>("pixelDigiSimLinkInputTag"))),
      tpToken_(consumes<std::vector<TrackingParticle>>(iConfig.getParameter<edm::InputTag>("trackingParticleInputTag"))),
      simTrackToken_(consumes<edm::SimTrackContainer>(iConfig.getParameter<edm::InputTag>("simTrackInputTag"))),
      topoToken_(esConsumes()),
      geomToken_(esConsumes()),
      ptThresholds_(iConfig.getParameter<std::vector<double>>("ptThresholds")) {
  nPass_.resize(ptThresholds_.size());
}

void SmartPixelsClusterCensusAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  const auto& topo = iSetup.getData(topoToken_);
  if (!geomDone_) {
    geomDone_ = true;
    // What IS a "module" here: one PixelGeomDetUnit, i.e. one DetId, i.e. the unit
    // the candidate loop scans. Recorded with hard numbers so "40 clusters per
    // module" can be turned into an occupancy per pixel and per ROC.
    const auto& geom = iSetup.getData(geomToken_);
    for (const auto* det : geom.detUnits()) {
      const DetId did = det->geographicalId();
      if (did.subdetId() != PixelSubdetector::PixelBarrel)
        continue;
      const unsigned lay = topo.pxbLayer(did);
      if (lay < 1 || lay > kNLayers)
        continue;
      const auto* pdu = dynamic_cast<const PixelGeomDetUnit*>(det);
      if (pdu == nullptr)
        continue;
      const PixelTopology& pt = pdu->specificTopology();
      ++geomModules_[lay - 1];
      geomRows_[lay - 1] = pt.nrows();
      geomCols_[lay - 1] = pt.ncolumns();
      geomRocsX_[lay - 1] = pt.rocsX();
      geomRocsY_[lay - 1] = pt.rocsY();
      geomRowsPerRoc_[lay - 1] = pt.rowsperroc();
      geomColsPerRoc_[lay - 1] = pt.colsperroc();
      geomPitchX_[lay - 1] = pt.pitch().first;
      geomPitchY_[lay - 1] = pt.pitch().second;
    }
  }

  edm::Handle<SiPixelRecHitCollection> recHits;
  iEvent.getByToken(recHitToken_, recHits);
  edm::Handle<edm::DetSetVector<PixelDigiSimLink>> simLinks;
  iEvent.getByToken(simLinkToken_, simLinks);
  edm::Handle<std::vector<TrackingParticle>> tps;
  iEvent.getByToken(tpToken_, tps);
  edm::Handle<edm::SimTrackContainer> simTracks;
  iEvent.getByToken(simTrackToken_, simTracks);

  const smartpixels::ParentMomentumMap parentMom =
      smartpixels::buildParentMomentumMap(*tps, simTracks.product());
  ++nEvents_;

  for (const auto& dsv : *recHits) {
    const DetId detId(dsv.detId());
    if (detId.subdetId() != PixelSubdetector::PixelBarrel)
      continue;
    const unsigned layer = topo.pxbLayer(detId);
    if (layer < 1 || layer > kNLayers)
      continue;
    // channel -> best simlink for this module (highest fraction wins), exactly as
    // the producer builds it.
    std::map<unsigned int, const PixelDigiSimLink*> linkByChannel;
    const auto dsl = simLinks->find(detId);
    if (dsl != simLinks->end()) {
      for (const auto& lk : *dsl) {
        auto it = linkByChannel.find(lk.channel());
        if (it == linkByChannel.end() || it->second->fraction() < lk.fraction())
          linkByChannel[lk.channel()] = &lk;
      }
    }

    int nOnThisModule = 0;
    for (const auto& rh : dsv) {
      const SiPixelCluster* cl = rh.cluster().isNonnull() ? &(*rh.cluster()) : nullptr;
      if (cl == nullptr)
        continue;
      ++nTotal_[layer - 1];
      ++nOnThisModule;

      std::map<std::pair<uint32_t, unsigned int>, double> qByTp;
      double qTot = 0.;
      for (const auto& px : cl->pixels()) {
        qTot += px.adc;
        const auto lit =
            linkByChannel.find(PixelDigi::pixelToChannel(static_cast<int>(px.x), static_cast<int>(px.y)));
        if (lit == linkByChannel.end())
          continue;
        qByTp[std::make_pair(lit->second->eventId().rawId(), lit->second->SimTrackId())] += px.adc;
      }
      if (qByTp.empty()) {
        ++nNoLink_[layer - 1];
        continue;
      }
      double qDom = 0.;
      std::pair<uint32_t, unsigned int> domKey;
      for (const auto& kv : qByTp)
        if (kv.second > qDom) {
          qDom = kv.second;
          domKey = kv.first;
        }
      if (qTot > 0. && qDom / qTot < 0.9)
        ++nShared_[layer - 1];

      const auto mit = parentMom.find(domKey);
      if (mit == parentMom.end()) {
        ++nNoParent_[layer - 1];
        continue;
      }
      const double pt = std::hypot(mit->second.px(), mit->second.py());
      for (size_t t = 0; t < ptThresholds_.size(); ++t)
        if (pt > ptThresholds_[t])
          ++nPass_[t][layer - 1];
    }
    if (nOnThisModule > 0) {
      ++nModulesOcc_[layer - 1];
      nClusOnOcc_[layer - 1] += nOnThisModule;
      if (nOnThisModule >= 4)
        ++nModulesGe4_[layer - 1];
      perModule_[layer - 1].push_back(nOnThisModule);
    }
  }
}

void SmartPixelsClusterCensusAnalyzer::endJob() {
  if (nEvents_ == 0)
    return;
  const double inv = 1.0 / static_cast<double>(nEvents_);
  std::ostringstream os;
  os << "\n=== SmartPixels cluster census: " << nEvents_ << " events ===\n";
  os << "per-event mean cluster count per TBPX layer, and survival under a truth-pT\n"
        "requirement on the cluster's DOMINANT charge contributor\n\n";
  os << "  layer      all    noLink   noParent    shared";
  for (double th : ptThresholds_)
    os << "     >" << th << "GeV";
  os << "\n";

  std::array<long, kNLayers> dummy{};
  long allTot = 0, allNoLink = 0, allShared = 0;
  std::vector<long> passTot(ptThresholds_.size(), 0);
  for (int l = 0; l < kNLayers; ++l) {
    os << "    L" << (l + 1) << "  " << std::fixed << std::setprecision(1) << std::setw(8) << nTotal_[l] * inv
       << std::setw(10) << nNoLink_[l] * inv << std::setw(11) << nNoParent_[l] * inv << std::setw(10)
       << nShared_[l] * inv;
    for (size_t t = 0; t < ptThresholds_.size(); ++t) {
      os << std::setw(12) << nPass_[t][l] * inv;
      passTot[t] += nPass_[t][l];
    }
    os << "\n";
    allTot += nTotal_[l];
    allNoLink += nNoLink_[l];
    allShared += nShared_[l];
    (void)dummy;
  }
  os << "  total  " << std::setw(8) << allTot * inv << std::setw(10) << allNoLink * inv << std::setw(11) << " "
     << std::setw(10) << allShared * inv;
  for (size_t t = 0; t < ptThresholds_.size(); ++t)
    os << std::setw(12) << passTot[t] * inv;
  os << "\n\n  retained fraction of ALL clusters:";
  for (size_t t = 0; t < ptThresholds_.size(); ++t)
    os << "   >" << ptThresholds_[t] << "GeV " << std::setprecision(3)
       << (allTot ? double(passTot[t]) / double(allTot) : 0.0);
  os << "\n\n  per-OCCUPIED-MODULE occupancy (what the refit candidate loop actually faces,\n"
        "  since it only examines clusters on the crossed module):\n";
  os << "  layer   modules/ev   clus/occ.module   median   p95   max   frac modules >=4\n";
  for (int l = 0; l < kNLayers; ++l) {
    auto v = perModule_[l];
    if (v.empty())
      continue;
    std::sort(v.begin(), v.end());
    const double med = v[v.size() / 2];
    const double p95 = v[static_cast<size_t>(0.95 * (v.size() - 1))];
    os << "    L" << (l + 1) << std::setw(11) << std::setprecision(1) << nModulesOcc_[l] * inv
       << std::setw(18) << std::setprecision(2)
       << (nModulesOcc_[l] ? double(nClusOnOcc_[l]) / double(nModulesOcc_[l]) : 0.0)
       << std::setw(9) << std::setprecision(0) << med << std::setw(6) << p95 << std::setw(6) << v.back()
       << std::setw(18) << std::setprecision(3)
       << (nModulesOcc_[l] ? double(nModulesGe4_[l]) / double(nModulesOcc_[l]) : 0.0) << "\n";
  }
  os << "\n  GEOMETRY: what a \"module\" is. One module == one PixelGeomDetUnit == one\n"
        "  DetId == the unit the refit candidate loop scans. A module is subdivided into\n"
        "  ROCs (readout chips); there is no separate 'sensor' DetId in CMSSW -- the\n"
        "  sensor and its ROC array are one detUnit, addressed as a single pixel matrix.\n";
  os << "  layer  modules   rows x cols     pixels/module   ROCs(x,y)   rows,cols per ROC"
        "   pitch x,y [um]   occupancy\n";
  for (int l = 0; l < kNLayers; ++l) {
    if (geomModules_[l] == 0)
      continue;
    const long pix = static_cast<long>(geomRows_[l]) * geomCols_[l];
    const double occ = (nModulesOcc_[l] && pix)
                           ? (double(nClusOnOcc_[l]) / double(nModulesOcc_[l])) / double(pix) : 0.0;
    os << "    L" << (l + 1) << std::setw(9) << geomModules_[l] << std::setw(10) << geomRows_[l] << " x "
       << geomCols_[l] << std::setw(14) << pix << std::setw(12) << geomRocsX_[l] << "," << geomRocsY_[l]
       << std::setw(14) << geomRowsPerRoc_[l] << "," << geomColsPerRoc_[l] << std::setw(14)
       << std::setprecision(1) << geomPitchX_[l] * 1e4 << "," << geomPitchY_[l] * 1e4
       << std::setw(14) << std::scientific << std::setprecision(2) << occ << std::fixed << "\n";
  }
  os << "  (occupancy = clusters per occupied module / pixels per module)\n";
  os << "\nNOTE 'noLink' clusters carry no simlink on any pixel (noise-like) and can\n"
        "never pass a truth-pT cut; 'noParent' are linked but their parent momentum is\n"
        "absent from the TP+SimTrack map. Both are counted in 'all' and excluded from\n"
        "every threshold, so a threshold column is a LOWER bound on what a real\n"
        "(non-truth) preselection would keep.\n";
  edm::LogSystem("SmartPixelsClusterCensus") << os.str();
}

void SmartPixelsClusterCensusAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("pixelRecHitInputTag", edm::InputTag("spixPixelRecHits"));
  desc.add<edm::InputTag>("pixelDigiSimLinkInputTag", edm::InputTag("simSiPixelDigis", "Pixel"));
  desc.add<edm::InputTag>("trackingParticleInputTag", edm::InputTag("mix", "MergedTrackTruth"));
  desc.add<edm::InputTag>("simTrackInputTag", edm::InputTag("g4SimHits"));
  desc.add<std::vector<double>>("ptThresholds", {0.5, 1.0, 1.5, 2.0})
      ->setComment("truth-pT thresholds applied to the dominant charge contributor's parent");
  descriptions.add("smartPixelsClusterCensusAnalyzer", desc);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(SmartPixelsClusterCensusAnalyzer);
