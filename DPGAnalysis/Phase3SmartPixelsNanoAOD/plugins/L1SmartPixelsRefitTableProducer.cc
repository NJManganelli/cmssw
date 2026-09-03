// Nano adapter for the SmartPixels refit sidecar (spec §4.4;
// L1Trigger/Phase3SmartPixels/doc/RefitSidecarSpec.md §4). Consumes a refit
// TTTrack collection and its 1:1 row-synced smartpixels::SmartPixelsRefitSidecar
// (same module, same instance label) and emits two nanoaod::FlatTables:
//
//   "hit"  : per-crossing LINK table (one row per layer-crossing record across
//            ALL tracks), carrying trackIdx (index of the owning track in the
//            variant track table -- the L1SC4NGJetCands link-table pattern) plus
//            the per-hit sidecar fields. Residuals and pulls stay at full float
//            precision; convenience angles/chi2 are reduced precision.
//   "trk"  : EXTENSION table on the variant track table (extension=true, SAME
//            name and length as that table) carrying the per-track trackInfo
//            summary + the materialized compact word (packCompactWord).
//
// Row-count invariants (asserted, throw SmartPixelsSyncBroken):
//   sidecar.trackInfo.size() == sidecar.hitInfo.size() == tracks.size()
//   sum_i sidecar.hitInfo[i].size() == nHitRows
// Truth (withGen tiers): no work here -- the reference-track truth table rows
// align 1:1 with every variant's rows by the output-sync invariant (spec §1),
// so analysis reuses trackIdx against the reference truth table.

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"

#include "L1Trigger/Phase3SmartPixels/interface/SmartPixelsRefitSidecar.h"
#include "L1Trigger/Phase3SmartPixels/interface/SmartPixelsTransmittedSubset.h"

#include <cstdint>
#include <vector>

class L1SmartPixelsRefitTableProducer : public edm::stream::EDProducer<> {
public:
  using L1Track = TTTrack<Ref_Phase2TrackerDigi_>;
  using TTTrackCollection = std::vector<L1Track>;

  explicit L1SmartPixelsRefitTableProducer(const edm::ParameterSet& cfg)
      : tracksToken_(consumes<TTTrackCollection>(cfg.getParameter<edm::InputTag>("tracks"))),
        sidecarToken_(consumes<smartpixels::SmartPixelsRefitSidecar>(cfg.getParameter<edm::InputTag>("sidecar"))),
        trackTableName_(cfg.getParameter<std::string>("trackTableName")),
        hitTableName_(cfg.getParameter<std::string>("hitTableName")) {
    produces<nanoaod::FlatTable>("hit");
    produces<nanoaod::FlatTable>("trk");
  }

  void produce(edm::Event& iEvent, const edm::EventSetup&) override {
    const auto& tracks = iEvent.get(tracksToken_);
    const auto& sidecar = iEvent.get(sidecarToken_);

    const size_t nTracks = tracks.size();
    if (sidecar.trackInfo.size() != nTracks || sidecar.hitInfo.size() != nTracks) {
      throw cms::Exception("SmartPixelsSyncBroken")
          << "SmartPixels refit sidecar not 1:1 with track collection: tracks=" << nTracks
          << " trackInfo=" << sidecar.trackInfo.size() << " hitInfo=" << sidecar.hitInfo.size()
          << " (spec RefitSidecarSpec.md §1 output-sync invariant).";
    }

    // ---- EXTENSION table on the variant track table (extension=true, same name/length) ----
    std::vector<uint8_t> status, nCrossings, nAcceptedHits, nKFUpdates, layerHitMask, maxWindowMult;
    std::vector<bool> refitPerformed, seedCovOK, parametrizedSeed, anyWindowTruncated;
    std::vector<float> chi2IncXTot, chi2IncYTot, chi2IncAlphaTot, chi2IncBetaTot;
    std::vector<int32_t> compactWord;
    status.reserve(nTracks);
    nCrossings.reserve(nTracks);
    nAcceptedHits.reserve(nTracks);
    nKFUpdates.reserve(nTracks);
    layerHitMask.reserve(nTracks);
    maxWindowMult.reserve(nTracks);
    refitPerformed.reserve(nTracks);
    seedCovOK.reserve(nTracks);
    parametrizedSeed.reserve(nTracks);
    anyWindowTruncated.reserve(nTracks);
    chi2IncXTot.reserve(nTracks);
    chi2IncYTot.reserve(nTracks);
    chi2IncAlphaTot.reserve(nTracks);
    chi2IncBetaTot.reserve(nTracks);
    compactWord.reserve(nTracks);

    // ---- per-hit LINK table (one row per crossing across all tracks) ----
    std::vector<int32_t> hitTrackIdx;
    std::vector<uint8_t> hitLayer;
    std::vector<uint32_t> hitDetId;
    std::vector<uint16_t> hitWindowMult;
    std::vector<int32_t> hitFlags;
    std::vector<bool> hitAccepted, windowTruncated, hasAlpha, hasBeta;
    std::vector<float> projResX, projResY, recoCotAlpha, recoCotBeta, sigAlpha, sigBeta;
    std::vector<float> recoLocalX, recoLocalY, sigX, sigY, recoCharge;
    std::vector<uint8_t> recoSizeX, recoSizeY;
    std::vector<bool> clusterMerged;
    std::vector<float> truthChargeFrac;
    std::vector<float> pullX, pullY, pullAlpha, pullBeta, selChi2Margin;
    std::vector<float> chi2IncX, chi2IncY, chi2IncAlpha, chi2IncBeta;
    std::vector<int32_t> selHitClass;
    std::vector<float> truthCotAlpha, truthCotBeta;

    for (size_t it = 0; it < nTracks; ++it) {
      const auto& ti = sidecar.trackInfo[it];
      status.push_back(ti.status);
      nCrossings.push_back(ti.nCrossings);
      nAcceptedHits.push_back(ti.nAcceptedHits);
      nKFUpdates.push_back(ti.nKFUpdates);
      layerHitMask.push_back(ti.layerHitMask);
      maxWindowMult.push_back(ti.maxWindowMult);
      refitPerformed.push_back(ti.status & smartpixels::trackstatus::kRefitPerformed);
      seedCovOK.push_back(ti.status & smartpixels::trackstatus::kSeedCovOK);
      parametrizedSeed.push_back(ti.status & smartpixels::trackstatus::kParametrizedSeed);
      anyWindowTruncated.push_back(ti.status & smartpixels::trackstatus::kAnyWindowTruncated);
      chi2IncXTot.push_back(ti.chi2IncXTot);
      chi2IncYTot.push_back(ti.chi2IncYTot);
      chi2IncAlphaTot.push_back(ti.chi2IncAlphaTot);
      chi2IncBetaTot.push_back(ti.chi2IncBetaTot);
      compactWord.push_back(static_cast<int32_t>(smartpixels::packCompactWord(ti)));

      for (const auto& hi : sidecar.hitInfo[it]) {
        hitTrackIdx.push_back(static_cast<int32_t>(it));
        hitLayer.push_back(hi.layer);
        hitDetId.push_back(hi.detId);
        hitWindowMult.push_back(hi.windowMult);
        hitFlags.push_back(hi.flags);
        hitAccepted.push_back(hi.flags & smartpixels::hitflag::kHitAccepted);
        windowTruncated.push_back(hi.flags & smartpixels::hitflag::kWindowTruncated);
        hasAlpha.push_back(hi.flags & smartpixels::hitflag::kHasAlpha);
        hasBeta.push_back(hi.flags & smartpixels::hitflag::kHasBeta);
        recoLocalX.push_back(hi.recoLocalX);
        recoLocalY.push_back(hi.recoLocalY);
        sigX.push_back(hi.sigX);
        sigY.push_back(hi.sigY);
        recoSizeX.push_back(hi.recoSizeX);
        recoSizeY.push_back(hi.recoSizeY);
        recoCharge.push_back(hi.recoCharge);
        clusterMerged.push_back(hi.flags & smartpixels::hitflag::kClusterMerged);
        truthChargeFrac.push_back(hi.truthChargeFrac);
        projResX.push_back(hi.projResX);
        projResY.push_back(hi.projResY);
        recoCotAlpha.push_back(hi.recoCotAlpha);
        recoCotBeta.push_back(hi.recoCotBeta);
        sigAlpha.push_back(hi.sigAlpha);
        sigBeta.push_back(hi.sigBeta);
        pullX.push_back(hi.pullX);
        pullY.push_back(hi.pullY);
        pullAlpha.push_back(hi.pullAlpha);
        pullBeta.push_back(hi.pullBeta);
        chi2IncX.push_back(hi.chi2IncX);
        chi2IncY.push_back(hi.chi2IncY);
        chi2IncAlpha.push_back(hi.chi2IncAlpha);
        chi2IncBeta.push_back(hi.chi2IncBeta);
        selChi2Margin.push_back(hi.selChi2Margin);
        selHitClass.push_back(hi.selHitClass);
        truthCotAlpha.push_back(hi.truthCotAlpha);
        truthCotBeta.push_back(hi.truthCotBeta);
      }
    }

    const size_t nHitRows = hitTrackIdx.size();

    // per-hit LINK table (standalone, not an extension)
    auto hitTable = std::make_unique<nanoaod::FlatTable>(nHitRows, hitTableName_, false, false);
    hitTable->addColumn<int32_t>(
        "trackIdx", hitTrackIdx, "index of the owning track in the " + trackTableName_ + " table");
    hitTable->addColumn<uint8_t>("layer", hitLayer, "TBPX layer 1..4");
    hitTable->addColumn<uint32_t>("detId", hitDetId, "module rawId of the crossing");
    hitTable->addColumn<uint16_t>("windowMult", hitWindowMult, "digis collected in window (post truncation)");
    hitTable->addColumn<int32_t>("flags", hitFlags, "packed flags (bit0 hitAccepted, bit1 windowTruncated, bit2 hasAlpha, bit3 hasBeta)");
    hitTable->addColumn<bool>("hitAccepted", hitAccepted, "selected hit accepted into the KF");
    hitTable->addColumn<bool>("windowTruncated", windowTruncated, "window hit the maxHitsPerWindow truncation");
    hitTable->addColumn<bool>("hasAlpha", hasAlpha, "synthesized cotAlpha available");
    hitTable->addColumn<bool>("hasBeta", hasBeta, "synthesized cotBeta available");
    hitTable->addColumn<float>("recoLocalX", recoLocalX, "reco cluster local x [cm] (-999 if none)");
    hitTable->addColumn<float>("recoLocalY", recoLocalY, "reco cluster local y [cm] (-999 if none)");
    hitTable->addColumn<float>("sigX", sigX, "reco local-x uncertainty from the pixel CPE [cm] (-999 if none)");
    hitTable->addColumn<float>("sigY", sigY, "reco local-y uncertainty from the pixel CPE [cm] (-999 if none)");
    hitTable->addColumn<uint8_t>("recoSizeX", recoSizeX, "selected cluster extent in pixels along local x");
    hitTable->addColumn<uint8_t>("recoSizeY", recoSizeY, "selected cluster extent in pixels along local y");
    hitTable->addColumn<float>("recoCharge", recoCharge, "selected cluster charge [ADC] (-999 if none)", /*mantissaBits=*/12);
    hitTable->addColumn<bool>("clusterMerged", clusterMerged, "TRUTH-ONLY: a second TP contributes more than clusterMergeFrac of the cluster charge");
    hitTable->addColumn<float>("truthChargeFrac", truthChargeFrac, "TRUTH-ONLY dominant contributor share of the cluster charge (-999 if none)", /*mantissaBits=*/12);
    hitTable->addColumn<float>("projResX", projResX, "reco minus projected-crossing local x [cm] (-999 if none)");
    hitTable->addColumn<float>("projResY", projResY, "reco minus projected-crossing local y [cm] (-999 if none)");
    hitTable->addColumn<float>("recoCotAlpha", recoCotAlpha, "reco cotAlpha of the selected hit (-999 if none)", /*mantissaBits=*/12);
    hitTable->addColumn<float>("recoCotBeta", recoCotBeta, "reco cotBeta of the selected hit (-999 if none)", /*mantissaBits=*/12);
    hitTable->addColumn<float>("sigAlpha", sigAlpha, "per-hit cotAlpha sigma from the PixelAV payload (-999 if none)", /*mantissaBits=*/12);
    hitTable->addColumn<float>("sigBeta", sigBeta, "per-hit cotBeta sigma from the PixelAV payload (-999 if none)", /*mantissaBits=*/12);
    hitTable->addColumn<float>("pullX", pullX, "KF pull x = r/sqrt(S) (-999 if none)");
    hitTable->addColumn<float>("pullY", pullY, "KF pull y = r/sqrt(S) (-999 if none)");
    hitTable->addColumn<float>("pullAlpha", pullAlpha, "KF pull cotAlpha = r/sqrt(S) (-999 if none)");
    hitTable->addColumn<float>("pullBeta", pullBeta, "KF pull cotBeta = r/sqrt(S) (-999 if none)");
    hitTable->addColumn<float>("chi2IncX", chi2IncX, "crossing chi2 increment, local-x position term (-999 if none; 0 if not applied)", /*mantissaBits=*/12);
    hitTable->addColumn<float>("chi2IncY", chi2IncY, "crossing chi2 increment, local-y position term (-999 if none; 0 if not applied)", /*mantissaBits=*/12);
    hitTable->addColumn<float>("chi2IncAlpha", chi2IncAlpha, "crossing chi2 increment, cotAlpha angle term (-999 if none; 0 if not applied)", /*mantissaBits=*/12);
    hitTable->addColumn<float>("chi2IncBeta", chi2IncBeta, "crossing chi2 increment, cotBeta angle term (-999 if none; 0 if not applied)", /*mantissaBits=*/12);
    hitTable->addColumn<float>("selChi2Margin", selChi2Margin, "runner-up minus best selection chi2 (>=0; -999 if <2 candidates or no accepted hit)", /*mantissaBits=*/12);
    hitTable->addColumn<int32_t>("selHitClass", selHitClass, "TRUTH-ONLY simlink class of selected hit: 0 sameTP, 1 otherTP, 2 noise, -1 none");
    hitTable->addColumn<float>("truthCotAlpha", truthCotAlpha, "TRUTH-ONLY unsmeared parent local cotAlpha of the selected hit (-999 if none)", /*mantissaBits=*/12);
    hitTable->addColumn<float>("truthCotBeta", truthCotBeta, "TRUTH-ONLY unsmeared parent local cotBeta of the selected hit (-999 if none)", /*mantissaBits=*/12);
    hitTable->setDoc("SmartPixels refit per-crossing records for the " + trackTableName_ +
                     " tracks (one row per layer crossing; trackIdx links to that track table)");
    iEvent.put(std::move(hitTable), "hit");

    // EXTENSION table on the variant track table (extension=true -> merged by name)
    auto trkTable = std::make_unique<nanoaod::FlatTable>(nTracks, trackTableName_, false, true);
    trkTable->addColumn<uint8_t>("spxStatus", status, "packed refit status (bit0 refitPerformed, bit1 seedCovOK, bit2 parametrizedSeed, bit3 anyWindowTruncated)");
    trkTable->addColumn<bool>("spxRefitPerformed", refitPerformed, "refit performed (else passthrough copy of input track)");
    trkTable->addColumn<bool>("spxSeedCovOK", seedCovOK, "seed covariance was usable");
    trkTable->addColumn<bool>("spxParametrizedSeed", parametrizedSeed, "seed covariance came from the parametrized model");
    trkTable->addColumn<bool>("spxAnyWindowTruncated", anyWindowTruncated, "at least one crossing window hit maxHitsPerWindow");
    trkTable->addColumn<uint8_t>("spxNCrossings", nCrossings, "valid layer crossings attempted");
    trkTable->addColumn<uint8_t>("spxNAcceptedHits", nAcceptedHits, "hits accepted into the KF");
    trkTable->addColumn<uint8_t>("spxNKFUpdates", nKFUpdates, "scalar-update groups applied (layers updated)");
    trkTable->addColumn<uint8_t>("spxLayerHitMask", layerHitMask, "accepted-hit bitmask bit0=L1..bit3=L4; popcount == nAcceptedHits");
    trkTable->addColumn<uint16_t>("spxMaxWindowMult", maxWindowMult, "max window multiplicity over this track's crossings");
    trkTable->addColumn<float>("spxChi2IncXTot", chi2IncXTot, "sum of local-x chi2 increments over crossings");
    trkTable->addColumn<float>("spxChi2IncYTot", chi2IncYTot, "sum of local-y chi2 increments over crossings");
    trkTable->addColumn<float>("spxChi2IncAlphaTot", chi2IncAlphaTot, "sum of cotAlpha chi2 increments over crossings (r-phi total = X + Alpha)");
    trkTable->addColumn<float>("spxChi2IncBetaTot", chi2IncBetaTot, "sum of cotBeta chi2 increments over crossings (r-z total = Y + Beta)");
    trkTable->addColumn<int32_t>("spxCompactWord", compactWord, "16-bit transmitted-subset compact word (packCompactWord; see SmartPixelsTransmittedSubset.h)");
    iEvent.put(std::move(trkTable), "trk");
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("tracks", edm::InputTag("l1tSmartPixelsTrackProducer", "Level1TTTracks"));
    desc.add<edm::InputTag>("sidecar", edm::InputTag("l1tSmartPixelsTrackProducer", "Level1TTTracks"));
    desc.add<std::string>("trackTableName", "L1TSmartPixelsTrack");
    desc.add<std::string>("hitTableName", "L1TSmartPixelsRefitHit");
    descriptions.addWithDefaultLabel(desc);
  }

private:
  const edm::EDGetTokenT<TTTrackCollection> tracksToken_;
  const edm::EDGetTokenT<smartpixels::SmartPixelsRefitSidecar> sidecarToken_;
  const std::string trackTableName_;
  const std::string hitTableName_;
};

DEFINE_FWK_MODULE(L1SmartPixelsRefitTableProducer);
