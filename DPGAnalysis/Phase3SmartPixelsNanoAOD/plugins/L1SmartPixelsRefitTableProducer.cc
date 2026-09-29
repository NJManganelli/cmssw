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
// Links are written only when their target table is in the nano: an empty
// trackTableName drops the "trk" extension and the hit trackIdx column, an empty
// clusterTableName drops selClusterIdx.
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

#include "DPGAnalysis/Phase2L1TNanoAOD/interface/L1TPRowLink.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"

#include "L1Trigger/Phase3SmartPixels/interface/SmartPixelsRefitSidecar.h"
#include "L1Trigger/Phase3SmartPixels/interface/SmartPixelsTransmittedSubset.h"

#include <cstdint>
#include <string>
#include <vector>

class L1SmartPixelsRefitTableProducer : public edm::stream::EDProducer<> {
public:
  using L1Track = TTTrack<Ref_Phase2TrackerDigi_>;
  using TTTrackCollection = std::vector<L1Track>;

  explicit L1SmartPixelsRefitTableProducer(const edm::ParameterSet& cfg)
      : tracksToken_(consumes<TTTrackCollection>(cfg.getParameter<edm::InputTag>("tracks"))),
        sidecarToken_(consumes<smartpixels::SmartPixelsRefitSidecar>(cfg.getParameter<edm::InputTag>("sidecar"))),
        tpRowMapToken_(consumes<edm::ValueMap<int>>(cfg.getParameter<edm::InputTag>("tpRowMap"))),
        trackTableName_(cfg.getParameter<std::string>("trackTableName")),
        clusterTableName_(cfg.getParameter<std::string>("clusterTableName")),
        hitTableName_(cfg.getParameter<std::string>("hitTableName")) {
    produces<nanoaod::FlatTable>("hit");
    if (!trackTableName_.empty())
      produces<nanoaod::FlatTable>("trk");
  }

  void produce(edm::Event& iEvent, const edm::EventSetup&) override {
    const auto& tracks = iEvent.get(tracksToken_);
    const auto& sidecar = iEvent.get(sidecarToken_);
    // The sidecar holds a bare key into the refit's TrackingParticle input, the same
    // mix:MergedTrackTruth the row map covers.
    const auto tpRowMapHandle = iEvent.getHandle(tpRowMapToken_);
    const edm::ValueMap<int>* tpRowMap = tpRowMapHandle.isValid() ? tpRowMapHandle.product() : nullptr;
    const edm::ProductID tpId = tpRowMap ? l1tnano::tpProductOf(*tpRowMap) : edm::ProductID();

    const size_t nTracks = tracks.size();
    if (sidecar.trackInfo.size() != nTracks || sidecar.hitInfo.size() != nTracks) {
      throw cms::Exception("SmartPixelsSyncBroken")
          << "SmartPixels refit sidecar not 1:1 with track collection: tracks=" << nTracks
          << " trackInfo=" << sidecar.trackInfo.size() << " hitInfo=" << sidecar.hitInfo.size()
          << " (spec RefitSidecarSpec.md §1 output-sync invariant).";
    }

    // ---- EXTENSION table on the variant track table (extension=true, same name/length) ----
    std::vector<uint8_t> status, nCrossings, nAcceptedHits, nKFUpdates, layerHitMask, maxWindowMult, nClustersOwnedByMatchedTp,
        nITClustersOwnedByMatchedTp;
    std::vector<bool> refitPerformed, otCovOK, anyWindowTruncated;
    std::vector<float> chi2IncXTot, chi2IncYTot, chi2IncAlphaTot, chi2IncBetaTot;
    std::vector<float> chi2ITAtOT, chi2ITAtRefit, shiftChi2, logDetRatio;
    std::vector<int32_t> matchedTpIdx;
    std::vector<int32_t> compactWord;
    status.reserve(nTracks);
    nCrossings.reserve(nTracks);
    nAcceptedHits.reserve(nTracks);
    nKFUpdates.reserve(nTracks);
    layerHitMask.reserve(nTracks);
    nClustersOwnedByMatchedTp.reserve(nTracks);
    nITClustersOwnedByMatchedTp.reserve(nTracks);
    maxWindowMult.reserve(nTracks);
    refitPerformed.reserve(nTracks);
    otCovOK.reserve(nTracks);
    anyWindowTruncated.reserve(nTracks);
    chi2IncXTot.reserve(nTracks);
    chi2IncYTot.reserve(nTracks);
    chi2IncAlphaTot.reserve(nTracks);
    chi2IncBetaTot.reserve(nTracks);
    chi2ITAtOT.reserve(nTracks);
    chi2ITAtRefit.reserve(nTracks);
    shiftChi2.reserve(nTracks);
    logDetRatio.reserve(nTracks);
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
    std::vector<int32_t> selClusterIdx;
    std::vector<float> projLocalX, projLocalY, projCotAlpha, projCotBeta,
        projSigX, projSigY,
        projOTLocalX, projOTLocalY, projOTSigX, projOTSigY,
        projOTCotAlpha, projOTCotBeta;
    std::vector<uint8_t> recoSizeX, recoSizeY;
    std::vector<bool> clusterMerged;
    std::vector<float> tpChargeFrac;
    std::vector<float> pullX, pullY, pullAlpha, pullBeta, selChi2Margin;
    std::vector<float> chi2IncX, chi2IncY, chi2IncAlpha, chi2IncBeta;
    std::vector<int32_t> selHitClass;
    std::vector<float> tpLocalCotAlpha, tpLocalCotBeta;

    for (size_t it = 0; it < nTracks; ++it) {
      const auto& ti = sidecar.trackInfo[it];
      status.push_back(ti.status);
      nCrossings.push_back(ti.nCrossings);
      nAcceptedHits.push_back(ti.nAcceptedHits);
      nKFUpdates.push_back(ti.nKFUpdates);
      layerHitMask.push_back(ti.layerHitMask);
      maxWindowMult.push_back(ti.maxWindowMult);
      refitPerformed.push_back(ti.status & smartpixels::trackstatus::kRefitPerformed);
      otCovOK.push_back(ti.status & smartpixels::trackstatus::kOTCovOK);
      anyWindowTruncated.push_back(ti.status & smartpixels::trackstatus::kAnyWindowTruncated);
      chi2IncXTot.push_back(ti.chi2IncXTot);
      chi2IncYTot.push_back(ti.chi2IncYTot);
      chi2IncAlphaTot.push_back(ti.chi2IncAlphaTot);
      chi2IncBetaTot.push_back(ti.chi2IncBetaTot);
      matchedTpIdx.push_back(l1tnano::tpIdxOfKey(tpRowMap, tpId, ti.matchedTpIdx));
      chi2ITAtOT.push_back(ti.chi2ITAtOT);
      chi2ITAtRefit.push_back(ti.chi2ITAtRefit);
      shiftChi2.push_back(ti.shiftChi2);
      logDetRatio.push_back(ti.logDetRatio);
      compactWord.push_back(static_cast<int32_t>(smartpixels::packCompactWord(ti)));
      nClustersOwnedByMatchedTp.push_back(ti.nClustersOwnedByMatchedTp);
      nITClustersOwnedByMatchedTp.push_back(ti.nITClustersOwnedByMatchedTp);

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
        selClusterIdx.push_back(hi.selClusterIdx);
        projLocalX.push_back(hi.projLocalX);
        projLocalY.push_back(hi.projLocalY);
        projCotAlpha.push_back(hi.projCotAlpha);
        projCotBeta.push_back(hi.projCotBeta);
        projOTLocalX.push_back(hi.projOTLocalX);
        projOTLocalY.push_back(hi.projOTLocalY);
        projSigX.push_back(hi.projSigX);
        projSigY.push_back(hi.projSigY);
        projOTSigX.push_back(hi.projOTSigX);
        projOTSigY.push_back(hi.projOTSigY);
        projOTCotAlpha.push_back(hi.projOTCotAlpha);
        projOTCotBeta.push_back(hi.projOTCotBeta);
        recoLocalX.push_back(hi.recoLocalX);
        recoLocalY.push_back(hi.recoLocalY);
        sigX.push_back(hi.sigX);
        sigY.push_back(hi.sigY);
        recoSizeX.push_back(hi.recoSizeX);
        recoSizeY.push_back(hi.recoSizeY);
        recoCharge.push_back(hi.recoCharge);
        clusterMerged.push_back(hi.flags & smartpixels::hitflag::kClusterMerged);
        tpChargeFrac.push_back(hi.tpChargeFrac);
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
        tpLocalCotAlpha.push_back(hi.tpLocalCotAlpha);
        tpLocalCotBeta.push_back(hi.tpLocalCotBeta);
      }
    }

    const size_t nHitRows = hitTrackIdx.size();

    // per-hit LINK table (standalone, not an extension)
    auto hitTable = std::make_unique<nanoaod::FlatTable>(nHitRows, hitTableName_, false, false);
    if (!trackTableName_.empty())
      hitTable->addColumn<int32_t>(
          "trackIdx", hitTrackIdx, "index of the owning track in the " + trackTableName_ + " table");
    hitTable->addColumn<uint8_t>("layer", hitLayer,
        "IT surface of the crossing: 1-4 = TBPX IL1-IL4, 11-14 = ID1-ID4 (L1TSmartPixelsCluster_surface)");
    hitTable->addColumn<uint32_t>("detId", hitDetId, "module rawId of the crossing");
    hitTable->addColumn<uint16_t>("windowMult", hitWindowMult, "digis collected in window (post truncation)");
    hitTable->addColumn<int32_t>("flags", hitFlags, "packed flags (bit0 hitAccepted, bit1 windowTruncated, bit2 hasAlpha, bit3 hasBeta)");
    hitTable->addColumn<bool>("hitAccepted", hitAccepted, "selected hit accepted into the KF");
    hitTable->addColumn<bool>("windowTruncated", windowTruncated, "window hit the maxHitsPerWindow truncation");
    hitTable->addColumn<bool>("hasAlpha", hasAlpha, "synthesized cotAlpha available");
    hitTable->addColumn<bool>("hasBeta", hasBeta, "synthesized cotBeta available");
    if (!clusterTableName_.empty())
      hitTable->addColumn<int32_t>(
          "selClusterIdx", selClusterIdx,
          "row of the SELECTED cluster in " + clusterTableName_ +
              " (-1 if none). EXACT link: both "
              "tables walk the same SiPixelRecHitCollection with the same filter and order. Consumers "
              "MUST assert cluster[selClusterIdx].detId == detId -- matching on position instead fails, "
              "since both store coordinates at ~5 um nano precision against a 25 um pitch");
    hitTable->addColumn<float>("projLocalX", projLocalX,
        "RUNNING projection: track state as of this layer's visit, local x [cm]");
    hitTable->addColumn<float>("projLocalY", projLocalY, "RUNNING projection, local y [cm]");
    hitTable->addColumn<float>("projCotAlpha", projCotAlpha, "RUNNING projection, expected cotAlpha", 12);
    hitTable->addColumn<float>("projCotBeta", projCotBeta, "RUNNING projection, expected cotBeta", 12);
    hitTable->addColumn<float>("projOTLocalX", projOTLocalX,
        "OT-ONLY projection (the OT track, no Kalman updates), local x [cm]");
    hitTable->addColumn<float>("projOTLocalY", projOTLocalY, "OT-ONLY projection, local y [cm]");
    hitTable->addColumn<float>("projSigX", projSigX,
        "sqrt((H C H^T)_xx) [cm]: REFIT-ORDER projection cone from the RUNNING covariance, "
        "recorded after multiple-scattering Q and before this layer's update, EXCLUDING the "
        "measurement term. Compare against projOTSigX to isolate what the refit buys. "
        "Present on every crossing, including the ~10% whose window is empty");
    hitTable->addColumn<float>("projSigY", projSigY, "sqrt((H C H^T)_yy) [cm]");
    hitTable->addColumn<float>("projOTSigX", projOTSigX,
        "sqrt((H C_OT H^T)_xx) [cm]: single-shot projection cone from the OT covariance, "
        "EXCLUDING the measurement term");
    hitTable->addColumn<float>("projOTSigY", projOTSigY, "sqrt((H C_OT H^T)_yy) [cm]");
    hitTable->addColumn<float>("projOTCotAlpha", projOTCotAlpha, "OT-ONLY expected cotAlpha", 12);
    hitTable->addColumn<float>("projOTCotBeta", projOTCotBeta, "OT-ONLY expected cotBeta", 12);
    hitTable->addColumn<float>("recoLocalX", recoLocalX, "reco cluster local x [cm] (-999 if none)");
    hitTable->addColumn<float>("recoLocalY", recoLocalY, "reco cluster local y [cm] (-999 if none)");
    hitTable->addColumn<float>("sigX", sigX, "reco local-x uncertainty from the pixel CPE [cm] (-999 if none)");
    hitTable->addColumn<float>("sigY", sigY, "reco local-y uncertainty from the pixel CPE [cm] (-999 if none)");
    hitTable->addColumn<uint8_t>("recoSizeX", recoSizeX, "selected cluster extent in pixels along local x");
    hitTable->addColumn<uint8_t>("recoSizeY", recoSizeY, "selected cluster extent in pixels along local y");
    hitTable->addColumn<float>("recoCharge", recoCharge, "selected cluster charge [ADC] (-999 if none)", /*mantissaBits=*/12);
    hitTable->addColumn<bool>("clusterMerged", clusterMerged, "TRUTH-ONLY: a second TP contributes more than clusterMergeFrac of the cluster charge");
    hitTable->addColumn<float>("tpChargeFrac", tpChargeFrac, "TRUTH-ONLY dominant contributor share of the cluster charge (-999 if none)", /*mantissaBits=*/12);
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
    hitTable->addColumn<int32_t>("selHitClass", selHitClass, "TRUTH-ONLY simlink class of the selected hit by its dominant charge contributor: 0 the OT track's "
        "loosely-genuine TP (L1TTrack_tpIdx, not the SmartPixels track's tpIdx), 1 another TP, 2 noise, -1 no hit accepted");
    hitTable->addColumn<float>("tpLocalCotAlpha", tpLocalCotAlpha, "TRUTH-ONLY unsmeared parent local cotAlpha of the selected hit (-999 if none)", /*mantissaBits=*/12);
    hitTable->addColumn<float>("tpLocalCotBeta", tpLocalCotBeta, "TRUTH-ONLY unsmeared parent local cotBeta of the selected hit (-999 if none)", /*mantissaBits=*/12);
    hitTable->setDoc(trackTableName_.empty()
                         ? std::string("SmartPixels refit per-crossing records (one row per layer crossing)")
                         : "SmartPixels refit per-crossing records for the " + trackTableName_ +
                               " tracks (one row per layer crossing; trackIdx links to that track table)");
    iEvent.put(std::move(hitTable), "hit");
    if (trackTableName_.empty())
      return;

    // EXTENSION table on the variant track table (extension=true -> merged by name)
    auto trkTable = std::make_unique<nanoaod::FlatTable>(nTracks, trackTableName_, false, true);
    trkTable->addColumn<uint8_t>("spixStatus", status, "packed refit status (bit0 refitPerformed, bit1 otCovOK, bit2 RETIRED, bit3 anyWindowTruncated)");
    trkTable->addColumn<bool>("spixRefitPerformed", refitPerformed, "refit performed (else passthrough copy of input track)");
    trkTable->addColumn<bool>("spixOTCovOK", otCovOK, "the OT track covariance was usable");
    trkTable->addColumn<bool>("spixAnyWindowTruncated", anyWindowTruncated, "at least one crossing window hit maxHitsPerWindow");
    trkTable->addColumn<uint8_t>("spixNCrossings", nCrossings, "valid layer crossings attempted");
    trkTable->addColumn<uint8_t>("spixNAcceptedHits", nAcceptedHits, "hits accepted into the KF");
    trkTable->addColumn<uint8_t>("spixNKFUpdates", nKFUpdates, "scalar-update groups applied (layers updated)");
    trkTable->addColumn<uint8_t>("spixLayerHitMask", layerHitMask, "accepted-hit bitmask, bits 0-3 = IL1-IL4, bits 4-7 = ID1-ID4; popcount == nAcceptedHits");
    trkTable->addColumn<uint16_t>("spixMaxWindowMult", maxWindowMult, "max window multiplicity over this track's crossings");
    trkTable->addColumn<float>("spixChi2IncXTot", chi2IncXTot, "sum of local-x chi2 increments over crossings");
    trkTable->addColumn<float>("spixChi2IncYTot", chi2IncYTot, "sum of local-y chi2 increments over crossings");
    trkTable->addColumn<float>("spixChi2IncAlphaTot", chi2IncAlphaTot, "sum of cotAlpha chi2 increments over crossings (r-phi total = X + Alpha)");
    trkTable->addColumn<float>("spixChi2IncBetaTot", chi2IncBetaTot, "sum of cotBeta chi2 increments over crossings (r-z total = Y + Beta)");
    trkTable->addColumn<int32_t>("tpIdx", matchedTpIdx,
        std::string("TRUTH-ONLY, SmartPixels cluster-majority match: the TP owning a strict majority "
                    "(> 50%) of this track's clusters, i.e. the 2 clusters of each OT stub (L1TOTStub_tpIdxInner, "
                    "L1TOTStub_tpIdxOuter) plus, when refitted, its accepted IT clusters "
                    "(L1TSmartPixelsCluster_tpIdx); passthrough tracks: the stub clusters only. A cluster's owner "
                    "is its dominant contributor (OT: TTClusterAssociationMap; IT: charge); ownerless clusters "
                    "count in the denominator. NOT the loosely-genuine TTTrackAssociator match of L1TTrack_tpIdx "
                    "(same name, different definition; compare via l1TrackIdx): ") +
            l1tnano::tpIdxDoc());
    trkTable->addColumn<uint8_t>("tpNClusters", nClustersOwnedByMatchedTp,
        "TRUTH-ONLY, clusters (OT stub clusters + accepted IT clusters) owned by the tpIdx TP; 0 if unmatched. "
        "Cluster purity = this / (2 x L1TTrack_nStubs[l1TrackIdx] + spixNAcceptedHits)");
    trkTable->addColumn<uint8_t>("tpNITClusters", nITClustersOwnedByMatchedTp,
        "TRUTH-ONLY, accepted IT clusters owned by the tpIdx TP; 0 if unmatched. IT purity = this / "
        "spixNAcceptedHits; OT stub clusters owned = tpNClusters - tpNITClusters");
    trkTable->addColumn<float>("spixChi2ITAtOT", chi2ITAtOT,
        "chi2 of the ACCEPTED IT hits against the UNMODIFIED OT track helix. Needs no stubs, so "
        "it travels with the track. CAVEAT: the refit selected these hits by minimising this, so "
        "a large OT-to-refit drop is equally the signature of chasing wrong hits");
    trkTable->addColumn<float>("spixChi2ITAtRefit", chi2ITAtRefit,
        "same accepted IT hits against the FINAL refit helix");
    trkTable->addColumn<float>("spixShiftChi2", shiftChi2,
        "da^T (C_OT - C_refit)^-1 da: how far the fit moved, in units of the information that "
        "moving it required (-999 if the covariance difference is not invertible)");
    trkTable->addColumn<float>("spixLogDetRatio", logDetRatio,
        "ln(det C_OT / det C_refit) >= 0: total information gained, parametrisation-independent");
    trkTable->addColumn<int32_t>("spixCompactWord", compactWord, "16-bit transmitted-subset compact word (packCompactWord; see SmartPixelsTransmittedSubset.h)");
    iEvent.put(std::move(trkTable), "trk");
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("tracks", edm::InputTag("l1tSmartPixelsTrackProducer", "Level1TTTracks"));
    desc.add<edm::InputTag>("sidecar", edm::InputTag("l1tSmartPixelsTrackProducer", "Level1TTTracks"));
    desc.add<edm::InputTag>("tpRowMap", edm::InputTag("l1tTrackingParticleTable"))
        ->setComment("L1TTP key -> row map from L1TrackingParticleTableProducer; tpIdx = -2 when absent");
    desc.add<std::string>("trackTableName", "L1TSmartPixelsTrack")
        ->setComment("track table the spix*/tp* extension and the hit trackIdx belong to; empty = neither is written");
    desc.add<std::string>("clusterTableName", "L1TSmartPixelsCluster")
        ->setComment("cluster table selClusterIdx indexes; empty = no selClusterIdx column");
    desc.add<std::string>("hitTableName", "L1TSmartPixelsRefitHit");
    descriptions.addWithDefaultLabel(desc);
  }

private:
  const edm::EDGetTokenT<TTTrackCollection> tracksToken_;
  const edm::EDGetTokenT<smartpixels::SmartPixelsRefitSidecar> sidecarToken_;
  const edm::EDGetTokenT<edm::ValueMap<int>> tpRowMapToken_;
  const std::string trackTableName_;
  const std::string clusterTableName_;
  const std::string hitTableName_;
};

DEFINE_FWK_MODULE(L1SmartPixelsRefitTableProducer);
