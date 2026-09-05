#ifndef L1Trigger_Phase3SmartPixels_SmartPixelsRefitSidecar_h
#define L1Trigger_Phase3SmartPixels_SmartPixelsRefitSidecar_h

// -*- C++ -*-
//
// Persistent SmartPixels refit sidecar data model. Authoritative contract:
// L1Trigger/Phase3SmartPixels/doc/RefitSidecarSpec.md §2. This header
// implements §2 of that spec EXACTLY (field names, types, sentinel rules). Any
// change requires a version bump in the spec + classes_def.xml first.
//
// One SmartPixelsRefitSidecar is emitted per refit track collection (same
// module/event as the tracks). It is the single persistent home for every NEW
// SmartPixels fact that does not exist in the OT-only TTTrack (per-layer hit
// residuals, synthesized angles, KF pulls, chi2 increments, window occupancy).
// The 1:1 output-sync invariant (spec §1) means row i of every vector derives
// from input track i, so downstream views are cheap index-aligned adapters.
//
// Sentinel rules (spec §2): any unavailable float is -999.f (consumers MUST test
// > -900.f); integers default to 0. Structs default-initialize accordingly so a
// freshly constructed record is already "empty / not applicable".

#include <cstdint>
#include <vector>

namespace smartpixels {

  // One entry per LAYER CROSSING attempted (a valid crossing, whether or not a
  // hit was accepted) -- NOT one entry per accepted hit.
  struct SmartPixelsRefitHitInfo {
    // --- crossing identity / window occupancy (always valid) ---
    uint8_t layer = 0;         // TBPX layer 1..4
    uint32_t detId = 0;        // module rawId of the crossing (diagnostics)
    uint16_t windowMult = 0;   // digis collected in window (post readout-order truncation)
    uint8_t flags = 0;         // bit0 hitAccepted; bit1 windowTruncated (maxHitsPerWindow hit);
                               // bit2 hasAlpha; bit3 hasBeta; bits4-7 reserved

    // --- selected hit, valid only when hitAccepted (else sentinel -999.f) ---
    // VOCABULARY (spec §2). Four prefixes, never mixed:
    //   reco*     what the sensor produced (a reconstructed quantity; today the
    //             position comes from the Phase-2 pixel CPE and the angles from a
    //             parametrized throw on truth -- both are reco, not truth)
    //   proj*     the OT-only track projected to this layer crossing
    //   projRes*  reco - proj  (the KF innovation numerator)
    //   truth*    unsmeared generator/parent quantity, TRUTH-ONLY
    // "residual" alone is ambiguous -- it means reco-vs-truth for a hit and
    // fitted-vs-truth for a track -- so it never appears unprefixed.
    // NOTE the payload ntuple (SmartPixelsPayloadAnalyzer) keeps its own historical
    // digi_parCot* branch names; they are a different product, not this contract.
    // Reco payload: what the sensor produced for the selected cluster. This is the
    // ONLY group a real sensor could encode into its per-cluster word, since it
    // knows nothing about any track -- hence it is stored in full, not just as a
    // track-relative residual.
    float recoLocalX = -999.f, recoLocalY = -999.f;    // cluster position, module-local [cm]
    float sigX = -999.f, sigY = -999.f;                // position uncertainty from the pixel CPE [cm]
    uint8_t recoSizeX = 0, recoSizeY = 0;              // cluster extent in pixels (shape: an angle handle)
    float recoCharge = -999.f;                         // cluster charge [ADC]
    float projResX = -999.f, projResY = -999.f;        // reco - proj, module-local [cm]
    float recoCotAlpha = -999.f, recoCotBeta = -999.f; // reco incidence angles
    float sigAlpha = -999.f, sigBeta = -999.f;         // per-hit angle sigmas from the PixelAV payload
    float pullX = -999.f, pullY = -999.f;              // KF pulls r_k/sqrt(S_k) from the scalar updates
    float pullAlpha = -999.f, pullBeta = -999.f;
    // Per-dimension scalar-update chi2 increments r^2/S. 0 when the hit was
    // accepted but that update was not applied (angle absent or numerics-gated);
    // -999.f when no hit was accepted.
    float chi2IncX = -999.f, chi2IncY = -999.f;
    float chi2IncAlpha = -999.f, chi2IncBeta = -999.f;
    // --- PROJECTION, stored for EVERY valid crossing (not only where a hit was
    // accepted). Crossings with no accepted hit are exactly where combinatorics
    // matter most, so making these conditional would hide the interesting cases.
    //   proj*      : the state as it stands when this layer is visited, i.e. the
    //                seed updated by whichever layers were visited BEFORE it under
    //                digiRefitLayerOrder. Order-dependent by construction.
    //   projSeed*  : the UNMODIFIED OT-only seed helix projected to this layer, with
    //                no Kalman updates at all. Order-independent, and the honest
    //                "cold start" a system matching all layers in parallel faces.
    // projSeedSig{X,Y} = sqrt(H C_seed H^T) in the module-local frame: the
    // single-shot projection cone from the OT fit's own covariance. It excludes the
    // measurement term, so it is the TRACK's uncertainty, not an innovation sigma.
    float projLocalX = -999.f, projLocalY = -999.f;
    float projCotAlpha = -999.f, projCotBeta = -999.f;
    float projSeedLocalX = -999.f, projSeedLocalY = -999.f;
    float projSeedSigX = -999.f, projSeedSigY = -999.f;
    float projSeedCotAlpha = -999.f, projSeedCotBeta = -999.f;

    // Row index of the SELECTED cluster in the untruncated cluster nano table
    // (L1TSmartPixelsCluster), or -1. Both are produced from the same
    // SiPixelRecHitCollection with the same filter and iteration order, so the
    // index is exact -- which matters because matching on POSITION does not work:
    // both tables store coordinates at 10-bit nano mantissa precision, ~5 um on a
    // 0.5 cm coordinate, against a 25 um pitch.
    // Consumers MUST assert cluster[selClusterIdx].detId == hitInfo.detId; that is a
    // cheap check that the two orderings have not silently diverged.
    int32_t selClusterIdx = -1;

    float selChi2Margin = -999.f;                      // runner-up minus best selection chi2 (>=0); how unambiguous
                                                       // the hit choice was. Sentinel -999.f when no hit accepted or the
                                                       // window held fewer than 2 candidates. Hardware-plausible.

    // --- TRUTH-ONLY (never hardware-available; excluded from every transmitted subset) ---
    int8_t selHitClass = -1;                           // selected-cluster class by DOMINANT charge contributor:
                                                       // 0 sameTP, 1 otherTP, 2 noise (no simlink), -1 none
    // TRUE incidence angles of the selected cluster's dominant TrackingParticle,
    // in the MODULE frame, unsmeared (-999.f if no parent). The tp prefix matches
    // the convention used everywhere else in this nano (L1TTrack_tpPt, the cluster
    // table's tpIdx/tpPt/...); "truth" was a second word for the same thing.
    // local* because cotAlpha/cotBeta are module-frame BY DEFINITION (PixelAV).
    float tpLocalCotAlpha = -999.f, tpLocalCotBeta = -999.f;
    float tpChargeFrac = -999.f;                    // dominant contributor's share of the cluster charge. < 1 means
                                                       // the cluster is shared; see hitflag::kClusterMerged. A cluster
                                                       // can be class 0 and still carry another TP's charge, which
                                                       // biases its position and makes its angle ill-defined -- an
                                                       // effect the pre-cluster (per-digi) hit model could not express.
  };

  // One entry per track (refit or passthrough).
  struct SmartPixelsRefitTrackInfo {
    uint8_t status = 0;         // bit0 refit performed (else passthrough); bit1 seedCovOK;
                                // bit2 RETIRED (was seedCovMode==parametrized, mode removed);
                                // bit3 anyWindowTruncated; bits4-7 reserved
    uint8_t nCrossings = 0;     // valid layer crossings attempted
    uint8_t nAcceptedHits = 0;  // hits accepted into the KF
    uint8_t nKFUpdates = 0;     // scalar-update groups applied (== layers updated)
    uint8_t layerHitMask = 0;   // accepted-hit bitmask, bit0=L1 .. bit3=L4;
                                // popcount(layerHitMask) == nAcceptedHits (exact)
    uint16_t maxWindowMult = 0; // max windowMult over this track's crossings
    // Per-dimension chi2-increment sums over crossings.
    float chi2IncXTot = -999.f, chi2IncYTot = -999.f;
    float chi2IncAlphaTot = -999.f, chi2IncBetaTot = -999.f;

    // --- seed-vs-refit diagnostics, computable WITHOUT stubs -----------------
    // A full joint chi2 of stubs + IT hits against the refit needs the stubs, which
    // a downstream consumer may not have. These use only the IT hits the refit
    // itself accepted, so they travel with the track.
    //
    // chi2ITAtSeed  : chi2 of THOSE hits against the UNMODIFIED OT seed helix.
    // chi2ITAtRefit : the same hits against the final refit helix.
    // Their DIFFERENCE is how much the refit improved its own IT description.
    // Beware the obvious trap: the refit chose these hits by minimising exactly
    // this quantity, so a large improvement is not by itself evidence of a better
    // track -- it is equally the signature of having picked wrong hits that the
    // fit then chased. It discriminates only in combination with the pulls.
    float chi2ITAtSeed = -999.f, chi2ITAtRefit = -999.f;
    // shiftChi2 = da^T (C_seed - C_refit)^-1 da, da = a_refit - a_seed: how far the
    // fit moved measured in units of the information that moving it required.
    // Sentinel if (C_seed - C_refit) is not invertible.
    float shiftChi2 = -999.f;
    // logDetRatio = ln(det C_seed / det C_refit) >= 0: total information gained,
    // parametrisation-independent and insensitive to WHICH direction shrank.
    float logDetRatio = -999.f;
    // Index of the TrackingParticle this track is truth-matched to, or -1.
    // Join key against the cluster table's truthTpIdx: equality means the cluster
    // came from this track's own particle. TRUTH-ONLY.
    int32_t matchedTpIdx = -1;
  };

  struct SmartPixelsRefitSidecar {
    std::vector<SmartPixelsRefitTrackInfo> trackInfo;               // size == N tracks
    std::vector<std::vector<SmartPixelsRefitHitInfo>> hitInfo;      // outer size == N tracks
  };

  // Bit accessors for SmartPixelsRefitHitInfo::flags (spec §2).
  namespace hitflag {
    inline constexpr uint8_t kHitAccepted = 0x1;     // bit0
    inline constexpr uint8_t kWindowTruncated = 0x2; // bit1
    inline constexpr uint8_t kHasAlpha = 0x4;        // bit2
    inline constexpr uint8_t kHasBeta = 0x8;         // bit3
    inline constexpr uint8_t kClusterMerged = 0x10;  // bit4: a second TP contributes
                                                     // more than clusterMergeFrac of the charge (TRUTH-ONLY)
  }  // namespace hitflag

  // Bit accessors for SmartPixelsRefitTrackInfo::status (spec §2).
  namespace trackstatus {
    inline constexpr uint8_t kRefitPerformed = 0x1;   // bit0
    inline constexpr uint8_t kSeedCovOK = 0x2;        // bit1
    // 0x4 (bit2) is RETIRED: it flagged the removed 'parametrized' seedCovMode.
    // Do not reuse it -- pre-removal files may have it set, and a new meaning
    // would silently mis-read them.
    inline constexpr uint8_t kAnyWindowTruncated = 0x8; // bit3
  }  // namespace trackstatus

}  // namespace smartpixels

#endif
