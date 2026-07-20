#ifndef L1Trigger_Phase3SmartPixels_SmartPixelsRefitSidecar_h
#define L1Trigger_Phase3SmartPixels_SmartPixelsRefitSidecar_h

// -*- C++ -*-
//
// Persistent SmartPixels refit sidecar data model. Authoritative contract:
// L1Trigger/Phase3SmartPixels/doc/RefitSidecarSpec.md (spec v0). This header
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
    float resX = -999.f, resY = -999.f;                // selected-hit local residual vs predicted crossing [cm]
    float cotAlphaMeas = -999.f, cotBetaMeas = -999.f; // synthesized measured angles
    float sigAlpha = -999.f, sigBeta = -999.f;         // per-hit angle sigmas from the PixelAV payload
    float pullX = -999.f, pullY = -999.f;              // KF pulls r_k/sqrt(S_k) from the scalar updates
    float pullAlpha = -999.f, pullBeta = -999.f;
    float chi2IncRPhi = -999.f;                        // sum over this crossing's scalar updates of r^2/S, x + alpha terms
    float chi2IncRZ = -999.f;                          //                                              y + beta terms
    float selChi2Margin = -999.f;                      // v0.3: runner-up minus best selection chi2 (>=0); how unambiguous
                                                       // the hit choice was. Sentinel -999.f when no hit accepted or the
                                                       // window held fewer than 2 candidates. Hardware-plausible.

    // --- TRUTH-ONLY (never hardware-available; excluded from every transmitted subset) ---
    int8_t selHitClass = -1;                           // selected-hit simlink class: 0 sameTP, 1 otherTP, 2 noise, -1 none
    float parCotAlpha = -999.f, parCotBeta = -999.f;   // selected-hit parent local angles (-999.f if no parent)
  };

  // One entry per track (refit or passthrough).
  struct SmartPixelsRefitTrackInfo {
    uint8_t status = 0;         // bit0 refit performed (else passthrough); bit1 seedCovOK;
                                // bit2 seedCovMode==parametrized; bit3 anyWindowTruncated; bits4-7 reserved
    uint8_t nCrossings = 0;     // valid layer crossings attempted
    uint8_t nAcceptedHits = 0;  // hits accepted into the KF
    uint8_t nKFUpdates = 0;     // scalar-update groups applied (== layers updated)
    uint8_t layerHitMask = 0;   // accepted-hit bitmask, bit0=L1 .. bit3=L4;
                                // popcount(layerHitMask) == nAcceptedHits (exact, v0.1)
    uint16_t maxWindowMult = 0; // max windowMult over this track's crossings
    float chi2IncRPhiTot = -999.f, chi2IncRZTot = -999.f;  // sums over crossings
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
  }  // namespace hitflag

  // Bit accessors for SmartPixelsRefitTrackInfo::status (spec §2).
  namespace trackstatus {
    inline constexpr uint8_t kRefitPerformed = 0x1;   // bit0
    inline constexpr uint8_t kSeedCovOK = 0x2;        // bit1
    inline constexpr uint8_t kParametrizedSeed = 0x4; // bit2
    inline constexpr uint8_t kAnyWindowTruncated = 0x8; // bit3
  }  // namespace trackstatus

}  // namespace smartpixels

#endif
