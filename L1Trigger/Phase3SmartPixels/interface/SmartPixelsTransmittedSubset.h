#ifndef L1Trigger_Phase3SmartPixels_SmartPixelsTransmittedSubset_h
#define L1Trigger_Phase3SmartPixels_SmartPixelsTransmittedSubset_h

// -*- C++ -*-
//
// SmartPixels transmitted-subset (hardware-boundary) model. Authoritative
// contract: L1Trigger/Phase3SmartPixels/doc/RefitSidecarSpec.md (spec v0) §3.
//
// This header is the SINGLE SOURCE OF TRUTH for the quantizers and the 16-bit
// "compact" (TS1) summary word. ngtagger-train MUST reproduce these bit-exactly
// in python when training on "compact"-tier features (same discipline as the
// trkquality two's-complement track-word decode). Header-only, constexpr where
// possible, no CMSSW dependencies -- so it can be mirrored trivially.
//
// v0 compact-word bit layout (spec §3):
//   bits 0-3   : per-layer accepted-hit bitmask (L1..L4); popcount == nAcceptedHits
//   bits 4-7   : q(chi2IncRPhiTot)
//   bits 8-11  : q(chi2IncRZTot)
//   bits 12-14 : occ = clamp(floor(log2(1 + maxWindowMult)), 0, 7)
//   bit  15    : reserved (0)
// with q(c) = clamp(round(2 * log2(1 + c)), 0, 15).
//
// The canonical entry point is packCompactWord(trackInfo) (spec v0.1): the
// SmartPixelsRefitTrackInfo now carries the exact layerHitMask and maxWindowMult,
// so no per-crossing state is needed. The explicit 4-arg form is the low-level
// primitive.

#include <cmath>
#include <cstdint>

#include "L1Trigger/Phase3SmartPixels/interface/SmartPixelsRefitSidecar.h"

namespace smartpixels {

  // 4-bit chi2 quantizer: q(c) = clamp(round(2 * log2(1 + c)), 0, 15).
  // Negative / sentinel inputs clamp to 0 (log2 domain guard).
  inline uint8_t quantizeChi2(double c) {
    if (!(c > 0.))
      return 0;
    double v = std::round(2.0 * std::log2(1.0 + c));
    if (v < 0.)
      v = 0.;
    if (v > 15.)
      v = 15.;
    return static_cast<uint8_t>(v);
  }

  // 3-bit window-occupancy quantizer: occ(m) = clamp(floor(log2(1 + m)), 0, 7).
  inline uint8_t quantizeOccupancy(unsigned m) {
    double v = std::floor(std::log2(1.0 + static_cast<double>(m)));
    if (v < 0.)
      v = 0.;
    if (v > 7.)
      v = 7.;
    return static_cast<uint8_t>(v);
  }

  // Assemble the 16-bit compact word from an explicit L1..L4 accepted-hit
  // bitmask (bit0=L1 .. bit3=L4), the two per-track chi2 totals, and the
  // per-track maximum window multiplicity. This is the exact-fidelity form: the
  // caller (the producer, which owns per-crossing layer identity) supplies the
  // true bitmask so popcount(mask) == nAcceptedHits holds physically.
  inline uint16_t packCompactWord(uint8_t layerHitMask,
                                  double chi2IncRPhiTot,
                                  double chi2IncRZTot,
                                  unsigned maxWindowMult) {
    uint16_t w = 0;
    w |= static_cast<uint16_t>(layerHitMask & 0xF);              // bits 0-3
    w |= static_cast<uint16_t>(quantizeChi2(chi2IncRPhiTot)) << 4;  // bits 4-7
    w |= static_cast<uint16_t>(quantizeChi2(chi2IncRZTot)) << 8;    // bits 8-11
    w |= static_cast<uint16_t>(quantizeOccupancy(maxWindowMult) & 0x7) << 12;  // bits 12-14
    // bit 15 reserved (0)
    return w;
  }

  // Canonical entry point (spec v0.1): assemble the compact word directly from a
  // SmartPixelsRefitTrackInfo, which carries the exact layerHitMask and
  // maxWindowMult. popcount(layerHitMask) == nAcceptedHits holds by construction.
  inline uint16_t packCompactWord(const SmartPixelsRefitTrackInfo& t) {
    return packCompactWord(t.layerHitMask, t.chi2IncRPhiTot, t.chi2IncRZTot, t.maxWindowMult);
  }

}  // namespace smartpixels

#endif
