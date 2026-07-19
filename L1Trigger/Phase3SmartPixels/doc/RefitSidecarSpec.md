# SmartPixels Refit Sidecar — data-model and adapter contract (spec v0)

Status: DRAFT v0 (2026-07-19). Companion to `PixelAVAngleResponseSpec.md`; same contract
discipline: consumers and producers MUST match this document exactly, and any change
requires a version bump here first.

Scope: defines where all NEW SmartPixels information (per-layer hit residuals, angles,
pulls, chi2 increments, window occupancy — information that does not exist in the
OT-only TTTrack) lives, and how every downstream consumer (refit TkQuality BDT, NG jet
tagger, E2E vertexing, DispVertexFinder, NanoAOD/ngtagger-train) reaches it. The design
goals, in priority order:

1. NO new track class. `TTTrack<Ref_Phase2TrackerDigi_>` is the only track type any
   official module ever sees.
2. One persistent home per fact. Downstream views are adapters over the sidecar, never
   copies that can drift.
3. Zero edits to official (non-SmartPixels) files. All new code lives in
   `L1Trigger/Phase3SmartPixels` (+ our nano package).
4. FPGA fidelity from day one: full-precision fields for training/diagnostics are
   separated from a quantized "transmitted subset" that models what could actually
   cross the hardware boundary, selectable by knob.

---

## 1. The 1:1 output-sync invariant (load-bearing contract)

The digiRefit/refit producer emits EXACTLY one output track per input track, in input
order; entry `i` of every output product is derived from input track `i`. Passthrough
fallback entries are byte-identical copies of the input track.

- This invariant is what makes "external synced collections" free: PF/Puppi candidate
  track refs (`l1TrackIdx`), `L1DispVertex.firstIndexTrk/secondIndexTrk`, and the
  reference-collection truth table all remain valid ROW-WISE against the refit
  collection and its sidecar, in both coopt and coexist workflows.
- REQUIREMENT: the producer MUST assert `outputTracks->size() == inputTracks->size()`
  and `sidecar sizes == inputTracks->size()` before `put()`. Violation throws
  `cms::Exception("SmartPixelsSyncBroken")`.
- Truth flows through the index: truth association is looked up on the INPUT track and
  applied to output row `i`. No truth product keyed to the refit collection is created.

## 2. Persistent product: `SmartPixelsRefitSidecar`

One product per refit track collection (same module, same event), instance label equal
to the producing module's variant suffix conventions (default instance `""` on the
prompt producer, `"Extended"` semantics follow the track collections). ROOT dictionary
via `src/classes.h` + `src/classes_def.xml` in this package.

```cpp
struct SmartPixelsRefitHitInfo {   // one entry per LAYER CROSSING attempted (not per accepted hit)
  uint8_t  layer;          // TBPX layer 1..4
  uint32_t detId;          // module rawId of the crossing (diagnostics)
  uint16_t windowMult;     // digis collected in window (post readout-order truncation)
  uint8_t  flags;          // bit0 hitAccepted; bit1 windowTruncated (maxHitsPerWindow hit);
                           // bit2 hasAlpha; bit3 hasBeta; bits4-7 reserved
  // --- selected hit, valid only when hitAccepted (else sentinel -999.f) ---
  float resX, resY;        // selected-hit local residual vs predicted crossing [cm]
  float cotAlphaMeas, cotBetaMeas;   // synthesized measured angles
  float sigAlpha, sigBeta; // per-hit angle sigmas from the PixelAV payload
  float pullX, pullY, pullAlpha, pullBeta;  // KF pulls r_k/sqrt(S_k) from the scalar updates
  float chi2IncRPhi;       // sum over this crossing's scalar updates of r^2/S, x + alpha terms
  float chi2IncRZ;         //                                              y + beta  terms
  // --- TRUTH-ONLY fields (never hardware-available; excluded from every transmitted subset) ---
  int8_t selHitClass;      // simlink class of the selected hit: 0 sameTP, 1 otherTP, 2 noise, -1 none
  float  parCotAlpha, parCotBeta;   // selected-hit parent local angles (-999.f if no parent)
};

struct SmartPixelsRefitTrackInfo { // one entry per track (refit or passthrough)
  uint8_t  status;         // bit0 refit performed (else passthrough); bit1 seedCovOK;
                           // bit2 seedCovMode==parametrized; bit3 anyWindowTruncated;
                           // bits4-7 reserved
  uint8_t  nCrossings;     // valid layer crossings attempted
  uint8_t  nAcceptedHits;  // hits accepted into the KF
  uint8_t  nKFUpdates;     // scalar-update groups applied (== layers updated)
  float chi2IncRPhiTot, chi2IncRZTot;  // sums over crossings
};

struct SmartPixelsRefitSidecar {
  std::vector<SmartPixelsRefitTrackInfo>            trackInfo;  // size == N tracks
  std::vector<std::vector<SmartPixelsRefitHitInfo>> hitInfo;    // outer size == N tracks
};
```

Rules:
- Sentinel for any unavailable float is `-999.f`; consumers MUST test `> -900.f`.
- Passthrough tracks have `status` bit0 unset, empty `hitInfo[i]`, zeros elsewhere.
- Parameter/covariance deltas (Δd0, Δpt, ...) are NOT stored: they are derivable from
  (output track i, input track i) via the invariant. Adapters compute them on demand.
- Fields marked TRUTH-ONLY exist for training labels and diagnostics; no inference
  adapter may forward them (enforced by the transmitted-subset definitions below).

## 3. Transmitted subsets (hardware-boundary model)

What crosses the track->correlator boundary in a real system is bits, not floats. Every
inference-side adapter MUST implement the `transmittedSubset` knob:

- `"score"` (TS0): nothing beyond the refit-quality BDT score already embedded in the
  track word MVA field / `trkMVA1`. 0 extra bits.
- `"compact"` (TS1): a 16-bit summary word (provisional v0 layout, pending a spare-bit
  audit of `TTTrack_TrackWord`):
  - bits 0-3   : per-layer accepted-hit bitmask (L1..L4); popcount == nAcceptedHits
  - bits 4-7   : `q(chi2IncRPhiTot)`
  - bits 8-11  : `q(chi2IncRZTot)`
  - bits 12-14 : `occ = clamp(floor(log2(1 + maxWindowMult)), 0, 7)`
  - bit  15    : reserved (0)
  with the shared quantizer `q(c) = clamp(round(2 * log2(1 + c)), 0, 15)`.
- `"full"` (TS2, studies-only upper bound): all non-truth fields at float precision.

REQUIREMENTS:
- The quantizer functions live in ONE header in this package
  (`interface/SmartPixelsTransmittedSubset.h`) and nowhere else; ngtagger-train MUST
  reproduce them bit-exactly in python when training on `"compact"`-tier features
  (same discipline as the trkquality two's-complement track-word decode).
- NanoAOD always writes FULL fidelity plus the materialized compact word, so offline
  training can use either tier and verify the quantizer round-trip.
- Physics gain curves quoted for hardware MUST state the subset tier; `"full"`-tier
  results are upper bounds by definition.

## 4. Adapters (the only consumers of the sidecar)

1. Refit TkQuality BDT (in-producer or adjacent module): features from
   `trackInfo` + `hitInfo` + (output − input) parameter deltas, gated by
   `transmittedSubset`; score written into the `trkMVA1` ctor slot + track-word MVA
   bits of OUR output collection. No new product.
2. Candidate glue producer `l1tSmartPixCandExtraProducer`: input = final Puppi
   candidate collection + refit tracks + sidecar (+ reference tracks in coexist).
   Resolution: `cand -> trackRef -> index -> sidecar row` (refs point at the refit
   collection in coopt; at the reference collection in coexist, mapped by the
   invariant). Output: one parallel per-candidate record collection, FULL fidelity,
   aligned index-wise to the candidate collection (neutral/unmatched candidates get
   sentinel rows). Subset selection + quantization happen in the consumers, via the
   shared header — one adapter serves all studies.
3. Cloned NG tagger producer (our package): consumes standard Puppi candidates + the
   glue product, applies `transmittedSubset`, packs extended constituent features,
   evaluates the retrained model, emits the STANDARD NGJet output type. Certification
   gate: with zero extra features and the stock model it MUST reproduce stock NGJets
   bitwise (clustering is untouched by construction).
4. NanoAOD (our nano package): `L1TSmartPixelsRefitHit` link-style table (one row per
   crossing, `trackIdx` column, the `L1SC4NGJetCands` pattern) + extension columns on
   the refit track table from `trackInfo` + the compact word + truth labels applied by
   index from the reference truth table.

## 5. Versioning and provenance

- This spec carries a version (v0). The structs carry no schema version field; the
  EDM/ROOT class version in `classes_def.xml` is bumped in lockstep with this doc.
- Any producer filling the sidecar MUST already record its payload provenance
  (fitSmartHitPayloads / PixelAV set descriptions); the sidecar itself adds none.
- Reserved for future versions (do NOT improvise): endcap layers (TEPX/TFPX ids),
  per-hit KF gain snapshots for `gainMode="lut"` studies, PixelAV high-res position
  fields (`spx_pos_*` era), per-candidate subset words materialized in the track word.
