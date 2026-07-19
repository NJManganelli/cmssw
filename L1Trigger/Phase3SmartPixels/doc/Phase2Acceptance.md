# digiRefit (Tier-2) Phase 2 acceptance criteria

Authoritative plan: `mem:smartpixels-tier2-refit-plan`. These criteria are fixed
BEFORE implementation so the validation is not fit to the result. "OT-only fit"
means the input TTTrack as produced by the tracklet chain (== `passthrough`).

## 1. Reproducibility (hard gate)

- Two identical `digiRefit` runs (same input, same config, same seed policy)
  produce **bitwise-identical** refit TTTrack collections. Randomness comes only
  from `RandomNumberGeneratorService`: each producer label gets a deterministic,
  label-derived seed (`_ensureDigiRefitRNG` in `customizeSmartPixels_cff.py`,
  `TRandom3`), and draws come from the service's per-stream engine
  (`getEngine(streamID)`) in `produce()`.
- Scope of the guarantee: reproducibility holds for identical jobs (same input
  file(s), same event set and order, same `nStreams`). It is **not**
  event-order-independent: splitting the input, skipping events, or multi-stream
  scheduling changes the draw sequence. Per-event deterministic seeding
  (hash of label/run/lumi/event) is planned alongside the Phase-3 producer
  changes, where split-job training productions first need it.

## 2. Passthrough invariance (hard gate)

- A `("passthrough", None)` variant run next to a `digiRefit` variant emits
  tracks **identical** to the standard tracklet collection (rInv/phi/tanL/z0/d0
  exactly equal), as already verified in the WF1 smoke. digiRefit must not
  perturb the passthrough path.

## 3. Matched-track resolution (primary physics gate)

For TP-matched tracks with `tp.pt() > 2 GeV`, over a fixed sample, comparing the
refit track parameters to the matched TP truth:

- (a) **Improvement over OT-only.** With `useAngles=alphaBeta` and the loosest
  truncation (`maxHitsPerWindow` and `maxKFUpdates` at their configured caps),
  the RMS of `Δd0`, `Δz0`, and `Δpt/pt` versus the TP is **<=** the OT-only RMS
  in the same sample (improvement or at worst parity; d0/pt are where the IT
  angle measurements carry information). Reported as a before/after table split
  by `tp_track_match` class, reusing the producer's existing `track_summary_`
  accumulators (extended with an `_new` set already present).
- (b) **Monotone degradation in `useAngles`.** RMS(Δd0), RMS(Δpt/pt) satisfy
  `none >= alpha >= alphaBeta` (position-only refit is the weakest; adding alpha
  then beta only helps or is neutral). `none` must be approximately a
  position-only refit (angles disabled), distinctly different from the input yet
  weaker than `alpha`.
- (c) **Monotone degradation as truncation tightens.** Reducing
  `maxHitsPerWindow` and `maxKFUpdates` from the caps toward 1 does not improve
  resolution (RMS is non-decreasing). `gainMode="lut"` is a reserved placeholder
  and MUST raise (NotImplemented) rather than silently degrade.

## 4. Knob liveness (functional gate)

Each FPGA-fidelity knob, toggled individually, changes the output in the
expected direction without crashing:

- `useAngles`: `none` != `alpha` != `alphaBeta` (refit collections differ).
- `maxHitsPerWindow`: 1 vs 8 changes which hits attach (collection differs when
  windows are multiply-occupied — expected in TTbar jets).
- `maxKFUpdates`: capping below the active-layer count changes outputs.
- `gainMode="lut"`: raises `cms::Exception` (NotImplemented placeholder).

## 5. Loud-failure guard (hard gate, folds old task #14)

- A `digiRefit` (or any truth-required mode) job fed a **broken truth config**
  (e.g. wrong `MCTruthTrackInputTag` label so the association map resolves to an
  empty/mismatched product, or `truthSource=fromFile` with maps absent) must
  throw `cms::Exception` after a reasonable window of events with **zero** truth
  matches, rather than silently emitting unmodified passthrough tracks. The
  guard counts matched tracks across the stream and throws in `endStream()` (or
  once a per-job event threshold is crossed with zero matches).

## 6. Fake separation (BDT-input gate; PU-dependent)

- With a PU200 sample, the digiRefit per-hit / per-track BDT-input distributions
  (n IT hits attached, per-layer position + angle pulls, IT chi2 increment, Δd0,
  ΔpT, window occupancy) for TP-matched vs fake tracks **separate visibly**
  (populations distinguishable by eye in 1-D projections of the strongest
  features).
- **noPU caveat:** `/work/testfiles` currently has only
  `RelValTTbar_D121_noPU.root`, which contains **zero** unmatched (fake) L1
  tracks (recorded in the plan memory). Therefore criterion 6 is **deferred**
  until the PU200 sample (user-fetching) is available; Phase-2 validation
  demonstrates criteria 1-5 on noPU and the BDT-input columns are produced and
  sane, but the matched-vs-fake separation is asserted only once PU200 lands.

## Sample & config for validation

- Input: `/work/testfiles/RelValTTbar_D121_noPU.root`, posture B (`truthSource=inJob`,
  `DIGI:pdigi_valid,L1TrackTrigger,...`) so digis/simlinks/maps are self-consistent.
- Angle payload: PixelAV stand-in from
  `python/validatePixelAVAngleSet.py --write-example` until the real fit arrives.
- Seed collection: both `nPar=4` (standard prompt) and `nPar=5` (extended prompt
  + covariance, default) are run for the A/B comparison.
