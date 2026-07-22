# PixelAV → SmartPixels angle-response payload specification (v0)

Interface contract between the PixelAV simulation + regression-NN team and the CMSSW
`digiRefit` consumer (`L1Trigger/Phase3SmartPixels`, Tier-2 interim refit). The deliverable
described here is the `pixelavAngleSet` payload referenced by `DIGIREFIT_DEFAULTS` in
`python/customizeSmartPixels_cff.py`.

Reference implementation / self-test: `python/validatePixelAVAngleSet.py`
(`--write-example` emits a structurally complete example payload; the same script validates
any candidate payload against this spec).

## 1. Role in the refit

For every in-window pixel digi with a truth-linked parent, the digiRefit producer computes
the parent particle's TRUE local incidence angles at the module and synthesizes the
"SmartPixels-measured" angles as

```
cot(alpha)_meas = cot(alpha)_true + bias_alpha(inputs) + N(0, sigma_alpha(inputs))
cot(beta)_meas  = cot(beta)_true  + bias_beta(inputs)  + N(0, sigma_beta(inputs))   [variant-gated]
```

with a per-hit validity throw `p_valid(inputs)` (an invalid throw degrades the hit to
position-only in the Kalman update; interacts with the `useAngles` fidelity handle).
The payload encodes the PixelAV + NN **response**: resolution (sigma), bias, and validity
probability of the NN angle estimate as functions of the true angles, layer, and local
magnetic field.

### 1a. Relationship to Stack A (`smarthit_true`) — CHARACTERIZES, does not DRIVE

The angle response above is the one part of the broader Stack A "smarthit_true" family
(true-hit efficiency, position residuals, angle sigma/bias — see
`python/fitSmartHitPayloads.py`) that the Tier-2 refit actually consumes, and it is consumed
via this dedicated `pixelavAngleSet` payload, NOT via a `smarthitTrueSet`. The rest of Stack A
**characterizes** true hits (it measures what a true hit looks like) but does **not drive**
Tier-2 synthesis: position comes from the real CMSSW pixel digis (digi/cluster fidelity), and
the angle comes from this PixelAV response applied to the truth-linked parent's true incidence
angle. There is therefore no `smarthitTrueSet` input to the digiRefit producer.

The `smarthitTrueSet` config key is nonetheless kept (present + validated) in
`DIGIREFIT_DEFAULTS`, **RESERVED** for a future use: a SmartPixels ASIC on-chip readout
inefficiency model (`smarthit_true_eff`) — a hit-loss mechanism internal to the SmartPixels
chip that the CMSSW digitizer cannot express, and that would gate/weight otherwise-present
true hits. Wiring it in requires (a) a hardware-derived payload and (b) a collaboration-agreed
semantics decision — hard gate (drop the hit) vs soft weight (down-weight in selection).
Until then the producer loads no Stack A payload; if `smarthitTrueSet` is set it emits a loud
`edm::LogWarning` (category `SmartPixelsStackAUnused`) and ignores it.

## 2. Conventions (MUST match exactly)

- **Local frame**: CMSSW module-local frame of the GeomDet surface (`surface().toLocal`).
  Local z is the sensor normal (depth), local x and y span the sensor plane. For TBPX
  planks, local x ≈ global r-phi direction, local y ≈ global z direction.
- **Angles (PixelAV convention)**: for a particle with local momentum (p_x, p_y, p_z):
  `cotAlpha = p_x / p_z`, `cotBeta = p_y / p_z`. All payload angle quantities are in these
  dimensionless cot units.
- **Local B**: `MagneticField->inTesla(globalModuleCenter)` rotated into the local frame.
  Signed components in Tesla.
- Charge sign is NOT an input: its effect enters through the signed cotAlpha and signed
  bLocalY (Lorentz drift direction).
- **SOURCE-CONVENTION WARNING (verified 2026-07-19; analytically CONFIRMED 2026-07-22 by
  the PixelAV author's coordinate transformation)**: PixelAV uses test-beam coordinates;
  the CMSSW frame relates to them by `x_cms = -y_pix, y_cms = -x_pix, z_cms = -z_pix`
  (proper rotation, det = +1). Since momenta transform like coordinates,
  `cotAlpha_cms = (-p_y_pix)/(-p_z_pix) = +cotBeta_pix` and
  `cotBeta_cms = +cotAlpha_pix`: for the COT angles this is a PURE SWAP with NO sign
  flip (the z negation cancels in both ratios). Non-ratio quantities do NOT enjoy the
  cancellation: local positions/residuals swap AND negate (`dx_cms = -dy_pix`), as do
  local-B components (`B_y_cms = -B_x_pix`) — relevant if position residuals or a real
  bLocalY axis are ever sourced from PixelAV-side dumps. The PixelAV
  simulation files and the SmartPixels regression-NN eval dumps label the angles
  OPPOSITE to this spec. In those files the fine 50 um pitch runs along PixelAV
  local-y (`y-midplane` span 50 um vs `x-midplane` 200 um), so their
  `cotBeta = n_y/n_z` is the r-phi BENDING angle. Mapping any PixelAV-side source
  into this spec therefore requires the swap:
  spec `cotAlpha` (bending) <- source `cotBeta`/`cotB*`;
  spec `cotBeta` (non-bending) <- source `cotAlpha`/`cotA*`.
  Additionally the NN dump residual columns are (true - pred); this spec's bias is
  median(pred - true), so the sign is negated on extraction. Extractors MUST apply
  this mapping in exactly one clearly-marked place (see
  ngtagger-train eval_spixel_angles/extract_pixelav_angle_payload.py,
  `SWAP_ALPHA_BETA`). Payload JSONs themselves are ALWAYS in spec convention;
  consumers are unaffected. Verified consistently across the Mlp_Slim-2bit,
  Conv2D_Max-2bit, and Conv1D_Full-2bit variants. UPDATE (2026-07-22 adversarial
  convention audit, ngtagger-train eval_spixel_angles/convention_audit/): the Lorentz
  ORIENTATION is now settled empirically without a two-sign bLocalY scan — the CMSSW
  cluster x-extent minimum sits at cotAlpha ~ −0.14..−0.19, the same SIGNED location as
  the PixelAV bending-response peak (−0.15) under the pure-swap mapping, excluding the
  mirrored (swap+sign-flip) alternative. Known coverage limitation from the same audit:
  ~44% of real TBPX crossings exceed the source beta coverage (|cotBeta| > 1.07) and
  land on clamped edge sigmas (beta pull widths grow to ~5 at high |eta|); the remedy is
  extended PixelAV simulation coverage at high |cotBeta|, not a mapping change.

## 3. Payload format

One **correctionlib schemav2** JSON (optionally .json.gz) `CorrectionSet` containing exactly
these corrections (names are load-bearing):

| name | output | constraints |
|---|---|---|
| `spx_angle_alpha_sigma` | sigma of (cotAlpha_NN − cotAlpha_true) | > 0, finite |
| `spx_angle_alpha_bias`  | mean of (cotAlpha_NN − cotAlpha_true)  | \|bias\| < 1.0 |
| `spx_angle_beta_sigma`  | sigma of (cotBeta_NN − cotBeta_true)   | > 0, finite |
| `spx_angle_beta_bias`   | mean of (cotBeta_NN − cotBeta_true)    | \|bias\| < 1.0 |
| `spx_angle_valid_prob`  | P(NN emits a usable angle estimate)     | in [0, 1] |
| `spx_angle_prng`        | deterministic N(0,1) deviate for the ALPHA throw: `hashprng` `stdnormal`, entropy order (`layer,cotAlpha,cotBeta,bLocalY`) | same 4-input signature |
| `spx_angle_prng_beta`   | deterministic N(0,1) deviate for the BETA throw: `hashprng` `stdnormal`, DISTINCT entropy permutation (`cotBeta,bLocalY,layer,cotAlpha`) ⇒ independent of alpha (matches the legacy producer's two independent draws) | same 4-input signature |
| `spx_angle_valid_flat`  | deterministic U(0,1) gate variate: `hashprng` `stdflat`, entropy order REVERSED (`bLocalY,cotBeta,cotAlpha,layer`) to decorrelate from both throws | same 4-input signature |
| `spx_angle_alpha_final` | fused-shift terminal: bias-table structure whose cells are Formula `"<bias> + prngAcc"` (same rounded bias constants as `spx_angle_alpha_bias`) | 5-input signature (4 + `prngAcc`) |
| `spx_angle_beta_final`  | ditto for beta | 5-input signature |

plus **exactly these CompoundCorrections** (the synthesis throw, factorized out of
compiled code):

| name | stack | ops | inputs | meaning |
|---|---|---|---|---|
| `spx_angle_alpha_smear` | [`spx_angle_alpha_sigma`, `spx_angle_prng`] | output_op `*` | 4 | sigma × N(0,1) (two-piece term; visualization + cross-validation) |
| `spx_angle_beta_smear`  | [`spx_angle_beta_sigma`, `spx_angle_prng_beta`] | output_op `*` | 4 | sigma × N(0,1) |
| `spx_angle_alpha_shift` | [`spx_angle_alpha_sigma`, `spx_angle_prng`, `spx_angle_alpha_final`] | `inputs_update=["prngAcc"]`, input_op `*`, output_op `last` | 4 + `prngAcc` | **fused** bias + sigma × N(0,1) |
| `spx_angle_beta_shift`  | [`spx_angle_beta_sigma`, `spx_angle_prng_beta`, `spx_angle_beta_final`] | ditto | 4 + `prngAcc` | **fused** bias + sigma × N(0,1) |

Fused-shift mechanism: the consumer passes `prngAcc = 1.0`; `inputs_update` folds each
stack output into `prngAcc` via `input_op "*"` (1 → sigma → sigma·z); the terminal
`spx_angle_X_final` returns bias + prngAcc; `output_op "last"` emits it.

**CONSUMER CONTRACT** (PRIMARY — one evaluate per angle; the producer implements exactly
this, no in-code RNG for angle synthesis):

```
cot(X)_meas = cot(X)_true
              + spx_angle_X_shift(layer, cotAlpha, cotBeta, bLocalY, 1.0)
              for X in {alpha, beta}      (the trailing 1.0 is prngAcc — REQUIRED)
accept the synthesized angles iff
    spx_angle_valid_flat(inputs) < spx_angle_valid_prob(inputs)
```

Equivalent two-piece form (identical bit-for-bit — same rounded bias decimals, same prng
nodes; kept for visualization and cross-validation):
`cot(X)_meas = cot(X)_true + spx_angle_X_bias(inputs) + spx_angle_X_smear(inputs)`.
Either form reproduces a throw ~ N(bias, sigma) per bin.

HashPRNG semantics (deliberate, documented): the throw is a pure hash of the input
tuple — **identical inputs give identical throws** (a reproducibility feature:
replay/debug gets bit-identical synthesis; float-distinct angles give independent
throws). Alpha and beta throws are INDEPENDENT via distinct entropy permutations,
matching the legacy producer's two independent engine draws.

Sensor **variants** (alpha-only NN vs alpha+beta NN, pitch/thickness options, orthogonal-angle
"beta" sensor) are delivered as **separate CorrectionSet files**, not extra axes — mirroring
the existing SmartPixels correction-set pattern. An alpha-only variant still ships the beta
corrections; set `spx_angle_valid_prob`'s beta usability via the consumer's `useAngles`
config (i.e., beta entries may be filled with any positive sigma; they will not be evaluated
when `useAngles="alpha"`).

### Inputs — identical names, types, and ORDER for all five corrections

| # | name | type | valid domain | notes |
|---|---|---|---|---|
| 1 | `layer`    | int (category) | 1, 2, 3, 4 | TBPX layer |
| 2 | `cotAlpha` | real | [−0.6, +0.6] | see measured ranges below |
| 3 | `cotBeta`  | real | [−6.0, +6.0] | grazing tracks at high \|eta\| |
| 4 | `bLocalY`  | real | [−4.0, +4.0] T | signed local-y B component |

- Binning nodes MUST use `"flow": "clamp"` so out-of-domain queries (tails, sentinel-free
  but extreme tracks) return the edge value instead of erroring.
- The consumer guarantees it never evaluates at sentinel values (unlinked digis carry
  cotAlpha/cotBeta = −999 in the derivation ntuple and are excluded before evaluation).
- `bLocalX`, `bLocalZ` are intentionally NOT inputs in v0 (measured < 1e−3 T in TBPX, see
  below). The endcap extension (TEPX/TFPX, B mostly along local z) will revise this input
  list under new correction names — do not overload these.

## 4. Measured input ranges (D121 RelVal, 50-event OT-track sample; quantiles per layer)

From `SmartPixelsPayloadAnalyzer` (q01/q99 quoted; simulate with margin to the domain edges):

| layer | cotAlpha q01..q99 (max extent) | cotBeta q01..q99 (max extent) |
|---|---|---|
| 1 | −0.27 .. +0.52 (−0.29, +0.55) | −5.0 .. +5.1 (−6.1, +5.9) |
| 2 | −0.13 .. +0.39 (−0.15, +0.44) | −3.2 .. +2.9 (−4.2, +3.8) |
| 3 | −0.17 .. +0.27 (−0.18, +0.31) | −1.9 .. +2.0 (−2.4, +2.5) |
| 4 | −0.13 .. +0.22 (−0.15, +0.27) | −1.5 .. +1.5 (−1.7, +1.8) |

Local B (all four TBPX layers, 4343 crossings): `bLocalX` and `bLocalZ` within ±1e−3 T;
`bLocalY = −3.811 T` with spread < 1.2e−3 T. (Sign is our frame convention for these
modules; modules with flipped local frames would show +3.811 — hence the signed axis.)

## 5. Required PixelAV simulation grid

- `cotAlpha`: −0.6 → +0.6, step ≤ 0.05 (finer near the Lorentz-peak region if the response
  varies rapidly there).
- `cotBeta`: −6 → +6; step ≤ 0.1 for |cotBeta| ≤ 1, ≤ 0.25 for 1 < |cotBeta| ≤ 3,
  ≤ 0.5 beyond. (Cluster length grows ~linearly with |cotBeta|; the tails matter for L1.)
- `bLocalY`: v0 minimum = both signs of the nominal field, i.e. {−3.81, +3.81} T.
  Recommended: {±3.6, ±3.7, ±3.81, ±3.9} T to cover field nonuniformity and future
  geometries; single-point fits are acceptable v0 (encode as one clamp bin centred on
  ±3.81 — the validator accepts this).
- Per layer if the sensor/electronics differ by layer; otherwise duplicate one fit across
  the four layer categories.

## 6. Provenance metadata (REQUIRED)

Top-level CorrectionSet `description` and each correction's `description` must record:
PixelAV version/tag, NN training identifier, sensor variant name, assumed E-field/bias and
temperature, the B grid actually simulated, date, and the git revision of the fitting code.
The payload is refused in review if provenance is missing.

## 7. Consumer-side usage (as it will appear in the digiRefit producer)

```cpp
#include "correction.h"   // correctionlib

auto cset = correction::CorrectionSet::from_file(
    edm::FileInPath(cfg.getParameter<std::string>("pixelavAngleSet")).fullPath());
corrAlphaShift_  = cset->compound().at("spx_angle_alpha_shift");   // fused bias + sigma×N(0,1)
corrValidProb_   = cset->at("spx_angle_valid_prob");
corrValidFlat_   = cset->at("spx_angle_valid_flat");               // U(0,1), HashPRNG
// per truth-linked digi, all inputs from the module-local frame — NO in-code RNG:
const std::vector<correction::Variable::Type> in{layer, cotAlphaTrue, cotBetaTrue, bLocalY};
if (corrValidFlat_->evaluate(in) < corrValidProb_->evaluate(in)) {
  const double meas = cotAlphaTrue
      + corrAlphaShift_->evaluate({layer, cotAlphaTrue, cotBetaTrue, bLocalY, 1.0});
  // deterministic per input tuple; the trailing 1.0 is the REQUIRED prngAcc seed
  // -> angle measurement enters the Kalman update (subject to useAngles)
}
```

## 8. Validation protocol

```bash
# emit + self-check the reference example
python3 validatePixelAVAngleSet.py --write-example spx_angle_response_example.json
# validate a candidate payload
python3 validatePixelAVAngleSet.py YOUR_PAYLOAD.json
```

The validator asserts: schema v2; the ten correction names AND the four compound names
(all load-bearing); exact input names/order/types on every node (including the trailing
`prngAcc` on the finals/shifts); `hashprng` nodetype with the required distributions and
PAIRWISE-DISTINCT entropy orders; compound stack membership + ops (`output_op "*"` for
the smears; `inputs_update=["prngAcc"]`, `input_op "*"`, `output_op "last"` for the
shifts); the fused ≡ two-piece identity `shift(in,1.0) == bias(in) + smear(in)` to 1e-12
over a probe grid; determinism across repeated evaluates and fresh loads; evaluation over
the full documented grid (including out-of-domain clamp probes at cotBeta = ±10,
bLocalY = ±5) with sigma > 0, |bias| < 1, p_valid ∈ [0,1], all finite.
CI-style: exit code 0 = spec-compliant.

**Back-compat is intentionally hard**: payloads without the prng/compound nodes FAIL
validation (and the producer refuses them at load). We control all payloads; regenerate
with `ngtagger-train/eval_spixel_angles/extract_pixelav_angle_payload.py` — the additions
are purely additive, so regenerated payloads keep the plain corrections bit-identical.

## 9. Reserved for future revisions (do not use these names yet)

- `spx_pos_x_sigma` / `spx_pos_x_bias` / `spx_pos_y_sigma` / `spx_pos_y_bias` — the
  PixelAV high-resolution *position* regression (the "mini-regression working back from the
  embedded hit + simlink"; deferred by plan).
- `spx_angle_*_endcap` — TEPX/TFPX extension with revised B inputs.
