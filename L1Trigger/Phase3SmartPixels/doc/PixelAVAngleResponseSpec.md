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
corrAlphaSigma_ = cset->at("spx_angle_alpha_sigma");
corrAlphaBias_  = cset->at("spx_angle_alpha_bias");
corrValidProb_  = cset->at("spx_angle_valid_prob");
// per truth-linked digi, all inputs from the module-local frame:
const std::vector<correction::Variable::Type> in{layer, cotAlphaTrue, cotBetaTrue, bLocalY};
if (rng.flat() < corrValidProb_->evaluate(in)) {
  const double meas = cotAlphaTrue + corrAlphaBias_->evaluate(in)
                      + rng.gauss(0., corrAlphaSigma_->evaluate(in));
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

The validator asserts: schema v2; the five correction names; exact input names/order/types;
evaluation over the full documented grid (including out-of-domain clamp probes at
cotBeta = ±10, bLocalY = ±5) with sigma > 0, |bias| < 1, p_valid ∈ [0,1], all finite.
CI-style: exit code 0 = spec-compliant.

## 9. Reserved for future revisions (do not use these names yet)

- `spx_pos_x_sigma` / `spx_pos_x_bias` / `spx_pos_y_sigma` / `spx_pos_y_bias` — the
  PixelAV high-resolution *position* regression (the "mini-regression working back from the
  embedded hit + simlink"; deferred by plan).
- `spx_angle_*_endcap` — TEPX/TFPX extension with revised B inputs.
