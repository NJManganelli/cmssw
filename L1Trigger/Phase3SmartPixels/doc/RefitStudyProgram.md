# v2.6 study program: sensor bit allocation, window sizing, and the two architectures

Status: **plan agreed, no producer code written for it yet.** Canonical home for
this program; the copy that briefly lived in the (unversioned) host
`planningAndPatches/` directory is retired. Companion contract:
`RefitSidecarSpec.md`. Study tools and results live in the `ngtagger-train`
repo (`docs/`, `eval_refitq/`) and are cross-referenced per section. Written 2026-09-02 after the
latency budget and the sensor payload constraint were made explicit; several
conclusions recorded earlier in the v2.6 work are corrected below.

## Version discipline: v2.6 IS NOT MINTED YET

Everything in this program lands **inside v2.6**. The commits already on the
branch that carry "v2.6" in their subject (the 4-way chi2 split and the spec
changelog) are **work in progress**, not a shipped schema: nothing is pushed and
no production has ever written a v2.6 file. Consequences, all load-bearing:

1. **ONE ClassVersion bump for the whole of v2.6.** HitInfo 4->5, TrackInfo
   4->5, Sidecar 5->6 already cover it. Adding the reco payload, the innovation
   sigmas and the projected chi2 fields must NOT bump again — re-pin the
   ClassVersion 5/5/6 checksums instead. There is no intermediate schema to stay
   compatible with, and inventing one would recreate exactly the version sprawl
   the `spec v0.x` retirement just cleaned up.
2. **The changelog carries a single consolidated v2.6 entry**, not one entry per
   commit. It is rewritten as the content settles, and only frozen when v2.6 is
   minted (= first production written with it).
3. **v2.6 is minted when a production is written with it**, and not before.
   After that point any schema change is v2.7 and the normal
   bump-and-keep-old-versions discipline applies.

---

## 0. Two corrections to earlier conclusions

**0a. There is no unquantized measurement anywhere in the real chain.** The
sensor emits ~16 bits per cluster of unique information (NOT counting fixed
sensor/module addressing), nominally x, y, cotAlpha and possibly cotBeta. That
budget follows from expected bandwidth of a design that does not exist yet, so
the split between fields is exactly what this program must determine.

The v2.6-era recommendation "keep the angle chi2 at full precision because the
refit BDT sits inside the producer, ahead of any transmission boundary" is
therefore **wrong**. The BDT never sees an unquantized sensor output. Worse, the
bit-width study behind it (`ngtagger-train docs/refit-chi2-bitwidth-study.md`)
measured float x/y/alpha/beta taken from simulation, so its headline result —
"the angle chi2 is information-rich and needs ~8 bits to encode" — is an upper
bound conditioned on a sensor that will not be built. With 2-3 bit sensor
angles the derived chi2 carries much less information and needs far fewer bits.

The two questions are ORDERED and were done backwards: sensor allocation is
upstream and dominant; encoding of derived quantities is downstream and cheap to
redo. `ngtagger calibrate-chi2-quant` already takes the allocation as input, so
the downstream answer is a rerun once an allocation is fixed, not new work.

**0b. The value of alpha/beta was measured only through the refit chi2, which
undervalues them.** In the architecture below, the pixel-only preprocessing stage
runs with no OT track and no OT stub information at all. Cluster incidence
angles are then one of the only handles available for track-agnostic filtering
(e.g. rejecting clusters whose angle is inconsistent with origin in the luminous
region) and for building pixel-only candidates. Bits spent on angles may pay off
mostly THERE, not in the refit's chi2. Any allocation study that scores only
refit-stage metrics will therefore systematically under-value beta.

---

## 1. The binding constraint: latency

| stage | budget | information available |
|---|---|---|
| total L1 latency | 12.5 us (11 us already consumed) | — |
| OT-only track building | ~5 us | OT stubs |
| **pixel-only preprocessing** | **~5 us, concurrent with the above** | **SmartPixels hits ONLY — no OT tracks, no OT stubs** |
| **match + refit** | **~1 us** | OT tracks + the preprocessed pixel product |

Three consequences that reshape every design question below.

1. **Combinatorial work is affordable only in the 5 us pixel-only phase.** That
   is where branching, candidate building and clustering can live. The 1 us
   match/refit phase can afford scoring, not searching.
2. **The primary data reduction cannot be track-guided.** The current producer
   model — project the track, open a window, collect digis, truncate at
   `maxHitsPerWindow=8` — is not the real architecture. In the real system the
   pixel side must be reduced BEFORE any track exists; the window is then a
   property of the *match* stage acting on an already-reduced product.
3. **`maxHitsPerWindow` is a latency knob, not just a fidelity knob.** Every
   candidate that survives into the match stage costs time in the 1 us budget.
   Window tightening therefore buys latency AND purity simultaneously, which
   raises the value of the window work in section 3 considerably.

---

## 2. The two major architectures

**Architecture A — refit on OT-only tracks (the current line).**
Pixel-only preprocessing during OT track building, then match and refit.
Two sub-variants for the preprocessing product:

- **A1, reduced hit lists.** Preprocessing organizes and truncates clusters
  (per module / per region), optionally filtering on angle consistency with the
  luminous region. The match stage opens a window per track per layer and scores
  candidates hit-by-hit. This is closest to what exists today.
- **A2, pixel-only candidates ("tracklets"/"tracksters") with branching.**
  Preprocessing builds multi-layer pixel candidates — e.g. 4 candidates from
  1xL1 + 1xL2 + 2xL3 + 2xL4 hits — and the refit is then evaluated against each
  candidate AS A SET rather than hit-by-hit.

  A2 is architecturally attractive for a reason worth stating explicitly: it
  moves the branching that the TMTT KF performs *during* fitting (full state
  branching over stubs per layer, chi2-pruned) into the phase that has time for
  it. The spec's own Q1 comparison already identified greedy per-layer selection
  as digiRefit's deliberate simplification versus TMTT; A2 recovers branching
  without spending the 1 us budget on it.

**Architecture B — inclusive SmartPixels + OT seeding and building.**
Drop OT-only track building; spend the ~5 us + ~1 us headroom on joint pixel+OT
seeding and fitting from scratch. More capable (no seed/refit split, no
extrapolation penalty, IT hits inform the seed rather than correct it) and more
expensive.

Both need physics-performance comparison. That comparison is the terminal
deliverable of this program, and it must be made at matched latency, not at
matched algorithmic ambition.

---

## 3. Measured evidence in hand (2026-09-01/02)

All from the only POST-guard v2.5 sidecar production (`clamp_on_f{1,2}`, 200
events PU200 TT, config AAAA, 32 799 refit tracks), plus existing v2.5 nano.

### 3a. The identity that makes v2.5 nano reusable

`pull[k] = r/sqrt(S)` and `chi2inc[k] = r*r/S` come from the same r, S behind the
same gate, so `chi2Inc<D>Tot == sum(pull<D>^2)` exactly — verified to 3.7e-8
median relative deviation against the stored joint columns. So the v2.6 split
totals are recoverable exactly from v2.5 files, and they are the same quantity
as `REFIT_BDT_FEATURES` features 5-8.

Also: because `aLin = a` is captured before the scalar-update loop, the
relinearization term vanishes for the first dimension, so
**`sqrt(S_x) = resX/pullX` exactly** — the innovation sigma is measurable from
existing nano with no new production.

### 3b. The static search windows are mis-sized in both directions

| layer | sqrt(S_x) median | static r-phi window | window/sigma median | at p90 |
|---|---|---|---|---|
| L1 | **307 um** | 500 um | **1.6x** | **0.8x** |
| L2 | 14.4 um | 1700 um | 118x | 7.8x |
| L3 | 15.4 um | 5000 um | 324x | 263x |
| L4 | 12.5 um | 9000 um | **722x** | 624x |

z/local-y: L1 1419 um vs 4500 um window (3.2x); L2-L4 ~50 um vs 2000-3500 um
windows (40-70x).

Two structural facts behind this:

- **The covariance collapses ~20x after the first pixel hit.** At L1 the
  prediction still carries the OT extrapolation; once one pixel measurement is
  in, S falls to approximately sigma_meas and stays there. The correct window
  changes by 20x within a single track between its first and second layer, so no
  per-layer constant can be right for both.
- **The pT dependence runs opposite to what a constant provides.** At L1,
  sqrt(S_x) is 342 um at 2-3 GeV and 75 um above 20 GeV, so the fixed window is
  1.5 sigma for soft tracks and 6.6 sigma for stiff ones — tightest exactly
  where the prediction is worst.

Note these per-layer S values are **order-dependent**: L1 shows seed-level
uncertainty only because it is processed first (the loop is hardcoded
inside-out, `L1SmartPixelsTrackProducer.cc:1651`). Seed-level S at each layer
requires single-layer runs; see 4b.

### 3c. Consequence, and the dominant failure mode

| layer | <windowMult> | wrong (otherTP) fraction of accepted hits |
|---|---|---|
| L1 | 6.1 | **47.7%** |
| L2 | 6.1 | 25.4% |
| L3 | 5.6 | 19.6% |
| L4 | 5.1 | 17.8% |

**54.7% of refits include at least one wrong hit** (0.82 wrong of 2.72 accepted
per track). Window area scales as sigma_x*sigma_y, so correctly sized windows at
L2-L4 would shrink area by ~1e4 and take contamination there towards zero. This
is a larger effect than any plausible bit-allocation gain.

### 3d. The error model is optimistic, and probably for an identifiable reason

Pull widths for hits truthfully matched to the track's own TP (robust sigma;
should be 1.0):

| dim | correct hits | wrong hits | ratio |
|---|---|---|---|
| pullX | 2.04 | 2.27 | **1.1** |
| pullY | 1.61 | 6.29 | 3.9 |
| pullAlpha | 2.20 | 5.81 | 2.6 |
| pullBeta | 1.30 | 46.3 | **35.6** |

29.8% of *correct* x-hits sit beyond 3 sigma. digiRefit has **no multiple
scattering term**, which the spec's Q1 table already flagged against TMTT
(`sigmaScat = KalmanMultiScattTerm/pT`, 0.00075). Its absence would inflate L1
pulls (extrapolation-dominated) while barely touching L2-L4
(measurement-dominated) — which is the observed pattern. The hand-tuned L1
window being 1.6x the modeled sigma looks like silent compensation for the same
gap.

**Ordering consequence: the MS term must be fixed before windows are derived
from the covariance**, or an optimistic S will produce windows that lose real
hits at exactly the rate the optimism implies.

Note also the position error model is `sigX = sigY = pitch/sqrt(12)` — digital
single-pixel resolution with no cluster interpolation. Bits spent on x/y ARE
that interpolation, so section 4a and the error model are coupled.

### 3e. Wrong-hit discrimination is where the four-way split pays

Single-feature AUC for "every accepted IT hit came from the track's own TP":
Beta 0.879, Y 0.768, X 0.753, Alpha 0.706, all-four-summed 0.941. Versus the
deployed `genuine` objective: 0.717 / 0.768 / 0.617 / 0.649 / 0.728.

The seed chi2 and the refit deltas are **complementary, not redundant**: seed
chi2 adds +0.0434 on `genuine` over the deltas but +0.0001 on `clean`; the
deltas add +0.0033 on `genuine` but **+0.3139** on `clean`. So the BDT wants
both, and the deployed `genuine` label is nearly blind to the failure mode the
refit actually controls.

Correction recorded in the spec: the original v2.6 rationale — that summing
angle into position chi2 dilutes the discriminant — is false. The plain 4-way
sum is the best single scalar (0.941). The split is still right, because the
angle terms are the strongest wrong-hit discriminants and the old (r-phi, r-z)
pairing hid them by mixing each with a position term.

### 3f. What wrong hits cost, and why a trust gate is the headline result

Measured 2026-09-02 (`eval_refitq/wronghits/param_vs_wronghits.py`), 31 775
refit tracks with a genuine TP match. Population by number of wrong hits used:

| n_wrong | 0 | 1 | 2 | 3 | 4+ |
|---|---|---|---|---|---|
| share | 46.7% | 35.4% | 12.3% | 4.1% | 1.5% |

Median |param - truth|, seed -> refit, by category:

| n_wrong | d0 [um] | z0 [um] | phi [mrad] | tanL [1e-3] | pt [GeV] |
|---|---|---|---|---|---|
| **0** | 247 -> **30.5** (+88%) | 1005 -> **87** (+91%) | 1.60 -> 0.81 (+50%) | 3.34 -> 1.41 (+58%) | +7% |
| 1 | 529 -> 576 (-9%) | 1683 -> 1744 (-4%) | 3.16 -> 5.78 (-83%) | -148% | -127% |
| 2 | 567 -> 840 (-48%) | 1517 -> 2556 (-69%) | -298% | -302% | -549% |
| 3 | 532 -> 1104 (-108%) | 1298 -> 2920 (-125%) | -633% | -426% | -1312% |
| 4+ | 430 -> 1328 (-209%) | 1116 -> 2812 (-152%) | -978% | -399% | -2055% |

**When the refit picks only correct hits it improves d0 by 8x and z0 by 11.5x,
and a single wrong hit annihilates the entire gain.** 91.9% / 94.7% of clean
refits improve d0 / z0; that falls to 45.5% / 48.5% at one wrong hit.

Note the categories are not random subsets: the seed error is itself ~2x worse
for n_wrong >= 1 (d0 529 vs 247 um), i.e. poor seeds attract wrong hits. That is
the expected consequence of 3b - a worse seed means larger S, hence a wider
correct window, hence more candidates - and it means window sizing and hit
purity are the same problem.

Population-level policy comparison (same tracks, three policies):

| | d0 median | d0 q95 | phi median | phi q95 | pt median |
|---|---|---|---|---|---|
| always seed | 350 um | 1547 | 2.245 mrad | 16.9 | 0.0295 GeV |
| always refit | 132 um | **2365** | 2.329 | **46.5** | **0.0504** |
| **oracle gate** | **118 um** | **1513** | **1.581** | 16.9 | 0.0286 |

As a blanket replacement the refit is a **mixed-to-negative** deal: it improves
the d0/z0 median but degrades every tail, and degrades phi/tanL/pt at nearly all
quantiles. With a perfect cleanliness gate it becomes a **strict improvement on
all five parameters at every quantile** (d0 and z0 ~3x better median, phi and
tanL ~30% better, pt neutral-to-better).

**Design consequence.** The refit-quality classifier's highest-value job is not
predicting seed genuineness — it is predicting HIT CLEANLINESS, so a consumer
can take refit parameters per track when clean and fall back to the seed
otherwise. And that is achievable: the `clean` label reaches AUC 0.94 from a
single summed chi2 scalar and 0.996 from a BDT on the float 4-way chi2 (3e), so
a realistic gate should recover much of the oracle column. This also explains
the long-standing "refit quality has limited impact" impression: measured
against `genuine`, the classifier was being asked the wrong question.

---

## 4. The study program

Ordered so that each stage removes a confound from the next.

### 4a. Sensor bit allocation (16 bits/cluster)

The actuator: quantize the MEASUREMENTS (x, y, cotAlpha, cotBeta) at configurable
widths and ranges before they enter the KF. Grid over allocations summing to ~16,
including the beta-or-not question:

```
(x7, y7, a2, b0)   (x6, y6, a2, b2)   (x5, y5, a3, b3)
(x6, y7, a3, b0)   (x7, y6, a2, b1)   (x8, y8, a0, b0)   <- position-only control
```

Scored on FOUR axes, because scoring only the refit undervalues angles (0b):

1. hit-selection purity (wrong-hit rate, per layer);
2. track parameter resolution vs truth (the physics deliverable);
3. refit-quality discrimination (both labels);
4. **preprocessing value**: can pixel-only filtering / candidate building use
   these angles without any OT information? Measured in the A2 study (4d).

Then rerun `calibrate-chi2-quant` on the winning allocation to fix the
downstream chi2 encoding for the BDT.

### 4b. Error model and window sizing

1. Add the MS inflation (TMTT-style `measTerm/pT`), re-measure pull widths per
   layer and per dimension; target robust sigma ~1.0 for correct hits.
2. Single-layer runs (`smartPixelsActiveLayers = 1000/0100/0010/0001` — the knob
   already exists, no code needed) to obtain **seed-level** S per layer. This is
   the missing input for the layer-order question and costs four short jobs.
3. Window policies, three arms:
   - static per-layer constants (today);
   - **LUT(pT bin, layer)** — captures the dominant 1/pT scaling and the
     L1-vs-rest split with no matrix algebra; likely most of the gain;
   - covariance-derived `n*sqrt(S)`, or equivalently a chi2-space cut
     `dx^2/S_x + dy^2/S_y < cut` which avoids the square root. `H` depends only
     on the track state, layer and geometry — not on hits — so it can be hoisted
     above candidate collection at zero physics cost (it is currently built at
     line 1833, after selection). Needs floor/ceiling guards in the style of the
     existing `jacobianMaxAbs` / `chi2UpdateGate`.

   Metrics per policy: truth-hit capture efficiency per layer, windowMult,
   wrong-hit inclusion, truncation rate, and downstream parameter resolution.

### 4c. Layer crossing order

Knob `digiRefitLayerOrder = insideOut | outsideIn | byAmbiguity`. The loop is
currently hardcoded inside-out, and because the KF updates state as it goes, a
wrong hit at the first layer corrupts the prediction for the rest — and the first
layer is currently L1, which has both the largest S and the worst contamination
(47.7%).

Prior after 3b, stated so it can be falsified: whichever layer goes first pays
the full seed uncertainty, and in r-phi that is cheapest at L1 (the extrapolation
converges as r->0; the static windows encode 18x growth outward). So inside-out
may well be correct in r-phi while z favours the opposite, and the ordering may
be second-order next to window sizing. The single-layer runs in 4b.2 settle it
empirically. `digiRefitMaxKFUpdates_` interacts: order decides WHICH layers are
used when the cap bites.

### 4d. Architecture: window matching vs pixel-only candidate building

Three arms, aligned to section 2:

- **A1** static/LUT/covariance window matching on reduced hit lists (4b).
- **A2** pixel-only candidate building with branching (e.g. 1+1+2+2 -> 4
  candidates), refit scored against each candidate as a set. Study questions:
  what branching factor fits 5 us; how often does a truth-consistent candidate
  exist; does set-level scoring beat greedy hit-by-hit selection at equal
  latency; and how much do angles contribute to building candidates WITHOUT OT
  information (feeds 4a axis 4).
- **B** inclusive pixel+OT seeding and building, compared at matched latency.

Phase-0 sizing for A2/B is possible on existing nano: pixel hit multiplicity per
module/region, and how often a truth-consistent multi-layer pixel combination
exists at all.

### 4e. Refit fidelity vs truth-matched hits (categorized by wrong-hit count)

Three hit-selection modes, `digiRefitHitSelection` (truth modes marked
studies-only, per the spec's TRUTH-ONLY discipline):

- `window` — production behaviour;
- `tpOracle` — use the hits of the TP matched to the OT-only track, ignoring the
  window. **This is the ceiling for refit performance** and the reference every
  other arm should be quoted against;
- `windowDropWrong` — window selection but wrong hits discarded, which isolates
  the damage wrong hits do from the damage window coverage does.

Break every result out by number of wrong hits included: 0, 1, 2, 3, 4+.
Deliverable: parameter-resolution-vs-truth per category, and the population of
each category. The observational half needs no new production — see 5.

---

## 5. What is measurable today, and what needs the rest of v2.6

**Available now, no schema change, no new production:**
- innovation sigma per layer/pT (3b) — `sqrt(S_x) = resX/pullX`;
- per-layer ambiguity and wrong-hit rates (3c);
- pull widths by truth class (3d);
- **4e's observational half in full**: the nano carries refit params
  (`rInv/phi/tanL/z0/d0` on the variant table), seed params (reference table),
  truth (`tp_pt/tp_d0/tp_z0/tp_eta/tp_phi/tp_tanL`) and per-hit `selHitClass`,
  so the n_wrong = 0..4+ categorization of parameter change vs truth can be done
  immediately. Only the ORACLE arm needs a producer knob.
- Phase-0 sizing for A2/B (4d).

**Needs the rest of the v2.6 sidecar payload and knobs:**

Per-hit additions (studies): `predLocalX/Y`, `predCotAlpha/Beta` (currently only
the residual and the measured angle are stored, so prediction and measurement
cannot be separated); `sigX/sigY` and `sqrt(S)` per dimension as actually used
(removes the resX/pullX reconstruction and its relinearization caveat); the
truth-matched hit's residuals even when another hit was selected
(`truthResX/Y`, `truthInWindow`) for 4e.

Per-track float chi2 additions: `chi2StubsAtSeed`, `chi2StubsAtRefit`,
`chi2ITAtRefit`, `chi2JointAtRefit`. These are FINAL-STATE chi2 values, not the
sequential KF increments (each of which is evaluated at its own predicted
state); "the chi2 of the refit track over stubs and IT hits" means the former,
and it is what a genuine joint fit would minimize. Carry both.

Nano precision: raise the per-hit chi2/pull columns off `mantissaBits=12` for
study productions. At 2^-12 they are too coarse for the fine-structure questions
in 4a.

Knobs: `digiRefitLayerOrder`, `digiRefitHitSelection`, `digiRefitMeasBits`
(nx, ny, nalpha, nbeta + ranges), MS-term parameter, window policy selector.

---

## 6. SETTLED (2026-09-02): the stub chi2 uses the float approximation

Decision: build the **float approximation** for the study phase. Rationale and
the caveat that comes with it are below; revisit only if a conclusion about
large parameter shifts turns out to matter.

Also settled in the same pass:
- SmartPixels hits carry **float x, y, cotAlpha, cotBeta** in the sidecar
  (alongside the predicted values, so measurement and prediction separate).
- The sidecar stores the float chi2 terms from the IT hits **and** the projected
  chi2 of stubs + IT hits re-evaluated against the refit track.
- Study-3 oracle is the **TP matched to the OT-only track** (the performance
  ceiling), with `windowDropWrong` retained as the second arm.

The two candidates that were weighed:

- **TMTT-style**: re-project each OT stub onto the refit helix, residuals in
  (r, phi, z) with per-module-type sigma (PS vs 2S), summed. Physically correct,
  matches what the L1 track finder already does, and is the form a firmware
  implementation would eventually need. More code, and it imports the OT error
  model with its own conventions.
- **Float approximation**: reuse the seed's per-stub residuals and apply a linear
  correction from the parameter shift (`dchi2 ~ 2 r^T W H da + da^T H^T W H da`).
  Much less code, trivially re-derivable when the fit changes, and adequate for
  ranking/first-order studies — but it is only valid for small parameter shifts
  and will misestimate exactly the pathological large-shift cases that matter
  most for wrong-hit diagnosis.

Chosen: the float approximation, because it unblocks 4e immediately and is cheap
to discard when the fit design changes. **Standing caveat to carry with every
number it produces:** it is a small-shift expansion, so it will misestimate
exactly the large-parameter-shift cases that wrong-hit diagnosis cares most
about — and 3f shows those are common (a single wrong hit moves d0 by ~290 um
median and phi by ~2.3 mrad, with 4+ wrong hits moving pt by 0.56 GeV). So the
approximation is trustworthy for the clean and 1-wrong categories and
progressively less so beyond; any conclusion drawn in the 2+ categories must be
flagged provisional until a TMTT-style re-projection exists.

Still open: whether per-stub records go in the sidecar (bigger, enables residual
studies) or only the totals.

---

## 7. Tooling status

| tool | state | purpose |
|---|---|---|
| `ngtagger calibrate-chi2-quant` | committed (`ed59c51`) | per-field (bits, k) for chi2 delta encoding; per-field constants, LUT vs log, BDT-verified selection. Rerun per allocation in 4a. |
| `eval_refitq/quantstudy/chi2_bitwidth_study.py` | committed (`63ff147`) | the bit-width/label study of record |
| `eval_refitq/quantstudy/chi2_code_efficiency.py` | committed | label-free code entropy/saturation diagnostics |
| `eval_refitq/quantstudy/chi2_perfield_k.py` | committed | per-field BDT-isolated k scan |
| window-sigma measurement (3b) | **written, not yet committed** | fold into `eval_refitq/windows/` |
| n_wrong categorization (4e observational) | **not written** | next, no new production needed |
| measurement-quantization emulator (4a) | **not written** | offline first (quantize nano values, re-derive), producer knob second |

Deliberately still pending: the new `REFIT_BDT_FEATURES` version, the
migration of `_dataio.py` / `refit_replay.py` / two tests / five
`modelspace/*` scripts off pre-v2.6 column names, and propagation of the three
v2.6 commits to the arm mirror and the CMSSW_17 backports — all held because the
producer will change again for this program.
