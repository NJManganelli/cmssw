# SmartPixels IT seeding and fitting — design note

Status: **design study, no code and no firmware.** Nothing here has been
implemented, simulated, or synthesised. The measured inputs are real and cited;
everything the design itself concludes is arithmetic on top of them and is
marked. Companion documents: `RefitStudyProgram.md` (the v2.6 study program that
produced most of the measured inputs), `RefitSidecarSpec.md` (the refit data
contract), `PixelAVAngleResponseSpec.md` (the angle payload). Study tooling
lives in the `ngtagger-train` repo (`eval_refitq/combinatorics/`).

Written 2026-09-05, after the SmartPixelsRecHit angle work landed and the
`globalCluster{Phi,CotTheta}` nano columns were specified.

## Reading the number tags

Every quantity in this document carries one of three tags. They are not
decoration — the design's conclusions are only as good as the weakest input
feeding them, and several load-bearing numbers are currently in the third
category.

- **[measured]** — read out of a production or a study in this workspace, with a
  citation. Trustworthy.
- **[derived]** — arithmetic in this document, on top of measured or assumed
  inputs. Reproducible but unreviewed; the algebra has not been independently
  checked.
- **[assumed]** — not established anywhere in this workspace. These are the
  ones that can move a conclusion. All are collected in §12.

---

## 1. Scope, and the two physics targets

An inner-tracker-only track seeder built from SmartPixels clusters, running
**in parallel to** the Phase-2 Outer Tracker track finder (which sees only OT
stubs), inside the same 4–5 µs latency envelope. Seeds only; the fit is §10.

Two targets, which want different things and should not be conflated:

**A — prompt, pT ≥ 2 GeV, IT-seeded then fit jointly with OT stubs.**
The OT already finds these tracks; the motive for seeding them from the IT is
the seed's impact-parameter quality, not its pT. The refit measurements show
where that pays: the 4-parameter baseline pins d0 ≡ 0, while the refit *creates*
a d0 with a 0.026–0.051 cm core and sharpens z0 by ~20%. A seed built on IT
clusters starts from that quality rather than extrapolating inward to it.
Because the OT stubs are available for confirmation, **A can tolerate a high
seed fake rate** — the combined fit's χ² does the cleaning.

**B — soft, 0.5–2 GeV, not found by the OT at all.**
Below the OT track finder's pT floor these tracks have no L1 representation.
There is no OT confirmation available by construction, so **B's seeder must be
self-sufficient on fakes.** B is the demanding case and should drive the design;
A falls out as the relaxed special case.

**η reach: |η| ≲ 1.4 today — a BLOCKED DEPENDENCY, not a design choice.**

This document analyses a TBPX-only seeder because the disc sensors cannot yet be
included, **not** because forward coverage is unwanted. The blocker is the angle
payload: **there is no PixelAV angle parametrisation for the disc sensors in
their B-field configuration.** A disc sensor sits perpendicular to B where a
barrel sensor sits parallel to it, so the charge-transport and Lorentz response
that maps cluster shape to incidence angle is a different problem, not a
rescaling of the barrel one *(mechanism is my inference; the missing
parametrisation is the stated fact)*. Until it exists, a cluster on a disc has no
usable angle, and every mechanism in §3 and §4 depends on having one.

Within TBPX alone the geometric reach is **[measured]** on `spix_postq_500.root`
by study (8), cross-tabulated over the whole (pT, η) grid:

| \|η\| | 0–0.8 | 1.2 | 1.6 | 2.0–2.4 |
|---|---|---|---|---|
| layers lit by the *typical* track | **4** | 3 | 2 | **1** |
| fraction of tracks reaching ≥3 layers | 91–100% | ~73% | ~16% | **~2%** |

The numbers are the same at 1 GeV and at 20 GeV — geometry, not rate — so within
a barrel-only configuration a 3-layer seed effectively does not exist above
|η| ≈ 1.6.

Three consequences. The genuine-track counts in §7 must be read as **restricted
to |η| < 1.4**. The 4–5 µs and board-count budgets in §5 are likewise for a
barrel-only system and will grow when the discs enter. And **the η reach of this
design is bounded by a payload gap that is expected to close**, so nothing here
should be written as though forward coverage were permanently excluded — see
§11 for the dependency.

**Not in scope here:** displaced tracks. Every transform in this document
assumes d0 ≈ 0, which is what makes a 2-D Hough affordable. **[derived]** at
r = 3 cm a d0 of only 100 µm shifts the position azimuth by 3.3 mrad — roughly
5× the multiple-scattering blur at 2 GeV — so L1 is exquisitely sensitive to
displacement and the prompt assumption is not a mild one. A displaced seeder
(the disappearing-track case in `DarkSectorL1Scoping.md` §2.5) needs a third
Hough dimension or a different method, and is deliberately excluded.

---

## 2. What the sensor gives you

**[measured]**, all from `mem:smartpixels-v2p6-state` unless noted:

| quantity | value |
|---|---|
| clusters/event, PU200 TBPX | 26 479 |
| per occupied module, L1–L4 | 40.6 / 32.7 / 31.6 / 19.7 (p95 up to 72) |
| TBPX modules | 864 (216 / 216 / 180 / 252) |
| pixel pitch | 25 × 100 µm |
| sensor angle resolution, cot α | 0.0225 |
| cluster fraction from pT above 0.5 / 1 / 1.5 / 2 GeV | 33.4 / 10.4 / 4.6 / 2.2 % |
| L1 tracks/event at PU200 | 189 (`planningAndPatches/trackInputMode-matrix.md`) |

The per-module means times the module counts reproduce the event total to 0.02%,
so essentially **every module is hit every event** and those four numbers are
also the per-layer totals: ≈ 8.8k / 7.1k / 5.7k / 5.0k on L1–L4.

Two findings set the whole problem:

1. **pT is not a local observable.** No raw cluster feature separates pT > 1 or
   pT > 2: all AUC 0.44–0.55, best working point keeps 90.9% of clusters. pT
   lives in curvature, which one sensor cannot see. **You cannot filter the soft
   flood at the sensor on cluster shape alone.**
2. **The angle is the exception, and only at the outer layers.** |cot α| gives
   AUC 0.538 / 0.651 / 0.733 / 0.830 for pT > 2 on L1–L4, against raw sizeX at
   0.528 / 0.534 / 0.570 / 0.615. The gap is the headroom an on-sensor angle
   regressor competes for, and it is largest at L3/L4.

**Nano schema** (in flight, being added for this data tier). Per cluster:
`localX, localY, localCotAlpha, localCotBeta` (module frame);
`globalR, globalPhi, globalZ` (cluster **position**, cylindrical);
`globalClusterPhi, globalClusterCotTheta` (cluster **direction**, global frame);
`charge, size, sizeX, sizeY, sigX, sigY, sigAlpha, sigBeta, hasAlpha, hasBeta,
layer, detId`. Truth counterparts under a `tp` prefix: `tpIdx, tpPt,
tpChargeFrac, tpLocalCotAlpha, tpLocalCotBeta, tpGlobalClusterPhi,
tpGlobalClusterCotTheta`.

Note that `globalPhi` (where the cluster is) and `globalClusterPhi` (where the
track is going) are a confusable pair on which the entire transform depends.

**Rotated uncertainties are present**: `sigGlobalClusterPhi` and
`sigGlobalClusterCotTheta` were added, which is what makes the segment lengths of
§3 settable from stored columns. `sigGlobalClusterCotTheta` is **[measured]** 21%
optimistic (pull width 1.210) — §11 item 4.

**Superseded:** earlier versions of this analysis carried a blocking caveat that
noise clusters received no angle, making `hasAlpha` a perfect proxy for "is a
real cluster" (99.5% vs 0.0%). After the `smarthit_noise_*` re-derivation this
is 99.5% / 100% and the confound is gone. The `*** RESULTS ARE CURRENTLY
CONFOUNDED ***` block in `spix_combinatorics_omnibus.py` study (6) should be
retired accordingly.

---

## 3. The transform: each cluster is a segment, not a line

For a prompt helix, write κ = q/pT and c = 0.3·B/2 ≈ 0.0057 for B = 3.8 T with
r in cm and κ in GeV⁻¹. The position azimuth and the direction azimuth turn at
different rates — the tangent rotates by r/R, the position azimuth by half that:

    φ_p  =  φ₀ − c·r·κ           (position, `globalPhi`)
    φ_d  =  φ₀ − 2·c·r·κ         (direction, `globalClusterPhi`)

**[derived]**, three consequences:

1. **The classic Hough line** is φ₀ = φ_p + c·r·κ — one line per cluster in
   (φ₀, κ), as in any position-only transform.
2. **Each cluster measures curvature by itself**: κ̂ = sin(φ_p − φ_d)/(c·r) —
   the exact form, since φ_p lags φ₀ by asin(c·r·κ) and φ_d by twice that. The
   small-angle version is only valid for stiff tracks. Two
   measurements at one radius over-determine two parameters, so a single cluster
   nominally fixes the track. This is the gradient-informed Hough of classical
   image processing, where an edge's gradient direction collapses its vote from
   a whole curve to a point.
3. **The error ellipse is a sliver lying along the Hough line.** Its transverse
   width is σ(φ_p) (small: ~0.33 mrad at L1, ~0.06 mrad at L4) and its extent
   *along* the line is σ(κ̂) = σ(φ_d)/(c·r) (large, because the direction is the
   poorly measured quantity).

So the vote footprint is a **line segment**, centred on κ̂, of half-length
n·σ(κ̂), lying along the cluster's Hough line. That is the object this design
votes with.

**Implementation note.** Do not implement literal geometric segment
intersection: pairwise intersection is O(N²), 9 M pairs/event at N = 3000. The
segment *is* the vote footprint. The FPGA-native form is an ordinary binned
accumulator in which the angle supplies a per-cluster column range
[κ_lo, κ_hi] and a weight ramp — a small modification to a standard Hough
transform, not a different algorithm.

**Weighting.** Hard truncation at 1σ discards 32% of real clusters. Graded
weights across the σ bands (3 / 2 / 1 for < 1σ / 1–2σ / 2–3σ) let the peak
threshold do the rejecting while every real cluster still votes somewhere.
Cost is ~2 extra bits per accumulator cell.

---

## 4. What the angle actually buys — the veto, not the segment

**[derived]**, using cot α ≈ c·r·κ and σ(cot α) = 0.0225 **[measured]**:

| layer | r [cm] **[measured]** | σ(κ) from the angle | ±1σ segment, as fraction of a pT>2 window | of a pT>1 window |
|---|---|---|---|---|
| L1 | 2.9 | 1.36 | 100% | 100% |
| L2 | 6.0 | 0.66 | 100% | 66% |
| L3 | 10.3 | 0.38 | 77% | 38% |
| L4 | 14.5 | 0.27 | 54% | 27% |

This independently reproduces the measured AUC ladder: the angle informs
curvature only at large radius, because the lever arm to the beamline is what
converts an angle into a pT.

**The segment shortening is modest.** In a pT > 2 window at the ±2σ needed for
efficiency, segments are full-length at every layer — no saving at all. In a
pT > 1 window at ±1σ the saving is ~1.8×; at ±2σ, ~1.25×. **This does not move
a feasibility argument, and the design should not be sold on it.**

**The veto is where the angle pays.** 97.8% of IT clusters come from tracks
below 2 GeV **[measured]**. A 2σ angle cut keeps |cot α| < c·r/2 + 2σ, which is
a per-cluster statement that a cluster *cannot* come from a track above some pT.
**[derived]**, at the measured radii, the cut rejects below 0.96 / 0.79 / 0.55 /
0.31 GeV on L4 / L3 / L2 / L1, giving ~9.6× / 5.9× / 3.2× / 2× per layer and
**~3.3× overall**. (An earlier revision quoted ~2.9× from a dropped factor of
two in the L2/L3 thresholds; the overall figure barely moves because L1, where
the veto is weakest, carries the most clusters.)

Two caveats on that 2.9×. The retention fractions are integrated over all four
layers, and soft tracks curl, so the pT spectrum at L4 is harder than at L1 —
which makes the L4 factor optimistic. The per-layer spectrum must be measured.
And critically, an angle veto tuned for target A rejects everything below 1 GeV,
which is **half of target B's signal**; see §8.

**The η angle is the more valuable one, and it has now been measured.**
`globalClusterCotTheta` gives z₀ directly per cluster: z₀ = z − r·cot θ, so
σ(z₀) = r·σ(cot θ). The radial dependence is *opposite* to the bending case —
the bending angle informs pT only at large r, the η angle informs z₀ best at
small r — so L1 and L4 do complementary jobs rather than L1 being pure
combinatorial noise.

**[measured]** by `spix_combinatorics_omnibus.py` study (7) on
`spix_postq_500.root` (100 events, PU200, post module-flip fix), as the robust
width of (`globalClusterCotTheta` − `tpGlobalClusterCotTheta`):

| | L1 | L2 | L3 | L4 |
|---|---|---|---|---|
| median cluster r [cm] | 2.9 | 6.0 | 10.3 | 14.5 |
| σ(cot θ) | 0.0265 | 0.0260 | 0.0246 | 0.0232 |
| σ(z₀) [µm] | **799** | 1599 | 2585 | **3405** |
| σ(κ̂) measured | 1.277 | 0.621 | 0.369 | 0.278 |
| σ(κ̂) predicted in the table above | 1.36 | 0.66 | 0.38 | 0.27 |

Two closure results worth recording. **[measured]** each cluster's Hough line
passes through the true (κ, φ₀) to a median of **0.51 mrad** (p90 2.5 mrad),
which validates the position-side construction end to end. And the measured
σ(κ̂) agrees with the §4 prediction σ(κ) = σ(cot α)/(c·r) to within 6% at every
layer — so the segment lengths this design votes with are now **measured**, not
merely derived.

The σ(z₀) row is numerically identical to a set of values retracted in an
earlier revision. Those were unjustified at the time (they multiplied a *local*
angle by r); the conversion is legitimate now that cot θ genuinely is dz/dr.
They were right by accident, and are now right on purpose.

> ### RESOLVED — the module-flip bug, and what it broke
>
> **Root cause (fixed 2026-09-05):** modules are physically **flipped** within a
> ladder — a real feature of the detector layout — and the direction-propagation
> code did not account for it. On flipped modules the propagated direction came
> out inverted. One bug, two symptoms, both now fixed:
>
> | symptom | before | after |
> |---|---|---|
> | `globalClusterCotTheta` behaved as `−localCotBeta`, not dz/dr | per-TP RMS of `z − r·cot θ` = 5.95 cm, *worse* than the 3.10 cm do-nothing baseline | **0.0087 cm** against a 3.18 cm baseline — a **365× collapse** |
> | `globalClusterPhi` direction sign scrambled per module | κ sign flipped cluster-to-cluster on one track (−0.609, +0.603, −0.602, +0.467 for |κ|=0.591); π branch present on some modules | π branch **gone** (median \|φ_p−φ_d\| = 0.098, zero clusters beyond 2.0 rad); per-cluster κ sign agrees with truth **98.1%** on clean ≥3-layer topologies, flat across layers (98.3/98.3/98.1/97.5) |
>
> **[measured]** on `spix_postq_500.root` (100 events), |κ|·pT = 1.0035 — the
> bend relation κ = sin(φ_p − φ_d)/(c·r) is exact.
>
> **One of my acceptance tests was wrong and is withdrawn.** I claimed a true
> dz/dr must give a layer-independent |cot θ| across the population, and that its
> 1.78 → 0.48 fall across L1–L4 was "the decisive test". That reasoning conflates
> per-track invariance with a population median: TBPX is a barrel, so high-|η|
> tracks only reach the inner layers and L1 samples a wider η range than L4.
> **[measured]** restricted to TPs that actually reach L4 the medians are flat —
> 0.577 / 0.532 / 0.494 / 0.477. The population trend was acceptance all along.
> The valid test is the self-consistency one in the table above.
>
> **Residual, not blocking:** ~2% of clusters on clean topologies still carry the
> wrong κ sign, and agreement falls to 82.5% when messy topologies (curlers,
> multiple clusters per layer, secondaries) are included. Segment shading is
> re-enabled on that basis, but a seeder would need to handle a few-percent sign
> error rather than assume none.
>
> ### `globalClusterPhi` has a second, independent problem: the sign
>
> Its *magnitude* validates — **[measured]** |κ| = |sin(φ_p − φ_d)|/(c·r)
> reproduces the truth at |κ|·pT = 1.008 (p25–p75 0.979–1.044) for pT > 2 GeV.
> But the **orientation convention varies per module**, so the sign does not.
> On one 1.69 GeV track (expected |κ| = 0.591):
>
> ```
> L1 r= 3.25  phi_p-phi_d = -0.0113   kappa = -0.609
> L3 r=10.21  phi_p-phi_d = +3.1065   kappa = +0.603   <- pi offset present
> L4 r=14.91  phi_p-phi_d = -0.0511   kappa = -0.602
> L4 r=14.93  phi_p-phi_d = +0.0397   kappa = +0.467   <- no offset, sign flipped
> ```
>
> Some modules carry a π offset and some do not, and the residual's sign flips
> independently. `sin()` absorbs the π but not the sign. **A seeder cannot use a
> per-cluster bending angle whose sign is ambiguous**, so §3's segment
> shortening — the entire reason for preferring this over a plain Hough — is
> blocked until the convention is fixed or the module orientation is published
> alongside it.
>
> Both problems live in `smartpixels::toGlobalDirection`. Suggested acceptance
> tests: (i) for one track, `z − r·cot θ` constant across its clusters; (ii)
> |cot θ| independent of layer; (iii) sign(φ_p − φ_d) constant along a track.
>
> **What still works, and what study (8) therefore draws:** the position-only
> Hough line φ₀ = φ_p + c·r·κ needs no direction column at all and is unaffected.
> The example panels fall back to it, with the truth curvature taken as
> |κ| = 1/tpPt and its sign from how the position azimuth turns with radius.

Three things the σ(cot θ) row settles:

- **σ(cot θ) ≈ 0.023–0.026, flat in layer, and comparable to σ(cot α) = 0.0225.**
  The pessimistic estimate — that β, read from cluster length against a 100 µm
  pitch, would land near 0.7 and be useless — is wrong by a factor ~30.
- **The predicted η ≈ 0 degradation does not appear.** cot θ is read from cluster
  length, which is shortest at η ≈ 0, so the resolution was expected to be worst
  exactly where most tracks are. **[measured]** σ(cot θ) by |η| bin is
  0.0243 / 0.0213 / 0.0264 / 0.0278 / 0.0281 / 0.0288 across |η| 0–2.4 — mildly
  *better* near the middle and degrading toward the endcap. The prediction in an
  earlier revision of this document was wrong and the η binning is retained only
  because it was the right thing to check.
- **The two targets are indistinguishable**: σ(z₀) = 1.489 mm for pT > 2 GeV and
  1.485 mm for 0.5–2 GeV. The z₀ lever works just as well for the soft target,
  which is the one that needs it.

**Caveat — the stored uncertainty is optimistic.** **[measured]** the pull
(`globalClusterCotTheta` − truth) / `sigGlobalClusterCotTheta` has width **1.208**
with median −0.017: unbiased, but the claimed σ is ~21% too small. Segment
lengths and vote weights (§3) must not be set from that column until it is
recalibrated, or every window will be ~20% too tight.

### 4a. THE ABOVE MEASUREMENT IS CIRCULAR — read before using it

The σ in §4 is **not a physical resolution**. It is the PixelAV payload's own
smear, recovered:

- `SmartPixelsRecHitProducer.cc:292-319` builds `trueCotA/trueCotB` as the
  dominant TP's helix **propagated to the hit**, with no multiple scattering.
  The file says so at `:32-40`: *"WHAT THIS STILL NEGLECTS: multiple scattering
  between the vertex and the module. The helix is the no-scattering limit."*
- `:330-331` then sets `cotB = trueCotB + corrBetaShift_->evaluate(...)`.
- `L1SmartPixelsClusterTableProducer.cc:199-202` fills
  `tpGlobalClusterCotTheta` from that same `trueCotBeta()`, rotated.

So (reco − truth) is *identically* the payload shift. **[measured]** the
signature is visible in the study's own output: σ(z₀) is 1.489 mm for pT > 2 GeV
and 1.485 mm for 0.5–2 GeV, flat to 0.3%, where a physical resolution would
degrade toward low pT. `globalClusterPhi` shares the identical construction, so
§3's segment lengths rest on the same artifact.

**How much the neglected scattering adds.** A kink at radius r_s displaces the
*extrapolated* z₀ by r_s·δ(cot θ), not r·δ — material outside the measurement
radius does not bias z₀. **[derived]** at 1.5% X₀/layer **[assumed]**, added in
quadrature to the values above:

| | η = 0, 1 GeV | η = 1.6, 1 GeV | η = 2.0, 0.5 GeV |
|---|---|---|---|
| L1 (only the beam pipe inside) | +0.0% | +1% | +3% |
| L4 (beam pipe + L1 + L2 + L3 inside) | +0.1% | +2% | +27% |

**The §6 conclusion survives, and survives where it matters**: L1 dominates the
z₀ lever (93.7 slices against L4's 22.0) and is nearly immune, because almost
nothing sits inside r = 2.9 cm. The worst corner — L4, high η, soft — is +27%,
inside the ±2× bracket §6 already quotes. What does **not** survive is the
**[measured]** tag: the numbers in §4 and §6 should be read as an estimate
resting on an assumed material budget until the closure test below is done.

**A second, separate fidelity problem.** The cluster's *shape* comes from the
Geant4 digitisation, which does scatter, while its *reported angle* is
synthesised from the helix, which does not. Shape and angle are therefore
mutually inconsistent for scattered tracks, and the noise inverse-CDF's
conditioning on `sizeY` (`:353-357`) is drawn from a slightly wrong joint
distribution. This affects anything using shape and angle together, not just z₀.

**What a `vz` closure test would and would not reach.** If the cluster truth
block carried the dominant TP's production `vz`, one could compare
z₀_pred = `globalZ` − `globalR`·`globalClusterCotTheta` against the TP's actual
z₀. That is **non-circular** — its real value — but it is a *bound*, not the
answer, and on the pessimistic side. **[derived]**, for one scatter of size δ at
radius r_s measured at radius r:

    simulation:      z0_pred - vz  =  +(r - r_s)*delta  -  r*(draw)
    real detector:   z0_pred - vz  =  -r_s*delta        -  r*(draw)

because the simulated cluster's POSITION is scattered (Geant4 made it) while its
ANGLE is not (the helix made it) — the two no longer describe the same
trajectory. For a cluster on L4 scattered at L1 that is +11.6δ where the physical
term is −2.9δ: four times too large, and opposite in sign.

So the two numbers bracket the truth rather than pinning it:

| | optimistic (payload draw only, §4) | pessimistic (`vz` closure) |
|---|---|---|
| L1, η ≈ 0 | 800 µm | +0.0% |
| L1, η = 2, 0.5 GeV | 800 µm | +0.3% |
| L4, η ≈ 0 | 3404 µm | +0.2% |
| L4, η = 2, 0.5 GeV | 3404 µm | +42% |

**The bracket is tight enough to design against.** The bounds agree to under 1%
everywhere except the high-η, low-pT corner, where they span 27–42%. L1 — which
carries the z₀ lever at 93.7 slices against L4's 22.0 — is robust under either,
because almost nothing sits inside r = 2.9 cm. If the `vz` closure confirms this,
the scattering question is **closed for design purposes** without being resolved
exactly.

**PSimHit is NOT the way to resolve it.** `PSimHit::localDirection()` is the only
scattering-aware angle available today, but it is g4SimHits/"SIM", pre-mixing,
and exists for **~1.7%** of PU200 clusters — signal only. Using it as a
validation reference would characterise the scattering term on a population whose
pT and η spectra differ from pileup, and any number so derived would carry that
asymmetry into every study that consumed it. This is the same rule
`SmartPixelsRecHitProducer.cc:36-39` already enforces when it refuses PSimHit as
a production input: *"manufacturing discrimination out of a truth-access
asymmetry."* The rule applies to validation as well. By contrast `vz` is uniform
— TrackingParticles are post-mixing and cover PU and signal alike (**[measured]**
86.5% of clusters, 573 815 of 662 993, carry a TP link) — which is why it is
priority 1 despite being only a bound.

**The exact answer requires better cluster simulation, not better truth access.**
The angle must be derived from the SIMULATED CLUSTER SHAPE for every cluster —
which is what a real smart pixel sensor does, and what the on-sensor regressor is
ultimately for — instead of being injected as helix-truth plus a payload draw.
Until that exists the scattering term stays unquantified, and this document
carries it as a bracket rather than a number.

**Overlap with the `Q` process-noise work.** Same physics, different path. `Q` is
calibrated for scattering along **OT seed → IT** (targets k = 167/121/114
µm·GeV on L2/L3/L4), which crosses far more material than the **vertex → IT**
path at issue here; **[derived]** the vertex→L2 displacement is ~3× smaller,
consistent with the path difference rather than with either number being wrong.
The material integral should be derived **once** and serve both.

---

## 5. System scale — bandwidth is not the constraint

Envelope: up to the Phase-2 Track Finder's own board count (Apollo boards, dual
VU13P), time-multiplex period up to 18, silicon at Versal Gen 2 or better at the
tier VU13P occupied in its generation. Take 9 φ sectors × TMP 18 = 162 nodes.

**[derived]:**

| quantity | value |
|---|---|
| full event payload | 26 479 × 88.6 bits = 2.35 Mb |
| aggregate rate | ~94 Tb/s |
| per node (one sector-event per 450 ns) | **~580 Gb/s** |
| as a fraction of one VU13P's transceivers | ~20% |
| clusters per node-event | 2 942 |
| cycles available (450 ns at 500 MHz) | 225 |
| required input width | **~13 clusters/cycle** |
| accumulator replication to sustain it | ~13 copies ⇒ ~416 BRAM |

Against 2 688 BRAM36 on a VU13P this is comfortable, before any generational
improvement. **Retention f = 1.0 — every cluster, no readout gate at all — is
affordable on the I/O and logic axes.** Latency is likewise not binding: voting
is ~225 cycles plus pipeline depth, against 4–5 µs.

An earlier version of this analysis concluded the opposite by computing the
aggregate (~32 VU13P-equivalents) and implicitly treating a dozen boards as the
budget. At 162 nodes there is roughly 10× margin on ingest. The error is
recorded because it inverts the design's priorities: it is not a data-movement
problem.

---

## 6. What does bind: the accumulator saturates

**[derived]**, at f = 1 with a plain 2-D (φ₀, κ) transform and a pT > 2 window:

- 735 clusters per layer per sector, spread over ~698 φ₀ bins of 1 mrad (set by
  the 0.68 mrad multiple-scattering blur at 2 GeV) ⇒ **~1.05 clusters per cell
  per layer**.
- Expected random 3-of-4 coincidences per cell = 4 × 1.05³ ≈ 4.6.
- Over ~22k cells ⇒ **~10⁵ fake triples per sector per event**, against ~190
  real tracks.

A threshold of 3 is meaningless when every cell already holds ~4 clusters. The
transform needs roughly 10³–10⁴ of combined suppression. Three multiplicative
levers:

| lever | factor on fakes | cost |
|---|---|---|
| z₀ slicing from a global direction | ~30 **[assumed — BLOCKED, see §4a retraction]** | none — a free extra dimension |
| retention f | f³ | ε³ on efficiency (§9) |
| finer φ₀ bins | 2–4 | floored by multiple scattering |

**[derived]** product ≈ 2 400, which is about what is needed with nothing to
spare. The load-bearing lever is z₀ slicing; §4 now measures the σ(cot θ) it
depends on.

**How the slicing factor is actually defined.** The suppression is
**Z_range / (4σ)**, not Z_range / (slice width): a cluster must vote into every
slice its z₀ could belong to, so it occupies ~4σ/w of the w-wide slices and w
cancels. Choosing a finer slicing buys nothing on its own — only a smaller σ
does. An earlier revision of this document paired "σ(z₀) ≈ 1 cm" with "30
slices", which is internally inconsistent by 4×: 30 slices over ±15 cm requires
**σ(z₀) ≤ 2500 µm**.

**[measured]** per-layer slice counts are **93.9 / 46.9 / 29.0 / 22.0** on L1–L4.
Naively quoting the cluster-weighted overall figure (50.9) would overstate the
benefit, because a seed's z₀ agreement window is set by the *quadrature sum* of
the layers it uses, i.e. dominated by the outermost. **[derived]** the honest
per-triplet suppression, P ≈ (4σ_ij/Z)(4σ_ik/Z) over Z = 30 cm:

| triplet | suppression | vs the assumed 1/30² |
|---|---|---|
| L1-L2-L3 | 8.6 × 10⁻⁴ | 1.29× **better** |
| L2-L3-L4 | 2.0 × 10⁻³ | 1.83× **worse** |

So the design's assumed factor of 30 is **validated to within ±2×**, with inner
triplets beating it and outer triplets falling short. The §7 fake numbers are
retained unchanged; the two triplet cases bracket them. It also sharpens the
complementarity of §4 — **the inner layers are where z₀ is won**, the outer
layers where pT is.

Directly visible in study (8): on the `pt2_eta0.4` example the z₀ slice cuts a
159-cluster cone to **24**.

One caveat carries over from §4a and is not removed by the module-flip fix:
σ(cot θ) is still the PixelAV payload's own smear, so every number here remains
an **optimistic bound** on the physical resolution.

---

## 7. Four design points, and a cell-count invariance

The four thresholds are a **design scan**, not a menu of deliverables — one
survives to implementation. Sizing rules: φ₀ bin width from the scattering blur
at threshold, θ_MS ≈ 1.36/pT mrad per ~1% X₀ layer **[assumed]**; κ column width
held fixed at 0.031 GeV⁻¹ so pT resolution does not degrade as the window widens.

| | 2.0 GeV | 1.5 | 1.0 | 0.5 |
|---|---|---|---|---|
| κ window | ±0.5 | ±0.67 | ±1.0 | ±2.0 |
| φ₀ bins/sector | 1024 | 768 | 512 | 256 |
| κ columns | 32 | 43 | 64 | 128 |
| **cells** | **32.8k** | **33.0k** | **32.8k** | **32.8k** |
| fakes, relative | 1 | 2.4 | 8 | 64 |
| fakes/event (30 z₀ slices) | ~490 | ~1 160 | ~3 900 | ~31 000 |
| genuine tracks/event **[assumed]** | ~190 | ~300 | ~700 | ~2 400 |
| fake : real | ~2.6 : 1 | ~3.9 : 1 | ~5.6 : 1 | ~13 : 1 |

**The accumulator size is invariant across all four points.** This is exact, not
coincidental: φ₀ bins scale as pT_min (scattering) and κ columns as 1/pT_min
(window at fixed resolution), so the product 513·pT × 64/pT cancels. The wobble
at 1.5 GeV is integer rounding only.

The practical consequence is worth protecting as a design constraint: **one
accumulator geometry, one firmware, four bin-map constant sets.** Different TMP
slices could run different working points against the same fabric. If a later
choice breaks the invariance it costs four firmwares instead of one.

Fakes scale as (2/pT_min)³ exactly, for the same reason the retention lever is
cubic — both act through per-cell density.

**Where the knee is.** 1.5 and 1.0 GeV are affordable at f = 1, needing no
readout gate at all. **0.5 GeV is the outlier** — 10× worse than 1.0, and the
only point that requires the retention lever to work (**[derived]** f ≈ 0.6 with
a better-than-random ranking to reach ~3 : 1). The scan's real question is
therefore narrow: *is 0.5 reachable, or does the system stop at 1.0?*

Two things could move the knee against 0.5 GeV specifically. cot θ is measured
from cluster length, which is shortest near η ≈ 0, so if the z₀ slicing degrades
there it degrades hardest for the softest working point. And at 0.5 GeV the
helix model is strained by scattering and energy loss over the IT path, so the
φ₀ bin width may need to exceed the single-layer θ_MS used above.

---

## 8. Two-stage architecture

The plausible implementation is not a single seeder but two: **unified IT+OT
seeding at > 2.0 GeV**, plus a **recovery stage on clusters the first stage did
not use**, reaching down to 1.5, 1.0 or 0.5 GeV.

**The cluster-removal benefit is not real, and should not be in the rationale.**
**[derived]** stage 1 finds ~190 tracks, each claiming ~4 IT clusters: ~760 of
26 479, or **2.9%**. The recovery stage faces essentially the full PU200 density.

**What the architecture does buy is window partitioning**, which is worth more
than the removal would have been. The recovery stage searches a κ *annulus*
(0.5 < |κ| < 1/pT_min) rather than a full window, and stage 1 keeps fine φ₀ bins
instead of being coarsened to accommodate soft tracks. **[derived]:**

| recovery down to | 1.5 GeV | 1.0 | 0.5 |
|---|---|---|---|
| recovery cells | 8.4k (0.26×) | 16.4k (0.50×) | 24.6k (0.75×) |
| recovery fakes/event | ~296 | ~1 942 | ~23 300 |
| **two-stage total** | **~786** | **~2 432** | **~23 790** |
| single-stage equivalent | ~1 160 | ~3 900 | ~31 000 |
| advantage | 0.68× | 0.62× | 0.77× |

The two-stage design is therefore **both cheaper and lower-fake than the
equivalent single wide-window design**, at 1.26–1.75× the stage-1 accumulator
rather than 2×, across all three recovery choices.

**Run the two stages in parallel, not serially.** "Recovery on unutilised
clusters" implies masking, which forces stage 2 to wait for stage 1 to find
peaks, fit them and claim hits — serialising two full pipelines inside 4–5 µs.
Since masking removes only 2.9%, that is real latency for almost nothing.
Running both accumulators concurrently on disjoint κ annuli and de-duplicating
at the end captures the entire window-partitioning benefit with no
serialisation; duplicate suppression is needed regardless.

**The readout gate cannot be shared the way the accumulator can.** Target A's
best gate (angle veto below 1 GeV, ~9.6× at L4) removes *all* of the 0.5 GeV
signal band and half of the 1.0 GeV band. Options: a union gate tuned to the
softest threshold (~3×, giving up most of A's reduction — which A does not need,
having OT confirmation); a pT-band-aware gate; or per-TMP-slice gates. The union
gate looks like the cheap answer and should be checked rather than assumed.

---

## 9. Retention as a scan axis, not a number

Which ranking function a front-end would use is unresolved and under study, so
the design must not depend on it. The robust formulation uses **two** numbers:

- **f** — fraction of *all* clusters kept. Sets bandwidth and combinatorics.
- **ε** — fraction of *needed* clusters kept. Sets efficiency.

A ranking function's entire value is the gap it opens between them. Total charge
has ε < f — worse than random: **[measured]** top-N-by-charge per module keeps
0.7 / 1.7 / 4.4 / **12.7** / 36.8 / 75.4 % of needed clusters at N = 1 / 2 / 4 /
8 / 16 / 32 (`spix_combinatorics_omnibus.json`, `charge_gate.survival_vs_topN`,
n = 9647), because high-pT tracks cross near normal incidence and make short,
low-charge clusters. Random has ε = f. A useful ranking has ε ≫ f.

Scaling, **[derived]**, for a k-layer seed:

| quantity | scaling |
|---|---|
| bandwidth, votes, logic | ∝ f |
| fake seeds | ∝ f^k |
| seed efficiency | ∝ ε^k |
| **net figure of merit** | **(ε/f)^k** |

Cutting retention is only worth doing to the extent the ranking beats random,
and the payoff is that advantage cubed. At k = 3, ε/f = 3 buys 27× purity at
fixed efficiency.

**This is why 2-layer seeds matter.** Requiring k-of-k is what makes ε^k brutal;
relaxing topology buys it back. **[derived]** at ε = 0.5: 3-of-3 = 12.5%,
3-of-4 = 31%, 2-of-4 = 69%. Low retention and multi-layer requirements are in
direct tension, and 2-layer seeds are the release valve — at the cost of being
the fake-dominated ones. The useful output is the crossing point: at what f does
the 2-layer efficiency gain overtake its fake cost?

**Study construction.** Every result becomes a curve over f, bracketed by two
reference gates that require no ranking to be known:

- **random gate** (ε = f): uniform subsample per module — pessimistic bound.
- **oracle gate** (ε = 1): keep truth-linked clusters first, fill the rest at
  random — optimistic bound.

Real rankings are points between the curves, droppable on later without
redesigning the study. Note the ε^k arithmetic assumes independent per-cluster
survival, which is false — a track's clusters are correlated in whatever the
ranking keys on, so tracks tend to lose all or none, making real efficiency
*better* than ε^k. Another reason to measure the bracket rather than model it.

---

## 10. Fitting

Deliberately thin, because seeding is what this note analyses. What is settled:

- **Target A fits jointly with OT stubs.** The IT seed supplies d0/z0 quality,
  the OT supplies the lever arm. The fit's χ² is also what makes A tolerant of a
  2.6 : 1 seed fake rate.
- **Target B fits IT-only**, over a 13 cm lever arm (r = 3 → 16 cm). **[derived]**
  this is not resolution-limited: the sagitta at 2 GeV is ~1200 µm against ~10 µm
  hit resolution, and multiple scattering dominates, giving roughly constant
  ~1% relative pT resolution. IT-only curvature measurement is physically sound;
  the difficulty is entirely in finding the right clusters, not in fitting them.
- The existing `digiRefit` is **not** this fit: it refits *existing* OT-seeded
  tracks and its layer order (`outsideIn`) is chosen for an OT-only seed
  projecting inward. An IT-seeded fit inverts that geometry and must not inherit
  the ordering conclusion without re-deriving it.

Open: whether an IT seed feeds the existing refit as a new `trackInputMode`, or
whether seeding and fitting become one producer. Not decided here.

---

## 11. What must be measured before this hardens

In priority order. The first two can move the design's conclusions; the rest
tighten it.

**PARTLY DONE — σ(cot θ) → σ(z₀), per layer and per η bin.** Was priority 1;
study (7) ran and is reported in §4 and §6, but §4a shows the residual it
measures is circular — reco and truth differ by construction only by the PixelAV
draw. The lever still looks real and the assumed factor of 30 still holds to
±2× once the neglected scattering is estimated, but that estimate rests on an
assumed material budget. Superseded by item 1 below. By-products that do stand:
the layer radii are genuinely measured, the feared η ≈ 0 degradation does not
exist, and `sigGlobalClusterCotTheta` is 21% optimistic.

**DONE — the module-flip bug.** Was priority 1. Modules are physically flipped
within a ladder and the direction propagation did not account for it; both
`globalClusterCotTheta` (behaving as `−localCotBeta`) and `globalClusterPhi`
(scrambled sign) were symptoms of that one cause. Fixed and verified — see the
§4a box. Residual: ~2% wrong-sign clusters on clean topologies.

1. **Add the dominant TP's production `vz` to the cluster truth block**, and
   redo study (7) as a z₀ closure test (§4a). It removes the circularity and is
   **uniform across PU and signal** — TrackingParticles are post-mixing, so
   **[measured]** 86.5% of clusters carry one. It yields a *pessimistic bound*,
   not the exact resolution, because the simulated cluster pairs a scattered
   position with an unscattered angle. That is acceptable: the bound plus §4's
   optimistic figure bracket the answer to under 1% except in the high-η, low-pT
   corner, which is enough to design against. One column. Highest value here.
   **Do NOT use `PSimHit` for this** — signal-only at ~1.7%, it would characterise
   scattering on a biased population and carry that asymmetry into everything
   downstream, which is the same failure `SmartPixelsRecHitProducer.cc:36-39`
   already refuses for production angles.
2. **PixelAV angle parametrisation for the DISC sensors, in the disc B-field
   configuration.** This is the single dependency standing between the design and
   |η| > 1.4 (§1). A disc sensor is perpendicular to B where a barrel sensor is
   parallel to it, so the barrel payload cannot be reused. Until it lands, disc
   clusters carry no usable angle and the seeder is barrel-only by necessity.
   Everything else in this note is written for TBPX and will need its geometry,
   occupancy and board-count numbers revisited when the discs enter — the
   analysis structure carries over, the numbers do not.
3. **Derive the cluster angle from the simulated cluster SHAPE**, for every
   cluster, instead of injecting helix-truth plus a payload draw. This is the
   only thing that resolves §4a exactly rather than bracketing it, it is what a
   real sensor does, and it is where the on-sensor regressor is headed anyway.
   Long-horizon; until it exists the scattering term stays a caveat, not a number.
4. **Recalibrate `sigGlobalClusterCotTheta`.** **[measured]** pull width 1.208,
   median −0.017 — unbiased but too small. Until fixed, no segment length or
   vote weight in §3 may be set from that column. Check
   `sigGlobalClusterPhi` the same way; §3's segment *lengths* depend on it and it
   has not been looked at.
5. **A single IT material model**, shared with the `Q` process-noise work.
   §4a's scattering estimate and `Q`'s k = 167/121/114 µm·GeV calibration are the
   same physics integrated over different paths (vertex→IT vs OT seed→IT). Derive
   the material budget once and let both consume it, rather than each carrying
   its own assumption.
6. **Per-layer pT spectrum of clusters.** The §4 angle-veto factors used
   all-layer retention, which is optimistic at L4 where the spectrum is harder.
7. **Efficiency denominator.** ε is undefined until "needed cluster" is defined,
   and the choice is not neutral: `mem:smartpixels-v2p6-state` warns that "TP
   lights ≥3 IT layers" measures IT traversal, not OT reach — a 0.3 GeV track has
   R ≈ 26 cm and still crosses r = 16 cm. Targets A and B need different
   denominators.
8. **Statistics.** Studies (7) and (8) now run on 100 events / 2.7 M clusters,
   ample for resolution widths but still not for anything rate-like; the
   fake-rate and efficiency studies of §7 and §9 need a larger production.
9. **Whether a front-end readout gate exists at all**, i.e. clusters per module
   per BX shippable at 40 MHz. §9 is constructed so the design does not depend
   on the answer, but the answer decides whether the ranking question is urgent
   or academic.
10. **Whether TBPX modules are tilted** in this geometry, which decides whether
   the local→global angle rotation is per-module or a per-layer constant.

**Study home.** These belong in
`ngtagger-train/eval_refitq/combinatorics/spix_combinatorics_omnibus.py` as
additional `study_*` functions — the script is explicitly designed to grow this
way ("add a function, add it to STUDIES, done — do not fork this script"). Two
adaptations are needed: the seeding studies want per-event *all-cluster* arrays
(`K`) rather than the per-crossing join (`P`), which the existing signature
already permits by ignoring `X`/`P`; and `_discover()` currently hard-requires a
refit hit table, which a clusters-only input would not have. Outputs become
arrays over f rather than scalars.

The script can honestly produce votes/event, accumulator occupancy, peak
multiplicity, candidates/event, efficiency and fake rate — **and their tails**,
which matter more than means for firmware: per-module occupancy is already
p95 = 72 against a mean of 40.6 **[measured]**, and buffers size to p99.9. It
cannot produce LUT or BRAM counts; any resource model must be labelled
separately in the JSON and not blended with measured quantities.

---

## 12. Assumptions ledger

Everything the design rests on that is **[assumed]**. Each is a place where a
conclusion could move.

| assumption | value used | consequence if wrong |
|---|---|---|
| ~~TBPX layer radii~~ | **[measured]** 2.9 / 6.0 / 10.3 / 14.5 cm (median cluster r) — retired; the §4 table's assumed 3.0/6.8/10.9/16.0 was ~10% high at L2/L4 | — |
| σ(cot θ) ⇒ 30 z₀ slices | **[measured]** 0.023–0.026 post-flip-fix; per-triplet suppression within ±2× of the assumed 30. Still an OPTIMISTIC BOUND: §4a's circularity is unaffected by the flip fix | **the design's largest single factor**; §11 items 1 and 3 close it |
| material per IT layer | ~1.5% X₀ | sets θ_MS (φ₀ bin counts in §7) **and** the §4a scattering estimate — now load-bearing in two places |
| genuine tracks/event above threshold | 190 / 300 / 700 / 2400 | moves the fake : real ratios, not the fake counts |
| PU z spread | ~5 cm | sets how much the z₀ slicing actually buys |
| clock | 500 MHz | §5 cycle budget |
| φ sectors × TMP | 9 × 18 = 162 nodes | §5 per-node numbers scale inversely |
| a front-end readout gate exists | unresolved | §9 is built to survive either answer |
| σ(φ_d) ≈ `sigAlpha` | untilted-barrel approximation | segment lengths in §3 |

Two further cautions that are not assumptions but modelling limits: the ε^k
efficiency scaling ignores correlated cluster loss (§9), and the whole document
assumes d0 ≈ 0 (§1).
