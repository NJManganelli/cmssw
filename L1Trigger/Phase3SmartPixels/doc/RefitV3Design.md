# Seedless refit v3: design note

v3 is the port of the Python replay's calibrated arm INN-D with the two-layer hit choice
(ngtagger-train `eval_refitq/replay/`, the reference; its specification is exported with the code
snapshot, the payloads and reference outputs). Select it with the producer's
`seedlessRefitPreset = "v3"` (vocabulary name `seedlessRefitV3`, nano tag `@SPSeedlessRefitV3`); the
v2 path (`seedlessRefitPreset = "v2"`, the default) is unchanged.

**Why v3 was adopted.** The physicist declared GO on 2026-10-20 over the replay scan's pre-declared
NO-GO: the results are good enough with placeholder angles, which are to be replaced by the real
sensor's angles; v3 is built so that this replacement is a payload change (doc/RefitV3Calibration.md).

## Invariants

- **Angles are used only locally:** each angle enters as the measured minus the predicted module-frame
  cot at the crossing. There is no beamline or vertex pointing term, and no consumer may add one.
- **The refit never reads a cluster's class** (T or C). The class exists only in the emulated sensor
  (`SmartPixelsRecHitProducer`) and in the TRUTH-ONLY sidecar labels the producer fills after the
  engine has chosen.
- **Every calibrated number is in the versioned payload set** `data/refit_v3/<set>/` (manifest with
  sources, registry ids and md5, checked at load). The angle use per consumer, angle, region and visit
  step is configuration (`seedlessRefitV3.angleUse`, `angleSteps`).

## Code

- `interface/SmartPixelsRefitV3.h`, `src/SmartPixelsRefitV3.cc`: the engine, framework-free. Geometry
  from the reconstruction geometry (`src/SmartPixelsRefitV3Geometry.cc`), in double precision: v3 does
  not use the float32 projector of v2.
- `plugins/L1SmartPixelsTrackProducer.cc`: `seedlessRefitPreset` / `seedlessRefitV3`; per event the
  clusters per module in cluster-table order; per OT track `Engine::run`; the sidecar.

## Per track (the replay's functions in brackets)

1. **Visit plan** (`arms.visit_plan`): the surfaces on the OT helix's path (barrel layers always; discs
   met within their radial extent +- 1 cm), outside-in, restricted to the scenario's active surfaces.
   A track is **disc-going** when a disc is on that path (before the mask): it selects the material
   model and the cut bins.
2. **Seed** (`arms.seed`): the OT covariance, the d0 floor, plus the seed model dC(|eta|, p).
3. **Process noise** before each visit (`arms._q_step`): Highland kinks on the material planes between
   the previously *visited* surface and this one: barrel-only tracks by radius (the calibration planes),
   disc-going tracks along their path within each plane's extent, discs at the fitted scale.
4. **Window** (`ot_projection.evaluate_surface`, geometry only): every module the 4-sigma position
   ellipse may overlap, by a distance bound from the nearest module's prediction, then the overlap
   test; every cluster with m^2 < nsig^2 in the prediction's (u, v) covariance. A failed nearest-module
   projection takes every module of the surface (side).
5. **Measurements** (`arms.measurements`): position errors from the R tables; the angle tables'
   centre and width at the predicted incidence key.
6. **Ranking** (`arms.rank_values`): the innovation chi2 over the positions and the usable angles the
   rank consumer uses, +1 per unusable dimension.
7. **Two-layer hit choice** (`arms._tree2_step`, `_tree2_commit`, `tree2_scores`): at each visit a window
   opens with three copies of the track (best applied, runner-up applied, skip); they are committed at
   the next visit by score = first-layer chi2 (skip: P) + the best second-layer chi2 capped at P (an
   empty layer costs P; at most one chosen skip), ties best, runner-up, skip; a window still open at the
   end commits on its first layer alone. The copies keep their own material bookkeeping (the previously
   visited surface) and their own 8-update budget. P is the 99 % gate cut of the window's surface and
   first / later visit step, in the angle dims the rank consumer uses at that step.
8. **Update** (`engine.kf_update`): sequential scalar updates x, y, cotAlpha, cotBeta with
   relinearisation; guarded Jacobian columns zeroed; a scalar with chi2 > chi2UpdateGate is skipped; a
   prediction with |cot| > 12 is not applied.
9. **Inner material** (`arms.inner_material`, rule `path`): after the last visit, C only.

Not ported: drop-and-refit (and the keep-one-cluster rule that acts only inside it), the forward hit
gate, the window and gate angle consumers. The configuration refuses them.

## Arithmetic

v3 crosses in double; v2 used float32 (with fma on aarch64), and the replay's entry gate showed that
float32 differences amplify after large-pull updates. A future FPGA-fidelity path starts from v2's
arithmetic findings, not from v3's.

## Known limitations (carried from the replay; v3cmssw/PLAN.md §6)

The disc material scale 2.67 is physically suspect; the IL1 / IL2 later-visit errors are 9-19 % too
small; 1.6-5.8 % of the pulls lie beyond 3 sigma; the disc u pulls of particles born beyond IL1 are 1.6x
too wide (prompt-fitted seed); AAAA barrel heavy-flavour attachment is -2.3 pt against v2; the beta
9.5-11.3 degree disc bin; the ID1 u errors at 10-20 GeV are too large; the seed model above |eta| 2.25 at
10-15 GeV; the two-layer choice's skip loss; the replay's GO numbers were made with truth-owner modules in
the window (slightly optimistic); the angles are placeholders.
