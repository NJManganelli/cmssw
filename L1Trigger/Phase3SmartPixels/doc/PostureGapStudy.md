# Posture-A/B gap study: rebuilding new-layout L1 tracks from persisted stubs ("posture C")

2026-07-19. Design study + executed proof-of-concept. Artifacts in `/work/spxsmoke/poc_postureC_*`
(config, validation scripts, logs, 10-event output). Nothing committed.

## 1. The gap

| | tracks | covariance (PR #51503) | PU fidelity | truth maps |
|---|---|---|---|---|
| **A** `fromFile` | file (old layout) | **all-zero** → digiRefit forced to `seedCovMode=parametrized` | real PU everywhere | file, consistent |
| **B** `inJob` | fresh (new layout) | real, per-track | **destroyed** — DIGI re-runs on signal-only `g4SimHits`; the PU exists only in the file's HLT-process digis | in-job, consistent |

Brute-force fix (re-run pileup mixing from minbias libraries) is prohibitive for our
resources and connectivity. Question: new-layout tracks with real covariance for BOTH
PV and PU tracks, without re-deriving anything from PU premix inputs?

## 2. Key fact (verified): the PU RelVal persists the entire stub tier WITH pileup

`edmDumpEventContent` of
`RelValTTbar_14TeV_PU_150X_mcRun4_realistic_v1_STD_D121_RegeneratedGS_PU-v1.root`
(full dump: `/work/spxsmoke/poc_postureC_fulldump.txt`), all HLT-process unless noted:

- `TTStubsFromPhase2TrackerDigis` : `StubAccepted` + `StubRejected` (TTStub DetSetVector)
- TTClusters: `TTStubsFromPhase2TrackerDigis:ClusterAccepted/ClusterRejected`,
  `TTClustersFromPhase2TrackerDigis:ClusterInclusive`
- All three truth-map sets: `TTClusterAssociatorFromPixelDigis` (Accepted/Inclusive/Rejected),
  `TTStubAssociatorFromPixelDigis` (Accepted/Rejected),
  `TTTrackAssociatorFromPixelDigis[Extended]:Level1TTTracks`
- `mix:MergedTrackTruth` TrackingParticles + TrackingVertices; `prunedTrackingParticles`
- `offlineBeamSpot` (and `hltOnlineBeamSpot`)
- Even the OT digis themselves: `mix:Tracker` Phase2TrackerDigis **and** the OT simlinks
  `simSiPixelDigis:Tracker` PixelDigiSimLinks (so re-running cluster/stub building from
  file digis is *also* possible — not needed, see §3)
- IT side (digiRefit inputs): `simSiPixelDigis:Pixel` digis + PixelDigiSimLinks; `g4SimHits`
  SimTracks/SimVertices (SIM process)

No gaps for the minimal chain. Everything the stub tier's PU content feeds is in the file.

## 3. Posture C: file stubs → in-job tracklet chain → fresh tracks + fresh truth map

The tracklet track finder consumes **stubs, not digis** (verified in the work-area release,
post-PR#51503):

- `trackerDTC::ProducerDTC` ← `TTStubsFromPhase2TrackerDigis:StubAccepted` → `TTDTC` (not
  persisted; must run in-job)
- `L1FPGATrackProducer` (`l1tTTTracksFromTrackletEmulation`) ← same TTStubs +
  `ProducerDTC:StubAccepted` + `offlineBeamSpot` (file). `readMoreMcTruth=False` → needs no
  truth inputs. ES: `tt::Setup`, `trklet::Setup`, `hph::Setup` (`trackerDTC::SetupRcd`) — the
  posture-A ES-only chain (`TrackTrigger_cff` + `TrackerDTC` setup + `ProducerHPH_cff`) plus the
  setups auto-attached by `DTC_cff` / the tracklet cfi.
- `TTTrackAssociatorFromPixelDigis` re-run in-job pointed at the new tracks, consuming the
  **file's** `ClusterAccepted`/`StubAccepted` truth maps. The new tracks' stub Refs point into
  the file's stub product (same ProductID the file maps are keyed by), so the lookups are
  coherent by construction. Output: fresh `TTTrackAssociationMap` keyed to the new collection.

No DIGI, no mixing, no L1TrackTrigger cluster/stub remaking. Same-label in-process production
shadows the file's HLT branches for all downstream consumers configured without a process name
— identical to the posture-B labeling model, so the digiRefit producer defaults need no changes.
Only the prompt chain is re-run in the PoC; the extended (displaced) chain works identically at
extra CPU if instantiated (its SmartPixels variants are deliberately left unscheduled in the PoC:
their inputs would resolve to old-layout file tracks and trip the trackCov guard, by design).

## 4. PoC results (executed, ARM native, `poc_postureC_cfg.py`)

Job: 10 PU200 events, `ProducerDTC * l1tTTTracksFromTrackletEmulation *
TTTrackAssociatorFromPixelDigis` + SmartPixels `passthrough` + `digiRefit`
(`seedCovMode=trackCov`, `useAngles=alpha`, example PixelAV payload). **rc=0**,
37.7 s wall total (incl. init), 2.9 GB RSS — ~2-3 s/event marginal. Output 315 MB/10 evt
(with TP+stub+TrackingVertex keeps for validation; production keeps would be far slimmer).

**(a) Production**: new track collection produced from file stubs. PASS (rc=0, branches present).

**(b) Fresh layout**: `helixCovMat` non-zero for **1738/1738 (100%)** new tracks;
`sigma(rInv)/rInv` q25/50/75 = **0.204% / 0.264% / 0.364%** — squarely in the 0.2–0.6% band
measured on posture-B in-job tracks. File (HLT) tracks confirmed 0/1633 non-zero (old layout).

**(c) Counts + PU tracks**: 1738 new vs 1633 file tracks (ratio 1.064 — current-algorithm
drift vs the file's older release, expected). Re-run associator: **1440 PU-TP-matched,
273 signal-TP-matched, 25 unmatched** (file map on file tracks: 1359/257/17 — same
composition). PU tracks are fully present and truth-flagged.

**(d) Parameter agreement** (stub-content matching, 1610/1738 pairs, 128 new-unpaired /
23 old-unpaired):

| param | q50 | q90 | max |
|---|---|---|---|
| \|ΔrInv\| | 0 (exact) | 4.2e-6 | 9.9e-4 |
| \|Δphi0\| | 0 (exact) | 9.4e-5 | 1.5e-2 |
| \|ΔtanL\| | 4.9e-4 | 1.5e-3 | 0.27 |
| \|Δz0\| [cm] | 1.5e-2 | 4.4e-2 | 7.1 |
| \|Δd0\| | 0 | 0 | 0 (4-par prompt) |
| \|Δpt\|/pt | 1.1e-7 | 1.0e-3 | 0.26 |

Median rInv/phi agreement is exact; the small tanL/z0 offsets and rare tails are consistent
with PR #51503's fit changes plus release drift — the algorithms are near-identical, as
hypothesized.

**(e) digiRefit smoke with `seedCovMode=trackCov`**: rc=0 and the
`SmartPixelsSeedCovMissing` endStream guard did **not** fire (it would have failed the job).
Engagement 88.2% (1438/1631 tracks pt>2); kick scales physical
(|Δpt|/pt q50/q90 = 0.027/0.149, |Δz0| = 0.12/0.34 cm, |Δphi| = 2.7e-3/2.0e-2).
TP-resolution sanity (matched, pt>2): signal q68 dpt/pt 0.020→0.067, |Δz0| 0.44→0.48 cm;
PU 0.017→0.077, 0.26→0.30 cm. The refit currently *degrades* aggregate resolution — this is
the **known example-payload limit** (structurally valid but physics-empty PixelAV angles,
recorded as Phase-2 known-limit #1), not a posture-C defect; quantitative gates await the real
payload. The point here: the full trackCov-seeded KF path runs end-to-end on real-PU tracks.

**Verdict: posture C WORKS.** It dominates A (real covariance) and B (real PU) simultaneously,
at ~seconds/event with zero re-mixing.

Caveats, honestly stated:
1. New tracks are the *current* algorithm, not the file's (+6.4% count drift). For SmartPixels
   development that is what we want; for validating against the file's own downstream objects
   (e.g. file L1 jets built from file tracks), the drift matters — those studies stay posture A.
2. Anything keyed to the *file track collection* (file track-truth map, file L1 objects built
   on file tracks) does not transfer to the new collection; the fresh map replaces it.
3. FWLite validation gotcha: `TrackingParticle::vz()/vertex()` deref the parentVertex ref —
   keep `TrackingVertexs_mix_MergedTrackTruth_*` alongside TPs in small validation outputs.
4. Extended/displaced chain not exercised (optional re-run, more CPU).

## 5. Alternatives assessed

### 5a. Track-level cross-collection stitching (new in-job signal tracks ⊕ file PU tracks)
Inferior on every axis; rejected with specifics:
- **Fails the primary goal**: the PU-track half still has zero covariance (they remain file
  tracks). Only the signal half gains covariance.
- **Mixed provenance**: posture-B signal tracks come from *re-digitized signal-only* events —
  their stubs are a different product with different content than the file's (no PU in the
  digitization), so the "signal" half is not even the file's signal.
- **Map incoherence**: a merged collection Refs two different stub products;
  `TTTrackAssociator` takes one cluster-map/stub-map pair per instance, so no single coherent
  truth map can cover the merge. Two maps keyed to two sub-collections breaks every downstream
  consumer expecting one map per collection.
- **Duplicates**: signal tracks exist in both halves; dedup requires fuzzy kinematic matching
  (the stub Refs are incomparable across products) — exactly the silent-inconsistency class the
  two-posture rule exists to prevent.

### 5b. Covariance backfill (kinematics-parametrized cov(pt, eta, nstub, chi2, …))
A correctionlib payload fitted on posture-B in-job tracks, applied to file tracks as an
upgraded seed mode (`seedCovMode="parametrizedKinematic"`). Cheap, no re-running, and strictly
better than the current 5-scalar `parametrized` mode. But: no per-track fluctuations (real cov
reflects the actual hit pattern and residuals), and off-diagonal correlations would have to be
parametrized too or dropped. Good enough for: aggregate-resolution studies, BDT feature
generation where the seed cov enters weakly, and any collection that cannot be re-derived.
Not good enough for: chi2-sensitive per-track selections and tails — the very regime the
refit KF cares about. With posture C working at ~2-3 s/event, backfill is a documented
fallback, not a needed path.

### 5c. Timeline dissolution
Once cms-L1TK PR #51503 reaches cms-sw master and a release used for a Phase-2 RelVal
campaign, new PU RelVals carry new-layout tracks with real covariance and posture A becomes
perfect for them (merge/campaign timing not under our control; the current PU RelVal is a
150X_mcRun4 campaign product and will not be regenerated for us). Posture C is the bridge:
it needs nothing from future campaigns and degrades gracefully into posture A the day the
files catch up.

## 6. Pileup-remix contingency (NOT needed — posture C works; kept for the record)
Classic mixing at PU200 needs O(600) minbias GEN-SIM draws per event (O(10k)-event minbias
pool ≥ 5–10 GB download, marginal on a constrained link); premix would need the PU200
premixed neutrino library tier (O(3–4 MB/event) but library access + DataMixer DIGI re-run).
CPU O(1–2 min/event) for DIGI+L1TT at PU200. If ever forced: a **one-time custom-keeps
production** (keep `mix` products + fresh tracks + maps) so downstream studies stay cheap.
All of this is dominated by posture C, which costs none of it.

## 7. Recommendation and productionization
Adopt posture C as the standard PU posture for SmartPixels development. Concretely:
- Add `truthSource="fromFileStubs"` as a third posture in `customizeSmartPixels_cff`:
  attach `ProducerDTC` + `l1tTTTracksFromTrackletEmulation` (+Extended, optional flag) +
  `TTTrackAssociatorFromPixelDigis[Extended]` in-process with default tags (file cluster/stub
  maps + file stubs + file beamspot resolve by shadowing); do NOT schedule DIGI or the
  cluster/stub producers/associators; digiRefit default `seedCovMode=trackCov` becomes valid
  on PU files.
- Loud guard: at config time, this posture must document (and at run time the existing
  `SmartPixelsTruthMissing`/`SmartPixelsSeedCovMissing` guards already enforce) that the file
  carries `TTStubsFromPhase2TrackerDigis:StubAccepted`, the Accepted cluster/stub maps,
  `mix:MergedTrackTruth`, and `offlineBeamSpot`.
- File-requirement summary for any future input: stub tier + maps + TPs + beamspot
  (posture-C core), plus `simSiPixelDigis:Pixel` digis+links and `g4SimHits` SimTracks
  (digiRefit), plus TrackingVertices for FWLite validation convenience.
- PoC artifacts to promote when productionizing: `/work/spxsmoke/poc_postureC_cfg.py`
  (config), `poc_postureC_validate.py` (criteria b–d), `poc_postureC_refitcheck.py`
  (criterion e).
