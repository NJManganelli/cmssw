# DPGAnalysis/Phase3SmartPixelsNanoAOD

NanoAOD tables and autoNANO flavours for the Phase-3 SmartPixels L1 tracks
(companion to `L1Trigger/Phase3SmartPixels`, which provides the
`L1SmartPixelsTrackProducer` collections and the workflow customisations).

- `l1tPh3SmartPixelsNano_cff.py`: `l1tPh3{,Ext}SmartPixelsTracksTable` (clones of the
  general `l1tTracksTable`) and `addPh3L1SmartPixelsTracks(process, srcLabel=...,
  srcExtendedLabel=..., tableSuffix=...)`, which can target ANY SmartPixels variant
  collection and be called once per variant.
- `addPh3L1SmartPixelsRefitTables(process, ...)` (same module) + the
  `L1SmartPixelsRefitTableProducer` plugin: the digiRefit sidecar nano adapter (spec
  `L1Trigger/Phase3SmartPixels/doc/RefitSidecarSpec.md` v0.1 §4.4). ONLY for digiRefit
  variants (the sidecar exists only in that mode); wired automatically from
  `smartPixelsCoexist`/`smartPixelsCoopt` when a digiRefit variant is present.
- autoNANO flavours (in `PhysicsTools/NanoAOD/autoNANO.py`): one per base L1 DPG tier
  x `{,withGen}` -- `L1TrkNanoSmartPix{,withGen}`, `L1PFNanoSmartPix{,withGen}`,
  `L1PFTrkNanoSmartPix{,withGen}` = `@L1TrkNano`/`@L1PFNano`/`@L1PFTrkNano`
  (`{,withGen}`) + the default SmartPixels track tables. The base tiers live in
  `DPGAnalysis/Phase2L1TNanoAOD` (not edited here).

## digiRefit sidecar tables (spec §4.4)

The digiRefit producer emits a `smartpixels::SmartPixelsRefitSidecar` 1:1 row-synced
with its refit TTTrack collection under the SAME instance label (`Level1TTTracks`).
`L1SmartPixelsRefitTableProducer` turns `(refit tracks + sidecar)` into two tables per
producer instance (prompt and extended cloned separately):

- **per-hit LINK table** `L1TSmartPixelsRefitHit<Suffix>` (extended:
  `L1TSmartPixelsExtRefitHit<Suffix>`): one row per layer-crossing record across all
  tracks. Columns: `trackIdx` (index of the owning track in the variant track table,
  the `L1SC4NGJetCands` link pattern), `layer`, `detId`, `windowMult`, `flags` +
  unpacked bools (`hitAccepted`, `windowTruncated`, `hasAlpha`, `hasBeta`), `resX`,
  `resY` (full float), `cotAlphaMeas`, `cotBetaMeas`, `sigAlpha`, `sigBeta` (float12),
  `pullX`, `pullY`, `pullAlpha`, `pullBeta` (full float), `chi2IncRPhi`, `chi2IncRZ`
  (float12), and TRUTH-ONLY `selHitClass`, `parCotAlpha`, `parCotBeta`.
- **track EXTENSION table** (`extension=True`, SAME name+length as the variant track
  table `L1TSmartPixelsTrack<Suffix>` / `L1TSmartPixelsExtTrack<Suffix>`): `spxStatus`
  + unpacked bools (`spxRefitPerformed`, `spxSeedCovOK`, `spxParametrizedSeed`,
  `spxAnyWindowTruncated`), `spxNCrossings`, `spxNAcceptedHits`, `spxNKFUpdates`,
  `spxLayerHitMask`, `spxMaxWindowMult`, `spxChi2IncRPhiTot`, `spxChi2IncRZTot`, and
  `spxCompactWord` (= `packCompactWord(trackInfo)`, the 16-bit transmitted-subset word).

Sentinel floats pass through as `-999.f` (consumers test `> -900`); passthrough tracks
have `spxRefitPerformed=false`, empty per-hit rows, and sentinel chi2 totals.

## Truth index-reuse (withGen tiers)

There is NO refit-specific truth table, by design. The 1:1 output-sync invariant
(spec §1) guarantees that row `i` of every variant's track table (and its extension
columns) derives from input reference track `i`. The `withGen` tiers' reference-track
truth table (`L1TTrackTruth` from `addPh2L1TrackTruth`) therefore aligns 1:1 with every
variant's rows: analysis reads the truth label for a variant track by reusing that same
row index against the reference truth table. Likewise, per-hit `trackIdx` -> variant
track row -> reference truth row. No truth is ever re-keyed to a refit collection.

## Candidate->track linking sentinels: -1 vs -2 (ratified convention)

Candidate->track cross-reference columns in the L1 nano tables distinguish two failure
modes with two sentinels. Do not conflate them:

- **`-1` = ELEMENT-level no-match.** The referenced product IS available and matching
  ran successfully for the collection, but this specific candidate genuinely has no
  track partner (e.g. a neutral candidate with no `pfTrack`, or a constituent `Ptr`
  that does not resolve into an otherwise-present candidate collection).
- **`-2` = PRODUCT-level failure.** The whole referenced collection is not available in
  this event (`isAvailable()` is false / a bare deref would be `ProductNotFound`), so
  linking could not run for ANY candidate. Every row of that column is `-2`.

Sites (all guarded, cf. the tolerance-family fix):

- `L1JetCandLinkTableProducer` (`DPGAnalysis/Phase2L1TNanoAOD`):
  - link table `candIdx` (jet-constituent -> `L1PuppiCand`): `-2` for all rows if the
    `L1PuppiCand` collection is unavailable; else `-1` for a constituent `Ptr` that does
    not resolve.
  - candidate extension `l1TrackIdx` (`L1PuppiCand` -> `L1TTrack`): `-2` for all
    candidates if the PFTrack (or TTTrack) product is unavailable — detected by probing
    `pfTrack().isAvailable()`, since a dropped-but-referenced PFTrack collection dangles;
    else `-1` for a candidate with no resolvable track ref. On posture-C PU RelVals (the
    file's PFTrack collection is dropped) this reads `-2` for all candidates.
- `L1PFCandTrackTruthTableProducer` (`DPGAnalysis/Phase2L1TNanoAOD`): the truth-status
  column `trkTruthStatus` is `0` when the candidate's underlying track resolved (the
  `trkGenuine`/`trkLooselyGenuine`/`trkCombinatoric`/`trkUnknown`/`trkTpPt` columns then
  carry the association result), `-1` for an available product with no resolvable track
  ref, and `-2` when the PFTrack / TTTrack / `TTTrackAssociationMap` product is
  unavailable (the truth columns stay at their unknown defaults). This producer now
  self-guards the PFCandidate->PFTrack->TTTrack deref, so the table is scheduled again
  under posture C (previously it was dropped because it hard-dereffed unstored refs).

Not a linking sentinel: the per-hit refit `trackIdx` in `L1SmartPixelsRefitTableProducer`
is a within-event index into the same-producer, length-checked, 1:1-row-synced refit
track collection (spec §1; a broken sync throws `SmartPixelsSyncBroken`). It has no
product-availability or no-match failure mode and therefore uses NO -1/-2 sentinel.

Both workflows keep `--procModifiers nano_l1_hlt` in the recipes: it removes the
HLT-side `dstTriggerAcceptFilter` and the nano-output `SelectEvents`, so ALL events
are stored (the denominator for trigger-rate/efficiency studies). L1-only and
L1+Reco chains are filter-free by construction.

## WF1 — coexist (in-file track-to-track comparison)

One file with the standard L1 objects AND SmartPixels track tables per variant:

```bash
cmsDriver.py stepWF1 -s RAW2DIGI,L1TrackTrigger,L1,L1P2GT,NANO:@L1PFTrkNanoSmartPixwithGen \
  --conditions auto:phase2_realistic_T35 --geometry ExtendedRun4D121 --era Phase2C22I13M9 \
  --procModifiers nano_l1_hlt --eventcontent NANOAODSIM --datatier NANOAODSIM \
  --customise_commands 'from L1Trigger.Phase3SmartPixels.customizeSmartPixels_cff import smartPixelsCoexist; process = smartPixelsCoexist(process, variants=[("passthrough", None), ("correctionlibRegression", "1100")])' \
  --filein <GEN-SIM-DIGI-RAW> --mc -n -1
```

Notes:
- The `@L1PFTrkNanoSmartPix*` flavour adds the DEFAULT-label tables; `smartPixelsCoexist`
  adds one table pair per requested variant (suffixed). For variant-only tables, use the
  plain `@L1PFTrkNano*` flavour + the coexist customise.
- Downstream objects (vertices, PF/Puppi, jets, MET, GT) still reflect the STANDARD
  tracks in this workflow.

## WF2 — co-opt ("pure view", file-to-file comparison)

One job per variant; every downstream consumer of the standard tracklet tracks reads
the chosen SmartPixels collection instead. Any nano flavour then reflects that single
interpretation:

```bash
# baseline (standard tracks): same command WITHOUT the customise_commands line
cmsDriver.py stepWF2 -s RAW2DIGI,L1TrackTrigger,L1,L1P2GT,NANO:@L1PFTrkNanowithGen \
  --conditions auto:phase2_realistic_T35 --geometry ExtendedRun4D121 --era Phase2C22I13M9 \
  --procModifiers nano_l1_hlt --eventcontent NANOAODSIM --datatier NANOAODSIM \
  --customise_commands 'from L1Trigger.Phase3SmartPixels.customizeSmartPixels_cff import smartPixelsCoopt; process = smartPixelsCoopt(process, mode="correctionlibRegression", activeSP="1100")' \
  --filein <GEN-SIM-DIGI-RAW> --mc -n -1
```

Notes:
- The injector log prints every `l1tTTTracksFromTrackletEmulation --> <module>`
  replacement (GTT input, vertex finder, correlator layers, ...); verify with
  `--no_exec` + `python3 stepWF2_*.py` before large runs.
- The L1TTrack nano table itself is protected by the injector skip-list, so it keeps
  the ORIGINAL tracks as an in-file reference; pass `addPh3Table=True` to also write
  the injected variant as an explicit SmartPixels table.
- Loop `activeSP` over 0000..1111 (or modes) across jobs; comparisons are file-to-file.
- For the full L1+HLT+Reco no-filter study, extend the steps per
  the recipes in the L1Nano documentation, e.g.
  `-s ...,HLT:NGTScouting,RECO,RECOSIM,PAT,NANO:@NGTScoutingVal+@L1PFTrkNanowithGen+@PHYS
  --procModifiers ngtScouting,nano_l1_hlt`.

## Requirements

The SmartPixels collections must exist: run within the `smartpixels-phase3` branch and,
for correctionlib modes, place the correction JSON under
`L1Trigger/Phase3SmartPixels/data/` (see that package's README).
