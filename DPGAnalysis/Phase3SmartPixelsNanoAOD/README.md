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
  `L1Trigger/Phase3SmartPixels/doc/RefitSidecarSpec.md` §4.4). ONLY for digiRefit
  variants (the sidecar exists only in that mode); wired automatically from
  `smartPixelsCoexist`/`smartPixelsCoopt` when a digiRefit variant is present.
- `l1tPh3SmartPixelsNanoFlavours_cff.py`: the autoNANO add-on flavours (below).

## Add-on flavours

Customize-only autoNANO entries, composed after a base tier `@Phase2L1DPG[withGen]` in any
order: `-s NANO:@Phase2L1DPGwithGen+@A+@B+...`. Each add-on only adds tables. Index columns
linking two tables exist only when both tables are in the nano; every add-on re-derives
them (`_reconcile`), so whichever end is added second creates the link.

| flavour | adds | links (when the other end exists) |
|---|---|---|
| `@L1Trk` | `L1TTrack`, `L1TExtTrack` (15 `cov_*`), GTT columns; with gen: track truth | `L1DispVertex` track indices |
| `@L1PF` | `L1PFCand`, `L1PuppiCand`, `L1ExtPuppiCand`, `L1HGCCluster`, SC4/SC8/SC4NG constituent links | candidate `l1TrackIdx` -> `L1TTrack`/`L1TExtTrack`, or the SmartPixels tables when a variant is wired |
| `@TrkStubs` | `L1TOTStub` (every OT stub, with truth) | OT track -> stub |
| `@SPTrk` | `L1TSmartPixelsTrack`, `L1TSmartPixelsExtTrack` (15 `cov_*`) of the wired variant | `l1TrackIdx` -> the OT track it was built from |
| `@SPClustersStubs` | `L1TOTStub`, `L1TSmartPixelsCluster` (every IT cluster) | as above |
| `@SPSeedlessRefit` | mode tag, digiRefit: `L1TSmartPixels[Ext]RefitHit` + `spix*` columns on the SmartPixels tables | hit -> SmartPixels track, hit -> cluster |
| `@SPParallelSeedRefit`, `@SPJointSeedFit` | reserved (raise) | |
| `@SPActiveLayers<mask>` | mask tag (16 masks; character k = ILk and IDk, A/I); `IIII` = passthrough control | |
| `@TPTruth1GeV` | `L1TTP`: charged TrackingParticles with \|p\| > 1 GeV | target of every `tpIdx` |

- One mode tag and one mask tag select the variant; the second of the two wires it through
  the L1 emulation (`smartPixelsCoopt`). A mode tag without a mask tag, a mask other than
  `IIII` without a mode tag, `@SPTrk` without a variant, two mode or mask tags, or two
  TrackingParticle flavours are errors (the completeness checks run when cmsRun builds the
  process description). `@SPActiveLayersIIII` wires the passthrough variant (SmartPixels
  tracks equal to their OT tracks); a mode tag next to it adds no mode tables.
- Settings not in a tag name: `smartPixelsNanoParams(process, ...)` in
  `--customise_commands` (`trackInputMode`, `extendedTracks`, `promptHnpar`, any
  `digiRefitConfig` key; the angle source defaults to `emu`). cmsDriver runs it after the
  flavours, so it re-applies the settings to what the tags wired; the result equals passing
  them first.
- The composition is recorded in `process.smartPixelsNano`.

## digiRefit sidecar tables (spec §4.4)

The digiRefit producer emits a `smartpixels::SmartPixelsRefitSidecar` 1:1 row-synced
with its refit TTTrack collection under the SAME instance label (`Level1TTTracks`).
`L1SmartPixelsRefitTableProducer` turns `(refit tracks + sidecar)` into two tables per
producer instance (prompt and extended cloned separately):

- **per-hit LINK table** `L1TSmartPixelsRefitHit<Suffix>` (extended:
  `L1TSmartPixelsExtRefitHit<Suffix>`): one row per layer-crossing record across all
  tracks. Columns: `trackIdx` (index of the owning track in the variant track table,
  the `L1SC4NGJetCands` link pattern), `layer`, `detId`, `windowMult`, `flags` +
  unpacked bools (`hitAccepted`, `windowTruncated`, `hasAlpha`, `hasBeta`), `projResX`,
  `projResY` (full float), `recoCotAlpha`, `recoCotBeta`, `sigAlpha`, `sigBeta` (float12),
  `pullX`, `pullY`, `pullAlpha`, `pullBeta` (full float), `chi2IncX`, `chi2IncY`,
  `chi2IncAlpha`, `chi2IncBeta` (float12), and TRUTH-ONLY `selHitClass`,
  `tpLocalCotAlpha`, `tpLocalCotBeta`.
- **track EXTENSION table** (`extension=True`, SAME name+length as the variant track
  table `L1TSmartPixelsTrack<Suffix>` / `L1TSmartPixelsExtTrack<Suffix>`): `spixStatus`
  + unpacked bools (`spixRefitPerformed`, `spixSeedCovOK`,
  `spixAnyWindowTruncated`; `spxParametrizedSeed` was REMOVED with the
  `parametrized` seedCovMode, and status bit2 is retired, not reused), `spixNCrossings`, `spixNAcceptedHits`, `spixNKFUpdates`,
  `spixLayerHitMask`, `spixMaxWindowMult`, `spixChi2IncXTot`, `spixChi2IncYTot`,
  `spixChi2IncAlphaTot`, `spixChi2IncBetaTot`, and `spixCompactWord`
  (= `packCompactWord(trackInfo)`, the 16-bit transmitted-subset word). The r-phi and
  r-z totals are `XTot + AlphaTot` and `YTot + BetaTot`.

Sentinel floats pass through as `-999.f` (consumers test `> -900`); passthrough tracks
have `spixRefitPerformed=false`, empty per-hit rows, and zeroed chi2 totals.

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
    else `-1` for a candidate with no resolvable track ref. On rebuildTracksFromStubs PU RelVals (the
    file's PFTrack collection is dropped) this reads `-2` for all candidates.
- `L1PFCandTrackTruthTableProducer` (`DPGAnalysis/Phase2L1TNanoAOD`): the truth-status
  column `trkTruthStatus` is `0` when the candidate's underlying track resolved (the
  `trkGenuine`/`trkLooselyGenuine`/`trkCombinatoric`/`trkUnknown`/`trkTpPt` columns then
  carry the association result), `-1` for an available product with no resolvable track
  ref, and `-2` when the PFTrack / TTTrack / `TTTrackAssociationMap` product is
  unavailable (the truth columns stay at their unknown defaults). This producer now
  self-guards the PFCandidate->PFTrack->TTTrack deref, so the table is scheduled again
  under rebuildTracksFromStubs (previously it was dropped because it hard-dereffed unstored refs).

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
cmsDriver.py stepWF1 -s RAW2DIGI,L1TrackTrigger,L1,L1P2GT,NANO:@L1PFTrkNanowithGen \
  --conditions auto:phase2_realistic_T35 --geometry ExtendedRun4D121 --era Phase2C22I13M9 \
  --procModifiers nano_l1_hlt --eventcontent NANOAODSIM --datatier NANOAODSIM \
  --customise_commands 'from L1Trigger.Phase3SmartPixels.customizeSmartPixels_cff import smartPixelsCoexist; process = smartPixelsCoexist(process, variants=[("passthrough", None), ("correctionlibRegression", "1100")])' \
  --filein <GEN-SIM-DIGI-RAW> --mc -n -1
```

Notes:
- Coexist runs on a BASE flavour (`@L1TrkNano*`/`@L1PFNano*`/`@L1PFTrkNano*`) and adds one
  SUFFIXED table pair per requested variant; it never writes the unsuffixed
  `L1TSmartPixelsTrack`/`L1TSmartPixelsExtTrack`. Those belong to `@SPTrk`, so coexist
  raises when it is composed. For the cluster tables add
  `addPh3L1SmartPixelsClusters` in `--customise_commands`.
- Downstream objects (vertices, PF/Puppi, jets, MET, GT) reflect the STANDARD tracks; in a
  NANO-only job they are re-emulated on demand from them (standard L1 tasks).

## WF2 — co-opt (canonical single-view nano, file-to-file comparison)

One job per variant. The standard L1 emulation (SimL1Emulator + P2GT, on demand) is
rewired by the generic injector to the chosen SmartPixels collection, so vertices,
PF/Puppi, jets and the NG tagger, GMT tk muons, track jets/HT/MET, displaced vertices
and the GT all follow it. The nano carries the pure OT tables `L1TTrack`/`L1TExtTrack`
and the variant as the unsuffixed `L1TSmartPixelsTrack`/`L1TSmartPixelsExtTrack`
(+ `L1TSmartPixels[Ext]RefitHit` for digiRefit):

```bash
cmsDriver.py stepWF2 -s NANO:@Phase2L1DPGwithGen+@L1PF+@L1Trk+@SPTrk+@SPClustersStubs+@TPTruth1GeV+@SPSeedlessRefit+@SPActiveLayersAAAA \
  --conditions auto:phase2_realistic_T35 --geometry ExtendedRun4D121 --era Phase2C22I13M9 \
  --procModifiers nano_l1_hlt --eventcontent NANOAODSIM --datatier NANOAODSIM \
  --customise_commands 'from DPGAnalysis.Phase3SmartPixelsNanoAOD.l1tPh3SmartPixelsNanoFlavours_cff import smartPixelsNanoParams; from DPGAnalysis.Phase3SmartPixelsNanoAOD.l1tPh3SmartPixelsNano_cff import produceMiniAODGenInputs; smartPixelsNanoParams(process, trackInputMode="rebuildTracksFromStubs"); process = produceMiniAODGenInputs(process)' \
  --filein <GEN-SIM-DIGI-RAW> --mc -n -1
```

Notes:
- The mode and mask tags are the single selector; `smartPixelsCoopt(process, variant=...)`
  on a stock flavour (`addNanoTables=True`) writes the same tables (`addCooptNanoTables`).
- The injector logs every rewired InputTag (`SmartPixels injector: <module> ...`); it
  never touches nano tables. Verify with `--no_exec` + `python3 stepWF2_*.py`.
- Table wiring with a variant wired: `L1PuppiCand`/`L1ExtPuppiCand` `l1TrackIdx` index
  `L1TSmartPixelsTrack`/`L1TSmartPixelsExtTrack` (the branch doc names the table);
  candidate track truth uses the associator of the SmartPixels collection; `L1DispVertex`
  track indices are `L1TSmartPixelsExtTrack` rows; `L1TTrack` keeps its own GTT hardware
  columns from a GTT conversion of the OT tracks.
- Every SmartPixels track table has `l1TrackIdx`: the row of the OT track it was seeded from.
- `produceMiniAODGenInputs` only for inputs without MINIAOD gen products (withGen tiers).
- Loop `activeSP` over 0000..1111 (or modes) across jobs; comparisons are file-to-file.
- For the full L1+HLT+Reco no-filter study, extend the steps per
  the recipes in the L1Nano documentation, e.g.
  `-s ...,HLT:NGTScouting,RECO,RECOSIM,PAT,NANO:@NGTScoutingVal+@L1PFTrkNanowithGen+@PHYS
  --procModifiers ngtScouting,nano_l1_hlt`.

## Requirements

The SmartPixels collections must exist: run within the `smartpixels-phase3` branch and,
for correctionlib modes, place the correction JSON under
`L1Trigger/Phase3SmartPixels/data/` (see that package's README).
