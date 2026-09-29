# DPGAnalysis/Phase3SmartPixelsNanoAOD

NanoAOD tables and autoNANO flavours for the Phase-3 SmartPixels L1 tracks
(companion to `L1Trigger/Phase3SmartPixels`, which provides the
`L1SmartPixelsTrackProducer` collections and the workflow customisations).

- `l1tPh3SmartPixelsNano_cff.py`: `l1tPh3{,Ext}SmartPixelsTracksTable` (clones of the
  general `l1tTracksTable`) and `addPh3L1SmartPixelsTracks(process, srcLabel=...,
  srcExtendedLabel=..., tableSuffix=...)`, which can target ANY SmartPixels variant
  collection and be called once per variant.
- `addPh3L1SmartPixelsRefitTables(process, ...)` (same module) + the
  `L1SmartPixelsRefitTableProducer` plugin: the seedlessRefit sidecar nano adapter (spec
  `L1Trigger/Phase3SmartPixels/doc/RefitSidecarSpec.md` §4.4). ONLY for seedlessRefit
  variants (the sidecar exists only in that mode); wired automatically from
  `smartPixelsCoexist`/`smartPixelsCoopt` when a seedlessRefit variant is present.
- `l1tPh3SmartPixelsNanoFlavours_cff.py`: the autoNANO add-on flavours, i.e. the nano tiers (below).

## Nano tiers: add-on flavours

The SmartPixels nano is composed on the cmsDriver line from customize-only autoNANO add-ons
after a base tier: `-s NANO:@Phase2L1DPG[withGen]+@A+@B+...` (`l1tPh3SmartPixelsNanoFlavours_cff.py`).

### Standard compositions

| alias | composition | purpose |
|---|---|---|
| `L1TSPStudyNano` | `@Phase2L1DPGwithGen+@L1PF+@L1Trk+@SPTrk+@SPClustersStubs+@TPTruth1GeVLinked+@SPSeedlessRefit`, plus a mask tag | the L1-trigger SmartPixels study tier |
| `SPStudyNano` | `@L1TSPStudyNano+@PHYS`, plus a mask tag | the full study content plus the regular nano (HLT, PUPPI MET, AK4 jets, muons, taus, ...); not L1-only |

```bash
COMMON="--conditions auto:phase2_realistic_T35 --geometry ExtendedRun4D121 --era Phase2C22I13M9 \
  --procModifiers nano_l1_hlt --datatier NANOAODSIM --eventcontent NANOAODSIM --mc --nThreads 4"
PARAMS='from DPGAnalysis.Phase3SmartPixelsNanoAOD.l1tPh3SmartPixelsNanoFlavours_cff import smartPixelsNanoParams; smartPixelsNanoParams(process, trackInputMode="rebuildTracksFromStubs")'
GEN='from DPGAnalysis.Phase3SmartPixelsNanoAOD.l1tPh3SmartPixelsNano_cff import produceMiniAODGenInputs; process = produceMiniAODGenInputs(process)'

# L1TSPStudyNano, seedlessRefit with all four IT layers, on a GEN-SIM-DIGI-RAW RelVal
cmsDriver.py spx -s NANO:@L1TSPStudyNano+@SPActiveLayersAAAA $COMMON \
  --filein file:<GEN-SIM-DIGI-RAW> --customise_commands "$PARAMS; $GEN"

# the same with the passthrough control (SmartPixels tracks = OT tracks, no mode tables)
cmsDriver.py spx -s NANO:@L1TSPStudyNano+@SPActiveLayersIIII $COMMON \
  --filein file:<GEN-SIM-DIGI-RAW> --customise_commands "$PARAMS; $GEN"

# SPStudyNano on an input WITH MINIAOD products (GEN-SIM-DIGI-RAW-MINIAOD): NANO only
cmsDriver.py spx -s NANO:@SPStudyNano+@SPActiveLayersAAAA $COMMON \
  --filein file:<...-MINIAOD> --customise_commands "$PARAMS"

# SPStudyNano on an input WITHOUT MINIAOD (today's RelVals): RECO+PAT in the job, with the
# step-3 options of the RelVal workflow; makeSpixConfig.py --reco builds this line itself
cmsDriver.py spx -s RAW2DIGI,RECO,RECOSIM,PAT,NANO:@SPStudyNano+@SPActiveLayersAAAA $COMMON \
  --filein file:<GEN-SIM-DIGI-RAW> --customise_commands "$PARAMS"
```

`produceMiniAODGenInputs` makes the MINIAOD gen products the withGen tables read; it is only
for inputs without them and jobs without PAT. `L1Trigger/Phase3SmartPixels/test/makeSpixConfig.py
--tier study --coopt` / `--reco` emit these commands. `@Phase2L1DPG` labels its L1 jet template and L1 PUPPI MET table
`l1jetTable` / `l1puppiMetTable` (as in cms-sw/cmssw#49477), so it composes with `@PHYS`, whose
`jetTable` / `puppiMetTable` would otherwise clash at configuration time.

### Flavours

| flavour | adds | links it creates (when the other end exists) |
|---|---|---|
| `@Phase2L1DPG[withGen]` (base) | L1 objects, GT; gen objects | – |
| `@L1Trk` | `L1TTrack`, `L1TExtTrack` (15 `cov_*`), GTT columns; with gen: track truth | `tpIdx` -> `L1TTP`; `L1DispVertex` track indices |
| `@L1PF` | `L1PFCand`, `L1PuppiCand`, `L1ExtPuppiCand`, `L1HGCCluster`, SC4/SC8/SC4NG constituent links | candidate `l1TrackIdx` -> `L1TTrack`/`L1TExtTrack`, or -> `L1TSmartPixels(Ext)Track` when a variant is wired |
| `@TrkStubs` | `L1TOTStub` (every OT stub, with truth) | track -> stub tables `<track table>Stub` (`trackIdx`, `otStubIdx`) for the OT and SmartPixels tracks; stub `tpIdx` |
| `@SPTrk` | `L1TSmartPixelsTrack`, `L1TSmartPixelsExtTrack` (15 `cov_*`) of the wired variant | `l1TrackIdx` -> the OT track it was built from |
| `@SPClustersStubs` | `L1TOTStub`, `L1TSmartPixelsCluster` (every IT cluster, TBPX + ID1-ID4) | track -> stub as above; SmartPixels track -> cluster `<track table>Cluster` (`trackIdx`, `clusterIdx`; seedlessRefit: the accepted IT hits, no rows for passthrough); cluster/stub `tpIdx` |
| `@SPSeedlessRefit` | mode tag: the seedlessRefit variant; its `L1TSmartPixels[Ext]RefitHit` and `spix*` columns on the SmartPixels tables | hit -> SmartPixels track (`trackIdx`), hit -> cluster (`selClusterIdx`) |
| `@SPParallelSeedRefit`, `@SPJointSeedFit` | reserved mode tags (raise `NotImplementedError`) | – |
| `@SPActiveLayers<mask>` | mask tag, 16 masks; character k = ILk and IDk, `A` active, `I` inactive; `IIII` = passthrough control | – |
| `@TPTruth1GeV`, `@TPTruthLinked`, `@TPTruth1GeVLinked`, `@TPTruthAll` | `L1TTP`, see below | target of every `tpIdx` |

### Tables and sizes

`L1TSPStudyNano+@SPActiveLayersAAAA`, 3 PU200 events of RelValTTbar D121 (file3),
`rebuildTracksFromStubs`, emu angles. Logical types as in coffea/awkward
(`var` = per-event list); a link table row is one (parent, child) pair.

| table | logical type | rows/event | kB/event |
|---|---|---:|---:|
| `L1TSmartPixelsCluster` | var * cluster | 73922 | 3123 |
| `L1TTP` (`@TPTruth1GeVLinked`) | var * TrackingParticle | 125110 | 3029 |
| `L1TOTStub` | var * stub | 15343 | 228 |
| `L1TSmartPixelsExtRefitHit` | var * (ext SP track, crossing) | 639 | 62 |
| `L1TSmartPixelsRefitHit` | var * (SP track, crossing) | 485 | 50 |
| `L1TSmartPixelsExtTrack` (+ `spix*`, `l1TrackIdx`) | var * track | 252 | 41 |
| `L1TExtTrack` (+ truth, GTT) | var * track | 252 | 34 |
| `L1TSmartPixelsTrack` (+ `spix*`, `l1TrackIdx`) | var * track | 188 | 31 |
| `L1TTrack` (+ truth, GTT) | var * track | 188 | 28 |
| `L1TExtTrackStub`, `L1TSmartPixelsExtTrackStub` | var * (track, stub) | 1295 each | 2.5 each |
| `L1TTrackStub`, `L1TSmartPixelsTrackStub` | var * (track, stub) | 1019 each | 2.1 each |
| `L1TSmartPixelsExtTrackCluster` | var * (track, cluster) | 581 | 1.4 |
| `L1TSmartPixelsTrackCluster` | var * (track, cluster) | 465 | 1.2 |
| `L1PFCand` | var * candidate | 420 | 10.6 |
| `L1HGCCluster` | var * cluster | 213 | 9.0 |
| `L1PuppiCand`, `L1ExtPuppiCand` | var * candidate | 66, 69 | 3.3, 3.4 |
| `L1SC4JetCands`, `L1SC8JetCands`, `L1SC4NGJetCands` | var * (jet, candidate) | 57, 66, 60 | 0.4 each |
| `L1DispVertex` | var * vertex | 75 | 3.5 |
| `GenPart` | var * particle | 50 | 1.7 |
| 60 other base-tier tables | | | 23 |
| **total** | | | **6693** |

The IIII control writes the same tables (no refit hit tables, empty `...Cluster` link tables),
6552 kB/event. Two tables dominate: every IT cluster, and the TrackingParticles linked to them
(89% of the `L1TTP` rows have pT < 1 GeV).

### TrackingParticle flavours (`L1TTP`)

Mutually exclusive; with none composed every `tpIdx` is -2.

| flavour | `L1TTP` rows |
|---|---|
| `@TPTruth1GeV` | charged TrackingParticles with 3-momentum \|p\| > 1 GeV |
| `@TPTruthLinked` | TrackingParticles linked by the truth associations to an object the nano writes: OT stubs (stub association), OT tracks and SmartPixels tracks (track associations), IT clusters (dominant charge contributor); neutral ones included if linked. Computed from the association products of the objects present |
| `@TPTruth1GeVLinked` | the union of the two above |
| `@TPTruthAll` | every TrackingParticle of the input (after its own pruning), neutrals included; d0/z0 of a neutral are a straight line |

`tpIdx` (and `spixMatchedTpIdx`): >= 0 row in `L1TTP`; -1 not truth-matched; -2 TrackingParticle
product or `L1TTP` unavailable; <= -3 matched TrackingParticle not selected into `L1TTP`,
-3 - tpIdx is its index in the unwritten TrackingParticle collection (cannot occur with
`@TPTruthAll`, nor for the linked objects with the linked selections). Two objects share a
TrackingParticle iff their `tpIdx` are equal. The `L1TTP` table doc states its selection.
`L1TTP_genPartIdx`: >= 0 GenPart row; -1 no generator particle; -2 no kept ancestor; <= -3
nearest kept ancestor is GenPart row -3 - genPartIdx.

### Rules

- **No drops.** A flavour defines the content; no table is ever dropped. A table whose input
  is missing gets its input produced (e.g. `produceMiniAODGenInputs`, the L1 emulation on demand).
- **Order-independent links.** A link column or link table exists only when both of its tables
  are in the nano. Every add-on re-derives all links, the variant wiring and the `L1TTP`
  inputs from what the job contains, so whichever end is added second creates the link and
  the order of the add-ons does not matter (the base tier comes first).
- **Variant selection.** One mode tag and one mask tag; the second of the two wires the variant
  through the L1 emulation (`smartPixelsCoopt`: GTT, correlator, GMT, ... all run on it).
  `@SPActiveLayersIIII` alone or with a mode tag wires the passthrough variant (no mode tables).
  Errors: two mode or mask tags, a mode tag without a mask tag, a mask other than `IIII` without
  a mode tag, `@SPTrk` without a variant, two TrackingParticle flavours, an add-on before the
  base tier. Completeness is checked when cmsRun builds the process description.
- **Settings not in a tag name** go through `smartPixelsNanoParams(process, ...)` in
  `--customise_commands`: `trackInputMode`, `extendedTracks`, `promptHnpar`, and any
  `seedlessRefitConfig` key (`angleSource` defaults to `emu`; `pixelavAngleSet`/`smarthitFakeSet`
  for `angleSource="payload"`). cmsDriver runs it after the flavours, so it re-applies the
  settings to what the tags wired; the result equals passing them first. The composition is
  recorded in `process.smartPixelsNano`.
- **`@PHYS` guard.** With the add-ons composed, the configuration fails unless the input has the
  MINIAOD products `@PHYS` reads or the job runs RECO and PAT; the message gives the `-s` string.
- **OT naming.** The OT track is never called a seed: quantities of the OT track alone carry
  `OT` (`projOT*`, `spixChi2ITAtOT`, `spixOTCovOK`, `C_OT`, `l1TrackIdx` = the OT track a
  SmartPixels track was built from). "Seed" means a real seed (tracklet seed type, IT/OT seeds).
- Coexist (several variants side by side, suffixed tables) stays a `--customise_commands` call
  (`smartPixelsCoexist`), never a canonical tier.
- **Output compression: ZSTD level 5.** Every SmartPixels nano recipe (any add-on flavour,
  `smartPixelsCoexist`, `smartPixelsCoopt`, hence `makeSpixConfig.py` and the CRAB pset
  `runUnifiedNano_cfg.py`) sets it on each `NanoAODOutputModule` through
  `useSmartPixelsNanoCompression` (`l1tPh3SmartPixelsNano_cff`), overriding the cmsDriver
  default LZMA 9. PU200 L1TSPStudyNano: 2.9 vs 3.8 s/event, file 1.18x larger, read ~13x
  faster (LZMA 9 output compression is serial and dominated the event time). Check a file with
  `ROOT.TFile(path).GetCompressionSettings()` (505 = ZSTD 5, 209 = LZMA 9).

## seedlessRefit sidecar tables (spec §4.4)

The seedlessRefit producer emits a `smartpixels::SmartPixelsRefitSidecar` 1:1 row-synced
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
  + unpacked bools (`spixRefitPerformed`, `spixOTCovOK`,
  `spixAnyWindowTruncated`; `spxParametrizedSeed` was REMOVED with the
  `parametrized` otCovMode, and status bit2 is retired, not reused), `spixNCrossings`, `spixNAcceptedHits`, `spixNKFUpdates`,
  `spixLayerHitMask`, `spixMaxWindowMult`, `spixChi2IncXTot`, `spixChi2IncYTot`,
  `spixChi2IncAlphaTot`, `spixChi2IncBetaTot`, and `spixCompactWord`
  (= `packCompactWord(trackInfo)`, the 16-bit transmitted-subset word). The r-phi and
  r-z totals are `XTot + AlphaTot` and `YTot + BetaTot`.

Sentinel floats pass through as `-999.f` (consumers test `> -900`); passthrough tracks
have `spixRefitPerformed=false`, empty per-hit rows, and zeroed chi2 totals.

`spixMatchedTpIdx` is the OT track's truth match and nothing else:
`TTTrackAssociationMap::findTrackingParticlePtr` of the OT track the SmartPixels track was
built from (non-null iff loosely genuine: one TP on every stub except at most one wrong stub,
PS or 2S), so it equals `L1TTrack_tpIdx[l1TrackIdx]` for refitted and passthrough tracks.
IT hits are not used: the one-wrong-stub tolerance is an OT-stub-count rule. IT-hit purity is
`spixNITHitsFromOTTp / spixNAcceptedHits` (accepted IT hits whose dominant charge contributor
is that TP).

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
(+ `L1TSmartPixels[Ext]RefitHit` for seedlessRefit):

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
- Every SmartPixels track table has `l1TrackIdx`: the row of the OT track it was built from.
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
