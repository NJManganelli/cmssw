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

# SPStudyNano on an input WITH MINIAOD products (GEN-SIM-DIGI-RAW-MINIAOD): NANO only. The configuration never
# reads the input, so state its kind (inputHasMiniAOD=True); the TrackingParticle table's generator link needs the
# gen pruner maps rebuilt (see "The generator link on inputs with MINIAOD")
PARAMS_MINIAOD='from DPGAnalysis.Phase3SmartPixelsNanoAOD.l1tPh3SmartPixelsNanoFlavours_cff import smartPixelsNanoParams; smartPixelsNanoParams(process, trackInputMode="rebuildTracksFromStubs", inputHasMiniAOD=True)'
LINK='from DPGAnalysis.Phase3SmartPixelsNanoAOD.l1tPh3SmartPixelsNano_cff import rebuildGenPruningMaps; process = rebuildGenPruningMaps(process)'
cmsDriver.py spx -s NANO:@SPStudyNano+@SPActiveLayersAAAA $COMMON \
  --filein file:<...-MINIAOD> --customise_commands "$PARAMS_MINIAOD; $LINK"

# SPStudyNano on an input WITHOUT MINIAOD (today's RelVals): RECO+PAT in the job, with the
# step-3 options of the RelVal workflow; makeSpixConfig.py --reco builds this line itself
cmsDriver.py spx -s RAW2DIGI,RECO,RECOSIM,PAT,NANO:@SPStudyNano+@SPActiveLayersAAAA $COMMON \
  --filein file:<GEN-SIM-DIGI-RAW> --customise_commands "$PARAMS"
```

`produceMiniAODGenInputs` makes the MINIAOD gen products the withGen tables read; it is only
for inputs without them and jobs without PAT (see "Gen content without MINIAOD" below). `L1Trigger/Phase3SmartPixels/test/makeSpixConfig.py
--tier study --coopt` / `--reco` emit these commands. `@Phase2L1DPG` labels its L1 jet template and L1 PUPPI MET table
`l1jetTable` / `l1puppiMetTable` (as in cms-sw/cmssw#49477), so it composes with `@PHYS`, whose
`jetTable` / `puppiMetTable` would otherwise clash at configuration time.

### Gen content without MINIAOD (`produceMiniAODGenInputs`)

`--procModifiers nano_l1_hlt` (cms-sw/cmssw#51004) lets the L1 nano run without RECO and PAT. In
the flavours used here it changes two things: it removes the stock gen tables that read PAT
products (`puTable`, `metMCTable`, `genJetAK8Table` + flavour), and it removes `hpsTauTable`,
whose input `l1tHPSPFTauProducerPuppi` is no longer in the L1 menu. On a GEN-SIM-DIGI-RAW input
(today's RelVals), `produceMiniAODGenInputs` rebuilds the MINIAOD gen products from the
persisted GEN-SIM products. It uses the stock MINIAOD modules with stock parameters and puts the
removed gen tables back. The reference is the same flavours run without `nano_l1_hlt` on an input
WITH MINIAOD. Their configurations differ only in the rows below; the column tables and all
other modules match.

| nano content | MINIAOD input (reference) | GEN-SIM-DIGI-RAW + `produceMiniAODGenInputs` | match |
|---|---|---|---|
| `GenPart` (`finalGenParticles` <- `prunedGenParticles`) | from the input | stock `prunedGenParticles(WithStatusOne)` on `genParticles` | faithful |
| `GenPart_iso` (`genIso` <- `packedGenParticles`) | from the input | stock `packedGenParticles` | faithful |
| `GenJet` | `slimmedGenJets` from the input | stock `slimmedGenJets` on the persisted `ak4GenJetsNoNu` | faithful |
| `GenJet` flavour | `slimmedGenJetsFlavourInfos` from the input | stock chain: `selectedHadronsAndPartonsForGenJetsFlavourInfos` (on `prunedGenParticles`), `ak4GenJetFlavourInfos`, `slimmedGenJetsFlavourInfos` | faithful |
| `GenJetAK8` + flavour | `slimmedGenJetsAK8` from the input, stock tables | stock `slimmedGenJetsAK8` on the persisted `ak8GenJetsNoNu`, the same stock tables | faithful |
| `GenMET` | stock `metMCTable`: `slimmedMETs.genMET` (a copy of `genMetTrue`) | `genMetTrueTable`: the NANOGEN construction (singleton candidate table on the persisted `genMetTrue`), stock name, columns and precision | same values, different table module |
| `Pileup` | stock `puTable` on `slimmedAddPileupInfo` + `offlineSlimmedPrimaryVertices` | `genPileupTable` = stock `puTable` with `savePUDensityVars=False`, on the stock `slimmedAddPileupInfo` (from `addPileupInfo`) | approximate: `nTrueInt`, `nPU`, `sumEOOT`, `sumLOOT`, `pthatmax` as in the reference; `pudensity`, `gpudensity` ABSENT (counted around the offline primary vertex, which needs RECO) |
| `GenVisTau`, `GenVtx` | same modules on the same inputs | same | identical |
| `L1hpsTau` | `hpsTauTable` reads `l1tHPSPFTauProducerPuppi` | not written (`nano_l1_hlt`) | missing in both on CMSSW_20 inputs: the producer is not in the L1 menu (cms-sw/cmssw#51803 removes the table) |

`SPStudyNano` on a GEN-SIM-DIGI-RAW input runs RECO and PAT in the job, so `@PHYS` writes
`GenMET`, `GenJetAK8` and `Pileup` with the stock modules. Against the MINIAOD-input reference,
only `L1hpsTau` differs. `produceMiniAODGenInputs` refuses that job, because PAT already defines
the products.

`produceMiniAODGenInputs` is only for inputs WITHOUT MINIAOD products: on an input that carries
them, the in-job copies would replace the stored ones. **No configuration reads its input file**
(no configuration-time file access: the same configuration must load on the grid, where the
input is only known to CRAB). The caller states the input kind instead: the generators choose
the customisation from their `--campaign` (d121: `produceMiniAODGenInputs`; spring24:
`customiseSpring24Inputs`), `runUnifiedNano_cfg.py` from its `hasMiniAOD` option.
`produceMiniAODGenInputs` refuses at configuration time when:
- the job defines these products already (PAT in the job);
- the tier has no gen tables;
- `nano_l1_hlt` is off.

On an input WITH MINIAOD under `nano_l1_hlt`, `L1TSPStudyNano` has no `GenMET`, `GenJetAK8` or
`Pileup`, because the modifier removes the stock tables; run such inputs without `nano_l1_hlt`.

### The generator link on inputs with MINIAOD (`rebuildGenPruningMaps`)

The TrackingParticle table's generator link (`L1TTP` `genPartIdx`, `genFromB`, `genFromC`; `@TPTruth1GeVLinked`)
follows the GenParticlePruner maps genParticles -> prunedGenParticlesWithStatusOne -> prunedGenParticles ->
finalGenParticles. MINIAOD stores the pruned collection but not the maps, so on such an input the job stops at the
first event (ProductNotFound, `Association<GenParticle>` `prunedGenParticlesWithStatusOne`).
`rebuildGenPruningMaps` reruns the two stock pruners on the stored `genParticles`, under their stock labels, so
`GenPart` reads the rebuilt collection. This is a documented exception to `produceMiniAODGenInputs`' rule (the
physicist, 2026-10-25), safe while the rebuilt collection equals the stored one. The pruner code is identical in
CMSSW_14_0_9 and CMSSW_17_0_0_pre4, and its configuration differs only under the heavy-ion, ppRef, UPC and oxygen
modifiers. Checked on 10 events of the Phase2Spring24 TTbar: every column bit-identical to the job without the
link, and the link columns filled. Repeat that check on a new sample family. Use it only on inputs with MINIAOD
products: the configuration does not inspect the input.

### Phase2Spring24 inputs (CMSSW_17; `customiseSpring24Inputs`)

The Phase2Spring24 samples (CMSSW_14_0_9 GEN-SIM-DIGI-RAW-MINIAOD; D110 = the T35 tracker of D121) are read only by
releases before CMSSW_20, which stores `io_v1` types: run the backport branch in a CMSSW_17 area. Use their own
geometry and era, the same conditions, and no `nano_l1_hlt`, because the input carries MINIAOD and the stock gen
tables read it. `customiseSpring24Inputs` applies the input fixes:
- the stored `l1tTrackerMuons_l1tTkMuonsGmt*` products are dropped on input, keeping their descendants. Their class
  changed (`TrackerMuon::muRef_` now refers to `l1t::SAMuon`), so reading them fails; the L1 emulation remakes them;
- `rebuildGenPruningMaps` runs when the job has the generator link.

```bash
SPRING24="--conditions auto:phase2_realistic_T35 --geometry ExtendedRun4D110 --era Phase2C17I13M9 \
  --datatier NANOAODSIM --eventcontent NANOAODSIM --mc --nThreads 4"
FIX='from DPGAnalysis.Phase3SmartPixelsNanoAOD.l1tPh3SmartPixelsNano_cff import customiseSpring24Inputs; process = customiseSpring24Inputs(process)'
cmsDriver.py spx -s NANO:@L1TSPStudyNanoV3+@SPActiveLayersAAAA $SPRING24 \
  --filein file:<Phase2Spring24 ...-MINIAOD> --customise_commands "$PARAMS; $FIX"
```

`makeSpixConfig.py --campaign spring24` emits these jobs. `makeCrabTaggerTier.py` infers the campaign from the
`--inputDataset` name: a `Phase2Spring24...` processed dataset with a MINIAOD tier, with no file access and no DAS
query. Compared with the CMSSW_20 tier on D121, the tier is complete: no column is missing or retyped. The extras are
the stock MINIAOD tables `nano_l1_hlt` removes (`Flag`, `Pileup` `pudensity` / `gpudensity`, `L1hpsTau`) and the
input's own HLT and L1 bits.

### Flavours

| flavour | adds | links it creates (when the other end exists) |
|---|---|---|
| `@Phase2L1DPG[withGen]` (base) | L1 objects, GT; gen objects | – |
| `@L1Trk` | `L1TTrack`, `L1TExtTrack` (15 `cov_*`), GTT columns; with gen: track truth | `tpIdx` -> `L1TTP`; `L1DispVertex` track indices |
| `@L1PF` | `L1PFCand`, `L1PuppiCand`, `L1ExtPuppiCand`, `L1HGCCluster`, SC4/SC8/SC4NG constituent links | candidate `l1TrackIdx` -> `L1TTrack`/`L1TExtTrack`, or -> `L1TSmartPixels(Ext)Track` when a variant is wired |
| `@TrkStubs` | `L1TOTStub` (every OT stub, with truth and the `isPS` / `isTilted` module flags) | track -> stub tables `<track table>Stub` (`trackIdx`, `otStubIdx`) for the OT and SmartPixels tracks; stub `tpIdx`, `tpIdxInner`, `tpIdxOuter` |
| `@SPTrk` | `L1TSmartPixelsTrack`, `L1TSmartPixelsExtTrack` (15 `cov_*`) of the wired variant | `l1TrackIdx` -> the OT track it was built from |
| `@SPClustersStubs` | `L1TOTStub`, `L1TSmartPixelsCluster` (every IT cluster, TBPX + ID1-ID4) | track -> stub as above; SmartPixels track -> cluster `<track table>Cluster` (`trackIdx`, `clusterIdx`; seedlessRefit: the accepted IT hits, no rows for passthrough); cluster/stub `tpIdx` |
| `@SPSeedlessRefit` | mode tag: the seedlessRefit variant; its `L1TSmartPixels[Ext]RefitHit` and `spix*` columns on the SmartPixels tables | hit -> SmartPixels track (`trackIdx`), hit -> cluster (`selClusterIdx`) |
| `@SPParallelSeedRefit`, `@SPJointSeedFit` | reserved mode tags (raise `NotImplementedError`) | – |
| `@SPActiveLayers<mask>` | mask tag, 16 masks; character k = ILk and IDk, `A` active, `I` inactive; `IIII` = passthrough control | – |
| `@TPTruth1GeV`, `@TPTruthLinked`, `@TPTruth1GeVLinked`, `@TPTruthAll` | `L1TTP`, see below | target of every `tpIdx` |

### Tables and sizes

`L1TSPStudyNano+@SPActiveLayersAAAA`, 100 PU200 events of RelValTTbar D121 (file3),
`rebuildTracksFromStubs`, emu angles, ZSTD level 5 (the canon_v2 production, build
92bad91b2ba). Logical types as in coffea/awkward (`var` = per-event list); a link table row is
one (parent, child) pair.

| table | logical type | rows/event | kB/event |
|---|---|---:|---:|
| `L1TSmartPixelsCluster` | var * cluster | 72469 | 3984 |
| `L1TTP` (`@TPTruth1GeVLinked`) | var * TrackingParticle | 130527 | 3360 |
| `L1TOTStub` | var * stub | 15250 | 345 |
| `L1TSmartPixelsExtRefitHit` | var * (ext SP track, crossing) | 980 | 96 |
| `L1TSmartPixelsRefitHit` | var * (SP track, crossing) | 770 | 76 |
| `L1TSmartPixelsExtTrack` (+ `spix*`, `l1TrackIdx`) | var * track | 259 | 39 |
| `L1TExtTrack` (+ truth, GTT) | var * track | 259 | 32 |
| `L1TSmartPixelsTrack` (+ `spix*`, `l1TrackIdx`) | var * track | 186 | 27 |
| `L1TTrack` (+ truth, GTT) | var * track | 186 | 24 |
| `L1TExtTrackStub`, `L1TSmartPixelsExtTrackStub` | var * (track, stub) | 1321 each | 3.2 each |
| `L1TTrackStub`, `L1TSmartPixelsTrackStub` | var * (track, stub) | 1008 each | 2.5 each |
| `L1TSmartPixelsExtTrackCluster` | var * (track, cluster) | 896 | 2.4 |
| `L1TSmartPixelsTrackCluster` | var * (track, cluster) | 743 | 2.0 |
| `L1PFCand` | var * candidate | 410 | 9.3 |
| `L1HGCCluster` | var * cluster | 200 | 7.9 |
| `L1PuppiCand`, `L1ExtPuppiCand` | var * candidate | 53, 56 | 1.4, 1.4 |
| `L1SC4JetCands`, `L1SC8JetCands`, `L1SC4NGJetCands` | var * (jet, candidate) | 47, 53, 49 | 0.1 each |
| `L1DispVertex` | var * vertex | 108 | 3.7 |
| `GenPart` | var * particle | 49 | 0.8 |
| 60 other base-tier tables | | | 3.6 |
| **total** | | | **8025** |

Later additions, measured on 3 events of the same file: `GenJetAK8`, `GenMET` and `Pileup`
(gen content without MINIAOD) add 0.6 kB/event, and the `L1TOTStub` `isPS` / `isTilted` flags
0.24 kB/event.

The IIII control writes the same tables (no refit hit tables, empty `...Cluster` link tables),
7824 kB/event on the same events. Two tables dominate: every IT cluster, and the
TrackingParticles linked to them (89% of the `L1TTP` rows have pT < 1 GeV).

### OT stub module flags (`L1TOTStub`)

Two bool columns carry the module flags the L1 track finder attaches to each stub, derived as it
derives them, from the stub's lower-sensor detId (the `detId` column):

- `isPS`: PS module (macro-pixel + strip sensors); false = 2S (two strip sensors).
  `TrackerGeometry::getDetectorType(detId) == Ph2PSP`, the test `trackerDTC::Setup` uses for the
  `psModule` that `L1FPGATrackProducer` hands to the tracklet stub as `isPSmodule`. The tracklet
  chain seeds on disks with PS stubs only, has separate PS and 2S match windows, and keeps a 2S
  disk stub's r only as a ring-half index (ring centre +- 2.5 cm). D121: barrel L1-L3 PS, L4-L6
  2S; disks D1-D2 rings 1-10 and D3-D5 rings 1-7 PS, the five outer rings 2S.
- `isTilted`: barrel module in a tilted ring (`TrackerTopology::barrelTiltTypeP2` = `tiltedZminus`
  or `tiltedZplus`), the tracklet stub's `tiltedBarrel`; the hybrid KF scales the z variance of
  these stubs only. D121: L1-L3 outside their flat central section, 12 rings per side; false for
  flat rings and on disks.

Cost: 0.25 kB/event at PU200 (ZSTD-5), 1 byte/row uncompressed.

### TrackingParticle flavours (`L1TTP`)

Mutually exclusive; with none composed every `tpIdx` is -2.

| flavour | `L1TTP` rows |
|---|---|
| `@TPTruth1GeV` | charged TrackingParticles with 3-momentum \|p\| > 1 GeV |
| `@TPTruthLinked` | TrackingParticles linked by the truth associations to an object the nano writes: OT stubs (stub association, and the owners of both stub clusters), OT tracks and SmartPixels tracks (track associations), IT clusters (dominant charge contributor); neutral ones included if linked. Computed from the association products of the objects present |
| `@TPTruth1GeVLinked` | the union of the two above |
| `@TPTruthAll` | every TrackingParticle of the input (after its own pruning), neutrals included; d0/z0 of a neutral are a straight line |

`tpIdx` (every `tpIdx*` column): >= 0 row in `L1TTP`; -1 not truth-matched; -2 TrackingParticle
product or `L1TTP` unavailable; <= -3 matched TrackingParticle not selected into `L1TTP`,
-3 - tpIdx is its index in the unwritten TrackingParticle collection (cannot occur with
`@TPTruthAll`, nor for the linked objects with the linked selections). Two objects share a
TrackingParticle iff their `tpIdx` are equal. The `L1TTP` table doc states its selection.
`L1TTP_genPartIdx`: >= 0 GenPart row; -1 no generator particle; -2 no kept ancestor; <= -3
nearest kept ancestor is GenPart row -3 - genPartIdx.

### Rules

- **No drops.** A flavour defines the content; no table is ever dropped. A table whose input
  is missing gets its input produced (e.g. `produceMiniAODGenInputs`, the L1 emulation on demand).
  The one exception is the named pseudo-nano tier `@SPTaggerTier` (below), signed off by the
  physicist; no canonical tier drops anything.
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
- **`@PHYS` guard.** With the add-ons composed and no RECO+PAT in the job, `@PHYS` reads MINIAOD products
  from the input. The configuration never reads the input, so building the process fails unless
  `smartPixelsNanoParams(process, inputHasMiniAOD=True)` states it; the message also gives the `-s` string
  with RECO and PAT.
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

## Tagger pseudo-nano tier (`@SPTaggerTier`, 2026-10-21)

A slim, retrainable tier for the jet tagger and its track heads. It is a named pseudo-nano tier: the physicist
signed off on its drops, and the canonical tiers stay complete.

- **One job per mask:** `L1TSPTaggerNanoV3` + a mask tag, as canon_v5 (the physicist, 2026-10-22).
  - Every table of a file belongs to that file's mask. The L1 downstream (PF / Puppi candidates, jets and their NG
    scores, e/gamma, muons, the GT bits, vertices, track jets and sums) is built from that mask's refit tracks,
    exactly as in the canonical nano.
  - Nothing is shared between masks, and no table carries a mask suffix.
  - **Why no merged multi-mask files:** with four variants in production, files that merge several masks'
    collections (one shared downstream, suffixed per-mask tables) confuse more than they save. A shared downstream
    would also be IIII's, built from the OT tracks, and could not give the SmartPixels-minus-IIII results that
    come through the PF candidates and jets.
- **Dropped from the output only:** `L1TSmartPixelsCluster`, `L1TTP`, `L1TOTStub`, `L1TOTStubDTC`.
  - Their producers still run where another table needs them.
  - nanoMetadata `smartPixelsTier` names the tier.
  - Every other table is bitwise identical to the canonical nano's.
- **Truth that travels with the tracks,** in every tier, canonical included (`_configureLinkedTruth`):
  - the SmartPixels track tables copy their OT track's matched TrackingParticle truth (`otTpIdx`, `otTpPdgId`,
    `otGenuine`, `otLooselyGenuine`, `otCombinatoric`, `otUnknown`, `otTpFromHardInteraction`,
    `otTp_pt/eta/phi0/tanL/charge/rInv/d0/z0`), the "otlink" convention;
  - the refit hit tables carry `selTpIdx`, the selected cluster's `tpIdx`;
  - the OT track tables carry `tp_rInv` = tp_charge × c_light/1e5 × Bz(0) / tp_pt.
- **Kept on purpose: the per-hit OT projection columns** (`projOT*`, full float). That was the physicist's decision
  of 2026-10-21, at about 0.36 MB/event (projected from 10 events; projOT* about 0.03 of it) against a
  0.30-0.35 target.
- **Track heads** (optional, `addSmartPixelsTrackHeads(process, modelDir)`):
  - the per-track quality heads (conifer, one model per head and window) on L1TExtTrack and L1TSmartPixelsExtTrack,
    with the job's mask's models;
  - the features are computed from the stored nano values;
  - as built, they need the OT track truth extension (`nStubs`), i.e. a gen tier.
- **Grid:** `L1Trigger/Phase3SmartPixels/test/makeCrabTaggerTier.py` writes one pset and CRAB config per mask. The
  campaign comes from the `--inputDataset` name (`..._STD_D121_...` RelVals or `Phase2Spring24...`, checked against
  the data tier); an unknown name, or no dataset (the TODO marker), needs `--campaign`.
- **A possible future option (not started): one cmsRun, separate files per mask.** The L1 emulation and the per-mask
  refits would share the upstream, with one NanoAOD output module per mask whose per-mask table producers write the
  canonical table names, so each file looks exactly like a per-mask job's. See the v3 port's RESUME.md for the
  effort and the validation it would need.

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
  (= `packCompactWord(trackInfo)`, the 16-bit transmitted-subset word), and TRUTH-ONLY
  `tpIdx`, `tpNOTClusters`, `tpNITClusters` (see "Truth links" below). The r-phi and
  r-z totals are `XTot + AlphaTot` and `YTot + BetaTot`.

Sentinel floats pass through as `-999.f` (consumers test `> -900`); passthrough tracks
have `spixRefitPerformed=false`, empty per-hit rows, and zeroed chi2 totals.

### Truth links of the track and stub tables

The same column name `tpIdx` carries each table's own documented definition:

- `L1TTrack_tpIdx`, `L1TExtTrack_tpIdx` (OT tracks): the TTTrackAssociator match,
  `TTTrackAssociationMap::findTrackingParticlePtr`, non-null iff loosely genuine (one TP contributes
  to every stub except at most one, PS or 2S).
- `L1TSmartPixelsTrack_tpIdx`, `L1TSmartPixelsExtTrack_tpIdx` (extension columns): the
  weighted cluster-majority match, the TP owning a strict majority (> 50%) of the track's cluster
  weight. Each of the two clusters of an OT stub weighs 1/2 (a stub = 1 point) and, for a refitted
  track, each accepted IT cluster weighs 1 (passthrough: the stub clusters only): an IT cluster
  with its angles tags a particle about as tightly as a stub's requirement of two correlated
  clusters. Each cluster has at most one owner, its dominant contributor (below); ownerless
  clusters count in the denominator. The integer counts `tpNOTClusters` (OT stub clusters owned by
  that TP, 0..2 x nStubs) and `tpNITClusters` (accepted IT clusters owned) recompute the
  weighted purity
  `(0.5 * tpNOTClusters + tpNITClusters) / (L1TTrack_nStubs[l1TrackIdx] + spixNAcceptedHits)`,
  and `tpIdx >= 0` iff it exceeds 0.5 (both counts 0 when unmatched). Compare with the OT match
  through `L1TTrack_tpIdx[l1TrackIdx]`; `selHitClass` keeps the OT match as its reference.
- `L1TOTStub_tpIdx`: the stub's principal TP, `TTStubAssociationMap::findTrackingParticlePtr`
  (non-null iff `tpGenuine`). `L1TOTStub_tpIdxInner` / `tpIdxOuter`: the owners of the stub's
  inner- and outer-sensor clusters (`TTStub::clusterRef(0)`, the lower sensor of the stack, and
  `clusterRef(1)`, the upper one).
- `L1TSmartPixelsCluster_tpIdx`: the IT cluster's owner, its dominant charge contributor.

An OT cluster's owner is its dominant contributor in the `TTClusterAssociationMap`
(`TTClusterAssociatorFromPixelDigis:ClusterAccepted`): the TP with the largest share of the
map's own weight, TP pT summed over the cluster's strip x simlink entries (the weight
`isGenuine` cuts on). It equals `findTrackingParticlePtr` for a genuine cluster and is also set
for a combinatoric (merged) one. From the nano, the SmartPixels match is reproduced by the
weighted majority over `L1TOTStub_tpIdxInner/Outer[otStubIdx]` (track -> stub links) and
`L1TSmartPixelsCluster_tpIdx[clusterIdx]` (track -> cluster links).

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

Both workflows keep `--procModifiers nano_l1_hlt` in the recipes. In these flavours it removes
`hpsTauTable` (its input is not in the L1 menu) and the PAT-reading gen tables, which
`produceMiniAODGenInputs` puts back (see "Gen content without MINIAOD"). With `@NGTScoutingVal`
it also removes the HLT-side `dstTriggerAcceptFilter` and the nano-output `SelectEvents`, so ALL
events are stored (the denominator for trigger-rate and efficiency studies). L1-only and L1+Reco
chains are filter-free by construction.

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
- Coexist runs on a BASE flavour (`@L1TrkNano*`/`@L1PFNano*`/`@L1PFTrkNano*`) or on the add-on
  flavours without `@SPTrk`, and adds one SUFFIXED table pair per requested variant; it never
  writes the unsuffixed `L1TSmartPixelsTrack`/`L1TSmartPixelsExtTrack`. Those belong to
  `@SPTrk`, so coexist raises when it is composed. For the cluster tables on a base flavour add
  `addPh3L1SmartPixelsClusters` in `--customise_commands`.
- Each variant also gets the links the add-on flavours give the canonical tables,
  `L1TSmartPixels[Ext]Track<Suffix>Stub` (-> `L1TOTStub`) and `...Cluster`
  (-> `L1TSmartPixelsCluster`), whenever those tables are in the nano, in either order.
- The multi-variant omnibus nano (all 15 non-trivial seedlessRefit masks) is one cmsDriver line:
  `makeSpixConfig.py --tier coexist --variant seedlessRefit:all` composes
  `NANO:@Phase2L1DPGwithGen+@L1Trk+@SPClustersStubs+@TPTruth1GeVLinked` with a single
  `smartPixelsCoexist(...)` call (+ `produceMiniAODGenInputs`).
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
