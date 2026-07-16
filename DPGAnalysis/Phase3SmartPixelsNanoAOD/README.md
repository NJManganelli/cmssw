# DPGAnalysis/Phase3SmartPixelsNanoAOD

NanoAOD tables and autoNANO flavours for the Phase-3 SmartPixels L1 tracks
(companion to `L1Trigger/Phase3SmartPixels`, which provides the
`L1SmartPixelsTrackProducer` collections and the workflow customisations).

- `l1tPh3SmartPixelsNano_cff.py`: `l1tPh3{,Ext}SmartPixelsTracksTable` (clones of the
  general `l1tTracksTable`) and `addPh3L1SmartPixelsTracks(process, srcLabel=...,
  srcExtendedLabel=..., tableSuffix=...)`, which can target ANY SmartPixels variant
  collection and be called once per variant.
- autoNANO flavours (in `PhysicsTools/NanoAOD/autoNANO.py`): `L1PFTrkNanoSmartPix`,
  `L1PFTrkNanoSmartPixwithGen` = `@L1PFTrkNano{,withGen}` + the default SmartPixels
  track tables.

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
