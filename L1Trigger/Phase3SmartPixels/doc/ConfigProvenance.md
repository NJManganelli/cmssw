# Config provenance: the cmsDriver commands behind the development configs

Development ran for months out of `cmssw/work/` and `cmssw/work/spxsmoke/`, which
accumulated **163 hand-kept `*_cfg.py` files**. They are being replaced by a
generator, `test/makeSpxConfig.py`. This document is what has to survive the
deletion: the commands that produced them, and the commit range each command is
still valid at.

## Why they are going away

The 163 files collapse to **15 distinct cmsDriver shapes**. Everything else was a
per-input-file dump (`..._file1` … `..._file10`), an architecture twin (`_x86`,
`_arm`), or a rename. Two structural problems made them actively harmful:

1. **Names encoded the episode, not the configuration.**
   `stepCOOPT_nano_fat_1100_coopt_file7_cfg.py` tells you which debugging session
   produced it. It does not tell you the pileup, the seed covariance mode, whether
   the extended chain was rebuilt, or whether the TrackingParticle truth columns
   survived — which are the four things that decide whether a file can answer a
   given question. Measured across all 163: `seedCovMode`, `trackInputMode`,
   `extendedTracks` and truth-column retention appear in **zero** filenames.

2. **They rot silently.** Every one of the 163 uses the retired `truthSource=`
   keyword with the pre-rename vocabulary (`fromFileStubs` / `fromFile` /
   `inJob`), which `_resolveTrackInputMode` still remaps with a `print`. And
   every config written before 2026-08 fails at import on
   `addPh2L1DisplacedVertices` → `extendPh2L1DisplacedVertices`.

This cost real time and produced two confounded results that had to be withdrawn:
a correct-hit pull comparison taken from a `parametrized` seed, and a first
layer-order attempt staged on a config whose `pruneAbsentSimpleTables` list
silently dropped `L1TTrack_genuine` / `L1TTrack_tp_*`.

## Compatibility watermarks

Use these to decide whether an archived command can be replayed. Commits are on
`smartpixels-phase3-tier2refit`.

| watermark | commit | what changes across it |
|---|---|---|
| `parametrized` seedCovMode **removed** | first commit after `ca9dc92` | Any command whose `digiRefitConfig` sets `"seedCovMode"` **fails at config time** from this point on. Replay such commands at `ca9dc92` or earlier. Nothing is lost: the mode only ever substituted a fixed diagonal for a real covariance. |
| `digiRefitLayerOrder` added, default `outsideIn` | same commit | Before it, the layer loop was hardcoded **inside-out**. Every refit number measured at or before `ca9dc92` is an `insideOut` number, and `insideOut` is measurably worse (§ below). |
| hit candidates become `SiPixelRecHit`s | `ca9dc92` | Before: one candidate per fired `PixelDigi`, position a pixel centre, error `pitch/sqrt(12)`. All combinatorics, pull and window numbers measured earlier are **digi-era** and not comparable. |
| sidecar vocabulary `reco*`/`proj*`/`projRes*`/`truth*` | `1c265e9` | `resX`/`resY` became `projResX`/`projResY`; several nano columns renamed. Analysis scripts targeting v2.6 names reject earlier files outright, by design. |
| 4-way chi2 split | `c521303` | Before: `chi2RPhi`/`chi2RZ` only. `chi2Inc{X,Y,Alpha,Beta}` do not exist earlier. |
| `trackInputMode` vocabulary | `9b8c48f`, `0d122f6`, `b890ce3` | `inJob`→`reemulateL1TrackFinding`, `fromFileStubs`→`rebuildTracksFromStubs`, `fromFile`→`useStoredTracks`. Old spellings still honored with a warning. |
| `extendPh2L1DisplacedVertices` rename | upstream, ~2026-08 | Any config generated before this **cannot be imported** under 20_1. Patch with `sed -e 's/addPh2L1DisplacedVertices/extendPh2L1DisplacedVertices/g'`. |
| input-file staging | n/a | The `NJM256GBSD` SD card is no longer mounted; D121 RelVals live under `/host_volumes/WDMac/smartpixels-cmssw-testfiles/`. Archived commands hardcode the old path. |

## The 15 shapes

Common to all: `--conditions auto:phase2_realistic_T35 --mc --no_exec`, plus
`--customise_commands` carrying a single `smartPixelsCoexist(...)` call. `D121` /
`Phase2C22I13M9` is the 20_1 geometry/era; `D110` / `Phase2C17I13M9` is the older
17_0-era pair and appears only in the pre-August `spx_*` files.

| n | `-s` steps | geometry/era | example file | notes |
|---|---|---|---|---|
| 64 | `NANO:@L1PFTrkNanoSmartPixwithGen` | D121 / C22 | `stepBASE_nano_fat_0000_baseline_file*` | the "fat" PF+track+truth tier; the 10-file fan-outs live here |
| 24 | `DIGI:pdigi_valid,L1TrackTrigger` | D121 / C22 | `dr_clampoff_cfg.py`, `stepDR_*` | digi-era refit development; **all `insideOut`, most `parametrized`** |
| 21 | `NANO:@L1TrkNanoSmartPixwithGen` | D121 / C22 | `stepCLAMP_clamp_on_f1_cfg.py` | track+truth tier; the grazing-clamp and wrong-hit studies |
| 7 | `NANO:@L1TrkNanoSmartPix` | D121 / C22 | `stepAB_abnew_run1_cfg.py` | **no truth columns** — cannot support resolution studies |
| 4 | `DIGI:pdigi_valid,L1TrackTrigger,L1,L1P2GT,NANO:@L1PFTrkNano` | D121 / C22 | `stepNANO_nano_stage1_cfg.py` | single-job DIGI→nano; `_arm`/`_x86` twins |
| 4 | `NANO:@L1PFTrkNanoSmartPix` | D121 / C22 | `stepS2D_nano_puD_PFTrkSmartPix_cfg.py` | |
| 3 | `NANO:@L1PFNanoSmartPix` | D121 / C22 | `stepS2D_nano_puD_PFSmartPix_cfg.py` | |
| 3 | `NANO:@L1PFNanoSmartPixwithGen` | D121 / C22 | `stepS2D_nano_puD_PFSmartPix_withGen_cfg.py` | |
| 3 | `DIGI:...,NANO:@L1PFTrkNano` (`NANOAOD` tier) | D121 / C22 | `stepWF1_cfg.py` | WF1 coexist proof |
| 3 | `L1TrackTrigger,L1,L1P2GT,NANO:@L1PFTrkNanowithGen` | D110 / C17 | `spx_fatflavor_cfg.py` | pre-August, 17_0-era |
| 2 | `L1TrackTrigger,L1,L1P2GT,NANO:@Phase2L1DPGwithGen` + 5 `--customise` | D110 / C17 | `spx_fat_nano_cfg.py` | uses `addPh2L1DisplacedVertices` (broken under 20_1) |
| 1 | same, no `--customise` | D110 / C17 | `spx_baseline_nano_cfg.py` | |
| 1 | `L1TrackTrigger,L1` | D110 / C17 | `spx_reemul17_cfg.py` | re-emulation only, no nano |
| 1 | same as row 11, `--nThreads 4` | D110 / C17 | `spx_valnano_newTT_cfg.py` | |
| 1 | same, with `extendPh2L1DisplacedVertices` | D110 / C17 | `v2p6_base_cfg.py` | the post-rename fix of row 11 |

Full per-file commands remain recoverable from each config's own header, which
cmsDriver writes as `# with command line options: ...`, until the files are
deleted. **Regenerate this table's detail before deleting them** if any number
still cited in `RefitStudyProgram.md` traces to a shape not listed above.

## Replacement

```bash
# the standard PU200 refit-development config (replaces the 21-file CLAMP family)
test/makeSpxConfig.py --pu 200 --tier trk-truth --variant digiRefit:1111 \
    --events 100 --needs-truth -o /work/spx_pu200.py

# an A/B on one axis: both arms identical except the knob
test/makeSpxConfig.py --pu 200 --tier trk-truth --variant digiRefit:1111 \
    --scan layerOrder=outsideIn,insideOut --needs-truth -o /work/spx_order
```

`--dry-run` prints the cmsDriver command, so the generator can also be used to
reconstruct an archived shape rather than replay a stale file.

Two guards exist because their absence caused the two withdrawn results:
`--needs-truth` refuses a tier that would drop `L1TTrack_genuine`/`tp_*`, and
`seedCovMode` is not an axis at all.

## Measured consequence of the `insideOut` watermark

PU200 D121, 100 events, `trackCov`, paired on event ID
(`eval_refitq/ordering/layer_order_ab.py`). This is why pre-`ca9dc92` refit
numbers should not be quoted forward:

| metric | `outsideIn` | `insideOut` (old default) |
|---|---|---|
| wrong-hit fraction L1/L2/L3/L4 | 0.225 / 0.136 / 0.077 / 0.052 | 0.482 / 0.205 / 0.131 / 0.094 |
| refits with zero wrong hits | **75.3%** | 53.9% |
| paired d0 error | **34.1 µm** | 77.1 µm |
| paired z0 error | **62.1 µm** | 182.8 µm |
