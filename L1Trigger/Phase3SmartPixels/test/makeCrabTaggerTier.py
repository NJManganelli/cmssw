#!/usr/bin/env python3
"""Generate ready-to-submit CRAB tasks for the SmartPixels tagger pseudo-nano tier (CMSSW_20_1_0_pre1, D121; or, with
--campaign spring24, the old Phase2Spring24 samples in a CMSSW_17 area, which alone can read them):
ONE TASK PER MASK, as canon_v5 (the physicist, 2026-10-22: no merged multi-mask collections), so every file's
L1 downstream (PF / Puppi, jets and their NG scores, e/gamma, muons, GT bits) comes from that mask's tracks.

Per mask, writes (run inside the CMSSW environment, which `crab submit` needs anyway):
  1. the pset, made by cmsDriver.py from NANO:@L1TSPTaggerNanoV3+@SPActiveLayers<mask>, with the canonical study
     posture (rebuildTracksFromStubs, extended tracks, the 5-parameter prompt fit, the emulator angles) and the
     generator inputs made in-job, as canon_v5; optionally the track heads (--heads DIR: the conifer models,
     shipped flat in the sandbox, SPEC.md);
  2. the CRAB config (the July unified-nano generator's conventions, makeCrabUnifiedNano.py), with the pset and,
     with --heads, the model files as JobType.inputFiles.

NOTHING is queried or submitted: the dataset is left blank (TODO) unless --inputDataset is given. A DAS query and
`crab submit` need the physicist's grid proxy and explicit go. The d121 inputs must be post-cms-sw/cmssw#51503
PU200 D121 GEN-SIM-DIGI-RAW (TTbar, H->tautau, ...), as the local test files; the spring24 inputs are the
Phase2Spring24 GEN-SIM-DIGI-RAW-MINIAOD samples (customiseSpring24Inputs: the input fixes, without nano_l1_hlt).

  python3 makeCrabTaggerTier.py [--campaign d121|spring24] [--masks AAAA AAAI IAAA IIII]
                                [--heads /path/to/heads_cmssw/models] [--inputDataset /A/B/C]
"""
import argparse
import os
import shutil
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from makeCrabUnifiedNano import DATASET_UNSET, primary_name, sanitize  # noqa: E402
from makeSpixConfig import CAMPAIGNS, INPUTS  # noqa: E402

POSTURE = ("import DPGAnalysis.Phase3SmartPixelsNanoAOD.l1tPh3SmartPixelsNano_cff as p3nano; "
           "from DPGAnalysis.Phase3SmartPixelsNanoAOD.l1tPh3SmartPixelsNanoFlavours_cff import smartPixelsNanoParams; "
           "smartPixelsNanoParams(process, trackInputMode='rebuildTracksFromStubs', extendedTracks=True, "
           "promptHnpar=5, angleSource='emu')")
# the generator inputs / input fixes per campaign. A CRAB pset is loaded before CRAB sets its input files, so the
# spring24 customisation cannot inspect the input: the dataset choice vouches for its MINIAOD products
GEN = {"d121": "process = p3nano.produceMiniAODGenInputs(process)",
       "spring24": "process = p3nano.customiseSpring24Inputs(process, assumeMiniAOD=True)"}
HEADS = ("from DPGAnalysis.Phase3SmartPixelsNanoAOD.l1tPh3SmartPixelsNanoFlavours_cff import addSmartPixelsTrackHeads; "
         "addSmartPixelsTrackHeads(process, '.', flat=True)")
# a placeholder input: CRAB replaces process.source.fileNames with the dataset's files
PLACEHOLDER = {"d121": ("file:/host_volumes/NJM256GBSD/smartpixels-cmssw-testfiles/"
                        "RelValTTbar_14TeV_PU_150X_mcRun4_realistic_v1_STD_D121_RegeneratedGS_PU-v1_file1.root"),
               "spring24": "file:" + INPUTS[("200", "17_0")]}


def model_files(models, mask):
    """The model files a mask's job reads: the OT heads, and the SmartPixels heads of a refit mask."""
    from DPGAnalysis.Phase3SmartPixelsNanoAOD.l1tPh3SmartPixelsNanoFlavours_cff import (
        TRACK_HEADS, trackHeadModelFile)
    out = []
    for head, _ in TRACK_HEADS:
        out.append(("OT", "IIII", "all", head))
        if mask != "IIII":
            out += [("SmartPixels", mask, w, head) for w in ("prompt", "displaced")]
    return [(trackHeadModelFile(models, *k), trackHeadModelFile(".", *k, flat=True)) for k in out]


def task(a, mask):
    """One mask's pset and CRAB config in <workDir>/<mask>/; returns the CRAB config path."""
    wd = os.path.join(a.workDir, mask)
    os.makedirs(wd, exist_ok=True)
    tag = f"SP{mask}"
    pset = os.path.join(wd, f"taggertier_{tag}_cfg.py")
    custom = POSTURE + "; " + GEN[a.campaign] + ("; " + HEADS if a.heads else "")
    camp = CAMPAIGNS[a.campaign]
    cmd = ["cmsDriver.py", f"taggertier_{tag}_cfg", "-s", f"NANO:@L1TSPTaggerNanoV3+@SPActiveLayers{mask}",
           "--conditions", "auto:phase2_realistic_T35", "--geometry", camp["geometry"], "--era", camp["era"]] + \
        (["--procModifiers", ",".join(camp["modifiers"])] if camp["modifiers"] else []) + \
        ["--datatier", "NANOAODSIM", "--eventcontent", "NANOAODSIM",
         "--filein", PLACEHOLDER[a.campaign], "--fileout", "file:taggertier.root", "--python_filename", pset,
         "-n", "-1", "--mc", "--nThreads", str(a.numCores), "--no_exec", "--customise_commands", custom]
    subprocess.run(cmd, check=True)

    inputs = []
    if a.heads:
        for src, flat in model_files(a.heads, mask):
            dst = os.path.join(wd, os.path.basename(flat))
            shutil.copyfile(src, dst)
            inputs.append(dst)

    primary = primary_name(a.inputDataset)
    request = sanitize(f"{a.requestName or 'taggertier'}_{a.outLFNDirBase.rstrip('/').split('/')[-1]}_{tag}_{primary}")[:100]
    wanted = {"d121": "a post-cms-sw/cmssw#51503 PU200 D121 RelVal (GEN-SIM-DIGI-RAW)",
              "spring24": "a Phase2Spring24 sample (GEN-SIM-DIGI-RAW-MINIAOD)"}[a.campaign]
    lines = ([] if a.inputDataset else
             [f"# TODO: config.Data.inputDataset is EMPTY ({DATASET_UNSET}): fill in {wanted}",
              "#       after a DAS query with the physicist's go.", ""]) + [
        f"# Auto-generated by makeCrabTaggerTier.py -- SmartPixels tagger pseudo-nano tier, one mask per task, "
        f"campaign {a.campaign}",
        f"#   composition: NANO:@L1TSPTaggerNanoV3+@SPActiveLayers{mask}{'  + track heads' if a.heads else ''}",
        "from CRABClient.UserUtilities import config",
        "config = config()",
        "",
        f"config.General.requestName = {request!r}",
        "config.General.workArea = 'crab_projects'",
        "config.General.transferOutputs = True",
        "config.General.transferLogs = True",
        "",
        "config.JobType.pluginName = 'Analysis'",
        f"config.JobType.psetName = {os.path.abspath(pset)!r}",
        f"config.JobType.numCores = {a.numCores}",
        f"config.JobType.maxMemoryMB = {a.maxMemoryMB}",
        "config.JobType.allowUndistributedCMSSW = True",
        f"config.JobType.inputFiles = {[os.path.abspath(f) for f in inputs]!r}",
        "",
        f"config.Data.inputDBS = {a.inputDBS!r}",
        f"config.Data.inputDataset = {a.inputDataset!r}",
        "config.Data.splitting = 'FileBased'",
        f"config.Data.unitsPerJob = {a.unitsPerJob}" + (f"\nconfig.Data.totalUnits = {a.totalUnits}" if a.totalUnits else ""),
        f"config.Data.outLFNDirBase = {a.outLFNDirBase!r}",
        "config.Data.publication = False",
        f"config.Data.outputDatasetTag = {request!r}",
        "",
        f"config.Site.storageSite = {a.storageSite!r}",
        "",
    ]
    crab = os.path.join(wd, f"crab_{request}.py")
    with open(crab, "w") as fh:
        fh.write("\n".join(lines))
    sys.stderr.write(f"{mask}: wrote {pset} and {crab} ({len(inputs)} model files staged)\n")
    return crab


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--campaign", default="d121", choices=sorted(CAMPAIGNS),
                    help="d121 (CMSSW_20 RelVals) or spring24 (the old Phase2Spring24 samples, CMSSW_17 area)")
    ap.add_argument("--masks", nargs="+", default=["AAAA", "AAAI", "IAAA", "IIII"],
                    help="one task per mask (IIII: the passthrough control, OT heads only)")
    ap.add_argument("--heads", default=None, help="the heads' model directory (SPEC.md layout); default: no heads")
    ap.add_argument("--inputDataset", default="", help="DBS dataset; default empty (TODO marker)")
    ap.add_argument("--outLFNDirBase", default="/store/user/lpcsmartpixelscms/taggertier/v1")
    ap.add_argument("--storageSite", default="T3_US_FNALLPC")
    ap.add_argument("--inputDBS", default="global")
    ap.add_argument("--unitsPerJob", type=int, default=1)
    ap.add_argument("--totalUnits", type=int, default=0)
    ap.add_argument("--numCores", type=int, default=4)
    ap.add_argument("--maxMemoryMB", type=int, default=10000)
    ap.add_argument("--requestName", default=None, help="request name prefix (default taggertier)")
    ap.add_argument("--workDir", default="crab_taggertier", help="where the psets, staged models and configs go")
    a = ap.parse_args(argv)
    for mask in a.masks:
        task(a, mask)
    sys.stderr.write("nothing queried or submitted\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
