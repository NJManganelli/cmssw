#!/usr/bin/env python3
"""Generate a ready-to-submit CRAB task for the SmartPixels tagger pseudo-nano tier (CMSSW_20_1_0_pre1, D121).

Writes two files (run inside the CMSSW environment, which `crab submit` needs anyway):
  1. the pset, made by cmsDriver.py from the autoNANO composition
       multi-mask (default):  NANO:@L1TSPTaggerNanoV3Multi   (the L1 emulation once, the seedless refit v3 for
                                                              AAAA / AAAI / IAAA, mask-suffixed tables)
       one mask (--mask M):   NANO:@L1TSPTaggerNanoV3+@SPActiveLayersM
     with the canonical study posture (rebuildTracksFromStubs, extended tracks, the 5-parameter prompt fit, the
     emulator angles) and the generator inputs made in-job, as canon_v5; optionally the track heads (--heads DIR:
     the conifer models, shipped flat in the sandbox, SPEC.md);
  2. the CRAB config (the July unified-nano generator's conventions, makeCrabUnifiedNano.py), with the pset and,
     with --heads, the model files as JobType.inputFiles.

NOTHING is queried or submitted: the dataset is left blank (TODO) unless --inputDataset is given. A DAS query and
`crab submit` need the physicist's grid proxy and explicit go. The RelVal inputs must be post-cms-sw/cmssw#51503
PU200 D121 GEN-SIM-DIGI-RAW (TTbar, H->tautau, ...), as the local test files.

  python3 makeCrabTaggerTier.py [--mask AAAA] [--heads /path/to/heads_cmssw/models] [--inputDataset /A/B/C]
"""
import argparse
import os
import shutil
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from makeCrabUnifiedNano import DATASET_UNSET, primary_name, sanitize  # noqa: E402

POSTURE = ("import DPGAnalysis.Phase3SmartPixelsNanoAOD.l1tPh3SmartPixelsNano_cff as p3nano; "
           "from DPGAnalysis.Phase3SmartPixelsNanoAOD.l1tPh3SmartPixelsNanoFlavours_cff import smartPixelsNanoParams; "
           "smartPixelsNanoParams(process, trackInputMode='rebuildTracksFromStubs', extendedTracks=True, "
           "promptHnpar=5, angleSource='emu'); process = p3nano.produceMiniAODGenInputs(process)")
HEADS = ("from DPGAnalysis.Phase3SmartPixelsNanoAOD.l1tPh3SmartPixelsNanoFlavours_cff import addSmartPixelsTrackHeads; "
         "addSmartPixelsTrackHeads(process, '.', flat=True)")
# a placeholder input: CRAB replaces process.source.fileNames with the dataset's files
PLACEHOLDER = ("file:/host_volumes/NJM256GBSD/smartpixels-cmssw-testfiles/"
               "RelValTTbar_14TeV_PU_150X_mcRun4_realistic_v1_STD_D121_RegeneratedGS_PU-v1_file1.root")


def model_files(models, masks):
    from DPGAnalysis.Phase3SmartPixelsNanoAOD.l1tPh3SmartPixelsNanoFlavours_cff import (
        TRACK_HEADS, trackHeadModelFile)
    out = []
    for head, _ in TRACK_HEADS:
        out.append(("OT", "IIII", "all", head))
        out += [("SmartPixels", m, w, head) for m in masks for w in ("prompt", "displaced")]
    return [(trackHeadModelFile(models, *k), trackHeadModelFile(".", *k, flat=True)) for k in out]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mask", default=None, help="one SmartPixels mask (default: the multi-mask job)")
    ap.add_argument("--heads", default=None, help="the heads' model directory (SPEC.md layout); default: no heads")
    ap.add_argument("--inputDataset", default="", help="DBS dataset; default empty (TODO marker)")
    ap.add_argument("--outLFNDirBase", default="/store/user/lpcsmartpixelscms/taggertier/v1")
    ap.add_argument("--storageSite", default="T3_US_FNALLPC")
    ap.add_argument("--inputDBS", default="global")
    ap.add_argument("--unitsPerJob", type=int, default=1)
    ap.add_argument("--totalUnits", type=int, default=0)
    ap.add_argument("--numCores", type=int, default=4)
    ap.add_argument("--maxMemoryMB", type=int, default=10000)
    ap.add_argument("--requestName", default=None)
    ap.add_argument("--workDir", default="crab_taggertier", help="where the pset, staged models and config go")
    a = ap.parse_args(argv)

    os.makedirs(a.workDir, exist_ok=True)
    tag = f"SP{a.mask}" if a.mask else "multimask"
    step = f"NANO:@L1TSPTaggerNanoV3+@SPActiveLayers{a.mask}" if a.mask else "NANO:@L1TSPTaggerNanoV3Multi"
    pset = os.path.join(a.workDir, f"taggertier_{tag}_cfg.py")
    custom = POSTURE + ("; " + HEADS if a.heads else "")
    cmd = ["cmsDriver.py", f"taggertier_{tag}_cfg", "-s", step, "--conditions", "auto:phase2_realistic_T35",
           "--geometry", "ExtendedRun4D121", "--era", "Phase2C22I13M9", "--procModifiers", "nano_l1_hlt",
           "--datatier", "NANOAODSIM", "--eventcontent", "NANOAODSIM", "--filein", PLACEHOLDER,
           "--fileout", "file:taggertier.root", "--python_filename", pset, "-n", "-1", "--mc",
           "--nThreads", str(a.numCores), "--no_exec", "--customise_commands", custom]
    subprocess.run(cmd, check=True)

    inputs = []
    if a.heads:
        masks = (a.mask,) if a.mask else ("AAAA", "AAAI", "IAAA")
        for src, flat in model_files(a.heads, masks):
            dst = os.path.join(a.workDir, os.path.basename(flat))
            shutil.copyfile(src, dst)
            inputs.append(dst)

    primary = primary_name(a.inputDataset)
    request = sanitize(a.requestName or f"taggertier_{a.outLFNDirBase.rstrip('/').split('/')[-1]}_{tag}_{primary}")[:100]
    lines = ([] if a.inputDataset else
             [f"# TODO: config.Data.inputDataset is EMPTY ({DATASET_UNSET}): fill in a post-cms-sw/cmssw#51503",
              "#       PU200 D121 RelVal (GEN-SIM-DIGI-RAW) after a DAS query with the physicist's go.", ""]) + [
        "# Auto-generated by makeCrabTaggerTier.py -- SmartPixels tagger pseudo-nano tier",
        f"#   composition: {step}{'  + track heads' if a.heads else ''}",
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
    crab = os.path.join(a.workDir, f"crab_{request}.py")
    with open(crab, "w") as fh:
        fh.write("\n".join(lines))
    sys.stderr.write(f"wrote {pset} and {crab} ({len(inputs)} model files staged); nothing submitted\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
