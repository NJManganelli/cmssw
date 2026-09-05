#!/usr/bin/env python3
"""Emit a SmartPixels cmsRun config from the axes that actually distinguish one.

WHY THIS EXISTS. The development area accumulated 163 hand-kept `*_cfg.py`
files. They collapse to FIFTEEN distinct cmsDriver shapes (see
doc/ConfigProvenance.md): everything else was a per-input-file dump, an
architecture twin, or a one-off rename. Worse, every name encoded the debugging
episode that produced it (`stepCOOPT_nano_fat_1100_coopt_file7`) and none
encoded the settings that decide whether the file is fit for a given
measurement. That cost real time and produced at least two confounded results:
a pull comparison taken from a `parametrized` seed, and a layer-order attempt
staged on a config that silently dropped the TrackingParticle truth columns.

So configs are GENERATED here, from named axes, and are not kept.

    # the standard PU200 refit-development config
    ./makeSpixConfig.py --pu 200 --tier trk-truth --variant digiRefit:1111 \
        --events 100 -o /work/spix_pu200.py

    # noPU smoke, two inner layers only, alpha-only angles
    ./makeSpixConfig.py --pu 0 --tier trk --variant digiRefit:1100 \
        --use-angles alpha --events 3 -o /work/spix_smoke.py

    # A/B on one axis: emit both arms, identical in every other respect
    ./makeSpixConfig.py --pu 200 --tier trk-truth --variant digiRefit:1111 \
        --scan layerOrder=outsideIn,insideOut -o /work/spix_order

`--dry-run` prints the cmsDriver command instead of running it, which is what
doc/ConfigProvenance.md records.

DELIBERATELY NOT AN AXIS: seedCovMode. digiRefit always seeds from the track's
own helixCovMat. The removed `parametrized` mode substituted a fixed diagonal
for a missing covariance and made every seed-covariance-dependent measurement
meaningless while still producing plausible numbers.
"""
from __future__ import annotations

import argparse
import os
import shlex
import subprocess
import sys

from DPGAnalysis.Phase3SmartPixelsNanoAOD.l1tPh3SmartPixelsNano_cff import (
    _ABSENT_MENU_TABLES_DEFAULT,
)

# --- input files, by (pileup, release-compat). See mem:smartpixels-testfile-release-compat.
#     Paths are the WDMac mount; the NJM256GBSD SD card is no longer attached.
TESTFILE_DIR = "/host_volumes/WDMac/smartpixels-cmssw-testfiles"
INPUTS = {
    ("200", "20_1"): f"{TESTFILE_DIR}/RelValTTbar_14TeV_PU_150X_mcRun4_realistic_v1_STD_D121_RegeneratedGS_PU-v1_file1.root",
    ("0", "20_1"): "/work/testfiles/RelValTTbar_D121_noPU_regen.root",
}

# --- nano tiers. "truth" variants keep the TrackingParticle extension
#     (L1TTrack_genuine, L1TTrack_tp_*), WITHOUT WHICH no resolution-vs-truth
#     or hit-purity study can run at all. That omission is exactly what made
#     the first layer-order attempt unusable, so the tier names say it.
TIERS = {
    "trk":        ("NANO:@L1TrkNanoSmartPix",              False),
    "trk-truth":  ("NANO:@L1TrkNanoSmartPixwithGen",       True),
    "pftrk":      ("NANO:@L1PFTrkNanoSmartPix",            False),
    "pftrk-truth": ("NANO:@L1PFTrkNanoSmartPixwithGen",    True),
    "pf":         ("NANO:@L1PFNanoSmartPix",               False),
    "pf-truth":   ("NANO:@L1PFNanoSmartPixwithGen",        True),
    # Payload tiers. "clusters" carries the UNTRUNCATED IT cluster table: ~26.5k
    # rows/event at PU200, measured 0.29 MB/event (88.6 stored bits/cluster).
    # "reco" is RESERVED and raises until its content is defined.
    "clusters":       ("NANO:@L1PFTrkNanoSmartPixClusters",        False),
    "clusters-truth": ("NANO:@L1PFTrkNanoSmartPixClusterswithGen", True),
    "reco":           ("NANO:@L1PFTrkNanoSmartPixReco",            False),
    "reco-truth":     ("NANO:@L1PFTrkNanoSmartPixRecowithGen",     True),
}


GEOMETRY, ERA, CONDITIONS = "ExtendedRun4D121", "Phase2C22I13M9", "auto:phase2_realistic_T35"

# digiRefitConfig keys this CLI exposes. Everything else takes DIGIREFIT_DEFAULTS.
REFIT_AXES = {
    "smarthitFakeSet": str,   # noise-angle inverse CDF; without it "has an angle" is a truth proxy
    "layerOrder": str,
    "useAngles": str,
    "maxHitsPerWindow": int,
    "maxKFUpdates": int,
    "measAngleMaxAbs": float,
    "predAngleMaxAbs": float,
    "clusterMergeFrac": float,
    "bdtModel": str,
}

# Tables keyed to the BASE (non-variant) SmartPixels producers. smartPixelsCoexist
# schedules only the requested VARIANT producers, so these resolve to a missing
# product and abort at the output module. This is a different reason from
# _ABSENT_MENU_TABLES_DEFAULT (which is about the input's reduced L1 menu), hence a
# separate list rather than a change upstream.
BASE_PRODUCER_TABLES = ("l1tPh3SmartPixelsTracksTable", "l1tPh3ExtSmartPixelsTracksTable")

# Tables that are absent on the D121 RelVals, whether from the reduced L1 menu or
# from objects those samples never persisted. This list is EMPIRICAL: it is the
# union arrived at over months of archived configs (see doc/ConfigProvenance.md),
# and it is reproduced here rather than rediscovered because each missing entry
# costs a full multi-minute job to find -- the failure is a ProductNotFound abort
# at the output module on the first event, one table at a time.
# For a PF-carrying tier on a PU sample, stitchPFTierForStubRebuild() un-prunes the
# Puppi/SC4/SC8 family; pass --keep-table to override individual entries here.
RELVAL_ABSENT_TABLES = (
    "gttTracksTable", "gttExtTracksTable",
    "dispVtxTable", "l1tDisplacedVertexTable",
    "l1tPuppiCandsTable", "l1tExtPuppiCandsTable", "l1tPFCandsTable",
    "l1tSC4JetCandsTable", "l1tSC4NGJetCandsTable",
    "l1tHGCClusterTable",
    "l1tPuppiCandHGCClusterLinkTable", "l1tExtPuppiCandHGCClusterLinkTable",
    "l1tPuppiCandTrackTruthTable", "l1tExtPuppiCandTrackTruthTable",
)

PAYLOAD_DIR = "/work/spxsmoke"
DEFAULT_ANGLE_SET = f"{PAYLOAD_DIR}/spix_angle_response_Conv1D_Full-2bit_v4fixed.json"


def build_customise(args, overrides):
    """The --customise_commands payload: one smartPixelsCoexist call plus pruning."""
    mode, _, activeSP = args.variant.partition(":")
    if mode == "digiRefit" and not activeSP:
        raise SystemExit("--variant digiRefit requires an activeSP, e.g. digiRefit:1111")

    refit = {"pixelavAngleSet": args.pixelav_angle_set}
    refit.update(overrides)
    variant = f'("{mode}", "{activeSP}")' if activeSP else f'("{mode}", None)'

    parts = [
        "import FWCore.ParameterSet.Config as cms",
        "from L1Trigger.Phase3SmartPixels.customizeSmartPixels_cff import smartPixelsCoexist",
        "from DPGAnalysis.Phase3SmartPixelsNanoAOD.l1tPh3SmartPixelsNano_cff import "
        "dropAbsentMenuTables, pruneAbsentSimpleTables, useGenParticlesFromFile, dropOrphanExtensionTables",
        f"process = smartPixelsCoexist(process, variants=[{variant}], addNanoTables=True, "
        f"trackInputMode='{args.track_input_mode}', extendedTracks={args.extended_tracks}, "
        f"promptHnpar={args.prompt_hnpar}, digiRefitConfig={refit!r})",
    ]
    # RelVals carry a reduced L1 menu, so some menu tables would resolve to a
    # missing product and abort with ProductNotFound at the output module. The
    # curated default (_ABSENT_MENU_TABLES_DEFAULT) covers the ones seen so far;
    # --drop-tables extends it. Not optional: every archived config that ran on a
    # RelVal needed this, and omitting it is a run-time abort, not a warning.
    # withGen tiers build a MINIAOD-shaped gen chain (finalGenParticles <-
    # prunedGenParticles). A GEN-SIM-DIGI-RAW RelVal has `genParticles` instead, so
    # the pruner must be repointed or genParticleTable aborts with ProductNotFound.
    if TIERS[args.tier][1]:
        parts.append("process = useGenParticlesFromFile(process)")
    drop = [t for t in (list(_ABSENT_MENU_TABLES_DEFAULT) + list(BASE_PRODUCER_TABLES)
                        + list(RELVAL_ABSENT_TABLES) + args.drop_tables)
            if t not in args.keep_table]
    parts.append(f"process = dropAbsentMenuTables(process, {tuple(drop)!r})")
    # AUTO-PRUNE is not optional. A RelVal lacks a long and sample-dependent tail
    # of menu objects, and each missing one is a ProductNotFound abort at the
    # output module on the first event -- discovered one table per multi-minute
    # job. The explicit drop lists above cannot keep up (OMTFpromptMuTable,
    # gttTracksTable, ... were each found that way), so the branch list is taken
    # from the INPUT ITSELF via edmDumpEventContent and everything absent is
    # pruned in one shot.
    if args.labels_file:
        parts.append(
            f"_avail = set(l.strip() for l in open({args.labels_file!r}) if l.strip())")
        parts.append("process = pruneAbsentSimpleTables(process, _avail)")
    # LAST: whatever the two prunes above removed may have orphaned an extension
    # table, which aborts the output module (with a segfault alongside).
    parts.append("process = dropOrphanExtensionTables(process)")
    return "; ".join(parts)


def ensure_labels_file(infile, explicit=None, cache_dir="/work/.spix_labels"):
    """Branch labels present in `infile`, cached. Mirrors the archived extract_labels.py."""
    if explicit:
        return explicit
    os.makedirs(cache_dir, exist_ok=True)
    tag = os.path.basename(infile).replace(".root", "")
    out = os.path.join(cache_dir, f"{tag}.labels.txt")
    if os.path.exists(out) and os.path.getsize(out) > 0:
        return out
    print(f"[labels] edmDumpEventContent {infile} -> {out}", file=sys.stderr)
    r = subprocess.run(["edmDumpEventContent", f"file:{infile}"],
                       capture_output=True, text=True, timeout=1800)
    if r.returncode != 0:
        print(f"[labels] WARNING: edmDumpEventContent failed; auto-prune disabled. "
              f"Expect ProductNotFound aborts.\n{r.stderr[-500:]}", file=sys.stderr)
        return None
    labels = sorted({p.split('"')[1].strip()
                     for p in r.stdout.splitlines() if len(p.split('"')) >= 2
                     if p.split('"')[1].strip()})
    with open(out, "w") as fh:
        fh.write("\n".join(labels) + "\n")
    print(f"[labels] {len(labels)} labels cached", file=sys.stderr)
    return out


def cmsdriver_cmd(args, out_py, overrides):
    steps, keeps_truth = TIERS[args.tier]
    if args.needs_truth and not keeps_truth:
        raise SystemExit(
            f"--tier {args.tier} does not keep the TrackingParticle truth extension, but "
            "--needs-truth was requested. Use a '-truth' tier; otherwise L1TTrack_genuine "
            "and L1TTrack_tp_* are absent and no resolution or hit-purity study can run.")
    infile = args.filein or INPUTS.get((args.pu, args.release_compat))
    if not infile:
        raise SystemExit(f"no known input for pu={args.pu} release={args.release_compat}; "
                         "pass --filein explicitly")
    args.labels_file = ensure_labels_file(infile, args.labels_file)
    return [
        "cmsDriver.py", os.path.splitext(os.path.basename(out_py))[0],
        "-s", steps,
        "--conditions", CONDITIONS,
        "--geometry", GEOMETRY,
        "--era", ERA,
        "--procModifiers", "nano_l1_hlt",
        "--datatier", "NANOAODSIM",
        "--eventcontent", "NANOAODSIM",
        "--filein", f"file:{infile}",
        "--fileout", f"file:{args.fileout or os.path.splitext(out_py)[0] + '.root'}",
        "--python_filename", out_py,
        "-n", str(args.events),
        "--mc", "--nThreads", str(args.threads), "--no_exec",
        "--customise_commands", build_customise(args, overrides),
    ]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pu", choices=["0", "200"], default="200")
    ap.add_argument("--release-compat", default="20_1", choices=["20_1"],
                    help="which release the input must be readable by "
                         "(see mem:smartpixels-testfile-release-compat)")
    ap.add_argument("--tier", choices=sorted(TIERS), default="trk-truth")
    ap.add_argument("--variant", default="digiRefit:1111",
                    help="mode[:activeSP], e.g. digiRefit:1111, passthrough")
    ap.add_argument("--events", type=int, default=100)
    ap.add_argument("--threads", type=int, default=4,
                    help="NOTE multi-threaded runs write events in COMPLETION order, so two "
                         "runs of the same job emit the same events at different positions. "
                         "Paired A/B analyses must sort by (run, lumi, event).")
    ap.add_argument("--track-input-mode", default="rebuildTracksFromStubs",
                    choices=["reemulateL1TrackFinding", "rebuildTracksFromStubs",
                             "useStoredTracks"])
    ap.add_argument("--extended-tracks", default="True", choices=["True", "False"],
                    help="rebuild the displaced chain too. False leaves the Extended "
                         "producer reading old-layout stored tracks, whose all-zero "
                         "helixCovMat makes every extended seed unusable.")
    ap.add_argument("--prompt-hnpar", type=int, default=5, choices=[4, 5])
    ap.add_argument("--use-angles", default=None, choices=["none", "alpha", "alphaBeta"])
    ap.add_argument("--pixelav-angle-set", default=DEFAULT_ANGLE_SET)
    ap.add_argument("--drop-tables", action="append", default=[], metavar="LABEL",
                    help="extra nano table module label to drop (repeatable), on top of "
                         "the curated _ABSENT_MENU_TABLES_DEFAULT")
    ap.add_argument("--keep-table", action="append", default=[], metavar="LABEL",
                    help="do NOT drop this table even though it is in the default "
                         "absent-on-RelVal list (repeatable)")
    ap.add_argument("--labels-file", default=None,
                    help="newline-separated available-branch list for pruneAbsentSimpleTables. "
                         "Derived from the input via edmDumpEventContent and cached if omitted.")
    ap.add_argument("--filein", default=None)
    ap.add_argument("--fileout", default=None)
    ap.add_argument("--needs-truth", action="store_true",
                    help="fail unless the tier keeps L1TTrack_genuine / tp_*")
    ap.add_argument("--set", action="append", default=[], metavar="KEY=VALUE",
                    help=f"digiRefitConfig override; one of {sorted(REFIT_AXES)}")
    ap.add_argument("--scan", default=None, metavar="KEY=V1,V2[,...]",
                    help="emit one config per value, identical otherwise (for A/B)")
    ap.add_argument("-o", "--output", required=True,
                    help="output .py path; with --scan, a prefix (_<value>.py appended)")
    ap.add_argument("--dry-run", action="store_true",
                    help="print the cmsDriver command instead of running it")
    args = ap.parse_args()

    def parse_kv(s):
        k, _, v = s.partition("=")
        if k not in REFIT_AXES:
            raise SystemExit(f"unknown digiRefit axis {k!r}; known: {sorted(REFIT_AXES)}")
        return k, REFIT_AXES[k](v)

    base = dict(parse_kv(s) for s in args.set)
    if args.use_angles:
        base["useAngles"] = args.use_angles

    jobs = []
    if args.scan:
        key, _, vals = args.scan.partition("=")
        if key not in REFIT_AXES:
            raise SystemExit(f"unknown scan axis {key!r}; known: {sorted(REFIT_AXES)}")
        for v in vals.split(","):
            ov = dict(base)
            ov[key] = REFIT_AXES[key](v)
            jobs.append((f"{args.output}_{v}.py", ov))
    else:
        jobs.append((args.output if args.output.endswith(".py") else args.output + ".py", base))

    for out_py, ov in jobs:
        cmd = cmsdriver_cmd(args, out_py, ov)
        if args.dry_run:
            print(" ".join(shlex.quote(c) for c in cmd))
            print()
            continue
        print(f"--> {out_py}", file=sys.stderr)
        subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()
