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

    # the canonical single-view nano: ONE variant drives the L1 emulation
    ./makeSpixConfig.py --pu 200 --tier clusters-truth --coopt --variant digiRefit:1111 \
        --events 100 -o /work/spix_coopt.py

No nano table is ever dropped: every table of the tier is produced (the L1
emulation runs on demand, the MINIAOD gen inputs of a truth tier are made in-job).

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
#     (base L1 tier, SmartPix payload, withGen). --coopt runs the SmartPix flavor
#     @<base>SmartPix<payload>[withGen]; coexist never writes its canonical unsuffixed
#     track tables, so it runs the base flavor @<base>[withGen] and adds the payload
#     tables with addPh3L1SmartPixels<payload>.
TIERS = {
    "trk":         ("L1TrkNano",   "", False),
    "trk-truth":   ("L1TrkNano",   "", True),
    "pftrk":       ("L1PFTrkNano", "", False),
    "pftrk-truth": ("L1PFTrkNano", "", True),
    "pf":          ("L1PFNano",    "", False),
    "pf-truth":    ("L1PFNano",    "", True),
    # Payload tiers. "clusters" carries the UNTRUNCATED IT cluster table: ~26.5k
    # rows/event at PU200, measured 0.29 MB/event (88.6 stored bits/cluster).
    # "reco" is RESERVED and raises until its content is defined.
    "clusters":       ("L1PFTrkNano", "Clusters", False),
    "clusters-truth": ("L1PFTrkNano", "Clusters", True),
    "reco":           ("L1PFTrkNano", "Reco",     False),
    "reco-truth":     ("L1PFTrkNano", "Reco",     True),
}


def nano_step(args):
    base, payload, truth = TIERS[args.tier]
    gen = "withGen" if truth else ""
    return f"NANO:@{base}SmartPix{payload}{gen}" if args.coopt else f"NANO:@{base}{gen}"


GEOMETRY, ERA, CONDITIONS = "ExtendedRun4D121", "Phase2C22I13M9", "auto:phase2_realistic_T35"

# digiRefitConfig keys this CLI exposes. Everything else takes DIGIREFIT_DEFAULTS.
REFIT_AXES = {
    "applyProcessNoise": bool,   # MS state-covariance inflation
    "multScattTerm": float,      # projected scattering angle per layer = term/pT [rad*GeV]
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

# v6 = v5 (physical units; v4 and earlier were x1.88 (alpha) / x6.53 (beta) too narrow) with a fine beta axis
DEFAULT_ANGLE_SET = "L1Trigger/Phase3SmartPixels/data/spix_angle_response_Conv1D_Full-2bit_v6.json"


def build_customise(args, overrides):
    """The --customise_commands payload: smartPixelsCoexist (or smartPixelsCoopt with
    --coopt), plus the in-job MINIAOD gen inputs for a truth tier."""
    # --variant is repeatable. smartPixelsCoexist takes a LIST of variants and gives
    # each its own producers and its own nano tables, so N configurations cost one
    # pass over the input instead of N: the expensive upstream work (unpacking,
    # track finding) is shared and only the refit is repeated. Needed for the
    # activeSP sweep, where 15 non-trivial masks as 15 separate jobs would be
    # 15 full reconstructions to vary one cheap downstream step.
    specs = args.variant if isinstance(args.variant, list) else [args.variant]
    if args.coopt and len(specs) != 1:
        raise SystemExit("--coopt takes exactly one --variant: it is the collection the whole "
                         "L1 emulation is built from")
    variants = []
    seen = set()
    for spec in specs:
        mode, _, activeSP = spec.partition(":")
        if mode == "digiRefit" and not activeSP:
            raise SystemExit(f"--variant digiRefit requires an activeSP, e.g. digiRefit:1111 (got '{spec}')")
        key = (mode, activeSP)
        if key in seen:
            raise SystemExit(f"--variant {spec} given more than once; each variant makes an "
                             "identically-named set of tables and the duplicate would collide")
        seen.add(key)
        variants.append(f'("{mode}", "{activeSP}")' if activeSP else f'("{mode}", None)')

    refit = {"pixelavAngleSet": args.pixelav_angle_set}
    refit.update(overrides)
    common = (f"trackInputMode='{args.track_input_mode}', extendedTracks={args.extended_tracks}, "
              f"promptHnpar={args.prompt_hnpar}, digiRefitConfig={refit!r}")
    _, payload, truth = TIERS[args.tier]
    parts = [
        "from L1Trigger.Phase3SmartPixels.customizeSmartPixels_cff import smartPixelsCoexist, smartPixelsCoopt",
        "import DPGAnalysis.Phase3SmartPixelsNanoAOD.l1tPh3SmartPixelsNano_cff as p3nano",
    ]
    if args.coopt:
        parts.append(f"process = smartPixelsCoopt(process, variant={variants[0]}, {common})")
    else:
        if payload:
            parts.append(f"process = p3nano.addPh3L1SmartPixels{payload}(process)")
        parts.append(f"process = smartPixelsCoexist(process, variants=[{', '.join(variants)}], "
                     f"addNanoTables=True, {common})")
    # withGen tiers read MINIAOD gen products (prunedGenParticles, slimmedGenJets, ...);
    # the GEN-SIM-DIGI-RAW RelVals carry only genParticles/ak4GenJetsNoNu, so make them.
    if truth:
        parts.append("process = p3nano.produceMiniAODGenInputs(process)")
    return "; ".join(parts)


def cmsdriver_cmd(args, out_py, overrides):
    if args.needs_truth and not TIERS[args.tier][2]:
        raise SystemExit(
            f"--tier {args.tier} does not keep the TrackingParticle truth extension, but "
            "--needs-truth was requested. Use a '-truth' tier; otherwise L1TTrack_genuine "
            "and L1TTrack_tp_* are absent and no resolution or hit-purity study can run.")
    infile = args.filein or INPUTS.get((args.pu, args.release_compat))
    if not infile:
        raise SystemExit(f"no known input for pu={args.pu} release={args.release_compat}; "
                         "pass --filein explicitly")
    return [
        "cmsDriver.py", os.path.splitext(os.path.basename(out_py))[0],
        "-s", nano_step(args),
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
    ap.add_argument("--variant", action="append", default=None,
                    help="mode[:activeSP], e.g. digiRefit:1111, passthrough. REPEATABLE: give it "
                         "several times to emit one table set per variant from a single pass over "
                         "the input, which is how the activeSP sweep avoids 15 reconstructions")
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
    ap.add_argument("--coopt", action="store_true",
                    help="canonical single-view nano (smartPixelsCoopt): the one --variant "
                         "replaces the tracklet tracks for the whole L1 emulation")
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
    if not args.variant:
        args.variant = ["digiRefit:1111"]

    def _tobool(v):
        # bool("False") is True in Python, so --set applyProcessNoise=False would
        # silently pass True and make an A/B degenerate. Parse the WORD.
        if v.lower() in ("false", "0", "no", "off"):
            return False
        if v.lower() in ("true", "1", "yes", "on"):
            return True
        raise SystemExit(f"cannot parse {v!r} as a boolean; use true/false")

    def parse_kv(s):
        k, _, v = s.partition("=")
        if k not in REFIT_AXES:
            raise SystemExit(f"unknown digiRefit axis {k!r}; known: {sorted(REFIT_AXES)}")
        conv = REFIT_AXES[k]
        return k, (_tobool(v) if conv is bool else conv(v))

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
