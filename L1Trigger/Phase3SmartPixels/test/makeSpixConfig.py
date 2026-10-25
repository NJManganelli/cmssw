#!/usr/bin/env python3
"""Emit a SmartPixels cmsRun config from the axes that actually distinguish one.

WHY THIS EXISTS. The development area accumulated 163 hand-kept `*_cfg.py`
files. They collapse to FIFTEEN distinct cmsDriver shapes (see
doc/ConfigProvenance.md): everything else was a per-input-file dump, an
architecture twin, or a one-off rename. Worse, every name encoded the debugging
episode that produced it (`stepCOOPT_nano_fat_1100_coopt_file7`) and none
encoded the settings that decide whether the file is fit for a given
measurement. That cost real time and produced at least two confounded results:
a pull comparison taken from a `parametrized` OT covariance, and a layer-order attempt
staged on a config that silently dropped the TrackingParticle truth columns.

So configs are GENERATED here, from named axes, and are not kept.

    # the standard PU200 refit-development config
    ./makeSpixConfig.py --pu 200 --tier trk-truth --variant seedlessRefit:1111 \
        --events 100 -o spix_pu200.py

    # noPU smoke, two inner layers only, alpha-only angles
    ./makeSpixConfig.py --pu 0 --tier trk --variant seedlessRefit:1100 \
        --use-angles alpha --events 3 -o spix_smoke.py

    # A/B on one axis: emit both arms, identical in every other respect
    ./makeSpixConfig.py --pu 200 --tier trk-truth --variant seedlessRefit:1111 \
        --scan layerOrder=outsideIn,insideOut -o spix_order

    # the canonical single-view nano (add-on flavours): ONE variant drives the L1 emulation
    ./makeSpixConfig.py --pu 200 --tier study --coopt --variant seedlessRefit:1111 \
        --events 100 -o spix_coopt.py

    # the multi-variant omnibus nano: all 15 non-trivial seedlessRefit masks side by side
    ./makeSpixConfig.py --pu 200 --tier coexist --variant seedlessRefit:all \
        --events 100 -o coexist15.py

    # the canonical nano on the old Phase2Spring24 samples, in a CMSSW_17 area (CMSSW_20 cannot read them)
    ./makeSpixConfig.py --campaign spring24 --tier study --coopt --variant seedlessRefitV3:1111 \
        --events 100 -o spix_spring24.py

    # SPStudyNano on a D121 RelVal: the study tier plus the regular nano (@PHYS); the job also runs
    # RAW2DIGI,RECO,RECOSIM,PAT with the step-3 options of the RelVal workflow (runTheMatrix --dryRun)
    ./makeSpixConfig.py --pu 200 --reco --variant seedlessRefit:1111 --events 2 -o spix_reco.py

No nano table is ever dropped: every table of the tier is produced (the L1
emulation runs on demand, the MINIAOD gen inputs of a truth tier are made in-job).

`--dry-run` prints the cmsDriver command instead of running it, which is what
doc/ConfigProvenance.md records.

DELIBERATELY NOT AN AXIS: otCovMode. seedlessRefit always starts from the OT track's
own helixCovMat. The removed `parametrized` mode substituted a fixed diagonal
for a missing covariance and made every OT-covariance-dependent measurement
meaningless while still producing plausible numbers.
"""
from __future__ import annotations

import argparse
import os
import shlex
import subprocess
import sys

# --- default test inputs, by (pileup, release-compat): file names in the directory $SPIX_TESTFILE_DIR (no
#     default: the directory is the operator's); otherwise pass --filein. See mem:smartpixels-testfile-release-compat.
TESTFILE_ENV = "SPIX_TESTFILE_DIR"
INPUTS = {
    ("200", "20_1"): "RelValTTbar_14TeV_PU_150X_mcRun4_realistic_v1_STD_D121_RegeneratedGS_PU-v1_file1.root",
    ("0", "20_1"): "RelValTTbar_D121_noPU_regen.root",
    # the old Phase2Spring24 TTbar (CMSSW_14_0_9, GEN-SIM-DIGI-RAW-MINIAOD, TrackingParticles above 1 GeV):
    # readable by CMSSW_17, not by CMSSW_20 (io_v1)
    ("200", "17_0"): "TT_TuneCP5_14TeV-powheg-pythia8_Phase2Spring24DIGIRECOMiniAOD_PU200_Trk1GeV_"
                     "140X_mcRun4_realistic_v4-v2.root",
}

# --- campaigns: the input samples' geometry and era, and how the job gets its generator inputs.
#     d121: the CMSSW_20 D121 RelVals (GEN-SIM-DIGI-RAW, no MINIAOD): nano_l1_hlt and produceMiniAODGenInputs.
#     spring24: the Phase2Spring24 samples (D110 = the same T35 tracker; MINIAOD in the input), in a CMSSW_17
#     area: no nano_l1_hlt, and customiseSpring24Inputs (the input fixes and the generator link's pruner maps).
CAMPAIGNS = {
    "d121": dict(geometry="ExtendedRun4D121", era="Phase2C22I13M9", modifiers=["nano_l1_hlt"], release="20_1",
                 gen="process = p3nano.produceMiniAODGenInputs(process)"),
    "spring24": dict(geometry="ExtendedRun4D110", era="Phase2C17I13M9", modifiers=[], release="17_0",
                     gen="process = p3nano.customiseSpring24Inputs(process)"),
}

# --- nano tiers. "truth" variants keep the TrackingParticle extension
#     (L1TTrack_genuine, L1TTrack_tp_*), WITHOUT WHICH no resolution-vs-truth
#     or hit-purity study can run at all. That omission is exactly what made
#     the first layer-order attempt unusable, so the tier names say it.
#     (base L1 tier, payload, withGen). coexist runs the base flavor @<base>[withGen]
#     and adds the payload tables with addPh3L1SmartPixels<payload>; --coopt composes the
#     SmartPixels add-on flavours (COOPT_ADDONS) with the variant's mode and mask tags.
TIERS = {
    "trk":         ("L1TrkNano",   "", False),
    "trk-truth":   ("L1TrkNano",   "", True),
    "pftrk":       ("L1PFTrkNano", "", False),
    "pftrk-truth": ("L1PFTrkNano", "", True),
    "pf":          ("L1PFNano",    "", False),
    "pf-truth":    ("L1PFNano",    "", True),
    # Payload tiers. "clusters" carries the UNTRUNCATED IT cluster table: ~74k
    # rows/event at PU200 (TBPX + ID1-ID4), measured 3.1 MB/event.
    "clusters":       ("L1PFTrkNano", "Clusters", False),
    "clusters-truth": ("L1PFTrkNano", "Clusters", True),
    # --coopt only: the standard L1 SmartPixels study tier (L1TSPStudyNano), and the same
    # plus the regular nano (SPStudyNano, --reco)
    "study":          (None, None, True),
    "reco":           (None, "PHYS", True),
    # coexist only: every --variant as suffixed tables next to the OT tracks, on the add-on
    # flavours (COEXIST_ADDONS); the L1 emulation stays on the OT tracks
    "coexist":        (None, None, True),
}
COEXIST_ADDONS = ["@Phase2L1DPGwithGen", "@L1Trk", "@SPClustersStubs", "@TPTruth1GeVLinked"]
COOPT_ADDONS = {
    "L1TrkNano": ["@L1Trk", "@SPTrk"],
    "L1PFNano": ["@L1PF", "@SPTrk"],
    "L1PFTrkNano": ["@L1PF", "@L1Trk", "@SPTrk"],
}


# seedlessRefit:all = the 15 non-trivial masks (IIII, the passthrough, has no refit)
ALL_SEEDLESSREFIT_MASKS = [a + b + c + d for a in "01" for b in "01" for c in "01" for d in "01"
                           if "1" in a + b + c + d]


def expand_variants(specs):
    return [f"seedlessRefit:{m}" for m in ALL_SEEDLESSREFIT_MASKS] if specs == ["seedlessRefit:all"] else specs


def variant_tags(spec):
    """--variant of a --coopt job -> its mode and mask tags."""
    mode, _, activeSP = spec.partition(":")
    if mode == "passthrough":
        return ["@SPActiveLayersIIII"]
    if mode in ("seedlessRefit", "seedlessRefitV3"):
        tag = "@SPSeedlessRefitV3" if mode == "seedlessRefitV3" else "@SPSeedlessRefit"
        return [tag, "@SPActiveLayers" + activeSP.replace("1", "A").replace("0", "I")]
    raise SystemExit(f"--coopt composes the nano flavours, whose tags select seedlessRefit[V3]:<mask> or "
                     f"passthrough; {spec!r} runs only with coexist")


def nano_step(args):
    base, payload, truth = TIERS[args.tier]
    gen = "withGen" if truth else ""
    if args.tier == "coexist":
        if args.coopt:
            raise SystemExit("--tier coexist runs without --coopt")
        return "NANO:" + "+".join(COEXIST_ADDONS)
    if not args.coopt:
        if base is None:
            raise SystemExit(f"--tier {args.tier} is a --coopt tier")
        return f"NANO:@{base}{gen}"
    tags = variant_tags(args.variant[0])
    if base is None:
        # the aliases carry the mode tag (L1TSPStudyNano: seedlessRefit; L1TSPStudyNanoV3: seedlessRefitV3, the
        # canonical flavour); the variant adds its mask tag
        v3 = "@SPSeedlessRefitV3" in tags
        if payload == "PHYS" and v3:
            raise SystemExit("--reco: SPStudyNano carries seedlessRefit; there is no seedlessRefitV3 alias")
        flavours = ["@SPStudyNano" if payload == "PHYS" else ("@L1TSPStudyNanoV3" if v3 else "@L1TSPStudyNano")]
        tags = [tag for tag in tags if tag.startswith("@SPActiveLayers")]
    else:
        flavours = [f"@Phase2L1DPG{gen}"] + COOPT_ADDONS[base] + (["@SPClustersStubs"] if payload else [])
    return "NANO:" + "+".join(flavours + tags)


CONDITIONS = "auto:phase2_realistic_T35"

# seedlessRefitConfig keys this CLI exposes. Everything else takes SEEDLESSREFIT_DEFAULTS.
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
    """The --customise_commands payload: smartPixelsCoexist (or, with --coopt, the nano
    flavours' parameter block), plus the in-job MINIAOD gen inputs for a truth tier."""
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
        if mode in ("seedlessRefit", "seedlessRefitV3") and not activeSP:
            raise SystemExit(f"--variant {mode} requires an activeSP, e.g. {mode}:1111 (got '{spec}')")
        key = (mode, activeSP)
        if key in seen:
            raise SystemExit(f"--variant {spec} given more than once; each variant makes an "
                             "identically-named set of tables and the duplicate would collide")
        seen.add(key)
        variants.append(f'("{mode}", "{activeSP}")' if activeSP else f'("{mode}", None)')

    # the nano flavours default to the emu angle source, coexist to the PixelAV payload
    angle_source = args.angle_source or ("emu" if args.coopt or args.tier == "coexist" else "payload")
    refit = {"angleSource": angle_source}
    if angle_source == "payload":
        refit["pixelavAngleSet"] = args.pixelav_angle_set
    refit.update(overrides)
    posture = (f"trackInputMode='{args.track_input_mode}', extendedTracks={args.extended_tracks}, "
               f"promptHnpar={args.prompt_hnpar}")
    _, payload, truth = TIERS[args.tier]
    parts = ["import DPGAnalysis.Phase3SmartPixelsNanoAOD.l1tPh3SmartPixelsNano_cff as p3nano"]
    if args.coopt:
        # the tags in the NANO step wire the variant; the rest goes through the parameter block
        params = ", ".join(f"{k}={v!r}" for k, v in sorted(refit.items()))
        parts += ["from DPGAnalysis.Phase3SmartPixelsNanoAOD.l1tPh3SmartPixelsNanoFlavours_cff import smartPixelsNanoParams",
                  f"smartPixelsNanoParams(process, {posture}, {params})"]
    else:
        parts.append("from L1Trigger.Phase3SmartPixels.customizeSmartPixels_cff import smartPixelsCoexist")
        if payload:
            parts.append(f"process = p3nano.addPh3L1SmartPixels{payload}(process)")
        parts.append(f"process = smartPixelsCoexist(process, variants=[{', '.join(variants)}], "
                     f"addNanoTables=True, {posture}, seedlessRefitConfig={refit!r})")
    # withGen tiers read MINIAOD gen products (prunedGenParticles, slimmedGenJets, ...);
    # the GEN-SIM-DIGI-RAW RelVals carry only genParticles/ak4GenJetsNoNu, so make them.
    # SPStudyNano gets them from its MINIAOD input or from the PAT step it runs.
    # The Spring24 inputs carry MINIAOD; their input fixes apply to every tier.
    if args.campaign == "spring24" or (truth and args.tier != "reco"):
        parts.append(CAMPAIGNS[args.campaign]["gen"])
    return "; ".join(parts)


def cmsdriver_cmd(args, out_py, overrides):
    if args.needs_truth and not TIERS[args.tier][2]:
        raise SystemExit(
            f"--tier {args.tier} does not keep the TrackingParticle truth extension, but "
            "--needs-truth was requested. Use a '-truth' tier; otherwise L1TTrack_genuine "
            "and L1TTrack_tp_* are absent and no resolution or hit-purity study can run.")
    infile = args.filein
    if not infile:
        name, testdir = INPUTS.get((args.pu, args.release_compat)), os.environ.get(TESTFILE_ENV)
        if not (name and testdir):
            raise SystemExit(f"no input: pass --filein, or set ${TESTFILE_ENV} to the directory of the default test "
                             f"input for pu={args.pu} release={args.release_compat}" + (f" ({name})" if name else ""))
        infile = os.path.join(testdir, name)
    if args.tier == "reco":
        return reco_cmd(args, out_py, overrides, infile)
    camp = CAMPAIGNS[args.campaign]
    return [
        "cmsDriver.py", os.path.splitext(os.path.basename(out_py))[0],
        "-s", nano_step(args),
        "--conditions", CONDITIONS,
        "--geometry", camp["geometry"],
        "--era", camp["era"],
    ] + (["--procModifiers", ",".join(camp["modifiers"])] if camp["modifiers"] else []) + [
        "--datatier", "NANOAODSIM",
        "--eventcontent", "NANOAODSIM",
        "--filein", f"file:{infile}",
        "--fileout", f"file:{args.fileout or os.path.splitext(out_py)[0] + '.root'}",
        "--python_filename", out_py,
        "-n", str(args.events),
        "--mc", "--nThreads", str(args.threads), "--no_exec",
        "--customise_commands", build_customise(args, overrides),
    ]


# The RelVal workflow the test inputs come from: TTbar_14TeV, Run4D121, PU200.
RELVAL_WORKFLOW = "34634.0"
# step-3 options the SmartPixels job sets itself; --pileup/--pileup_input also go, since the
# job reads the mixed TrackingParticles from the input and the pileup files are on DAS
STEP3_OWN_OPTIONS = {"-s", "--step", "-n", "--filein", "--fileout", "--datatier", "--eventcontent",
                     "--pileup", "--pileup_input", "--python_filename", "--no_exec", "--nThreads"}


def relval_step3(workflow):
    """(steps up to PAT, [other cmsDriver options]) of the workflow's step 3, from
    runTheMatrix --dryRun: the RECO the RelVal itself was made with."""
    import glob
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        subprocess.run(["runTheMatrix.py", "-w", "upgrade", "-l", workflow, "--dryRun"], cwd=tmp,
                       check=True, capture_output=True)
        logs = glob.glob(os.path.join(tmp, "*", "cmdLog"))
        lines = [line for log in logs for line in open(log) if "cmsDriver.py step3" in line]
    if not lines:
        raise SystemExit(f"runTheMatrix --dryRun of workflow {workflow} gave no step 3")
    tokens = shlex.split(lines[0].split(">")[0])
    tokens = tokens[tokens.index("cmsDriver.py") + 2:]
    steps, options, i = None, [], 0
    while i < len(tokens):
        opt = tokens[i]
        value = tokens[i + 1] if i + 1 < len(tokens) and not tokens[i + 1].startswith("-") else None
        i += 2 if value is not None else 1
        if opt in ("-s", "--step"):
            steps = value.split(",")
        elif opt not in STEP3_OWN_OPTIONS:
            options.append((opt, value))
    if steps is None or "PAT" not in steps:
        raise SystemExit(f"step 3 of workflow {workflow} runs no PAT: {steps}")
    return steps[:steps.index("PAT") + 1], options


def reco_cmd(args, out_py, overrides, infile):
    """SPStudyNano on a D121 RelVal (GEN-SIM-DIGI-RAW, no MINIAOD products): RECO and PAT in the job, with the
    RelVal workflow's step-3 options. (On a MINIAOD input the job is NANO only, with
    smartPixelsNanoParams(process, inputHasMiniAOD=True): see the README.)"""
    steps, options = relval_step3(args.reco_workflow)
    nano = ",".join(steps + [nano_step(args)])
    modifiers = ["nano_l1_hlt"]
    cmd = ["cmsDriver.py", os.path.splitext(os.path.basename(out_py))[0], "-s", nano]
    for opt, value in options:
        if opt == "--procModifiers":
            modifiers = value.split(",") + modifiers
        else:
            cmd += [opt] + ([value] if value is not None else [])
    return cmd + [
        "--procModifiers", ",".join(modifiers),
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
    ap.add_argument("--campaign", default="d121", choices=sorted(CAMPAIGNS),
                    help="the input samples: d121 (the CMSSW_20 RelVals) or spring24 (the old Phase2Spring24 "
                         "samples, read by a CMSSW_17 area)")
    ap.add_argument("--release-compat", default=None, choices=["20_1", "17_0"],
                    help="which release the input must be readable by (default: the campaign's; "
                         "see mem:smartpixels-testfile-release-compat)")
    ap.add_argument("--tier", choices=sorted(TIERS), default="trk-truth")
    ap.add_argument("--variant", action="append", default=None,
                    help="mode[:activeSP], e.g. seedlessRefit:1111, passthrough. REPEATABLE: give it "
                         "several times to emit one table set per variant from a single pass over "
                         "the input, which is how the activeSP sweep avoids 15 reconstructions; "
                         "seedlessRefit:all = the 15 non-trivial masks")
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
                         "helixCovMat makes every extended OT track unusable for the refit.")
    ap.add_argument("--prompt-hnpar", type=int, default=5, choices=[4, 5])
    ap.add_argument("--use-angles", default=None, choices=["none", "alpha", "alphaBeta"])
    ap.add_argument("--angle-source", default=None, choices=["emu", "payload"],
                    help="rec-hit angle source (default: emu with the nano flavours, i.e. --coopt "
                         "or --tier coexist; payload otherwise)")
    ap.add_argument("--pixelav-angle-set", default=DEFAULT_ANGLE_SET)
    ap.add_argument("--reco", action="store_true",
                    help="SPStudyNano (--tier reco --coopt) on a D121 RelVal: the study tier plus the regular nano, "
                         "with RECO and PAT in the job")
    ap.add_argument("--reco-workflow", default=RELVAL_WORKFLOW,
                    help="RelVal workflow whose step-3 options the RECO+PAT steps take (%(default)s)")
    ap.add_argument("--coopt", action="store_true",
                    help="canonical single-view nano (the SmartPixels nano flavours): the one "
                         "--variant replaces the tracklet tracks for the whole L1 emulation")
    ap.add_argument("--filein", default=None)
    ap.add_argument("--fileout", default=None)
    ap.add_argument("--needs-truth", action="store_true",
                    help="fail unless the tier keeps L1TTrack_genuine / tp_*")
    ap.add_argument("--set", action="append", default=[], metavar="KEY=VALUE",
                    help=f"seedlessRefitConfig override; one of {sorted(REFIT_AXES)}")
    ap.add_argument("--scan", default=None, metavar="KEY=V1,V2[,...]",
                    help="emit one config per value, identical otherwise (for A/B)")
    ap.add_argument("-o", "--output", required=True,
                    help="output .py path; with --scan, a prefix (_<value>.py appended)")
    ap.add_argument("--dry-run", action="store_true",
                    help="print the cmsDriver command instead of running it")
    args = ap.parse_args()
    release = CAMPAIGNS[args.campaign]["release"]
    if args.release_compat not in (None, release):
        raise SystemExit(f"--campaign {args.campaign} inputs are read by {release}, not {args.release_compat}")
    args.release_compat = release
    if args.campaign == "spring24" and args.reco:
        raise SystemExit("--reco takes the D121 RelVal workflow's step-3 options; it has no spring24 form")
    if not args.variant:
        args.variant = ["seedlessRefit:1111"]
    args.variant = expand_variants(args.variant)
    if args.reco:
        args.tier, args.coopt = "reco", True

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
            raise SystemExit(f"unknown seedlessRefit axis {k!r}; known: {sorted(REFIT_AXES)}")
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
