#!/usr/bin/env python3
"""Pure-python (no cmsRun) checks of the SmartPixels refit-mode config surface.

Exercises Phase-0 of the tier-2 refit plan in
L1Trigger/Phase3SmartPixels/python/customizeSmartPixels_cff.py:

 (a) suffix/label generation for digiRefit variants works and labels still end
     so that derived nano tables end in 'Table' (keep-pattern gotcha);
 (b) a digiRefit variant wires the producer (Phase 2): params land on the module
     and no RandomNumberGeneratorService is required (the producer seeds a local
     per-event RNG from hash(label,run,lumi,event)); two identical builds produce
     identical configs; an empty pixelavAngleSet raises ValueError;
 (c) mode 'refit' (Tier 3) raises the RESERVED NotImplementedError;
 (d) invalid digiRefitConfig keys / enum values raise ValueError;
 (e) existing modes (passthrough, correctionlibRegression) still normalize and
     suffix EXACTLY as before (byte-identical regression guard).

Run inside the container:
  source /cvmfs/cms.cern.ch/cmsset_default.sh
  cd $CMSSW_BASE/src && eval $(scramv1 runtime -sh)
  python3 L1Trigger/Phase3SmartPixels/test/testRefitModesConfig.py
"""
import sys

import FWCore.ParameterSet.Config as cms

from L1Trigger.Phase3SmartPixels.customizeSmartPixels_cff import (
    ACTIVE_LAYER_MODES,
    ALL_MODES,
    DIGIREFIT_DEFAULTS,
    PASSTHROUGH_MODES,
    RESERVED_MODES,
    TRUTHSOURCE_CHOICES,
    TRUTH_REQUIRED_MODES,
    _applyTruthSource,
    _normalizeVariants,
    _resolveDigiRefitConfig,
    addSmartPixelsTrackProducerVariants,
    attachFromFileStubsChain,
    smartPixelsCoexist,
    smartPixelsVariantLabels,
    smartPixelsVariantSuffix,
)

_failures = []


def check(cond, msg):
    if cond:
        print("  ok:   " + msg)
    else:
        print("  FAIL: " + msg)
        _failures.append(msg)


def expect_raises(exc_type, fn, msg, contains=None):
    try:
        fn()
    except exc_type as e:
        if contains is not None and contains not in str(e):
            print(f"  FAIL: {msg} -- raised {exc_type.__name__} but message "
                  f"{str(e)!r} missing {contains!r}")
            _failures.append(msg)
        else:
            print(f"  ok:   {msg} (raised {exc_type.__name__})")
    except Exception as e:  # noqa: BLE001
        print(f"  FAIL: {msg} -- expected {exc_type.__name__}, got "
              f"{type(e).__name__}: {e}")
        _failures.append(msg)
    else:
        print(f"  FAIL: {msg} -- no exception raised")
        _failures.append(msg)


# --- (a) digiRefit suffix / label generation -------------------------------
print("[a] digiRefit suffix/label generation")
check("digiRefit" in ACTIVE_LAYER_MODES, "digiRefit is an active-layer mode")
check("digiRefit" in TRUTH_REQUIRED_MODES, "digiRefit is truth-required")
check("refit" in RESERVED_MODES, "refit is a reserved mode")
check("digiRefit" in ALL_MODES and "refit" in ALL_MODES,
      "both digiRefit and refit are in the mode vocabulary")

suffix = smartPixelsVariantSuffix("digiRefit", "1100")
check(suffix == "digiRefitAAII",
      f"smartPixelsVariantSuffix('digiRefit','1100') == 'digiRefitAAII' (got {suffix!r})")

prompt, extended = smartPixelsVariantLabels("digiRefit", "1100")
check(prompt == "l1tSmartPixelsTrackProducerWdigiRefitAAII",
      f"prompt label (got {prompt!r})")
check(extended == "l1tSmartPixelsTrackProducerExtendedWdigiRefitAAII",
      f"extended label (got {extended!r})")
# keep-pattern gotcha: nano table labels are formed as <src>+suffix+'Table';
# emulate that the derived table label ends in 'Table'.
table_label = "l1tPh3SmartPixelsTracks" + suffix + "Table"
check(table_label.endswith("Table"),
      f"derived nano table label ends in 'Table' (got {table_label!r})")

expect_raises(ValueError, lambda: smartPixelsVariantSuffix("digiRefit", None),
              "digiRefit without activeSP raises ValueError")

# --- (b) digiRefit variant wires the producer (Phase 2) ----------------------
print("[b] digiRefit -> wired producer variants (Phase 2)")
_DR_CFG = {"pixelavAngleSet": "dummy/pixelav_angle_example.json"}

p_dr = cms.Process("TEST")
p_dr, _dr_mods = addSmartPixelsTrackProducerVariants(p_dr, variants=[("digiRefit", "1100")],
                                                     digiRefitConfig=_DR_CFG)
_dr_prompt = "l1tSmartPixelsTrackProducerWdigiRefitAAII"
check(hasattr(p_dr, _dr_prompt), "digiRefit prompt producer variant exists")
_m = getattr(p_dr, _dr_prompt)
check(_m.smartPixelsEmulatorMode.value() == "digiRefit", "mode set to digiRefit")
check(_m.digiRefitSeedCovMode.value() == "trackCov", "seedCovMode default is trackCov")
check(_m.digiRefitSeedNPar.value() == 5, "seedNPar default is 5")
check(_m.digiRefitUseAngles.value() == "alpha", "useAngles default is alpha")
check(_m.digiRefitPixelavAngleSet.value() == _DR_CFG["pixelavAngleSet"],
      "pixelavAngleSet carried onto the module")
check(list(_m.digiRefitParamSigmas) == list(DIGIREFIT_DEFAULTS["paramSigmas"]),
      "paramSigmas carried onto the module")
# KF numerical guards (spec §6b) default onto the module and take the
# investigation-chosen values.
check(_m.digiRefitJacobianMaxAbs.value() == DIGIREFIT_DEFAULTS["jacobianMaxAbs"] == 1.0e4,
      "jacobianMaxAbs default (1e4) carried onto the module")
check(_m.digiRefitChi2UpdateGate.value() == DIGIREFIT_DEFAULTS["chi2UpdateGate"] == 2.0e6,
      "chi2UpdateGate default (2e6) carried onto the module")
# and overrides flow through
_m_ovr = getattr(
    addSmartPixelsTrackProducerVariants(
        cms.Process("TEST"), variants=[("digiRefit", "1100")],
        digiRefitConfig={"pixelavAngleSet": "dummy/x.json", "jacobianMaxAbs": 5.0e5,
                         "chi2UpdateGate": 3.0e6})[0],
    _dr_prompt)
check(_m_ovr.digiRefitJacobianMaxAbs.value() == 5.0e5
      and _m_ovr.digiRefitChi2UpdateGate.value() == 3.0e6,
      "jacobianMaxAbs / chi2UpdateGate overrides flow onto the module")
# Refit-quality BDT model (spec §6a): default empty (no scoring), override flows.
check(_m.digiRefitBdtModel.value() == DIGIREFIT_DEFAULTS["bdtModel"] == "",
      "bdtModel default is empty (no BDT scoring) and carried onto the module")
_m_bdt = getattr(
    addSmartPixelsTrackProducerVariants(
        cms.Process("TEST"), variants=[("digiRefit", "1100")],
        digiRefitConfig={"pixelavAngleSet": "dummy/x.json",
                         "bdtModel": "L1Trigger/Phase3SmartPixels/data/refitq_model.json"})[0],
    _dr_prompt)
check(_m_bdt.digiRefitBdtModel.value() == "L1Trigger/Phase3SmartPixels/data/refitq_model.json",
      "bdtModel override flows onto the module")
# The RNG scheme is now producer-side (local per-event engine seeded from
# hash(label,run,lumi,event)); NO RandomNumberGeneratorService is wired.
check(not hasattr(p_dr, "RandomNumberGeneratorService"),
      "no RandomNumberGeneratorService is added (RNG is per-event, producer-side)")
# Determinism at config level: two identical builds produce identical module
# configs (the dump strings match byte-for-byte).
p_dr2 = cms.Process("TEST")
p_dr2, _ = addSmartPixelsTrackProducerVariants(p_dr2, variants=[("digiRefit", "1100")],
                                               digiRefitConfig=_DR_CFG)
check(getattr(p_dr2, _dr_prompt).dumpPython() == _m.dumpPython(),
      "two identical digiRefit builds produce identical module configs")

# a digiRefit variant WITHOUT a payload must fail loudly at config time
def _build_digirefit_nopayload():
    p = cms.Process("TEST")
    addSmartPixelsTrackProducerVariants(p, variants=[("digiRefit", "1100")])
expect_raises(ValueError, _build_digirefit_nopayload,
              "digiRefit without pixelavAngleSet raises ValueError",
              contains="pixelavAngleSet")

# --- (c) mode 'refit' raises the RESERVED error -----------------------------
print("[c] refit (Tier 3) -> reserved NotImplementedError")
expect_raises(NotImplementedError, lambda: smartPixelsVariantSuffix("refit"),
              "smartPixelsVariantSuffix('refit') raises reserved NotImplementedError",
              contains="RESERVED")
def _coexist_refit():
    p = cms.Process("TEST")
    smartPixelsCoexist(p, variants=[("refit", None)])
expect_raises(NotImplementedError, _coexist_refit,
              "smartPixelsCoexist(refit) raises reserved NotImplementedError",
              contains="Tier 3")

# --- (d) invalid digiRefitConfig -------------------------------------------
print("[d] digiRefitConfig validation")
expect_raises(ValueError,
              lambda: _resolveDigiRefitConfig({"notAKey": 1}),
              "unknown digiRefitConfig key raises ValueError",
              contains="unknown key")
expect_raises(ValueError,
              lambda: _resolveDigiRefitConfig({"useAngles": "bogus"}),
              "invalid useAngles enum raises ValueError",
              contains="useAngles")
expect_raises(ValueError,
              lambda: _resolveDigiRefitConfig({"gainMode": "bogus"}),
              "invalid gainMode enum raises ValueError",
              contains="gainMode")
expect_raises(ValueError,
              lambda: _resolveDigiRefitConfig(["not", "a", "dict"]),
              "non-dict digiRefitConfig raises ValueError",
              contains="must be a dict")
# KF numerical guards must be positive doubles (spec §6b).
expect_raises(ValueError,
              lambda: _resolveDigiRefitConfig({"jacobianMaxAbs": -1.0}),
              "non-positive jacobianMaxAbs raises ValueError",
              contains="jacobianMaxAbs")
expect_raises(ValueError,
              lambda: _resolveDigiRefitConfig({"chi2UpdateGate": 0}),
              "zero chi2UpdateGate raises ValueError",
              contains="chi2UpdateGate")
expect_raises(ValueError,
              lambda: _resolveDigiRefitConfig({"chi2UpdateGate": "big"}),
              "non-numeric chi2UpdateGate raises ValueError",
              contains="chi2UpdateGate")
# Refit-quality BDT model must be a string path (spec §6a).
expect_raises(ValueError,
              lambda: _resolveDigiRefitConfig({"bdtModel": 17}),
              "non-string bdtModel raises ValueError",
              contains="bdtModel")

# valid config merges over defaults and does NOT mutate the defaults
merged = _resolveDigiRefitConfig({"minHits": 3, "useAngles": "alphaBeta"})
check(merged["minHits"] == 3 and merged["useAngles"] == "alphaBeta",
      "valid overrides applied")
check(merged["windowRPhi"] == DIGIREFIT_DEFAULTS["windowRPhi"],
      "untouched keys keep defaults")
check(DIGIREFIT_DEFAULTS["minHits"] != 3,
      "DIGIREFIT_DEFAULTS not mutated by merge")
check(set(merged) == set(DIGIREFIT_DEFAULTS),
      "merged config has exactly the default key set")

# an invalid digiRefitConfig must fail loudly even before the digiRefit
# variant hits its (Phase-0) NotImplementedError
def _coexist_bad_cfg():
    p = cms.Process("TEST")
    smartPixelsCoexist(p, variants=[("digiRefit", "1100")],
                       digiRefitConfig={"notAKey": 1})
expect_raises(ValueError, _coexist_bad_cfg,
              "invalid digiRefitConfig raises ValueError via coexist (before NotImplemented)",
              contains="unknown key")

# --- (e) existing modes unchanged (byte-identical regression guard) ---------
print("[e] existing modes normalize/suffix exactly as before")
check(smartPixelsVariantSuffix("passthrough") == "passthrough",
      "passthrough suffix unchanged")
check(smartPixelsVariantSuffix("passthroughFloat") == "passthroughFloat",
      "passthroughFloat suffix unchanged")
check(smartPixelsVariantSuffix("trackingParticleTruth") == "trackingParticleTruth",
      "trackingParticleTruth suffix unchanged")
check(smartPixelsVariantSuffix("correctionlibRegression", "1100") == "correctionlibRegressionAAII",
      "correctionlibRegression '1100' suffix unchanged")
check(smartPixelsVariantSuffix("correctionlibRegression", "0000") == "correctionlibRegressionIIII",
      "correctionlibRegression '0000' suffix unchanged")
lbl = smartPixelsVariantLabels("correctionlibRegression", "1100")
check(lbl == ("l1tSmartPixelsTrackProducerWcorrectionlibRegressionAAII",
              "l1tSmartPixelsTrackProducerExtendedWcorrectionlibRegressionAAII"),
      "correctionlibRegression labels unchanged")
check(_normalizeVariants(["passthrough", ("correctionlibRegression", "1100")])
      == [("passthrough", None), ("correctionlibRegression", "1100")],
      "_normalizeVariants unchanged for existing inputs")
# unknown mode still rejected
expect_raises(ValueError, lambda: smartPixelsVariantSuffix("totallyBogusMode"),
              "unknown mode raises ValueError", contains="unknown SmartPixels mode")

# --- (f) truthSource="fromFileStubs" (posture C) wiring ---------------------
print("[f] truthSource='fromFileStubs' (posture C) wiring")
check("fromFileStubs" in TRUTHSOURCE_CHOICES,
      "fromFileStubs is in the truthSource vocabulary")
check(set(TRUTHSOURCE_CHOICES) == {"inJob", "fromFile", "fromFileStubs"},
      "truthSource vocabulary is exactly {inJob, fromFile, fromFileStubs}")


def _build_fromfilestubs(extendedTracks, seedCovMode="trackCov"):
    """Full coexist build with truthSource=fromFileStubs on a process carrying a
    dummy Path (so the posture-C Task has somewhere to associate)."""
    p = cms.Process("TEST")
    p.pdummy = cms.Path()  # somewhere for the posture-C Task to associate
    return smartPixelsCoexist(
        p, variants=[("digiRefit", "1100")], addNanoTables=False,
        truthSource="fromFileStubs", extendedTracks=extendedTracks,
        digiRefitConfig={"pixelavAngleSet": "dummy/x.json", "seedCovMode": seedCovMode})


# extended ON (default): full prompt+extended posture-C chain present
p_ffs = _build_fromfilestubs(extendedTracks=True)
_PROMPT_CHAIN = ["ProducerDTC", "l1tTTTracksFromTrackletEmulation",
                 "TTTrackAssociatorFromPixelDigis"]
_EXT_CHAIN = ["l1tTTTracksFromExtendedTrackletEmulation",
              "TTTrackAssociatorFromPixelDigisExtended"]
for m in _PROMPT_CHAIN:
    check(hasattr(p_ffs, m), f"posture-C chain module '{m}' present (default labels)")
for m in _EXT_CHAIN:
    check(hasattr(p_ffs, m), f"posture-C extended chain module '{m}' present (extendedTracks=True)")
# the re-run track associator points at the NEW tracklet tracks (default label)
check(list(p_ffs.TTTrackAssociatorFromPixelDigis.TTTracks)[0]
      == cms.InputTag("l1tTTTracksFromTrackletEmulation", "Level1TTTracks"),
      "prompt TTTrackAssociatorFromPixelDigis re-run against the new tracklet tracks")
check(list(p_ffs.TTTrackAssociatorFromPixelDigisExtended.TTTracks)[0]
      == cms.InputTag("l1tTTTracksFromExtendedTrackletEmulation", "Level1TTTracks"),
      "extended associator re-run against the new extended tracklet tracks")
# the cluster/stub map tags carry NO process name -> resolve to the file's HLT maps
check(p_ffs.TTTrackAssociatorFromPixelDigis.TTClusterTruth.getProcessName() == "",
      "TTClusterTruth has no process name (resolves to file HLT maps)")
check(p_ffs.TTTrackAssociatorFromPixelDigis.TTStubTruth.getProcessName() == "",
      "TTStubTruth has no process name (resolves to file HLT maps)")
# DIGI-tier + cluster/stub associators are ABSENT (never scheduled in-process)
for m in ("TTClusterAssociatorFromPixelDigis", "TTStubAssociatorFromPixelDigis",
          "simSiPixelDigis", "mix", "offlineBeamSpot"):
    check(not hasattr(p_ffs, m),
          f"DIGI-tier/cluster-stub module '{m}' NOT scheduled in-process (posture C)")
# the posture-C Task exists and contains exactly the chain modules
check(hasattr(p_ffs, "l1tSmartPixelsFromFileStubsTask"),
      "posture-C Task l1tSmartPixelsFromFileStubsTask exists")
_task_names = {mod.label_() for mod in p_ffs.l1tSmartPixelsFromFileStubsTask._collection}
check(_task_names == set(_PROMPT_CHAIN + _EXT_CHAIN),
      f"posture-C Task holds exactly prompt+extended chain (got {sorted(_task_names)})")

# extended OFF: prompt chain only. NOTE: the extended EMULATOR module
# (l1tTTTracksFromExtendedTrackletEmulation) is still defined on the process --
# the l1tTTTracksFromTrackletEmulation_cfi import brings both prompt and extended
# emulator objects regardless -- but it is NOT on the posture-C Task and its
# associator is not cloned. The knob's load-bearing effect is Task membership.
p_ffs_p = _build_fromfilestubs(extendedTracks=False)
for m in _PROMPT_CHAIN:
    check(hasattr(p_ffs_p, m), f"prompt-only build still has '{m}'")
check(not hasattr(p_ffs_p, "TTTrackAssociatorFromPixelDigisExtended"),
      "extendedTracks=False -> extended track associator NOT cloned")
_task_names_p = {mod.label_() for mod in p_ffs_p.l1tSmartPixelsFromFileStubsTask._collection}
check(_task_names_p == set(_PROMPT_CHAIN),
      f"prompt-only posture-C Task holds exactly the prompt chain (got {sorted(_task_names_p)})")
check("l1tTTTracksFromExtendedTrackletEmulation" not in _task_names_p,
      "extendedTracks=False -> extended emulator NOT on the posture-C Task")
# extended ON vs OFF changes the scheduled module set (knob is live)
_task_names_on = {mod.label_() for mod in p_ffs.l1tSmartPixelsFromFileStubsTask._collection}
check(_task_names_on != _task_names_p and _EXT_CHAIN[0] in _task_names_on,
      "extendedTracks on/off changes the scheduled posture-C module set")

# unknown truthSource still raises loudly
def _bad_truthsource():
    p = cms.Process("TEST")
    _applyTruthSource(p, "bogusPosture", [("passthrough", None)])
expect_raises(ValueError, _bad_truthsource,
              "unknown truthSource raises ValueError", contains="truthSource must be one of")

# attachFromFileStubsChain returns the chain module list (prompt-only when off)
p_direct = cms.Process("TEST")
p_direct.pdummy = cms.Path()
_p, _chain = attachFromFileStubsChain(p_direct, extendedTracks=False)
check(_chain == _PROMPT_CHAIN,
      f"attachFromFileStubsChain(extendedTracks=False) returns prompt chain (got {_chain})")

# --- summary ---------------------------------------------------------------
print()
if _failures:
    print(f"RESULT: FAILED ({len(_failures)} failing check(s))")
    for f in _failures:
        print("  - " + f)
    sys.exit(1)
print("RESULT: PASSED (all checks green)")
sys.exit(0)
