#!/usr/bin/env python3
"""Pure-python (no cmsRun) checks of the SmartPixels refit-mode config surface.

Exercises Phase-0 of the tier-2 refit plan in
L1Trigger/Phase3SmartPixels/python/customizeSmartPixels_cff.py:

 (a) suffix/label generation for seedlessRefit variants works and labels still end
     so that derived nano tables end in 'Table' (keep-pattern gotcha);
 (b) a seedlessRefit variant wires the producer (Phase 2): params land on the module
     and no RandomNumberGeneratorService is required (the producer seeds a local
     per-event RNG from hash(label,run,lumi,event)); two identical builds produce
     identical configs; an empty pixelavAngleSet raises ValueError;
 (c) mode 'refit' (Tier 3) raises the RESERVED NotImplementedError;
 (d) invalid seedlessRefitConfig keys / enum values raise ValueError;
 (e) existing modes (passthrough, correctionlibRegression) still normalize and
     suffix EXACTLY as before (byte-identical regression guard);
 (i) one cluster table for barrel + discs: discScope validated, rec hits and
     table agree in either customise order, no separate disc table; stub tables
     carry the OT stub collection for otStubIdx; L1TTP cuts on |p|; a cut on a
     row-aligned track table is refused;
 (k) coexist variants get their track -> stub / cluster link tables in either order.

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
    SEEDLESSREFIT_DEFAULTS,
    PASSTHROUGH_MODES,
    PROMPT_HNPAR_CHOICES,
    RESERVED_MODES,
    TRACKINPUTMODE_CHOICES,
    TRUTH_REQUIRED_MODES,
    _applyTrackInputMode,
    _normalizeVariants,
    _resolveSeedlessRefitConfig,
    addSmartPixelsTrackProducerVariants,
    attachStubRebuildChain,
    smartPixelsCoexist,
    smartPixelsCoopt,
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


# --- (a) seedlessRefit suffix / label generation -------------------------------
print("[a] seedlessRefit suffix/label generation")
check("seedlessRefit" in ACTIVE_LAYER_MODES, "seedlessRefit is an active-layer mode")
check("seedlessRefit" in TRUTH_REQUIRED_MODES, "seedlessRefit is truth-required")
check("refit" in RESERVED_MODES, "refit is a reserved mode")
check("seedlessRefit" in ALL_MODES and "refit" in ALL_MODES,
      "both seedlessRefit and refit are in the mode vocabulary")

suffix = smartPixelsVariantSuffix("seedlessRefit", "1100")
check(suffix == "seedlessRefitAAII",
      f"smartPixelsVariantSuffix('seedlessRefit','1100') == 'seedlessRefitAAII' (got {suffix!r})")

prompt, extended = smartPixelsVariantLabels("seedlessRefit", "1100")
check(prompt == "l1tSmartPixelsTrackProducerWseedlessRefitAAII",
      f"prompt label (got {prompt!r})")
check(extended == "l1tSmartPixelsTrackProducerExtendedWseedlessRefitAAII",
      f"extended label (got {extended!r})")
# keep-pattern gotcha: nano table labels are formed as <src>+suffix+'Table';
# emulate that the derived table label ends in 'Table'.
table_label = "l1tPh3SmartPixelsTracks" + suffix + "Table"
check(table_label.endswith("Table"),
      f"derived nano table label ends in 'Table' (got {table_label!r})")

expect_raises(ValueError, lambda: smartPixelsVariantSuffix("seedlessRefit", None),
              "seedlessRefit without activeSP raises ValueError")

# --- (b) seedlessRefit variant wires the producer (Phase 2) ----------------------
print("[b] seedlessRefit -> wired producer variants (Phase 2)")
_DR_CFG = {"pixelavAngleSet": "dummy/pixelav_angle_example.json"}

p_dr = cms.Process("TEST")
p_dr, _dr_mods = addSmartPixelsTrackProducerVariants(p_dr, variants=[("seedlessRefit", "1100")],
                                                     seedlessRefitConfig=_DR_CFG)
_dr_prompt = "l1tSmartPixelsTrackProducerWseedlessRefitAAII"
check(hasattr(p_dr, _dr_prompt), "seedlessRefit prompt producer variant exists")
_m = getattr(p_dr, _dr_prompt)
check(_m.smartPixelsEmulatorMode.value() == "seedlessRefit", "mode set to seedlessRefit")
check(_m.seedlessRefitOTNPar.value() == 5, "otNPar default is 5")
check(_m.seedlessRefitUseAngles.value() == "alphaBeta", "useAngles default is alphaBeta (the production setting)")
check(_m.smartPixelsRecHitInputTag.value() == "spixSmartPixelsRecHits",
      "refit reads the SmartPixelsRecHit product, not raw pixel rec hits")
check(hasattr(p_dr, "spixSmartPixelsRecHits"),
      "the upstream SmartPixelsRecHitProducer is scheduled")
check(p_dr.spixSmartPixelsRecHits.angleSet.value() == _DR_CFG["pixelavAngleSet"],
      "pixelavAngleSet configures the RECHIT producer (single source of the angle)")
# KF numerical guards (spec §6b) default onto the module and take the
# investigation-chosen values.
check(_m.seedlessRefitJacobianMaxAbs.value() == SEEDLESSREFIT_DEFAULTS["jacobianMaxAbs"] == 1.0e4,
      "jacobianMaxAbs default (1e4) carried onto the module")
check(_m.seedlessRefitChi2UpdateGate.value() == SEEDLESSREFIT_DEFAULTS["chi2UpdateGate"] == 2.0e6,
      "chi2UpdateGate default (2e6) carried onto the module")
# and overrides flow through
_m_ovr = getattr(
    addSmartPixelsTrackProducerVariants(
        cms.Process("TEST"), variants=[("seedlessRefit", "1100")],
        seedlessRefitConfig={"pixelavAngleSet": "dummy/x.json", "jacobianMaxAbs": 5.0e5,
                         "chi2UpdateGate": 3.0e6})[0],
    _dr_prompt)
check(_m_ovr.seedlessRefitJacobianMaxAbs.value() == 5.0e5
      and _m_ovr.seedlessRefitChi2UpdateGate.value() == 3.0e6,
      "jacobianMaxAbs / chi2UpdateGate overrides flow onto the module")
# Grazing-angle clamps (spec §6b): defaults 12.0 carried, overrides flow.
check(_m.seedlessRefitMeasAngleMaxAbs.value() == SEEDLESSREFIT_DEFAULTS["measAngleMaxAbs"] == 12.0,
      "measAngleMaxAbs default (12.0) carried onto the module")
check(_m.seedlessRefitPredAngleMaxAbs.value() == SEEDLESSREFIT_DEFAULTS["predAngleMaxAbs"] == 12.0,
      "predAngleMaxAbs default (12.0) carried onto the module")
_m_clamp = getattr(
    addSmartPixelsTrackProducerVariants(
        cms.Process("TEST"), variants=[("seedlessRefit", "1100")],
        seedlessRefitConfig={"pixelavAngleSet": "dummy/x.json", "measAngleMaxAbs": 8.0,
                         "predAngleMaxAbs": 15.0})[0],
    _dr_prompt)
check(_m_clamp.seedlessRefitMeasAngleMaxAbs.value() == 8.0
      and _m_clamp.seedlessRefitPredAngleMaxAbs.value() == 15.0,
      "measAngleMaxAbs / predAngleMaxAbs overrides flow onto the module")
# Kalman layer visit order: the DEFAULT is outsideIn (L4->L1), because the
# projection to a layer carries the unmodelled scattering between that layer and
# the OT, which grows inward. insideOut remains selectable for the comparison.
check(_m.seedlessRefitLayerOrder.value() == SEEDLESSREFIT_DEFAULTS["layerOrder"] == "outsideIn",
      "layerOrder default (outsideIn) carried onto the module")
_m_order = getattr(
    addSmartPixelsTrackProducerVariants(
        cms.Process("TEST"), variants=[("seedlessRefit", "1100")],
        seedlessRefitConfig={"pixelavAngleSet": "dummy/x.json", "layerOrder": "insideOut"})[0],
    _dr_prompt)
check(_m_order.seedlessRefitLayerOrder.value() == "insideOut",
      "layerOrder override (insideOut) flows onto the module")
try:
  addSmartPixelsTrackProducerVariants(
      cms.Process("TEST"), variants=[("seedlessRefit", "1100")],
      seedlessRefitConfig={"layerOrder": "outsidein"})
  check(False, "bad layerOrder must raise")
except ValueError:
  check(True, "bad layerOrder raises ValueError at config time")
# Refit-quality BDT model (spec §6a): default empty (no scoring), override flows.
check(_m.seedlessRefitBdtModel.value() == SEEDLESSREFIT_DEFAULTS["bdtModel"] == "",
      "bdtModel default is empty (no BDT scoring) and carried onto the module")
_m_bdt = getattr(
    addSmartPixelsTrackProducerVariants(
        cms.Process("TEST"), variants=[("seedlessRefit", "1100")],
        seedlessRefitConfig={"pixelavAngleSet": "dummy/x.json",
                         "bdtModel": "L1Trigger/Phase3SmartPixels/data/refitq_model.json"})[0],
    _dr_prompt)
check(_m_bdt.seedlessRefitBdtModel.value() == "L1Trigger/Phase3SmartPixels/data/refitq_model.json",
      "bdtModel override flows onto the module")
# disc surfaces: per-disc windows and the disc scope reach the module, bad windows raise
check(tuple(_m.seedlessRefitWindowDiscRPhi) == SEEDLESSREFIT_DEFAULTS["windowDiscRPhi"]
      and tuple(_m.seedlessRefitWindowDiscR) == SEEDLESSREFIT_DEFAULTS["windowDiscR"]
      and _m.seedlessRefitDiscScope.value() == SEEDLESSREFIT_DEFAULTS["discScope"] == "sp",
      "disc windows and discScope (sp: the refit visits ID1-ID4) are carried onto the module")
check(_m.seedlessRefitMaxKFUpdates.value() == 8, "maxKFUpdates default 8 (every IT surface) is carried onto the module")
expect_raises(ValueError, lambda: _resolveSeedlessRefitConfig({"windowDiscR": (0.3, 0.3)}),
              "a disc window with the wrong length raises ValueError", contains="windowDiscR")
check(_resolveSeedlessRefitConfig({"windowDiscRPhi": 0.1})["windowDiscRPhi"] == (0.1,) * 4,
      "a scalar disc window is broadcast to ID1-ID4")
# smarthitTrueSet is RESERVED: NOT consumed by Tier-2 but kept as a
# validated key for a future ASIC-efficiency model. Shipped/default recipes must
# pass NO real path (empty), and setting one must still validate + wire cleanly
# (the producer will warn at construction -- cmsRun-gated, DEFERRED).
check("smarthitTrueSet" in SEEDLESSREFIT_DEFAULTS and SEEDLESSREFIT_DEFAULTS["smarthitTrueSet"] == "",
      "smarthitTrueSet is a RESERVED key defaulting to empty (no shipped path)")
check(not hasattr(_m, "seedlessRefitSmarthitTrueSet"),
      "smarthitTrueSet no longer reaches the refit module: the payload surface moved upstream "
      "to SmartPixelsRecHitProducer with the angle synthesis")
# No RandomNumberGeneratorService: the refit is now engine-free entirely (its old
# local engine existed only for the noise-angle draw, which moved upstream with
# the synthesis).
check(not hasattr(p_dr, "RandomNumberGeneratorService"),
      "no RandomNumberGeneratorService is added (RNG is per-event, producer-side)")
# Determinism at config level: two identical builds produce identical module
# configs (the dump strings match byte-for-byte).
p_dr2 = cms.Process("TEST")
p_dr2, _ = addSmartPixelsTrackProducerVariants(p_dr2, variants=[("seedlessRefit", "1100")],
                                               seedlessRefitConfig=_DR_CFG)
check(getattr(p_dr2, _dr_prompt).dumpPython() == _m.dumpPython(),
      "two identical seedlessRefit builds produce identical module configs")

# a seedlessRefit variant WITHOUT a payload must fail loudly at config time
def _build_seedlessrefit_nopayload():
    p = cms.Process("TEST")
    addSmartPixelsTrackProducerVariants(p, variants=[("seedlessRefit", "1100")])
expect_raises(ValueError, _build_seedlessrefit_nopayload,
              "seedlessRefit without pixelavAngleSet raises ValueError",
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

# --- (d) invalid seedlessRefitConfig -------------------------------------------
print("[d] seedlessRefitConfig validation")
expect_raises(ValueError,
              lambda: _resolveSeedlessRefitConfig({"notAKey": 1}),
              "unknown seedlessRefitConfig key raises ValueError",
              contains="unknown key")
expect_raises(ValueError,
              lambda: _resolveSeedlessRefitConfig({"useAngles": "bogus"}),
              "invalid useAngles enum raises ValueError",
              contains="useAngles")
expect_raises(ValueError,
              lambda: _resolveSeedlessRefitConfig({"gainMode": "bogus"}),
              "invalid gainMode enum raises ValueError",
              contains="gainMode")
expect_raises(ValueError,
              lambda: _resolveSeedlessRefitConfig(["not", "a", "dict"]),
              "non-dict seedlessRefitConfig raises ValueError",
              contains="must be a dict")
# KF numerical guards must be positive doubles (spec §6b).
expect_raises(ValueError,
              lambda: _resolveSeedlessRefitConfig({"jacobianMaxAbs": -1.0}),
              "non-positive jacobianMaxAbs raises ValueError",
              contains="jacobianMaxAbs")
expect_raises(ValueError,
              lambda: _resolveSeedlessRefitConfig({"chi2UpdateGate": 0}),
              "zero chi2UpdateGate raises ValueError",
              contains="chi2UpdateGate")
expect_raises(ValueError,
              lambda: _resolveSeedlessRefitConfig({"chi2UpdateGate": "big"}),
              "non-numeric chi2UpdateGate raises ValueError",
              contains="chi2UpdateGate")
# Grazing-angle clamps must be positive doubles (spec §6b).
expect_raises(ValueError,
              lambda: _resolveSeedlessRefitConfig({"measAngleMaxAbs": -1.0}),
              "non-positive measAngleMaxAbs raises ValueError",
              contains="measAngleMaxAbs")
expect_raises(ValueError,
              lambda: _resolveSeedlessRefitConfig({"predAngleMaxAbs": 0}),
              "zero predAngleMaxAbs raises ValueError",
              contains="predAngleMaxAbs")
# Refit-quality BDT model must be a string path (spec §6a).
expect_raises(ValueError,
              lambda: _resolveSeedlessRefitConfig({"bdtModel": 17}),
              "non-string bdtModel raises ValueError",
              contains="bdtModel")

# valid config merges over defaults and does NOT mutate the defaults
merged = _resolveSeedlessRefitConfig({"minHits": 3, "useAngles": "alpha"})
check(merged["minHits"] == 3 and merged["useAngles"] == "alpha",
      "valid overrides applied")
check(merged["windowRPhi"] == SEEDLESSREFIT_DEFAULTS["windowRPhi"],
      "untouched keys keep defaults")
check(SEEDLESSREFIT_DEFAULTS["minHits"] != 3,
      "SEEDLESSREFIT_DEFAULTS not mutated by merge")
check(set(merged) == set(SEEDLESSREFIT_DEFAULTS) | {"v3"},
      "merged config has exactly the default key set (+ the v3 element switches)")

# an invalid seedlessRefitConfig must fail loudly even before the seedlessRefit
# variant hits its (Phase-0) NotImplementedError
def _coexist_bad_cfg():
    p = cms.Process("TEST")
    smartPixelsCoexist(p, variants=[("seedlessRefit", "1100")],
                       seedlessRefitConfig={"notAKey": 1})
expect_raises(ValueError, _coexist_bad_cfg,
              "invalid seedlessRefitConfig raises ValueError via coexist (before NotImplemented)",
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

# --- (f) trackInputMode="rebuildTracksFromStubs" wiring ---------------------
print("[f] trackInputMode='rebuildTracksFromStubs' wiring")
check("rebuildTracksFromStubs" in TRACKINPUTMODE_CHOICES,
      "rebuildTracksFromStubs is in the trackInputMode vocabulary")
check(set(TRACKINPUTMODE_CHOICES) == {"reemulateL1TrackFinding", "useStoredTracks",
                                      "rebuildTracksFromStubs"},
      "trackInputMode vocabulary is exactly the three renamed modes")


def _build_fromfilestubs(extendedTracks):
    """Full coexist build with trackInputMode=rebuildTracksFromStubs on a process
    carrying a dummy Path (so the stub-rebuild Task has somewhere to associate)."""
    p = cms.Process("TEST")
    p.pdummy = cms.Path()  # somewhere for the stub-rebuild Task to associate
    return smartPixelsCoexist(
        p, variants=[("seedlessRefit", "1100")], addNanoTables=False,
        trackInputMode="rebuildTracksFromStubs", extendedTracks=extendedTracks,
        seedlessRefitConfig={"pixelavAngleSet": "dummy/x.json"})


# extended ON (default): full prompt+extended stub-rebuild chain present
p_ffs = _build_fromfilestubs(extendedTracks=True)
_PROMPT_CHAIN = ["ProducerDTC", "l1tTTTracksFromTrackletEmulation",
                 "TTTrackAssociatorFromPixelDigis"]
_EXT_CHAIN = ["l1tTTTracksFromExtendedTrackletEmulation",
              "TTTrackAssociatorFromPixelDigisExtended"]
for m in _PROMPT_CHAIN:
    check(hasattr(p_ffs, m), f"stub-rebuild chain module '{m}' present (default labels)")
for m in _EXT_CHAIN:
    check(hasattr(p_ffs, m), f"stub-rebuild extended chain module '{m}' present (extendedTracks=True)")
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
          f"DIGI-tier/cluster-stub module '{m}' NOT scheduled in-process (rebuildTracksFromStubs)")
# the stub-rebuild Task exists and contains exactly the chain modules
check(hasattr(p_ffs, "l1tSmartPixelsStubRebuildTask"),
      "stub-rebuild Task l1tSmartPixelsStubRebuildTask exists")
_task_names = {mod.label_() for mod in p_ffs.l1tSmartPixelsStubRebuildTask._collection}
check(_task_names == set(_PROMPT_CHAIN + _EXT_CHAIN),
      f"stub-rebuild Task holds exactly prompt+extended chain (got {sorted(_task_names)})")

# extended OFF: prompt chain only. NOTE: the extended EMULATOR module
# (l1tTTTracksFromExtendedTrackletEmulation) is still defined on the process --
# the l1tTTTracksFromTrackletEmulation_cfi import brings both prompt and extended
# emulator objects regardless -- but it is NOT on the stub-rebuild Task and its
# associator is not cloned. The knob's load-bearing effect is Task membership.
p_ffs_p = _build_fromfilestubs(extendedTracks=False)
for m in _PROMPT_CHAIN:
    check(hasattr(p_ffs_p, m), f"prompt-only build still has '{m}'")
check(not hasattr(p_ffs_p, "TTTrackAssociatorFromPixelDigisExtended"),
      "extendedTracks=False -> extended track associator NOT cloned")
_task_names_p = {mod.label_() for mod in p_ffs_p.l1tSmartPixelsStubRebuildTask._collection}
check(_task_names_p == set(_PROMPT_CHAIN),
      f"prompt-only stub-rebuild Task holds exactly the prompt chain (got {sorted(_task_names_p)})")
check("l1tTTTracksFromExtendedTrackletEmulation" not in _task_names_p,
      "extendedTracks=False -> extended emulator NOT on the stub-rebuild Task")
# extended ON vs OFF changes the scheduled module set (knob is live)
_task_names_on = {mod.label_() for mod in p_ffs.l1tSmartPixelsStubRebuildTask._collection}
check(_task_names_on != _task_names_p and _EXT_CHAIN[0] in _task_names_on,
      "extendedTracks on/off changes the scheduled stub-rebuild module set")

# unknown trackInputMode still raises loudly
def _bad_trackinputmode():
    p = cms.Process("TEST")
    _applyTrackInputMode(p, "bogusMode", [("passthrough", None)])
expect_raises(ValueError, _bad_trackinputmode,
              "unknown trackInputMode raises ValueError", contains="trackInputMode must be one of")

# RETIRED pre-rename spellings now RAISE, and the error names the replacement.
# They were accepted with a warning for a transition cycle in which nothing
# migrated (106 of 222 archived configs used truthSource=, zero used
# trackInputMode=), so the alias was a second vocabulary rather than a path.
from L1Trigger.Phase3SmartPixels.customizeSmartPixels_cff import _resolveTrackInputMode
for _old, _new in (("inJob", "reemulateL1TrackFinding"),
                   ("redigitizePVignorePU", "reemulateL1TrackFinding"),
                   ("fromFile", "useStoredTracks"),
                   ("fromFileStubs", "rebuildTracksFromStubs")):
    expect_raises(ValueError, lambda o=_old: _resolveTrackInputMode(o),
                  f"retired trackInputMode spelling {_old!r} raises",
                  contains="RETIRED")
    expect_raises(ValueError, lambda o=_old: _resolveTrackInputMode(o),
                  f"the error for {_old!r} names {_new!r} as the replacement",
                  contains=_new)
# the truthSource= keyword is gone from the signatures entirely -> TypeError
expect_raises(TypeError,
              lambda: smartPixelsCoexist(cms.Process("TEST"), truthSource="fromFile"),
              "retired truthSource= keyword is a TypeError on smartPixelsCoexist",
              contains="truthSource")
expect_raises(TypeError,
              lambda: smartPixelsCoopt(cms.Process("TEST"), truthSource="fromFile"),
              "retired truthSource= keyword is a TypeError on smartPixelsCoopt",
              contains="truthSource")

# PU-safety guard: reemulateL1TrackFinding + a scheduled digitisation step is
# refused (signal-only re-digitization would destroy pileup and its truth).
def _reemulate_with_digi_step():
    p = cms.Process("TEST")
    p.digitisation_step = cms.Path()
    _applyTrackInputMode(p, "reemulateL1TrackFinding", [("passthrough", None)])
expect_raises(RuntimeError, _reemulate_with_digi_step,
              "reemulateL1TrackFinding + digitisation step raises (PU guard)",
              contains="DESTROYS pileup")

# attachStubRebuildChain returns the chain module list (prompt-only when off)
p_direct = cms.Process("TEST")
p_direct.pdummy = cms.Path()
_p, _chain = attachStubRebuildChain(p_direct, extendedTracks=False)
check(_chain == _PROMPT_CHAIN,
      f"attachStubRebuildChain(extendedTracks=False) returns prompt chain (got {_chain})")

# --- (g) promptHnpar (PRIME target: 5-par prompt OT track) ------------------
print("[g] promptHnpar wiring (5-par prompt OT track -> real 5x5 seedlessRefit starting cov)")
check(set(PROMPT_HNPAR_CHOICES) == {4, 5},
      "promptHnpar vocabulary is exactly {4, 5}")


def _build_ffs_hnpar(promptHnpar):
    """coexist + rebuildTracksFromStubs seedlessRefit build with a given promptHnpar."""
    p = cms.Process("TEST")
    p.pdummy = cms.Path()
    return smartPixelsCoexist(
        p, variants=[("seedlessRefit", "1100")], addNanoTables=False,
        trackInputMode="rebuildTracksFromStubs", extendedTracks=True, promptHnpar=promptHnpar,
        seedlessRefitConfig={"pixelavAngleSet": "dummy/x.json"})


# promptHnpar=5 -> PRIME target: prompt tracklet producer set to Hnpar=5, extended
# stays 5, and the prompt seedlessRefit variant uses otNPar=5.
p_h5 = _build_ffs_hnpar(5)
check(p_h5.l1tTTTracksFromTrackletEmulation.Hnpar.value() == 5,
      "promptHnpar=5 sets prompt l1tTTTracksFromTrackletEmulation.Hnpar = 5")
check(p_h5.l1tTTTracksFromTrackletEmulation.Extended.value() is False,
      "promptHnpar=5 keeps prompt Extended = False (5-par PROMPT, not displaced)")
check(p_h5.l1tTTTracksFromExtendedTrackletEmulation.Hnpar.value() == 5,
      "extended chain stays Hnpar = 5")
_dr5 = getattr(p_h5, "l1tSmartPixelsTrackProducerWseedlessRefitAAII")
check(_dr5.seedlessRefitOTNPar.value() == 5,
      "prompt seedlessRefit variant uses otNPar=5 (prime target)")

# promptHnpar=4 (ablation only, still selectable) -> prompt Hnpar=4, d0 pinned.
p_h4 = _build_ffs_hnpar(4)
check(p_h4.l1tTTTracksFromTrackletEmulation.Hnpar.value() == 4,
      "promptHnpar=4 (ablation) leaves prompt Hnpar = 4 (d0 pinned; weak-d0-prior refit fallback)")
check(p_h4.l1tTTTracksFromExtendedTrackletEmulation.Hnpar.value() == 5,
      "promptHnpar=4 still has extended chain at Hnpar = 5")
# 5-par-only framing: default coexist (no promptHnpar passed) == promptHnpar=5.
p_hdef = smartPixelsCoexist(
    (lambda: (lambda pp: (setattr(pp, "pdummy", cms.Path()) or pp))(cms.Process("TEST")))(),
    variants=[("seedlessRefit", "1100")], addNanoTables=False, trackInputMode="rebuildTracksFromStubs",
    seedlessRefitConfig={"pixelavAngleSet": "dummy/x.json"})
check(p_hdef.l1tTTTracksFromTrackletEmulation.Hnpar.value() == 5,
      "default promptHnpar (unspecified) == 5 (5-par-only framing: real prompt d0 + 5x5 cov)")

# out-of-vocabulary promptHnpar raises loudly
expect_raises(ValueError, lambda: attachStubRebuildChain(
    (lambda pp: (setattr(pp, "pdummy", cms.Path()) or pp))(cms.Process("TEST")), promptHnpar=6),
    "promptHnpar=6 raises ValueError", contains="promptHnpar must be one of")

# --- (h) OT track -> stub link tables ---------------------------------------
print("[h] OT track -> stub link tables")
from DPGAnalysis.Phase3SmartPixelsNanoAOD.l1tPh3SmartPixelsNano_cff import (
    addPh3L1SmartPixelsStubTables, addPh3L1SmartPixelsTracks)


def _proc_with_nanotask():
    p = cms.Process("TEST")
    p.l1tPh2NanoTask = cms.Task()
    return p


# prompt + extended link tables, keyed to the REFERENCE collections (not variants)
p_stub = _proc_with_nanotask()
addPh3L1SmartPixelsStubTables(p_stub, doExtended=True)
check(hasattr(p_stub, "l1tPh3TrackStubLinkTable") and hasattr(p_stub, "l1tPh3ExtTrackStubLinkTable"),
      "prompt and extended track -> stub link modules exist (labels end in 'Table')")
check(p_stub.l1tPh3TrackStubLinkTable.tracks
      == cms.InputTag("l1tTTTracksFromTrackletEmulation", "Level1TTTracks")
      and p_stub.l1tPh3ExtTrackStubLinkTable.tracks
      == cms.InputTag("l1tTTTracksFromExtendedTrackletEmulation", "Level1TTTracks"),
      "link tables keyed to the reference tracklet collections")
check(p_stub.l1tPh3TrackStubLinkTable.trackTableName.value() == "L1TTrack"
      and p_stub.l1tPh3ExtTrackStubLinkTable.trackTableName.value() == "L1TExtTrack",
      "trackIdx links to L1TTrack / L1TExtTrack")
check(p_stub.l1tPh3TrackStubLinkTable.linkTableName.value() == "L1TTrackStub"
      and p_stub.l1tPh3ExtTrackStubLinkTable.linkTableName.value() == "L1TExtTrackStub",
      "link table names are L1TTrackStub / L1TExtTrackStub")
_nano_mods = set(p_stub.l1tPh2NanoTask.moduleNames())
check({"l1tPh3TrackStubLinkTable", "l1tPh3ExtTrackStubLinkTable", "l1tPh3SmartPixelsAllStubTable"} <= _nano_mods,
      "both link tables and L1TOTStub (their otStubIdx target) added to l1tPh2NanoTask")
_n_before = len(p_stub.l1tPh2NanoTask.moduleNames())
addPh3L1SmartPixelsStubTables(p_stub, doExtended=True)
check(len(p_stub.l1tPh2NanoTask.moduleNames()) == _n_before,
      "addPh3L1SmartPixelsStubTables is idempotent (no duplicate modules)")

p_stub_p = _proc_with_nanotask()
addPh3L1SmartPixelsStubTables(p_stub_p, doExtended=False)
check(hasattr(p_stub_p, "l1tPh3TrackStubLinkTable") and not hasattr(p_stub_p, "l1tPh3ExtTrackStubLinkTable"),
      "doExtended=False adds the prompt link table only")

# coexist(addNanoTables=True) wires the link tables alongside the track tables
p_cox = cms.Process("TEST")
p_cox.pdummy = cms.Path()
p_cox.l1tPh2NanoTask = cms.Task()
smartPixelsCoexist(p_cox, variants=[("passthrough", None)], addNanoTables=True,
                   trackInputMode="reemulateL1TrackFinding")
check(hasattr(p_cox, "l1tPh3TrackStubLinkTable") and hasattr(p_cox, "l1tPh3SmartPixelsAllStubTable"),
      "smartPixelsCoexist(addNanoTables=True) emits the track -> stub links and L1TOTStub")
_cox_tables = {m.name.value() for m in p_cox.producers_().values()
               if isinstance(getattr(m, "name", None), cms.string)}
check("L1TSmartPixelsTrackPassthrough" in _cox_tables
      and not _cox_tables & {"L1TSmartPixelsTrack", "L1TSmartPixelsExtTrack"},
      "smartPixelsCoexist writes suffixed variant tables only, never the canonical unsuffixed ones")
check(not hasattr(p_cox, "l1tPh3SmartPixelsOTTrackLinkTable"),
      "smartPixelsCoexist adds no unsuffixed OT-track link table")

# coexist with the canonical SmartPixels tables (@SPTrk) raises
p_coxtier = cms.Process("TEST")
p_coxtier.pdummy = cms.Path()
p_coxtier.l1tPh2NanoTask = cms.Task()
addPh3L1SmartPixelsTracks(p_coxtier)
expect_raises(RuntimeError, lambda: smartPixelsCoexist(p_coxtier, variants=[("passthrough", None)]),
              "smartPixelsCoexist refuses a job with the canonical SmartPixels tables",
              contains="canonical SmartPixels track tables")

# --- (i) unified cluster table, disc scope, stub index, TP table, cut guard ---
print("[i] unified cluster table / discScope / otStubIdx / L1TTP / track-table cut guard")
from DPGAnalysis.Phase3SmartPixelsNanoAOD.l1tPh3SmartPixelsNano_cff import (
    addPh3L1SmartPixelsClusters, assertTrackTablesUncut)
import DPGAnalysis.Phase3SmartPixelsNanoAOD.l1tPh3SmartPixelsNano_cff as _p3nano
from DPGAnalysis.Phase2L1TNanoAOD.l1tPh2TrkNanotables_cff import l1tTracksTable, l1tTrackingParticleTable

check(SEEDLESSREFIT_DEFAULTS["discScope"] == "sp", "default discScope is sp (ID1-ID4)")
expect_raises(ValueError, lambda: _resolveSeedlessRefitConfig({"discScope": "tepx"}),
              "invalid discScope raises ValueError", contains="discScope")
check(not hasattr(_p3nano, "addPh3L1SmartPixelsDiscClusters")
      and not hasattr(_p3nano, "l1tPh3SmartPixelsDiscClusterTable"),
      "the separate disc cluster table is gone")


def _emu_refit(scope, tableFirst):
    p = cms.Process("TEST")
    p.pdummy = cms.Path()
    p.l1tPh2NanoTask = cms.Task()
    if tableFirst:
        addPh3L1SmartPixelsClusters(p)
    cfg = {"angleSource": "emu"}
    if scope is not None:
        cfg["discScope"] = scope
    smartPixelsCoexist(p, variants=[("seedlessRefit", "1111")], addNanoTables=True,
                       trackInputMode="reemulateL1TrackFinding", seedlessRefitConfig=cfg)
    if not tableFirst:
        addPh3L1SmartPixelsClusters(p)
    return p


for scope in (None, "none"):
    want = scope or "sp"
    for tableFirst in (True, False):
        pc = _emu_refit(scope, tableFirst)
        check(pc.spixSmartPixelsRecHits.discScope.value() == want
              and pc.l1tPh3SmartPixelsClusterTable.discScope.value() == want,
              f"rec hits and cluster table share discScope={want} (table added "
              f"{'before' if tableFirst else 'after'} the refit)")
        check(not hasattr(pc.spixSmartPixelsRecHits, "includeDiscs"),
              "rec hits have no includeDiscs parameter")
        check(pc.spixSmartPixelsRecHits.angleSource.value() == "emu"
              and pc.spixSmartPixelsRecHits.emuLibrary.value().endswith("spix_emu_library_v2.root")
              and pc.spixSmartPixelsRecHits.emuBarrelA1.value().endswith("spix_emu_a1_barrel_v2.onnx")
              and pc.spixSmartPixelsRecHits.emuDiscA1.value().endswith("spix_emu_a1_disc_v2.onnx")
              and len(pc.spixSmartPixelsRecHits.emuBarrelScalars) == 23
              and len(pc.spixSmartPixelsRecHits.emuDiscScalars) == 21,
              "emu source wired with its models, scalars and library")
        check(not hasattr(pc.spixSmartPixelsRecHits, "emuTruthPtMin")
              and not hasattr(pc.spixSmartPixelsRecHits, "emuDonorK"),
              "emu class T has no pT cut and no donor library")
expect_raises(ValueError, lambda: _resolveSeedlessRefitConfig({"angleSource": "emu", "discScope": "all"}),
              "emu with discScope all raises ValueError", contains="discScope")
expect_raises(ValueError, lambda: _resolveSeedlessRefitConfig({"angleSource": "emu", "pixelavAngleSet": "x.json"}),
              "emu with a pixelavAngleSet raises ValueError", contains="pixelavAngleSet")
expect_raises(ValueError, lambda: _resolveSeedlessRefitConfig({"angleSource": "cnn"}),
              "the removed cnn angle source raises ValueError", contains="angleSource")
check(not hasattr(_p3nano.l1tPh3SmartPixelsClusterTable, "subdetector"),
      "cluster table has no subdetector option")
check(not hasattr(_p3nano.l1tPh3SmartPixelsAllStubTable, "barrelOnly"),
      "L1TOTStub is unfiltered (row == stub key), no barrelOnly option")
check(p_stub.l1tPh3TrackStubLinkTable.stubs
      == cms.InputTag("TTStubsFromPhase2TrackerDigis", "StubAccepted")
      and p_stub.l1tPh3TrackStubLinkTable.stubs == _p3nano.l1tPh3SmartPixelsAllStubTable.stubs,
      "track -> stub link tables index the same stub collection L1TOTStub is built from")
check(hasattr(l1tTrackingParticleTable, "minP") and not hasattr(l1tTrackingParticleTable, "minPt")
      and l1tTrackingParticleTable.minP.value() == 1.0,
      "L1TTP selects on 3-momentum |p| > 1 GeV")


def _cut_track_table():
    p = cms.Process("TEST")
    p.l1tTracksTable = l1tTracksTable.clone(cut=cms.string("pt > 2"))
    assertTrackTablesUncut(p)
expect_raises(RuntimeError, _cut_track_table, "a cut on the L1TTrack table is refused",
              contains="cut")
p_ok = cms.Process("TEST")
p_ok.l1tTracksTable = l1tTracksTable.clone()
assertTrackTablesUncut(p_ok)
check(True, "an uncut L1TTrack table passes the guard")

# --- (j) nano add-on flavours: registry and the @PHYS input guard -----------
print("[j] nano add-on flavours")
from PhysicsTools.NanoAOD.autoNANO import autoNANO
import DPGAnalysis.Phase3SmartPixelsNanoAOD.l1tPh3SmartPixelsNanoFlavours_cff as _flav

check(sum(k.startswith("SPActiveLayers") for k in autoNANO) == 16
      and all(hasattr(_flav, autoNANO[k]["customize"].split(".")[1]) for k in autoNANO
              if autoNANO[k].get("customize", "").startswith("DPGAnalysis/Phase3SmartPixelsNanoAOD/")),
      "16 mask tags registered, every SmartPixels flavour resolves to a customisation")
check(not any("SmartPix" in k for k in autoNANO), "the combined SmartPix flavours are gone")
expect_raises(NotImplementedError, lambda: _flav.useSPJointSeedFit(cms.Process("TEST")),
              "@SPJointSeedFit is reserved", contains="RESERVED")


def _phys_job(recoInJob):
    p = cms.Process("NANO")
    p.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring())
    p.tab = cms.EDProducer("X")
    p.nanoSequenceMC = cms.Sequence(p.tab)
    p.nanoAOD_step = cms.Path(p.nanoSequenceMC)
    p.l1tPh2NanoTask = cms.Task()
    if recoInJob:
        p.reco = cms.EDProducer("Y")
        p.reconstruction_step = cms.Path(p.reco)
        p.patTask = cms.Task(p.reco)
        p.schedule = cms.Schedule(p.reconstruction_step, p.nanoAOD_step)
        p.schedule.associate(p.patTask)
    return p


expect_raises(RuntimeError, lambda: _flav._checkRecoInputs(_phys_job(False), _flav._state(_phys_job(False))),
              "@PHYS without RECO+PAT and without an input to check is refused", contains="MINIAOD")
_flav._checkRecoInputs(_phys_job(True), _flav._state(_phys_job(True)))
check(True, "@PHYS with RECO+PAT in the job passes the guard")

# coexist variants get the canonical tables' links, whichever of coexist and the tables comes first
print("[k] coexist variant track -> stub / cluster links, either order")


def _coexist_links(tablesFirst):
    p = cms.Process("TEST")
    p.pdummy = cms.Path()
    p.l1tDummyTable = cms.EDProducer("X")
    p.p2L1TablesTask = cms.Task(p.l1tDummyTable)
    p.l1tPh2NanoTask = cms.Task(p.p2L1TablesTask)
    if tablesFirst:
        _flav.addSPClusterStubTables(p)
    smartPixelsCoexist(p, variants=[("seedlessRefit", "1111"), ("seedlessRefit", "1100")], addNanoTables=True,
                       trackInputMode="rebuildTracksFromStubs", seedlessRefitConfig={"angleSource": "emu"})
    if not tablesFirst:
        _flav.addSPClusterStubTables(p)
    return p


_links = {}
for tablesFirst in (True, False):
    pl = _coexist_links(tablesFirst)
    for sfx in ("SeedlessRefitAAAA", "SeedlessRefitAAII"):
        for ext in ("", "Ext"):
            for kind, target in (("Stub", "L1TOTStub"), ("Cluster", "L1TSmartPixelsCluster")):
                label = f"l1tPh3{ext}SmartPixelsTrack{kind}Link{sfx}Table"
                mod = getattr(pl, label, None)
                track = f"L1TSmartPixels{ext}Track{sfx}"
                check(mod is not None and label in pl.l1tPh2NanoTask.moduleNames()
                      and mod.linkTableName.value() == f"{track}{kind}" and mod.trackTableName.value() == track
                      and mod.tracks == cms.InputTag(smartPixelsVariantLabels(
                          "seedlessRefit", "1111" if sfx.endswith("AAAA") else "1100")[ext == "Ext"], "Level1TTTracks"),
                      f"{track}{kind} -> {target} (tables {'before' if tablesFirst else 'after'} coexist)")
                if mod is not None:
                    _links.setdefault(label, set()).add(mod.dumpPython())
    check(pl.smartPixelsNano.trackInputModeApplied.value() == "rebuildTracksFromStubs",
          "coexist records its trackInputMode for later reconciles")
check(all(len(v) == 1 for v in _links.values()), "the link modules do not depend on the order")

print("[l] refit v3 preset (seedlessRefitV3)")
from L1Trigger.Phase3SmartPixels.customizeSmartPixels_cff import SEEDLESSREFIT_V3_DEFAULTS, SEEDLESSREFIT_PRESETS
check(smartPixelsVariantSuffix("seedlessRefitV3", "1111") == "seedlessRefitV3AAAA", "v3 suffix")
check("seedlessRefitV3" in TRUTH_REQUIRED_MODES and SEEDLESSREFIT_PRESETS["seedlessRefitV3"] == "v3",
      "v3 is truth-required and maps to preset v3")
_r3 = _resolveSeedlessRefitConfig({"v3": {"hitChoice": "best", "angleSteps": "first"}})
check(_r3["v3"]["hitChoice"] == "best" and _r3["v3"]["angleSteps"] == "first"
      and _r3["v3"]["nsig"] == SEEDLESSREFIT_V3_DEFAULTS["nsig"], "v3 overrides merge over the v3 defaults")
expect_raises(ValueError, lambda: _resolveSeedlessRefitConfig({"v3": {"nope": 1}}), "unknown v3 key raises")
expect_raises(ValueError, lambda: _resolveSeedlessRefitConfig({"v3": {"angleUse": ["window:alpha:barrel"]}}),
              "a not-implemented angle consumer raises")
expect_raises(ValueError, lambda: _resolveSeedlessRefitConfig({"v3": {"scattering": "producer"}}),
              "a not-implemented scattering model raises")
_p3 = cms.Process("TESTV3")
_p3, _mods3 = addSmartPixelsTrackProducerVariants(_p3, variants=[("seedlessRefitV3", "1111"), ("seedlessRefit", "1111")],
                                                   seedlessRefitConfig={"angleSource": "emu"})
_m3 = getattr(_p3, smartPixelsVariantLabels("seedlessRefitV3", "1111")[0])
_m2 = getattr(_p3, smartPixelsVariantLabels("seedlessRefit", "1111")[0])
check(_m3.smartPixelsEmulatorMode.value() == "seedlessRefit" and _m3.seedlessRefitPreset.value() == "v3",
      "the v3 producer runs mode seedlessRefit with preset v3")
check(_m2.seedlessRefitPreset.value() == "v2", "the v2 producer keeps preset v2")
check(list(_m3.seedlessRefitV3.angleUse.value()) == list(SEEDLESSREFIT_V3_DEFAULTS["angleUse"])
      and _m3.seedlessRefitV3.knownBad.value() == "" and _m3.seedlessRefitV3.closureDumpFile.value() == "",
      "the v3 PSet carries the defaults, no known-bad, no dump")

# --- summary ---------------------------------------------------------------
print()
if _failures:
    print(f"RESULT: FAILED ({len(_failures)} failing check(s))")
    for f in _failures:
        print("  - " + f)
    sys.exit(1)
print("RESULT: PASSED (all checks green)")
sys.exit(0)
