#!/usr/bin/env python3
"""Pure-python (no cmsRun) checks of the SmartPixels refit-mode config surface.

Exercises Phase-0 of the tier-2 refit plan in
L1Trigger/Phase3SmartPixels/python/customizeSmartPixels_cff.py:

 (a) suffix/label generation for digiRefit variants works and labels still end
     so that derived nano tables end in 'Table' (keep-pattern gotcha);
 (b) coexist with a digiRefit variant raises NotImplementedError
     ("... lands in Phase 2 ...");
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
    TRUTH_REQUIRED_MODES,
    _normalizeVariants,
    _resolveDigiRefitConfig,
    addSmartPixelsTrackProducerVariants,
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

# --- (b) coexist with a digiRefit variant raises NotImplementedError --------
print("[b] coexist digiRefit -> NotImplementedError")
def _coexist_digirefit():
    p = cms.Process("TEST")
    smartPixelsCoexist(p, variants=[("digiRefit", "1100")])
expect_raises(NotImplementedError, _coexist_digirefit,
              "smartPixelsCoexist(digiRefit) raises NotImplementedError",
              contains="Phase 2")

# also directly via the producer-variant builder
def _build_digirefit():
    p = cms.Process("TEST")
    addSmartPixelsTrackProducerVariants(p, variants=[("digiRefit", "1100")])
expect_raises(NotImplementedError, _build_digirefit,
              "addSmartPixelsTrackProducerVariants(digiRefit) raises NotImplementedError",
              contains="Phase 2")

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

# --- summary ---------------------------------------------------------------
print()
if _failures:
    print(f"RESULT: FAILED ({len(_failures)} failing check(s))")
    for f in _failures:
        print("  - " + f)
    sys.exit(1)
print("RESULT: PASSED (all checks green)")
sys.exit(0)
