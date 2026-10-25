#!/usr/bin/env python3
"""makeCrabTaggerTier.campaign_of (scram unit test, testCrabDatasetCampaign): the campaign of a CRAB task from the
dataset NAME alone (no file access, no DAS query).

The Phase2Spring24 name below is an INFERRED form of the TTbar sample (from the local file's name), not confirmed in
DAS; the RelVal name is the physicist's (2026-10-25).
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from makeCrabTaggerTier import campaign_of  # noqa: E402

RELVAL = ("/RelValTTbar_14TeV/CMSSW_20_1_0_pre2-PU_150X_mcRun4_realistic_v1_STD_D121_RegeneratedGS_PU_20260917_163358_"
          "RV24-v1/GEN-SIM-DIGI-RAW")
SPRING24 = ("/TT_TuneCP5_14TeV-powheg-pythia8/Phase2Spring24DIGIRECOMiniAOD-PU200_Trk1GeV_140X_mcRun4_realistic_v4-v2/"
            "GEN-SIM-DIGI-RAW-MINIAOD")
failures = []


def expect(dataset, campaign, want):
    try:
        got = campaign_of(dataset, campaign)
    except SystemExit as e:
        got = SystemExit
        msg = str(e)
    ok = got == want
    print(f"[{'ok' if ok else 'FAIL'}] campaign_of({dataset!r}, {campaign!r}) -> "
          f"{'refused: ' + msg if got is SystemExit else got}")
    if not ok:
        failures.append((dataset, campaign))


# inference
expect(RELVAL, None, "d121")
expect(SPRING24, None, "spring24")
# an explicit campaign that agrees, and one that overrides an unknown name (checked against the tier)
expect(RELVAL, "d121", "d121")
expect("/X/SomeOtherProduction-v1/GEN-SIM-DIGI-RAW", "d121", "d121")
expect("/X/SomeOtherProduction-v1/GEN-SIM-DIGI-RAW-MINIAOD", "spring24", "spring24")
# unknown names without --campaign
expect("/X/SomeOtherProduction-v1/GEN-SIM-DIGI-RAW", None, SystemExit)
# contradictions between the campaign and the tier
expect("/TT/Phase2Spring24DIGIRECOMiniAOD-PU200-v2/GEN-SIM-DIGI-RAW", None, SystemExit)
expect(RELVAL.replace("/GEN-SIM-DIGI-RAW", "/GEN-SIM-DIGI-RAW-MINIAOD"), None, SystemExit)
expect(RELVAL, "spring24", SystemExit)
expect(SPRING24, "d121", SystemExit)
# empty or malformed dataset paths
expect("", None, SystemExit)
expect("", "d121", "d121")
for bad in ("/A/B", "A/B/C", "/A//C", "/A/B/C/D", "/A/B/"):
    expect(bad, None, SystemExit)

print(f"testCrabDatasetCampaign: {'FAILED' if failures else 'passed'}")
sys.exit(1 if failures else 0)
