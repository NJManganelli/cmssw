#!/usr/bin/env python3
"""Configuration hygiene of the SmartPixels packages (scram unit test, testConfigHygiene):

1. No local-machine paths in any Python file under python/ or test/: no laptop mounts or home directories
   (/host_volumes, /Volumes/, /Users/, the drive names, the development container's /work/). Configurations and
   generators must run anywhere, the grid included; inputs come from arguments or placeholders.
2. No configuration-time file access in configuration code (python/*_cff.py, *_cfi.py, *_cfg.py, customise
   modules, and test/*_cfg.py): no uproot, no ROOT file opening. A configuration must load without its input:
   on the grid CRAB sets the input after loading the pset. Input kinds are stated explicitly (--campaign,
   inputHasMiniAOD, hasMiniAOD).
"""
import os
import re
import sys

PACKAGES = ("L1Trigger/Phase3SmartPixels", "DPGAnalysis/Phase3SmartPixelsNanoAOD", "DPGAnalysis/Phase2L1TNanoAOD",
            "DataFormats/Phase3SmartPixels")
LOCAL_PATHS = re.compile(r"/host_volumes|/Volumes/|/Users/|NJM256GBSD|WDMac|(?<![\w.])/work/")
FILE_ACCESS = re.compile(r"^\s*(import|from)\s+uproot\b|\buproot\.\w|TFile\.Open\s*\(|ROOT\.TFile\s*\(|inputEventBranches")
SELF = os.path.abspath(__file__)


def is_configuration(path):
    name = os.path.basename(path)
    return name.endswith(("_cff.py", "_cfi.py", "_cfg.py")) or name.startswith(("customize", "customise"))


def main():
    src = os.path.abspath(os.path.join(os.path.dirname(SELF), "..", "..", ".."))
    problems, n = [], 0
    for pkg in PACKAGES:
        for sub in ("python", "test"):
            root = os.path.join(src, pkg, sub)
            for dirpath, _, files in os.walk(root):
                for f in files:
                    path = os.path.join(dirpath, f)
                    if not f.endswith(".py") or os.path.abspath(path) == SELF:
                        continue
                    n += 1
                    conf = is_configuration(path)
                    for i, line in enumerate(open(path, encoding="utf-8", errors="replace"), 1):
                        if LOCAL_PATHS.search(line):
                            problems.append(f"{os.path.relpath(path, src)}:{i}: local path: {line.strip()[:120]}")
                        if conf and FILE_ACCESS.search(line):
                            problems.append(f"{os.path.relpath(path, src)}:{i}: configuration-time file access: "
                                            f"{line.strip()[:120]}")
    if problems:
        print("\n".join(problems))
        print(f"testConfigHygiene: FAILED ({len(problems)} problems in {n} files)")
        return 1
    print(f"testConfigHygiene: passed ({n} Python files)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
