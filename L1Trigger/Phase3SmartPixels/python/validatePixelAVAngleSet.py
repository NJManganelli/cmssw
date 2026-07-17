#!/usr/bin/env python3
"""Validator + reference example for the PixelAV angle-response payload.

Spec: L1Trigger/Phase3SmartPixels/doc/PixelAVAngleResponseSpec.md

  validatePixelAVAngleSet.py PAYLOAD.json            # validate a candidate payload
  validatePixelAVAngleSet.py --write-example X.json  # emit + self-validate the example
"""
import argparse
import gzip
import json
import math
import sys

CORRECTIONS = {
    "spx_angle_alpha_sigma": dict(kind="sigma"),
    "spx_angle_alpha_bias": dict(kind="bias"),
    "spx_angle_beta_sigma": dict(kind="sigma"),
    "spx_angle_beta_bias": dict(kind="bias"),
    "spx_angle_valid_prob": dict(kind="prob"),
}
INPUTS = [("layer", "int"), ("cotAlpha", "real"), ("cotBeta", "real"), ("bLocalY", "real")]
LAYERS = [1, 2, 3, 4]
COTALPHA_GRID = [round(-0.6 + 0.1 * i, 3) for i in range(13)]
COTBETA_GRID = [round(-6.0 + 0.5 * i, 3) for i in range(25)]
BLOCALY_GRID = [-3.81, 3.81]
# out-of-domain probes: must CLAMP, not raise
CLAMP_PROBES = [
    (1, 0.0, 10.0, -3.81),
    (1, 0.0, -10.0, -3.81),
    (4, 0.9, 0.0, -3.81),
    (2, 0.0, 0.0, -5.0),
    (2, 0.0, 0.0, 5.0),
]


def _example_multibinning(kind):
    """Coarse but structurally complete MultiBinning over (cotAlpha, cotBeta, bLocalY)."""
    ca_edges = [-0.6, -0.2, 0.0, 0.2, 0.6]
    cb_edges = [-6.0, -3.0, -1.0, 0.0, 1.0, 3.0, 6.0]
    by_edges = [-4.0, 0.0, 4.0]
    ncell = (len(ca_edges) - 1) * (len(cb_edges) - 1) * (len(by_edges) - 1)
    content = []
    for ica in range(len(ca_edges) - 1):
        for icb in range(len(cb_edges) - 1):
            cb_mid = 0.5 * (cb_edges[icb] + cb_edges[icb + 1])
            for iby in range(len(by_edges) - 1):
                if kind == "sigma":
                    # toy: resolution grows with |cotBeta| (longer clusters)
                    content.append(round(0.02 + 0.01 * abs(cb_mid), 5))
                elif kind == "bias":
                    content.append(round(0.001 * cb_mid, 6))
                else:  # prob
                    content.append(round(max(0.5, 0.98 - 0.05 * abs(cb_mid)), 4))
    assert len(content) == ncell
    return {
        "nodetype": "multibinning",
        "inputs": ["cotAlpha", "cotBeta", "bLocalY"],
        "edges": [ca_edges, cb_edges, by_edges],
        "content": content,
        "flow": "clamp",
    }


def write_example(path):
    corrections = []
    for name, meta in CORRECTIONS.items():
        corrections.append(
            {
                "name": name,
                "description": (
                    f"EXAMPLE {meta['kind']} payload. Provenance: pixelav=EXAMPLE-vX, "
                    "nn=EXAMPLE-training-id, sensor=EXAMPLE-variant, bias/temp=EXAMPLE, "
                    "bGrid=[-3.81,+3.81]T, date=1970-01-01, fitrev=deadbeef"
                ),
                "version": 0,
                "inputs": [
                    {"name": "layer", "type": "int", "description": "TBPX layer 1-4"},
                    {"name": "cotAlpha", "type": "real", "description": "true local px/pz"},
                    {"name": "cotBeta", "type": "real", "description": "true local py/pz"},
                    {"name": "bLocalY", "type": "real", "description": "local-y B [T], signed"},
                ],
                "output": {"name": name, "type": "real"},
                "data": {
                    "nodetype": "category",
                    "input": "layer",
                    "content": [
                        {"key": l, "value": _example_multibinning(meta["kind"])} for l in LAYERS
                    ],
                },
            }
        )
    cset = {
        "schema_version": 2,
        "description": (
            "EXAMPLE PixelAV SmartPixels angle-response payload; see "
            "doc/PixelAVAngleResponseSpec.md. NOT physics — structural reference only."
        ),
        "corrections": corrections,
    }
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "wt") as f:
        json.dump(cset, f, indent=1)
    print(f"wrote example payload: {path}")


def validate(path):
    import correctionlib

    errors = []
    raw = json.load(gzip.open(path, "rt") if path.endswith(".gz") else open(path))
    if raw.get("schema_version") != 2:
        errors.append(f"schema_version={raw.get('schema_version')} != 2")
    if "EXAMPLE" not in raw.get("description", "") and not raw.get("description"):
        errors.append("CorrectionSet.description missing (provenance required)")
    cset = correctionlib.CorrectionSet.from_file(path)
    present = set(cset)
    missing = set(CORRECTIONS) - present
    if missing:
        errors.append(f"missing corrections: {sorted(missing)}")
    for cjson in raw.get("corrections", []):
        name = cjson.get("name")
        if name not in CORRECTIONS:
            continue
        got = [(i["name"], i["type"]) for i in cjson.get("inputs", [])]
        if got != INPUTS:
            errors.append(f"{name}: inputs {got} != required {INPUTS}")
        if not cjson.get("description"):
            errors.append(f"{name}: description/provenance missing")
    if errors:
        return errors

    for name, meta in CORRECTIONS.items():
        corr = cset[name]
        nchecked = 0
        for layer in LAYERS:
            for ca in COTALPHA_GRID:
                for cb in COTBETA_GRID:
                    for by in BLOCALY_GRID:
                        try:
                            v = corr.evaluate(layer, ca, cb, by)
                        except Exception as e:  # noqa: BLE001
                            errors.append(f"{name}({layer},{ca},{cb},{by}) raised: {e}")
                            return errors
                        if not math.isfinite(v):
                            errors.append(f"{name}({layer},{ca},{cb},{by}) not finite")
                        elif meta["kind"] == "sigma" and v <= 0:
                            errors.append(f"{name}({layer},{ca},{cb},{by})={v} (sigma<=0)")
                        elif meta["kind"] == "bias" and abs(v) >= 1.0:
                            errors.append(f"{name}({layer},{ca},{cb},{by})={v} (|bias|>=1)")
                        elif meta["kind"] == "prob" and not (0.0 <= v <= 1.0):
                            errors.append(f"{name}({layer},{ca},{cb},{by})={v} (not in [0,1])")
                        nchecked += 1
                        if errors and len(errors) > 10:
                            return errors
        for probe in CLAMP_PROBES:
            try:
                corr.evaluate(*probe)
            except Exception as e:  # noqa: BLE001
                errors.append(f"{name}{probe} must clamp, raised: {e}")
        print(f"  {name}: {nchecked} grid points OK, clamp probes OK")
    return errors


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("payload", nargs="?", help="candidate payload JSON(.gz)")
    ap.add_argument("--write-example", metavar="OUT", help="emit the reference example, then validate it")
    args = ap.parse_args()
    path = args.write_example or args.payload
    if not path:
        ap.error("give a payload to validate or --write-example OUT")
    if args.write_example:
        write_example(path)
    errors = validate(path)
    if errors:
        print("SPEC VIOLATIONS:")
        for e in errors:
            print("  -", e)
        sys.exit(1)
    print(f"RESULT: {path} is spec-compliant")


if __name__ == "__main__":
    main()
