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
# HashPRNG nodes (REQUIRED since the synthesis-throw factorization):
# name -> (distribution, exact entropy hash order). Orders are PAIRWISE DISTINCT --
# that is the decorrelation mechanism (alpha/beta throws independent; gate independent).
PRNG_CORRECTIONS = {
    "spx_angle_prng": ("stdnormal", ["layer", "cotAlpha", "cotBeta", "bLocalY"]),
    "spx_angle_prng_beta": ("stdnormal", ["cotBeta", "bLocalY", "layer", "cotAlpha"]),
    "spx_angle_valid_flat": ("stdflat", ["bLocalY", "cotBeta", "cotAlpha", "layer"]),
}
# Fused-shift terminal corrections (REQUIRED): bias-table structure, Formula cells,
# 5-input signature (the 4 plain inputs + prngAcc).
FINAL_CORRECTIONS = ["spx_angle_alpha_final", "spx_angle_beta_final"]
# Two-piece smear CompoundCorrections (REQUIRED): name -> exact stack; output_op "*".
COMPOUNDS = {
    "spx_angle_alpha_smear": ["spx_angle_alpha_sigma", "spx_angle_prng"],
    "spx_angle_beta_smear": ["spx_angle_beta_sigma", "spx_angle_prng_beta"],
}
# Fused shift CompoundCorrections (REQUIRED, the PRIMARY consumer contract):
# name -> exact stack; inputs_update=["prngAcc"], input_op "*", output_op "last";
# consumer passes prngAcc=1.0.
SHIFT_COMPOUNDS = {
    "spx_angle_alpha_shift": ["spx_angle_alpha_sigma", "spx_angle_prng", "spx_angle_alpha_final"],
    "spx_angle_beta_shift": ["spx_angle_beta_sigma", "spx_angle_prng_beta", "spx_angle_beta_final"],
}
MISSING_PRNG_MSG = (
    "payload predates the HashPRNG synthesis-throw factorization (missing {names}). "
    "Regenerate it with ngtagger-train/eval_spixel_angles/extract_pixelav_angle_payload.py "
    "-- the new nodes are purely additive (plain corrections stay bit-identical)."
)
INPUTS = [("layer", "int"), ("cotAlpha", "real"), ("cotBeta", "real"), ("bLocalY", "real")]
INPUTS_SHIFT = INPUTS + [("prngAcc", "real")]
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
    input_block = [
        {"name": "layer", "type": "int", "description": "TBPX layer 1-4 (entropy)"},
        {"name": "cotAlpha", "type": "real", "description": "true local px/pz (entropy)"},
        {"name": "cotBeta", "type": "real", "description": "true local py/pz (entropy)"},
        {"name": "bLocalY", "type": "real", "description": "local-y B [T] (entropy)"},
    ]
    input_block_shift = input_block + [
        {"name": "prngAcc", "type": "real",
         "description": "accumulator seed -- consumer MUST pass 1.0"},
    ]
    for name, (dist, hash_order) in PRNG_CORRECTIONS.items():
        corrections.append(
            {
                "name": name,
                "description": f"EXAMPLE deterministic {dist} variate (HashPRNG over the input tuple).",
                "version": 0,
                "inputs": input_block,
                "output": {"name": name, "type": "real"},
                "data": {"nodetype": "hashprng", "inputs": hash_order, "distribution": dist},
            }
        )
    for name in FINAL_CORRECTIONS:
        mb = _example_multibinning("bias")
        mb["content"] = [
            {"nodetype": "formula", "expression": f"{b} + x", "parser": "TFormula",
             "variables": ["prngAcc"]}
            for b in mb["content"]
        ]
        corrections.append(
            {
                "name": name,
                "description": "EXAMPLE fused-shift terminal: per-bin Formula '<bias> + prngAcc'.",
                "version": 0,
                "inputs": input_block_shift,
                "output": {"name": name, "type": "real"},
                "data": {
                    "nodetype": "category",
                    "input": "layer",
                    "content": [{"key": l, "value": mb} for l in LAYERS],
                },
            }
        )
    compound_corrections = [
        {
            "name": name,
            "description": ("EXAMPLE stochastic synthesis term: sigma * N(0,1). "
                            "cotX_meas = cotX_true + bias + smear."),
            "inputs": input_block,
            "output": {"name": name, "type": "real"},
            "inputs_update": [],
            "input_op": "*",
            "output_op": "*",
            "stack": stack,
        }
        for name, stack in COMPOUNDS.items()
    ] + [
        {
            "name": name,
            "description": ("EXAMPLE FUSED synthesis shift: bias + sigma*N(0,1) in one "
                            "evaluate; consumer passes prngAcc=1.0."),
            "inputs": input_block_shift,
            "output": {"name": name, "type": "real"},
            "inputs_update": ["prngAcc"],
            "input_op": "*",
            "output_op": "last",
            "stack": stack,
        }
        for name, stack in SHIFT_COMPOUNDS.items()
    ]
    cset = {
        "schema_version": 2,
        "description": (
            "EXAMPLE PixelAV SmartPixels angle-response payload; see "
            "doc/PixelAVAngleResponseSpec.md. NOT physics — structural reference only. "
            "CONSUMER CONTRACT: cotX_meas = cotX_true + spx_angle_X_bias(inputs) + "
            "spx_angle_X_smear(inputs); accept iff spx_angle_valid_flat < spx_angle_valid_prob."
        ),
        "corrections": corrections,
        "compound_corrections": compound_corrections,
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
    # HashPRNG factorization nodes are a HARD requirement (no back-compat: we control
    # all payloads; see spec section 3/8).
    missing_prng = (set(PRNG_CORRECTIONS) | set(FINAL_CORRECTIONS)) - present
    raw_compounds = {c.get("name"): c for c in raw.get("compound_corrections", [])}
    missing_comp = (set(COMPOUNDS) | set(SHIFT_COMPOUNDS)) - set(raw_compounds)
    if missing_prng or missing_comp:
        errors.append(MISSING_PRNG_MSG.format(names=sorted(missing_prng | missing_comp)))
    for cjson in raw.get("corrections", []):
        name = cjson.get("name")
        if name not in CORRECTIONS and name not in PRNG_CORRECTIONS and name not in FINAL_CORRECTIONS:
            continue
        got = [(i["name"], i["type"]) for i in cjson.get("inputs", [])]
        want = INPUTS_SHIFT if name in FINAL_CORRECTIONS else INPUTS
        if got != want:
            errors.append(f"{name}: inputs {got} != required {want}")
        if not cjson.get("description"):
            errors.append(f"{name}: description/provenance missing")
        if name in PRNG_CORRECTIONS:
            dist, order = PRNG_CORRECTIONS[name]
            data = cjson.get("data", {})
            if data.get("nodetype") != "hashprng":
                errors.append(f"{name}: nodetype {data.get('nodetype')} != hashprng")
            if data.get("distribution") != dist:
                errors.append(f"{name}: distribution {data.get('distribution')} != {dist}")
            if data.get("inputs") != order:
                errors.append(f"{name}: entropy order {data.get('inputs')} != required {order} "
                              "(pairwise-distinct orders are the decorrelation mechanism)")
    for name, stack in COMPOUNDS.items():
        cjson = raw_compounds.get(name)
        if cjson is None:
            continue
        got = [(i["name"], i["type"]) for i in cjson.get("inputs", [])]
        if got != INPUTS:
            errors.append(f"{name}: inputs {got} != required {INPUTS}")
        if cjson.get("stack") != stack:
            errors.append(f"{name}: stack {cjson.get('stack')} != required {stack}")
        if cjson.get("output_op") != "*":
            errors.append(f"{name}: output_op {cjson.get('output_op')} != '*'")
    for name, stack in SHIFT_COMPOUNDS.items():
        cjson = raw_compounds.get(name)
        if cjson is None:
            continue
        got = [(i["name"], i["type"]) for i in cjson.get("inputs", [])]
        if got != INPUTS_SHIFT:
            errors.append(f"{name}: inputs {got} != required {INPUTS_SHIFT}")
        if cjson.get("stack") != stack:
            errors.append(f"{name}: stack {cjson.get('stack')} != required {stack}")
        if cjson.get("inputs_update") != ["prngAcc"]:
            errors.append(f"{name}: inputs_update {cjson.get('inputs_update')} != ['prngAcc']")
        if cjson.get("input_op") != "*":
            errors.append(f"{name}: input_op {cjson.get('input_op')} != '*'")
        if cjson.get("output_op") != "last":
            errors.append(f"{name}: output_op {cjson.get('output_op')} != 'last'")
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

    # ---- HashPRNG factorization: functional checks -------------------------------
    probes = [(l, ca, cb, by) for l in LAYERS for ca in (-0.15, 0.02, 0.31)
              for cb in (-2.0, 0.6) for by in BLOCALY_GRID]
    for pr in probes:
        u = cset["spx_angle_valid_flat"].evaluate(*pr)
        if not (0.0 <= u <= 1.0):
            errors.append(f"spx_angle_valid_flat{pr}={u} not in [0,1]")
        for cname, stack in COMPOUNDS.items():
            z = cset[stack[1]].evaluate(*pr)
            if not math.isfinite(z):
                errors.append(f"{stack[1]}{pr} not finite")
            v = cset.compound[cname].evaluate(*pr)
            want = cset[stack[0]].evaluate(*pr) * z
            if abs(v - want) > 1e-12 * max(1.0, abs(want)):
                errors.append(f"{cname}{pr}={v} != sigma*prng={want}")
            if v != cset.compound[cname].evaluate(*pr):
                errors.append(f"{cname}{pr} non-deterministic across repeated evals")
        # fused == two-piece: shift(in, 1.0) == bias(in) + smear(in), to 1e-12
        for x in ("alpha", "beta"):
            fused = cset.compound[f"spx_angle_{x}_shift"].evaluate(*pr, 1.0)
            two = (cset[f"spx_angle_{x}_bias"].evaluate(*pr)
                   + cset.compound[f"spx_angle_{x}_smear"].evaluate(*pr))
            if abs(fused - two) > 1e-12 * max(1.0, abs(two)):
                errors.append(f"spx_angle_{x}_shift{pr} fused={fused} != bias+smear={two}")
    # alpha/beta throws must NOT be identical (independent entropy permutations)
    pr0 = probes[0]
    if cset["spx_angle_prng"].evaluate(*pr0) == cset["spx_angle_prng_beta"].evaluate(*pr0):
        errors.append("spx_angle_prng and spx_angle_prng_beta return identical deviates "
                      "(throws must be independent)")
    cset2 = correctionlib.CorrectionSet.from_file(path)
    for cname in ("spx_angle_alpha_smear", "spx_angle_beta_smear"):
        if cset.compound[cname].evaluate(*pr0) != cset2.compound[cname].evaluate(*pr0):
            errors.append(f"{cname} non-deterministic across fresh loads")
    for cname in ("spx_angle_alpha_shift", "spx_angle_beta_shift"):
        if cset.compound[cname].evaluate(*pr0, 1.0) != cset2.compound[cname].evaluate(*pr0, 1.0):
            errors.append(f"{cname} non-deterministic across fresh loads")
    if not errors:
        print(f"  hashprng/compounds: {len(probes)} probes OK (stack product exact, "
              "fused==bias+smear, deterministic, throws independent, gate in [0,1])")
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
