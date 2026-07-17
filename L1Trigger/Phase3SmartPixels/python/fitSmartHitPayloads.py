"""Offline fit: SmartPixelsPayloadAnalyzer ntuple -> Stack A/B correctionlib JSONs.

Reads the per-crossing ntuple written by SmartPixelsPayloadAnalyzer and derives,
per (TBPX layer, |eta| bin), the v0 "smarthit_true" (Stack A) and
"smarthit_fake" (Stack B) payloads used by the Tier-2 digiRefit producer
(Phase 2).

Stack A "smarthit_true" (per layer, |eta|):
  - eff            : true-hit efficiency = P(>=1 class-0 digi in window | TP-matched crossing)
  - res_x_sigma    : robust position residual width in local x [cm] (class-0 digis)
  - res_y_sigma    : robust position residual width in local y [cm] (class-0 digis)
  - ang_alpha_sigma: robust width of (parent cotAlpha - track cotAlpha) [class-0]
  - ang_beta_sigma : robust width of (parent cotBeta  - track cotBeta ) [class-0]

Stack B "smarthit_fake" (per layer, |eta|):
  - mult_otherTP   : mean number of class-1 (other-TP) digis per window
  - mult_noise     : mean number of class-2 (noise) digis per window
  - ang_alpha_width: robust width of the inclusive class-1/2 cotAlpha distribution
  - ang_beta_width : robust width of the inclusive class-1/2 cotBeta distribution

Robust sigma = (q84 - q16)/2 (v0). Double-Gaussian core+tail is the documented
upgrade path. Angle residuals use the parent's ORIGIN-momentum local angle
(origin-momentum approximation, matching the analyzer).

Run (cmsenv):
  python3 L1Trigger/Phase3SmartPixels/python/fitSmartHitPayloads.py \
      --in_file payload_ntuple.root \
      --out_true smarthit_true_v0.json \
      --out_fake smarthit_fake_v0.json
"""
import argparse
import numpy as np
import uproot
import correctionlib
import correctionlib.schemav2 as cs

# |eta| binning shared by both stacks (barrel coverage; clamp at the ends).
ABS_ETA_EDGES = [0.0, 0.3, 0.6, 0.9, 1.2, 1.6, 2.4]


def robust_sigma(x):
    """(q84-q16)/2; NaN-safe, returns 0.0 for empty/degenerate input."""
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size < 2:
        return 0.0
    q16, q84 = np.quantile(x, [0.16, 0.84])
    return float((q84 - q16) / 2.0)


def _binning(edges, content, name, flow="clamp"):
    return cs.Binning(
        nodetype="binning",
        input=name,
        edges=list(map(float, edges)),
        content=[float(c) for c in content],
        flow=flow,
    )


def _layer_category(layers, per_layer_binnings, output_default=0.0):
    """cs.Category over integer 'layer' -> a per-layer abs_eta cs.Binning."""
    return cs.Category(
        nodetype="category",
        input="layer",
        content=[cs.CategoryItem(key=int(l), value=b) for l, b in zip(layers, per_layer_binnings)],
        default=cs.Binning(
            nodetype="binning",
            input="abs_eta",
            edges=[ABS_ETA_EDGES[0], ABS_ETA_EDGES[-1]],
            content=[output_default],
            flow="clamp",
        ),
    )


def build_payloads(in_file, out_true, out_fake, tree="smartPixelsPayloadAnalyzer/crossings"):
    f = uproot.open(in_file)
    t = f[tree]
    arr = t.arrays(
        [
            "trk_eta", "trk_tpMatched", "layer", "trk_cotAlpha", "trk_cotBeta",
            "nWinDigi", "digi_class", "digi_resx", "digi_resy",
            "digi_parCotAlpha", "digi_parCotBeta",
        ],
        library="np",
    )
    n = len(arr["layer"])
    print(f"Loaded {n} crossings from {in_file}")
    if n == 0:
        raise RuntimeError("empty ntuple: no crossings to fit")

    abs_eta = np.abs(arr["trk_eta"])
    layers = sorted(set(int(x) for x in np.unique(arr["layer"])))
    eta_lo = np.array(ABS_ETA_EDGES[:-1])
    eta_hi = np.array(ABS_ETA_EDGES[1:])
    nbins = len(eta_lo)

    # per-(layer, eta-bin) accumulators
    A = {l: {"eff": [0.0] * nbins, "rx": [0.0] * nbins, "ry": [0.0] * nbins,
             "aa": [0.0] * nbins, "ab": [0.0] * nbins} for l in layers}
    B = {l: {"mo": [0.0] * nbins, "mn": [0.0] * nbins,
             "aa": [0.0] * nbins, "ab": [0.0] * nbins} for l in layers}

    summary_rows = []
    for l in layers:
        for b in range(nbins):
            sel = (arr["layer"] == l) & (abs_eta >= eta_lo[b]) & (abs_eta < eta_hi[b])
            idx = np.nonzero(sel)[0]
            n_cross = idx.size
            n_matched = 0
            n_matched_hit = 0
            resx, resy = [], []
            ang_a_true, ang_b_true = [], []
            ang_a_fake, ang_b_fake = [], []
            mult_o, mult_n = [], []
            for i in idx:
                cls = np.asarray(arr["digi_class"][i])
                trkA = arr["trk_cotAlpha"][i]
                trkB = arr["trk_cotBeta"][i]
                matched = arr["trk_tpMatched"][i] == 1
                # Stack A: efficiency + true-hit residuals (class-0)
                m0 = cls == 0
                if matched:
                    n_matched += 1
                    if np.any(m0):
                        n_matched_hit += 1
                        resx.extend(np.asarray(arr["digi_resx"][i])[m0].tolist())
                        resy.extend(np.asarray(arr["digi_resy"][i])[m0].tolist())
                        pA = np.asarray(arr["digi_parCotAlpha"][i])[m0]
                        pB = np.asarray(arr["digi_parCotBeta"][i])[m0]
                        good = pA > -900
                        ang_a_true.extend((pA[good] - trkA).tolist())
                        ang_b_true.extend((pB[good] - trkB).tolist())
                # Stack B: window multiplicity + inclusive fake-angle spread
                m1 = cls == 1
                m2 = cls == 2
                mult_o.append(int(np.count_nonzero(m1)))
                mult_n.append(int(np.count_nonzero(m2)))
                for mm in (m1, m2):
                    pA = np.asarray(arr["digi_parCotAlpha"][i])[mm]
                    pB = np.asarray(arr["digi_parCotBeta"][i])[mm]
                    ga = pA > -900
                    ang_a_fake.extend(pA[ga].tolist())
                    ang_b_fake.extend(pB[ga].tolist())

            eff = (n_matched_hit / n_matched) if n_matched else 0.0
            A[l]["eff"][b] = eff
            A[l]["rx"][b] = robust_sigma(resx)
            A[l]["ry"][b] = robust_sigma(resy)
            A[l]["aa"][b] = robust_sigma(ang_a_true)
            A[l]["ab"][b] = robust_sigma(ang_b_true)
            B[l]["mo"][b] = float(np.mean(mult_o)) if mult_o else 0.0
            B[l]["mn"][b] = float(np.mean(mult_n)) if mult_n else 0.0
            B[l]["aa"][b] = robust_sigma(ang_a_fake)
            B[l]["ab"][b] = robust_sigma(ang_b_fake)
            summary_rows.append(
                (l, f"{eta_lo[b]:.1f}-{eta_hi[b]:.1f}", n_cross, n_matched,
                 round(eff, 3), round(A[l]["rx"][b], 5), round(B[l]["mo"][b], 2),
                 round(B[l]["mn"][b], 2)))

    # ---- Stack A correctionlib set ----
    def stackA_corr(key, desc, per_layer_key):
        return cs.Correction(
            name=key, version=1, description=desc,
            inputs=[cs.Variable(name="layer", type="int", description="TBPX layer (1..4)"),
                    cs.Variable(name="abs_eta", type="real", description="|eta| of the L1 track")],
            output=cs.Variable(name=key, type="real", description=desc),
            data=_layer_category(layers, [_binning(ABS_ETA_EDGES, A[l][per_layer_key], "abs_eta") for l in layers]),
        )

    cset_true = cs.CorrectionSet(
        schema_version=2,
        corrections=[
            stackA_corr("smarthit_true_eff", "true-hit efficiency P(class-0 in window | TP-matched)", "eff"),
            stackA_corr("smarthit_true_res_x_sigma", "true-hit local-x residual robust sigma [cm]", "rx"),
            stackA_corr("smarthit_true_res_y_sigma", "true-hit local-y residual robust sigma [cm]", "ry"),
            stackA_corr("smarthit_true_ang_alpha_sigma", "true-hit (parent-track) cotAlpha robust sigma", "aa"),
            stackA_corr("smarthit_true_ang_beta_sigma", "true-hit (parent-track) cotBeta robust sigma", "ab"),
        ],
    )

    def stackB_corr(key, desc, per_layer_key):
        return cs.Correction(
            name=key, version=1, description=desc,
            inputs=[cs.Variable(name="layer", type="int", description="TBPX layer (1..4)"),
                    cs.Variable(name="abs_eta", type="real", description="|eta| of the L1 track")],
            output=cs.Variable(name=key, type="real", description=desc),
            data=_layer_category(layers, [_binning(ABS_ETA_EDGES, B[l][per_layer_key], "abs_eta") for l in layers]),
        )

    cset_fake = cs.CorrectionSet(
        schema_version=2,
        corrections=[
            stackB_corr("smarthit_fake_mult_otherTP", "mean class-1 (other-TP) digis per window", "mo"),
            stackB_corr("smarthit_fake_mult_noise", "mean class-2 (noise) digis per window", "mn"),
            stackB_corr("smarthit_fake_ang_alpha_width", "inclusive fake cotAlpha robust width", "aa"),
            stackB_corr("smarthit_fake_ang_beta_width", "inclusive fake cotBeta robust width", "ab"),
        ],
    )

    with open(out_true, "w") as fo:
        fo.write(cset_true.model_dump_json(exclude_unset=True))
    with open(out_fake, "w") as fo:
        fo.write(cset_fake.model_dump_json(exclude_unset=True))
    print(f"Wrote {out_true} and {out_fake}")

    # ---- validation: reload via correctionlib and evaluate a few points ----
    tset = correctionlib.CorrectionSet.from_file(out_true)
    fset = correctionlib.CorrectionSet.from_file(out_fake)
    test_layer = layers[0]
    for probe_eta in (0.1, 0.5, 1.0):
        eff = tset["smarthit_true_eff"].evaluate(test_layer, probe_eta)
        rx = tset["smarthit_true_res_x_sigma"].evaluate(test_layer, probe_eta)
        mo = fset["smarthit_fake_mult_otherTP"].evaluate(test_layer, probe_eta)
        print(f"  eval layer={test_layer} |eta|={probe_eta}: eff={eff:.3f} res_x_sigma={rx:.5f} mult_otherTP={mo:.2f}")

    # ---- summary table ----
    print("\nSummary (layer, |eta|, nCross, nMatched, eff, res_x_sigma[cm], mult_otherTP, mult_noise):")
    print(f"  {'lay':>3} {'|eta|':>9} {'nCr':>6} {'nMat':>6} {'eff':>6} {'resx':>8} {'multO':>6} {'multN':>6}")
    for r in summary_rows:
        print(f"  {r[0]:>3} {r[1]:>9} {r[2]:>6} {r[3]:>6} {r[4]:>6} {r[5]:>8} {r[6]:>6} {r[7]:>6}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--in_file", default="payload_ntuple.root")
    ap.add_argument("--out_true", default="smarthit_true_v0.json")
    ap.add_argument("--out_fake", default="smarthit_fake_v0.json")
    ap.add_argument("--tree", default="smartPixelsPayloadAnalyzer/crossings")
    opts = ap.parse_args()
    build_payloads(opts.in_file, opts.out_true, opts.out_fake, tree=opts.tree)
