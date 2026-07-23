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
  - noise_cotAlpha : per-layer inverse CDF (layer, quantile) of the inclusive
  - noise_cotBeta    wrong-hit cotAlpha/cotBeta population, consumed by the
                     digiRefit producer to draw angles for no-link (noise)
                     digis: cot = F^-1(layer, u), u ~ U(0,1)

True noise digis carry no simlink and hence no measurable parent angle, so the
noise inverse CDFs use the inclusive in-window wrong-hit population (class-1,
plus any class-2 with a valid parent) as the stand-in: a noise cluster's
reconstructed angle should look like "some random in-detector cluster's angle".
The producer draws alpha and beta with independent uniforms (correlation
neglected, v1). Layers with too few entries fall back to the pooled all-layer
distribution; if the pool itself is empty the noise corrections are NOT written
and the producer falls back to position-only noise digis (its designed default).

Robust sigma = (q84 - q16)/2 (v0). Double-Gaussian core+tail is the documented
upgrade path. Angle residuals use the parent's ORIGIN-momentum local angle
(origin-momentum approximation, matching the analyzer).

Provenance (input file, tree, crossings/events, git rev, timestamp, optional
--sample/--provenance free text) is recorded in the CorrectionSet description
and every correction description, per the payload-provenance requirement in
doc/PixelAVAngleResponseSpec.md.

Run (cmsenv):
  python3 L1Trigger/Phase3SmartPixels/python/fitSmartHitPayloads.py \
      --in_file payload_ntuple.root \
      --out_true smarthit_true_v0.json \
      --out_fake smarthit_fake_v0.json \
      --sample RelValTTbar_D121_noPU --version 1
"""
import argparse
import datetime
import os
import subprocess
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


def _git_rev():
    """Best-effort git revision of the tree containing this script ('unknown' offline)."""
    here = os.path.dirname(os.path.abspath(__file__))
    try:
        out = subprocess.run(["git", "-C", here, "rev-parse", "HEAD"],
                             capture_output=True, text=True, timeout=10)
        rev = out.stdout.strip()[:12]
        if out.returncode == 0 and rev:
            dirty = subprocess.run(["git", "-C", here, "diff", "--quiet", "HEAD"],
                                   capture_output=True, timeout=10)
            if dirty.returncode == 1:
                rev += "+dirty"
            return rev
    except Exception:
        pass
    return "unknown"


def _inverse_cdf_binning(values, n_quantile_bins, name="quantile"):
    """Step-function inverse CDF as a Binning over the uniform quantile u in [0,1].

    content[i] = empirical quantile of `values` at the i-th bin midpoint, so
    evaluate(u ~ U(0,1)) reproduces the empirical distribution to a resolution
    of 1/n_quantile_bins. Content is non-decreasing by construction; clamp flow
    covers u == 0 and u == 1 exactly.
    """
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    edges = np.linspace(0.0, 1.0, n_quantile_bins + 1)
    mids = 0.5 * (edges[:-1] + edges[1:])
    content = np.quantile(v, mids)
    return _binning(edges, content, name)


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


def build_payloads(in_file, out_true, out_fake, tree="smartPixelsPayloadAnalyzer/crossings",
                   sample="", extra_provenance="", version=1, noise_quantile_bins=40):
    f = uproot.open(in_file)
    t = f[tree]
    arr = t.arrays(
        [
            "event", "trk_eta", "trk_tpMatched", "layer", "trk_cotAlpha", "trk_cotBeta",
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

    # ---- provenance: recorded in the CorrectionSet description and every
    # correction description (payload-provenance requirement, see
    # doc/PixelAVAngleResponseSpec.md section 6) ----
    n_events = int(np.unique(arr["event"]).size)
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    provenance = (
        f"fitSmartHitPayloads.py rev={_git_rev()} generated={stamp}; "
        f"input={os.path.basename(in_file)} tree={tree} crossings={n} events={n_events} "
        f"layers={layers} abs_eta_edges={ABS_ETA_EDGES}"
    )
    if sample:
        provenance += f"; sample={sample}"
    if extra_provenance:
        provenance += f"; {extra_provenance}"
    print(f"Provenance: {provenance}")

    def _desc(desc):
        return f"{desc} | {provenance}"

    # per-(layer, eta-bin) accumulators
    A = {l: {"eff": [0.0] * nbins, "rx": [0.0] * nbins, "ry": [0.0] * nbins,
             "aa": [0.0] * nbins, "ab": [0.0] * nbins} for l in layers}
    B = {l: {"mo": [0.0] * nbins, "mn": [0.0] * nbins,
             "aa": [0.0] * nbins, "ab": [0.0] * nbins} for l in layers}
    # per-layer inclusive wrong-hit angle pools (eta-inclusive) for the
    # smarthit_noise_* inverse CDFs
    noise_a = {l: [] for l in layers}
    noise_b = {l: [] for l in layers}

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
                    noise_a[l].extend(pA[ga].tolist())
                    noise_b[l].extend(pB[ga].tolist())

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
            name=key, version=version, description=_desc(desc),
            inputs=[cs.Variable(name="layer", type="int", description="TBPX layer (1..4)"),
                    cs.Variable(name="abs_eta", type="real", description="|eta| of the L1 track")],
            output=cs.Variable(name=key, type="real", description=desc),
            data=_layer_category(layers, [_binning(ABS_ETA_EDGES, A[l][per_layer_key], "abs_eta") for l in layers]),
        )

    cset_true = cs.CorrectionSet(
        schema_version=2,
        description=_desc("SmartPixels Tier-2 Stack A smarthit_true payload"),
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
            name=key, version=version, description=_desc(desc),
            inputs=[cs.Variable(name="layer", type="int", description="TBPX layer (1..4)"),
                    cs.Variable(name="abs_eta", type="real", description="|eta| of the L1 track")],
            output=cs.Variable(name=key, type="real", description=desc),
            data=_layer_category(layers, [_binning(ABS_ETA_EDGES, B[l][per_layer_key], "abs_eta") for l in layers]),
        )

    fake_corrs = [
        stackB_corr("smarthit_fake_mult_otherTP", "mean class-1 (other-TP) digis per window", "mo"),
        stackB_corr("smarthit_fake_mult_noise", "mean class-2 (noise) digis per window", "mn"),
        stackB_corr("smarthit_fake_ang_alpha_width", "inclusive fake cotAlpha robust width", "aa"),
        stackB_corr("smarthit_fake_ang_beta_width", "inclusive fake cotBeta robust width", "ab"),
    ]

    # ---- smarthit_noise_* inverse CDFs (consumed by the digiRefit producer
    # for no-link digis: cot = F^-1(layer, u), u ~ U(0,1), independent draws
    # for alpha and beta). Stand-in population: inclusive in-window wrong-hit
    # parent angles (see module docstring). ----
    kMinNoiseEntries = 10
    pooled_a = [x for l in layers for x in noise_a[l]]
    pooled_b = [x for l in layers for x in noise_b[l]]
    if len(pooled_a) >= kMinNoiseEntries and len(pooled_b) >= kMinNoiseEntries:
        def noise_corr(key, desc, pools, pooled):
            per_layer = [pools[l] if len(pools[l]) >= kMinNoiseEntries else pooled for l in layers]
            return cs.Correction(
                name=key, version=version, description=_desc(desc),
                inputs=[cs.Variable(name="layer", type="int", description="TBPX layer (1..4)"),
                        cs.Variable(name="quantile", type="real",
                                    description="uniform random quantile u in [0,1]")],
                output=cs.Variable(name=key, type="real", description=desc),
                data=cs.Category(
                    nodetype="category",
                    input="layer",
                    content=[cs.CategoryItem(key=int(l),
                                             value=_inverse_cdf_binning(p, noise_quantile_bins))
                             for l, p in zip(layers, per_layer)],
                    default=_inverse_cdf_binning(pooled, noise_quantile_bins),
                ),
            )
        fake_corrs.append(noise_corr(
            "smarthit_noise_cotAlpha",
            "inverse CDF of the inclusive wrong-hit cotAlpha population (noise-angle stand-in)",
            noise_a, pooled_a))
        fake_corrs.append(noise_corr(
            "smarthit_noise_cotBeta",
            "inverse CDF of the inclusive wrong-hit cotBeta population (noise-angle stand-in)",
            noise_b, pooled_b))
        print(f"Noise inverse CDFs written from {len(pooled_a)} inclusive wrong-hit angles "
              f"(per-layer: {[len(noise_a[l]) for l in layers]})")
    else:
        print("WARNING: too few inclusive wrong-hit angles "
              f"({len(pooled_a)} < {kMinNoiseEntries}); smarthit_noise_* corrections NOT written "
              "(producer falls back to position-only noise digis)")

    cset_fake = cs.CorrectionSet(
        schema_version=2,
        description=_desc("SmartPixels Tier-2 Stack B smarthit_fake payload"),
        corrections=fake_corrs,
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
    if any(c.name == "smarthit_noise_cotAlpha" for c in cset_fake.corrections):
        for u in (0.05, 0.5, 0.95):
            na = fset["smarthit_noise_cotAlpha"].evaluate(test_layer, u)
            nb = fset["smarthit_noise_cotBeta"].evaluate(test_layer, u)
            print(f"  noise draw layer={test_layer} u={u}: cotAlpha={na:.3f} cotBeta={nb:.3f}")
        # an inverse CDF must be non-decreasing in u, for every layer and the default
        probe_layers = layers + [max(layers) + 1]  # +1 exercises the Category default
        for nm in ("smarthit_noise_cotAlpha", "smarthit_noise_cotBeta"):
            for l in probe_layers:
                vals = [fset[nm].evaluate(int(l), float(u)) for u in np.linspace(0.0, 1.0, 101)]
                if not all(y2 >= y1 - 1e-12 for y1, y2 in zip(vals, vals[1:])):
                    raise RuntimeError(f"{nm} layer {l}: inverse CDF not monotone")
        print("  noise inverse-CDF monotonicity (all layers + default): OK")

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
    ap.add_argument("--sample", default="", help="sample name recorded in the payload provenance")
    ap.add_argument("--provenance", default="", help="extra free-text provenance to record")
    ap.add_argument("--version", type=int, default=1, help="version stamped on every correction")
    ap.add_argument("--noise_quantile_bins", type=int, default=40,
                    help="quantile bins for the smarthit_noise_* inverse CDFs")
    opts = ap.parse_args()
    build_payloads(opts.in_file, opts.out_true, opts.out_fake, tree=opts.tree,
                   sample=opts.sample, extra_provenance=opts.provenance,
                   version=opts.version, noise_quantile_bins=opts.noise_quantile_bins)
