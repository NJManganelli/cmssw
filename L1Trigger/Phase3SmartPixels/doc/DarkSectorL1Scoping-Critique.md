# Critique + Contributions — DarkSectorL1Scoping.md (v0)

Reviewer role: friendly-but-rigorous critic/collaborator. Research/writing only.
Companion to `DarkSectorL1Scoping.md`; the author revises against this — do not merge inline.
Date: 2026-07-20. Web research was available; all new citations verified against source text
(HTML/ar5iv mirrors where PDFs would not text-extract). Verdict headers flag confidence.

**Headline (read this first).** The document is well-structured, honestly hedged, and the
infrastructure grounding is accurate. But the single most load-bearing claim — Rank 1, the SUEP
n_tracks-at-PV handle — is **contradicted by the two papers the draft itself cites**, in the exact
regime the draft leans on. The whitepaper (2203.07314) and its companion (2211.05720) state
plainly that at the CMS Phase-2 baseline pT>2 GeV, *a large-track-multiplicity trigger has
negligible efficiency for essentially all SUEP models*, not merely "marginal for low-mass." Rank 1
as written (multiplicity-led, sphericity as backup) is **inverted**: sphericity/isotropy and track
pT-sum (HT_trk) are the survivors; raw multiplicity is the variable that dies. This is fixable and
the fix strengthens the doc, but it is not cosmetic. Details in (A.1) and (F.1).

---

## (A) Load-bearing-claim audit

### A.1 — SUEP flagship (make-or-break). VERDICT: claim is directionally right that the physics
is track/vertex-shaped, but the SPECIFIC Rank-1 handle (n_tracks at the hardest vertex) is
**refuted at the pT>2 GeV baseline** by the draft's own sources. Re-anchor on sphericity + track-HT.

What the sources actually say (now extracted — the draft flagged these as un-extracted):

- **arXiv:2211.05720 (Optimizing Trigger-Level Track Reconstruction), verified text:**
  "For a track-trigger similar to the CMS HL-LHC upgrade, where pT>2 GeV is the baseline threshold,
  a trigger requiring a large track multiplicity is inefficient for most SUEP signal models."
  Quantified: at **pT>2 GeV there is negligible efficiency for all signals with any nTrack
  threshold considered**; at **pT>0.5 GeV** you reach 50% efficiency only for **mediators ≥200 GeV
  with nTrack>100**; at **pT>1 GeV**, 50% only for **mediators ≥600 GeV**. The SUEP charged-pT
  spectrum "peaks sharply in the lowest pT bins," most tracks **below 1 GeV**, and the spectrum
  depends on the dark-meson mass/temperature, **not** the mediator mass. Their recommended
  discriminants are **track pT-sum (HT_trk)** and **event-shape/isotropy**, explicitly *not* plain
  multiplicity.
- **arXiv:2203.07314 (whitepaper), verified text:** same conclusion; for SUEP, **d0 is not
  discriminating (all prompt)** — efficiency is entirely pT-threshold-driven; and even at low pT
  "the overall efficiency is still likely to be low, because high-energy ISR is required" for their
  benchmark. Recommend prompt-track reconstruction down to **pT≈1 GeV** (lower than the displaced-
  track threshold).
- **arXiv:2607.09621 (CMS scouting), verified text:** this search does **NOT** use a track-
  multiplicity trigger at all — it uses a **hadronic HT trigger (HT lowered 1 TeV→410 GeV via
  scouting)**; per-particle reco efficiency is **45–95% depending on pT**; offline SR selection is
  **n_constituent_SUEP > 50 AND boosted-sphericity S > 0.5**. The efficiency "can fall below 5% for
  low m_S or high T_D/m_φ" and "exceeds 90% for large m_S." NB the improvement it claims (~1 order
  of magnitude, especially for T_D and m_φ ≲ 2 GeV) is a *HT-threshold* gain, not a track-
  multiplicity gain. **The draft mis-attributes to 2607.09621 a "pT>2 GeV multiplicity is
  inefficient" statement that lives in 2203.07314 / 2211.05720, not in the scouting paper.** Fix the
  citation (A.4).

Per-regime verdict for the author (this is what the brief asked for):

| SUEP regime | L1 tractability with a track handle |
|---|---|
| High m_S (≳600 GeV–1 TeV), moderate T_D/m_φ | **Tractable** — hardest tracks poke above pT>1–2 GeV; sphericity + track-HT separable; this is the regime the doc can honestly claim. |
| Mid m_S (200–600 GeV) | **Marginal** — needs pT floor pushed to ≤1 GeV (below GTT baseline) to keep 50% eff; multiplicity alone insufficient. |
| Low m_S and/or high T_D, m_φ ≲ 2 GeV (soft, "more unclustered") | **Not L1-tractable via GTT tracks at pT>2 GeV.** This is the SUEP the scouting search bought back with a HT threshold, and it is precisely the regime where the GTT pT floor guts the signal. Honest ceiling: L1 does not recover this without either a sub-2-GeV track path or an HT/anomaly path. |

Net: the draft's own "honest prior" (high-mass separable, low-mass marginal) is **too optimistic on
the mechanism**. The correct statement is: high-mass SUEP is separable *if and only if the
discriminant is sphericity/track-HT, not multiplicity*, and low-mass SUEP is out of reach at the
GTT baseline threshold regardless of discriminant. The physics survives; the specific "n_tracks at
PV" flagship needs to be demoted under "sphericity + track-HT at PV."

Second, an under-examined premise: **"all SUEP tracks share one vertex, PU spreads over 200."**
True, but the discriminating power of *that* fact is weakest exactly where you need it — a soft SUEP
contributes few tracks *above threshold* at its PV, and at 200 PU the hardest PU vertex already
carries a real hard-scatter's tracks. The per-vertex n_tracks *excess* is a tail competition between
(few surviving hard SUEP tracks) and (a genuine hard-scatter PV). Sphericity is what breaks the
degeneracy — a hard-scatter PV is jetty (S→0), a SUEP PV is isotropic (S→1) — so **sphericity is
load-bearing and multiplicity is the assist**, the reverse of the draft's ordering.

### A.2 — Quirk GNN efficiency numbers. VERDICT: citation is ACCURATE; add three caveats the draft omits.

Verified against arXiv:2410.00269 (ar5iv): SM tracks **97.9%**, quirks with the SM-trained
(helical) pipeline **10.2%**, in the **8-layer** geometry. GNN (Exa.TrkX) quirk efficiency:
**91.5%** (well-behaved quirks, *no background*), **56.3%** (all quirks, no background), **71.9%**
(all quirks *with* SM background, benchmark), **79.0%** (25-layer geometry). So the draft's
"56–91%" band is real but the endpoints are **not comparable operating points** — 91.5% is the
easy, background-free, well-behaved subset; 56.3% is all-quirks-no-background; the realistic
with-background number is ~72%. The draft should cite **~72% (with background)** as the honest
headline, not the 91% ceiling.

Caveats the draft must add:
- **Offline only, no trigger study.** The paper does no L1/trigger implementation. Draft's §4 use of
  it as motivation for an *L1* subsystem is an extrapolation — flag it.
- **Pileup is 5–10 events, not 200.** "inclusion of 5–10 pileup events degrades performance by
  ~5–10%." Extrapolating to 200 PU is unwarranted and probably optimistic; the 10.2%/72% gap could
  shrink at 200 PU. Do not let the numbers read as PU-200-validated.
- **Quirk mass 100–5000 GeV, DY production, σ≈1 fb at 500 GeV; oscillation length ~0.1 cm to >600 cm.**
  The efficiency depends strongly on oscillation length; short-oscillation quirks look nearly helical
  (easy), long ones are the hard/novel case. "56–91%" hides this axis.

### A.3 — Disappearing-track handle: is layerHitMask expressible at L1? VERDICT: the draft's own
§2.5 already concedes the core problem, but §5/table placement UNDER-states it. CHALLENGE UPHELD and
sharpened.

The critic's framing is correct and the draft half-acknowledges it: the OT-stub-seeded L1 track flow
**never forms a track object with no OT continuation**, so a disappearing track is invisible to the
*finder*, not merely unobserved by a downstream tagger. `layerHitMask` is a property of a *formed
refit track* — if no track is formed, there is no row to carry a mask. So this is a **track-FINDER
problem**, same category as quirks, and it belongs adjacent to §2.6/§4, not as a "bucket b/c
observable." The draft says "(b)/(c) — the sidecar has the hit-pattern info, but ... producing it
needs the §4 subsystem" — that sentence is right, but the summary table (§2.7) lists the handle as
"inner hits, no OT / SmartPixels hit patterns, refit sidecar layerHitMask" which reads as if the
sidecar *provides* the observable. It does not, unless the object is formed upstream. Recommend:
explicitly reclassify disappearing tracks as a **track-formation problem co-located with quirks**,
both gated on §4, and stop implying `layerHitMask` is the enabling piece (it is downstream of the
missing piece). One nuance in the draft's favor: SmartPixels *inner* hits are exactly what a
short-stub inner-seeded finder would need — so §4 is genuinely the right home; just don't call it an
observable.

### A.4 — Numerics/citations that need checking or fixing.

- **MIS-ATTRIBUTION (fix):** the "pT>2 GeV plain-multiplicity inefficient for most SUEP" claim in
  §2.1 is cited to 2607.09621; it belongs to **2211.05720 / 2203.07314**. 2607.09621 is a *HT*-based
  scouting search. (A.1)
- **"10.2% → 56–91%"** — accurate but present the honest operating point (~72% w/ background); note
  offline / low-PU. (A.2)
- **CICADA/AXOL1TL latencies** ("sub-200 ns", "<100 ns"): plausible and consistent with public
  numbers; not independently re-verified here — low risk, leave but do not sharpen.
- **L1 track path "~4 µs + ~1 µs transfer"** (§3): the total **Phase-2 L1 latency is 12.5 µs**
  (verified), of which the track trigger is an allocation, not the whole budget. The "~4 µs" figure
  is in the right ballpark for the L1 tracking allocation but the draft should cite the **12.5 µs
  total** as the frame and note GTT sits inside it, else a reader thinks 4 µs is the whole L1.
- **|d0| ~ 10 cm displaced reach** (§2.2): consistent with the LST/HLT displaced literature
  (r_vertex ≳ 5 cm acceptance is the quoted new HLT feature); the "~10 cm" L1 extended-track number
  is plausible but I could not pin the exact GTT extended-track d0 ceiling to a primary firmware
  source — mark **needs-checking against the GTT/DisplacedVertexFinder firmware note**, not a review
  talk.
- **"Belt of fire" / isotropy as HLT discriminant** — supported (2211.05720 explicitly). Good to add.

---

## (B) Feasibility-bucket audit

Re-derived buckets, and the hidden dependencies the draft glosses:

- **Rank 1 (SUEP): draft says "(a) offline." PARTIALLY WRONG — it is (a) only after 3 non-trivial
  build items.** The "vertex tooling exists" claim is true for *per-z-bin n_tracks* (e87538e) but the
  offline reach study additionally needs, none of which the nano carries cleanly today:
  1. **A sphericity/isotropy implementation over PV-associated tracks** — new code (small, but the
     draft's own A.1 correction makes this the *primary* discriminant, so it is not optional).
  2. **A per-vertex track-to-vertex association** — the FastHisto histogram gives per-z-bin *counts*,
     but computing a sphericity tensor needs the **(px,py,pz) of the individual tracks assigned to the
     hardest bin**, i.e. an explicit track→vertex membership list, not just a bin count. Confirm the
     nano exposes track z0 + the vertex z so you can reconstruct membership; if only the histogram
     count is persisted, you cannot form the tensor. **This is a real hidden dependency — flag it.**
  3. **Truth-vertex matching** to define "the SUEP PV" for signal and to build the "hardest PU vertex"
     background envelope — ties directly to the truth-posture discipline (the draft's §6 cross-cutting
     point). At 200 PU, "hardest vertex" is an algorithmic choice (FastHisto argmax) that may not be
     the SUEP vertex; the study must decide truth-PV vs reco-hardest-PV and they diverge for soft SUEP.
  So Rank 1 is **(a) offline once you write sphericity + confirm track→vertex membership in nano**;
  honestly it is "(a) with a small (b)-flavored producer/nano dependency." Say so.
- **Rank 2 (anomaly score): (a/b) stands** — but it inherits Rank 1's track→vertex-membership
  dependency for its per-vertex features. Correct that they share the dependency, not just the
  extractor.
- **Rank 3 (emerging/displaced): (a/b) stands**, contingent on the `setTrackWordBits` nPar=5 refit-d0
  re-encode (the draft's own open item #3 / spec §7-D). Good that it is flagged; elevate it to a
  *precondition*, not a footnote — if d0 is silently zeroed the whole study is null.
- **Rank 4 (disappearing): (b/c) — agree, but per A.3 it is (c)-gated on §4, closer to Rank 5 than
  the table implies.**
- **Rank 5 (quirks): (c) stands.**

---

## (C) Filled literature gaps

### C.1 — SUEP soft-pT vs GTT pT>2 GeV (the Rank-1 make-or-break). CLOSED. See A.1.
Bottom line for the doc: **at pT>2 GeV, multiplicity is dead for SUEP; sphericity/isotropy + track-HT
are the survivors; sub-2-GeV or low-m_S SUEP is out of L1 reach at the GTT baseline.** Numbers in A.1.

### C.2 — Track-Based Triggers whitepaper (2203.07314) + companion (2211.05720). CLOSED (HTML mirror).
Parameterization recovered: efficiency is a function of **(pT threshold, d0 threshold, nTrack
threshold)**; for SUEP d0 is irrelevant (prompt), pT threshold dominates, and the recommended
future-trigger posture is **prompt tracking down to pT≈1 GeV, lower than the displaced-track pT
threshold**, plus **track-HT and event-shape** discriminants rather than bare multiplicity. For the
LLP/displaced benchmarks (long-lived staus, Higgs-portal scalars) the recommendation is **increasing
d0 coverage with pT** — i.e. displaced tracks can tolerate a higher pT floor. This directly supports
Rank 3's extended-d0 handle and simultaneously undercuts Rank 1's multiplicity framing.

### C.3 — Phase-2 GTT vertex + DisplacedVertexFinder firmware constraints. PARTIALLY CLOSED.
Confirmed: **Phase-2 L1 total latency 12.5 µs**; GTT runs on **Xilinx VU13P** (VU13P: ~1.3×10^4
DSP48E2, large URAM/BRAM); the flow is FastHisto/vertex → track-quality/vertex-based selections;
28 Gb/s links; hls4ml NN blocks demonstrated within ~50 ns / 2 clock ticks in related L1 ML work.
**Not closed:** I could not pull a primary GTT firmware note giving the *per-vertex accumulator*
budget (whether a sphericity tensor — 6 sums of products per track, accumulated per z-bin — fits the
existing FastHisto pass, and at what DSP/latency cost). This is the specific number the bucket-(b)
producer scope needs. Recommend the author route this as a **direct question to the GTT firmware
group** (the draft's §6 already has the right instinct; make it a concrete ask: "can the FastHisto
pass accumulate 6 additional per-z-bin sum-of-products for a sphericity tensor within budget?").
A sphericity tensor is cheap arithmetically (multiply-accumulate) so the honest prior is "plausible,"
but assert nothing without the firmware group.

---

## (D) Strategic remit — ADJUDICATED (coordinator/user decision folded in) + convergence pressure-test

**Remit decision (settled by the user via coordinator, not re-litigated here): BOTH tracks, in
explicit parallel, with an upstream-payoff frame.** Write §D of the doc as two named tracks:

- **Track 1 — near-term observables on the EXISTING track/vertex collections.** Sphericity/track-HT
  at the PV (Rank 1, corrected), extended-d0 + dispvtx (Rank 3), the unsupervised track/vertex
  anomaly score (Rank 2). These are sample-gated, prototypable in ngtagger now, and — crucially —
  **directly upstreamable to the baseline Phase-2 GTT/correlator**: they are observables on
  collections that exist in the standard (non-SmartPixels) Phase-2 track trigger. Their value does
  **not** depend on SmartPixels being adopted. Frame them as **"payoffs of the R&D line that survive
  SmartPixels non-adoption."**
- **Track 2 — motivate the NEW subsystem (graph-track condensation / reseeding).** Longer-horizon,
  design-gated, track-finder-group-owned. Frame explicitly as a **standalone contribution to the
  baseline Phase-2 track trigger** whose value is **decoupled from SmartPixels' own adoption risk** —
  a de-risking argument for the whole speculative R&D line: even if SmartPixels never enters CMS
  Phase-3 (or ATLAS), a non-helical/condensation reseed stage is a contribution to Phase-2 tracking
  in its own right. This is the strongest strategic sentence available and the doc should lead §4/§D
  with it.

**Convergence claim pressure-test (this is where I earn my keep — the "both + upstream" decision does
NOT soften this).** The claim: one hit-graph-condensation layer serves THREE signatures at once
(quirk recovery + SUEP/displaced reseeding + SmartPixels-OT merge). **Verdict: partially real, and
weaker than the draft implies — the three want DIFFERENT things from the subsystem, and only two of
the three are genuinely the same mechanism.**

- **Quirk recovery** needs a **non-helical track-formation model** — a hit-graph segmentation with
  *no trajectory prior* (the whole point of 2410.00269 / the 2509.08878 follow-up). This is the hard,
  novel capability. VERIFIED as a distinct, published need.
- **SUEP/displaced reseeding** needs **regional re-running of an (still essentially helical) finder
  at a lower pT floor / looser prompt assumption** in hit-dense pockets. A SUEP track *is* helical —
  it is just soft and dropped by the pT>2 GeV / nstub cut, not non-helical. So SUEP reseeding does
  **not** need the non-helical machinery at all; it needs a *lower-threshold regional refinder*. That
  is a real capability but a **different** one from quirk condensation. **The draft conflates "the
  layer has the hits so it can reseed" with "the same non-helical model serves both" — it does not.**
  Displaced reseeding is intermediate (helical but high-d0), closer to SUEP's need than quirk's.
- **SmartPixels-OT merge** needs a **spatial association / reconciliation** of inner smarttracks with
  OT tracks — an instance-matching / linking problem. Object-condensation latent clustering *could*
  express this, so there is genuine architectural overlap with quirk condensation (both are
  "instance-segment hits into objects without a trajectory fit"). This is the **real** convergence
  axis: **quirk recovery and smarttrack/OT merge both want object-condensation-style hit clustering.**
  SUEP/displaced reseeding is the odd one out (it wants a cheaper regional refinder, not condensation).

So the honest convergence statement is: **the subsystem serves three signatures, but via TWO
distinct capabilities — (i) object-condensation hit clustering (quirks + smarttrack/OT merge +
disappearing-track formation), and (ii) regional low-threshold reseeding (SUEP-soft + displaced).**
Both live in the same *placement* (the layer that holds all hits + partial tracks), which is the true
unifying claim — **it is a placement/architecture convergence, not a single-model convergence.** That
is still a strong argument (one subsystem, one home, two capabilities, four+ signatures) and it is
*more* defensible than "one condensation layer does everything." Rewrite the §4 convergence bullet to
this two-capability framing; it will survive a track-finder-group reading, whereas "one layer serves
three" invites the objection "reseeding isn't condensation."

---

## (E) Additive contributions (signatures/handles the draft under-developed or missed)

**E.1 — The quirk↔SUEP convergence is a GIFT the draft leaves on the table (highest-value add).**
arXiv:2506.11192 ("Soft-unclustered-energy patterns from quirks", PRD) shows quirk-antiquirk
de-excitation radiates **hundreds of soft isotropic pions** — i.e. **quirks produce a SUEP-like
signature** (high multiplicity, transverse-sphericity isotropy, tracks spatially closer than QCD).
This unifies §2.1 (SUEP) and §2.6 (quirks) into one physics story and makes the §4 subsystem serve a
*single coherent* target rather than two disconnected ones. **BUT** — and this is the pressure-test —
that same paper's *recommended* trigger strategy is to use the **hard resonant dijet from the final
quirk-antiquirk annihilation** to trigger, with the soft tracks as *offline* discovery enhancement.
That is exactly the **associated-production crutch the whole draft is trying to escape.** So the
quirk-SUEP link is a double-edged addition: it strengthens the "one subsystem, one signature family"
argument, but it is *also* a literature data-point that the community's own answer for these events is
"trigger on the hard object." The doc should cite it *and* confront it: our differentiator is
precisely triggering on the soft structure *without* the annihilation dijet — state that as the
explicit gap we fill, or concede the hard-object path is the realistic near-term trigger and we buy
back the ISR-less/annihilation-less phase space. Either way, engaging 2506.11192 sharpens the thesis.

**E.2 — vertex dx/dy transverse-position significance (e87538e / the just-prototyped kernel tool) as
a SECOND SUEP/dark-vertex handle beyond n_tracks.** The draft mentions the dx/dy tooling only for
"sharpening the hardest-vertex definition." It is more than that: a SUEP PV, being isotropic and
high-multiplicity, gives a **well-constrained, high-significance transverse (dx,dy) fit** (many tracks
→ small covariance), whereas a soft PU vertex gives a poorly-constrained one. The **transverse-position
significance itself** (from the least-squares (dx,dy) normal equations already prototyped) is a cheap
per-vertex discriminant orthogonal to n_tracks and complementary to sphericity — it encodes "this
vertex is real and track-rich" without needing the pT-sum. Add it as a **third feature in the Rank-1
(and Rank-2) feature vector**: (sphericity S, track-HT, dx/dy significance [, n_tracks as assist]).
This is a genuine, infra-backed addition the draft missed.

**E.3 — Dark-photon / dilepton soft-resonance handle (missing signature).** The draft covers dark
*showers* but omits the **dark-photon → soft dilepton** signature (low-mass A' → μμ/ee), which is a
mainstream dark-sector target. L1 handle: a **low-pT, small-opening-angle track pair from a common
displaced (or prompt) vertex**, potentially with the muon system for the μμ mode. This rides the same
extended-track + vertex infra (Rank 3-adjacent) and the L1 muon track finder. Worth a short §2.x
entry: it is a lower-multiplicity, higher-purity companion to the high-multiplicity SUEP handle and is
arguably *more* L1-tractable (two objects with a resonance constraint beats a soft pile). Feasibility
(a/b): dilepton-vertex association is expressible; a soft-dimuon L1 seed may need GMT+GTT correlation.

**E.4 — Does the refit d0/z0 improvement materially help any of these? HONEST ANSWER: marginally, and
only for Rank 3/4, not Rank 1/2.** From `mem:smartpixels-phase3-progress`: refit gives ~20% z0
resolution improvement (0.192→0.153 cm) and *creates* a d0 measurement (core 0.026–0.051 cm) where
the 4-par baseline pins d0≡0. Implications:
  - **Rank 1/2 (SUEP, prompt, isotropy-driven):** refit d0/z0 helps **negligibly** — the discriminant
    is event-shape and multiplicity of *prompt* tracks; a 20% z0 sharpening tightens vertex
    membership slightly (helps the track→vertex association in B.1) but is not the physics lever.
    Do not oversell refit for SUEP.
  - **Rank 3/4 (displaced/disappearing):** refit d0 is **material** — the baseline has *no* d0, so the
    refit d0 (or the extended-track nPar=5 d0) is the enabling measurement for any displacement
    handle. This is where the refit investment pays into the dark-sector program. State it precisely:
    **refit helps the displaced signatures and barely touches the isotropic ones.**

**E.5 — Timing / MTD as a displaced handle (scope question).** The draft omits track/vertex **timing**.
If MTD or L1 track-timing is in scope, a displaced/long-lived vertex carries a **time offset** vs the
prompt PV — a powerful, low-background displaced discriminant orthogonal to d0. Honest caveat: MTD
timing at L1 in Phase-2 is itself a design-gated capability and may be out of scope for a
track-trigger-anchored study. Recommend a one-line "out of current scope but flagged as the natural
next orthogonal handle for Rank 3/4" — do not build it into a rank, but do not silently omit it.

**E.6 — Sample-availability triage (the honest "what won't get produced soon" question).** The draft's
§6 lists the samples but does not *triage by realism*. My read:
  - **Realistically producible / partially in hand:** 200-PU min-bias/Zero-Bias background envelope
    (most reusable, standard) — this is the one procurement that unblocks Rank 1/2 *and* anomaly
    training. Push it first; it is not exotic.
  - **Producible with effort, standard generators:** SUEP signal MC (public SUEP generator exists;
    used by 2403.05311/2607.09621) and emerging-jet cτ scans (Pythia HV module). Medium lift.
  - **Hard / unlikely soon:** quirk/folded-SUSY MC (specialized generators, non-helical propagation
    in the sim — 2410.00269 built custom samples; this is a real production project, not a config).
    Compressed-EWK/disappearing-chargino is producible but the *disappearing* reconstruction needs the
    §4 object anyway, so the sample is not the binding constraint.
  - **Triage verdict:** Ranks 1–3 are sample-limited but the samples are *feasible*; Ranks 4–5 are
    **double-gated** (need both an exotic sample *and* the §4 subsystem), so they are honestly
    "feasibility-study / white-paper" items, not near-term prototypes. The draft should say this
    explicitly so nobody scopes a quirk trigger before the offline condensation study exists.

---

## (F) Ranked change-list (most-impactful first; mechanical for the author)

1. **Fix Rank 1 (the make-or-break).** Demote raw n_tracks; promote **sphericity/isotropy + track-HT
   (scalar sum of PV-track pT)** to the *primary* discriminant, n_tracks to *assist*. Insert the A.1
   per-regime table (high-m_S tractable / mid marginal / low-m_S & high-T_D out of reach at pT>2 GeV).
   State the pT-floor reality: separability requires pushing to ≤1 GeV, below the GTT baseline. This
   is the single most important edit.
2. **Fix the citation mis-attribution (A.4):** move the "pT>2 GeV multiplicity is inefficient for
   most SUEP" claim from 2607.09621 to **2211.05720 / 2203.07314**; re-describe 2607.09621 as the
   *HT-based scouting* search (1 TeV→410 GeV) with offline (n_const>50, S>0.5) selection.
3. **Add the quirk↔SUEP unification (E.1) and confront its trigger caveat.** Cite 2506.11192; state
   that quirks *are* a SUEP source (unifying §2.1/§2.6/§4), and explicitly address that the paper's
   own recommended trigger is the annihilation dijet — declaring the soft-structure trigger as *our*
   differentiator (or conceding the hard-object path). Add the 2509.08878 non-helical follow-up as
   further §4 prior art.
4. **Correct the convergence claim in §4/§D (D above):** reframe from "one condensation layer serves
   three signatures" to **"one subsystem (placement), two capabilities — object-condensation hit
   clustering (quirks + smarttrack/OT merge + disappearing-track formation) and regional low-threshold
   reseeding (SUEP-soft + displaced)."** This is the version that survives a track-finder-group read.
5. **Rewrite §D remit to the settled BOTH-tracks + upstream-payoff frame (D above):** Track 1 =
   near-term observables upstreamable to baseline Phase-2 GTT now (value survives SmartPixels
   non-adoption); Track 2 = the subsystem as a standalone Phase-2 track-trigger contribution decoupled
   from SmartPixels adoption risk (de-risking argument). Lead §4 with the non-adoption-survival
   sentence.
6. **Add the hidden Rank-1 dependencies (B.1):** state that the offline study needs (i) a sphericity
   implementation, (ii) confirmed **track→vertex membership in nano** (histogram counts alone are
   insufficient for the tensor — verify track z0 + vertex z are both persisted), (iii) a truth-PV vs
   reco-hardest-PV decision tied to the truth-posture discipline. Re-label Rank 1 "(a) offline with a
   small nano/producer dependency," not clean "(a)."
7. **Reclassify disappearing tracks (A.3)** as a track-FORMATION problem co-located with quirks
   (both §4-gated); stop implying `layerHitMask` is the enabling observable (it is downstream of the
   object that does not yet exist). Update the §2.7 table row accordingly.
8. **Add vertex dx/dy significance (E.2)** as an explicit third Rank-1/Rank-2 feature (isotropy S,
   track-HT, dx/dy significance, n_tracks-assist), citing the just-prototyped kernel tooling.
9. **Sharpen the quirk-GNN caveats (A.2):** headline **~72% (with background, 8-layer)** not 91%;
   add "offline only, no trigger study, PU 5–10 not 200, efficiency depends on oscillation length."
10. **Elevate the `setTrackWordBits` nPar=5 refit-d0 re-encode (B, Rank 3)** from open-question #3 to
    a stated **precondition** of Ranks 3/4 (silent d0 zeroing nulls the study).
11. **Add E.4 honesty on refit reach:** refit d0/z0 materially helps only displaced/disappearing
    (Rank 3/4); negligible for isotropic SUEP (Rank 1/2). Don't oversell refit for the flagship.
12. **Add the dark-photon soft-dilepton signature (E.3)** as a short §2.x — lower-multiplicity,
    higher-purity, possibly more L1-tractable companion handle.
13. **Add sample-triage realism (E.6):** mark 200-PU min-bias as the reusable unblocker to procure
    first; label Ranks 4–5 as double-gated feasibility items, not near-term prototypes.
14. **Fix the latency frame (A.4):** cite **12.5 µs total Phase-2 L1 latency**, with GTT/track as an
    allocation inside it (the "~4 µs" reads as the whole budget otherwise).
15. **(Optional) Add timing/MTD (E.5)** as a flagged out-of-scope orthogonal displaced handle.

---

## Reviewer's residual gaps (could not fully close)

- **GTT per-vertex accumulator firmware budget (C.3):** could not locate a primary GTT firmware note
  quantifying whether a per-z-bin sphericity tensor (6 extra sum-of-products) fits the FastHisto pass.
  Routed as a concrete question to the GTT firmware group in §6. Arithmetic prior is "cheap/plausible"
  but assert nothing.
- **Exact GTT extended-track d0 ceiling** ("~10 cm", §2.2): consistent with LST/HLT displaced
  literature but not pinned to a GTT primary source; marked needs-checking.
- **CICADA/AXOL1TL exact latency figures:** taken as cited (not independently re-verified); low risk.

Sources newly verified in this pass: arXiv:2211.05720, 2203.07314, 2607.09621, 2410.00269,
2506.11192, 2509.08878, 2601.17544; Phase-2 L1 12.5 µs / VU13P GTT platform facts.
