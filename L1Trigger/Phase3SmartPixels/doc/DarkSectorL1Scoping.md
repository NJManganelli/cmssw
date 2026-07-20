# Dark-Sector & Soft-BSM L1 Sensitivity — Physics Scoping + Feasibility (v1)

Status: research/scoping design doc, uncommitted. No code, no build, no samples yet.
Version note: v1 incorporates the friendly-critic pass (`DarkSectorL1Scoping-Critique.md`). The
load-bearing change from v0 is the **SUEP handle inversion** (§2.1): at the GTT pT>2 GeV baseline,
track *multiplicity* is near-useless for SUEP; **sphericity/isotropy + track-HT** are the survivors.
Other major v1 changes: quirk↔SUEP unification (§2.6), the two-capability subsystem reframe (§4), the
BOTH-tracks upstream-payoff remit (§D), and the additive handles (vertex dx/dy significance,
dark-photon dileptons, timing) in §2 and §5.

Grounds every recommendation in the SmartPixels/GTT/refit/nano/ngtagger stack (see
`mem:smartpixels-refit-sidecar-architecture`, `mem:smartpixels-phase3-progress`,
`mem:smartpixels-tq-mva-usage-map`, `mem:smartpixels-modelspace-findings`,
`mem:smartpixels-darksector-study-plan`, `mem:smartpixels-next-production-plan`).

Feasibility buckets: **(a)** computable from existing/near-term nano now; **(b)** needs a new
GTT/correlator observable or producer; **(c)** fundamentally hard.

---

## 1. Framing — the sensitivity gap and why it is track/vertex-shaped

### 1.1 The unifying problem

Dark-sector and soft-BSM signatures share one failure mode against the current trigger: they are
**soft and/or lack a hard trigger object**. Their visible energy falls under HT/MET/jet/lepton
thresholds, or — worse for the isotropic cases — it is smeared across the whole detector so no single
object clears a threshold and the event is indistinguishable from pileup at the object level.

The community's workaround is the **associated-production tag**: require the dark sector to recoil
against a hard SM object (ISR jet, leptonic Z/W, VBF jets) so *something* fires a conventional
trigger. This works but is expensive: it throws away the bulk of the cross-section (inclusive
production dwarfs associated production) and biases the kinematics toward the boosted tail, where
several of these signatures are least characteristic (SUEP isotropy is a rest-frame property;
boosting it costs the very handle you want). The existing CMS SUEP search is the textbook example —
it relies on high-threshold hadronic triggers for the boosted topology
([CMS PRL 133, 191902](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.133.191902);
[arXiv:2403.05311](https://arxiv.org/abs/2403.05311)); the follow-up **data-scouting** search
([arXiv:2607.09621](https://arxiv.org/html/2607.09621)) exists precisely to claw back the hadronic
**HT** threshold (1 TeV → 410 GeV) that the associated-activity requirement imposes — note it is an
HT trigger, not a track-multiplicity trigger (this distinction is load-bearing in §2.1).

The **opportunity**: for Phase-2, per-event L1 tracks and vertices exist at 40 MHz. The information
these signatures carry — event shape, track-HT, impact parameter, per-vertex structure, inner-layer
hit patterns — lives exactly at the layer this SmartPixels/GTT effort instruments. The right place to
catch soft/isotropic/displaced physics is a **track- or vertex-based, signature-agnostic observable
in GTT / the L1 or L2 correlator / the Global Trigger**, *before* the calorimeter-driven thresholds
decide the event is boring.

### 1.2 The existing CMS L1 anomaly triggers — what they cover, what they miss

Two model-independent anomaly triggers are already deployed/prototyped, and any track/vertex proposal
must be positioned as **complementary**:

- **CICADA** (Calorimeter Image Convolutional Anomaly Detection): convolutional autoencoder on the
  **calorimeter region energy map** (18×14), trained on Zero Bias, scoring reconstruction error.
  Teacher (~300k params) distilled to ~10k for sub-200 ns FPGA inference
  ([arXiv:2411.19506](https://arxiv.org/html/2411.19506v1);
  [arXiv:2510.15672](https://arxiv.org/abs/2510.15672)).
- **AXOL1TL**: VAE on **standard L1 objects** (pT^miss, 4 e/γ, 4 μ, 10 jets as (pT,η,φ)),
  latency < 100 ns, deployed since 2024 ([arXiv:2411.19506](https://arxiv.org/pdf/2411.19506);
  [CDS 2942560](https://cds.cern.ch/record/2942560)).

Both are **calorimeter/object-level**. Their common blind spot is exactly our target class:

- **Isotropic soft physics (SUEP)** produces no anomalous *object* — energy spreads into many
  sub-threshold tracks that never become jets/e/γ/μ. To a calo image it resembles diffuse pileup; to
  the object VAE there is nothing to feed. The discriminating structure (isotropy of PV-associated
  tracks) is **only visible in the track/vertex layer**.
- **Displaced / disappearing** signatures are defined by tracker geometry (impact parameter,
  inner-layer-only hit patterns) that never reaches the calorimeter representation.
- **Non-helical (quirk)** tracks are a tracker-level pathology invisible to calo/object VAEs.

A **track/vertex-based anomaly or novel-observable trigger is therefore a genuine third axis**:
CICADA covers calo images, AXOL1TL covers object kinematics, and a GTT-level track/vertex handle
covers *charged-particle event-shape, track-HT, displacement, and hit-pattern* structure the other
two are architecturally blind to. That is the strategic pitch.

---

## 2. Signature-by-signature map (core)

### 2.1 SUEP (Soft Unclustered Energy Patterns) — the flagship, re-anchored on isotropy + track-HT

**Physics.** Strongly-coupled hidden valley with large 't Hooft coupling: a scalar mediator (mass
m_S) decays into a dark shower that hadronizes into **many low-pT, near-isotropic** dark hadrons
decaying back to SM ([CMS PRL 133, 191902](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.133.191902);
Hidden Valley, Strassler-Zurek). The charged-pT spectrum **peaks sharply in the lowest pT bins**,
with most tracks **below 1 GeV**, and the spectrum is set by the dark-meson mass/temperature T_D,
**not** the mediator mass ([arXiv:2211.05720](https://arxiv.org/abs/2211.05720)).

**Why L1 misses it.** No hard jet/MET/lepton; visible pT smeared across many sub-threshold tracks.
Boosted-topology triggers (the current handle) sacrifice the isotropy and most of the cross-section.

**The refuted handle, and the correction.** The v0 draft led with *track multiplicity at the hardest
vertex*. This is **wrong at the GTT baseline** and is contradicted by the very papers v0 cited.
Verbatim ([arXiv:2211.05720](https://arxiv.org/abs/2211.05720)): "for a track-trigger similar to the
CMS HL-LHC upgrade, where pT>2 GeV is the baseline threshold, a trigger requiring a large track
multiplicity is inefficient for most SUEP signal models." Quantified:
- at **pT>2 GeV**: negligible efficiency for all signals at any nTrack threshold;
- at **pT>1 GeV**: 50% efficiency only for m_S ≳ 600 GeV;
- at **pT>0.5 GeV**: 50% only for m_S ≳ 200 GeV with nTrack>100.

The reason multiplicity dies is mechanical: the SUEP charged spectrum sits mostly below the GTT
2 GeV floor, so *few SUEP tracks survive at the PV to be counted*, while the hardest PU vertex
already carries a genuine hard-scatter's tracks — the per-vertex n_tracks *excess* is a losing tail
competition. What breaks the degeneracy is **event shape**: a hard-scatter PV is jetty (sphericity
S → 0), a SUEP PV is isotropic (S → 1). So the ordering inverts:

**Primary discriminants: (1) sphericity/isotropy S of PV-associated tracks; (2) track-HT (scalar
sum of PV-track pT, HT_trk). Assist: n_tracks.** The whitepaper
([arXiv:2203.07314](https://arxiv.org/pdf/2203.07314)) reaches the same conclusion (for SUEP, d0 is
irrelevant — all prompt — and event-shape/track-HT beat bare multiplicity), and the CMS scouting
search uses boosted-sphericity S > 0.5 (with n_const > 50) as its offline SR discriminant
([arXiv:2607.09621](https://arxiv.org/html/2607.09621)).

**Per-regime L1-reach verdict (the honest ceiling):**

| SUEP regime | L1 tractability with a GTT track handle |
|---|---|
| High m_S (≳600 GeV–1 TeV), moderate T_D/m_φ | **Tractable** — hardest tracks poke above pT>1–2 GeV; sphericity + track-HT separable. This is the regime we can honestly claim. |
| Mid m_S (200–600 GeV) | **Marginal** — needs the pT floor pushed to ≤1 GeV (below the GTT baseline) to keep 50% efficiency; multiplicity alone insufficient. |
| Low m_S and/or high T_D, m_φ ≲ 2 GeV (soft, "more unclustered") | **Not L1-tractable via GTT tracks at pT>2 GeV.** This is exactly the SUEP the scouting search bought back with an HT threshold; the GTT pT floor guts it. L1 does not recover it without a sub-2-GeV track path or an HT/anomaly path. |

Net: the physics *is* track/vertex-shaped, but only via **isotropy + track-HT**, and only the
high-m_S regime is honestly in reach at the GTT baseline. Where multiplicity does contribute, it is
an assist to the shape variables, not the lead.

**Where it computes.** GTT, in/adjacent to the vertex producer (per-z-bin structure is built there).

**Our infra.** FastHisto/NNVtx per-bin `n_tracks` (ngtagger `e87538e`); the vertex dx/dy + kernel
tooling (peak-finder kernels, two-close-vertex scan) both for sharpening the hardest-vertex
definition *and* as a third discriminant (see §2.1a). ngtagger as the observable/MVA vehicle.

**Feasibility: (a) offline, but with a small (b)-flavored producer/nano dependency** — *not* clean
(a). The reach study additionally needs three items the nano does not carry cleanly today:
1. **A sphericity/isotropy implementation over PV-associated tracks** — new code. Since the
   correction above makes sphericity the *primary* discriminant, this is not optional.
2. **Track→vertex membership** — the FastHisto histogram gives per-z-bin *counts*, but the sphericity
   tensor needs the **(px,py,pz) of the individual tracks in the hardest bin**. Confirm the nano
   exposes track z0 *and* vertex z so membership is reconstructable; if only bin counts are
   persisted, the tensor cannot be formed. **Rank 2 (anomaly score) inherits this same dependency.**
3. **Truth-PV vs reco-hardest-PV decision** — at 200 PU "hardest vertex" is an algorithmic argmax
   that may not be the SUEP vertex; for soft SUEP they diverge. This ties to the truth-posture
   discipline (`mem:smartpixels-*`: everything-from-file OR everything-in-job, never mix — a stale
   association map silently passes tracks through unmodified and would invalidate the per-vertex
   truth match).

*This is the flagship and the recommended first prototype (§5), in its sphericity + track-HT form.*

### 2.1a Vertex transverse-position (dx/dy) significance — a third SUEP/dark-vertex handle

The just-prototyped least-squares (dx,dy) from d0/φ normal equations in the FastHisto parallel
histograms (`e87538e` + kernel tooling, `mem:smartpixels-tq-mva-usage-map`) yields a per-vertex
**transverse-position significance**. A SUEP PV, being high-multiplicity, gives a **well-constrained
(dx,dy) fit** (many tracks → small covariance → high significance); a soft PU vertex gives a
poorly-constrained one. This is a cheap per-vertex discriminant **orthogonal to n_tracks and
complementary to sphericity** ("this vertex is real and track-rich" without needing the pT-sum).
Add it as a third feature: **(S, HT_trk, dx/dy-significance [, n_tracks-assist])** for Ranks 1 and 2.
Feasibility: (a/b), same tooling as the sphericity handle.

### 2.2 Emerging jets (displaced dark-pion decays)

**Physics.** Dark pions with cτ ~ mm–m decay along a jet, producing **many displaced vertices, one
per dark meson**, high-d0 tracks, prompt-poor cores ([CMS arXiv:2510.12347](https://arxiv.org/pdf/2510.12347);
[ATLAS IOP 2025](https://iopscience.iop.org/article/10.1088/1361-6633/adfe17); theory
Schwaller-Stolarski-Weiler). Offline taggers use **median 2D impact parameter** and **prompt-track
energy fraction**.

**Why L1 misses it.** Prompt-only finding drops the displaced tracks that *define* the jet; the core
looks empty; current searches require an accompanying hard jet.

**The observable.** (i) **Extended (5-parameter) L1 tracks with large |d0|** clustered in a jet
region; (ii) **displaced-vertex multiplicity along the jet** from the L1 DisplacedVertexFinder. The
Phase-2 track trigger reconstructs displaced tracks and a Displaced Vertex Track Trigger prototype
exists ([PoS LHCP2024 274](https://pos.sissa.it/478/274/pdf); [CDS 2901311](https://cds.cern.ch/record/2901311)).
The whitepaper's recommendation of **increasing d0 coverage with pT** directly supports this handle
([arXiv:2203.07314](https://arxiv.org/pdf/2203.07314)).

**Where it computes.** GTT (extended tracks + DisplacedVertexFinder) → L2 correlator.

**Our infra.** Extended-track d0 (`seedCovMode=trackCov`, `extendedTracks=True`,
`mem:smartpixels-next-production-plan`); **dispvtx retraining** in ngtagger (`train-dispvtx`).
**PRECONDITION (elevated from a footnote):** verify `setTrackWordBits` re-encodes refit d0 into the
extended word for nPar=5 (spec §7-D). If d0 is silently zeroed, the entire displaced study is null.
**Feasibility: (a/b)** contingent on the d0-word precondition.

*Exact GTT extended-track d0 ceiling: **needs-checking.** The "~10 cm" figure is consistent with the
LST/HLT displaced literature (r_vertex ≳ 5 cm acceptance) but I could not pin it to a primary GTT
firmware note — treat as plausible, not pinned; route to the GTT/DisplacedVertexFinder firmware note.*

### 2.3 Semi-visible jets

**Physics.** Dark shower with a partially-invisible hadron fraction → a jet **collinear with MET**:
small Δφ(MET, jet), the *inverse* of the standard QCD MET-cleaning cut
([Cohen-Lisanti-Lou arXiv:1503.00009](https://arxiv.org/abs/1503.00009v1)). These events are
**actively vetoed** by the QCD-style Δφ(MET,jet) cut; total MET is often modest.

**The observable.** **min Δφ(MET, jet)** as a *signal* variable, plus MT(jet, MET). A correlation of
existing L1 objects, not a new tracking primitive.

**Where it computes.** L2 correlator / Global Trigger (needs L1 PUPPI MET + jets together).

**Our infra.** Least direct — object-level, partly overlapping AXOL1TL's domain. Our contribution is
improving the *inputs* (Puppi MET quality via better vertex/track association), not the observable.
**Feasibility: (a)** as a GT topological condition; low priority *for us* (thin infra overlap).

### 2.4 QCD-like / prompt visible dark jets (the 4th quadrant)

**Physics.** Frame via the **dark-shower parameter plane**: (dark-hadron lifetime) × (invisible /
detector-stable fraction). Quadrants: prompt+visible (QCD-like), prompt+invisible (semi-visible/MET),
displaced+visible (emerging jets), displaced+invisible (mixed). The prompt-visible "4th quadrant"
dark jets look almost exactly like QCD
([Cohen-Lisanti-Lou arXiv:1707.05326](https://arxiv.org/pdf/1707.05326);
[Bernreuther et al. arXiv:2301.07732](https://arxiv.org/pdf/2301.07732);
[Snowmass dark-showers arXiv:2203.09503](https://arxiv.org/abs/2203.09503)).

**Why L1 misses it — and why it may be irreducible at L1.** They differ from QCD only in **subtle
substructure** (dark-hadron multiplicity, Lund-plane radiation), below L1 granularity. This is a
*jet-tagger* problem (our NG tagger's home), not a trigger-object problem. At L1 the realistic move
is **keep these events via an inclusive track/HT/anomaly path** and separate at HLT/offline.

**Feasibility: (c) as a dedicated L1 tag** (irreducible at L1 granularity); **(a)** only as
"route via anomaly/inclusive-HT, don't throw it away."

### 2.5 Compressed SUSY (small ΔM electroweakinos) — disappearing tracks (a track-FORMATION problem)

**Physics.** Near-degenerate binos/winos/higgsinos: decay products **ultra-soft** (ΔM ~ 0.3–few GeV).
The associated-production-free handle is the **disappearing track**: a long-lived chargino traverses
the inner tracker then decays to an invisible neutralino + an undetectable ultra-soft pion, leaving
**inner pixel hits with no outer-tracker continuation** ([CMS arXiv:1804.07321](https://arxiv.org/pdf/1804.07321);
[compressed-EWK low-pT-track CDS 2960051](https://cds.cern.ch/record/2960051)).

**Why L1 misses it — sharpened.** Standard L1 tracking is **OT-stub-seeded**, so it **never forms a
track object with no OT continuation**. A disappearing track is invisible to the *finder*, not merely
un-tagged downstream. Crucially, `layerHitMask` in the refit sidecar is a property of a *formed*
refit track — **if no track is formed, there is no row to carry a mask.** So the sidecar does **not**
provide this observable; it is downstream of the missing piece. This is a **track-FORMATION problem,
co-located with quirks (§2.6), both gated on the §4 subsystem** — not a "bucket-b observable."

**The would-be observable (once the object exists).** Inner IT/SmartPixels hits (a clean short stub),
no OT continuation, isolated (low surrounding track/calo activity). SmartPixels *inner* hits are
exactly what a short-stub inner-seeded finder would need — so §4 is genuinely the right home; just
not an "observable" today.

**Our infra.** SmartPixels inner-layer hit patterns are the substrate for the §4 formation stage.
**Feasibility: (c), §4-gated** (closer to Rank 5 than the v0 table implied). Also soft-lepton and
soft-isolated-track handles ride the low-pT extended-track path (less novel, more tractable).

### 2.6 Quirks / squirks — non-helical tracks, AND a SUEP source

**Physics.** A quirk–antiquirk pair bound by a macroscopic IR flux tube (confinement scale ~ 100 eV–
keV) undergoes **non-helical, oscillatory motion** (tracks curve toward each other) plus anomalous
ionization ([Kang-Luty arXiv:0805.4642](https://arxiv.org/abs/0805.4642); folded SUSY Burdman et al.;
[arXiv:1708.02243](https://ar5iv.labs.arxiv.org/html/1708.02243)).

**The quirk↔SUEP unification (new in v1).** Quirk-antiquirk **de-excitation radiates hundreds of soft
isotropic pions** — i.e. **quirks are a SUEP source** (high multiplicity, transverse-sphericity
isotropy, tracks spatially closer than QCD) ([arXiv:2506.11192](https://arxiv.org/abs/2506.11192)).
This unifies §2.1 and §2.6 into one physics family and makes the §4 subsystem serve a *coherent*
target. **But confront the caveat honestly:** that same paper's *recommended* trigger is the **hard
resonant dijet from the final quirk-antiquirk annihilation**, with the soft tracks used *offline*
(their searches: cut-based, supervised, and a CATHODE weakly-supervised anomaly search on the soft
low-pT tracks *accompanying hard-scatter events*). That is exactly the **associated-production crutch
this whole document is escaping.** So the link is double-edged: it strengthens the "one subsystem, one
signature family" argument, *and* it is a literature data-point that the community's own near-term
answer is "trigger on the hard object." **Our explicit differentiator:** trigger on the soft isotropic
structure *without* the annihilation dijet, buying back the annihilation-less / ISR-less phase space —
via the isotropy + track-HT handle (§2.1) at the GTT, and (for the non-helical tracks themselves) via
§4. State this as the gap we fill; if it proves infeasible at the GTT baseline, we concede the
hard-object path is the near-term trigger and we contribute only the offline-enhancement piece.

**Why L1 misses the tracks — plainly.** **L1 track finding assumes a helix.** Stub-based
road/Hough/KF finders are built on the helical hypothesis; a quirk track violates it at the
trajectory level, so the hits are **missed or mis-reconstructed *before any MVA or trigger observable
ever sees them***. No downstream tagging helps if the track was never formed. Two offline studies
quantify this and set the *honest* operating point (v1 corrects v0's optimistic "56–91%"):
- [arXiv:2410.00269](https://arxiv.org/html/2410.00269) (Exa.TrkX GNN, 8-layer): SM tracks 97.9%;
  quirks with the **SM-trained helical** pipeline **10.2%**. Quirk-aware GNN: 91.5% (well-behaved,
  **no background**), 56.3% (all quirks, no background), **71.9% (all quirks *with* SM background)** —
  cite **~72% (with background)** as the honest headline, not the 91% background-free ceiling.
- [arXiv:2509.08878](https://arxiv.org/html/2509.08878) (model-agnostic non-helical finder, defines
  trajectories implicitly in training, generalizes beyond it): **32% efficiency with 0 fakes** on
  quirks — a stricter, more honest operating point.

**Caveats the numbers hide:** **offline only** (no L1/trigger implementation — §4's L1 framing is an
extrapolation, flagged); **pileup 5–10, not 200** (the paper notes ~5–10% degradation at 5–10 PU;
extrapolating to 200 PU is unwarranted and the 10.2%/72% gap could shrink); **efficiency depends
strongly on oscillation length** (short-oscillation quirks look nearly helical/easy, long ones are the
hard novel case — "56–91%" hides this axis).

**Where it computes.** *Before* GTT, at track *formation* — not addressable by any GTT/correlator
observable. Needs a non-helical formation stage (§4) or is lost.

**Our infra.** None of the current stack helps directly (all downstream of helical finding). A
hit-graph condensation stage (§4) over SmartPixels + OT hits is the only path. **Feasibility: (c)** —
fundamentally hard, and a **track-finder-group** problem, not an ngtagger prototype (see §6).

### 2.7 Dark photon → soft dilepton (new in v1)

**Physics.** A' → μμ / ee at low mass — a mainstream dark-sector target, lower-multiplicity and
higher-purity than the dark-shower cases.

**The observable.** A **low-pT, small-opening-angle track pair from a common (prompt or displaced)
vertex**, with a resonance (mass) constraint; the μμ mode adds the L1 muon system. This is arguably
*more* L1-tractable than SUEP — two objects with a resonance constraint beats a soft pile — and it is
a natural high-purity companion to the high-multiplicity SUEP handle.

**Our infra.** Rides the extended-track + vertex infra (Rank 3-adjacent) and the L1 muon track finder.
**Feasibility: (a/b)** — dilepton-vertex association is expressible; a soft-dimuon L1 seed may need
GMT+GTT correlation.

### 2.8 Summary table

| Signature | L1-available handle | Where | Our infra | Bucket |
|---|---|---|---|---|
| SUEP | **sphericity/isotropy + track-HT** at PV (dx/dy-signif; n_tracks assist) | GTT | FastHisto per-bin n_tracks, vertex dx/dy+kernel, ngtagger | **a-with-small-b** |
| Emerging jets | extended high-d0 + dispvtx multiplicity in jet | GTT→L2 | extended tracks, DisplacedVertexFinder, train-dispvtx | **a/b** (d0-word precond.) |
| Semi-visible jets | small Δφ(MET,jet) + MT | L2/GT | (thin — object level) | a |
| QCD-like dark jets | substructure (irreducible at L1) | HLT/offline | NG tagger; L1 anomaly gate | **c** at L1 |
| Compressed SUSY | disappearing track = **track FORMATION** (inner hits, no OT) | §4→GTT | SmartPixels inner hits (NOT layerHitMask alone) | **c, §4-gated** |
| Quirks | non-helical **track formation**; also a SUEP source | §4 (pre-GTT) | needs §4 condensation | **c** |
| Dark photon | soft dilepton, resonance+vertex constraint | GTT+GMT | extended tracks, muon finder | **a/b** |

---

## 3. MVA vs traditional — where does an MVA genuinely help?

Per handle and location, under a hard L1 constraint. **Latency frame (corrected):** the total
**Phase-2 L1 latency is 12.5 µs** (verified); the L1 track path is an *allocation inside it*
(~4 µs + ~1 µs transfer is the tracking allocation, not the whole budget). GTT sits within the 12.5 µs.
"MVA" here means a *small*, quantized, hls4ml-deployable model; related L1 ML demonstrates NN blocks
in ~50 ns / 2 clock ticks on the GTT platform ([arXiv:2203.12852](https://ar5iv.labs.arxiv.org/html/2203.12852)).

- **SUEP (isotropy + track-HT).** A **well-chosen observable + threshold beats an MVA** for the first
  cut — S and HT_trk at the PV are cheap, interpretable, and physics-explicit (isotropy). An MVA earns
  its place only in the **joint tail**: a tiny BDT/MLP on (S, HT_trk, dx/dy-significance, n_tracks)
  shapes the multi-D boundary against the PU fluctuation envelope better than rectangular cuts — a
  *second-order* gain. **Verdict: observable-first; MVA as a boundary-shaper.** Consistent with the
  ngtagger finding that *features >> architecture* at limited stats (`mem:smartpixels-modelspace-findings`).
- **Emerging jets / displaced.** Traditional: |d0| threshold × dispvtx count. MVA helps more (offline
  emerging-jet taggers are GNNs because displacement is multi-track-correlated); at L1 a compact
  dispvtx-quality MVA is plausible, but extended-d0 + dispvtx-multiplicity cuts are the honest first
  step.
- **Disappearing / quirks.** Pattern-recognition/formation problems (§2.5/2.6); MVA is intrinsic but
  the bottleneck is *forming the object at all* (§4), so classification MVA is premature.
- **The signature-agnostic angle (the real MVA opportunity).** The ngtagger **embedding /
  anomaly-detection** machinery (SupCon pretrain, DeepSet-over-constituents, `deepset_contrastive`,
  EBOPs-as-constraint; `mem:smartpixels-modelspace-findings`) is the natural vehicle for a
  **track/vertex-level unsupervised anomaly score** — an AE/VAE over per-vertex summary features
  (S, HT_trk, dx/dy-significance, pT-profile moments, dispvtx count), trained on Zero-Bias/min-bias,
  scoring reconstruction error. This is the **track-based sibling of CICADA/AXOL1TL**, arguably the
  highest-leverage MVA in the program because it is *signature-agnostic*: SUEP, prompt dark jets, and
  anything with anomalous charged structure in one trigger bit, on the third axis. Latency is fine —
  an AE over ~O(10) *features* (not raw tracks) is small; EBOPs-as-constraint gives a firmware-budget
  training path. **Strongest MVA-specific recommendation.**

**Bottom line.** For the concrete handles, a well-chosen observable + threshold is the right first
move. The MVA that changes the game is the **unsupervised track/vertex anomaly score**, not a
per-signature classifier.

---

## 4. NEW SUBSYSTEM — one placement, TWO capabilities

**Lead with the de-risking sentence (§D remit).** A track-condensation / reseeding layer between the
finders and GTT/GMT/correlator is a **standalone contribution to the baseline Phase-2 track trigger
whose value is decoupled from SmartPixels' own adoption risk**: even if SmartPixels never enters CMS
Phase-3 (or ATLAS), a non-helical/condensation reseed stage is a Phase-2 tracking contribution in its
own right. Design it as Phase-2-baseline-useful, with SmartPixels as an *enhancement on top*.

**Placement.** Between `{L1TrackFinder + L1SmartPixelsTrackFinder}` and `{GTT / GMT / L1Correlator}`.
This is the unique point holding **all tracks + smarttracks** (and, expanded, **their hits**).
Everything downstream sees finished helical tracks; everything upstream commits to helices per-region.
A layer that sees *hits + partial tracks together* can do things neither end can.

**The convergence claim, corrected.** v0 said "one condensation layer serves three signatures." That
over-claims — the three want **different things**, and only two share a mechanism. The honest statement
is **one subsystem (placement), TWO capabilities**:

1. **Object-condensation hit clustering** — instance-segment hits into objects with **no trajectory
   prior**. Serves: **quirk recovery** (the hard novel need — [arXiv:2410.00269](https://arxiv.org/html/2410.00269),
   [arXiv:2509.08878](https://arxiv.org/html/2509.08878)); **disappearing-track formation** (inner
   hits, no OT — §2.5); **SmartPixels-OT merge** (spatial instance-matching of inner smarttracks with
   OT tracks — object-condensation latent clustering expresses this). These three share the same
   "segment hits into objects without a fit" mechanism — the *real* convergence axis.
2. **Regional low-threshold reseeding** — regionally re-run an **(still essentially helical)** finder
   at a lower pT floor / looser prompt assumption in hit-dense pockets. Serves: **SUEP-soft** tracks
   (a SUEP track *is* helical — just soft and dropped by the pT>2 GeV / nstub cut) and **displaced**
   tracks (helical but high-d0). This does **not** need the non-helical machinery — it needs a cheaper
   regional refinder. Distinct capability, same placement.

So: **one subsystem, one home, two capabilities, 4+ signatures.** This is a *placement/architecture*
convergence, not a single-model convergence — and it **survives a track-finder-group read**, whereas
"one condensation layer does everything" invites the objection "reseeding isn't condensation."

**Prior art (cited).** Exa.TrkX / GNN4ITk hit-graph edge-classification, no helix assumed
([exatrkx.github.io](https://exatrkx.github.io/); [arXiv:2203.12852](https://ar5iv.labs.arxiv.org/html/2203.12852));
**Object Condensation** (Kieseler), one-stage grid-free multi-object reconstruction
([EPJC 80, 886](https://link.springer.com/article/10.1140/epjc/s10052-020-08461-2)), demonstrated for
**high-pileup tracking** ([arXiv:2312.03823](https://arxiv.org/abs/2312.03823)); **GravNet /
DeepJetCore** distance-weighted GNNs on **FPGA via hls4ml** ([arXiv:2008.03601](https://arxiv.org/pdf/2008.03601);
Belle II GNN-ETM L1 EM-calo trigger [arXiv:2602.15118](https://arxiv.org/pdf/2602.15118)) — object
condensation *has* run in a real L1 trigger, on calo; real-time FPGA graph building exists
([arXiv:2307.07289](https://arxiv.org/pdf/2307.07289)); regional LLP hit-selection for the L1 track
trigger ([arXiv:1907.09846](https://arxiv.org/pdf/1907.09846)).

**Feasibility and latency reality (flag speculation).** **Speculative at L1 hardware today.** Full
GNN/condensation tracking at 40 MHz within the L1 budget is unsolved; the deployed real-time
condensation examples are calo clusters, not full tracker hit graphs. The realistic near-term form is
**regional / on-demand** (only where an anomaly flag fires), not a global replacement — bounded
latency, exactly how the quirk papers imagine trigger use (distillation + accelerators,
[arXiv:2410.00269](https://arxiv.org/html/2410.00269)). **The doable, non-speculative first step is
the offline feasibility study:** build the hit-graph / condensation model in ngtagger over the
SmartPixels + OT hit collections we persist (sidecar per-crossing hit records; nano link tables), on
a quirk / disappearing-track MC sample, and measure recovery efficiency vs the helical baseline. It
reuses the DeepSet/graph machinery the tagger has.

**Verdict:** **promising for the offline feasibility study; speculative as L1 hardware; needs-X**,
where X = (a) a quirk/disappearing MC sample, (b) an offline condensation study in ngtagger proving
recoverable efficiency, (c) track-finder-group engagement for any hardware path. **Prototype the
offline study; do not scope hardware yet.** The de-risking frame makes even the offline study a
Phase-2-baseline contribution regardless of SmartPixels adoption.

---

## D. Strategic remit — BOTH tracks, upstream-payoff frame (settled)

Two named tracks, in explicit parallel, framed by the **upstream-payoff / de-risking** argument:

- **Track 1 — near-term observables on the EXISTING track/vertex collections.** Sphericity/track-HT +
  dx/dy-significance at the PV (Rank 1), extended-d0 + dispvtx (Rank 3), the unsupervised track/vertex
  anomaly score (Rank 2), the dark-photon dilepton handle (Rank 6). Sample-gated, prototypable in
  ngtagger now, and — crucially — **directly upstreamable to the baseline Phase-2 GTT/correlator**:
  they are observables on collections that exist in the *standard, non-SmartPixels* Phase-2 track
  trigger. **Their value does not depend on SmartPixels being adopted** — payoffs of the R&D line that
  survive SmartPixels non-adoption.
- **Track 2 — the NEW subsystem (§4).** Longer-horizon, design-gated, track-finder-group-owned. A
  **standalone contribution to the baseline Phase-2 track trigger, decoupled from SmartPixels'
  adoption risk** — a de-risking argument for the whole speculative R&D line.

Design deliverables to be Phase-2-baseline-useful; SmartPixels as enhancement on top. This is the
strongest strategic framing available and both §4 and the roadmap lead with it.

**Where does the refit d0/z0 improvement materially help? Honest triage.** From
`mem:smartpixels-phase3-progress`: refit gives ~20% z0 sharpening (0.192→0.153 cm) and *creates* a d0
measurement (core 0.026–0.051 cm) where the 4-par baseline pins d0≡0.
- **Rank 1/2 (SUEP, prompt, isotropy-driven): negligible.** The discriminant is event-shape/track-HT
  of *prompt* tracks; a 20% z0 sharpening tightens vertex membership slightly (helps §2.1 dependency
  #2) but is not the physics lever. **Do not oversell refit for the isotropic flagship.**
- **Rank 3/4 (displaced/disappearing): material.** The baseline has *no* d0, so the refit/extended-d0
  is the enabling measurement for any displacement handle. **This is where the refit investment pays
  into the dark-sector program.**

**Timing / MTD (flagged, out of current scope).** A displaced/long-lived vertex carries a **time
offset** vs the prompt PV — a powerful displaced discriminant orthogonal to d0. MTD/L1 track-timing at
L1 in Phase-2 is itself design-gated and likely out of scope for a track-trigger-anchored study.
Flagged as the natural next orthogonal handle for Ranks 3/4; not built into a rank.

---

## 5. Priority ranking + roadmap

Ranked by (reach gain) × (feasibility) × (uses-our-infra).

### Rank 1 — GTT SUEP isotropy + track-HT tag  ★ recommended first prototype

- **Reach.** Recovers *inclusive high-m_S* SUEP (drops the associated-hadronic-activity requirement).
  Honest ceiling: mid-m_S marginal, low-m_S/high-T_D out of reach at the GTT pT>2 GeV floor (§2.1).
- **Feasibility: (a) offline with a small (b) dependency** — needs a sphericity implementation +
  confirmed track→vertex membership in nano + a truth-PV/reco-PV decision (§2.1). The vertex per-bin
  n_tracks and kernel tooling exist; the shape variables do not.
- **Uses-our-infra.** Maximal: FastHisto/NNVtx per-bin n_tracks (`e87538e`), vertex dx/dy+kernel
  tooling (also feeds the dx/dy-significance feature), ngtagger.
- **Needs.** **SUEP signal MC** (m_S scan + T_D variations) at 200 PU; **200-PU min-bias-like
  background** for the PU per-vertex envelope (the reusable unblocker — procure first).
- **First prototype (ngtagger).** (i) hardest FastHisto z-bin; (ii) over its tracks compute the
  sphericity tensor eigenvalue combination S, track-HT (Σ pT), and dx/dy-significance, with n_tracks
  as assist; (iii) joint density (S, HT_trk[, dx/dy-signif]) for SUEP vs the hardest-PU-vertex
  envelope; (iv) rectangular + small 2–4-feature BDT; (v) reach = **signal efficiency vs PU-driven
  fake-trigger rate** (rate–efficiency curve) — the trigger analogue of a ROC, computable **without
  full trigger simulation** from the per-vertex observable distributions.
- **Dependencies.** SUEP + min-bias samples; sphericity code; track→vertex membership confirmation;
  truth-posture decision.

### Rank 2 — Track/vertex unsupervised anomaly score (CICADA/AXOL1TL third axis)

- **Reach.** Signature-agnostic: one bit catching SUEP + prompt dark jets + anything with anomalous
  charged structure. Highest *strategic* value; complements the deployed anomaly triggers.
- **Feasibility: (a/b).** Offline AE over per-vertex features is straightforward in ngtagger;
  **inherits Rank 1's track→vertex-membership dependency** (not just the feature extractor).
- **Uses-our-infra.** ngtagger anomaly/embedding machinery, the same per-vertex features as Rank 1.
- **Needs.** Min-bias 200-PU for training; SUEP (+ prompt-dark-jet) for evaluation. No dedicated
  signal for training (unsupervised).
- **First prototype.** VAE/AE over (S, HT_trk, dx/dy-signif, pT-profile moments, dispvtx count)
  trained on PU-only; anomaly = reconstruction error; signal-eff-at-rate on SUEP. **Build together
  with Rank 1** (shared extractor + shared membership dependency).

### Rank 3 — Emerging-jet / displaced high-d0 + dispvtx tag

- **Reach.** Recovers emerging jets without the accompanying-hard-jet trigger. Refit d0 is *material*
  here (§D triage).
- **Feasibility: (a/b),** contingent on the **`setTrackWordBits` nPar=5 refit-d0 precondition** (spec
  §7-D) — a precondition, not a footnote (silent d0 zeroing nulls the study).
- **Uses-our-infra.** Extended-track d0, DisplacedVertexFinder, `train-dispvtx`.
- **Needs.** **Emerging-jet MC** (cτ scan) at 200 PU; d0-word verification; GTT d0-ceiling check.
- **First prototype.** Per-jet-region median |d0| + dispvtx multiplicity; efficiency-vs-rate on
  emerging-jet MC vs QCD/PU.

### Rank 4 — Disappearing-track L1 handle (§4-gated)

- **Reach.** Unique (compressed SUSY without ISR+MET); leverages SmartPixels' defining capability.
- **Feasibility: (c), §4-gated** — the *object* (inner hits, no OT) does not exist in the OT-seeded
  flow; needs §4 formation. Do the §4 offline study first. **Double-gated** (exotic sample + §4).
- **Uses-our-infra.** SmartPixels inner hits (via §4), refit sidecar per-layer info *once a track is
  formed*.
- **Needs.** Compressed-EWK/disappearing MC; the §4 offline condensation study.

### Rank 5 — Quirk non-helical recovery (offline feasibility only)

- **Reach.** Opens a class currently *lost at reconstruction* (10.2% helical baseline → ~72% / 32%
  offline). Highest novelty, lowest near-term feasibility; unifies with SUEP (§2.6).
- **Feasibility: (c),** track-finder-group problem; hardware speculative; **double-gated** (quirk MC +
  §4). Offline only, low-PU-validated so far.
- **First step.** Offline hit-graph condensation study in ngtagger on a quirk MC sample; recovery vs
  helical baseline. A *feasibility* deliverable, not a trigger.
- **Needs.** **Quirk MC** (Kang-Luty / folded-SUSY generators, custom non-helical propagation —
  2410.00269 built custom samples; a real production project). Track-finder-group engagement (§6).

### Rank 6 — Dark-photon soft-dilepton tag

- **Reach.** Lower-multiplicity, higher-purity companion to SUEP; possibly *more* L1-tractable (two
  objects + resonance constraint).
- **Feasibility: (a/b).** Dilepton-vertex association expressible; soft-dimuon may need GMT+GTT
  correlation.
- **Uses-our-infra.** Extended tracks + vertex; L1 muon finder.
- **Needs.** Dark-photon (A' → μμ/ee) MC; a GMT+GTT correlation prototype for the dimuon mode.

**Single best first prototype:** **Rank 1, the GTT SUEP isotropy + track-HT tag** (with dx/dy
significance and n_tracks as assist). The vertex tooling exists, it needs the reusable 200-PU min-bias
sample plus one SUEP sample, and it attacks the associated-production limiter with a fully offline,
sample-in / reach-curve-out study. Its **now-fuller dependency list**: (1) SUEP MC + 200-PU min-bias;
(2) a sphericity/isotropy implementation over PV tracks; (3) confirmed track→vertex membership in nano
(track z0 + vertex z persisted); (4) a truth-PV vs reco-hardest-PV decision under the truth-posture
discipline; (5) the dx/dy-significance feature from the kernel tooling. Build the feature extractor so
**Rank 2 reuses it** — one build, two deliverables.

**Sample-availability triage (realism).**
- **Realistically producible / most reusable:** **200-PU min-bias / Zero-Bias background** — unblocks
  Rank 1/2 *and* anomaly training; standard, not exotic. **Procure first.**
- **Producible with effort, standard generators:** SUEP (public generator, used by
  2403.05311/2607.09621) and emerging-jet cτ scans (Pythia HV module). Medium lift.
- **Hard / unlikely soon (double-gated feasibility items, not near-term prototypes):**
  quirk/folded-SUSY MC (custom non-helical sim propagation) and disappearing-chargino (its *reco*
  needs the §4 object regardless). **Ranks 4–5 are feasibility / white-paper items** — do not scope a
  quirk trigger before the offline condensation study exists.

---

## 6. Open questions for collaborators

**For the L1 track-finder group (not ngtagger prototypes):**
1. **Quirk non-helical tracks (§2.6, §4):** does the Phase-2 stub-based finder have *any* acceptance
   for non-helical trajectories, or are quirk hits discarded at stub-forming? Room for a regional
   non-helical / condensation re-seed stage, and its latency envelope? The load-bearing question for
   the quirk program.
2. **Disappearing-track object (§2.5):** can an "inner IT/SmartPixels hits, no OT" object be formed in
   the L1 flow at all? What would feeding inner-tracker stubs into a short GTT-adjacent stage take?
3. **Extended-track d0 fidelity (Rank 3 precondition):** confirm `setTrackWordBits` preserves refit d0
   in the nPar=5 extended word end-to-end (spec §7-D).
4. **§4 subsystem latency:** is a regional/on-demand condensation stage between the finders and
   GTT/correlator admissible in the 12.5 µs Phase-2 budget, or strictly HLT/offline?

**For the GTT firmware group (open gap):**
5. **Per-vertex accumulator budget:** can the FastHisto pass accumulate **6 additional per-z-bin
   sum-of-products** for a **sphericity tensor** within the DSP/latency budget (VU13P platform)? A
   sphericity tensor is arithmetically cheap (multiply-accumulate), so the prior is "plausible," but
   assert nothing without the firmware group. This is the specific number the bucket-(b) SUEP producer
   scope needs.
6. **Extended-track d0 ceiling:** pin the exact GTT/DisplacedVertexFinder d0 acceptance (the "~10 cm"
   is plausible but unpinned to a primary firmware source).

**Sample-procurement list (gating dependency), in priority order:**
- **200-PU min-bias / Zero-Bias-equivalent** — single most reusable; unblocks Ranks 1/2 + anomaly
  training. **First.**
- **SUEP** signal MC, 200 PU, m_S scan + T_D variations (Rank 1/2).
- **Emerging-jet** MC, cτ scan, 200 PU (Rank 3).
- **Dark-photon** (A' → μμ/ee) MC (Rank 6).
- **Compressed-EWK / disappearing-chargino** MC (Rank 4) — double-gated with §4.
- **Quirk / folded-SUSY** MC (Rank 5) — hardest to generate (custom non-helical sim); flag early.
- (Lower) **prompt dark-jet** for anomaly-score stress test (Rank 2).

**Cross-cutting:** settle the truth-posture per-sample *before* production (everything-from-file OR
everything-in-job, never mix) — mixing silently passes tracks through unmodified and would invalidate
any per-vertex-truth SUEP study (Rank 1 dependency #4).

---

## Appendix: literature-pass status (v1)

Web research was available and used; all cited claims carry links. Sources verified this pass and in
the critique: arXiv:2211.05720, 2203.07314, 2607.09621, 2403.05311, 2410.00269, 2509.08878,
2506.11192, 2312.03823, 1503.00009, 1707.05326, 2301.07732, 2203.09503, 2510.12347, 1804.07321,
2411.19506, 2510.15672, 2008.03601, 2602.15118, 2307.07289, 1907.09846; Object Condensation EPJC
80, 886; Phase-2 L1 12.5 µs / VU13P GTT platform facts.

**Residual open gaps (routed, not closed):**
- **GTT per-vertex-accumulator firmware budget** — whether a per-z-bin sphericity tensor (6 extra
  sum-of-products) fits the FastHisto pass; routed to the GTT firmware group (§6.5). Arithmetic prior
  "plausible"; assert nothing.
- **Exact GTT extended-track d0 ceiling** ("~10 cm") — consistent with LST/HLT displaced literature,
  not pinned to a primary GTT source; needs-checking (§2.2, §6.6).
- **CICADA/AXOL1TL exact latencies** (sub-200 ns / <100 ns) — taken as cited, not independently
  re-verified; low risk.
