import FWCore.ParameterSet.Config as cms

# ---------------------------------------------------------------------------
# Phase-3 SmartPixels customisations
#
# Variant model: a SmartPixels track variant is (mode, activeSP), where mode is
# an emulator mode of L1SmartPixelsTrackProducer and activeSP (for the modes
# that use it) is the 4-bit pixel-layer configuration, e.g. "1100".
# Collection labels encode the variant: activeSP maps 1->A(ctive), 0->I(nactive),
# e.g. mode="correctionlibRegression", activeSP="1100" ->
#   l1tSmartPixelsTrackProducerWcorrectionlibRegressionAAII
# Future modes (e.g. a genuinely fitted L1SmartPixelsTracks collection) only
# need to be appended to the mode lists below; nothing structural assumes
# "regression"/"correctionlib".
#
# Two workflow entry points:
#  - smartPixelsCoexist (WF1): standard tracks AND SmartPixels variant tables
#    in one L1Nano file (in-file track-to-track comparison).
#  - smartPixelsCoopt (WF2, "pure view"): rewire every downstream consumer of
#    the standard tracklet tracks to ONE variant; compare nano files job-to-job.
# ---------------------------------------------------------------------------

# Modes that produce one collection per producer instance, no activeSP config
PASSTHROUGH_MODES = ["passthrough", "passthroughFloat", "passthroughHW", "trackingParticleTruth"]
# Modes parametrized by the active pixel-layer configuration (need activeSP).
# digiRefit reuses the activeSP encoding to select the refit layer set (BPIX
# layers), e.g. "1100" -> suffix AAII selects the two inner barrel layers.
ACTIVE_LAYER_MODES = ["correctionlibRegression", "correctionlibTPToySmear", "digiRefit"]
# Modes for which the TP/truth association is LOAD-BEARING: without it the
# producer silently passes tracks through unmodified (per-track fallback),
# which invalidates any study. Truth availability must be guaranteed.
# digiRefit additionally needs the pixel digis + PixelDigiSimLinks (posture B),
# not just the TTTrack association maps; the requirement is flagged the same way.
TRUTH_REQUIRED_MODES = ["trackingParticleTruth", "correctionlibRegression",
                        "correctionlibTPToySmear", "digiRefit"]
# Refit "tier model" (see L1Trigger/Phase3SmartPixels/README.md and
# mem:smartpixels-tier2-refit-plan):
#   digiRefit = Tier 2 interim refit (real pixel digis + synthesized angles).
#               Active-layer + truth-required; implementation lands in Phase 2.
#   refit     = Tier 3 true system (ingests a real SmartTracklet collection from
#               an L1SmartTracksFinder). RESERVED: accepted by the mode
#               vocabulary but building a variant with it raises NotImplementedError.
RESERVED_MODES = ["refit"]
# The complete accepted mode vocabulary (order = documentation order, not priority).
ALL_MODES = PASSTHROUGH_MODES + ACTIVE_LAYER_MODES + RESERVED_MODES

REGRESSION_CONFIGS = ["0000",
                      "0001", "0010", "0100", "1000",
                      "0011", "0101", "0110", "1001", "1010", "1100",
                      "0111", "1011", "1101", "1110",
                      "1111",
                      ]

DEFAULT_CORRECTION_SET = "L1Trigger/Phase3SmartPixels/data/spixel_smear_all_configs_barrel_CalV1_v2p1_compound.json"


# ---------------------------------------------------------------------------
# digiRefit (Tier 2) config surface — defined in full now (Phase 0); the
# producer that consumes it lands in Phase 2. Every entry is documented with
# units where relevant. Enum-valued entries are validated against the *_CHOICES
# sets below; unknown keys or out-of-vocabulary values raise loudly at config
# time (ValueError), so the surface is exercised even though nothing runs yet.
# ---------------------------------------------------------------------------
DIGIREFIT_USEANGLES_CHOICES = ("none", "alpha", "alphaBeta")
DIGIREFIT_GAINMODE_CHOICES = ("full", "lut")
DIGIREFIT_SEEDCOVMODE_CHOICES = ("trackCov", "parametrized")

DIGIREFIT_DEFAULTS = {
    # --- search windows (module-local frame, PER-LAYER TBPX L1-L4) ---
    # rphi spread of the beamline-constrained extrapolation GROWS outward
    # (MS-compensation bulge between origin anchor and OT stubs; PSimHit-
    # verified q68 0.04/0.16/0.49/0.90 cm at 2-5 GeV); z shrinks outward.
    # A scalar is accepted and broadcast to all four layers.
    "windowRPhi": (0.05, 0.17, 0.5, 0.9),  # r-phi half-widths [cm]
    "windowZ": (0.45, 0.35, 0.25, 0.2),    # z half-widths [cm]
    # --- track emission ---
    "minHits": 1,             # min attached IT hits required to emit a refit track
    # --- fidelity handle: which synthesized angle(s) enter the refit ---
    "useAngles": "alpha",     # "none" | "alpha" | "alphaBeta"
    # --- FPGA-fidelity / truncation handles (float impl, switchable) ---
    "maxHitsPerWindow": 8,    # combinatorics truncation: max in-window digis kept per layer
    "maxKFUpdates": 4,        # max Kalman updates applied (layer/update cap)
    "gainMode": "full",       # "full" (exact gain) | "lut" (RESERVED table-driven placeholder)
    # KF numerical guards (spec §6b, FPGA-fidelity handles). Defaults from the PU
    # chi2-tail investigation: physical |H|<~200 and the physical (incl. wrong-hit)
    # chi2-increment ceiling is ~1.9e6, so 1e4/2e6 kill only the non-physical
    # grazing-crossing tail (chi2 up to ~1e10) and leave the physical bulk untouched.
    "jacobianMaxAbs": 1.0e4,  # |H[k][j]| above this (or non-finite) zeroes that Jacobian column
    "chi2UpdateGate": 2.0e6,  # scalar update with r^2/S above this is skipped entirely
    # Grazing-angle clamps (spec v0.4 §6b). measAngleMaxAbs is LOAD-BEARING (a
    # synthesized measured cotAlpha/cotBeta beyond it invalidates that angle only,
    # removing the measurement-driven chi2 pathology at source); predAngleMaxAbs is
    # secondary hygiene (a predicted crossing angle beyond it invalidates the
    # crossing at the projector). Both default 12.0.
    "measAngleMaxAbs": 12.0,  # |synthesized measured cot| above this clears hasAlpha/hasBeta
    "predAngleMaxAbs": 12.0,  # |predicted crossing cot| above this invalidates the crossing
    # --- Kalman seed (user decision 2026-07-18: config-switchable) ---
    # CMSSW_17 backport (Branch B): TTTrack has NO helix covariance (Ian's covMatrix
    # PR cms-sw/cmssw#51503 is Branch C, not applied here), so the default seed
    # covariance is "parametrized". "trackCov" raises loudly (see _applyDigiRefitConfig).
    "seedNPar": 5,            # 4 | 5: seed-track parametrization entering the KF
    "seedCovMode": "parametrized",  # "parametrized" (default here) | "trackCov" (needs #51503 backport, Branch C)
    "paramSigmas": (1e-4, 1e-3, 2e-3, 0.06, 0.05),  # parametrized-mode sigmas (rInv[cm^-1],phi0,tanL,z0[cm],d0[cm])
    # --- correctionlib payload paths (empty defaults acceptable for Phase 0) ---
    # RESERVED: Stack A "smarthit_true" payload. Tier-2 does NOT consume it (position
    # comes from the real digis, angle from the PixelAV response; Stack A only
    # CHARACTERIZES true hits). The key is kept (present + validated) for a future
    # SmartPixels ASIC on-chip readout-inefficiency model (smarthit_true_eff), so
    # wiring is drop-in. Shipped recipes pass no path; if set, the producer warns
    # (SmartPixelsStackAUnused) and loads nothing.
    "smarthitTrueSet": "",    # RESERVED (see note above) -- do not ship a real path
    "smarthitFakeSet": "",    # Stack B "smarthit_fake" payload (window multiplicity / fakes)
    "pixelavAngleSet": "",    # PixelAV angle sigma/bias response payload
    # --- optional refit-aware TkQuality BDT ---
    "bdtModel": "",           # optional path to the refit TkQuality BDT model (empty = none)
}


def _resolveDigiRefitConfig(digiRefitConfig=None):
  """Merge a user digiRefitConfig dict over DIGIREFIT_DEFAULTS with loud validation.

  Returns a NEW dict (defaults never mutated). Unknown keys and out-of-vocabulary
  enum values raise ValueError at config time (Phase 0 requirement: the surface
  is exercised and mistakes fail loudly, not silently). Passing None yields a
  fresh copy of the defaults.
  """
  resolved = dict(DIGIREFIT_DEFAULTS)
  if digiRefitConfig is None:
    return resolved
  if not isinstance(digiRefitConfig, dict):
    raise ValueError(f"digiRefitConfig must be a dict, got {type(digiRefitConfig).__name__}")
  unknown = set(digiRefitConfig) - set(DIGIREFIT_DEFAULTS)
  if unknown:
    raise ValueError(
        f"digiRefitConfig has unknown key(s) {sorted(unknown)}; "
        f"valid keys are {sorted(DIGIREFIT_DEFAULTS)}")
  resolved.update(digiRefitConfig)
  if resolved["useAngles"] not in DIGIREFIT_USEANGLES_CHOICES:
    raise ValueError(
        f"digiRefitConfig['useAngles']={resolved['useAngles']!r} invalid; "
        f"must be one of {DIGIREFIT_USEANGLES_CHOICES}")
  if resolved["gainMode"] not in DIGIREFIT_GAINMODE_CHOICES:
    raise ValueError(
        f"digiRefitConfig['gainMode']={resolved['gainMode']!r} invalid; "
        f"must be one of {DIGIREFIT_GAINMODE_CHOICES}")
  if resolved["seedCovMode"] not in DIGIREFIT_SEEDCOVMODE_CHOICES:
    raise ValueError(
        f"digiRefitConfig['seedCovMode']={resolved['seedCovMode']!r} invalid; "
        f"must be one of {DIGIREFIT_SEEDCOVMODE_CHOICES}")
  # CMSSW_17 backport (Branch B): TTTrack carries no helix covariance in this
  # release, so seedCovMode='trackCov' (which seeds the KF from helixCovMat())
  # cannot run. It is gated here at config time AND in the C++ constructor.
  if resolved["seedCovMode"] == "trackCov":
    raise ValueError(
        "digiRefitConfig['seedCovMode']='trackCov' seeds the refit from "
        "TTTrack::helixCovMat(), which does not exist in this release. The TTTrack "
        "helix covariance was added by cms-sw/cmssw PR #51503 (the SmartPixels "
        "'Branch C' backport), which is NOT applied here. Use seedCovMode='parametrized' "
        "(the default in CMSSW_17), or apply the #51503 covMatrix backport to enable trackCov.")
  if resolved["seedNPar"] not in (4, 5):
    raise ValueError(
        f"digiRefitConfig['seedNPar']={resolved['seedNPar']!r} invalid; must be 4 or 5")
  if len(tuple(resolved["paramSigmas"])) != 5:
    raise ValueError(
        f"digiRefitConfig['paramSigmas'] must have exactly 5 entries "
        f"(rInv, phi0, tanL, z0, d0), got {len(tuple(resolved['paramSigmas']))}")
  for wk in ("windowRPhi", "windowZ"):
    wv = resolved[wk]
    if isinstance(wv, (int, float)):
      wv = (float(wv),) * 4  # scalar convenience: broadcast to all layers
    wv = tuple(float(x) for x in wv)
    if len(wv) != 4 or any(x <= 0 for x in wv):
      raise ValueError(
          f"digiRefitConfig['{wk}'] must be a positive scalar or 4 per-layer "
          f"positive values (TBPX L1-L4), got {resolved[wk]!r}")
    resolved[wk] = wv
  # KF numerical guards + grazing-angle clamps must be positive doubles (spec §6b).
  for gk in ("jacobianMaxAbs", "chi2UpdateGate", "measAngleMaxAbs", "predAngleMaxAbs"):
    gv = resolved[gk]
    if not isinstance(gv, (int, float)) or isinstance(gv, bool) or gv <= 0:
      raise ValueError(
          f"digiRefitConfig['{gk}'] must be a positive number, got {resolved[gk]!r}")
    resolved[gk] = float(gv)
  # Refit-quality BDT model path (spec §6a): a string; "" disables scoring.
  # Feature-count validation (n_features == 17) happens producer-side at load.
  if not isinstance(resolved["bdtModel"], str):
    raise ValueError(
        f"digiRefitConfig['bdtModel'] must be a string path (\"\" = none), "
        f"got {resolved['bdtModel']!r}")
  return resolved


def _checkReservedMode(mode):
  """Raise NotImplementedError for RESERVED modes (accepted by the vocabulary
  but not buildable). Currently only 'refit' (Tier 3)."""
  if mode in RESERVED_MODES:
    raise NotImplementedError(
        f"SmartPixels mode '{mode}' is RESERVED for the true SmartTracklet-input "
        f"system (Tier 3): it will ingest a real SmartTracklet collection from an "
        f"L1SmartTracksFinder. Not implemented. Use mode 'digiRefit' for the "
        f"interim Tier 2 refit (real pixel digis + synthesized angles).")


def smartPixelsVariantSuffix(mode, activeSP=None):
  """Unique label/table suffix for a variant, e.g. 'passthrough' or 'correctionlibRegressionAAII'.

  digiRefit reuses the activeSP encoding (e.g. '1100' -> 'AAII'), so
  smartPixelsVariantSuffix('digiRefit', '1100') -> 'digiRefitAAII' and the
  derived table labels still end in 'Table' downstream (nano keep-pattern gotcha).
  """
  if mode not in ALL_MODES:
    raise ValueError(f"unknown SmartPixels mode '{mode}'; valid modes are {ALL_MODES}")
  _checkReservedMode(mode)
  if mode in ACTIVE_LAYER_MODES:
    if activeSP is None:
      raise ValueError(f"SmartPixels mode '{mode}' requires activeSP (e.g. '1100')")
    return mode + activeSP.replace("0", "I").replace("1", "A")
  return mode


def smartPixelsVariantLabels(mode, activeSP=None):
  """(promptLabel, extendedLabel) of the collections produced for a variant."""
  suffix = smartPixelsVariantSuffix(mode, activeSP)
  return (f"l1tSmartPixelsTrackProducerW{suffix}",
          f"l1tSmartPixelsTrackProducerExtendedW{suffix}")


def _normalizeVariants(variants):
  """Accept 'mode', ('mode', activeSP) or {'mode':..., 'activeSP':...}; return [(mode, activeSP)]."""
  normalized = []
  for variant in variants:
    if isinstance(variant, str):
      normalized.append((variant, None))
    elif isinstance(variant, dict):
      normalized.append((variant["mode"], variant.get("activeSP")))
    else:
      mode, activeSP = variant
      normalized.append((mode, activeSP))
  return normalized


def _allVariants():
  return ([(mode, None) for mode in PASSTHROUGH_MODES]
          + [("correctionlibRegression", rconf) for rconf in REGRESSION_CONFIGS])


def _applyDigiRefitConfig(module, resolved):
  """Push a resolved (validated) digiRefit config dict onto a producer module."""
  module.digiRefitWindowRPhi = cms.vdouble(*resolved["windowRPhi"])
  module.digiRefitWindowZ = cms.vdouble(*resolved["windowZ"])
  module.digiRefitMinHits = cms.int32(resolved["minHits"])
  module.digiRefitUseAngles = cms.string(resolved["useAngles"])
  module.digiRefitMaxHitsPerWindow = cms.int32(resolved["maxHitsPerWindow"])
  module.digiRefitMaxKFUpdates = cms.int32(resolved["maxKFUpdates"])
  module.digiRefitGainMode = cms.string(resolved["gainMode"])
  module.digiRefitJacobianMaxAbs = cms.double(resolved["jacobianMaxAbs"])
  module.digiRefitChi2UpdateGate = cms.double(resolved["chi2UpdateGate"])
  module.digiRefitMeasAngleMaxAbs = cms.double(resolved["measAngleMaxAbs"])
  module.digiRefitPredAngleMaxAbs = cms.double(resolved["predAngleMaxAbs"])
  module.digiRefitSeedNPar = cms.int32(resolved["seedNPar"])
  module.digiRefitSeedCovMode = cms.string(resolved["seedCovMode"])
  module.digiRefitParamSigmas = cms.vdouble(*resolved["paramSigmas"])
  module.digiRefitPixelavAngleSet = cms.string(resolved["pixelavAngleSet"])
  module.digiRefitSmarthitFakeSet = cms.string(resolved["smarthitFakeSet"])
  # RESERVED, not consumed by Tier-2: passed through so future ASIC-efficiency
  # wiring is drop-in. Non-empty triggers a LogWarning in the producer ctor.
  module.digiRefitSmarthitTrueSet = cms.string(resolved["smarthitTrueSet"])
  module.digiRefitBdtModel = cms.string(resolved["bdtModel"])


def addSmartPixelsTrackProducerVariants(process, variants=None, correctionSet=DEFAULT_CORRECTION_SET,
                                        digiRefitConfig=None):
  """STEP1: instantiate the SmartPixels track producer (prompt+extended) for each variant.

  variants=None reproduces the legacy behavior: all passthrough-family modes plus
  all 16 regression configs, produced concurrently (they are cheap smears/copies
  of the same base tracklet tracks).

  digiRefitConfig: a dict merged over DIGIREFIT_DEFAULTS (validated loudly; see
  _resolveDigiRefitConfig). Only relevant when a digiRefit variant is present.
  Phase 2: digiRefit variants are fully wired (params only); the producer seeds a
  local RNG per-event from hash(module label, run, lumi, event), so no
  RandomNumberGeneratorService is needed. pixelavAngleSet is REQUIRED (ValueError).

  Returns (process, moduleNames); scheduling the modules is up to the caller
  (see smartPixelsCoexist/smartPixelsCoopt for path association).
  """
  if variants is None:
    variants = _allVariants()
  variants = _normalizeVariants(variants)

  # Resolve+validate the digiRefit surface up front so a malformed dict fails
  # loudly at config time regardless of where the digiRefit variant sits.
  digiRefitResolved = _resolveDigiRefitConfig(digiRefitConfig)

  process.load('L1Trigger.Phase3SmartPixels.l1tSmartPixelsTrackProducer_cfi')
  modules = []
  for mode, activeSP in variants:
    # smartPixelsVariantLabels -> smartPixelsVariantSuffix rejects unknown modes
    # and raises the RESERVED (Tier 3 'refit') NotImplementedError here.
    prompt, extended = smartPixelsVariantLabels(mode, activeSP)

    if mode == "digiRefit" and not digiRefitResolved["pixelavAngleSet"]:
      raise ValueError(
          "digiRefit variant requested but digiRefitConfig['pixelavAngleSet'] is empty; "
          "the PixelAV angle-response payload is REQUIRED (see doc/PixelAVAngleResponseSpec.md; "
          "validatePixelAVAngleSet.py --write-example emits a structural stand-in for testing)")

    modules.append(prompt)
    modules.append(extended)

    setattr(process, prompt, process.l1tSmartPixelsTrackProducer.clone())
    setattr(process, extended, process.l1tSmartPixelsTrackProducerExtended.clone())
    for label in (prompt, extended):
      module = getattr(process, label)
      module.smartPixelsEmulatorMode = cms.string(mode)
      if mode in ACTIVE_LAYER_MODES:
        module.smartPixelsActiveLayers = activeSP
        module.smartPixelsCorrectionSet = cms.FileInPath(correctionSet)
      if mode == "digiRefit":
        _applyDigiRefitConfig(module, digiRefitResolved)
        # No RandomNumberGeneratorService: the digiRefit producer uses a LOCAL
        # engine seeded per-event from hash(module label, run, lumi, event), so
        # outputs are event-order-independent and split-job invariant (see
        # doc/Phase2Acceptance.md §1 and the producer's digiRefitSeed()).

  return process, modules


TRUTH_ASSOCIATOR_LABELS = ("TTClusterAssociatorFromPixelDigis",
                           "TTStubAssociatorFromPixelDigis",
                           "TTTrackAssociatorFromPixelDigis",
                           "TTTrackAssociatorFromPixelDigisExtended")


def useTruthAssociationFromFile(process, associatorLabels=TRUTH_ASSOCIATOR_LABELS):
  """Force the TT truth association maps to be read FROM THE INPUT FILE.

  The associators need mix:Tracker PixelDigiSimLinks, which are TRANSIENT DIGI
  products: they exist in-memory in a job that runs DIGI (central .774-style
  DIGI+L1+NANO jobs, or GEN-SIM-DIGI-RAW(-MINIAOD) inputs that retained them),
  but are dropped from most stripped tiers/RelVals. When the input file instead
  carries the READY-MADE association maps (STEP1-style output with
  'keep *_TTTrackAssociatorFromPixelDigis*_*_*' etc., as in the b-tagging
  runJetNTuple workflow), the in-process associator modules must be REMOVED so
  consumption resolves to the file products instead of re-running in-job.
  The file must also carry the TrackingParticles the maps point to
  (mix:MergedTrackTruth).
  """
  for label in associatorLabels:
    if not hasattr(process, label):
      continue
    module = getattr(process, label)
    for container in (list(process.tasks.values()) + list(process.sequences.values())
                      + list(process.paths.values()) + list(process.endpaths.values())):
      try:
        container.remove(module)
      except Exception:
        pass
    delattr(process, label)
    print(f"SmartPixels: removed in-process associator '{label}' (maps read from input file)")
  return process


# The three valid truth postures (see doc/PostureGapStudy.md and
# mem:smartpixels-pu-posture-note). Association maps are Ref/Ptr-keyed to specific
# collections; mixing postures yields a stale map + remade tracks = SILENT
# passthrough (study-breaking), so the vocabulary is closed and validated loudly.
TRUTHSOURCE_CHOICES = ("inJob", "fromFile", "fromFileStubs")


# FRAMING (user directive 2026-07-20): the ONLY story is 5-par OT-only vs 5-par
# OT+IT. promptHnpar=5 is the DEFAULT for all SmartPix productions henceforth so the
# PROMPT reference (OT-only L1TTrack) carries a REAL fitted d0 + 5x5 covariance, the
# right seed for b-tagging impact-parameter / vertexing resolution-vs-truth. 4-par
# (d0 pinned to 0) remains SELECTABLE for ablation but is off the default path.
#
# The prompt tracklet producer's shipped default is Extended=False, Hnpar=4. Setting
# Hnpar=5 (keeping Extended=False) makes the PROMPT collection 5-par. Static trace
# (2026-07-20): no guard couples Hnpar to Extended; the prompt USEHYBRID path runs
# TMTT KFParamsComb dimensioned purely by nHelixPar (KFbase.cc:81 pins d0 on
# nHelixPar==4, NOT on extended; L1fittedTrack.h fills the full 5x5 for nHelixPar==5).
# This is the PRIME-TARGET seed for the prompt digiRefit (seedCovMode=trackCov); at
# scale it VALIDATED (mem:smartpixels-prime-target-production: 100% nonzero prompt d0,
# spxSeedCovOK all-True across 4 variants x 3 PU files). The extended chain is always
# Hnpar=5 regardless.
PROMPT_HNPAR_CHOICES = (4, 5)


def attachFromFileStubsChain(process, extendedTracks=True, promptHnpar=5):
  """POSTURE C: rebuild NEW-layout L1 tracks from the FILE'S persisted stub tier.

  promptHnpar (4 | 5): helix-parameter count of the PROMPT tracklet fit. 5 (DEFAULT,
  the 5-par-only framing) makes the prompt collection carry a REAL fitted d0 + 5x5
  covariance, so the prompt digiRefit seeds from a genuine 5-par covariance
  (seedCovMode=trackCov) and the OT-only prompt reference has real d0 for b-tagging /
  resolution-vs-truth. 4 (ablation only) keeps d0 pinned to 0 with a 4x4 covariance --
  the prompt digiRefit then seeds with the weak-d0-prior fallback (the producer's
  nFitPars()==5 guard auto-selects it). The extended chain is always Hnpar=5 regardless.
  Only affects the SEED helix/covariance; the stub geometry (and hence the Part-B stub
  tables) is identical either way.

  The PU RelVal persists the entire stub tier WITH pileup (TTStubs + all three
  truth-association map sets + TrackingParticles + offlineBeamSpot, all HLT-process).
  The tracklet track finder consumes STUBS, not digis, so we re-run in-job:

    file TTStubsFromPhase2TrackerDigis:StubAccepted
      -> trackerDTC::ProducerDTC (TTDTC, transient)
      -> l1tTTTracksFromTrackletEmulation (readMoreMcTruth=False, file beamspot)
      -> TTTrackAssociatorFromPixelDigis re-run vs the FILE'S cluster/stub maps

  giving NEW TTTracks (post-PR#51503 layout, real helix covariance, PU tracks
  included) + a FRESH track-truth map keyed to the new collection. No DIGI, no
  mixing, no cluster/stub remaking -> so seedCovMode="trackCov" becomes VALID on
  PU files. Same-label in-process production shadows the file's HLT branches for
  every downstream consumer configured without a process name (identical to the
  posture-B labeling model), so the digiRefit producer defaults need no changes.

  Process-name handling of the file maps: the re-run associators consume
  TTClusterAssociatorFromPixelDigis:ClusterAccepted and
  TTStubAssociatorFromPixelDigis:StubAccepted with NO process name; because we do
  NOT schedule in-process producers with those labels, consumption resolves to the
  file's HLT products (the new tracks' stub Refs point into the file's stub
  product, the same ProductID the file maps are keyed by -> coherent by
  construction). Likewise offlineBeamSpot (no in-process producer scheduled) and
  the TrackingParticles (mix:MergedTrackTruth) resolve to the file.

  extendedTracks: also rebuild the extended (displaced) chain
  (l1tTTTracksFromExtendedTrackletEmulation + TTTrackAssociatorFromPixelDigisExtended,
  Hnpar=5) so the extended digiRefit variant's inputs are fresh new-layout tracks
  and trackCov is valid for them too. The prompt-only PoC left it off because its
  extended inputs would have resolved to old-layout file tracks (tripping the
  trackCov guard) -- with the chain rebuilt that no longer applies.

  NEVER schedules DIGI, the cluster/stub producers, the cluster/stub associators,
  or offlineBeamSpot: those need mix:Tracker PixelDigiSimLinks (transient DIGI
  products) and would recompute the beamspot signal-only. Only ProducerDTC, the
  tracklet emulator(s), and the track associator(s) are put on a Task.
  """
  # ES chain (mirrors the PoC's imports). DTC_cff pulls trackerDTC::ProducerSetup;
  # the tracklet cfi pulls trklet::ProducerSetup; ProducerHPH_cff pulls hph::Setup
  # (consumed by the digiRefit producer). TrackTrigger_cff provides the TTStub
  # algorithm ES records the setups reference.
  if promptHnpar not in PROMPT_HNPAR_CHOICES:
    raise ValueError(f"promptHnpar must be one of {PROMPT_HNPAR_CHOICES}, got {promptHnpar!r}")

  process.load('L1Trigger.TrackTrigger.TrackTrigger_cff')
  process.load('L1Trigger.TrackerDTC.DTC_cff')
  process.load('L1Trigger.TrackFindingTracklet.l1tTTTracksFromTrackletEmulation_cfi')
  process.load('L1Trigger.TrackFindingTracklet.ProducerHPH_cff')
  process.load('SimTracker.TrackTriggerAssociation.TTTrackAssociation_cfi')

  # PRIME TARGET: 5-par PROMPT seed. Keep Extended=False (the prompt seed geometry)
  # and only bump the helix-parameter count so the prompt fit floats a real d0 + 5x5
  # covariance. Extended stays Hnpar=5 unconditionally (its shipped clone default).
  process.l1tTTTracksFromTrackletEmulation.Hnpar = cms.uint32(promptHnpar)

  # Re-run the PROMPT track associator against the NEW tracks; cluster/stub maps
  # and TPs resolve to the file's HLT branches (no in-process producer for them).
  process.TTTrackAssociatorFromPixelDigis.TTTracks = cms.VInputTag(
      cms.InputTag("l1tTTTracksFromTrackletEmulation", "Level1TTTracks"))

  chainModules = ["ProducerDTC",
                  "l1tTTTracksFromTrackletEmulation",
                  "TTTrackAssociatorFromPixelDigis"]

  if extendedTracks:
    # Extended (displaced) chain: same stubs/beamspot, Hnpar=5, its own associator
    # (clone of the prompt one pointed at the extended tracks), mirroring
    # L1HybridEmulationTracks_cff.
    process.TTTrackAssociatorFromPixelDigisExtended = \
        process.TTTrackAssociatorFromPixelDigis.clone(
            TTTracks=cms.VInputTag(
                cms.InputTag("l1tTTTracksFromExtendedTrackletEmulation", "Level1TTTracks")))
    chainModules += ["l1tTTTracksFromExtendedTrackletEmulation",
                     "TTTrackAssociatorFromPixelDigisExtended"]

  # A Task so unscheduled mode auto-runs the chain; associate it to every path.
  process.l1tSmartPixelsFromFileStubsTask = cms.Task(
      *[getattr(process, m) for m in chainModules])
  for _, path in process.paths_().items():
    path.associate(process.l1tSmartPixelsFromFileStubsTask)
  for _, epath in process.endpaths_().items():
    epath.associate(process.l1tSmartPixelsFromFileStubsTask)
  print(f"SmartPixels posture-C: prompt tracklet Hnpar={promptHnpar} "
        + ("(5-par PRIME seed: real d0 + 5x5 cov -> prompt digiRefit trackCov)"
           if promptHnpar == 5
           else "(4-par: d0 pinned 0 -> prompt digiRefit weak-d0-prior fallback)"))
  return process, chainModules


def _applyTruthSource(process, truthSource, variants, extendedTracks=True, promptHnpar=5):
  if truthSource not in TRUTHSOURCE_CHOICES:
    raise ValueError(f"truthSource must be one of {TRUTHSOURCE_CHOICES}, got '{truthSource}'")
  # Truth is load-bearing for these modes (silent per-track passthrough otherwise):
  # make the requirement visible in the log either way.
  truthModes = [mode for mode, _ in variants if mode in TRUTH_REQUIRED_MODES]
  if truthModes:
    print(f"SmartPixels: modes {truthModes} REQUIRE the TP/truth association "
          f"(truthSource={truthSource}); without it tracks pass through unmodified.")
  removed = []
  attached = []
  if truthSource == "fromFile":
    process = useTruthAssociationFromFile(process)
    removed = list(TRUTH_ASSOCIATOR_LABELS)
  elif truthSource == "fromFileStubs":
    # Posture C: rebuild tracks + track-truth map(s) from the file's stubs;
    # remove ONLY the cluster/stub associators (they need mix:Tracker simlinks) --
    # the TRACK associator IS re-run (that is the difference from fromFile), and
    # the DIGI-tier + cluster/stub producers are never scheduled here either.
    process = useTruthAssociationFromFile(
        process, associatorLabels=("TTClusterAssociatorFromPixelDigis",
                                   "TTStubAssociatorFromPixelDigis"))
    process, attached = attachFromFileStubsChain(process, extendedTracks=extendedTracks,
                                                 promptHnpar=promptHnpar)
    removed = ["TTClusterAssociatorFromPixelDigis", "TTStubAssociatorFromPixelDigis"]
  # Load-bearing summary line (posture, chains attached, what was removed).
  if truthSource == "inJob":
    chainDesc = "none (in-process associators run unscheduled)"
  else:
    chainDesc = attached or "none"
  ext = extendedTracks if truthSource == "fromFileStubs" else "n/a"
  print(f"SmartPixels truthSource={truthSource}: chain attached={chainDesc}; "
        f"removed in-process modules={removed or 'none'}; extendedTracks={ext}")
  return process


def _scheduleVariantModules(process, modules, taskName):
  task = cms.Task(*[getattr(process, label) for label in modules])
  setattr(process, taskName, task)
  for _, path in process.paths_().items():
    path.associate(task)
  return process


def injectSmartPixelsTrackProducer(process,
                                   smartPixelsEmulatorMode="passthrough",
                                   addAssociation=True,
                                   l1tSmartPixelsTrackProducerLabel="l1tSmartPixelsTrackProducer",
                                   l1tSmartPixelsTrackProducerExtendedLabel="l1tSmartPixelsTrackProducerExtended",
                                   skipModuleTypes=None,
                                   printProcessInfo=False):
  """WF2 core: rewrite every producer InputTag that reads the standard tracklet
  tracks (l1tTTTracksFromTrackletEmulation[Extended] / Level1TTTracks) to the
  given SmartPixels collections, co-opting the full downstream chain (GTT input,
  vertex finder, correlator layers, track jets/HT/MET, GT, muon matching, ...).

  skipModuleTypes protects (by C++ type) the producers that must keep reading
  the original tracks: the SmartPixels producers themselves (their input), the
  TTTrack associators, the tracklet track finder, and the L1Nano track-table
  producer (so the L1TTrack nano table remains the untouched reference).
  """
  if skipModuleTypes is None:
    skipModuleTypes = ["L1SmartPixelsTrackProducer",
                       "TTTrackAssociator_Phase2TrackerDigi_",
                       "L1FPGATrackProducer",
                       "SimpleL1TTTrackCandidateFlatTableProducer",
                       # keep the truth table aligned with the (skipped) reference track
                       # table; Ref-based truth lookups only resolve on the original
                       # collection the associator ran on
                       "L1TrackTruthTableProducer"]

  def print_helper(item):
    if hasattr(item, "items"):
      for key, value in item.items():
        print("\t>>", key, "== ", value)
    else:
      print("\t==", item)

  if printProcessInfo:
    print("====Process Info====")
    for section in ["source_", "schedule_", "sequences_", "paths_", "endpaths_",
                    "producers_", "filters_", "analyzers_"]:
      print(section, "\n")
      print_helper(getattr(process, section)())

  for modname, module in process.producers_().items():
    if hasattr(module, '_TypedParameterizable__type') and module._TypedParameterizable__type in skipModuleTypes:
      print("Skipping module type ", module._TypedParameterizable__type, " with module name ", modname)
      continue
    for param_name in module.parameterNames_():
      param = getattr(module, param_name)
      if isinstance(param, cms.InputTag):
        if param.getModuleLabel() == "l1tTTTracksFromTrackletEmulation":
          print("module._TypedParameterizable__type --> ", module._TypedParameterizable__type)
          print("l1tTTTracksFromTrackletEmulation --> ", modname, param_name, param)
          setattr(module, param_name, cms.InputTag(l1tSmartPixelsTrackProducerLabel, "Level1TTTracks"))
          print("                                 --> ", modname, param_name, getattr(module, param_name))
        if param.getModuleLabel() == "l1tTTTracksFromExtendedTrackletEmulation":
          print("l1tTTTracksFromExtendedTrackletEmulation --> ", modname, param_name, param)
          setattr(module, param_name, cms.InputTag(l1tSmartPixelsTrackProducerExtendedLabel, "Level1TTTracks"))
          print("                                         --> ", modname, param_name, getattr(module, param_name))

  if addAssociation:
    process.load('L1Trigger.Phase3SmartPixels.l1tSmartPixelsTrackProducer_cfi')
    process.l1tSmartPixelsTrackProducer.smartPixelsEmulatorMode = cms.string(smartPixelsEmulatorMode)
    process.l1tSmartPixelsTrackProducerExtended.smartPixelsEmulatorMode = cms.string(smartPixelsEmulatorMode)
    process.l1tSmartPixelsTrackProducerTask = cms.Task(process.l1tSmartPixelsTrackProducer, process.l1tSmartPixelsTrackProducerExtended)

    print("Associating L1SmartPixelsTrackProducer into the process")
    for pathname, cmspath in process.paths_().items():
      cmspath.associate(process.l1tSmartPixelsTrackProducerTask)

  return process


# ---------------------------------------------------------------------------
# WF1: coexist — standard tracks AND SmartPixels variant tables in one L1Nano
# ---------------------------------------------------------------------------
def smartPixelsCoexist(process, variants=None, correctionSet=DEFAULT_CORRECTION_SET, addNanoTables=True,
                       truthSource="inJob", digiRefitConfig=None, extendedTracks=True, promptHnpar=5):
  """Add SmartPixels track collections (default: one passthrough variant)
  alongside the standard tracks, plus one pair of L1Nano track tables per
  variant. Nothing downstream is rewired: all other L1 objects still reflect
  the standard tracks, enabling in-file track-to-track comparisons.

  truthSource (three postures — see doc/PostureGapStudy.md):
   - 'inJob' (default): run the TT truth associators in-job — valid when the job
     also runs DIGI or the input retained mix:Tracker simlinks. Fresh new-layout
     tracks with real covariance; PU is destroyed by re-digitization of
     signal-only g4SimHits, so this is the no-PU / re-digitized posture.
   - 'fromFile': remove ALL in-process associators so STEP1-style association
     maps are read from the input file. Real PU, but the file's tracks are old
     layout -> helixCovMat all-zero -> digiRefit must use
     seedCovMode='parametrized' (the trackCov guard throws loudly otherwise).
   - 'fromFileStubs' (POSTURE C, PU + real covariance): rebuild NEW-layout tracks
     from the file's persisted stub tier (ProducerDTC -> tracklet emulator(s) ->
     re-run TRACK associator vs the file's cluster/stub maps); remove only the
     cluster/stub associators, never DIGI. Real PU AND real per-track covariance,
     so seedCovMode='trackCov' (the digiRefit default) is VALID here.
  Truth is REQUIRED for the regression/TP/digiRefit modes (silent per-track
  passthrough otherwise).

  extendedTracks (fromFileStubs only): also rebuild the extended (displaced)
  chain so the extended digiRefit variant's inputs are fresh new-layout tracks
  (trackCov valid for them too). Default True.

  promptHnpar (fromFileStubs only; 4 | 5): helix-parameter count of the PROMPT
  tracklet fit. DEFAULT 5 (the 5-par-only framing: 5-par OT-only vs 5-par OT+IT is
  the whole story): the prompt collection carries a REAL fitted d0 + 5x5 covariance,
  so the OT-only prompt reference has real d0 for b-tagging / resolution-vs-truth and
  the prompt digiRefit variant seeds from a genuine 5-par covariance
  (seedCovMode=trackCov, the digiRefit default). promptHnpar=4 (ablation only) keeps
  d0 pinned to 0; the prompt digiRefit then uses the weak-d0-prior fallback (the
  producer's nFitPars()==5 guard auto-selects it per-collection). The extended chain
  is always Hnpar=5. Only fromFileStubs honors this knob (inJob/fromFile use the
  shipped prompt cfi default of 4). The 5-par prompt path is VALIDATED at scale
  (mem:smartpixels-prime-target-production: 100% nonzero prompt d0, spxSeedCovOK
  all-True across 4 variants x 3 PU files).

  digiRefitConfig: dict merged over DIGIREFIT_DEFAULTS (validated loudly),
  relevant only for a digiRefit variant. NOTE: combining truthSource='fromFile'
  (old posture A) with seedCovMode='trackCov' is not special-cased at config
  time -- it fails loudly at runtime via the SmartPixelsSeedCovMissing guard, by
  design; use 'fromFileStubs' for trackCov on PU files.

  cmsDriver (defaults):  --customise L1Trigger/Phase3SmartPixels/customizeSmartPixels_cff.smartPixelsCoexist
  cmsDriver (explicit):  --customise_commands 'from L1Trigger.Phase3SmartPixels.customizeSmartPixels_cff import smartPixelsCoexist; process = smartPixelsCoexist(process, variants=[("correctionlibRegression", "1100")])'
  """
  if variants is None:
    variants = [("passthrough", None)]
  variants = _normalizeVariants(variants)

  process, modules = addSmartPixelsTrackProducerVariants(process, variants, correctionSet,
                                                         digiRefitConfig=digiRefitConfig)
  process = _scheduleVariantModules(process, modules, "l1tSmartPixelsCoexistTask")
  process = _applyTruthSource(process, truthSource, variants, extendedTracks=extendedTracks,
                              promptHnpar=promptHnpar)

  if addNanoTables:
    from DPGAnalysis.Phase3SmartPixelsNanoAOD.l1tPh3SmartPixelsNano_cff import (
        addPh3L1SmartPixelsTracks, addPh3L1SmartPixelsRefitTables,
        addPh3L1SmartPixelsStubTables)
    for mode, activeSP in variants:
      prompt, extended = smartPixelsVariantLabels(mode, activeSP)
      suffix = smartPixelsVariantSuffix(mode, activeSP)
      addPh3L1SmartPixelsTracks(process, srcLabel=prompt, srcExtendedLabel=extended, tableSuffix=suffix)
      # digiRefit variants ALSO emit a SmartPixelsRefitSidecar -> add the per-hit
      # link table + track extension table. Other modes produce no sidecar.
      if mode == "digiRefit":
        addPh3L1SmartPixelsRefitTables(process, srcLabel=prompt, srcExtendedLabel=extended, tableSuffix=suffix)
    # Reference-track OT stub-position tables (task Part B): added ONCE per job,
    # keyed to the reference L1TTrack / L1TExtTrack collections (not the refit
    # variants -- the refit reuses the same stub refs). Idempotent, so calling it
    # per coexist is safe. The prompt reference is 5-par when promptHnpar=5 (this
    # only affects the seed helix; the stub geometry is unchanged either way).
    addPh3L1SmartPixelsStubTables(process, doExtended=extendedTracks)
  return process


# ---------------------------------------------------------------------------
# WF2: co-opt — ONE variant feeds the entire downstream chain ("pure view")
# ---------------------------------------------------------------------------
def smartPixelsCoopt(process, mode="passthrough", activeSP=None,
                     correctionSet=DEFAULT_CORRECTION_SET,
                     addPh3Table=False, skipModuleTypes=None,
                     truthSource="inJob", digiRefitConfig=None, extendedTracks=True, promptHnpar=5):
  """Produce ONE SmartPixels variant in-job and inject it into every downstream
  consumer of the standard tracklet tracks. Any nano flavor run in this job then
  reflects that single track interpretation; comparisons are file-to-file
  against a baseline job without this customisation (one job per variant).

  addPh3Table=True additionally writes the injected collections as explicit
  Ph3 SmartPixels track tables (cross-check); default is the pure view, where
  the standard L1TTrack table (protected by the skip-list) remains the only
  extra reference.

  truthSource: see smartPixelsCoexist — 'inJob' (default), 'fromFile', or
  'fromFileStubs' (posture C: rebuild new-layout tracks from file stubs, real PU +
  real covariance, seedCovMode='trackCov' valid). Truth is REQUIRED for the
  regression/TP/digiRefit modes.

  extendedTracks (fromFileStubs only): also rebuild the extended chain. Default True.

  digiRefitConfig: dict merged over DIGIREFIT_DEFAULTS (validated loudly),
  relevant only when mode='digiRefit'. mode='refit' raises the RESERVED
  NotImplementedError. Combining truthSource='fromFile' with
  seedCovMode='trackCov' fails loudly at runtime (SmartPixelsSeedCovMissing), by
  design; use 'fromFileStubs' for trackCov on PU files.

  cmsDriver (defaults):  --customise L1Trigger/Phase3SmartPixels/customizeSmartPixels_cff.smartPixelsCoopt
  cmsDriver (explicit):  --customise_commands 'from L1Trigger.Phase3SmartPixels.customizeSmartPixels_cff import smartPixelsCoopt; process = smartPixelsCoopt(process, mode="correctionlibRegression", activeSP="1100")'
  """
  process, modules = addSmartPixelsTrackProducerVariants(process, [(mode, activeSP)], correctionSet,
                                                         digiRefitConfig=digiRefitConfig)
  process = _scheduleVariantModules(process, modules, "l1tSmartPixelsCooptTask")
  process = _applyTruthSource(process, truthSource, [(mode, activeSP)], extendedTracks=extendedTracks,
                              promptHnpar=promptHnpar)

  prompt, extended = smartPixelsVariantLabels(mode, activeSP)
  process = injectSmartPixelsTrackProducer(process,
                                           addAssociation=False,
                                           l1tSmartPixelsTrackProducerLabel=prompt,
                                           l1tSmartPixelsTrackProducerExtendedLabel=extended,
                                           skipModuleTypes=skipModuleTypes)

  if addPh3Table:
    from DPGAnalysis.Phase3SmartPixelsNanoAOD.l1tPh3SmartPixelsNano_cff import (
        addPh3L1SmartPixelsTracks, addPh3L1SmartPixelsRefitTables,
        addPh3L1SmartPixelsStubTables)
    suffix = smartPixelsVariantSuffix(mode, activeSP)
    addPh3L1SmartPixelsTracks(process, srcLabel=prompt, srcExtendedLabel=extended, tableSuffix=suffix)
    if mode == "digiRefit":
      addPh3L1SmartPixelsRefitTables(process, srcLabel=prompt, srcExtendedLabel=extended, tableSuffix=suffix)
    # Reference-track OT stub-position tables (task Part B), keyed to the reference
    # L1TTrack / L1TExtTrack collections. In coopt the injector rewires downstream
    # consumers to the variant, but the reference tracklet collections still exist
    # (they are the injector's protected input) so the stub tables resolve.
    addPh3L1SmartPixelsStubTables(process, doExtended=extendedTracks)
  return process


# ---------------------------------------------------------------------------
# Legacy no-argument wrappers (cmsDriver --customise calls func(process) only)
# ---------------------------------------------------------------------------
def injectSmartPixelsTrackProducer_passthrough(process):
  return injectSmartPixelsTrackProducer(process, smartPixelsEmulatorMode="passthrough")

def injectSmartPixelsTrackProducer_passthroughFloat(process):
  return injectSmartPixelsTrackProducer(process, smartPixelsEmulatorMode="passthroughFloat")

def injectSmartPixelsTrackProducer_passthroughHW(process):
  return injectSmartPixelsTrackProducer(process, smartPixelsEmulatorMode="passthroughHW")

def injectSmartPixelsTrackProducer_trackingParticleTruth(process):
  return injectSmartPixelsTrackProducer(process, smartPixelsEmulatorMode="trackingParticleTruth")

def injectSmartPixelsTrackProducer_correctionlibRegression(process):
  return injectSmartPixelsTrackProducer(process, smartPixelsEmulatorMode="correctionlibRegression")

def referenceSmartPixelsTrackProducer(process):
  return injectSmartPixelsTrackProducer(process, addAssociation=False)
