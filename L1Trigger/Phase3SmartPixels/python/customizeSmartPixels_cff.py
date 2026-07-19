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
    # --- Kalman seed (user decision 2026-07-18: config-switchable, trackCov default) ---
    "seedNPar": 5,            # 4 | 5: seed-track parametrization entering the KF
    "seedCovMode": "trackCov",  # "trackCov" (TTTrack helixCovMat) | "parametrized" (ablation/fallback)
    "paramSigmas": (1e-4, 1e-3, 2e-3, 0.06, 0.05),  # parametrized-mode sigmas (rInv[cm^-1],phi0,tanL,z0[cm],d0[cm])
    # --- correctionlib payload paths (empty defaults acceptable for Phase 0) ---
    "smarthitTrueSet": "",    # Stack A "smarthit_true" payload (eff/residual/angle response)
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
  module.digiRefitSeedNPar = cms.int32(resolved["seedNPar"])
  module.digiRefitSeedCovMode = cms.string(resolved["seedCovMode"])
  module.digiRefitParamSigmas = cms.vdouble(*resolved["paramSigmas"])
  module.digiRefitPixelavAngleSet = cms.string(resolved["pixelavAngleSet"])
  module.digiRefitSmarthitFakeSet = cms.string(resolved["smarthitFakeSet"])


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
    print(f"SmartPixels truthSource=fromFile: removed in-process associator '{label}' (maps read from input file)")
  return process


def _applyTruthSource(process, truthSource, variants):
  if truthSource not in ("inJob", "fromFile"):
    raise ValueError(f"truthSource must be 'inJob' or 'fromFile', got '{truthSource}'")
  # Truth is load-bearing for these modes (silent per-track passthrough otherwise):
  # make the requirement visible in the log either way.
  truthModes = [mode for mode, _ in variants if mode in TRUTH_REQUIRED_MODES]
  if truthModes:
    print(f"SmartPixels: modes {truthModes} REQUIRE the TP/truth association "
          f"(truthSource={truthSource}); without it tracks pass through unmodified.")
  if truthSource == "fromFile":
    process = useTruthAssociationFromFile(process)
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
                       truthSource="inJob", digiRefitConfig=None):
  """Add SmartPixels track collections (default: one passthrough variant)
  alongside the standard tracks, plus one pair of L1Nano track tables per
  variant. Nothing downstream is rewired: all other L1 objects still reflect
  the standard tracks, enabling in-file track-to-track comparisons.

  truthSource: 'inJob' (default) runs the TT truth associators in-job — valid
  when the job also runs DIGI or the input retained mix:Tracker simlinks;
  'fromFile' removes them so STEP1-style association maps are read from the
  input file. Truth is REQUIRED for the regression/TP modes (silent per-track
  passthrough otherwise).

  digiRefitConfig: dict merged over DIGIREFIT_DEFAULTS (validated loudly),
  relevant only for a digiRefit variant. Phase 0: a digiRefit variant raises
  NotImplementedError at producer instantiation (surface reserved, no run).

  cmsDriver (defaults):  --customise L1Trigger/Phase3SmartPixels/customizeSmartPixels_cff.smartPixelsCoexist
  cmsDriver (explicit):  --customise_commands 'from L1Trigger.Phase3SmartPixels.customizeSmartPixels_cff import smartPixelsCoexist; process = smartPixelsCoexist(process, variants=[("correctionlibRegression", "1100")])'
  """
  if variants is None:
    variants = [("passthrough", None)]
  variants = _normalizeVariants(variants)

  process, modules = addSmartPixelsTrackProducerVariants(process, variants, correctionSet,
                                                         digiRefitConfig=digiRefitConfig)
  process = _scheduleVariantModules(process, modules, "l1tSmartPixelsCoexistTask")
  process = _applyTruthSource(process, truthSource, variants)

  if addNanoTables:
    from DPGAnalysis.Phase3SmartPixelsNanoAOD.l1tPh3SmartPixelsNano_cff import (
        addPh3L1SmartPixelsTracks, addPh3L1SmartPixelsRefitTables)
    for mode, activeSP in variants:
      prompt, extended = smartPixelsVariantLabels(mode, activeSP)
      suffix = smartPixelsVariantSuffix(mode, activeSP)
      addPh3L1SmartPixelsTracks(process, srcLabel=prompt, srcExtendedLabel=extended, tableSuffix=suffix)
      # digiRefit variants ALSO emit a SmartPixelsRefitSidecar -> add the per-hit
      # link table + track extension table. Other modes produce no sidecar.
      if mode == "digiRefit":
        addPh3L1SmartPixelsRefitTables(process, srcLabel=prompt, srcExtendedLabel=extended, tableSuffix=suffix)
  return process


# ---------------------------------------------------------------------------
# WF2: co-opt — ONE variant feeds the entire downstream chain ("pure view")
# ---------------------------------------------------------------------------
def smartPixelsCoopt(process, mode="passthrough", activeSP=None,
                     correctionSet=DEFAULT_CORRECTION_SET,
                     addPh3Table=False, skipModuleTypes=None,
                     truthSource="inJob", digiRefitConfig=None):
  """Produce ONE SmartPixels variant in-job and inject it into every downstream
  consumer of the standard tracklet tracks. Any nano flavor run in this job then
  reflects that single track interpretation; comparisons are file-to-file
  against a baseline job without this customisation (one job per variant).

  addPh3Table=True additionally writes the injected collections as explicit
  Ph3 SmartPixels track tables (cross-check); default is the pure view, where
  the standard L1TTrack table (protected by the skip-list) remains the only
  extra reference.

  truthSource: see smartPixelsCoexist — 'inJob' (default) or 'fromFile'.
  Truth is REQUIRED for the regression/TP modes.

  digiRefitConfig: dict merged over DIGIREFIT_DEFAULTS (validated loudly),
  relevant only when mode='digiRefit'. Phase 0: mode='digiRefit' raises
  NotImplementedError at producer instantiation (surface reserved, no run);
  mode='refit' raises the RESERVED NotImplementedError.

  cmsDriver (defaults):  --customise L1Trigger/Phase3SmartPixels/customizeSmartPixels_cff.smartPixelsCoopt
  cmsDriver (explicit):  --customise_commands 'from L1Trigger.Phase3SmartPixels.customizeSmartPixels_cff import smartPixelsCoopt; process = smartPixelsCoopt(process, mode="correctionlibRegression", activeSP="1100")'
  """
  process, modules = addSmartPixelsTrackProducerVariants(process, [(mode, activeSP)], correctionSet,
                                                         digiRefitConfig=digiRefitConfig)
  process = _scheduleVariantModules(process, modules, "l1tSmartPixelsCooptTask")
  process = _applyTruthSource(process, truthSource, [(mode, activeSP)])

  prompt, extended = smartPixelsVariantLabels(mode, activeSP)
  process = injectSmartPixelsTrackProducer(process,
                                           addAssociation=False,
                                           l1tSmartPixelsTrackProducerLabel=prompt,
                                           l1tSmartPixelsTrackProducerExtendedLabel=extended,
                                           skipModuleTypes=skipModuleTypes)

  if addPh3Table:
    from DPGAnalysis.Phase3SmartPixelsNanoAOD.l1tPh3SmartPixelsNano_cff import (
        addPh3L1SmartPixelsTracks, addPh3L1SmartPixelsRefitTables)
    suffix = smartPixelsVariantSuffix(mode, activeSP)
    addPh3L1SmartPixelsTracks(process, srcLabel=prompt, srcExtendedLabel=extended, tableSuffix=suffix)
    if mode == "digiRefit":
      addPh3L1SmartPixelsRefitTables(process, srcLabel=prompt, srcExtendedLabel=extended, tableSuffix=suffix)
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
