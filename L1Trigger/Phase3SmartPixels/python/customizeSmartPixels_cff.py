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
# Modes parametrized by the active pixel-layer configuration (need activeSP)
ACTIVE_LAYER_MODES = ["correctionlibRegression", "correctionlibTPToySmear"]

REGRESSION_CONFIGS = ["0000",
                      "0001", "0010", "0100", "1000",
                      "0011", "0101", "0110", "1001", "1010", "1100",
                      "0111", "1011", "1101", "1110",
                      "1111",
                      ]

DEFAULT_CORRECTION_SET = "L1Trigger/Phase3SmartPixels/data/spixel_smear_all_configs_barrel_CalV1_v2p1_compound.json"


def smartPixelsVariantSuffix(mode, activeSP=None):
  """Unique label/table suffix for a variant, e.g. 'passthrough' or 'correctionlibRegressionAAII'."""
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


def addSmartPixelsTrackProducerVariants(process, variants=None, correctionSet=DEFAULT_CORRECTION_SET):
  """STEP1: instantiate the SmartPixels track producer (prompt+extended) for each variant.

  variants=None reproduces the legacy behavior: all passthrough-family modes plus
  all 16 regression configs, produced concurrently (they are cheap smears/copies
  of the same base tracklet tracks).
  Returns (process, moduleNames); scheduling the modules is up to the caller
  (see smartPixelsCoexist/smartPixelsCoopt for path association).
  """
  if variants is None:
    variants = _allVariants()
  variants = _normalizeVariants(variants)

  process.load('L1Trigger.Phase3SmartPixels.l1tSmartPixelsTrackProducer_cfi')
  modules = []
  for mode, activeSP in variants:
    prompt, extended = smartPixelsVariantLabels(mode, activeSP)
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

  return process, modules


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
def smartPixelsCoexist(process, variants=None, correctionSet=DEFAULT_CORRECTION_SET, addNanoTables=True):
  """Add SmartPixels track collections (default: one passthrough variant)
  alongside the standard tracks, plus one pair of L1Nano track tables per
  variant. Nothing downstream is rewired: all other L1 objects still reflect
  the standard tracks, enabling in-file track-to-track comparisons.

  cmsDriver (defaults):  --customise L1Trigger/Phase3SmartPixels/customizeSmartPixels_cff.smartPixelsCoexist
  cmsDriver (explicit):  --customise_commands 'from L1Trigger.Phase3SmartPixels.customizeSmartPixels_cff import smartPixelsCoexist; process = smartPixelsCoexist(process, variants=[("correctionlibRegression", "1100")])'
  """
  if variants is None:
    variants = [("passthrough", None)]
  variants = _normalizeVariants(variants)

  process, modules = addSmartPixelsTrackProducerVariants(process, variants, correctionSet)
  process = _scheduleVariantModules(process, modules, "l1tSmartPixelsCoexistTask")

  if addNanoTables:
    from DPGAnalysis.Phase3SmartPixelsNanoAOD.l1tPh3SmartPixelsNano_cff import addPh3L1SmartPixelsTracks
    for mode, activeSP in variants:
      prompt, extended = smartPixelsVariantLabels(mode, activeSP)
      addPh3L1SmartPixelsTracks(process,
                                srcLabel=prompt,
                                srcExtendedLabel=extended,
                                tableSuffix=smartPixelsVariantSuffix(mode, activeSP))
  return process


# ---------------------------------------------------------------------------
# WF2: co-opt — ONE variant feeds the entire downstream chain ("pure view")
# ---------------------------------------------------------------------------
def smartPixelsCoopt(process, mode="passthrough", activeSP=None,
                     correctionSet=DEFAULT_CORRECTION_SET,
                     addPh3Table=False, skipModuleTypes=None):
  """Produce ONE SmartPixels variant in-job and inject it into every downstream
  consumer of the standard tracklet tracks. Any nano flavor run in this job then
  reflects that single track interpretation; comparisons are file-to-file
  against a baseline job without this customisation (one job per variant).

  addPh3Table=True additionally writes the injected collections as explicit
  Ph3 SmartPixels track tables (cross-check); default is the pure view, where
  the standard L1TTrack table (protected by the skip-list) remains the only
  extra reference.

  cmsDriver (defaults):  --customise L1Trigger/Phase3SmartPixels/customizeSmartPixels_cff.smartPixelsCoopt
  cmsDriver (explicit):  --customise_commands 'from L1Trigger.Phase3SmartPixels.customizeSmartPixels_cff import smartPixelsCoopt; process = smartPixelsCoopt(process, mode="correctionlibRegression", activeSP="1100")'
  """
  process, modules = addSmartPixelsTrackProducerVariants(process, [(mode, activeSP)], correctionSet)
  process = _scheduleVariantModules(process, modules, "l1tSmartPixelsCooptTask")

  prompt, extended = smartPixelsVariantLabels(mode, activeSP)
  process = injectSmartPixelsTrackProducer(process,
                                           addAssociation=False,
                                           l1tSmartPixelsTrackProducerLabel=prompt,
                                           l1tSmartPixelsTrackProducerExtendedLabel=extended,
                                           skipModuleTypes=skipModuleTypes)

  if addPh3Table:
    from DPGAnalysis.Phase3SmartPixelsNanoAOD.l1tPh3SmartPixelsNano_cff import addPh3L1SmartPixelsTracks
    addPh3L1SmartPixelsTracks(process,
                              srcLabel=prompt,
                              srcExtendedLabel=extended,
                              tableSuffix=smartPixelsVariantSuffix(mode, activeSP))
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
