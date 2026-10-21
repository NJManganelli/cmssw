import json

import FWCore.ParameterSet.Config as cms

# SmartPixels L1 nano add-on flavours (PhysicsTools/NanoAOD/autoNANO.py), composed on the
# cmsDriver line after the base tier, in any order:
#
#   -s NANO:@Phase2L1DPGwithGen+@L1Trk+@SPTrk+@SPSeedlessRefit+@SPActiveLayersAAAA+...
#
# Every add-on only ADDS tables. Index columns that link two tables ("links") exist only
# when both tables are in the nano, so each add-on ends with _reconcile, which derives all
# links, the wiring of the selected SmartPixels variant and the TrackingParticle selection
# from what the job contains at that point. Whichever end is added second therefore
# creates a link, and the order of the add-ons does not matter.
#
# The SmartPixels variant is chosen by one mode tag (@SPSeedlessRefit = seedlessRefit) and one
# mask tag (@SPActiveLayers<mask>, character k = ILk and IDk, A active, I inactive); the
# second of the two wires it through the L1 emulation (smartPixelsCoopt). IIII is the
# passthrough control: SmartPixels tables and wiring, tracks equal to their OT tracks, and
# no mode tables (a passthrough variant has no refit sidecar). Settings that are not in a
# tag name go through smartPixelsNanoParams (see there).
#
# The composition is recorded in process.smartPixelsNano. Completeness (a mode tag without
# a mask tag, ...) can only be judged once every customisation has run, so it is checked
# when cmsRun builds the process description.

_STATE = "smartPixelsNano"

# seedlessRefitConfig overrides the add-on flavours apply over SEEDLESSREFIT_DEFAULTS.
_TIER_SEEDLESSREFIT_DEFAULTS = {"angleSource": "emu"}

_MODE_TAGS = {"SPSeedlessRefit": "seedlessRefit", "SPSeedlessRefitV3": "seedlessRefitV3"}
_TP_TRUTH_TAGS = {"TPTruth1GeV": "charged1GeV",
                  "TPTruthLinked": "linked",
                  "TPTruth1GeVLinked": "charged1GeVOrLinked",
                  "TPTruthAll": "all"}
MASKS = [a + b + c + d for a in "AI" for b in "AI" for c in "AI" for d in "AI"]
PASSTHROUGH_MASK = "IIII"

# OT (tracklet) and SmartPixels track tables: side -> (table module label, table name)
_OT_TABLES = {"prompt": ("l1tTracksTable", "L1TTrack"), "ext": ("l1tExtTracksTable", "L1TExtTrack")}
_SP_TABLES = {"prompt": ("l1tPh3SmartPixelsTracksTable", "L1TSmartPixelsTrack"),
              "ext": ("l1tPh3ExtSmartPixelsTracksTable", "L1TSmartPixelsExtTrack")}
_OT_TRACKS = {"prompt": "l1tTTTracksFromTrackletEmulation", "ext": "l1tTTTracksFromExtendedTrackletEmulation"}
_OT_ASSOCIATORS = {"prompt": "TTTrackAssociatorFromPixelDigis", "ext": "TTTrackAssociatorFromPixelDigisExtended"}
_CAND_SIDE = {"L1PuppiCand": "prompt", "L1ExtPuppiCand": "ext"}
_OT_STUB_TABLE = "l1tPh3SmartPixelsAllStubTable"
_OT_DTC_TABLE = "l1tPh3SmartPixelsDTCStubTable"
_CLUSTER_TABLE = "l1tPh3SmartPixelsClusterTable"
_TP_TABLE = "l1tTrackingParticleTable"


class _SmartPixelsNanoState(cms.PSet):
    """The composition record; validated when the process description is filled."""

    def pythonTypeName(self):
        return "cms.PSet"

    def insertInto(self, parameterSet, myname):
        _validateComposition(self)
        super().insertInto(parameterSet, myname)


def _state(process):
    if not hasattr(process, _STATE):
        setattr(process, _STATE, _SmartPixelsNanoState(
            addOns=cms.vstring(),
            modeTag=cms.string(""),
            mask=cms.string(""),
            tpTruth=cms.string(""),
            variant=cms.string(""),
            coexistVariants=cms.vstring(),
            trackInputMode=cms.string("reemulateL1TrackFinding"),
            trackInputModeApplied=cms.string(""),
            extendedTracks=cms.bool(True),
            promptHnpar=cms.uint32(5),
            seedlessRefitConfig=cms.string(json.dumps(_TIER_SEEDLESSREFIT_DEFAULTS, sort_keys=True)),
            tier=cms.string(""),
        ))
    return getattr(process, _STATE)


def _validateComposition(st):
    addOns = list(st.addOns)
    modeTag, mask = st.modeTag.value(), st.mask.value()
    if modeTag and not mask:
        raise RuntimeError(f"SmartPixels nano: the mode tag @{modeTag} needs a mask tag "
                           f"(@SPActiveLayers<mask>, e.g. @SPActiveLayersAAAA)")
    if mask and mask != PASSTHROUGH_MASK and not modeTag:
        raise RuntimeError(f"SmartPixels nano: the mask tag @SPActiveLayers{mask} needs a mode tag "
                           f"(e.g. @SPSeedlessRefit); only @SPActiveLayers{PASSTHROUGH_MASK} stands alone")
    if "SPTrk" in addOns and not st.variant.value():
        raise RuntimeError("SmartPixels nano: @SPTrk needs a SmartPixels variant: compose a mode tag "
                           f"and a mask tag, or @SPActiveLayers{PASSTHROUGH_MASK} for the passthrough control")


def _scheduled(process):
    return set(process.l1tPh2NanoTask.moduleNames())


def addCooptNanoTables(process, variant):
    """smartPixelsCoopt(addNanoTables=True) on a stock flavour (e.g. @L1PFTrkNanowithGen): the
    canonical SmartPixels track tables of the variant coopt wired, L1TOTStub, and every link,
    as the add-on flavours derive them."""
    st = _begin(process, "coopt")
    if st.variant.value():
        raise RuntimeError(f"addCooptNanoTables: variant {st.variant.value()} is already wired")
    st.variant = _variantString(variant)
    _addSPTrackTables(process)
    _addOTStubTable(process)
    return _reconcile(process)


def addCoexistNanoTables(process, variants, trackInputMode, seedlessRefitConfig):
    """smartPixelsCoexist(addNanoTables=True), after it made the suffixed variant tables: records
    the variants, their posture and refit settings, so that every reconcile, now or by a later
    add-on, gives each variant the links the add-on flavours give the canonical tables
    (<track table>Stub, <track table>Cluster)."""
    st = _state(process)
    applied = st.trackInputModeApplied.value()
    if applied and applied != trackInputMode:
        raise RuntimeError(f"smartPixelsCoexist: trackInputMode {trackInputMode!r}, but this job already "
                           f"applied {applied!r}")
    st.trackInputMode = trackInputMode
    st.trackInputModeApplied = trackInputMode
    st.seedlessRefitConfig = json.dumps(seedlessRefitConfig or {}, sort_keys=True)
    for variant in variants:
        if _variantString(variant) not in st.coexistVariants:
            st.coexistVariants.append(_variantString(variant))
    return _reconcile(process)


def _begin(process, tag):
    if not hasattr(process, "l1tPh2NanoTask") or not (
            set(process.p2L1TablesTask.moduleNames()) & _scheduled(process)):
        raise RuntimeError(f"SmartPixels nano: @{tag} is an add-on; compose it after the base tier, "
                           "e.g. NANO:@Phase2L1DPGwithGen+@" + tag)
    st = _state(process)
    st.addOns.append(tag)
    return st


def _withGen(process):
    return hasattr(process, "genNanoTask") and process.genNanoTask in process.l1tPh2NanoTask._collection


def _seedlessRefitConfig(st):
    from L1Trigger.Phase3SmartPixels.customizeSmartPixels_cff import _resolveSeedlessRefitConfig
    return _resolveSeedlessRefitConfig(json.loads(st.seedlessRefitConfig.value()))


def _variant(st):
    """(mode, activeSP) the tags select, or None while incomplete."""
    mask = st.mask.value()
    if not mask:
        return None
    if mask == PASSTHROUGH_MASK:
        return ("passthrough", None)
    if not st.modeTag.value():
        return None
    return (_MODE_TAGS[st.modeTag.value()], mask.replace("A", "1").replace("I", "0"))


def _variantString(variant):
    return f"{variant[0]}:{variant[1] or ''}"


def _wiredVariant(st):
    if not st.variant.value():
        return None
    mode, _, activeSP = st.variant.value().partition(":")
    return (mode, activeSP or None)


def _addToNano(process, taskName, labels):
    if not hasattr(process, taskName):
        setattr(process, taskName, cms.Task())
        process.l1tPh2NanoTask.add(getattr(process, taskName))
    task = getattr(process, taskName)
    for label in labels:
        task.add(getattr(process, label))


# ---------------------------------------------------------------------------
# Add-on entry points (autoNANO customisations)
# ---------------------------------------------------------------------------
def addL1TrkTables(process):
    """@L1Trk: L1TTrack / L1TExtTrack (with the helix covariance) and their GTT columns;
    with the gen tier also the track truth columns."""
    _begin(process, "L1Trk")
    from DPGAnalysis.Phase2L1TNanoAOD.l1tPh2Nano_cff import addPh2L1Tracks, addPh2L1TrackTruth
    addPh2L1Tracks(process)
    if _withGen(process):
        addPh2L1TrackTruth(process)
    return _reconcile(process)


def addL1PFTables(process):
    """@L1PF: L1 PF / Puppi candidates, HGCal clusters and the SC4/SC8/SC4NG constituent links."""
    _begin(process, "L1PF")
    from DPGAnalysis.Phase2L1TNanoAOD.l1tPh2Nano_cff import addPh2L1PFCandidates
    addPh2L1PFCandidates(process)
    return _reconcile(process)


def _addOTStubTable(process):
    if not hasattr(process, _OT_STUB_TABLE):
        from DPGAnalysis.Phase3SmartPixelsNanoAOD.l1tPh3SmartPixelsNano_cff import l1tPh3SmartPixelsAllStubTable
        setattr(process, _OT_STUB_TABLE, l1tPh3SmartPixelsAllStubTable.clone())
    _addToNano(process, "p3L1SmartPixelsStubTask", [_OT_STUB_TABLE])


def addTrkStubTables(process):
    """@TrkStubs: every OT stub (L1TOTStub, with truth)."""
    _begin(process, "TrkStubs")
    _addOTStubTable(process)
    return _reconcile(process)


def addSPTrkTables(process):
    """@SPTrk: L1TSmartPixelsTrack / L1TSmartPixelsExtTrack of the wired variant."""
    _begin(process, "SPTrk")
    _addSPTrackTables(process)
    return _reconcile(process)


def _addSPTrackTables(process):
    from DPGAnalysis.Phase3SmartPixelsNanoAOD.l1tPh3SmartPixelsNano_cff import (
        l1tPh3SmartPixelsTracksTable, l1tPh3ExtSmartPixelsTracksTable)
    for side, template in (("prompt", l1tPh3SmartPixelsTracksTable), ("ext", l1tPh3ExtSmartPixelsTracksTable)):
        label = _SP_TABLES[side][0]
        if not hasattr(process, label):
            setattr(process, label, template.clone())
    _addToNano(process, "p3L1SmartPixelsTracksTask", [label for label, _ in _SP_TABLES.values()])


def addSPClusterStubTables(process):
    """@SPClustersStubs: every OT stub (L1TOTStub) and every IT cluster (L1TSmartPixelsCluster)."""
    _begin(process, "SPClustersStubs")
    _addOTStubTable(process)
    if not hasattr(process, _CLUSTER_TABLE):
        from DPGAnalysis.Phase3SmartPixelsNanoAOD.l1tPh3SmartPixelsNano_cff import l1tPh3SmartPixelsClusterTable
        setattr(process, _CLUSTER_TABLE, l1tPh3SmartPixelsClusterTable.clone())
    _addToNano(process, "p3L1SmartPixelsClusterTask", [_CLUSTER_TABLE])
    return _reconcile(process)


def _useModeTag(process, tag):
    st = _begin(process, tag)
    if st.modeTag.value():
        raise RuntimeError(f"SmartPixels nano: two mode tags, @{st.modeTag.value()} and @{tag}; compose one")
    st.modeTag = tag
    return _reconcile(process)


def useSPSeedlessRefit(process):
    """@SPSeedlessRefit: the seedlessRefit variant (OT track + IT hits, Kalman refit); adds its
    per-crossing refit hit tables and spix* columns on the SmartPixels track tables."""
    return _useModeTag(process, "SPSeedlessRefit")


def useSPSeedlessRefitV3(process):
    """@SPSeedlessRefitV3: the seedless refit v3 (the replay's calibrated arm with the two-layer hit choice,
    doc/RefitV3Design.md); the same tables as @SPSeedlessRefit, with the v3 columns filled."""
    return _useModeTag(process, "SPSeedlessRefitV3")


# The tagger pseudo-nano tier (the physicist's sign-off, 2026-10-21): a named tier that keeps every track,
# refit, jet, PF/Puppi, GT, generator and remaining table and drops only the bulk truth/cluster/stub tables
# below. The canonical tiers keep everything; the truth the tagger needs travels with the tracks (otTp*,
# selTpIdx; _configureLinkedTruth). The full-float projOT* refit-hit columns stay (the physicist, 2026-10-21:
# about 0.36 MB/event for one mask). README.md, "Tagger pseudo-nano tier".
TAGGER_TIER = "SPTaggerTier"
TAGGER_TIER_DROPPED = ("L1TSmartPixelsCluster", "L1TTP", "L1TOTStub", "L1TOTStubDTC")


def useSPTaggerTier(process):
    """@SPTaggerTier: the tagger pseudo-nano tier: the output drops L1TSmartPixelsCluster, L1TTP, L1TOTStub
    and L1TOTStubDTC (their producers still run when another table needs them). Compose it last or anywhere:
    every reconcile re-derives the drops from the job's tables."""
    st = _begin(process, TAGGER_TIER)
    st.tier = TAGGER_TIER
    return _reconcile(process)


def _applyTier(process, st):
    if not st.tier.value():
        return
    labels = sorted(label for label, module in process.producers_().items()
                    if any(isinstance(getattr(module, p, None), cms.string)
                           and getattr(module, p).value() in TAGGER_TIER_DROPPED for p in ("name", "tableName")))
    drops = [f"drop nanoaodFlatTable_{label}_*_*" for label in labels]
    for module in process.outputModules_().values():
        if module.type_() != "NanoAODOutputModule":
            continue
        cmds = list(module.outputCommands)
        module.outputCommands = cms.untracked.vstring(*(cmds + [d for d in drops if d not in cmds]))
    if not hasattr(process, "nanoMetadata"):
        return
    process.nanoMetadata.strings.smartPixelsTier = cms.string(
        f"{st.tier.value()}: drops {', '.join(TAGGER_TIER_DROPPED)}")


# One job, the L1 emulation once (from the OT tracks, as the passthrough control IIII) and the seedless
# refit v3 once per mask, each with mask-suffixed SmartPixels tables (smartPixelsCoexist).
MULTI_MASKS = ("AAAA", "AAAI", "IAAA")


def useSPMultiMaskV3(process):
    """@SPMultiMaskV3: coexisting seedless refit v3 variants for AAAA, AAAI and IAAA with suffixed tables
    (L1TSmartPixelsTrackSeedlessRefitV3<mask>, ...RefitHitSeedlessRefitV3<mask>, links, otlink truth); the
    L1 emulation reads the OT tracks. Compose without @SPTrk and without mode or mask tags (coexist writes
    no canonical SmartPixels tables)."""
    st = _begin(process, "SPMultiMaskV3")
    if st.modeTag.value() or st.mask.value():
        raise RuntimeError("SmartPixels nano: @SPMultiMaskV3 replaces the mode and mask tags; compose neither")
    from L1Trigger.Phase3SmartPixels.customizeSmartPixels_cff import smartPixelsCoexist
    variants = [("seedlessRefitV3", m.replace("A", "1").replace("I", "0")) for m in MULTI_MASKS]
    return smartPixelsCoexist(process, variants=variants, addNanoTables=True,
                              trackInputMode=st.trackInputMode.value(),
                              seedlessRefitConfig=json.loads(st.seedlessRefitConfig.value()),
                              extendedTracks=st.extendedTracks.value(), promptHnpar=st.promptHnpar.value())


# The tagger's per-track heads (heads_cmssw/SPEC.md): six heads, one conifer model per head and window.
TRACK_HEADS = (("clean", "probability"), ("d0", "margin"), ("z0", "margin"), ("tanl", "margin"),
               ("rinv", "margin"), ("phi", "margin"))


def trackHeadModelFile(modelDir, coll, mask, window, head, flat=False):
    """<modelDir>/<OT|SmartPixels>/<mask>/<window>/<head>.json, or <modelDir>/<coll>_<mask>_<window>_<head>.json
    (flat: one directory, e.g. the files a grid job receives in its working directory)."""
    return f"{modelDir}/{coll}_{mask}_{window}_{head}.json" if flat else f"{modelDir}/{coll}/{mask}/{window}/{head}.json"


def addSmartPixelsTrackHeads(process, modelDir, heads=TRACK_HEADS, columnPrefix="head_", writeFeatures=False,
                             flat=False):
    """The per-track heads on the extended (5-parameter) track tables of the job: the OT models on L1TExtTrack
    and, on every SmartPixels extended table (canonical and coexist), the models of its variant's mask
    (trackHeadModelFile). Call after the add-on flavours, e.g. in --customise_commands. The prompt
    (4-parameter) family has no models."""
    from DPGAnalysis.Phase3SmartPixelsNanoAOD.l1tPh3SmartPixelsNano_cff import l1tSmartPixelsTrackHeadsTable
    st = _state(process)
    sched = _scheduled(process)

    def models(coll, mask, windows):
        return [cms.PSet(name=cms.string(h), output=cms.string(out), models=cms.VPSet(
            *[cms.PSet(window=cms.string(w), file=cms.string(trackHeadModelFile(modelDir, coll, mask, w, h, flat)))
              for w in windows])) for h, out in heads]

    otTables = cms.VInputTag(cms.InputTag("l1tExtTracksTable"), cms.InputTag("l1tExtTracksTruthTable"))
    added = []
    if "l1tExtTracksTable" in sched:
        setattr(process, "l1tExtTracksHeadsTable", l1tSmartPixelsTrackHeadsTable.clone(
            collection="OT", tableName="L1TExtTrack", columnPrefix=columnPrefix, otTables=otTables,
            heads=cms.VPSet(*models("OT", "IIII", ("all",))), writeFeatures=writeFeatures))
        added.append("l1tExtTracksHeadsTable")
    variants = [("", _wiredVariant(st))] if st.variant.value() else []
    from L1Trigger.Phase3SmartPixels.customizeSmartPixels_cff import smartPixelsVariantSuffix
    for spec in st.coexistVariants:
        mode, _, activeSP = spec.partition(":")
        suffix = smartPixelsVariantSuffix(mode, activeSP or None)
        variants.append((suffix[:1].upper() + suffix[1:], (mode, activeSP or None)))
    for suffix, (mode, activeSP) in variants:
        if mode not in _MODE_TAGS.values() or not activeSP:
            continue
        mask = activeSP.replace("1", "A").replace("0", "I")
        table, refit = f"l1tPh3ExtSmartPixelsTracks{suffix}Table", f"l1tPh3ExtSmartPixelsRefit{suffix}Table"
        if table not in sched or refit not in sched:
            continue
        label = f"l1tPh3ExtSmartPixelsHeads{suffix}Table"
        setattr(process, label, l1tSmartPixelsTrackHeadsTable.clone(
            collection="SmartPixels", tableName=f"L1TSmartPixelsExtTrack{suffix}", columnPrefix=columnPrefix,
            otTables=otTables, spTables=cms.VInputTag(cms.InputTag(table), cms.InputTag(refit, "trk")),
            hitTable=cms.InputTag(refit, "hit"), heads=cms.VPSet(*models("SmartPixels", mask, ("prompt", "displaced"))),
            writeFeatures=writeFeatures))
        added.append(label)
    _addToNano(process, "p3L1SmartPixelsHeadsTask", added)
    return process


def _reserved(tag):
    def customise(process):
        raise NotImplementedError(f"SmartPixels nano: @{tag} is RESERVED; the mode is not implemented. "
                                  "Use @SPSeedlessRefit.")
    customise.__name__ = f"use{tag}"
    return customise


useSPParallelSeedRefit = _reserved("SPParallelSeedRefit")
useSPJointSeedFit = _reserved("SPJointSeedFit")


def _maskTag(mask):
    def customise(process):
        st = _begin(process, f"SPActiveLayers{mask}")
        if st.mask.value():
            raise RuntimeError(f"SmartPixels nano: two mask tags, @SPActiveLayers{st.mask.value()} and "
                               f"@SPActiveLayers{mask}; compose one")
        st.mask = mask
        return _reconcile(process)
    customise.__name__ = f"useSPActiveLayers{mask}"
    customise.__doc__ = (f"@SPActiveLayers{mask}: active IT layers of the wired variant "
                         "(character k = ILk and IDk; A active, I inactive)"
                         + ("; IIII = passthrough control" if mask == PASSTHROUGH_MASK else ""))
    return customise


for _mask in MASKS:
    globals()[f"useSPActiveLayers{_mask}"] = _maskTag(_mask)


def _tpTruthTag(tag):
    def customise(process):
        st = _begin(process, tag)
        if st.tpTruth.value():
            raise RuntimeError(f"SmartPixels nano: TrackingParticle truth flavours are mutually exclusive; "
                               f"got @{st.tpTruth.value()} and @{tag}")
        st.tpTruth = tag
        from DPGAnalysis.Phase2L1TNanoAOD.l1tPh2Nano_cff import (
            addPh2L1TrackingParticles, addPh2L1TrackingParticlesWithGen)
        (addPh2L1TrackingParticlesWithGen if _withGen(process) else addPh2L1TrackingParticles)(process)
        return _reconcile(process)
    customise.__name__ = f"use{tag}"
    return customise


for _tag in _TP_TRUTH_TAGS:
    globals()[f"use{_tag}"] = _tpTruthTag(_tag)


# ---------------------------------------------------------------------------
# Settings that are not in a tag name
# ---------------------------------------------------------------------------
_POSTURE_KEYS = ("trackInputMode", "extendedTracks", "promptHnpar")


def smartPixelsNanoParams(process, **params):
    """Settings of the SmartPixels nano that no tag names:
      trackInputMode, extendedTracks, promptHnpar  (see smartPixelsCoopt), and
      any seedlessRefitConfig key (SEEDLESSREFIT_DEFAULTS), e.g. angleSource (default "emu"),
      pixelavAngleSet / smarthitFakeSet (the payload names for angleSource="payload").

    cmsDriver runs --customise_commands AFTER the autoNANO customisations, so the variant is
    usually wired before this runs. The settings are then applied to what the tags built,
    which yields the same configuration as passing them first (hand-written configs may
    call this before the add-ons). trackInputMode leaves the default only once.

      --customise_commands "from DPGAnalysis.Phase3SmartPixelsNanoAOD.l1tPh3SmartPixelsNanoFlavours_cff import smartPixelsNanoParams; smartPixelsNanoParams(process, trackInputMode='rebuildTracksFromStubs')"
    """
    from L1Trigger.Phase3SmartPixels.customizeSmartPixels_cff import (
        SEEDLESSREFIT_DEFAULTS, _resolveSeedlessRefitConfig, _resolveTrackInputMode)
    st = _state(process)
    unknown = set(params) - set(_POSTURE_KEYS) - set(SEEDLESSREFIT_DEFAULTS)
    if unknown:
        raise ValueError(f"smartPixelsNanoParams: unknown setting(s) {sorted(unknown)}")
    if "trackInputMode" in params:
        mode = _resolveTrackInputMode(params["trackInputMode"])
        current = st.trackInputMode.value()
        if mode != current and current != "reemulateL1TrackFinding":
            raise RuntimeError(f"smartPixelsNanoParams: trackInputMode is already {current!r}")
        st.trackInputMode = mode
    if "extendedTracks" in params:
        st.extendedTracks = bool(params["extendedTracks"])
    if "promptHnpar" in params:
        st.promptHnpar = int(params["promptHnpar"])
    overrides = json.loads(st.seedlessRefitConfig.value())
    overrides.update({k: v for k, v in params.items() if k in SEEDLESSREFIT_DEFAULTS})
    _resolveSeedlessRefitConfig(overrides)
    st.seedlessRefitConfig = json.dumps(overrides, sort_keys=True)
    _applySeedlessRefitSettings(process, st)
    if list(st.addOns):
        _reconcile(process)
    return process


def _applySeedlessRefitSettings(process, st):
    from L1Trigger.Phase3SmartPixels.customizeSmartPixels_cff import (
        _applySeedlessRefitConfig, configureSmartPixelsRecHits, smartPixelsVariantLabels)
    resolved = _seedlessRefitConfig(st)
    from L1Trigger.Phase3SmartPixels.customizeSmartPixels_cff import SEEDLESSREFIT_PRESETS
    for spec in [st.variant.value()] + list(st.coexistVariants):
        mode, _, activeSP = spec.partition(":")
        if mode not in SEEDLESSREFIT_PRESETS:
            continue
        if resolved["angleSource"] == "payload" and not resolved["pixelavAngleSet"]:
            raise ValueError("smartPixelsNanoParams: angleSource='payload' needs pixelavAngleSet")
        for label in smartPixelsVariantLabels(mode, activeSP or None):
            _applySeedlessRefitConfig(getattr(process, label), resolved, SEEDLESSREFIT_PRESETS[mode])
    if hasattr(process, "spixSmartPixelsRecHits"):
        configureSmartPixelsRecHits(process, resolved["pixelavAngleSet"], resolved["smarthitFakeSet"],
                                    angleSource=resolved["angleSource"], discScope=resolved["discScope"])


# ---------------------------------------------------------------------------
# Links and wiring, derived from the current content of the job
# ---------------------------------------------------------------------------
def _reconcile(process):
    from DPGAnalysis.Phase3SmartPixelsNanoAOD.l1tPh3SmartPixelsNano_cff import useSmartPixelsNanoCompression
    useSmartPixelsNanoCompression(process)
    st = _state(process)
    _checkRecoInputs(process, st)
    _wireVariant(process, st)
    _applyPosture(process, st)
    sched = _scheduled(process)
    variant = _wiredVariant(st)
    _configureSmartPixelsTracks(process, st, variant, sched)
    _configureModeTables(process, variant, sched)
    _configureStubLinks(process, st, variant, sched)
    if _CLUSTER_TABLE in sched:
        _ensureRecHits(process, st)
    _configureCandidateLinks(process, st, variant, sched)
    _configureDisplacedVertexLink(process, variant, sched)
    _configureGTTColumns(process, variant, sched)
    _configureTrackingParticles(process, st, variant, sched)
    _configureOTStubDTC(process, sched)
    _configureOTSeedType(process, sched)
    _configureLinkedTruth(process, sched)
    _applyTier(process, st)
    return process


_OT_TRUTH_TABLES = {"L1TTrack": "l1tTracksTruthTable", "L1TExtTrack": "l1tExtTracksTruthTable"}


def _configureLinkedTruth(process, sched):
    """Truth that travels with the tracks, so a tier without L1TTP / L1TSmartPixelsCluster keeps it: every
    SmartPixels track table (canonical and coexist) copies its OT track's TP truth (otTp*, the otlink
    convention) when the OT truth table is in the job, and every refit hit table carries the selected
    cluster's tpIdx (selTpIdx) when the cluster table has its truth columns."""
    for label, module in process.producers_().items():
        if label not in sched:
            continue
        if module.type_() == "L1SmartPixelsOTTrackLinkTableProducer":
            truth = _OT_TRUTH_TABLES.get(module.otTrackTableName.value())
            module.otTruthTable = cms.InputTag(truth if truth in sched else "")
        elif module.type_() == "L1SmartPixelsRefitTableProducer":
            withTruth = (_CLUSTER_TABLE in sched and module.clusterTableName.value()
                         and getattr(process, _CLUSTER_TABLE).doTruth.value())
            module.clusterTable = cms.InputTag(_CLUSTER_TABLE if withTruth else "")


def _wireVariant(process, st):
    selected, wired = _variant(st), _wiredVariant(st)
    if selected is not None and wired is None:
        from L1Trigger.Phase3SmartPixels.customizeSmartPixels_cff import smartPixelsCoopt
        smartPixelsCoopt(process, variant=selected, addNanoTables=False,
                         trackInputMode=st.trackInputMode.value(), extendedTracks=st.extendedTracks.value(),
                         promptHnpar=st.promptHnpar.value(), seedlessRefitConfig=_seedlessRefitConfig(st))
        st.variant = _variantString(selected)
        st.trackInputModeApplied = st.trackInputMode.value()
        return
    if selected is not None and selected != wired:
        raise RuntimeError(f"SmartPixels nano: variant {st.variant.value()} is already wired, "
                           f"the tags select {_variantString(selected)}")
    if wired is not None:
        # consumers added after the wiring (e.g. the displaced-vertex producer) are rewired too
        from L1Trigger.Phase3SmartPixels.customizeSmartPixels_cff import (
            injectSmartPixelsTrackProducer, smartPixelsVariantLabels)
        prompt, extended = smartPixelsVariantLabels(*wired)
        injectSmartPixelsTrackProducer(process, addAssociation=False, l1tSmartPixelsTrackProducerLabel=prompt,
                                       l1tSmartPixelsTrackProducerExtendedLabel=extended,
                                       skipModuleLabels=[conv for conv, _ in _OT_GTT_CONVERSIONS.values()])


_FILE_ASSOCIATORS = {"useStoredTracks": None,  # all of TRUTH_ASSOCIATOR_LABELS
                     "rebuildTracksFromStubs": ("TTClusterAssociatorFromPixelDigis", "TTStubAssociatorFromPixelDigis")}


def _applyPosture(process, st):
    """trackInputMode applies to the whole job, SmartPixels or not. After it is applied, every
    reconcile removes again the associators whose maps must come from the input, since a
    later add-on may have loaded them."""
    from L1Trigger.Phase3SmartPixels.customizeSmartPixels_cff import (
        TRUTH_ASSOCIATOR_LABELS, _applyTrackInputMode, useTruthAssociationFromFile)
    mode = st.trackInputMode.value()
    if mode == "reemulateL1TrackFinding":
        return
    if st.trackInputModeApplied.value() != mode:
        _applyTrackInputMode(process, mode, [], extendedTracks=st.extendedTracks.value(),
                             promptHnpar=st.promptHnpar.value())
        st.trackInputModeApplied = mode
    else:
        useTruthAssociationFromFile(process, _FILE_ASSOCIATORS[mode] or TRUTH_ASSOCIATOR_LABELS)


def _ensureRecHits(process, st):
    from L1Trigger.Phase3SmartPixels.customizeSmartPixels_cff import ensureSmartPixelsRecHits
    resolved = _seedlessRefitConfig(st)
    ensureSmartPixelsRecHits(process, resolved["pixelavAngleSet"], resolved["smarthitFakeSet"],
                             angleSource=resolved["angleSource"], discScope=resolved["discScope"])


def _spLabels(variant):
    from L1Trigger.Phase3SmartPixels.customizeSmartPixels_cff import smartPixelsVariantLabels
    prompt, extended = smartPixelsVariantLabels(*variant)
    return {"prompt": prompt, "ext": extended}


def _spTrackTables(st, variant, sched):
    """(variant, side, table label suffix, track table name) of every SmartPixels track table in
    the nano: the canonical pair of the wired variant and the suffixed pairs of coexist variants."""
    from L1Trigger.Phase3SmartPixels.customizeSmartPixels_cff import smartPixelsVariantSuffix
    tables = []
    if variant is not None:
        tables += [(variant, side, "", name) for side, (label, name) in _SP_TABLES.items() if label in sched]
    for spec in st.coexistVariants:
        mode, _, activeSP = spec.partition(":")
        coexist = (mode, activeSP or None)
        suffix = smartPixelsVariantSuffix(*coexist)
        suffix = suffix[:1].upper() + suffix[1:]
        for side, (label, name) in _SP_TABLES.items():
            ext = "Ext" if side == "ext" else ""
            if f"l1tPh3{ext}SmartPixelsTracks{suffix}Table" in sched:
                tables.append((coexist, side, suffix, f"{name}{suffix}"))
    return tables


def _configureSmartPixelsTracks(process, st, variant, sched):
    """SmartPixels track tables read the wired variant; l1TrackIdx links each one to the OT
    track it was built from when that OT table is in the nano."""
    from DPGAnalysis.Phase3SmartPixelsNanoAOD.l1tPh3SmartPixelsNano_cff import l1tPh3SmartPixelsOTTrackLinkTable
    for side, (label, name) in _SP_TABLES.items():
        if label not in sched:
            continue
        if variant is not None:
            getattr(process, label).src = cms.InputTag(_spLabels(variant)[side], "Level1TTTracks")
        linkLabel = f"l1tPh3{'Ext' if side == 'ext' else ''}SmartPixelsOTTrackLinkTable"
        if _OT_TABLES[side][0] in sched and variant is not None and not hasattr(process, linkLabel):
            setattr(process, linkLabel, l1tPh3SmartPixelsOTTrackLinkTable.clone(
                src=cms.InputTag(_spLabels(variant)[side], "Level1TTTracks"),
                otTracks=cms.InputTag(_OT_TRACKS[side], "Level1TTTracks"),
                trackTableName=name,
                otTrackTableName=_OT_TABLES[side][1]))
            _addToNano(process, "p3L1SmartPixelsLinkTask", [linkLabel])


def _configureModeTables(process, variant, sched):
    """The seedlessRefit variant's own tables: per-crossing refit hits (trackIdx -> SmartPixels
    track, selClusterIdx -> L1TSmartPixelsCluster when those tables are present) and the spix*
    columns on the SmartPixels track tables."""
    if variant is None or variant[0] not in ("seedlessRefit", "seedlessRefitV3"):
        return
    from DPGAnalysis.Phase3SmartPixelsNanoAOD.l1tPh3SmartPixelsNano_cff import l1tPh3SmartPixelsRefitTable
    for side, (spLabel, spName) in _SP_TABLES.items():
        ext = "Ext" if side == "ext" else ""
        label = f"l1tPh3{ext}SmartPixelsRefitTable"
        if not hasattr(process, label):
            setattr(process, label, l1tPh3SmartPixelsRefitTable.clone(
                tracks=cms.InputTag(_spLabels(variant)[side], "Level1TTTracks"),
                sidecar=cms.InputTag(_spLabels(variant)[side], "Level1TTTracks"),
                hitTableName=f"L1TSmartPixels{ext}RefitHit"))
            _addToNano(process, "p3L1SmartPixelsRefitTask", [label])
        module = getattr(process, label)
        module.trackTableName = spName if spLabel in sched else ""
        module.clusterTableName = "L1TSmartPixelsCluster" if _CLUSTER_TABLE in sched else ""


def _configureStubLinks(process, st, variant, sched):
    """Link-only track -> stub tables (<track table>Stub: trackIdx, otStubIdx -> L1TOTStub)
    for the OT and the SmartPixels track tables (canonical and coexist), and SmartPixels
    track -> cluster tables (<track table>Cluster: trackIdx, clusterIdx -> L1TSmartPixelsCluster;
    for seedlessRefit the hits accepted into the refit), each when both of its tables are in
    the nano."""
    from DPGAnalysis.Phase3SmartPixelsNanoAOD.l1tPh3SmartPixelsNano_cff import (
        addPh3L1SmartPixelsStubTables, l1tPh3TrackStubLinkTable, l1tPh3SmartPixelsTrackClusterLinkTable)
    if _OT_STUB_TABLE in sched:
        doPrompt, doExtended = (_OT_TABLES[s][0] in sched for s in ("prompt", "ext"))
        if doPrompt or doExtended:
            addPh3L1SmartPixelsStubTables(process, doPrompt=doPrompt, doExtended=doExtended)
    for spVariant, side, suffix, tableName in _spTrackTables(st, variant, sched):
        ext = "Ext" if side == "ext" else ""
        tracks = cms.InputTag(_spLabels(spVariant)[side], "Level1TTTracks")
        if _OT_STUB_TABLE in sched:
            label = f"l1tPh3{ext}SmartPixelsTrackStubLink{suffix}Table"
            if not hasattr(process, label):
                setattr(process, label, l1tPh3TrackStubLinkTable.clone(
                    tracks=tracks, linkTableName=f"{tableName}Stub", trackTableName=tableName))
                _addToNano(process, "p3L1SmartPixelsLinkTask", [label])
        if _CLUSTER_TABLE in sched:
            label = f"l1tPh3{ext}SmartPixelsTrackClusterLink{suffix}Table"
            if not hasattr(process, label):
                setattr(process, label, l1tPh3SmartPixelsTrackClusterLinkTable.clone(
                    tracks=tracks, linkTableName=f"{tableName}Cluster", trackTableName=tableName,
                    sidecar=tracks if spVariant[0] in ("seedlessRefit", "seedlessRefitV3") else cms.InputTag("")))
                _addToNano(process, "p3L1SmartPixelsLinkTask", [label])


def _configureCandidateLinks(process, st, variant, sched):
    """Candidate l1TrackIdx -> the tracks the correlator was built from: the SmartPixels
    tracks when a variant is wired, else the OT tracks; no column without that table.
    Candidate track truth (stock withGen PF flavours) uses the associator of those tracks."""
    from L1Trigger.Phase3SmartPixels.customizeSmartPixels_cff import smartPixelsTrackAssociator
    for label, module in process.producers_().items():
        if label not in sched or module.type_() not in ("L1JetCandLinkTableProducer",
                                                        "L1PFCandTrackTruthTableProducer"):
            continue
        side = _CAND_SIDE[module.candTableName.value()]
        if variant is not None:
            tracks, (tableLabel, tableName) = _spLabels(variant)[side], _SP_TABLES[side]
        else:
            tracks, (tableLabel, tableName) = _OT_TRACKS[side], _OT_TABLES[side]
        module.tracks = cms.InputTag(tracks, "Level1TTTracks")
        if module.type_() == "L1PFCandTrackTruthTableProducer":
            associator = (smartPixelsTrackAssociator(process, tracks, side == "ext") if variant is not None
                          else _OT_ASSOCIATORS[side])
            module.truthAssociation = cms.InputTag(associator, "Level1TTTracks")
            module.trackTableName = tableName
        else:
            module.trackTableName = tableName if tableLabel in sched else ""


def _configureDisplacedVertexLink(process, variant, sched):
    """L1DispVertex firstIndexTrk/secondIndexTrk -> the extended tracks the vertex finder
    ran on (SmartPixels when a variant is wired), when that table is in the nano."""
    if "dispVtxTable" not in sched:
        return
    tableLabel, tableName = _SP_TABLES["ext"] if variant is not None else _OT_TABLES["ext"]
    if tableLabel not in sched:
        # the link was made before a variant was wired and its target is not in the nano
        if process.p2L1DisplacedVertexTask in process.l1tPh2NanoTask._collection:
            process.l1tPh2NanoTask.remove(process.p2L1DisplacedVertexTask)
        return
    from DPGAnalysis.Phase2L1TNanoAOD.l1tPh2Nano_cff import extendPh2L1DisplacedVertices
    extendPh2L1DisplacedVertices(process)
    for col, which in (("firstIndexTrk", "higher"), ("secondIndexTrk", "lower")):
        getattr(process.l1tDisplacedVertexTable.variables, col).doc = cms.string(
            f"index of the {which}-pt track in the {tableName} table")


# GTT extension table -> (GTT input conversion of the OT tracks, its output instance)
_OT_GTT_CONVERSIONS = {"gttTracksTable": ("l1tGTTInputProducerOT", "Level1TTTracksConverted"),
                       "gttExtTracksTable": ("l1tGTTInputProducerExtendedOT", "Level1TTTracksExtendedConverted")}


def _configureGTTColumns(process, variant, sched):
    """The GTT hardware columns of L1TTrack / L1TExtTrack describe the OT tracks. With a
    variant wired the L1 GTT input converts the SmartPixels tracks, so the OT tables get a
    GTT input conversion of their own."""
    if variant is None:
        return
    from L1Trigger.L1TTrackMatch.l1tGTTInputProducer_cfi import l1tGTTInputProducer, l1tGTTInputProducerExtended
    for table, cfi, side in (("gttTracksTable", l1tGTTInputProducer, "prompt"),
                             ("gttExtTracksTable", l1tGTTInputProducerExtended, "ext")):
        conv, instance = _OT_GTT_CONVERSIONS[table]
        if table not in sched:
            continue
        if not hasattr(process, conv):
            setattr(process, conv, cfi.clone(l1TracksInputTag=cms.InputTag(_OT_TRACKS[side], "Level1TTTracks")))
            _addToNano(process, "smartPixelsOTGTTConversionTask", [conv])
        getattr(process, table).src = cms.InputTag(conv, instance)


def _configureOTStubDTC(process, sched):
    """L1TOTStubDTC next to L1TOTStub whenever the prompt tracklet chain runs in the job: it reads
    the chain's TTDTC, which is transient and exists only there."""
    if _OT_STUB_TABLE not in sched or not hasattr(process, _OT_TRACKS["prompt"]):
        return
    dtc = getattr(process, _OT_TRACKS["prompt"]).InputTagTTDTC
    if not hasattr(process, dtc.getModuleLabel()):
        return
    if not hasattr(process, _OT_DTC_TABLE):
        from DPGAnalysis.Phase3SmartPixelsNanoAOD.l1tPh3SmartPixelsNano_cff import l1tPh3SmartPixelsDTCStubTable
        setattr(process, _OT_DTC_TABLE, l1tPh3SmartPixelsDTCStubTable.clone(
            dtc=dtc, stubs=getattr(process, _OT_STUB_TABLE).stubs))
    _addToNano(process, "p3L1SmartPixelsStubTask", [_OT_DTC_TABLE])


def _configureOTSeedType(process, sched):
    """The tracklet seed type of every OT track (TTTrack::trackSeedType: L1L2 = 0, L2L3 = 1, L3L4 = 2,
    L5L6 = 3, D1D2 = 4, D3D4 = 5, L1D1 = 6, L2D1 = 7; the extended seeds above)."""
    from PhysicsTools.NanoAOD.common_cff import Var
    for label, _ in _OT_TABLES.values():
        if label in sched and not hasattr(getattr(process, label).variables, "seedType"):
            getattr(process, label).variables.seedType = Var("trackSeedType()", "int", doc="tracklet seed type")


def _configureTrackingParticles(process, st, variant, sched):
    """L1TTP selection of the TrackingParticle flavour. The linked selections take their
    TrackingParticles from the truth associations of the objects the nano writes: OT stubs
    (L1TOTStub; stub association, and the owners of both clusters of each stub), OT tracks,
    SmartPixels tracks (the associator of the wired variant) and IT clusters (their dominant
    TrackingParticle)."""
    if _TP_TABLE not in sched or not st.tpTruth.value():
        return
    selection = _TP_TRUTH_TAGS[st.tpTruth.value()]
    module = getattr(process, _TP_TABLE)
    module.selection = selection
    if selection not in ("linked", "charged1GeVOrLinked"):
        return
    from L1Trigger.Phase3SmartPixels.customizeSmartPixels_cff import smartPixelsTrackAssociator
    stubs, tracks, tps = [], [], []
    if _OT_STUB_TABLE in sched:
        stubTable = getattr(process, _OT_STUB_TABLE)
        stubs.append(stubTable.stubTruth)
        if not hasattr(process, "l1tPh3SmartPixelsStubClusterTPs"):
            process.l1tPh3SmartPixelsStubClusterTPs = cms.EDProducer(
                "L1SmartPixelsStubClusterTPLinkProducer",
                stubs=stubTable.stubs, clusterTruth=stubTable.clusterTruth,
                trackingParticles=module.trackingParticles)
            _addToNano(process, "p3L1SmartPixelsStubTask", ["l1tPh3SmartPixelsStubClusterTPs"])
        tps.append(cms.InputTag("l1tPh3SmartPixelsStubClusterTPs"))
    for side in ("prompt", "ext"):
        if _OT_TABLES[side][0] in sched:
            tracks.append(cms.InputTag(_OT_ASSOCIATORS[side], "Level1TTTracks"))
        if variant is not None and _SP_TABLES[side][0] in sched:
            tracks.append(cms.InputTag(smartPixelsTrackAssociator(process, _spLabels(variant)[side], side == "ext"),
                                       "Level1TTTracks"))
    if _CLUSTER_TABLE in sched:
        if not hasattr(process, "l1tPh3SmartPixelsClusterTPs"):
            process.l1tPh3SmartPixelsClusterTPs = cms.EDProducer(
                "L1SmartPixelsClusterTPLinkProducer",
                smartPixelsRecHits=getattr(process, _CLUSTER_TABLE).smartPixelsRecHits)
            _addToNano(process, "p3L1SmartPixelsClusterTask", ["l1tPh3SmartPixelsClusterTPs"])
        tps.append(cms.InputTag("l1tPh3SmartPixelsClusterTPs"))
    module.linkedStubAssociations = cms.VInputTag(*stubs)
    module.linkedTrackAssociations = cms.VInputTag(*tracks)
    module.linkedTrackingParticles = cms.VInputTag(*tps)


# ---------------------------------------------------------------------------
# @PHYS (the regular nano: HLT, PUPPI MET, AK4 jets, muons, taus, ...) reads MINIAOD
# products. They come from the input (GEN-SIM-DIGI-RAW-MINIAOD) or from RECO+PAT run in
# the same job; today's RelVals (GEN-SIM-DIGI-RAW) need the latter.
# ---------------------------------------------------------------------------
_MINIAOD_BRANCHES = ("patPackedCandidates_packedPFCandidates__", "patJets_slimmedJets__",
                     "patMuons_slimmedMuons__", "recoGenParticles_prunedGenParticles__")
_RECO_STEPS = "RAW2DIGI,RECO,RECOSIM,PAT"


def inputHasMiniAOD(fileName):
    """True if the first file of the input carries the MINIAOD products @PHYS reads."""
    from DPGAnalysis.Phase3SmartPixelsNanoAOD.l1tPh3SmartPixelsNano_cff import inputEventBranches
    branches = inputEventBranches(fileName)
    return all(any(b.startswith(p) for b in branches) for p in _MINIAOD_BRANCHES)


def _pathLabels(path):
    labels = set()

    class Visitor:
        def enter(self, x):
            if hasattr(x, "label_"):
                try:
                    labels.add(x.label_())
                except Exception:
                    pass

        def leave(self, x):
            pass
    path.visit(Visitor())
    return labels


def _checkRecoInputs(process, st):
    """With @PHYS composed, fail at configuration time unless the input has the MINIAOD
    products or the job runs RECO and PAT."""
    if not hasattr(process, "nanoAOD_step") or not (
            _pathLabels(process.nanoAOD_step) & {"nanoSequence", "nanoSequenceMC", "nanoSequenceFS"}):
        return
    recoInJob = hasattr(process, "reconstruction_step") and hasattr(process, "patTask") and (
        process.schedule is None or any(t is process.patTask for t in process.schedule._tasks))
    if recoInJob:
        return
    fileNames = list(getattr(process.source, "fileNames", []))
    if not fileNames:
        raise RuntimeError("SmartPixels nano: @PHYS needs MINIAOD products; the job has no input file to check "
                           "and runs no RECO+PAT")
    if inputHasMiniAOD(fileNames[0]):
        return
    # the NANO spec as given to cmsDriver; the add-ons seen so far otherwise (e.g. under cmsRun)
    import sys
    given = [arg.split("NANO:", 1)[1].split(",")[0] for arg in sys.argv if "NANO:" in arg]
    spec = given[0] if given else "+".join(["@Phase2L1DPGwithGen" if _withGen(process) else "@Phase2L1DPG"]
                                           + ["@" + tag for tag in st.addOns if tag != "coopt"] + ["@PHYS"])
    raise RuntimeError(
        f"SmartPixels nano: @PHYS reads MINIAOD products ({', '.join(b.split('_')[1] for b in _MINIAOD_BRANCHES)}) "
        f"that {fileNames[0]} does not carry and this job does not produce. Run RECO and PAT in the job, "
        f"with the RelVal workflow's own step-3 options: -s {_RECO_STEPS},NANO:{spec} "
        "(L1Trigger/Phase3SmartPixels/test/makeSpixConfig.py --reco emits the full command)")
