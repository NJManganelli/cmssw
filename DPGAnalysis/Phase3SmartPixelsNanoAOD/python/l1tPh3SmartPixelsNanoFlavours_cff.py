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
# The SmartPixels variant is chosen by one mode tag (@SPSeedlessRefit = digiRefit) and one
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

# digiRefitConfig overrides the add-on flavours apply over DIGIREFIT_DEFAULTS.
_TIER_DIGIREFIT_DEFAULTS = {"angleSource": "emu"}

_MODE_TAGS = {"SPSeedlessRefit": "digiRefit"}
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
            trackInputMode=cms.string("reemulateL1TrackFinding"),
            trackInputModeApplied=cms.string(""),
            extendedTracks=cms.bool(True),
            promptHnpar=cms.uint32(5),
            digiRefitConfig=cms.string(json.dumps(_TIER_DIGIREFIT_DEFAULTS, sort_keys=True)),
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


def _digiRefitConfig(st):
    from L1Trigger.Phase3SmartPixels.customizeSmartPixels_cff import _resolveDigiRefitConfig
    return _resolveDigiRefitConfig(json.loads(st.digiRefitConfig.value()))


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
    """@SPSeedlessRefit: the digiRefit variant (OT track + IT hits, Kalman refit); adds its
    per-crossing refit hit tables and spix* columns on the SmartPixels track tables."""
    return _useModeTag(process, "SPSeedlessRefit")


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
      any digiRefitConfig key (DIGIREFIT_DEFAULTS), e.g. angleSource (default "emu"),
      pixelavAngleSet / smarthitFakeSet (the payload names for angleSource="payload").

    cmsDriver runs --customise_commands AFTER the autoNANO customisations, so the variant is
    usually wired before this runs. The settings are then applied to what the tags built,
    which yields the same configuration as passing them first (hand-written configs may
    call this before the add-ons). trackInputMode leaves the default only once.

      --customise_commands "from DPGAnalysis.Phase3SmartPixelsNanoAOD.l1tPh3SmartPixelsNanoFlavours_cff import smartPixelsNanoParams; smartPixelsNanoParams(process, trackInputMode='rebuildTracksFromStubs')"
    """
    from L1Trigger.Phase3SmartPixels.customizeSmartPixels_cff import (
        DIGIREFIT_DEFAULTS, _resolveDigiRefitConfig, _resolveTrackInputMode)
    st = _state(process)
    unknown = set(params) - set(_POSTURE_KEYS) - set(DIGIREFIT_DEFAULTS)
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
    overrides = json.loads(st.digiRefitConfig.value())
    overrides.update({k: v for k, v in params.items() if k in DIGIREFIT_DEFAULTS})
    _resolveDigiRefitConfig(overrides)
    st.digiRefitConfig = json.dumps(overrides, sort_keys=True)
    _applyDigiRefitSettings(process, st)
    if list(st.addOns):
        _reconcile(process)
    return process


def _applyDigiRefitSettings(process, st):
    from L1Trigger.Phase3SmartPixels.customizeSmartPixels_cff import (
        _applyDigiRefitConfig, configureSmartPixelsRecHits, smartPixelsVariantLabels)
    resolved = _digiRefitConfig(st)
    variant = st.variant.value()
    if variant.startswith("digiRefit:"):
        if resolved["angleSource"] == "payload" and not resolved["pixelavAngleSet"]:
            raise ValueError("smartPixelsNanoParams: angleSource='payload' needs pixelavAngleSet")
        for label in smartPixelsVariantLabels("digiRefit", variant.split(":")[1]):
            _applyDigiRefitConfig(getattr(process, label), resolved)
    if hasattr(process, "spixSmartPixelsRecHits"):
        configureSmartPixelsRecHits(process, resolved["pixelavAngleSet"], resolved["smarthitFakeSet"],
                                    angleSource=resolved["angleSource"], discScope=resolved["discScope"])


# ---------------------------------------------------------------------------
# Links and wiring, derived from the current content of the job
# ---------------------------------------------------------------------------
def _reconcile(process):
    st = _state(process)
    _wireVariant(process, st)
    _applyPosture(process, st)
    sched = _scheduled(process)
    variant = _wiredVariant(st)
    _configureSmartPixelsTracks(process, st, variant, sched)
    _configureModeTables(process, variant, sched)
    _configureStubLinks(process, variant, sched)
    if _CLUSTER_TABLE in sched:
        _ensureRecHits(process, st)
    _configureCandidateLinks(process, st, variant, sched)
    _configureDisplacedVertexLink(process, variant, sched)
    _configureGTTColumns(process, variant, sched)
    _configureTrackingParticles(process, st, variant, sched)
    return process


def _wireVariant(process, st):
    selected, wired = _variant(st), _wiredVariant(st)
    if selected is not None and wired is None:
        from L1Trigger.Phase3SmartPixels.customizeSmartPixels_cff import smartPixelsCoopt
        smartPixelsCoopt(process, variant=selected, addNanoTables=False,
                         trackInputMode=st.trackInputMode.value(), extendedTracks=st.extendedTracks.value(),
                         promptHnpar=st.promptHnpar.value(), digiRefitConfig=_digiRefitConfig(st))
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
    resolved = _digiRefitConfig(st)
    ensureSmartPixelsRecHits(process, resolved["pixelavAngleSet"], resolved["smarthitFakeSet"],
                             angleSource=resolved["angleSource"], discScope=resolved["discScope"])


def _spLabels(variant):
    from L1Trigger.Phase3SmartPixels.customizeSmartPixels_cff import smartPixelsVariantLabels
    prompt, extended = smartPixelsVariantLabels(*variant)
    return {"prompt": prompt, "ext": extended}


def _configureSmartPixelsTracks(process, st, variant, sched):
    """SmartPixels track tables read the wired variant; l1TrackIdx links each one to the OT
    track it was built from when that OT table is in the nano."""
    from DPGAnalysis.Phase3SmartPixelsNanoAOD.l1tPh3SmartPixelsNano_cff import l1tPh3SmartPixelsSeedLinkTable
    for side, (label, name) in _SP_TABLES.items():
        if label not in sched:
            continue
        if variant is not None:
            getattr(process, label).src = cms.InputTag(_spLabels(variant)[side], "Level1TTTracks")
        linkLabel = f"l1tPh3{'Ext' if side == 'ext' else ''}SmartPixelsSeedLinkTable"
        if _OT_TABLES[side][0] in sched and variant is not None and not hasattr(process, linkLabel):
            setattr(process, linkLabel, l1tPh3SmartPixelsSeedLinkTable.clone(
                src=cms.InputTag(_spLabels(variant)[side], "Level1TTTracks"),
                seedTracks=cms.InputTag(_OT_TRACKS[side], "Level1TTTracks"),
                trackTableName=name,
                seedTableName=_OT_TABLES[side][1]))
            _addToNano(process, "p3L1SmartPixelsLinkTask", [linkLabel])


def _configureModeTables(process, variant, sched):
    """The digiRefit variant's own tables: per-crossing refit hits (trackIdx -> SmartPixels
    track, selClusterIdx -> L1TSmartPixelsCluster when those tables are present) and the spix*
    columns on the SmartPixels track tables."""
    if variant is None or variant[0] != "digiRefit":
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


def _configureStubLinks(process, variant, sched):
    """Link-only track -> stub tables (<track table>Stub: trackIdx, otStubIdx -> L1TOTStub)
    for the OT and the SmartPixels track tables, and SmartPixels track -> cluster tables
    (<track table>Cluster: trackIdx, clusterIdx -> L1TSmartPixelsCluster; for digiRefit the
    hits accepted into the refit), each when both of its tables are in the nano."""
    from DPGAnalysis.Phase3SmartPixelsNanoAOD.l1tPh3SmartPixelsNano_cff import (
        addPh3L1SmartPixelsStubTables, l1tPh3TrackStubLinkTable, l1tPh3SmartPixelsTrackClusterLinkTable)
    if _OT_STUB_TABLE in sched:
        doPrompt, doExtended = (_OT_TABLES[s][0] in sched for s in ("prompt", "ext"))
        if doPrompt or doExtended:
            addPh3L1SmartPixelsStubTables(process, doPrompt=doPrompt, doExtended=doExtended)
    if variant is None:
        return
    for side, (tableLabel, tableName) in _SP_TABLES.items():
        if tableLabel not in sched:
            continue
        ext = "Ext" if side == "ext" else ""
        tracks = cms.InputTag(_spLabels(variant)[side], "Level1TTTracks")
        if _OT_STUB_TABLE in sched:
            label = f"l1tPh3{ext}SmartPixelsTrackStubLinkTable"
            if not hasattr(process, label):
                setattr(process, label, l1tPh3TrackStubLinkTable.clone(
                    tracks=tracks, linkTableName=f"{tableName}Stub", trackTableName=tableName))
                _addToNano(process, "p3L1SmartPixelsLinkTask", [label])
        if _CLUSTER_TABLE in sched:
            label = f"l1tPh3{ext}SmartPixelsTrackClusterLinkTable"
            if not hasattr(process, label):
                setattr(process, label, l1tPh3SmartPixelsTrackClusterLinkTable.clone(
                    tracks=tracks, linkTableName=f"{tableName}Cluster", trackTableName=tableName,
                    sidecar=tracks if variant[0] == "digiRefit" else cms.InputTag("")))
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


def _configureTrackingParticles(process, st, variant, sched):
    """L1TTP selection of the TrackingParticle flavour. The linked selections take their
    TrackingParticles from the truth associations of the objects the nano writes: OT stubs
    (L1TOTStub), OT tracks, SmartPixels tracks (the associator of the wired variant) and IT
    clusters (their dominant TrackingParticle)."""
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
        stubs.append(getattr(process, _OT_STUB_TABLE).stubTruth)
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
