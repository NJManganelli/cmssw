import FWCore.ParameterSet.Config as cms

# Phase-3 SmartPixels L1 track nano tables.
#
# SmartPixels tracking is a beyond-Phase-2 (Phase-3) study, so it lives in its own
# package (companion to L1Trigger/Phase3SmartPixels) and is enabled independently
# via the autoNANO SmartPix flavours (see PhysicsTools/NanoAOD/autoNANO.py), the
# smartPixelsCoexist/smartPixelsCoopt customisations, or by calling
# addPh3L1SmartPixelsTracks directly. Requires the SmartPixels track collections
# from L1Trigger/Phase3SmartPixels.
from DPGAnalysis.Phase2L1TNanoAOD.l1tPh2TrkNanotables_cff import l1tTracksTable

l1tPh3SmartPixelsTracksTable = l1tTracksTable.clone(
    src = cms.InputTag("l1tSmartPixelsTrackProducer", "Level1TTTracks"),
    name = cms.string("L1TSmartPixelsTrack"),
    doc = cms.string("L1T Phase-3 SmartPixels Tracks"),
)

l1tPh3ExtSmartPixelsTracksTable = l1tTracksTable.clone(
    src = cms.InputTag("l1tSmartPixelsTrackProducerExtended", "Level1TTTracks"),
    name = cms.string("L1TSmartPixelsExtTrack"),
    doc = cms.string("L1T Phase-3 SmartPixels Extended Tracks"),
)

# l1TrackIdx on a SmartPixels track table: the row of the OT track (L1TTrack /
# L1TExtTrack) it was seeded from, read from the producer's explicit seed links.
l1tPh3SmartPixelsSeedLinkTable = cms.EDProducer(
    "L1SmartPixelsSeedLinkTableProducer",
    src = cms.InputTag("l1tSmartPixelsTrackProducer", "Level1TTTracks"),
    seedTracks = cms.InputTag("l1tTTTracksFromTrackletEmulation", "Level1TTTracks"),
    trackTableName = cms.string("L1TSmartPixelsTrack"),
    seedTableName = cms.string("L1TTrack"),
)

def addPh3L1SmartPixelsTracks(process,
                              srcLabel="l1tSmartPixelsTrackProducer",
                              srcExtendedLabel="l1tSmartPixelsTrackProducerExtended",
                              tableSuffix=""):
    """Add one pair (prompt+extended) of SmartPixels track tables to the L1Nano task.

    srcLabel/srcExtendedLabel select the SmartPixels collections (any variant
    produced by L1Trigger/Phase3SmartPixels, e.g.
    l1tSmartPixelsTrackProducerWdigiRefitAAAA); tableSuffix
    disambiguates the table/branch names when several variants coexist in one
    job (WF1). The defaults (the autoNANO SmartPix flavors) create the canonical
    unsuffixed pair, whose input only smartPixelsCoopt provides (it repoints the
    pair at its variant); without coopt the job stops with ProductNotFound.
    Each table carries l1TrackIdx -> L1TTrack / L1TExtTrack (the seed track row).
    """
    nameSuffix = tableSuffix[:1].upper() + tableSuffix[1:] if tableSuffix else ""
    # Module labels must END in "Table": the NANOAOD event content only keeps
    # "nanoaodFlatTable_*Table_*_*", and the output module's consumes (which
    # trigger the unscheduled table producers) derive from that pattern.
    promptTableName = f"l1tPh3SmartPixelsTracks{nameSuffix}Table"
    extendedTableName = f"l1tPh3ExtSmartPixelsTracks{nameSuffix}Table"

    setattr(process, promptTableName, l1tPh3SmartPixelsTracksTable.clone(
        src = cms.InputTag(srcLabel, "Level1TTTracks"),
        name = cms.string(f"L1TSmartPixelsTrack{nameSuffix}"),
        doc = cms.string(f"L1T Phase-3 SmartPixels Tracks ({srcLabel})"),
    ))
    setattr(process, extendedTableName, l1tPh3ExtSmartPixelsTracksTable.clone(
        src = cms.InputTag(srcExtendedLabel, "Level1TTTracks"),
        name = cms.string(f"L1TSmartPixelsExtTrack{nameSuffix}"),
        doc = cms.string(f"L1T Phase-3 SmartPixels Extended Tracks ({srcExtendedLabel})"),
    ))

    promptLinkName = f"l1tPh3SmartPixelsSeedLink{nameSuffix}Table"
    extendedLinkName = f"l1tPh3ExtSmartPixelsSeedLink{nameSuffix}Table"
    setattr(process, promptLinkName, l1tPh3SmartPixelsSeedLinkTable.clone(
        src = cms.InputTag(srcLabel, "Level1TTTracks"),
        trackTableName = cms.string(f"L1TSmartPixelsTrack{nameSuffix}"),
    ))
    setattr(process, extendedLinkName, l1tPh3SmartPixelsSeedLinkTable.clone(
        src = cms.InputTag(srcExtendedLabel, "Level1TTTracks"),
        seedTracks = cms.InputTag("l1tTTTracksFromExtendedTrackletEmulation", "Level1TTTracks"),
        trackTableName = cms.string(f"L1TSmartPixelsExtTrack{nameSuffix}"),
        seedTableName = cms.string("L1TExtTrack"),
    ))

    task = cms.Task(*[getattr(process, t) for t in (promptTableName, extendedTableName,
                                                   promptLinkName, extendedLinkName)])
    setattr(process, f"p3L1SmartPixelsTracksTask{nameSuffix}", task)
    process.l1tPh2NanoTask.add(task)
    return process


# ---------------------------------------------------------------------------
# digiRefit sidecar tables (spec RefitSidecarSpec.md §4.4)
# ---------------------------------------------------------------------------
# The digiRefit producer emits a smartpixels::SmartPixelsRefitSidecar 1:1 row-
# synced with its refit TTTrack collection, under the SAME instance label
# ("Level1TTTracks"). The L1SmartPixelsRefitTableProducer turns (tracks+sidecar)
# into a per-crossing LINK table ("hit", one row per crossing, trackIdx links to
# the variant track table) and an EXTENSION table ("trk", extension=True, same
# name+length as the variant track table). See the plugin doc for the schemas.
l1tPh3SmartPixelsRefitTable = cms.EDProducer(
    "L1SmartPixelsRefitTableProducer",
    tracks = cms.InputTag("l1tSmartPixelsTrackProducer", "Level1TTTracks"),
    sidecar = cms.InputTag("l1tSmartPixelsTrackProducer", "Level1TTTracks"),
    trackTableName = cms.string("L1TSmartPixelsTrack"),
    hitTableName = cms.string("L1TSmartPixelsRefitHit"),
)


def addPh3L1SmartPixelsRefitTables(process,
                                   srcLabel="l1tSmartPixelsTrackProducer",
                                   srcExtendedLabel="l1tSmartPixelsTrackProducerExtended",
                                   tableSuffix=""):
    """Add the digiRefit sidecar tables (per-hit link + track extension) for one
    variant (prompt+extended). ONLY valid for digiRefit variants -- the sidecar
    product only exists in digiRefit mode. Non-digiRefit variants get tracks-only
    tables via addPh3L1SmartPixelsTracks (no sidecar is produced, so consuming it
    would throw ProductNotFound). The extension table name MUST equal the variant
    track table name (from addPh3L1SmartPixelsTracks) so nano merges the columns.
    """
    nameSuffix = tableSuffix[:1].upper() + tableSuffix[1:] if tableSuffix else ""
    # Module labels must END in "Table" (keep-pattern gotcha): NANOAOD event
    # content keeps only nanoaodFlatTable_*Table_*_* and the output module's
    # consumes (which trigger the unscheduled producers) derive from that.
    promptRefitLabel = f"l1tPh3SmartPixelsRefit{nameSuffix}Table"
    extendedRefitLabel = f"l1tPh3ExtSmartPixelsRefit{nameSuffix}Table"

    setattr(process, promptRefitLabel, l1tPh3SmartPixelsRefitTable.clone(
        tracks = cms.InputTag(srcLabel, "Level1TTTracks"),
        sidecar = cms.InputTag(srcLabel, "Level1TTTracks"),
        trackTableName = cms.string(f"L1TSmartPixelsTrack{nameSuffix}"),
        hitTableName = cms.string(f"L1TSmartPixelsRefitHit{nameSuffix}"),
    ))
    setattr(process, extendedRefitLabel, l1tPh3SmartPixelsRefitTable.clone(
        tracks = cms.InputTag(srcExtendedLabel, "Level1TTTracks"),
        sidecar = cms.InputTag(srcExtendedLabel, "Level1TTTracks"),
        trackTableName = cms.string(f"L1TSmartPixelsExtTrack{nameSuffix}"),
        hitTableName = cms.string(f"L1TSmartPixelsExtRefitHit{nameSuffix}"),
    ))

    task = cms.Task(getattr(process, promptRefitLabel), getattr(process, extendedRefitLabel))
    setattr(process, f"p3L1SmartPixelsRefitTask{nameSuffix}", task)
    process.l1tPh2NanoTask.add(task)
    return process


# ---------------------------------------------------------------------------
# Reference-track OT stub position tables (task Part B; viz true-anchor + refit
# research). The L1SmartPixelsStubPosTableProducer walks getStubRefs() on a
# reference TTTrack collection and emits ONE per-stub LINK table (one row per stub
# across all tracks) carrying the stub's real global position (x,y,z,r,phi), its
# tracker layer/disk id, isBarrel, detId, and FE bend. trackIdx links back to the
# reference track table (the L1SC4NGJetCands / refit-hit link-table pattern).
#
# KEYED TO THE REFERENCE COLLECTIONS, NOT THE REFIT VARIANTS (task Part B.2):
#   prompt  -> l1tTTTracksFromTrackletEmulation:Level1TTTracks       (nano "L1TTrack")
#   extended-> l1tTTTracksFromExtendedTrackletEmulation:Level1TTTracks (nano "L1TExtTrack")
# The refit variants reuse the SAME stub refs (setStubRefs) so the stub geometry is
# identical -- but the SEED's stubs are the physics constraint being visualized, and
# the reference tables are the untouched anchor that trackIdx and the L1TrackTruth
# table both resolve against. Emitting once against the reference collections (rather
# than once per variant) is therefore both correct AND avoids redundant identical
# copies. Under trackInputMode=rebuildTracksFromStubs these reference labels are re-run in-job
# from the file's stubs, so the stub geometry is the fresh new-layout one.
l1tPh3SmartPixelsStubTable = cms.EDProducer(
    "L1SmartPixelsStubPosTableProducer",
    tracks = cms.InputTag("l1tTTTracksFromTrackletEmulation", "Level1TTTracks"),
    # the collection L1TOTStub is built from: otStubIdx is a row there (product ID asserted)
    stubs = cms.InputTag("TTStubsFromPhase2TrackerDigis", "StubAccepted"),
    stubTableName = cms.string("L1TTrackStub"),
    trackTableName = cms.string("L1TTrack"),
)


def addPh3L1SmartPixelsStubTables(process,
                                  promptSrc=("l1tTTTracksFromTrackletEmulation", "Level1TTTracks"),
                                  extendedSrc=("l1tTTTracksFromExtendedTrackletEmulation", "Level1TTTracks"),
                                  promptTrackTable="L1TTrack",
                                  extendedTrackTable="L1TExtTrack",
                                  doExtended=True):
    """Add the reference-track OT stub-position tables (prompt + extended).

    One per-stub LINK table per reference collection, so the viz can render the
    true OT anchors of whichever seed (prompt or extended) it visualizes. The
    stub table's trackIdx links into promptTrackTable / extendedTrackTable (the
    reference L1TTrack / L1TExtTrack nano tables), which the L1TrackTruth table
    also aligns to -- so stubs, kinematics and truth all share one track index.

    Idempotent per table name: skips if the label already exists (a job that
    calls this twice, e.g. once per SmartPix flavor, adds each table once).
    """
    # Module labels must END in "Table" (keep-pattern gotcha): NANOAOD event
    # content keeps only nanoaodFlatTable_*Table_*_* and the output module's
    # consumes (which trigger the unscheduled producers) derive from that.
    promptLabel = "l1tPh3SmartPixelsStubTable"
    extendedLabel = "l1tPh3ExtSmartPixelsStubTable"
    newModules = []
    if not hasattr(process, promptLabel):
        setattr(process, promptLabel, l1tPh3SmartPixelsStubTable.clone(
            tracks = cms.InputTag(*promptSrc),
            stubTableName = cms.string(f"{promptTrackTable}Stub"),
            trackTableName = cms.string(promptTrackTable),
        ))
        newModules.append(promptLabel)
    if doExtended and not hasattr(process, extendedLabel):
        setattr(process, extendedLabel, l1tPh3SmartPixelsStubTable.clone(
            tracks = cms.InputTag(*extendedSrc),
            stubTableName = cms.string(f"{extendedTrackTable}Stub"),
            trackTableName = cms.string(extendedTrackTable),
        ))
        newModules.append(extendedLabel)
    if newModules:
        task = cms.Task(*[getattr(process, m) for m in newModules])
        process.p3L1SmartPixelsStubTask = task
        process.l1tPh2NanoTask.add(task)
    assertTrackTablesUncut(process)
    return process


# Track tables whose row i MUST be collection element i: the truth / covariance /
# refit extensions, the stub tables' trackIdx and the refit hit tables' trackIdx
# all index by collection position, and a variant track row matches the L1TTrack
# row only because both are uncut.
_ROW_ALIGNED_TRACK_TABLES = ("L1TTrack", "L1TExtTrack", "L1TSmartPixelsTrack", "L1TSmartPixelsExtTrack")


def assertTrackTablesUncut(process):
    """Raise if any main L1 track table (reference or SmartPixels variant) has a cut.

    A cut silently shifts rows against every extension and link table keyed by
    collection position, so it must be applied in analysis, not in the table.
    """
    def _val(mod, attr, default=None):
        # declared-but-unset (_RequiredParameter) parameters carry no value
        try:
            return getattr(mod, attr).value()
        except AttributeError:
            return default

    for label in process.producers_():
        mod = getattr(process, label)
        name, cut = _val(mod, "name"), _val(mod, "cut")
        if not isinstance(name, str) or not isinstance(cut, str) or _val(mod, "extension", False):
            continue
        if not name.startswith(_ROW_ALIGNED_TRACK_TABLES):
            continue
        if cut.strip():
            raise RuntimeError(
                f"track table '{label}' ({name}) has cut={cut!r}: its rows would no "
                "longer be collection positions, breaking trackIdx/otStubIdx links and the "
                "truth/refit extensions. Cut in analysis instead.")
    return process


# ---------------------------------------------------------------------------
# Coopt table wiring: the L1 emulation consumes the SmartPixels tracks
# ---------------------------------------------------------------------------
def retargetTablesToSmartPixelsTracks(process, srcLabel, srcExtendedLabel):
    """Coopt (smartPixelsCoopt): the L1 emulation was built from the SmartPixels
    tracks, so every table that indexes the tracks an L1 object was built from must
    index the canonical SmartPixels tables, and say so.

      * candidate links (L1JetCandLinkTableProducer: L1PuppiCand / L1ExtPuppiCand
        l1TrackIdx) -> L1TSmartPixelsTrack / L1TSmartPixelsExtTrack rows;
      * candidate track truth (L1PFCandTrackTruthTableProducer) -> the SmartPixels
        tracks with the TTTrackAssociator OF THAT COLLECTION (never the tracklet
        associator: its maps are keyed by tracklet-track Ptrs);
      * L1DispVertex firstIndexTrk/secondIndexTrk: DisplacedVertexProducer ran on the
        extended SmartPixels tracks, so its indices are L1TSmartPixelsExtTrack rows.

    L1TTrack/L1TExtTrack stay the pure OT tables. Their GTT hardware extension
    (gttTracksTable/gttExtTracksTable: hwPt/hwEta after GTT input conversion) is fed
    by a GTT input conversion of the OT tracks themselves, because the L1 GTT input
    (l1tGTTInputProducer[Extended]) now converts the SmartPixels tracks.
    Explicit table wiring, not the injector: the injector never touches tables.
    """
    from L1Trigger.Phase3SmartPixels.customizeSmartPixels_cff import smartPixelsTrackAssociator
    sp = {"l1tTTTracksFromTrackletEmulation": (srcLabel, "L1TSmartPixelsTrack", False),
          "l1tTTTracksFromExtendedTrackletEmulation": (srcExtendedLabel, "L1TSmartPixelsExtTrack", True)}
    scheduled = process.l1tPh2NanoTask.moduleNames()
    for label, mod in process.producers_().items():
        if label not in scheduled or mod.type_() not in ("L1JetCandLinkTableProducer",
                                                         "L1PFCandTrackTruthTableProducer"):
            continue
        target = sp.get(mod.tracks.getModuleLabel())
        if target is None:
            raise RuntimeError(f"retargetTablesToSmartPixelsTracks: {label}.tracks = "
                               f"{mod.tracks.configValue()} is not a tracklet collection")
        spLabel, spTable, extended = target
        mod.tracks = cms.InputTag(spLabel, "Level1TTTracks")
        mod.trackTableName = cms.string(spTable)
        if mod.type_() == "L1PFCandTrackTruthTableProducer":
            mod.truthAssociation = cms.InputTag(smartPixelsTrackAssociator(process, spLabel, extended),
                                                "Level1TTTracks")
        print(f"SmartPixels coopt tables: {label} tracks -> {spLabel} ({spTable})")

    if hasattr(process, "l1tDisplacedVertexTable"):
        for col, which in (("firstIndexTrk", "higher"), ("secondIndexTrk", "lower")):
            getattr(process.l1tDisplacedVertexTable.variables, col).doc = cms.string(
                f"index of the {which}-pt track in the L1TSmartPixelsExtTrack table")

    from L1Trigger.L1TTrackMatch.l1tGTTInputProducer_cfi import (
        l1tGTTInputProducer, l1tGTTInputProducerExtended)
    gtt = []
    for table, conv, cfi, otLabel, instance in (
            ("gttTracksTable", "l1tGTTInputProducerOT", l1tGTTInputProducer,
             "l1tTTTracksFromTrackletEmulation", "Level1TTTracksConverted"),
            ("gttExtTracksTable", "l1tGTTInputProducerExtendedOT", l1tGTTInputProducerExtended,
             "l1tTTTracksFromExtendedTrackletEmulation", "Level1TTTracksExtendedConverted")):
        if not hasattr(process, table):
            continue
        setattr(process, conv, cfi.clone(l1TracksInputTag=cms.InputTag(otLabel, "Level1TTTracks")))
        getattr(process, table).src = cms.InputTag(conv, instance)
        gtt.append(getattr(process, conv))
    if gtt:
        process.smartPixelsOTGTTConversionTask = cms.Task(*gtt)
        process.l1tPh2NanoTask.add(process.smartPixelsOTGTTConversionTask)
    return process


def dropAbsentMenuTables(process, moduleLabels):
    """Remove nano table modules (by label). NOT for canonical tiers (a tier defines
    its content; produce a missing input instead): kept for runUnifiedNano_cfg's
    MINIAOD-only reco tables on DIGI-RAW inputs. Silent no-op for absent labels."""
    for label in moduleLabels:
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
        print(f"SmartPixels: removed nano table '{label}'")
    return process


def useGenParticlesFromFile(process, genSrc="genParticles"):
    """Adapt the withGen gen chain to a GEN-SIM(-DIGI-RAW) input that has
    `genParticles` but NOT the MINIAOD `prunedGenParticles` by repointing
    finalGenParticles and DROPPING the gen tables that need packedGenParticles or
    slimmed gen jets. Not for canonical tiers (it drops tables): use
    produceMiniAODGenInputs. Kept for runUnifiedNano_cfg.
    """
    if hasattr(process, "finalGenParticles"):
        process.finalGenParticles.src = cms.InputTag(genSrc)
    if hasattr(process, "genParticleTable") and hasattr(process.genParticleTable, "externalVariables"):
        if hasattr(process.genParticleTable.externalVariables, "iso"):
            del process.genParticleTable.externalVariables.iso
    process = dropAbsentMenuTables(process, ("genIso",))
    return dropAbsentMenuTables(process, ("genJetTable", "genJetFlavourTable", "patJetPartonsNano",
                                          "genJetAK8Table", "genJetAK8FlavourAssociation",
                                          "genJetAK8FlavourTable", "genJetFlavourAssociation"))


def produceMiniAODGenInputs(process):
    """withGen tiers on a GEN-SIM(-DIGI-RAW) input: produce the MINIAOD gen products
    the stock gen tables read (prunedGenParticles, packedGenParticles, slimmedGenJets,
    slimmedGenJetsFlavourInfos) from the persisted genParticles / ak4GenJetsNoNu with
    the stock PAT slimming modules, on demand. Only for inputs WITHOUT those products:
    on a MINIAOD-carrying input they would shadow the persisted ones.
    """
    from PhysicsTools.PatAlgos.slimming.genParticles_cff import (
        prunedGenParticlesWithStatusOne, prunedGenParticles, packedGenParticles)
    from PhysicsTools.PatAlgos.slimming.slimmedGenJets_cfi import slimmedGenJets
    from PhysicsTools.PatAlgos.slimming.slimmedGenJetsFlavourInfos_cfi import slimmedGenJetsFlavourInfos
    from PhysicsTools.JetMCAlgos.HadronAndPartonSelector_cfi import selectedHadronsAndPartons
    from PhysicsTools.JetMCAlgos.AK4GenJetFlavourInfos_cfi import ak4GenJetFlavourInfos
    modules = (("prunedGenParticlesWithStatusOne", prunedGenParticlesWithStatusOne.clone()),
               ("prunedGenParticles", prunedGenParticles.clone()),
               ("packedGenParticles", packedGenParticles.clone()),
               ("slimmedGenJets", slimmedGenJets.clone()),
               ("selectedHadronsAndPartonsForGenJetsFlavourInfos",
                selectedHadronsAndPartons.clone(particles="genParticles")),
               ("ak4GenJetFlavourInfos", ak4GenJetFlavourInfos.clone(jets="ak4GenJetsNoNu")),
               ("slimmedGenJetsFlavourInfos", slimmedGenJetsFlavourInfos.clone()))
    for label, module in modules:
        if hasattr(process, label):
            raise RuntimeError(f"produceMiniAODGenInputs: '{label}' is already defined in this job")
        setattr(process, label, module)
    process.smartPixelsMiniAODGenTask = cms.Task(*[getattr(process, l) for l, _ in modules])
    process.l1tPh2NanoTask.add(process.smartPixelsMiniAODGenTask)
    return process


# ---------------------------------------------------------------------------
# Untruncated IT cluster table  (tier: L1PFTrkNanoSmartPixClusters)
# ---------------------------------------------------------------------------
# EVERY IT pixel cluster, with no window, no truncation and no truth-pT filter:
# all TBPX clusters first, then the disc clusters in discScope (default sp =
# ID1-ID4), one table with a surface column (1-4 = IL1-IL4, 11-14 = ID1-ID4).
#
# ~26.5k barrel rows/event at PU200, measured at 88.6 stored bits/cluster = 0.29 MB/event
# for the barrel; ID1-ID4 add ~46k rows/event.
# Belongs ONLY to the Clusters tier. The size is the point: at PU200 the crossed module
# carries 40.6/32.7/31.6/19.7 clusters on L1-L4 (p95 up to 72), the static window
# admits about 2 of them, and maxHitsPerWindow truncates at 8. Every question
# about window sizing or candidate multiplicity is a question about the clusters
# the current window DISCARDED, so any filtered table is unable to answer it --
# and a truth-pT filter is actively misleading, since it removes the soft
# clusters that do the confusing (a >2 GeV cut keeps 2.2% of clusters).
l1tPh3SmartPixelsClusterTable = cms.EDProducer(
    "L1SmartPixelsClusterTableProducer",
    smartPixelsRecHits = cms.InputTag("spixSmartPixelsRecHits"),
    tableName = cms.string("L1TSmartPixelsCluster"),
    discScope = cms.string("sp"),  # MUST equal spixSmartPixelsRecHits.discScope (checked at run time)
    maxLayer = cms.uint32(4),
    doTruth = cms.bool(True),
)


# The COMPLETE OT stub collection, part of the Clusters tier (not a separate
# variant). ~15k rows/event at PU200 expected, roughly doubling the tier payload.
#
# It belongs here rather than in a tier of its own because the questions this tier
# exists to answer are comparative: an IT seeding design is only interesting
# relative to what the OT tracklet already does, and that comparison needs the OT
# INPUT multiplicity. The existing stub tables (L1TTrackStub, L1TExtTrackStub) are
# on-track only -- 1004 rows/event at PU200 against ~190 tracks -- so they cannot
# supply it, and the persisted edmNew::DetSetVector is unreadable by uproot. Any
# study wanting both halves would otherwise have to splice two data tiers.
l1tPh3SmartPixelsAllStubTable = cms.EDProducer(
    "L1SmartPixelsAllStubTableProducer",
    stubs = cms.InputTag("TTStubsFromPhase2TrackerDigis", "StubAccepted"),
    stubTruth = cms.InputTag("TTStubAssociatorFromPixelDigis", "StubAccepted"),
    tableName = cms.string("L1TOTStub"),
    doTruth = cms.bool(True),
)


def addPh3L1SmartPixelsClusters(process, recHitLabel="spixSmartPixelsRecHits", doTruth=True):
    """Add the untruncated IT cluster table + the complete OT stub table (Clusters tier).

    Requires the SmartPixels-owned cluster -> rec-hit chain (spixPixelClusters ->
    spixPixelRecHits). That chain is created by
    customizeSmartPixels_cff._addPixelRecHitChain whenever a digiRefit variant is
    requested; this function does NOT create it, because creating a second,
    differently-configured clusterizer here is exactly how the cluster table would
    silently stop describing the clusters the refit actually saw.
    """
    # The producer need NOT exist yet: autoNANO --customise functions run BEFORE
    # --customise_commands, so smartPixelsCoexist -- which builds
    # spixSmartPixelsRecHits and owns the PixelAV payload -- has not run at this
    # point. An InputTag is just a label and resolves at run time. This table
    # deliberately does NOT build its own rec hits: it must describe the SAME ones,
    # with the SAME angles, that the refit consumed. smartPixelsCoexist asserts
    # that pairing at the end, where the final process is visible.
    # Module label must END in 'Table' (NANOAOD keeps nanoaodFlatTable_*Table_*_*).
    process.l1tPh3SmartPixelsClusterTable = l1tPh3SmartPixelsClusterTable.clone(
        smartPixelsRecHits = cms.InputTag(recHitLabel),
        doTruth = cms.bool(doTruth),
    )
    # Follow the rec hits' disc scope when they already exist; if they are built
    # later, ensureSmartPixelsRecHits sets it on this table.
    if hasattr(process, recHitLabel) and hasattr(getattr(process, recHitLabel), "discScope"):
        process.l1tPh3SmartPixelsClusterTable.discScope = getattr(process, recHitLabel).discScope
    # The OT stub table reads TTStubsFromPhase2TrackerDigis:StubAccepted straight
    # from the input file with NO process name, so it resolves whether the stubs
    # were persisted by the RelVal or rebuilt in-job -- the same reason the
    # associator labels are left process-less in customizeSmartPixels_cff.
    # doTruth follows the cluster table: an OT side without truth could only be
    # compared on combinatorics while the IT side has efficiency and fake rate,
    # and that asymmetry silently biases any comparison drawn from it.
    process.l1tPh3SmartPixelsAllStubTable = l1tPh3SmartPixelsAllStubTable.clone(
        doTruth = cms.bool(doTruth),
    )
    task = cms.Task(process.l1tPh3SmartPixelsClusterTable,
                    process.l1tPh3SmartPixelsAllStubTable)
    process.p3L1SmartPixelsClusterTask = task
    process.l1tPh2NanoTask.add(task)
    return process


def addPh3L1SmartPixelsReco(process):
    """RESERVED (tier: L1PFTrkNanoSmartPixReco): full offline RECO quantities for
    physics studies, plus enough refit information to interpret them, and NO
    cluster table.

    The name is reserved so the tier can be referred to and scheduled for, but the
    content is not specified yet. It raises rather than silently producing a tier
    identical to L1PFTrkNanoSmartPix, which would be indistinguishable downstream
    from the real thing and would quietly invalidate any study that used it.
    """
    raise NotImplementedError(
        "L1PFTrkNanoSmartPixReco is RESERVED: the tier name exists but its content "
        "(offline RECO quantities + refit-interpretation columns, no clusters) is not "
        "defined yet. Use L1PFTrkNanoSmartPix or L1PFTrkNanoSmartPixClusters.")
