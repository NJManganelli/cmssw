import FWCore.ParameterSet.Config as cms

# Phase-3 SmartPixels L1 track nano tables.
#
# SmartPixels tracking is a beyond-Phase-2 (Phase-3) study, so it lives in its own
# package (companion to L1Trigger/Phase3SmartPixels). These are the table definitions
# and the coexist helpers; the canonical tiers are the add-on flavours of
# l1tPh3SmartPixelsNanoFlavours_cff (see PhysicsTools/NanoAOD/autoNANO.py).
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
# L1TExtTrack) it was built from, read from the producer's explicit OT-track links.
l1tPh3SmartPixelsOTTrackLinkTable = cms.EDProducer(
    "L1SmartPixelsOTTrackLinkTableProducer",
    src = cms.InputTag("l1tSmartPixelsTrackProducer", "Level1TTTracks"),
    otTracks = cms.InputTag("l1tTTTracksFromTrackletEmulation", "Level1TTTracks"),
    trackTableName = cms.string("L1TSmartPixelsTrack"),
    otTrackTableName = cms.string("L1TTrack"),
)

def addPh3L1SmartPixelsTracks(process,
                              srcLabel="l1tSmartPixelsTrackProducer",
                              srcExtendedLabel="l1tSmartPixelsTrackProducerExtended",
                              tableSuffix=""):
    """Add one pair (prompt+extended) of SmartPixels track tables to the L1Nano task.

    srcLabel/srcExtendedLabel select the SmartPixels collections (any variant
    produced by L1Trigger/Phase3SmartPixels, e.g.
    l1tSmartPixelsTrackProducerWseedlessRefitAAAA); tableSuffix
    disambiguates the table/branch names when several variants coexist in one
    job (WF1, smartPixelsCoexist). The canonical unsuffixed pair belongs to the @SPTrk
    flavour. Each table carries l1TrackIdx -> L1TTrack / L1TExtTrack (the row of the OT track it was built from).
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

    promptLinkName = f"l1tPh3SmartPixelsOTTrackLink{nameSuffix}Table"
    extendedLinkName = f"l1tPh3ExtSmartPixelsOTTrackLink{nameSuffix}Table"
    setattr(process, promptLinkName, l1tPh3SmartPixelsOTTrackLinkTable.clone(
        src = cms.InputTag(srcLabel, "Level1TTTracks"),
        trackTableName = cms.string(f"L1TSmartPixelsTrack{nameSuffix}"),
    ))
    setattr(process, extendedLinkName, l1tPh3SmartPixelsOTTrackLinkTable.clone(
        src = cms.InputTag(srcExtendedLabel, "Level1TTTracks"),
        otTracks = cms.InputTag("l1tTTTracksFromExtendedTrackletEmulation", "Level1TTTracks"),
        trackTableName = cms.string(f"L1TSmartPixelsExtTrack{nameSuffix}"),
        otTrackTableName = cms.string("L1TExtTrack"),
    ))

    task = cms.Task(*[getattr(process, t) for t in (promptTableName, extendedTableName,
                                                   promptLinkName, extendedLinkName)])
    setattr(process, f"p3L1SmartPixelsTracksTask{nameSuffix}", task)
    process.l1tPh2NanoTask.add(task)
    return process


# ---------------------------------------------------------------------------
# seedlessRefit sidecar tables (spec RefitSidecarSpec.md §4.4)
# ---------------------------------------------------------------------------
# The seedlessRefit producer emits a smartpixels::SmartPixelsRefitSidecar 1:1 row-
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
    clusterTableName = cms.string("L1TSmartPixelsCluster"),
    hitTableName = cms.string("L1TSmartPixelsRefitHit"),
)


def addPh3L1SmartPixelsRefitTables(process,
                                   srcLabel="l1tSmartPixelsTrackProducer",
                                   srcExtendedLabel="l1tSmartPixelsTrackProducerExtended",
                                   tableSuffix=""):
    """Add the seedlessRefit sidecar tables (per-hit link + track extension) for one
    variant (prompt+extended). ONLY valid for seedlessRefit variants -- the sidecar
    product only exists in seedlessRefit mode. Non-seedlessRefit variants get tracks-only
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
# Track -> OT stub link tables. Each is LINK-ONLY: one row per (track, stub) pair with
# trackIdx (row in the track table) and otStubIdx (row in L1TOTStub, the complete OT
# stub table, which carries the stub position, layer, bend and truth). The ref is
# asserted to point into the stub collection L1TOTStub is built from.
# The SmartPixels tracks carry the stub refs of the OT track they were built from, so
# the same producer serves both.
# ---------------------------------------------------------------------------
l1tPh3TrackStubLinkTable = cms.EDProducer(
    "L1SmartPixelsTrackStubLinkTableProducer",
    tracks = cms.InputTag("l1tTTTracksFromTrackletEmulation", "Level1TTTracks"),
    stubs = cms.InputTag("TTStubsFromPhase2TrackerDigis", "StubAccepted"),
    linkTableName = cms.string("L1TTrackStub"),
    trackTableName = cms.string("L1TTrack"),
    stubTableName = cms.string("L1TOTStub"),
)

# SmartPixels track -> IT cluster link table (trackIdx, clusterIdx -> L1TSmartPixelsCluster):
# the IT hits the variant used, from the seedlessRefit sidecar; no rows without a sidecar.
l1tPh3SmartPixelsTrackClusterLinkTable = cms.EDProducer(
    "L1SmartPixelsTrackClusterLinkTableProducer",
    tracks = cms.InputTag("l1tSmartPixelsTrackProducer", "Level1TTTracks"),
    sidecar = cms.InputTag(""),
    linkTableName = cms.string("L1TSmartPixelsTrackCluster"),
    trackTableName = cms.string("L1TSmartPixelsTrack"),
    clusterTableName = cms.string("L1TSmartPixelsCluster"),
)


def addPh3L1SmartPixelsStubTables(process,
                                  promptSrc=("l1tTTTracksFromTrackletEmulation", "Level1TTTracks"),
                                  extendedSrc=("l1tTTTracksFromExtendedTrackletEmulation", "Level1TTTracks"),
                                  promptTrackTable="L1TTrack",
                                  extendedTrackTable="L1TExtTrack",
                                  doPrompt=True,
                                  doExtended=True):
    """Add L1TOTStub and the link-only OT track -> stub tables (<track table>Stub:
    trackIdx, otStubIdx) of the reference collections (prompt + extended).

    Idempotent per module label: a job that calls this twice adds each table once.
    """
    # Module labels must END in "Table" (keep-pattern gotcha): NANOAOD event
    # content keeps only nanoaodFlatTable_*Table_*_* and the output module's
    # consumes (which trigger the unscheduled producers) derive from that.
    labels = ["l1tPh3SmartPixelsAllStubTable"]
    if not hasattr(process, labels[0]):
        process.l1tPh3SmartPixelsAllStubTable = l1tPh3SmartPixelsAllStubTable.clone()
    for do, label, src, trackTable in ((doPrompt, "l1tPh3TrackStubLinkTable", promptSrc, promptTrackTable),
                                       (doExtended, "l1tPh3ExtTrackStubLinkTable", extendedSrc, extendedTrackTable)):
        if not do:
            continue
        if not hasattr(process, label):
            setattr(process, label, l1tPh3TrackStubLinkTable.clone(
                tracks = cms.InputTag(*src),
                linkTableName = cms.string(f"{trackTable}Stub"),
                trackTableName = cms.string(trackTable),
            ))
        labels.append(label)
    if not hasattr(process, "p3L1SmartPixelsStubTask"):
        process.p3L1SmartPixelsStubTask = cms.Task()
        process.l1tPh2NanoTask.add(process.p3L1SmartPixelsStubTask)
    process.p3L1SmartPixelsStubTask.add(*[getattr(process, label) for label in labels])
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


# Output compression of every SmartPixels nano (the add-on flavours, smartPixelsCoexist,
# smartPixelsCoopt). ZSTD 5 instead of the cmsDriver default LZMA 9, whose serial output
# compression dominated the event time (PU200: 2.9 vs 3.8 s/event, file 1.18x larger,
# reads ~13x faster).
NANO_COMPRESSION_ALGORITHM = "ZSTD"
NANO_COMPRESSION_LEVEL = 5


def useSmartPixelsNanoCompression(process):
    """Set NANO_COMPRESSION_ALGORITHM/LEVEL on every NanoAODOutputModule defined so far
    (cmsDriver defines the output before the customisations run). Idempotent."""
    for out in process.outputModules_().values():
        if out.type_() == "NanoAODOutputModule":
            out.compressionAlgorithm = cms.untracked.string(NANO_COMPRESSION_ALGORITHM)
            out.compressionLevel = cms.untracked.int32(NANO_COMPRESSION_LEVEL)
    return process


_inputBranches = {}


def inputEventBranches(fileName):
    """Names of the Events branches of an EDM file (read with uproot, so no EDM dictionaries
    are needed); cached per file."""
    if fileName not in _inputBranches:
        import uproot
        try:
            with uproot.open(fileName[len("file:"):] if fileName.startswith("file:") else fileName) as f:
                _inputBranches[fileName] = f["Events"].keys(recursive=False)
        except Exception as e:
            raise RuntimeError(f"SmartPixels nano: cannot read {fileName} to inspect its products: {e}")
    return _inputBranches[fileName]


def produceMiniAODGenInputs(process, assumeNoMiniAOD=False):
    """withGen tiers on a GEN-SIM-DIGI-RAW input (no MINIAOD, no PAT in the job), run with
    --procModifiers nano_l1_hlt: rebuild the MINIAOD gen products the stock gen tables read
    from the persisted GEN-SIM products, with the stock MINIAOD slimming modules and their
    stock parameters (PhysicsTools/PatAlgos slimming + miniAOD_customizeMC), and put back the
    stock gen tables nano_l1_hlt excludes for lack of PAT inputs:
      prunedGenParticles(WithStatusOne), packedGenParticles <- genParticles
      slimmedGenJets <- ak4GenJetsNoNu;  slimmedGenJetsAK8 <- ak8GenJetsNoNu
      slimmedGenJetsFlavourInfos (selectedHadronsAndPartonsForGenJetsFlavourInfos,
        ak4GenJetFlavourInfos)
      GenJetAK8 (+ flavour): the stock genJetAK8Table / genJetAK8FlavourAssociation /
        genJetAK8FlavourTable, unchanged.
      GenMET: the stock table reads pat::MET (slimmedMETs, needs RECO); genMetTrueTable is
        the stock NANOGEN construction (customizeNanoGEN: singleton candidate table on
        genMetTrue) with the stock GenMET columns and precision. slimmedMETs carries a copy
        of genMetTrue as its genMET, so the values are the same.
      Pileup: slimmedAddPileupInfo <- addPileupInfo, and the stock puTable with
        savePUDensityVars=False (genPileupTable). pudensity/gpudensity are counted around the
        offline primary vertex (offlineSlimmedPrimaryVertices, needs RECO) and are absent;
        every other column is the stock one.
    Refused when the input carries any of these products (they would shadow the persisted
    ones; there the stock tables read them directly), when the job already defines them
    (PAT in the job), on a tier without the gen tables, or without nano_l1_hlt.
    assumeNoMiniAOD: skip the input inspection, for jobs whose input cannot be read at
    configuration time (CRAB replaces fileNames after loading the pset); the caller then
    vouches for an input without MINIAOD products.
    """
    from PhysicsTools.PatAlgos.slimming.genParticles_cff import (
        prunedGenParticlesWithStatusOne, prunedGenParticles, packedGenParticles)
    from PhysicsTools.PatAlgos.slimming.slimmedGenJets_cfi import slimmedGenJets, slimmedGenJetsAK8
    from PhysicsTools.PatAlgos.slimming.slimmedGenJetsFlavourInfos_cfi import slimmedGenJetsFlavourInfos
    from PhysicsTools.PatAlgos.slimming.slimmedAddPileupInfo_cfi import slimmedAddPileupInfo
    from PhysicsTools.JetMCAlgos.HadronAndPartonSelector_cfi import selectedHadronsAndPartonsForGenJetsFlavourInfos
    from PhysicsTools.JetMCAlgos.AK4GenJetFlavourInfos_cfi import ak4GenJetFlavourInfos
    from PhysicsTools.NanoAOD.simpleSingletonCandidateFlatTableProducer_cfi import (
        simpleSingletonCandidateFlatTableProducer)
    producers = {"prunedGenParticlesWithStatusOne": prunedGenParticlesWithStatusOne,
                 "prunedGenParticles": prunedGenParticles,
                 "packedGenParticles": packedGenParticles,
                 "slimmedGenJets": slimmedGenJets,
                 "slimmedGenJetsAK8": slimmedGenJetsAK8,
                 "selectedHadronsAndPartonsForGenJetsFlavourInfos": selectedHadronsAndPartonsForGenJetsFlavourInfos,
                 "ak4GenJetFlavourInfos": ak4GenJetFlavourInfos,
                 "slimmedGenJetsFlavourInfos": slimmedGenJetsFlavourInfos,
                 "slimmedAddPileupInfo": slimmedAddPileupInfo}
    # products whose presence marks a MINIAOD-carrying input, beyond the ones made here
    miniAODMarkers = set(producers) | {"slimmedMETs", "offlineSlimmedPrimaryVertices"}

    if not hasattr(process, "genNanoTask"):
        raise RuntimeError("produceMiniAODGenInputs: the job has no gen tables (a withGen tier)")
    if process.genNanoTask.contains(process.metMCTable):
        raise RuntimeError("produceMiniAODGenInputs: the stock gen tables that read PAT products are "
                           "scheduled; it is meant for --procModifiers nano_l1_hlt, which removes them")
    defined = sorted(label for label in miniAODMarkers if hasattr(process, label))
    if defined:
        raise RuntimeError(f"produceMiniAODGenInputs: {defined} already defined in this job (PAT in the job?)")
    if not assumeNoMiniAOD:
        fileNames = list(getattr(process.source, "fileNames", []))
        if not fileNames:
            raise RuntimeError("produceMiniAODGenInputs: the job has no input file to check for MINIAOD products")
        # product branches are <type>_<module label>_<instance>_<process>.
        labels = {b.split("_")[1] for b in inputEventBranches(fileNames[0]) if b.count("_") >= 3}
        present = sorted(labels & miniAODMarkers)
        if present:
            raise RuntimeError(f"produceMiniAODGenInputs: {fileNames[0]} carries MINIAOD products {present}; "
                               "they must not be shadowed -- run without produceMiniAODGenInputs (note: under "
                               "nano_l1_hlt the stock GenMET/GenJetAK8/Pileup tables are then absent; see "
                               "DPGAnalysis/Phase3SmartPixelsNanoAOD/README.md)")

    for label, module in producers.items():
        setattr(process, label, module.clone())
    process.genMetTrueTable = simpleSingletonCandidateFlatTableProducer.clone(
        src="genMetTrue",
        name=process.metMCTable.name,
        doc=process.metMCTable.doc,
        variables=cms.PSet(pt=process.metMCTable.variables.pt.clone(expr="pt"),
                           phi=process.metMCTable.variables.phi.clone(expr="phi")))
    process.genPileupTable = process.puTable.clone(savePUDensityVars=False)
    process.smartPixelsMiniAODGenTask = cms.Task(
        *[getattr(process, label) for label in producers],
        process.genJetAK8Table, process.genJetAK8FlavourAssociation, process.genJetAK8FlavourTable,
        process.genMetTrueTable, process.genPileupTable)
    process.l1tPh2NanoTask.add(process.smartPixelsMiniAODGenTask)
    return process


# ---------------------------------------------------------------------------
# Untruncated IT cluster table  (@SPClustersStubs)
# ---------------------------------------------------------------------------
# EVERY IT pixel cluster, with no window, no truncation and no truth-pT filter:
# all TBPX clusters first, then the disc clusters in discScope (default sp =
# ID1-ID4), one table with a surface column (1-4 = IL1-IL4, 11-14 = ID1-ID4).
#
# ~74k rows/event at PU200 (TBPX + ID1-ID4), measured 3.1 MB/event.
# The size is the point: at PU200 the crossed module
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


# The COMPLETE OT stub collection (@TrkStubs, @SPClustersStubs).
# ~15k rows/event at PU200 expected.
#
# The questions the SmartPixels tiers exist to answer are comparative: an IT seeding design is only interesting
# relative to what the OT tracklet already does, and that comparison needs the OT
# INPUT multiplicity. The track -> stub link tables (L1TTrackStub, ...) are on-track
# only -- 1004 rows/event at PU200 against ~190 tracks -- and read every stub from
# here (otStubIdx); the persisted edmNew::DetSetVector is unreadable by uproot.
l1tPh3SmartPixelsAllStubTable = cms.EDProducer(
    "L1SmartPixelsAllStubTableProducer",
    stubs = cms.InputTag("TTStubsFromPhase2TrackerDigis", "StubAccepted"),
    stubTruth = cms.InputTag("TTStubAssociatorFromPixelDigis", "StubAccepted"),
    clusterTruth = cms.InputTag("TTClusterAssociatorFromPixelDigis", "ClusterAccepted"),
    tableName = cms.string("L1TOTStub"),
    doTruth = cms.bool(True),
)

# The DTC stub frames the tracklet track finder reads (L1TOTStubDTC): the integer coordinates, the
# bend code, the nonant routing, DTC acceptance and the stream order, keyed to L1TOTStub rows. It
# needs the in-job TTDTC (transient), so the flavours schedule it only next to the tracklet chain,
# reading the same TTDTC tag the chain reads.
l1tPh3SmartPixelsDTCStubTable = cms.EDProducer(
    "L1SmartPixelsDTCStubTableProducer",
    dtc = cms.InputTag("ProducerDTC", "StubAccepted"),
    stubs = cms.InputTag("TTStubsFromPhase2TrackerDigis", "StubAccepted"),
    tableName = cms.string("L1TOTStubDTC"),
    withPositions = cms.bool(False),
)


def addPh3L1SmartPixelsClusters(process, recHitLabel="spixSmartPixelsRecHits", doTruth=True):
    """Add the untruncated IT cluster table + the complete OT stub table (coexist jobs;
    the @SPClustersStubs flavour adds the same tables).

    Requires the SmartPixels-owned cluster -> rec-hit chain (spixPixelClusters ->
    spixPixelRecHits). That chain is created by
    customizeSmartPixels_cff._addPixelRecHitChain whenever a seedlessRefit variant is
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

