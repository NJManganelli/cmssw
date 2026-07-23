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

p3L1SmartPixelsTracksTask = cms.Task(
    l1tPh3SmartPixelsTracksTable,
    l1tPh3ExtSmartPixelsTracksTable,
)

def addPh3L1SmartPixelsTracks(process,
                              srcLabel="l1tSmartPixelsTrackProducer",
                              srcExtendedLabel="l1tSmartPixelsTrackProducerExtended",
                              tableSuffix=""):
    """Add one pair (prompt+extended) of SmartPixels track tables to the L1Nano task.

    srcLabel/srcExtendedLabel select the SmartPixels collections (any variant
    produced by L1Trigger/Phase3SmartPixels, e.g.
    l1tSmartPixelsTrackProducerWcorrectionlibRegressionAAII); tableSuffix
    disambiguates the table/branch names when several variants coexist in one
    job (WF1). Defaults reproduce the plain flavour reading the default
    producer labels.
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

    task = cms.Task(getattr(process, promptTableName), getattr(process, extendedTableName))
    setattr(process, f"p3L1SmartPixelsTracksTask{nameSuffix}", task)
    process.l1tPh2NanoTask.add(task)
    return process


# ---------------------------------------------------------------------------
# digiRefit sidecar tables (spec RefitSidecarSpec.md v0.1 §4.4)
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
# copies. Under posture C (fromFileStubs) these reference labels are re-run in-job
# from the file's stubs, so the stub geometry is the fresh new-layout one.
l1tPh3SmartPixelsStubTable = cms.EDProducer(
    "L1SmartPixelsStubPosTableProducer",
    tracks = cms.InputTag("l1tTTTracksFromTrackletEmulation", "Level1TTTracks"),
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
    return process


# ---------------------------------------------------------------------------
# PF/Puppi/jet re-emulation for the reduced-menu PU RelVals (posture-C spirit)
# ---------------------------------------------------------------------------
# The 200PU RelVals persist l1tLayer1:PuppiRegional (the per-region Puppi PF
# candidates) but NOT the flat l1tLayer2Deregionizer:Puppi collection, and their
# pre-persisted SC4/SC8 jets have constituent Ptrs dangling into that absent
# deregionizer product -- so consuming the file jets yields candIdx == -1 in the
# L1JetCandLinkTableProducer. We therefore RE-EMULATE the flat Puppi candidates
# and the SeededCone jets in-job from the persisted regional Puppi, exactly as
# posture C rebuilds tracks from the file's stubs:
#
#   l1tLayer1:PuppiRegional (FILE) -> DeregionizerProducer (l1tLayer2Deregionizer)
#                                  -> L1SeedConePFJetProducer SC4 / SC8 (corrected)
#
# The in-job jets' constituents Ptr into the fresh in-job deregionizer product, so
# the L1PuppiCand candIdx crossref resolves. Same-label shadowing hides any file
# HLT branch for downstream default-configured consumers, like posture C's tracks.
# Only the PROMPT (non-extended) Puppi is re-run; L1PuppiCand.l1TrackIdx points at
# the posture-C re-run l1tTTTracksFromTrackletEmulation. This is the hook the
# PF-carrying SmartPix flavor wiring calls (see customize TODO in the report).
# ---------------------------------------------------------------------------
# COOPT single-coherent-view: re-run the WHOLE Layer-1 correlator + emulated PV
# FROM one config's (refit) tracks, so vertex->PF->PUPPI->jets follow that view.
# ---------------------------------------------------------------------------
# The posture-C PU RelVal persists l1tLayer1:PuppiRegional + l1tVertexFinderEmulator
# built from the FILE'S ORIGINAL tracklet tracks. reemulateJetSideForPFTier() reads
# that file Puppi -> the jets are the FILE-track PF view, IDENTICAL whether or not a
# refit variant was injected (empirically byte-identical coexist vs coopt). To make
# the downstream genuinely COHERENT with an injected (refit) track collection we must
# RE-RUN the correlator + emulated vertex in-job from those tracks. Both chains read
# tracks through exactly two InputTags that carry the reference label
# l1tTTTracksFromTrackletEmulation:Level1TTTracks:
#   correlator: l1tPFTracksFromL1Tracks (feeds every region of l1ctLayer1_cff)
#   vertex:     l1tGTTInputProducer -> l1tTrackSelectionProducer -> l1tVertexFinderEmulator
# We repoint those two to `trackSrc` (the injected refit collection). The re-run
# correlator emits l1tLayer1:PuppiRegional + l1tVertexFinderEmulator:L1VerticesEmulation
# with the SAME labels -> they shadow the file products, so reemulateJetSideForPFTier /
# stitchPFTierForPostureC / addNGJetTier consume the refit-driven Puppi/vertex with no
# further change. l1ctLayer1_cff imports clean on native aarch64 (no module-scope ONNX).
# Every non-track correlator input is persisted in the RelVal EXCEPT the barrel EM
# clusters (l1tPhase2GCTBarrelToCorrelatorLayer1Emulator:GCTEmDigiClusters); doBarrelEM
# re-emulates them from persisted GCT products (2 extra producers), else the barrel EM
# input is blanked (charged-hadron + HGCal jets are unaffected -- the tagger's
# track-driven content is intact either way).
def stitchCorrelatorFromTracks(process,
                               trackSrc="l1tTTTracksFromTrackletEmulation",
                               extendedTrackSrc="l1tTTTracksFromExtendedTrackletEmulation",
                               doBarrelEM=False, doExtended=True):
    """Re-run Layer-1 correlator + emulated primary vertex from `trackSrc` so the
    whole PF/Puppi/jet chain is coherent with that (refit) track view. Must run
    BEFORE reemulateJetSideForPFTier/stitchPFTierForPostureC/addNGJetTier (it
    produces the l1tLayer1:PuppiRegional + emulated vertex those consume). Returns
    process. Intended for the coopt (WF2) single-coherent-view productions."""
    import FWCore.ParameterSet.Config as cms
    process.load("L1Trigger.Phase2L1ParticleFlow.l1ctLayer1_cff")
    # Emulated PV chain (GTTInput -> selection -> vertex), track-driven.
    from L1Trigger.L1TTrackMatch.l1tGTTInputProducer_cfi import l1tGTTInputProducer
    from L1Trigger.L1TTrackMatch.l1tTrackSelectionProducer_cfi import l1tTrackSelectionProducer
    from L1Trigger.VertexFinder.l1tVertexProducer_cfi import l1tVertexProducer

    if not hasattr(process, "l1tGTTInputProducer"):
        process.l1tGTTInputProducer = l1tGTTInputProducer.clone()
    if not hasattr(process, "l1tTrackSelectionProducer"):
        process.l1tTrackSelectionProducer = l1tTrackSelectionProducer.clone()
    if not hasattr(process, "l1tVertexFinderEmulator"):
        process.l1tVertexFinderEmulator = l1tVertexProducer.clone(VertexReconstruction=dict(Algorithm="fastHistoEmulation"))

    # Repoint the ONLY two track tags carrying the reference label to `trackSrc`
    # (the correlator PFTrack decoder + the GTTInput vertex feeder). Everything
    # else in the correlator resolves to persisted calo/HGCal/muon products.
    process.l1tPFTracksFromL1Tracks.L1TrackTag = cms.InputTag(trackSrc, "Level1TTTracks")
    process.l1tGTTInputProducer.l1TracksInputTag = cms.InputTag(trackSrc, "Level1TTTracks")
    if doExtended and hasattr(process, "l1tPFTracksFromL1TracksExtended"):
        process.l1tPFTracksFromL1TracksExtended.L1TrackTag = cms.InputTag(extendedTrackSrc, "Level1TTTracks")

    if doBarrelEM:
        # Re-emulate the barrel EM clusters (the one non-persisted correlator input)
        # from persisted GCT products. Inputs verified persisted: GCTFullTowers,
        # GCTClusters, GCTDigitizedClusterToCorrelator, simHcalTriggerPrimitiveDigis.
        from L1Trigger.L1CaloTrigger.l1tPhase2CaloPFClusterEmulator_cfi import l1tPhase2CaloPFClusterEmulator
        from L1Trigger.L1CaloTrigger.l1tPhase2GCTBarrelToCorrelatorLayer1Emulator_cfi import l1tPhase2GCTBarrelToCorrelatorLayer1Emulator
        if not hasattr(process, "l1tPhase2CaloPFClusterEmulator"):
            process.l1tPhase2CaloPFClusterEmulator = l1tPhase2CaloPFClusterEmulator.clone()
        if not hasattr(process, "l1tPhase2GCTBarrelToCorrelatorLayer1Emulator"):
            process.l1tPhase2GCTBarrelToCorrelatorLayer1Emulator = \
                l1tPhase2GCTBarrelToCorrelatorLayer1Emulator.clone()
    else:
        # No barrel EM: blank the tag so the barrel correlator runs without it.
        # HGCal/HF hadronic + all charged (track) PF are unaffected.
        process.l1tLayer1Barrel.emClusters = cms.InputTag("")

    # Schedule the correlator regions + vertex trio unscheduled on the nano task
    # (data dependency runs them ahead of the deregionizer/PF-cand tables). We do
    # NOT re-run the calo-cluster producers (l1tPFClustersFromCombinedCaloHCal/HF,
    # l1tPFClustersFromL1EGClusters): their outputs are PERSISTED in the RelVal and
    # re-running them needs the CaloTPGRecord ES setup not loaded in a NANO job.
    # Only the correlator regions (which read TRACKS -> the injected refit view) and
    # the track PFTrack decoder are re-run; calo/HGCal/muon inputs resolve to the
    # file's persisted products.
    corrModules = ["l1tPFTracksFromL1Tracks",
                   "l1tLayer1Barrel", "l1tLayer1HGCal", "l1tLayer1HGCalNoTK",
                   "l1tLayer1HF", "l1tLayer1"]
    if doExtended:
        corrModules += ["l1tPFTracksFromL1TracksExtended", "l1tLayer1BarrelExtended",
                        "l1tLayer1HGCalExtended", "l1tLayer1Extended"]
    vtxModules = ["l1tGTTInputProducer", "l1tTrackSelectionProducer", "l1tVertexFinderEmulator"]
    emModules = (["l1tPhase2CaloPFClusterEmulator",
                  "l1tPhase2GCTBarrelToCorrelatorLayer1Emulator"] if doBarrelEM else [])
    members = [m for m in (corrModules + vtxModules + emModules) if hasattr(process, m)]
    process.smartPixelsCorrelatorTask = cms.Task(*[getattr(process, m) for m in members])
    if hasattr(process, "l1tPh2NanoTask"):
        process.l1tPh2NanoTask.add(process.smartPixelsCorrelatorTask)
    else:
        raise RuntimeError("stitchCorrelatorFromTracks: l1tPh2NanoTask not present.")
    print(f"SmartPixels coopt correlator: re-ran Layer-1 correlator + emulated PV from "
          f"'{trackSrc}' (barrelEM={doBarrelEM}); l1tLayer1:PuppiRegional + "
          f"l1tVertexFinderEmulator now refit-driven (shadow the file products).")
    return process


def reemulateJetSideForPFTier(process,
                              regionalPuppi=("l1tLayer1", "PuppiRegional"),
                              doSC4=True, doSC8=True):
    """Schedule the deregionizer + SC4/SC8 SeededCone jet producers in-job so the
    L1PuppiCand / L1SC4JetCands / L1SC8JetCands / L1puppiJetSC{4,8} nano tables are
    filled on a reduced-menu (posture-C) PU input. Idempotent: skips a producer if
    its default label already exists (e.g. a full-menu input already ran it)."""
    # Build from the BASE cfi modules, NOT l1pfJetMet_cff: that cff imports the
    # ONNX-backed NG-jet producer at module scope, which bus-errors under native
    # aarch64 cpuid probing and is not needed here.
    from L1Trigger.Phase2L1ParticleFlow.l1tDeregionizerProducer_cfi import l1tDeregionizerProducer
    from L1Trigger.Phase2L1ParticleFlow.l1SeedConePFJetEmulatorProducer_cfi import l1SeedConePFJetEmulatorProducer
    _JEC = "L1Trigger/Phase2L1ParticleFlow/data/jecs/jecs_20220308.root"
    seq = []
    if not hasattr(process, "l1tLayer2Deregionizer"):
        process.l1tLayer2Deregionizer = l1tDeregionizerProducer.clone(
            RegionalPuppiCands=cms.InputTag(*regionalPuppi))
        seq.append(process.l1tLayer2Deregionizer)
    if doSC4 and not hasattr(process, "l1tSC4PFL1PuppiCorrectedEmulator"):
        process.l1tSC4PFL1PuppiCorrectedEmulator = l1SeedConePFJetEmulatorProducer.clone(
            L1PFObjects="l1tLayer2Deregionizer:Puppi",
            doCorrections=cms.bool(True), correctorFile=cms.string(_JEC),
            correctorDir=cms.string('L1PuppiSC4EmuJets'))
        seq.append(process.l1tSC4PFL1PuppiCorrectedEmulator)
    if doSC8 and not hasattr(process, "l1tSC8PFL1PuppiCorrectedEmulator"):
        process.l1tSC8PFL1PuppiCorrectedEmulator = l1SeedConePFJetEmulatorProducer.clone(
            L1PFObjects="l1tLayer2Deregionizer:Puppi",
            coneSize=cms.double(0.8), wideConeJet=cms.bool(True),
            doCorrections=cms.bool(True), correctorFile=cms.string(_JEC),
            correctorDir=cms.string('L1PuppiSC4EmuJets'))
        seq.append(process.l1tSC8PFL1PuppiCorrectedEmulator)
    if seq:
        # Schedule the re-emulation modules UNSCHEDULED (on a Task), so the
        # framework runs them by data dependency ahead of the nano PF-cand tables
        # that consume l1tLayer2Deregionizer:Puppi / the corrected jets. A bare
        # cms.Path attribute is NOT added to process.schedule automatically and
        # would leave the products unproduced (ProductNotFound at output). Prefer
        # the nano task when present; else fall back to a scheduled Path.
        reemulTask = cms.Task(*seq)
        process.smartPixelsJetReemulTask = reemulTask
        if hasattr(process, "l1tPh2NanoTask"):
            process.l1tPh2NanoTask.add(reemulTask)
        else:
            process.smartPixelsJetReemulPath = cms.Path(reemulTask)
            if hasattr(process, "schedule") and process.schedule is not None:
                process.schedule.append(process.smartPixelsJetReemulPath)
    return process


# Nano table modules whose source objects are absent from some RelVal L1 menus
# (posture-A fromFile productions read the file's HLT-process objects, and older
# RelVals predate the NGJet producer / HPS PF taus). Removing the corresponding
# tables avoids a ProductNotFound at output time. hpsTauTable is already dropped
# by the nano_l1_hlt modifier; sc4NGJetTable + its cand link table are not.
# NOTE: l1tSC8JetCandsTable / l1tSC4JetCandsTable (plain seededcone link tables)
# are intentionally NOT in this default drop list -- when the PF-carrying tier
# calls reemulateJetSideForPFTier() their deregionizer Puppi source exists in-job.
_ABSENT_MENU_TABLES_DEFAULT = ("sc4NGJetTable", "l1tSC4NGJetCandsTable", "hpsTauTable")


def stitchPFTierForPostureC(process, doSC8=True):
    """Wire the PF/Puppi/jet nano tables for a posture-C (fromFileStubs) PU run of
    a PF-carrying SmartPix flavor (L1PFNano* / L1PFTrkNano*).

    The 200PU RelVals persist l1tLayer1:PuppiRegional but not the flat
    l1tLayer2Deregionizer:Puppi, so the PF-cand tables would otherwise resolve to
    dangling refs (candIdx == -1) or ProductNotFound. This helper:
      1. re-emulates the PROMPT deregionizer + SC4/SC8 corrected SeededCone jets
         in-job (reemulateJetSideForPFTier), so l1tLayer2Deregionizer:Puppi and the
         corrected jets exist as fresh in-job products;
      2. schedules p2L1PFCandsTask on the nano task and UN-PRUNES the L1PuppiCand /
         SC4 / SC8 link-table family (l1tPuppiCandsTable, l1tSC4JetCandsTable,
         l1tSC8JetCandsTable, l1tPuppiCandTrackTruthTable) so those columns write.
         l1tPuppiCandTrackTruthTable is now scheduled (no longer dropped): its
         producer isAvailable()-guards the PFTrack deref and self-reports the
         product-level failure via trkTruthStatus == -2 rather than throwing;
      3. drops the members that depend on products still absent under posture C
         (extended deregionizer, layer-1 PF, NG-tagged jets, HGCal clusters):
         l1tExtPuppiCandsTable, l1tPFCandsTable, l1tSC4NGJetCandsTable,
         l1tHGCClusterTable, l1tPuppiCandHGCClusterLink,
         l1tExtPuppiCandHGCClusterLink, l1tExtPuppiCandTrackTruthTable;
      4. leaves exactly ONE plain link table writing the L1PuppiCand extension
         (l1tSC4JetCandsTable: writeCandExtension=True; l1tSC8JetCandsTable:
         writeCandExtension=False) so the jetIdx/l1TrackIdx extension columns are
         written once. On posture-C PU the file's PFTrack refs dangle, so the
         candidate-extension deref is isAvailable()-guarded and l1TrackIdx == -2
         (product-level failure) for all candidates; -1 marks an element-level
         no-match when the product is present.
    """
    process = reemulateJetSideForPFTier(process, doSC4=True, doSC8=doSC8)

    from DPGAnalysis.Phase2L1TNanoAOD.l1tPh2PFCandsNanotables_cff import p2L1PFCandsTask
    if not hasattr(process, "l1tPh2NanoTask"):
        raise RuntimeError("stitchPFTierForPostureC: l1tPh2NanoTask not present; "
                           "the PF-carrying nano flavor did not schedule its base task.")
    # Schedule the whole PF-cand task, then drop the posture-C-unavailable members.
    process.l1tPh2NanoTask.add(p2L1PFCandsTask)
    # l1tPuppiCandTrackTruthTable (L1PFCandTrackTruthTableProducer) is NO LONGER dropped:
    # the producer now isAvailable()-guards the PFCandidate->PFTrack->TTTrack deref (the
    # tolerance-family fix, cf. the SC4 link-table candidate extension). On posture-C PU
    # RelVals the file's PFTrack collection is not stored, so the producer detects the
    # product-level failure and writes trkTruthStatus == -2 (with the truth columns at
    # their unknown defaults) instead of throwing ProductNotFound. It is scheduled again.
    # l1tExtPuppiCandTrackTruthTable stays dropped: it targets the extended deregionizer
    # product, which is genuinely absent under posture C (not merely dangling).
    _drop = ["l1tExtPuppiCandsTable", "l1tPFCandsTable", "l1tSC4NGJetCandsTable",
             "l1tHGCClusterTable", "l1tPuppiCandHGCClusterLink",
             "l1tExtPuppiCandHGCClusterLink", "l1tExtPuppiCandTrackTruthTable"]
    if not doSC8:
        _drop.append("l1tSC8JetCandsTable")
    process = dropAbsentMenuTables(process, _drop)
    print("SmartPixels posture-C PF tier: re-emulated deregionizer+SC4"
          + ("+SC8" if doSC8 else "")
          + "; kept L1PuppiCand + SC4"
          + ("/SC8" if doSC8 else "")
          + " link tables (SC4 writes the plain-cand extension) + PuppiCand track-truth.")
    return process


def dropAbsentMenuTables(process, moduleLabels=_ABSENT_MENU_TABLES_DEFAULT):
    """Remove nano table modules (by label) for L1 objects this input file lacks.

    Removes each module from every task/sequence/path/endpath then deletes the
    attribute, so the NanoAOD output module no longer consumes its product.
    Silent no-op for labels not present. Intended for posture-A SmartPixels
    productions on RelVals with a reduced L1 menu (see the l1nano workflow memory).
    """
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
        print(f"SmartPixels posture-A: removed nano table '{label}' (source object absent from this input menu)")
    return process


# FlatTable-producer C++ types in the Phase-2/3 L1 nano packages whose primary
# input is a single L1 object collection selected by an InputTag parameter. The
# value is the parameter name carrying that InputTag. The auto-prune below only
# considers these; link/extension/truth producers with multiple/optional inputs
# are left alone (and are handled by the explicit drop list when needed).
_SIMPLE_TABLE_SRC_PARAM = {
    "SimpleTriggerL1SAMuonFlatTableProducer": "src",
    "SimpleTriggerL1TrackerMuonFlatTableProducer": "src",
    "SimpleTriggerL1PFJetFlatTableProducer": "src",
    "SimpleTriggerL1PFTauFlatTableProducer": "src",
    "SimpleTriggerL1HPSPFTauFlatTableProducer": "src",
    "SimpleTriggerL1EGFlatTableProducer": "src",
    "SimpleTriggerL1PFCandidateFlatTableProducer": "src",
    "SimpleTriggerL1VertexWordFlatTableProducer": "src",
    "SimpleL1TTTrackCandidateFlatTableProducer": "src",
    "SimpleL1DisplacedVtxCandidateFlatTableProducer": "src",
    "SimpleTriggerL1EtSumFlatTableProducer": "src",
    "SimpleTriggerL1HGCalMulticlusterFlatTableProducer": "src",
}


def pruneAbsentSimpleTables(process, availableLabels, protectLabels=()):
    """Auto-drop simple L1 FlatTable producers whose single src collection is not
    available (posture-A on a reduced-menu RelVal).

    availableLabels: iterable of module labels the input FILE provides (from its
    provenance) -- a producer's src label is considered resolvable if it is a
    producer in the process OR appears in availableLabels. protectLabels are
    never dropped regardless (e.g. the SmartPixels + reference track tables).
    Only the _SIMPLE_TABLE_SRC_PARAM producer types are auto-pruned; everything
    else is untouched (drop those explicitly via dropAbsentMenuTables).
    """
    import FWCore.ParameterSet.Config as cms
    available = set(availableLabels)
    protect = set(protectLabels)
    to_drop = []
    for label, module in process.producers_().items():
        if label in protect:
            continue
        ctype = getattr(module, "_TypedParameterizable__type", None)
        param = _SIMPLE_TABLE_SRC_PARAM.get(ctype)
        if param is None or not hasattr(module, param):
            continue
        src = getattr(module, param)
        if not isinstance(src, cms.InputTag):
            continue
        srcLabel = src.getModuleLabel()
        if srcLabel in available or hasattr(process, srcLabel):
            continue
        to_drop.append(label)
    return dropAbsentMenuTables(process, to_drop)


def useGenParticlesFromFile(process, genSrc="genParticles"):
    """Adapt the withGen gen chain to a GEN-SIM(-DIGI-RAW) input that has
    `genParticles` but NOT the MINIAOD `prunedGenParticles` (posture-A RelVals).

    Repoints the finalGenParticles pruner at `genSrc` so genParticleTable and the
    gen-iso/tau tables resolve, then auto-prunes any remaining gen/PAT-derived
    tables whose source is still absent (e.g. slimmed collections). No-op if the
    gen chain was not scheduled (non-withGen tiers).
    """
    import FWCore.ParameterSet.Config as cms
    if hasattr(process, "finalGenParticles"):
        process.finalGenParticles.src = cms.InputTag(genSrc)
        print(f"SmartPixels posture-A: repointed finalGenParticles.src -> '{genSrc}' (no MINIAOD prunedGenParticles)")
    # genIso (GenPartIsoProducer) needs MINIAOD packedGenParticles; drop it and
    # the genParticleTable.iso external variable that consumes it. Core GenPart
    # columns (kinematics, pdgId, status, mother, statusFlags) are unaffected.
    if hasattr(process, "genParticleTable") and hasattr(process.genParticleTable, "externalVariables"):
        if hasattr(process.genParticleTable.externalVariables, "iso"):
            del process.genParticleTable.externalVariables.iso
            print("SmartPixels posture-A: dropped genParticleTable.externalVariables.iso (no packedGenParticles)")
    process = dropAbsentMenuTables(process, ("genIso",))
    # The gen-jet / gen-jet-flavour chain is MINIAOD-coupled (slimmedGenJets et
    # al.); it is not needed for SmartPixels track truth-matching (genParticles +
    # the L1 track truth table are). Drop it so the withGen tier runs on a
    # GEN-SIM RelVal. genParticles / gen taus / gen vertices stay.
    _GENJET_TABLES = ("genJetTable", "genJetFlavourTable", "patJetPartonsNano",
                      "genJetAK8Table", "genJetAK8FlavourAssociation", "genJetAK8FlavourTable",
                      "genJetFlavourAssociation")
    process = dropAbsentMenuTables(process, _GENJET_TABLES)
    return process


# ---------------------------------------------------------------------------
# GenJet flavour tier from GEN-SIM(-DIGI-RAW) content (unified-nano principle):
# jet-flavour labels for the tagger, WITHOUT MINIAOD.
# ---------------------------------------------------------------------------
# useGenParticlesFromFile() DROPS the whole AK4 genJet chain because the stock
# nano wiring reads MINIAOD (genJetTable.src=slimmedGenJets,
# genJetFlavourTable.jetFlavourInfos=slimmedGenJetsFlavourInfos,
# patJetPartonsNano.particles=prunedGenParticles). On a GEN-SIM RelVal those are
# absent BUT the file already persists reco::GenJet "ak4GenJets" and
# reco::GenParticle "genParticles" -- so the AK4 genJet + parton/hadron flavour
# tier is fully producible in-job by REPOINTING (not dropping): read the file's
# ak4GenJets and run HadronAndPartonSelector + JetFlavourClustering on genParticles.
# This is the (a) gap for the unified withGen nano: GenJet_partonFlavour/hadronFlavour.
def useGenJetsInJob(process, genJetSrc="ak4GenJets", genParticlesSrc="genParticles"):
    """Populate the AK4 GenJet + partonFlavour/hadronFlavour tables from GEN-SIM
    content (persisted reco::GenJet ak4GenJets + reco::GenParticle genParticles),
    replacing the MINIAOD-coupled default sources. Call AFTER useGenParticlesFromFile
    (which drops the MINIAOD genJet tables); this re-adds the producible ones.

    Wires, from PhysicsTools.NanoAOD.jetMC_cff:
      patJetPartonsNano (HadronAndPartonSelector) <- genParticlesSrc
      genJetFlavourAssociation (JetFlavourClustering) <- ak4GenJets + partons
      genJetFlavourTable (GenJetFlavourTableProducer) <- ak4GenJets + the assoc
      genJetTable (simpleGenJetFlatTableProducer, "GenJet") <- ak4GenJets
    The GenJet table then carries partonFlavour (int16) + hadronFlavour (uint8)
    aligned per index to the P4 columns -- the tagger's jet-flavour labels.
    """
    import FWCore.ParameterSet.Config as cms
    from PhysicsTools.NanoAOD.jetMC_cff import (
        genJetTable as _genJetTable, patJetPartonsNano as _patJetPartonsNano,
        genJetFlavourAssociation as _genJetFlavourAssociation,
        genJetFlavourTable as _genJetFlavourTable)

    # HadronAndPartonSelector reads genParticles (GEN-SIM) + the generator product.
    process.patJetPartonsNano = _patJetPartonsNano.clone(
        particles=cms.InputTag(genParticlesSrc))
    # JetFlavourClustering ghost-associates b/c hadrons + partons to the AK4 genJets.
    process.genJetFlavourAssociation = _genJetFlavourAssociation.clone(
        jets=cms.InputTag(genJetSrc))
    # GenJetFlavourTableProducer writes partonFlavour/hadronFlavour onto the GenJet
    # table (matched by deltaR to the assoc), keyed to the SAME ak4GenJets src+cut.
    process.genJetFlavourTable = _genJetFlavourTable.clone(
        src=cms.InputTag(genJetSrc),
        jetFlavourInfos=cms.InputTag("genJetFlavourAssociation"))
    # The GenJet flat table itself, from the file's ak4GenJets (P4 + flavour join).
    process.genJetTable = _genJetTable.clone(
        src=cms.InputTag(genJetSrc),
        doc=cms.string(f"AK4 GenJets ({genJetSrc}) from GEN-SIM visible genParticles"))

    task = cms.Task(process.patJetPartonsNano, process.genJetFlavourAssociation,
                    process.genJetFlavourTable, process.genJetTable)
    process.smartPixelsGenJetTask = task
    if hasattr(process, "l1tPh2NanoTask"):
        process.l1tPh2NanoTask.add(task)
    else:
        raise RuntimeError("useGenJetsInJob: l1tPh2NanoTask not present; the withGen "
                           "flavor did not schedule its base gen task.")
    print(f"SmartPixels GenJet tier: AK4 GenJet + partonFlavour/hadronFlavour from "
          f"'{genJetSrc}' + '{genParticlesSrc}' (in-job flavour clustering, no MINIAOD).")
    return process


# ---------------------------------------------------------------------------
# NG jet-tagger tier (unified-nano (b) gap): L1puppiJetSC4NG + L1SC4NGJetCands +
# L1ExtPuppiCand + L1HGCCluster on a reduced-menu (posture-C) PU input.
# ---------------------------------------------------------------------------
# The stock production DROPS sc4NGJetTable / l1tSC4NGJetCandsTable via
# dropAbsentMenuTables because the NG producer + its inputs are not scheduled.
# They ARE producible in this build: L1TSC4NGJetProducer is an hls4mlEmulator
# model (L1TSC4NGJetModel_v1_0_1.so, compiled for el9_aarch64_gcc13) loaded in the
# producer constructor -- NOT the module-scope ONNX import in l1pfJetMet_cff that
# bus-errors on native aarch64 (that is why we import the cfi STANDALONE here, as
# reemulateJetSideForPFTier does for the deregionizer/SeededCone cfis, never the
# full l1pfJetMet_cff). The NG producer reads corrected SeededCone l1t::PFJet jets
# (which reemulateJetSideForPFTier already produces as l1tSC4PFL1PuppiCorrectedEmulator)
# and emits per-jet b/c/uds/g/tau/e/mu tag scores. The link + cand + HGCcluster
# tables all resolve from products the RelVal persists (l1tLayer1[Extended]:PuppiRegional,
# l1tHGCalBackEndLayer2Producer 3D clusters) re-emulated in-job.
def addNGJetTier(process,
                 correctedJets=("l1tSC4PFL1PuppiCorrectedEmulator", ""),
                 promptDeregPuppi=("l1tLayer2Deregionizer", "Puppi"),
                 extRegionalPuppi=("l1tLayer1Extended", "PuppiRegional"),
                 hgcClusters=("l1tHGCalBackEndLayer2Producer",
                              "HGCalBackendLayer2Processor3DClustering"),
                 coherent=True):
    """Schedule the NG jet-tagger tier + un-prune its nano tables.

    coherent=True (single-coherent-view default): the NG jets are clustered from
    the SAME prompt corrected SeededCone jets that fill the base PF/Puppi/SC4 tier
    (correctedJets, produced by reemulateJetSideForPFTier), and the NG jet<->cand
    link table points its L1ExtPuppiCand candidate table at the PROMPT deregionizer
    Puppi -- so the whole jet tier follows one coherent track->PF->Puppi view.
    coherent=False reproduces the stock NG wiring (jets + cands from the EXTENDED
    deregionizer, re-emulated from extRegionalPuppi).

    Requires reemulateJetSideForPFTier() to have run (prompt deregionizer + SC4
    corrected jets) and the base p2L1PFCandsTask family scheduled (as
    stitchPFTierForPostureC does). Idempotent per producer/table label.
    """
    import FWCore.ParameterSet.Config as cms
    from L1Trigger.Phase2L1ParticleFlow.l1tSC4NGJetProducer_cfi import l1tSC4NGJetProducer
    from L1Trigger.Phase2L1ParticleFlow.l1tDeregionizerProducer_cfi import (
        l1tDeregionizerProducerExtended)

    seq = []
    ngJetsInput = cms.InputTag(*correctedJets) if correctedJets[1] else cms.InputTag(correctedJets[0])
    if not coherent:
        # Extended-deregionizer path (stock NG training input): re-emulate the
        # extended flat Puppi in-job from the file's extended regional Puppi.
        if not hasattr(process, "l1tLayer2DeregionizerExtended"):
            process.l1tLayer2DeregionizerExtended = l1tDeregionizerProducerExtended.clone(
                RegionalPuppiCands=cms.InputTag(*extRegionalPuppi))
            seq.append(process.l1tLayer2DeregionizerExtended)

    # The NG jet producer: hls4ml model in its ctor; reads corrected l1t::PFJet.
    if not hasattr(process, "l1tSC4NGJetProducer"):
        process.l1tSC4NGJetProducer = l1tSC4NGJetProducer.clone(jets=ngJetsInput)
        seq.append(process.l1tSC4NGJetProducer)

    if seq:
        task = cms.Task(*seq)
        process.smartPixelsNGJetTask = task
        if hasattr(process, "l1tPh2NanoTask"):
            process.l1tPh2NanoTask.add(task)
        else:
            raise RuntimeError("addNGJetTier: l1tPh2NanoTask not present.")

    # Un-prune the NG nano tables. They live in l1tPh2Nanotables_cff (jet-side
    # tables) and l1tPh2PFCandsNanotables_cff (cand/link/hgccluster). Import and
    # attach any not already present, repointing to the coherent products.
    from DPGAnalysis.Phase2L1TNanoAOD.l1tPh2Nanotables_cff import sc4NGJetTable as _sc4NG
    from DPGAnalysis.Phase2L1TNanoAOD.l1tPh2PFCandsNanotables_cff import (
        l1tSC4NGJetCandsTable as _ngCands, l1tExtPuppiCandsTable as _extPuppi,
        l1tHGCClusterTable as _hgc)

    addTables = []
    if not hasattr(process, "sc4NGJetTable"):
        _t = _sc4NG.clone()
        # The stock sc4NGJetTable inherits pfJetTable.externalVariables (btagScore/
        # llpTagScore from l1tBJetProducer*/l1tTOoLLiP* -- separate NN producers NOT
        # scheduled here). The multi-class NG scores come from `variables`
        # (getTagScore); drop the unresolvable external NN scores so it does not
        # ProductNotFound at output.
        if hasattr(_t, "externalVariables"):
            del _t.externalVariables
        process.sc4NGJetTable = _t
        addTables.append("sc4NGJetTable")

    if not hasattr(process, "l1tExtPuppiCandsTable"):
        # The L1ExtPuppiCand candidate table (constituents referenced by the NG
        # link table). In coherent mode point it at the PROMPT deregionizer Puppi
        # so the NG cands share the base-tier view; else the extended deregionizer.
        _e = _extPuppi.clone(
            src=cms.InputTag(*promptDeregPuppi) if coherent else _extPuppi.src)
        process.l1tExtPuppiCandsTable = _e
        addTables.append("l1tExtPuppiCandsTable")

    if not hasattr(process, "l1tSC4NGJetCandsTable"):
        _c = _ngCands.clone()
        if coherent:
            _c.cands = cms.InputTag(*promptDeregPuppi)
            _c.tracks = cms.InputTag("l1tTTTracksFromTrackletEmulation", "Level1TTTracks")
            _c.trackTableName = cms.string("L1TTrack")
        process.l1tSC4NGJetCandsTable = _c
        addTables.append("l1tSC4NGJetCandsTable")

    if not hasattr(process, "l1tHGCClusterTable"):
        process.l1tHGCClusterTable = _hgc.clone(src=cms.InputTag(*hgcClusters))
        addTables.append("l1tHGCClusterTable")

    if addTables:
        ngTablesTask = cms.Task(*[getattr(process, t) for t in addTables])
        process.smartPixelsNGTablesTask = ngTablesTask
        process.l1tPh2NanoTask.add(ngTablesTask)
    print("SmartPixels NG jet tier: L1puppiJetSC4NG + L1SC4NGJetCands + L1ExtPuppiCand"
          " + L1HGCCluster scheduled (coherent=%s); NG jets from %s."
          % (coherent, ngJetsInput))
    return process
