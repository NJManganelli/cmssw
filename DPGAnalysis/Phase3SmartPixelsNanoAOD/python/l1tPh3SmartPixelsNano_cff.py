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
