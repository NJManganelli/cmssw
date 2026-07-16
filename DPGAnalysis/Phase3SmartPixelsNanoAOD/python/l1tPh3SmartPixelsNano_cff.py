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
    promptTableName = f"l1tPh3SmartPixelsTracksTable{nameSuffix}"
    extendedTableName = f"l1tPh3ExtSmartPixelsTracksTable{nameSuffix}"

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
