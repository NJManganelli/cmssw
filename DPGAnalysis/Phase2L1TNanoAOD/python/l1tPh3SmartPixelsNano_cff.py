import FWCore.ParameterSet.Config as cms

# Phase-3 SmartPixels L1 track nano tables.
#
# SmartPixels tracking is a beyond-Phase-2 (Phase-3) study, so it is kept out of the
# general Phase-2 L1Nano cff files and enabled independently via the autoNANO
# SmartPix flavours (see PhysicsTools/NanoAOD/autoNANO.py) or by calling
# addPh3L1SmartPixelsTracks as a customise. Requires the
# l1tSmartPixelsTrackProducer{,Extended} from the smartpixels development branch.
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

def addPh3L1SmartPixelsTracks(process):
    process.l1tPh2NanoTask.add(p3L1SmartPixelsTracksTask)
    return process
