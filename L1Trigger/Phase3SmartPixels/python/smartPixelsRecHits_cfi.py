import FWCore.ParameterSet.Config as cms

# SmartPixels rec hits: IT pixel rec hits + the sensor angle estimate, produced
# ONCE so the refit and every nano table read the same object. See the plugin
# header for why the angle is helix-propagated rather than taken at the parent
# production vertex.
smartPixelsRecHits = cms.EDProducer(
    "SmartPixelsRecHitProducer",
    pixelRecHits=cms.InputTag("spixPixelRecHits"),
    pixelDigiSimLink=cms.InputTag("simSiPixelDigis", "Pixel"),
    trackingParticles=cms.InputTag("mix", "MergedTrackTruth"),
    maxLayer=cms.uint32(4),
    discScope=cms.string("sp"),  # disc clusters: none | sp (ID1-ID4) | tfpx | all; payload gives them no angle
    clusterMergeFrac=cms.double(0.1),
    angleSet=cms.string(""),   # REQUIRED: PixelAV angle-response payload
    noiseSet=cms.string(""),   # noise-angle inverse CDF; empty = unlinked clusters get NO
                               # angle, which makes "has an angle" a truth proxy
    angleSource=cms.string("payload"),  # "emu": SmartPixels-output emulator (see the plugin header)
    pixelDigis=cms.InputTag("simSiPixelDigis", "Pixel"),
    # emu only (customizeSmartPixels_cff.SPIX_EMU_ANGLE_SOURCE fills the files)
    emuBarrelModel=cms.string(""),
    emuDiscModel=cms.string(""),
    emuBarrelA1=cms.string(""),
    emuDiscA1=cms.string(""),
    emuBarrelScalars=cms.vstring(),
    emuDiscScalars=cms.vstring(),
    emuLibrary=cms.string(""),
    emuSigmaSet=cms.string(""),
    emuInflation=cms.double(1.5),
    emuDiscAlphaOffset=cms.double(-0.312),
    emuSatEpsDeg=cms.double(0.01),
    emuDumpK=cms.uint32(32),
    emuTruthMaxHelixDist=cms.double(1.0),
    emuHelixStep=cms.double(0.2),
    emuDebug=cms.bool(False),
)
