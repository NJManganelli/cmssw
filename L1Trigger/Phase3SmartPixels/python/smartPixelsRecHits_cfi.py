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
    measAngleMaxAbs=cms.double(12.0),
    angleSet=cms.string(""),   # REQUIRED: PixelAV angle-response payload
    noiseSet=cms.string(""),   # noise-angle inverse CDF; empty = unlinked clusters get NO
                               # angle, which makes "has an angle" a truth proxy
    angleSource=cms.string("payload"),  # "cnn": cluster CNN, no truth (see the plugin header)
    pixelDigis=cms.InputTag("simSiPixelDigis", "Pixel"),
    cnnBarrelModel=cms.string(""),
    cnnDiscModel=cms.string(""),
    cnnBarrelScalars=cms.vstring(),
    cnnDiscScalars=cms.vstring(),
    cnnFoldOffset=cms.vdouble(0., 0., 0.202, 0.202, 0.202),
    cnnSigmaSet=cms.string(""),
)
