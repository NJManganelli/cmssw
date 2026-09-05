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
    clusterMergeFrac=cms.double(0.1),
    measAngleMaxAbs=cms.double(12.0),
    angleSet=cms.string(""),   # REQUIRED: PixelAV angle-response payload
    noiseSet=cms.string(""),   # noise-angle inverse CDF; empty = unlinked clusters get NO
                               # angle, which makes "has an angle" a truth proxy
)
