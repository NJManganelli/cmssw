import FWCore.ParameterSet.Config as cms

smartPixelsClusterCensusAnalyzer = cms.EDAnalyzer(
    "SmartPixelsClusterCensusAnalyzer",
    pixelRecHitInputTag=cms.InputTag("spixPixelRecHits"),
    pixelDigiSimLinkInputTag=cms.InputTag("simSiPixelDigis", "Pixel"),
    trackingParticleInputTag=cms.InputTag("mix", "MergedTrackTruth"),
    simTrackInputTag=cms.InputTag("g4SimHits"),
    ptThresholds=cms.vdouble(0.5, 1.0, 1.5, 2.0),
)
