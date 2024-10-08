import FWCore.ParameterSet.Config as cms

hltTrackerClusterCheck = cms.EDProducer("ClusterCheckerEDProducer",
    ClusterCollectionLabel = cms.InputTag("siStripClusters"),
    MaxNumberOfPixelClusters = cms.uint32(40000),
    MaxNumberOfStripClusters = cms.uint32(400000),
    PixelClusterCollectionLabel = cms.InputTag("hltSiPixelClusters"),
    cut = cms.string('strip < 400000 && pixel < 40000 && (strip < 50000 + 10*pixel) && (pixel < 5000 + 0.1*strip)'),
    doClusterCheck = cms.bool(False),
    mightGet = cms.optional.untracked.vstring,
    silentClusterCheck = cms.untracked.bool(False)
)
