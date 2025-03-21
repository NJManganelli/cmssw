import FWCore.ParameterSet.Config as cms

l1tSmartPixelsTrackProducer = cms.EDProducer('L1SmartPixelsTrackProducer',
    # MyProcess is the (unsigned) PDGID corresponding to the process which is run
    # e.g. single electron/positron = 11
    #      single pion+/pion- = 211
    #      single muon+/muon- = 13
    #      pions in jets = 6
    #      taus = 15
    #      all TPs = 1 (pp collisions)
    MyProcess = cms.int32(1),
    DebugMode = cms.bool(False),      # printout lots of debug statements
    L1Tk_minNStub = cms.int32(4),     # L1 tracks with >= 4 stubs
    # 4Param Fit Tracks, HYBRID
    L1TrackInputTag = cms.InputTag("l1tTTTracksFromTrackletEmulation",  "Level1TTTracks"),         # TTTrack input
    MCTruthTrackInputTag = cms.InputTag( "TTTrackAssociatorFromPixelDigis",  "Level1TTTracks"),  # MCTruth input
    # other input collections
    L1StubInputTag = cms.InputTag("TTStubsFromPhase2TrackerDigis","StubAccepted"),
    MCTruthClusterInputTag = cms.InputTag("TTClusterAssociatorFromPixelDigis", "ClusterAccepted"),
    MCTruthStubInputTag = cms.InputTag("TTStubAssociatorFromPixelDigis", "StubAccepted"),
    TrackingParticleInputTag = cms.InputTag("mix", "MergedTrackTruth"),
    outputCollectionName = cms.string("Level1TTTracks"),
    smartPixelsEmulatorMode = cms.string("passthrough"),  # passthrough, passthroughFloat, passthroughHW, trackingParticleTruth
)

l1tSmartPixelsTrackProducerExtended = l1tSmartPixelsTrackProducer.clone(
    # 5Param Fit Tracks, HYBRID
    L1TrackInputTag = cms.InputTag("l1tTTTracksFromExtendedTrackletEmulation",  "Level1TTTracks"),         # TTTrack input
    MCTruthTrackInputTag = cms.InputTag( "TTTrackAssociatorFromPixelDigisExtended",  "Level1TTTracks"),  # MCTruth input
)

### https://github.com/cms-L1TK/cmssw/blob/L1TK-dev-14_2_0_pre4/L1Trigger/TrackFindingTracklet/python/Producer_cfi.py
#process.load( 'L1Trigger.TrackFindingTracklet.Producer_cff' )
#process.load( 'L1Trigger.TrackFindingTracklet.Analyzer_cff' )
#NHELIXPAR = 4
#L1TRK_NAME  = process.TrackFindingTrackletAnalyzer_params.OutputLabelTFP.value()
#L1TRK_LABEL = process.TrackFindingTrackletProducer_params.BranchTTTracks.value()
#L1TRUTH_NAME = "TTTrackAssociatorFromPixelDigis"
#process.TTTrackAssociatorFromPixelDigis.TTTracks = cms.VInputTag( cms.InputTag(L1TRK_NAME, L1TRK_LABEL) )
#process.HybridNewKF = cms.Sequence(process.L1THybridTracks + process.ProducerTM + process.ProducerDR + process.ProducerKF + process.ProducerTQ + process.ProducerTFP)
#process.TTTracksEmulation = cms.Path(process.HybridNewKF)
#process.TTTracksEmulationWithTruth = cms.Path(process.HybridNewKF +  process.TrackTriggerAssociatorTracks)