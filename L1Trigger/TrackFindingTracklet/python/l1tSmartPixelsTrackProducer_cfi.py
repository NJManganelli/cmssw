import FWCore.ParameterSet.Config as cms

l1tSmartPixelsTrackProducer = cms.EDAnalyzer('L1SmartPixelsTrackProducer',
      # MyProcess is the (unsigned) PDGID corresponding to the process which is run
      # e.g. single electron/positron = 11
      #      single pion+/pion- = 211
      #      single muon+/muon- = 13
      #      pions in jets = 6
      #      taus = 15
      #      all TPs = 1 (pp collisions)
       MyProcess = cms.int32(1),
       DebugMode = cms.bool(False),      # printout lots of debug statements
       SaveStubs = cms.bool(False),      # save some info for *all* stubs
       L1Tk_nPar = cms.int32(4),         # use 4 or 5-parameter L1 tracking?
       L1Tk_minNStub = cms.int32(4),     # L1 tracks with >= 4 stubs
       TP_minNStub = cms.int32(4),       # require TP to have >= X number of stubs associated with it
       TP_minNStubLayer = cms.int32(4),  # require TP to have stubs in >= X layers/disks
       TP_minPt = cms.double(1.9),       # only save TPs with pt > X GeV
       TP_maxEta = cms.double(2.5),      # only save TPs with |eta| < X
       TP_maxZ0 = cms.double(30.0),      # only save TPs with |z0| < X cm
       # 4Param Fit Tracks, HYBRID
       L1TrackInputTag = cms.InputTag("l1tTTTracksFromTrackletEmulation",  "Level1TTTracks"),         # TTTrack input
       MCTruthTrackInputTag = cms.InputTag( "TTTrackAssociatorFromPixelDigis",  "Level1TTTracks"),  # MCTruth input
       # 5Param Fit Tracks, HYBRID
       ##L1TrackInputTag = cms.InputTag("l1tTTTracksFromExtendedTrackletEmulation",  "Level1TTTracks"),         # TTTrack input
       ##MCTruthTrackInputTag = cms.InputTag( "TTTrackAssociatorFromPixelDigisExtended",  "Level1TTTracks"),  # MCTruth input
       # 4Param Fit Tracks, HYBRID_NEWKF / HYBRID_REDUCED
       ##L1TrackInputTag = cms.InputTag("???",  "???"),         # TTTrack input
       ##MCTruthTrackInputTag = cms.InputTag( "???",  "???"),  # MCTruth input
       # other input collections
       L1StubInputTag = cms.InputTag("TTStubsFromPhase2TrackerDigis","StubAccepted"),
       MCTruthClusterInputTag = cms.InputTag("TTClusterAssociatorFromPixelDigis", "ClusterAccepted"),
       MCTruthStubInputTag = cms.InputTag("TTStubAssociatorFromPixelDigis", "StubAccepted"),
       TrackingParticleInputTag = cms.InputTag("mix", "MergedTrackTruth"),
       ###TrackingVertexInputTag = cms.InputTag("mix", "MergedTrackTruth"),
       #### tracking in jets (--> requires AK4 genjet collection present!)
       ###TrackingInJets = cms.bool(False),
       ###GenJetInputTag = cms.InputTag("ak4GenJets", "")
)
### https://github.com/cms-L1TK/cmssw/blob/L1TK-dev-14_2_0_pre4/L1Trigger/TrackFindingTracklet/python/Producer_cfi.py
if (L1TRKALGO == 'HYBRID'):
    process.TTTracksEmulation = cms.Path(process.L1THybridTracks)
    process.TTTracksEmulationWithTruth = cms.Path(process.L1THybridTracksWithAssociators)
    NHELIXPAR = 4
    L1TRK_NAME  = "l1tTTTracksFromTrackletEmulation"
    L1TRK_LABEL = "Level1TTTracks"
    L1TRUTH_NAME = "TTTrackAssociatorFromPixelDigis"

# HYBRID: extended tracking
elif (L1TRKALGO == 'HYBRID_DISPLACED'):
    process.TTTracksEmulation = cms.Path(process.L1TExtendedHybridTracks)
    process.TTTracksEmulationWithTruth = cms.Path(process.L1TExtendedHybridTracksWithAssociators)
    NHELIXPAR = 5
    L1TRK_NAME  = "l1tTTTracksFromExtendedTrackletEmulation"
    L1TRK_LABEL = "Level1TTTracks"
    L1TRUTH_NAME = "TTTrackAssociatorFromPixelDigisExtended"

# HYBRID_NEWKF: prompt tracking or reduced
elif (L1TRKALGO == 'HYBRID_NEWKF' or L1TRKALGO == 'HYBRID_REDUCED'):
    process.load( 'L1Trigger.TrackFindingTracklet.Producer_cff' )
    process.load( 'L1Trigger.TrackFindingTracklet.Analyzer_cff' )
    NHELIXPAR = 4
    L1TRK_NAME  = process.TrackFindingTrackletAnalyzer_params.OutputLabelTFP.value()
    L1TRK_LABEL = process.TrackFindingTrackletProducer_params.BranchTTTracks.value()
    L1TRUTH_NAME = "TTTrackAssociatorFromPixelDigis"
    process.TTTrackAssociatorFromPixelDigis.TTTracks = cms.VInputTag( cms.InputTag(L1TRK_NAME, L1TRK_LABEL) )
    process.HybridNewKF = cms.Sequence(process.L1THybridTracks + process.ProducerTM + process.ProducerDR + process.ProducerKF + process.ProducerTQ + process.ProducerTFP)
    process.TTTracksEmulation = cms.Path(process.HybridNewKF)
    #process.TTTracksEmulationWithTruth = cms.Path(process.HybridNewKF +  process.TrackTriggerAssociatorTracks)
    # Optionally include code producing performance plots & end-of-job summary.
    process.load( 'SimTracker.TrackTriggerAssociation.StubAssociator_cff' )
    process.TTTracksEmulationWithTruth = cms.Path(process.HybridNewKF +  process.TrackTriggerAssociatorTracks + process.StubAssociator +  process.AnalyzerTracklet + process.AnalyzerTM + process.AnalyzerDR + process.AnalyzerKF + process.AnalyzerTQ + process.AnalyzerTFP )
    from L1Trigger.TrackFindingTracklet.Customize_cff import *
    if (L1TRKALGO == 'HYBRID_NEWKF'):
        fwConfig( process )
    if (L1TRKALGO == 'HYBRID_REDUCED'):
        reducedConfig( process )
    # Needed by L1TrackNtupleMaker
    process.HitPatternHelperSetup.useNewKF = True


process.L1TrackNtuple = L1TrackNtupleMaker.clone(
   L1Tk_nPar = NHELIXPAR, # use 4 or 5-parameter L1 tracking?
   L1TrackInputTag = (L1TRK_NAME, L1TRK_LABEL),         # TTTrack input
   MCTruthTrackInputTag = (L1TRUTH_NAME, L1TRK_LABEL),  # MCTruth input
)