import FWCore.ParameterSet.Config as cms

L1TrackSelectionProducer = cms.EDProducer('L1TrackSelectionProducer',
  l1TracksInputTag = cms.InputTag("L1GTTInputProducer","Level1TTTracksConverted"),
  # If no vertex collection is provided, then the DeltaZ cuts will not be run
  outputCollectionName = cms.string("Level1TTTracksSelected"),
  cutSet = cms.PSet(
                    ptMin = cms.double(2.0), # pt must be greater than this value, [GeV]
                    absEtaMax = cms.double(2.4), # absolute value of eta must be less than this value
                    absZ0Max = cms.double(15.0), # z0 must be less than this value, [cm]
                    nStubsMin = cms.int32(4), # number of stubs must be greater than or equal to this value
                    nPSStubsMin = cms.int32(0), # the number of stubs in the PS Modules must be greater than or equal to this value

                    reducedBendChi2Max = cms.double(2.25), # bend chi2 must be less than this value
                    reducedChi2RZMax = cms.double(5.0), # chi2rz/dof must be less than this value
                    reducedChi2RPhiMax = cms.double(20.0), # chi2rphi/dof must be less than this value
                    ),
  processSimulatedTracks = cms.bool(True), # return selected tracks after cutting on the floating point values
  processEmulatedTracks = cms.bool(True), # return selected tracks after cutting on the bitwise emulated values
  debug = cms.int32(0) # Verbosity levels: 0, 1, 2, 3, 4
)

L1TrackSelectionProducerExtended = L1TrackSelectionProducer.clone(
  l1TracksInputTag = cms.InputTag("L1GTTInputProducerExtended","Level1TTTracksExtendedConverted"),
  outputCollectionName = cms.string("Level1TTTracksExtendedSelected"),
  cutSet = cms.PSet(
                    ptMin = cms.double(3.0), # pt must be greater than this value, [GeV]
                    absEtaMax = cms.double(2.4), # absolute value of eta must be less than this value
                    absZ0Max = cms.double(15.0), # z0 must be less than this value, [cm]
                    nStubsMin = cms.int32(4), # number of stubs must be greater than or equal to this value
                    nPSStubsMin = cms.int32(0), # the number of stubs in the PS Modules must be greater than or equal to this value

                    reducedBendChi2Max = cms.double(2.4), # bend chi2 must be less than this value
                    reducedChi2RZMax = cms.double(10.0), # chi2rz/dof must be less than this value
                    reducedChi2RPhiMax = cms.double(40.0), # chi2rphi/dof must be less than this value
                    ),
)


