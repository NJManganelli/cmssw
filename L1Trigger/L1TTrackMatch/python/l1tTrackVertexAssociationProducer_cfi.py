import FWCore.ParameterSet.Config as cms

l1tTrackVertexAssociationProducer = cms.EDProducer('L1TrackVertexAssociationProducer',
  l1SelectedTracksInputTag = cms.InputTag("l1tTrackSelectionProducer", "Level1TTTracksSelected"),
  l1SelectedTracksEmulationInputTag = cms.InputTag("l1tTrackSelectionProducer", "Level1TTTracksSelectedEmulation"),
  # If no vertex collection is provided, then the DeltaZ cuts will not be run
  l1VerticesInputTag = cms.InputTag("l1tVertexFinderSim", "l1vertices"),
  l1VerticesEmulationInputTag = cms.InputTag("l1tVertexFinderEmu", "l1verticesEmulation"),
  outputCollectionName = cms.string("Level1TTTracksSelectedAssociated"),
  cutSet = cms.PSet(
                    #deltaZMaxEtaBounds = cms.vdouble(0.0, absEtaMax.value), # these values define the bin boundaries in |eta|
                    #deltaZMax = cms.vdouble(0.5), # delta z must be less than these values, there will be one less value here than in deltaZMaxEtaBounds, [cm]
                    deltaZMaxEtaBounds = cms.vdouble(0.0, 0.7, 1.0, 1.2, 1.6, 2.0, 2.4), # these values define the bin boundaries in |eta|
                    deltaZMax = cms.vdouble(0.37, 0.50, 0.60, 0.75, 1.00, 1.60), # delta z must be less than these values, there will be one less value here than in deltaZMaxEtaBounds, [cm]
                    ),
  useDisplacedTracksDeltaZOverride = cms.double(-1.0), # override the deltaZ cut value for displaced tracks
  processSimulatedTracks = cms.bool(True), # return selected tracks after cutting on the floating point values
  processEmulatedTracks = cms.bool(True), # return selected tracks after cutting on the bitwise emulated values
  debug = cms.int32(0) # Verbosity levels: 0, 1, 2, 3, 4
)

L1TrackVertexAssociationProducerExtended = l1tTrackVertexAssociationProducer.clone(
  l1SelectedTracksInputTag = cms.InputTag("l1tTrackSelectionProducer", "Level1TTTracksSelected"),
  l1SelectedTracksEmulationInputTag = cms.InputTag("l1tTrackSelectionProducer", "Level1TTTracksSelectedEmulation"),
  outputCollectionName = cms.string("Level1TTTracksExtendedSelectedAssociated"),
  cutSet = cms.PSet(
                    #deltaZMaxEtaBounds = cms.vdouble(0.0, absEtaMax.value), # these values define the bin boundaries in |eta|
                    #deltaZMax = cms.vdouble(0.5), # delta z must be less than these values, there will be one less value here than in deltaZMaxEtaBounds, [cm]
                    deltaZMaxEtaBounds = cms.vdouble(0.0, 0.7, 1.0, 1.2, 1.6, 2.0, 2.4), # these values define the bin boundaries in |eta|
                    deltaZMax = cms.vdouble(3.0, 3.0, 3.0, 3.0, 3.0, 3.0), # delta z must be less than these values, there will be one less value here than in deltaZMaxEtaBounds, [cm]
                    ),
  useDisplacedTracksDeltaZOverride = cms.double(3.0), # Use promt/displaced tracks
  processSimulatedTracks = cms.bool(True), # return selected tracks after cutting on the floating point values
  processEmulatedTracks = cms.bool(True), # return selected tracks after cutting on the bitwise emulated values
  debug = cms.int32(0) # Verbosity levels: 0, 1, 2, 3, 4
)


