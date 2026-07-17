import FWCore.ParameterSet.Config as cms

# Phase-1 payload-derivation analyzer (Tier-2 digiRefit plan).
# Pure file-reading: consumes persisted TTTracks, the TTTrack association map,
# TrackingParticles, the simSiPixelDigis:Pixel PixelDigi + PixelDigiSimLink
# collections, and g4SimHits SimTracks/SimVertices. Records, per L1Track x
# TBPX-layer crossing, the in-window pixel digis with per-digi truth class and
# module-LOCAL angle/B-field (PixelAV cotAlpha/cotBeta convention).
smartPixelsPayloadAnalyzer = cms.EDAnalyzer(
    'SmartPixelsPayloadAnalyzer',
    l1TracksInputTag           = cms.InputTag("l1tTTTracksFromTrackletEmulation", "Level1TTTracks"),
    mcTruthTrackInputTag       = cms.InputTag("TTTrackAssociatorFromPixelDigis", "Level1TTTracks"),
    trackingParticleInputTag   = cms.InputTag("mix", "MergedTrackTruth"),
    pixelDigiInputTag          = cms.InputTag("simSiPixelDigis", "Pixel"),
    pixelDigiSimLinkInputTag   = cms.InputTag("simSiPixelDigis", "Pixel"),
    simTrackInputTag           = cms.InputTag("g4SimHits"),
    simVertexInputTag          = cms.InputTag("g4SimHits"),
    # MEASUREMENT window half-widths [cm]. Deliberately much wider than the
    # refit windows in DIGIREFIT_DEFAULTS: the OT-track expected-position
    # resolution at TBPX radii is O(100um) in r-phi but O(mm) in z (L1 z0
    # resolution, no material model), and the payload must measure the FULL
    # residual distribution — truncation belongs to the refit, not here.
    windowRPhi                 = cms.double(0.30),
    windowZ                    = cms.double(0.50),
    nLayers                    = cms.int32(4),   # TBPX layers 1..4
    debug                      = cms.bool(False),
)
