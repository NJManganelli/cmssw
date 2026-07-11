import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import Var

#### L1T tracks from the tracklet emulation, with float and hardware (track word) values
l1tTracksTable = cms.EDProducer(
    "SimpleL1TTTrackCandidateFlatTableProducer",
    src = cms.InputTag("l1tTTTracksFromTrackletEmulation", "Level1TTTracks"),
    name = cms.string("L1TTrack"),
    doc = cms.string("L1T Prompt Tracks"),
    cut = cms.string(""),
    singleton = cms.bool(False), # the number of entries is variable
    variables = cms.PSet(
        rInv = Var("rInv()", float, doc="rInv"),
        pt = Var("momentum().perp()", float, doc="pt"),
        eta = Var("eta()", float, doc="eta"),
        phi = Var("phi()", float, doc="phi"),
        localPhi = Var("localPhi()", float, doc="local phi"),
        tanL = Var("tanL()", float, doc="tanL"),
        z0 = Var("z0()", float, doc="z0"),
        d0 = Var("d0()", float, doc="d0"),
        trkMVA1 = Var("trkMVA1()", float, doc="track MVA1"),
        trkMVA2 = Var("trkMVA2()", float, doc="track MVA2"),
        trkMVA3 = Var("trkMVA3()", float, doc="track MVA3"),
        hitPattern = Var("hitPattern()", "uint", doc="hit pattern"),
        chi2XYRed = Var("chi2XYRed()", float, doc="chi2XY per dof"),
        chi2ZRed = Var("chi2ZRed()", float, doc="chi2Z per dof"),
        chi2BendRed = Var("chi2BendRed()", float, doc="stub bend chi2 per stub (formerly stubPtConsistency)"),
        chi2Bend = Var("chi2Bend()", float, doc="stub bend chi2"),
        phiSector = Var("phiSector()", "uint", doc = "phi sector number"),
        etaSector = Var("etaSector()", "uint", doc = "eta sector number"),
        hwValid = Var("getValidBits()", "uint", doc = "hardware valid bit"),
        hwRinv = Var("getRinvBits()", "uint", doc = "hardware 1/R before GTT conversion"),
        hwTanl = Var("getTanlBits()", "uint", doc = "hardware tanL before GTT conversion"),
        hwPhi = Var("getPhiBits()", "uint", doc = "hardware phi"),
        hwZ0 = Var("getZ0Bits()", "uint", doc = "hardware z0"),
        hwD0 = Var("getD0Bits()", "uint", doc = "hardware d0"),
        hwChi2RPhi = Var("getChi2RPhiBits()", "uint", doc = "hardware Chi2RPhi"),
        hwChi2RZ = Var("getChi2RZBits()", "uint", doc = "hardware Chi2RZ"),
        hwBendChi2 = Var("getBendChi2Bits()", "uint", doc = "hardware BendChi2"),
        hwHitPattern = Var("getHitPatternBits()", "uint", doc = "hardware HitPattern"),
        hwMVAQuality = Var("getMVAQualityBits()", "uint", doc = "hardware MVA Quality, corresponding to trkMVA1"),
        hwMVAOther = Var("getMVAOtherBits()", "uint", doc = "hardware MVA Other, reserved for trkMVA2 and trkMVA3"),
     )
)

#### hardware values after GTT input conversion, as extension tables of the main track tables
gttTracksTable = l1tTracksTable.clone(
    src = cms.InputTag('l1tGTTInputProducer','Level1TTTracksConverted'),
    extension = cms.bool(True),
    variables = cms.PSet(
        hwPt = Var("getRinvBits()", "uint", doc = "hardware pt after GTT conversion"),
        hwEta = Var("getTanlBits()", "uint", doc = "hardware eta after GTT conversion"),
        )
)

l1tExtTracksTable = l1tTracksTable.clone(
    src = cms.InputTag("l1tTTTracksFromExtendedTrackletEmulation", "Level1TTTracks"),
    name = cms.string("L1TExtTrack"),
    doc = cms.string("L1T Extended Tracks"),
)

gttExtTracksTable = l1tExtTracksTable.clone(
    src = cms.InputTag('l1tGTTInputProducerExtended','Level1TTTracksExtendedConverted'),
    extension = cms.bool(True),
    variables = cms.PSet(
        hwPt = Var("getRinvBits()", "uint", doc = "hardware pt after GTT conversion"),
        hwEta = Var("getTanlBits()", "uint", doc = "hardware eta after GTT conversion"),
        )
)

#### SmartPixels-augmented track collections
#### (require the L1SmartPixelsTrackProducer from the smartpixels development branch)
l1tSmartPixelsTracksTable = l1tTracksTable.clone(
    src = cms.InputTag("l1tSmartPixelsTrackProducer", "Level1TTTracks"),
    name = cms.string("L1TSmartPixelsTrack"),
    doc = cms.string("L1T SmartPixels Tracks"),
)

l1tExtSmartPixelsTracksTable = l1tTracksTable.clone(
    src = cms.InputTag("l1tSmartPixelsTrackProducerExtended", "Level1TTTracks"),
    name = cms.string("L1TSmartPixelsExtTrack"),
    doc = cms.string("L1T SmartPixels Extended Tracks"),
)

#### MC truth extension (nStubs always; genuine/fake + TrackingParticle info
#### when the TTTrackAssociator ran in the workflow) -- inputs for the
#### track-quality GBDT training
l1tTracksTruthTable = cms.EDProducer(
    "L1TrackTruthTableProducer",
    tracks = cms.InputTag("l1tTTTracksFromTrackletEmulation", "Level1TTTracks"),
    truthAssociation = cms.InputTag("TTTrackAssociatorFromPixelDigis", "Level1TTTracks"),
    trackTableName = cms.string("L1TTrack"),
)

l1tExtTracksTruthTable = l1tTracksTruthTable.clone(
    tracks = cms.InputTag("l1tTTTracksFromExtendedTrackletEmulation", "Level1TTTracks"),
    truthAssociation = cms.InputTag("TTTrackAssociatorFromPixelDigisExtended", "Level1TTTracks"),
    trackTableName = cms.string("L1TExtTrack"),
)

p2L1TrackTruthTask = cms.Task(
    l1tTracksTruthTable,
    l1tExtTracksTruthTable,
)

p2L1TracksTask = cms.Task(
    l1tTracksTable,
    gttTracksTable,
    l1tExtTracksTable,
    gttExtTracksTable,
)

p2L1SmartPixelsTracksTask = cms.Task(
    l1tSmartPixelsTracksTable,
    l1tExtSmartPixelsTracksTable,
)
