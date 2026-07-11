import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import Var

#### L1 PF / Puppi candidate tables
#### The variable set is a superset of the per-candidate inputs consumed by the
#### L1TSC4NGJetID jet tagger and written by the FastPUPPI jet ntupler
#### (jet_pfcand_*): pt, eta, phi, mass, charge, id, z0, dxy, puppiWeight,
#### emID, quality plus the hardware words. Jet-relative quantities
#### (pt_rel, deta, dphi, ...) are recomputed offline from the candidate and
#### jet tables.
l1tPuppiCandsTable = cms.EDProducer(
    "SimpleTriggerL1PFCandidateFlatTableProducer",
    src = cms.InputTag("l1tLayer2Deregionizer", "Puppi"),
    name = cms.string("L1PuppiCand"),
    doc = cms.string("Deregionized Puppi candidates (input to the SC4 jets and the NGJet tagger)"),
    cut = cms.string(""),
    singleton = cms.bool(False), # the number of entries is variable
    variables = cms.PSet(
        pt = Var("pt()", float, doc="pt"),
        eta = Var("eta()", float, doc="eta"),
        phi = Var("phi()", float, doc="phi"),
        mass = Var("mass()", float, doc="mass"),
        charge = Var("charge()", "int16", doc="charge"),
        pdgId = Var("pdgId()", "int", doc="PDG id"),
        id = Var("id()", "uint8", doc="particle type: ChargedHadron=0, Electron=1, NeutralHadron=2, Photon=3, Muon=4"),
        z0 = Var("z0()", float, doc="z0 (vz)"),
        dxy = Var("dxy()", float, doc="dxy"),
        caloEta = Var("caloEta()", float, doc="eta at calorimeter surface"),
        caloPhi = Var("caloPhi()", float, doc="phi at calorimeter surface"),
        puppiWeight = Var("puppiWeight()", float, doc="puppi weight"),
        hwPt = Var("hwPt()", "int", doc="hardware pt"),
        hwEta = Var("hwEta()", "int", doc="hardware eta"),
        hwPhi = Var("hwPhi()", "int", doc="hardware phi"),
        hwQual = Var("hwQual()", "int", doc="hardware quality / particle type"),
        hwZ0 = Var("hwZ0()", "int16", doc="hardware z0"),
        hwDxy = Var("hwDxy()", "int16", doc="hardware dxy"),
        hwTkQuality = Var("hwTkQuality()", "uint16", doc="hardware track quality"),
        hwPuppiWeight = Var("hwPuppiWeight()", "uint16", doc="hardware puppi weight"),
        hwEmID = Var("hwEmID()", "uint16", doc="hardware EM id"),
        encodedPuppi64 = Var("encodedPuppi64()", "uint64", doc="64-bit hardware puppi word"),
     )
)

l1tExtPuppiCandsTable = l1tPuppiCandsTable.clone(
    src = cms.InputTag("l1tLayer2DeregionizerExtended", "Puppi"),
    name = cms.string("L1ExtPuppiCand"),
    doc = cms.string("Deregionized Puppi candidates from extended (displaced) tracking"),
)

l1tPFCandsTable = l1tPuppiCandsTable.clone(
    src = cms.InputTag("l1tLayer1", "PF"),
    name = cms.string("L1PFCand"),
    doc = cms.string("Layer-1 particle flow candidates (pre-Puppi)"),
)

#### Jet <-> constituent association tables (PFNano JetPFCands style), with
#### candidate-table extension columns (jetIdx backref, l1TrackIdx crossref)
l1tSC4NGJetCandsTable = cms.EDProducer(
    "L1JetCandLinkTableProducer",
    jets = cms.InputTag("l1tSC4NGJetProducer", "l1tSC4NGJets"),
    cands = cms.InputTag("l1tLayer2DeregionizerExtended", "Puppi"),
    tracks = cms.InputTag("l1tTTTracksFromExtendedTrackletEmulation", "Level1TTTracks"),
    linkTableName = cms.string("L1SC4NGJetCands"),
    candTableName = cms.string("L1ExtPuppiCand"),
    jetTableName = cms.string("L1puppiJetSC4NG"),
    trackTableName = cms.string("L1TExtTrack"),
    maxConstituentsTagger = cms.uint32(16),
    writeCandExtension = cms.bool(True),
)

l1tSC4JetCandsTable = l1tSC4NGJetCandsTable.clone(
    jets = cms.InputTag("l1tSC4PFL1PuppiCorrectedEmulator"),
    cands = cms.InputTag("l1tLayer2Deregionizer", "Puppi"),
    tracks = cms.InputTag("l1tTTTracksFromTrackletEmulation", "Level1TTTracks"),
    linkTableName = cms.string("L1SC4JetCands"),
    candTableName = cms.string("L1PuppiCand"),
    jetTableName = cms.string("L1puppiJetSC4"),
    trackTableName = cms.string("L1TTrack"),
)

#### HGCal multicluster table (placeholder for the jet ntupler cluster_*
#### inputs). Covers hOverE / sigmaRR / zBarycenter; the egVsPion / egVsPU
#### MVA outputs are not stored on any central dataformat (they were
#### fork-only PFCluster extensions) and need HGCal TPG support upstream.
l1tHGCClusterTable = cms.EDProducer(
    "SimpleTriggerL1HGCalMulticlusterFlatTableProducer",
    src = cms.InputTag("l1tHGCalBackEndLayer2Producer", "HGCalBackendLayer2Processor3DClustering"),
    name = cms.string("L1HGCCluster"),
    doc = cms.string("HGCal backend multiclusters (endcap calo input to correlator layer-1)"),
    cut = cms.string(""),
    singleton = cms.bool(False), # the number of entries is variable
    variables = cms.PSet(
        pt = Var("pt()", float, doc="pt"),
        eta = Var("eta()", float, doc="eta"),
        phi = Var("phi()", float, doc="phi"),
        hOverE = Var("hOverE()", float, doc="hadronic over electromagnetic energy"),
        sigmaRRTot = Var("sigmaRRTot()", float, doc="total sigmaRR shower shape"),
        sigmaRRMean = Var("sigmaRRMean()", float, doc="mean sigmaRR shower shape"),
        sigmaZZ = Var("sigmaZZ()", float, doc="sigmaZZ shower shape"),
        sigmaEtaEtaTot = Var("sigmaEtaEtaTot()", float, doc="total sigmaEtaEta shower shape"),
        sigmaPhiPhiTot = Var("sigmaPhiPhiTot()", float, doc="total sigmaPhiPhi shower shape"),
        zBarycenter = Var("zBarycenter()", float, doc="z of the shower barycenter (signed; abs for |z|)"),
        eMax = Var("eMax()", float, doc="energy of the most energetic layer"),
        showerLength = Var("showerLength()", "int16", doc="shower length in layers"),
        coreShowerLength = Var("coreShowerLength()", "int16", doc="core shower length in layers"),
        firstLayer = Var("firstLayer()", "int16", doc="first layer of the shower"),
        maxLayer = Var("maxLayer()", "int16", doc="layer with maximum energy"),
        hwQual = Var("hwQual()", "int", doc="hardware quality / ID bits"),
    )
)

#### caloPtr -> HGCal multicluster crossref (endcap candidates only; barrel
#### candidates point at GCT/PFCluster collections and get -1)
l1tPuppiCandHGCClusterLink = cms.EDProducer(
    "L1PFCandClusterLinkTableProducer",
    cands = cms.InputTag("l1tLayer2Deregionizer", "Puppi"),
    clusters = cms.InputTag("l1tHGCalBackEndLayer2Producer", "HGCalBackendLayer2Processor3DClustering"),
    candTableName = cms.string("L1PuppiCand"),
    columnName = cms.string("hgcClusterIdx"),
    clusterTableName = cms.string("L1HGCCluster"),
)

l1tExtPuppiCandHGCClusterLink = l1tPuppiCandHGCClusterLink.clone(
    cands = cms.InputTag("l1tLayer2DeregionizerExtended", "Puppi"),
    candTableName = cms.string("L1ExtPuppiCand"),
)

p2L1PFCandsTask = cms.Task(
    l1tPuppiCandsTable,
    l1tExtPuppiCandsTable,
    l1tPFCandsTable,
    l1tSC4NGJetCandsTable,
    l1tSC4JetCandsTable,
    l1tHGCClusterTable,
    l1tPuppiCandHGCClusterLink,
    l1tExtPuppiCandHGCClusterLink,
)
