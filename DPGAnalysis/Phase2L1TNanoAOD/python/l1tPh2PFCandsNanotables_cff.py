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

p2L1PFCandsTask = cms.Task(
    l1tPuppiCandsTable,
    l1tPFCandsTable,
)

p2L1ExtPFCandsTask = cms.Task(
    l1tExtPuppiCandsTable,
)
