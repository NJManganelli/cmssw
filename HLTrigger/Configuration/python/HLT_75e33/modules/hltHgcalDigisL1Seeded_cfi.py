import FWCore.ParameterSet.Config as cms

hltHgcalDigisL1Seeded = cms.EDProducer("HLTHGCalDigisInRegionsProducer",
    etaPhiRegions = cms.VPSet(cms.PSet(
        inputColl = cms.InputTag("hltL1TEGammaHGCFilteredCollectionProducer"),
        maxDEta = cms.double(0.0),
        maxDPhi = cms.double(0.0),
        maxDeltaR = cms.double(0.35),
        maxEt = cms.double(999999.0),
        minEt = cms.double(5.0),
        type = cms.string('L1P2GTCandidate')
    )),
    inputCollTags = cms.VInputTag("hltHgcalDigis:EE", "hltHgcalDigis:HEback", "hltHgcalDigis:HEfront"),
    outputProductNames = cms.vstring(
        'EE',
        'HEback',
        'HEfront'
    )
)
