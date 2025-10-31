import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
from Configuration.StandardSequences.Eras import eras
from PhysicsTools.NanoAOD.common_cff import Var, ExtVar
import os
import glob


options = VarParsing.VarParsing('analysis')
options.register ('emulatorMode',
                  'correctionlibRegression', # default value
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.string,
                  "SmartPixels mode: (passthrough, passthroughFloat, passthroughHW, trackingParticleTruth, correctionlibRegression)")
options.register ('activeSP',
                  '0000', # default value
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.string,
                  "Active layers as 4 digit code from innermost to outermost layer: '0000', ..., '0101', ..., '1111'")
options.register ('noGTT',
                  False, # default value
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.bool,
                  "Disable passing SmartPixels tracks to Global Track Trigger")
options.register ('noGMT',
                  False, # default value
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.bool,
                  "Disable passing SmartPixels tracks to Global Muon Trigger")
options.register ('noCLX',
                  False, # default value
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.bool,
                  "Disable passing SmartPixels tracks to Correlator Layer 1 and 2")
options.parseArguments()

# Smart Pixels configurations
if options.emulatorMode.lower() == "correctionlibregression":
    spixfourparam = f"l1tSmartPixelsTrackProducerWcorrectionlibRegression{options.activeSP.replace('0', 'I').replace('1', 'A')}"
    spixfiveparam = f"l1tSmartPixelsTrackProducerExtendedWcorrectionlibRegression{options.activeSP.replace('0', 'I').replace('1', 'A')}"
elif options.emulatorMode.lower() == "correctionlibtptoysmear":
    raise NotImplementedError("correctionlibTPToySmear not yet supported in this configuration, is it available in the l1tSmartPixelsTrackProducer?")
else:
    assert options.emulatorMode in ["passthrough", "passthroughFloat", "passthroughHW", "trackingParticleTruth"]
    spixfourparam = f"l1tSmartPixelsTrackProducerW{options.emulatorMode}"
    spixfiveparam = f"l1tSmartPixelsTrackProducerExtendedW{options.emulatorMode}"
spixpostfix = spixfourparam.replace("l1tSmartPixelsTrackProducerW", "injectedSP")

skipModulesForSP = []
if options.noGTT:
    spixpostfix += "_noGTT"
    skipModulesForSP.append("L1GTTInputProducer")
if options.noGMT:
    spixpostfix += "_noGMT"
    skipModulesForSP.append("Phase2L1TGMTTkMuonProducer")
if options.noCLX:
    spixpostfix += "_noCLX"
    skipModulesForSP.append("PFTrackProducerFromL1Tracks")
if len(skipModulesForSP) == 0:
    skipModulesForSP = None #uses defaults
else:
    # Add defaults too for explicit list of producer to skip co-opting inputs to
    skipModulesForSP += ["L1SmartPixelsTrackProducer",
                         "TTTrackAssociator_Phase2TrackerDigi_",
                         "L1FPGATrackProducer",
                         "SimpleL1TTTrackCandidateFlatTableProducer"
                         ]

process = cms.Process("RESP", eras.Phase2C17I13M9)

process.load('Configuration.StandardSequences.Services_cff')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False), allowUnscheduled = cms.untracked.bool(False) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))
process.MessageLogger.cerr.FwkReport.reportEvery = 250
inputMC = [f'root://cmseos.fnal.gov//eos/uscms/store/group/lpcsmartpixelscms/jettagging/step1/v2p1/Test/TT_TuneCP5_14TeV-powheg-pythia8/crab_jettagging_step1_v2p1_TestPathsOkay/251030_203121/0000/output_{x}.root' for x in range(1, 3)]
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(*inputMC),
    inputCommands = cms.untracked.vstring("keep *", 
            "drop l1tPFClusters_*_*_*",
            "drop l1tPFTracks_*_*_*",
            "drop l1tPFCandidates_*_*_*",
            "drop l1tTkPrimaryVertexs_*_*_*",
            "drop l1tKMTFTracks_*_*_*")
)

process.load('Configuration.Geometry.GeometryExtendedRun4D110Reco_cff')
process.load('Configuration.Geometry.GeometryExtendedRun4D110_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff') # needed to read HCal TPs
process.load('SimCalorimetry.HGCalSimProducers.hgcalDigitizer_cfi') # needed for HGCAL_noise_fC
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('RecoMET.Configuration.GenMETParticles_cff')
process.load('RecoMET.METProducers.genMetTrue_cfi')

from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
from RecoMET.METProducers.pfMet_cfi import pfMet

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '141X_mcRun4_realistic_v3', '')

# NOTE: we need this to avoid saving the stubs
process.l1tTrackSelectionProducer.processSimulatedTracks = False

from L1Trigger.L1CaloTrigger.l1tPhase2L1CaloEGammaEmulator_cfi import l1tPhase2L1CaloEGammaEmulator
process.l1tPhase2L1CaloEGammaEmulator = l1tPhase2L1CaloEGammaEmulator.clone()

process.extraPFStuff = cms.Task(
        process.l1tPhase2L1CaloEGammaEmulator,
        process.l1tSAMuonsGmt,
        process.l1tGTTInputProducer,
        process.l1tTrackSelectionProducer,
        process.l1tVertexFinderEmulator,
        process.L1TLayer1TaskInputsTask,
        process.L1TLayer1Task,
        process.L1TLayer2EGTask)


def addJetNTuple(trktype = "extended", nparam = 5, tagged = True):
    # create new jet tupler
    jetColl = "l1tSC4PFL1PuppiExtendedEmulator"
    jetCollCorr = "l1tSC4PFL1PuppiExtendedEmulator"
    if trktype == "baseline":
        jetColl = "l1tSC4PFL1PuppiEmulator"
        jetCollCorr = "l1tSC4PFL1PuppiCorrectedEmulator"
        process.l1tPFTracksFromL1Tracks.nParam = cms.uint32(nparam)
    else:
        process.l1tPFTracksFromL1TracksExtended.nParam = cms.uint32(nparam)
    if tagged:
        jetColl = ("l1tSC4NGJetProducer","l1tSC4NGJets")
        jetCollCorr = "l1tSC4PFL1PuppiCorrectedEmulator"

    process.outnano = cms.EDAnalyzer("JetNTuplizer",
        genJets = cms.InputTag("ak4GenJetsNoNu"),
        genParticles = cms.InputTag("genParticles"),
        scPuppiJets = cms.InputTag(jetColl),
        scPuppiJetsCorr = cms.InputTag(jetCollCorr),
        nnTaus = cms.InputTag("l1tNNTauProducerPuppi","L1PFTausNN"),
        genJetsFlavour = cms.InputTag("genFlavourInfo"),
        vtx = cms.InputTag("l1tVertexFinderEmulator","L1VerticesEmulation"),
        bjetIDs = cms.InputTag("l1tBJetProducerPuppiCorrectedEmulator", "L1PFBJets"),
        electrons = cms.InputTag("l1tLayer2EG","L1CtTkElectron"),
        muons = cms.InputTag("l1tSAMuonsGmt","promptSAMuons"),
    )
    process.endTuple = cms.EndPath(process.outnano)
    outName = "jetTuple_v2p1_"+trktype+"_"+str(nparam)+"_"+spixpostfix + ".root"
    process.TFileService = cms.Service("TFileService", fileName = cms.string(outName))

# to check available tags:
process.p = cms.Path()
process.p.associate(process.extraPFStuff)
process.p.associate(process.L1TPFJetsExtendedTask)
process.p.associate(process.L1TBJetsTask)
#process.p.associate(process.l1tSC4NGJetTask)
process.TFileService = cms.Service("TFileService", fileName = cms.string("jetTuple.root"))

def addNNPuppiTaus():
    process.load("L1Trigger.Phase2L1ParticleFlow.L1NNTauProducer_cff")
    process.l1tNNTauProducerPuppi.maxtaus = cms.int32(500)
    process.extraPFStuff.add(process.l1tNNTauProducerPuppi)

def addSeededConeJets():
    process.extraPFStuff.add(process.L1TPFJetsTask)
    process.extraPFStuff.add(process.L1TPFJetsExtendedTask)

def addMultitagging(trktype = "extended"):
    if trktype == "extended":
        process.l1tSC4NGJetProducer.jets = cms.InputTag("l1tSC4PFL1PuppiExtendedEmulator")
    else:
        process.l1tSC4NGJetProducer.jets = cms.InputTag("l1tSC4PFL1PuppiEmulator")
    process.l1tSC4NGJetProducer.maxJets = cms.int32(500)
    l1tSC4NGJetModelPath_so = "L1TSC4NGJetModel/L1TSC4NGJetModel_v0/L1TSC4NGJetModel_v0.so"
    if not os.path.isfile(l1tSC4NGJetModelPath_so):
        os.system('tar -xzf L1TSC4NGJetModels.tar.gz')
    process.l1tSC4NGJetProducer.l1tSC4NGJetModelPath = cms.string(l1tSC4NGJetModelPath_so.replace(".so", ""))
    process.extraPFStuff.add(process.l1tSC4NGJetProducer)

def addBtagging(jetColl): #extended TRK
    process.load("L1Trigger.Phase2L1ParticleFlow.L1BJetProducer_cff")
    process.l1tBJetProducerPuppiCorrectedEmulator.jets = cms.InputTag(jetColl)
    process.l1tBJetProducerPuppiCorrectedEmulator.maxJets = cms.int32(500)
    process.l1tBJetProducerPuppiCorrectedEmulator.useRawPt = cms.bool(True)
    process.extraPFStuff.add(process.L1TBJetsTask)
    #process.l1pfjetTable.jets.scPuppiBJet = cms.InputTag('l1tBJetProducerPuppiCorrectedEmulator')  

def addGenJetFlavourTable():
    process.load("PhysicsTools.JetMCAlgos.AK4PFJetsMCFlavourInfos_cfi")
    process.load("PhysicsTools.JetMCAlgos.HadronAndPartonSelector_cfi")
    process.selectedHadronsAndPartons.partonMode = cms.string("Pythia8")
    process.genFlavourInfo = process.ak4JetFlavourInfos.clone(jets = "ak4GenJetsNoNu")
    process.p += process.selectedHadronsAndPartons
    process.p += process.genFlavourInfo

def goMT(nthreads=1):
    process.options.numberOfThreads = cms.untracked.uint32(nthreads)
    process.options.numberOfStreams = cms.untracked.uint32(0)

if True:
    process.source.fileNames  = cms.untracked.vstring(*inputMC)
    goMT()
    trktype = "extended"
    nparam = 5
    addSeededConeJets()
    addMultitagging(trktype = trktype)
    addBtagging(("l1tSC4NGJetProducer","l1tSC4NGJets"))
    addNNPuppiTaus()
    addGenJetFlavourTable()
    addJetNTuple(trktype = trktype, nparam = nparam)
    if False:
        open("debug_dump_runJetNTuple.py", "w").write(process.dumpPython())

# from FastPUPPI.NtupleProducer.SmartPixelsCustomize_cff import referenceSmartPixelsTrackProducer
# process = referenceSmartPixelsTrackProducer(process)
from L1Trigger.TrackFindingTracklet.Customize_cff import injectSmartPixelsTrackProducer
process = injectSmartPixelsTrackProducer(process,
                                         # Only matters if addAssociation is on, which creates new SmartPixels producers; we only reference embedded collections here
                                         smartPixelsEmulatorMode=None,
                                         addAssociation=False,
                                         l1tSmartPixelsTrackProducerLabel=spixfourparam,
                                         l1tSmartPixelsTrackProducerExtendedLabel=spixfiveparam,
                                         skipModuleTypes=skipModulesForSP,
                                         printProcessInfo=False)
process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
