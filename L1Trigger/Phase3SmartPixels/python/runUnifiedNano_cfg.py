"""Unified SmartPixels NanoAOD production pset (CRAB-driven).

Produces, in ONE cmsRun job on a GEN-SIM-DIGI-RAW-MINIAOD input:
  (1) the L1 nano tier @L1PFTrkNanowithGen (fullest autoNANO flavor:
      L1TTrack + truth, L1PuppiCand/PFCand crossrefs, SC4/SC8 + NG jet tier,
      HGCCluster, GenJet parton/hadron flavour, ...), re-emulating
      L1TrackTrigger + SimL1Emulator (posture B) from the input, AND
  (2) the SmartPixels digiRefit track + refit-hit + stub tables for the chosen
      activeSP, layered via smartPixelsCoexist (truthSource='inJob' -- the
      re-emulated tracklet tracks carry a real helixCovMat, so seedCovMode
      'trackCov' is valid here; the covariance is native to this branch via
      the cms-L1TK merge, no covMatrix backport needed).

Input must be POST-cms-sw/cmssw#51503: that PR rewrote the persisted L1
track-trigger DataFormats, so pre-break samples (CMSSW_14/15-era Phase2Spring24
DIGI-RAW, the input the CMSSW_17_0 version of this pset targeted) are NOT valid
here. Local validation uses the D121 RelVals; the CRAB generator ships with an
empty default dataset because no post-break dataset path is recorded yet.

Both land in the SAME NanoAOD Events tree.

This pset is driven by CRAB (a generator script emits crab configs pointing at
it). CRAB (pluginName='Analysis') OVERRIDES process.source.fileNames from the
dataset, so the PoolSource below is a valid placeholder taken from
options.inputFiles.

Local test (D121 RelVal from the SD mount):
  cmsRun runUnifiedNano_cfg.py maxEvents=2 activeSP=1100 seedCovMode=trackCov \
    pixelavAngleSet=/work/spxsmoke/spx_angle_response_Conv1D_Full-2bit_v4fixed.json \
    inputFiles=file:/host_volumes/NJM256GBSD/smartpixels-cmssw-testfiles/RelValTTbar_14TeV_PU_150X_mcRun4_realistic_v1_STD_D121_RegeneratedGS_PU-v1_file1.root \
    outputFile=/work/spx_unifiednano_test.root
"""
import os

import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from Configuration.Eras.Era_Phase2C22I13M9_cff import Phase2C22I13M9

# ---------------------------------------------------------------------------
# Parameter contract (the CRAB generator passes these via JobType.pyCfgParams)
# ---------------------------------------------------------------------------
options = VarParsing('analysis')
options.register('activeSP', '1100', VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "SmartPixels digiRefit active-layer encoding")
options.register('tier', 'L1PFTrkNanowithGen', VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "autoNANO flavor (fullest); documented, wiring is explicit below")
options.register('seedCovMode', 'trackCov', VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "digiRefit seed covariance: trackCov | parametrized")
options.register('pixelavAngleSet', 'spx_angle_response.json',
                 VarParsing.multiplicity.singleton, VarParsing.varType.string,
                 "PixelAV angle-response payload (see resolution below)")
options.setDefault('maxEvents', -1)
options.setDefault('outputFile', 'file:spx_unifiednano.root')
options.setDefault('inputFiles',
                   'file:/host_volumes/NJM256GBSD/smartpixels-cmssw-testfiles/RelValTTbar_14TeV_PU_150X_mcRun4_realistic_v1_STD_D121_RegeneratedGS_PU-v1_file1.root')
options.parseArguments()


def _resolvePixelavAngleSet(payload):
  """CRAB payload resolution: if the given path exists use it; else try
  edm.FileInPath; else look in the job cwd (CRAB ships JobType.inputFiles
  there). For local testing pass an absolute path."""
  if os.path.isabs(payload) and os.path.exists(payload):
    return payload
  if os.path.exists(payload):
    return os.path.abspath(payload)
  try:
    from FWCore.ParameterSet.Types import FileInPath
    return FileInPath(payload).fullPath()
  except Exception:
    pass
  cwd_candidate = os.path.join(os.getcwd(), os.path.basename(payload))
  if os.path.exists(cwd_candidate):
    return cwd_candidate
  # Fall through with the basename in cwd -- the producer validates loudly at
  # ctor if it cannot be read (CRAB ships it via JobType.inputFiles).
  return os.path.basename(payload)


_pixelavAngleSet = _resolvePixelavAngleSet(options.pixelavAngleSet)
print("runUnifiedNano: resolved pixelavAngleSet -> %s" % _pixelavAngleSet)

process = cms.Process('NANO', Phase2C22I13M9)

# ---------------------------------------------------------------------------
# Standard configurations (unified-nano base, from the proven spx_int_fatflavor_cfg)
# ---------------------------------------------------------------------------
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtendedRun4D121Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.L1TrackTrigger_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.SimPhase2L1GlobalTriggerEmulator_cff')
process.load('L1Trigger.Configuration.Phase2GTMenus.SeedDefinitions.step1_2024.l1tGTMenu_cff')
process.load('DPGAnalysis.Phase2L1TNanoAOD.l1tPh2Nano_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# re-emul fix (a): Hcal trigger-tower geometry / calo TP ES for standalone L1 re-emulation
process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff')

# re-emul fix (b): OT PixelDigiSimLink lives at simSiPixelDigis:Tracker in the file;
# the current-process 'mix' shadows mix:Tracker, so repoint the cluster associator.
for _m in ["TTClusterAssociatorFromPixelDigis"]:
  if hasattr(process, _m):
    getattr(process, _m).digiSimLinks = cms.InputTag("simSiPixelDigis", "Tracker")

process.maxEvents = cms.untracked.PSet(
    input=cms.untracked.int32(options.maxEvents),
    output=cms.optional.untracked.allowed(cms.int32, cms.PSet),
)

# Input source (CRAB pluginName='Analysis' OVERRIDES fileNames from the dataset)
process.source = cms.Source("PoolSource",
    fileNames=cms.untracked.vstring(options.inputFiles),
    secondaryFileNames=cms.untracked.vstring(),
)
process.options = cms.untracked.PSet(
    numberOfThreads=cms.untracked.uint32(1),
    numberOfStreams=cms.untracked.uint32(0),
    wantSummary=cms.untracked.bool(False),
)

process.configurationMetadata = cms.untracked.PSet(
    annotation=cms.untracked.string('runUnifiedNano (unified @%s + digiRefit %s)'
                                    % (options.tier, options.activeSP)),
    name=cms.untracked.string('Applications'),
    version=cms.untracked.string('$Revision: 1.0 $'),
)

# ---------------------------------------------------------------------------
# NANOAOD output (one Events tree; carries both unified and digiRefit tables)
# ---------------------------------------------------------------------------
process.NANOAODoutput = cms.OutputModule("NanoAODOutputModule",
    compressionAlgorithm=cms.untracked.string('LZMA'),
    compressionLevel=cms.untracked.int32(9),
    dataset=cms.untracked.PSet(
        dataTier=cms.untracked.string('NANOAOD'),
        filterName=cms.untracked.string(''),
    ),
    fileName=cms.untracked.string(options.outputFile),
    outputCommands=process.NANOAODEventContent.outputCommands,
)

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T35', '')

# ---------------------------------------------------------------------------
# Paths / EndPaths (unified-nano re-emulation + GT menu + nano sequence)
# ---------------------------------------------------------------------------
process.L1TrackTrigger_step = cms.Path(process.L1TrackTrigger)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.Phase2L1GTProducer = cms.Path(process.l1tGTProducerSequence)
process.Phase2L1GTAlgoBlockProducer = cms.Path(process.l1tGTAlgoBlockProducerSequence)
# The GT AlgoBlock producer consumes the HLTPathStatus of every seed path, so the
# full GT menu seed-path set must be scheduled (as in the proven unified base config).
process.pDoubleEGEle37_24 = cms.Path(process.DoubleEGEle3724)
process.pDoubleIsoTkPho22_12 = cms.Path(process.DoubleIsoTkPho2212)
process.pDoublePuppiJet112_112 = cms.Path(process.DoublePuppiJet112112)
process.pDoublePuppiJet160_35_mass620 = cms.Path(process.DoublePuppiJet16035Mass620)
process.pDoublePuppiTau52_52 = cms.Path(process.DoublePuppiTau5252)
process.pDoubleTkEle25_12 = cms.Path(process.DoubleTkEle2512)
process.pDoubleTkElePuppiHT_8_8_390 = cms.Path(process.DoubleTkElePuppiHT)
process.pDoubleTkMuPuppiHT_3_3_300 = cms.Path(process.DoubleTkMuPuppiHT)
process.pDoubleTkMuPuppiJetPuppiMet_3_3_60_130 = cms.Path(process.DoubleTkMuPuppiJetPuppiMet)
process.pDoubleTkMuon15_7 = cms.Path(process.DoubleTkMuon157)
process.pDoubleTkMuonTkEle5_5_9 = cms.Path(process.DoubleTkMuonTkEle559)
process.pDoubleTkMuon_4_4_OS_Dr1p2 = cms.Path(process.DoubleTkMuon44OSDr1p2)
process.pDoubleTkMuon_4p5_4p5_OS_Er2_Mass7to18 = cms.Path(process.DoubleTkMuon4p5OSEr2Mass7to18)
process.pDoubleTkMuon_OS_Er1p5_Dr1p4 = cms.Path(process.DoubleTkMuonOSEr1p5Dr1p4)
process.pIsoTkEleEGEle22_12 = cms.Path(process.IsoTkEleEGEle2212)
process.pNNPuppiTauPuppiMet_55_190 = cms.Path(process.NNPuppiTauPuppiMet)
process.pPuppiHT400 = cms.Path(process.PuppiHT400)
process.pPuppiHT450 = cms.Path(process.PuppiHT450)
process.pPuppiMET200 = cms.Path(process.PuppiMET200)
process.pPuppiMHT140 = cms.Path(process.PuppiMHT140)
process.pPuppiTauTkIsoEle45_22 = cms.Path(process.PuppiTauTkIsoEle4522)
process.pPuppiTauTkMuon42_18 = cms.Path(process.PuppiTauTkMuon4218)
process.pQuadJet70_55_40_40 = cms.Path(process.QuadJet70554040)
process.pSingleEGEle51 = cms.Path(process.SingleEGEle51)
process.pSingleIsoTkEle28 = cms.Path(process.SingleIsoTkEle28)
process.pSingleIsoTkPho36 = cms.Path(process.SingleIsoTkPho36)
process.pSinglePuppiJet230 = cms.Path(process.SinglePuppiJet230)
process.pSingleTkEle36 = cms.Path(process.SingleTkEle36)
process.pSingleTkMuon22 = cms.Path(process.SingleTkMuon22)
process.pTkEleIsoPuppiHT_26_190 = cms.Path(process.TkEleIsoPuppiHT)
process.pTkElePuppiJet_28_40_MinDR = cms.Path(process.TkElePuppiJetMinDR)
process.pTkEleTkMuon10_20 = cms.Path(process.TkEleTkMuon1020)
process.pTkMuPuppiJetPuppiMet_3_110_120 = cms.Path(process.TkMuPuppiJetPuppiMet)
process.pTkMuTriPuppiJet_12_40_dRMax_DoubleJet_dEtaMax = cms.Path(process.TkMuTriPuppiJetdRMaxDoubleJetdEtaMax)
process.pTkMuonDoubleTkEle6_17_17 = cms.Path(process.TkMuonDoubleTkEle61717)
process.pTkMuonPuppiHT6_320 = cms.Path(process.TkMuonPuppiHT6320)
process.pTkMuonTkEle7_23 = cms.Path(process.TkMuonTkEle723)
process.pTkMuonTkIsoEle7_20 = cms.Path(process.TkMuonTkIsoEle720)
process.pTripleTkMuon5_3_3 = cms.Path(process.TripleTkMuon533)
process.pTripleTkMuon_5_3_0_DoubleTkMuon_5_3_OS_MassTo9 = cms.Path(process.TripleTkMuon530OSMassMax9)
process.pTripleTkMuon_5_3p5_2p5_OS_Mass5to17 = cms.Path(process.TripleTkMuon53p52p5OSMass5to17)
process.nanoAOD_step = cms.Path(process.l1tPh2NanoSequence)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.NANOAODoutput_step = cms.EndPath(process.NANOAODoutput)

process.schedule = cms.Schedule(
    process.L1TrackTrigger_step,
    process.L1simulation_step,
    process.Phase2L1GTProducer,
    process.Phase2L1GTAlgoBlockProducer,
    process.pDoubleEGEle37_24, process.pDoubleIsoTkPho22_12,
    process.pDoublePuppiJet112_112, process.pDoublePuppiJet160_35_mass620,
    process.pDoublePuppiTau52_52, process.pDoubleTkEle25_12,
    process.pDoubleTkElePuppiHT_8_8_390, process.pDoubleTkMuPuppiHT_3_3_300,
    process.pDoubleTkMuPuppiJetPuppiMet_3_3_60_130, process.pDoubleTkMuon15_7,
    process.pDoubleTkMuonTkEle5_5_9, process.pDoubleTkMuon_4_4_OS_Dr1p2,
    process.pDoubleTkMuon_4p5_4p5_OS_Er2_Mass7to18, process.pDoubleTkMuon_OS_Er1p5_Dr1p4,
    process.pIsoTkEleEGEle22_12, process.pNNPuppiTauPuppiMet_55_190,
    process.pPuppiHT400, process.pPuppiHT450, process.pPuppiMET200, process.pPuppiMHT140,
    process.pPuppiTauTkIsoEle45_22, process.pPuppiTauTkMuon42_18,
    process.pQuadJet70_55_40_40, process.pSingleEGEle51, process.pSingleIsoTkEle28,
    process.pSingleIsoTkPho36, process.pSinglePuppiJet230, process.pSingleTkEle36,
    process.pSingleTkMuon22, process.pTkEleIsoPuppiHT_26_190,
    process.pTkElePuppiJet_28_40_MinDR, process.pTkEleTkMuon10_20,
    process.pTkMuPuppiJetPuppiMet_3_110_120,
    process.pTkMuTriPuppiJet_12_40_dRMax_DoubleJet_dEtaMax,
    process.pTkMuonDoubleTkEle6_17_17, process.pTkMuonPuppiHT6_320,
    process.pTkMuonTkEle7_23, process.pTkMuonTkIsoEle7_20, process.pTripleTkMuon5_3_3,
    process.pTripleTkMuon_5_3_0_DoubleTkMuon_5_3_OS_MassTo9,
    process.pTripleTkMuon_5_3p5_2p5_OS_Mass5to17,
    process.nanoAOD_step,
    process.endjob_step,
    process.NANOAODoutput_step,
)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# ---------------------------------------------------------------------------
# Unified @L1PFTrkNanowithGen customise chain (what the autoNANO flavor expands to)
# ---------------------------------------------------------------------------
from DPGAnalysis.Phase2L1TNanoAOD.l1tPh2Nano_cff import (
    addPh2L1Objects, addPh2GTObjects, addPh2L1Tracks, addPh2L1PFCandidates,
    addPh2L1TrackTruth, addPh2L1DisplacedVertices, addGenObjects)

process = addPh2L1Objects(process)
process = addPh2GTObjects(process)
process = addPh2L1Tracks(process)
process = addPh2L1PFCandidates(process)
process = addPh2L1TrackTruth(process)
process = addPh2L1DisplacedVertices(process)
process = addGenObjects(process)

# ---------------------------------------------------------------------------
# Layer the SmartPixels digiRefit tables (WF1 coexist) into the SAME nano output.
# truthSource='inJob': this job re-emulates L1TrackTrigger + SimL1Emulator, so the
# in-process TT truth associators run and the re-emulated tracklet tracks carry a
# real helixCovMat -> seedCovMode='trackCov' is valid on the covMatrix build.
# ---------------------------------------------------------------------------
from L1Trigger.Phase3SmartPixels.customizeSmartPixels_cff import smartPixelsCoexist

process = smartPixelsCoexist(
    process,
    variants=[("digiRefit", options.activeSP)],
    truthSource="inJob",
    addNanoTables=True,
    digiRefitConfig={
        "pixelavAngleSet": _pixelavAngleSet,
        "seedCovMode": options.seedCovMode,
    },
)

# re-emul fix (b), reasserted: coexist(truthSource='inJob') relies on the in-process
# cluster associator; repoint it at the file's simSiPixelDigis:Tracker simlinks in
# case the sequence (re)created the module after the first pass above.
for _m in ["TTClusterAssociatorFromPixelDigis"]:
  if hasattr(process, _m):
    getattr(process, _m).digiSimLinks = cms.InputTag("simSiPixelDigis", "Tracker")

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
