"""Cluster ADC-map ntuple over the PU200 ttbar D121 RelVal (GEN-SIM-DIGI-RAW).

Rebuilds the SmartPixels-owned cluster chain exactly as the nano jobs do
(customizeSmartPixels_cff.ensureSmartPixelsRecHits: spixPixelClusters ->
spixPixelRecHits -> spixSmartPixelsRecHits, from simSiPixelDigis:Pixel) with the
v4fixed PixelAV payload + noise v6 (as spix_sweep15_500), then writes one TTree row
per TBPX cluster.

  cmsRun runClusterAdcNtuple_cfg.py maxEvents=5 outputFile=/work/clusterntuple/test5.root
"""
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from Configuration.Eras.Era_Phase2C22I13M9_cff import Phase2C22I13M9

BASE = "file:/host_volumes/WDMac/smartpixels-cmssw-testfiles/"
DS = "RelValTTbar_14TeV_PU_150X_mcRun4_realistic_v1_STD_D121_RegeneratedGS_PU-v1"

o = VarParsing("analysis")
o.register("fileIdx", "1", VarParsing.multiplicity.singleton, VarParsing.varType.string,
           "comma list of RelVal file indices")
o.register("nThreads", 4, VarParsing.multiplicity.singleton, VarParsing.varType.int, "threads")
o.register("skip", 0, VarParsing.multiplicity.singleton, VarParsing.varType.int, "skipEvents")
o.setDefault("maxEvents", 5)
o.setDefault("outputFile", "/work/clusterntuple/test5.root")
o.parseArguments()

process = cms.Process("SPIXCLNTUP", Phase2C22I13M9)
process.load("Configuration.StandardSequences.Services_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.Geometry.GeometryExtendedRun4D121Reco_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, "auto:phase2_realistic_T35", "")

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(o.maxEvents))
files = [BASE + "%s_file%s.root" % (DS, i.strip()) for i in o.fileIdx.split(",")]
process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(*files),
                            skipEvents=cms.untracked.uint32(o.skip))
process.options.numberOfThreads = o.nThreads
process.options.numberOfStreams = 0
process.MessageLogger.cerr.FwkReport.reportEvery = 10

from L1Trigger.Phase3SmartPixels.customizeSmartPixels_cff import ensureSmartPixelsRecHits
ensureSmartPixelsRecHits(process,
                         "/work/spxsmoke/spix_angle_response_Conv1D_Full-2bit_v4fixed.json",
                         "/work/spxsmoke/smarthit_noise_v6.json")

process.clusterAdcNtuple = cms.EDAnalyzer("SmartPixelsClusterAdcNtuplizer")
process.TFileService = cms.Service("TFileService", fileName=cms.string(o.outputFile))
process.p = cms.Path(process.clusterAdcNtuple, process.spixSmartPixelsRecHitTask)
