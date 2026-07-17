"""Pure file-reading job that runs SmartPixelsPayloadAnalyzer over a RelVal.

Nothing is remade: this reads the already-persisted TTTracks, association maps,
TrackingParticles, simSiPixelDigis:Pixel digis+simlinks and g4SimHits, and writes
a TFileService ntuple of L1Track x TBPX-layer crossings for payload derivation.

Run (in the container, cmsenv):
  cmsRun L1Trigger/Phase3SmartPixels/test/runSmartPixelsPayloadNtuple_cfg.py \
      maxEvents=2 \
      inputFiles=file:/work/testfiles/RelValTTbar_D121_noPU.root \
      outputFile=payload_ntuple.root
"""
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from Configuration.Eras.Era_Phase2C22I13M9_cff import Phase2C22I13M9

options = VarParsing('analysis')
options.register('nLayers', 4, VarParsing.multiplicity.singleton, VarParsing.varType.int,
                 "Number of TBPX layers to consider (1..4)")
options.register('windowRPhi', 0.30, VarParsing.multiplicity.singleton, VarParsing.varType.float,
                 "r-phi MEASUREMENT window half-width [cm] (wide on purpose; see cfi)")
options.register('windowZ', 0.50, VarParsing.multiplicity.singleton, VarParsing.varType.float,
                 "z MEASUREMENT window half-width [cm] (wide on purpose; see cfi)")
options.register('debug', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
                 "Verbose per-crossing logging")
options.setDefault('maxEvents', 2)
options.setDefault('inputFiles', 'file:/work/testfiles/RelValTTbar_D121_noPU.root')
options.setDefault('outputFile', 'payload_ntuple.root')
options.parseArguments()

process = cms.Process("SPXPAYLOAD", Phase2C22I13M9)

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.Geometry.GeometryExtendedRun4D121Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T35', '')

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))
process.source = cms.Source("PoolSource",
                            fileNames=cms.untracked.vstring(options.inputFiles))

process.load('L1Trigger.Phase3SmartPixels.smartPixelsPayloadAnalyzer_cfi')
process.smartPixelsPayloadAnalyzer.nLayers = cms.int32(options.nLayers)
process.smartPixelsPayloadAnalyzer.windowRPhi = cms.double(options.windowRPhi)
process.smartPixelsPayloadAnalyzer.windowZ = cms.double(options.windowZ)
process.smartPixelsPayloadAnalyzer.debug = cms.bool(options.debug)

process.TFileService = cms.Service("TFileService",
                                   fileName=cms.string(options.outputFile))

process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.p = cms.Path(process.smartPixelsPayloadAnalyzer)
