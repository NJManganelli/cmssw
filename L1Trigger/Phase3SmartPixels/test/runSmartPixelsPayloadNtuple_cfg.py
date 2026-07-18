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
options.register('windowRPhi', [], VarParsing.multiplicity.list, VarParsing.varType.float,
                 "Per-layer r-phi MEASUREMENT window half-widths [cm] (1 value = broadcast; "
                 "default captures the full L1-L4 extrapolation spread)")
options.register('windowZ', [], VarParsing.multiplicity.list, VarParsing.varType.float,
                 "Per-layer z MEASUREMENT window half-widths [cm] (1 value = broadcast)")
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

def _perLayerWindow(vals, default):
    # empty -> full-capture defaults; single value -> broadcast; else exactly 4
    if not vals:
        return list(default)
    if len(vals) == 1:
        return [float(vals[0])] * 4
    if len(vals) != 4:
        raise ValueError(f"window option needs 1 or 4 values, got {vals}")
    return [float(v) for v in vals]

process.smartPixelsPayloadAnalyzer.windowRPhi = cms.vdouble(*_perLayerWindow(options.windowRPhi, (0.15, 0.5, 1.5, 2.7)))
process.smartPixelsPayloadAnalyzer.windowZ = cms.vdouble(*_perLayerWindow(options.windowZ, (0.7, 0.6, 0.5, 0.4)))
process.smartPixelsPayloadAnalyzer.debug = cms.bool(options.debug)

process.TFileService = cms.Service("TFileService",
                                   fileName=cms.string(options.outputFile))

process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.p = cms.Path(process.smartPixelsPayloadAnalyzer)
