"""Posture-A (fromFile) digiRefit job for PU RelVals.

Nothing is remade: the digiRefit producer consumes the file's TTTracks,
simSiPixelDigis:Pixel digis+simlinks, association maps and TrackingParticles
(all HLT-process, mutually consistent, WITH pileup). Re-running DIGI here would
re-digitize the signal-only g4SimHits and silently lose the PU.

File tracks are old-layout (schema-evolved, all-zero helixCovMat), so
seedCovMode='parametrized' is REQUIRED — 'trackCov' throws
SmartPixelsSeedCovMissing, by design.

Run (in the container, cmsenv):
  cmsRun L1Trigger/Phase3SmartPixels/test/runDigiRefitFromFile_cfg.py \
      maxEvents=100 useAngles=alpha \
      inputFiles=file:/work/testfiles/RelValTTbar_14TeV_PU_150X_mcRun4_realistic_v1_STD_D121_RegeneratedGS_PU-v1.root \
      outputFile=digirefit_pu.root
"""
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
# CMSSW_17 backport (Branch B): the 140X-compatible MC uses the D110 layout /
# Phase2C17I13M9 era (see mem:smartpixels-cmssw17-backport-state), not D121/C22.
from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9

options = VarParsing('analysis')
options.register('useAngles', 'alpha', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                 "digiRefit angle usage: none | alpha | alphaBeta")
options.register('activeSP', '1100', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                 "active layer encoding for the digiRefit variant")
options.register('pixelavAngleSet', '/work/spxsmoke/pixelav_angle_example.json',
                 VarParsing.multiplicity.singleton, VarParsing.varType.string,
                 "PixelAV angle-response payload path")
options.setDefault('maxEvents', 100)
options.setDefault('inputFiles',
                   'file:/work/testfiles/TT_TuneCP5_14TeV-powheg-pythia8_PU200_Trk1GeV_140X_mcRun4_realistic_v4-v2_file1_100evts_numEvent100.root')
options.setDefault('outputFile', 'digirefit_pu.root')
options.parseArguments()

process = cms.Process("SPXREFIT", Phase2C17I13M9)
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.Geometry.GeometryExtendedRun4D110Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T35', '')
# ES-only provider for the hph::Setup the producer consumes (no event products
# remade). CMSSW_17 (Branch B) layout: hph::Setup is keyed on hph::SetupRcd, which
# depends on tt::SetupRcd (L1Trigger.TrackTrigger.Setup_cff) and
# trackerTFP::DataFormatsRcd (L1Trigger.TrackerTFP.DataFormats_cff). ProducerHPH_cff
# imports and provides all of those itself, so loading it alone suffices; the old
# L1Trigger.TrackerDTC.Setup_cff path does not exist in this release.
process.load('L1Trigger.TrackTrigger.TrackTrigger_cff')  # TTStub algorithm records
process.load('L1Trigger.TrackFindingTracklet.ProducerHPH_cff')

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))
process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(options.inputFiles))
process.MessageLogger = cms.Service("MessageLogger",
                                    cerr=cms.untracked.PSet(FwkReport=cms.untracked.PSet(reportEvery=cms.untracked.int32(10))))

from L1Trigger.Phase3SmartPixels.customizeSmartPixels_cff import addSmartPixelsTrackProducerVariants
process, spxModules = addSmartPixelsTrackProducerVariants(
    process,
    variants=[("passthrough", None), ("digiRefit", options.activeSP)],
    digiRefitConfig={"pixelavAngleSet": options.pixelavAngleSet,
                     "useAngles": options.useAngles,
                     "seedCovMode": "parametrized"})

process.p = cms.Path()
for m in spxModules:
    process.p += getattr(process, m)

process.out = cms.OutputModule("PoolOutputModule",
    fileName=cms.untracked.string(options.outputFile),
    outputCommands=cms.untracked.vstring(
        "drop *",
        "keep *_l1tSmartPixelsTrackProducer*_*_*",
        "keep *_l1tTTTracksFromTrackletEmulation_Level1TTTracks_*",
        "keep *_TTTrackAssociatorFromPixelDigis_Level1TTTracks_*",
        "keep *_mix_MergedTrackTruth_*",
    ))
process.e = cms.EndPath(process.out)
