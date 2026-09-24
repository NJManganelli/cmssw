# Dump the inner-tracker pixel module geometry (origin + local axes as global
# vectors, bounds, pitch) to JSON for offline projection studies.
#
# Era, geometry and GlobalTag MUST match the job that produced the nano being
# analysed; the defaults are the D121 PU200 SmartPixels production
# (Phase2C22I13M9, ExtendedRun4D121, auto:phase2_realistic_T35).
#
#   cmsRun dumpModuleGeometry_cfg.py jsonFile=/work/spix_module_geometry_D121.json
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from Configuration.Eras.Era_Phase2C22I13M9_cff import Phase2C22I13M9

opts = VarParsing()
opts.register("jsonFile", "spix_module_geometry_D121.json", VarParsing.multiplicity.singleton,
              VarParsing.varType.string, "JSON output path")
opts.register("globalTag", "auto:phase2_realistic_T35", VarParsing.multiplicity.singleton,
              VarParsing.varType.string, "GlobalTag (match the nano production)")
opts.parseArguments()

process = cms.Process("SPIXGEOM", Phase2C22I13M9)
process.load("Configuration.StandardSequences.Services_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.Geometry.GeometryExtendedRun4D121Reco_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, opts.globalTag, "")

process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(1))

process.dump = cms.EDAnalyzer("SmartPixelsModuleGeometryDumper",
                              outputFile=cms.string(opts.jsonFile))
process.p = cms.Path(process.dump)
