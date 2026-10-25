# Dump the inner-tracker pixel module geometry (origin + local axes as global
# vectors, bounds, pitch) to JSON for offline projection studies.
#
# Era, geometry and GlobalTag MUST match the job that produced the nano being
# analysed; the defaults are the D121 PU200 SmartPixels production
# (Phase2C22I13M9, ExtendedRun4D121, auto:phase2_realistic_T35). campaign=spring24
# is the Phase2Spring24 production (Phase2C17I13M9, ExtendedRun4D110; a CMSSW_17 area).
#
#   cmsRun dumpModuleGeometry_cfg.py jsonFile=spix_module_geometry_D121.json
#   cmsRun dumpModuleGeometry_cfg.py campaign=spring24 jsonFile=spix_module_geometry_D110.json
import importlib

import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

CAMPAIGNS = {"d121": ("Phase2C22I13M9", "ExtendedRun4D121"), "spring24": ("Phase2C17I13M9", "ExtendedRun4D110")}

opts = VarParsing()
opts.register("jsonFile", "spix_module_geometry_D121.json", VarParsing.multiplicity.singleton,
              VarParsing.varType.string, "JSON output path")
opts.register("globalTag", "auto:phase2_realistic_T35", VarParsing.multiplicity.singleton,
              VarParsing.varType.string, "GlobalTag (match the nano production)")
opts.register("campaign", "d121", VarParsing.multiplicity.singleton, VarParsing.varType.string,
              "d121 or spring24: the era and geometry of the nano production")
opts.parseArguments()
if opts.campaign not in CAMPAIGNS:
    raise SystemExit(f"campaign {opts.campaign!r} not in {sorted(CAMPAIGNS)}")
eraName, geometry = CAMPAIGNS[opts.campaign]
era = getattr(importlib.import_module(f"Configuration.Eras.Era_{eraName}_cff"), eraName)

process = cms.Process("SPIXGEOM", era)
process.load("Configuration.StandardSequences.Services_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load(f"Configuration.Geometry.Geometry{geometry}Reco_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, opts.globalTag, "")

process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(1))

process.dump = cms.EDAnalyzer("SmartPixelsModuleGeometryDumper",
                              outputFile=cms.string(opts.jsonFile))
process.p = cms.Path(process.dump)
