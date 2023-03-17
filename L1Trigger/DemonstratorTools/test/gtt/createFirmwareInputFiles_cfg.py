import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.ParameterSet.VarParsing as VarParsing

# PART 1 : PARSE ARGUMENTS

options = VarParsing.VarParsing ('analysis')
options.register('debug',
                 0,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "Print out additional debugging information")
options.register ('format',
                  'EMP', # default value
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.string,
                  "File format (APx, EMP or X20)")
options.register('threads',
                 1, # default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "Number of threads to run")
options.register('streams',
                 0, # default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "Number of streams to run")
options.parseArguments()

inputFiles = []
for filePath in options.inputFiles:
    if filePath.endswith(".root"):
        inputFiles.append(filePath)
    elif filePath.endswith("_cff.py"):
        filePath = filePath.replace("/python/","/")
        filePath = filePath.replace("/", ".")
        inputFilesImport = getattr(__import__(filePath.strip(".py"),fromlist=["readFiles"]),"readFiles")
        inputFiles.extend( inputFilesImport )
    else:
        inputFiles += FileUtils.loadListFromFile(filePath)

# PART 2: SETUP MAIN CMSSW PROCESS 

process = cms.Process("GTTFileWriter")

process.load('Configuration.Geometry.GeometryExtended2026D95Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D95_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')
process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(inputFiles),
    inputCommands = cms.untracked.vstring("keep *", "drop l1tTkPrimaryVertexs_L1TkPrimaryVertex__*")
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )
process.options = cms.untracked.PSet(
    numberOfThreads = cms.untracked.uint32(options.threads),
    numberOfStreams = cms.untracked.uint32(options.streams if options.streams>0 else 0)
)

process.load('L1Trigger.L1TTrackMatch.l1tGTTInputProducer_cfi')
process.load('L1Trigger.L1TTrackMatch.l1tTrackSelectionProducer_cfi')
process.load('L1Trigger.VertexFinder.l1tVertexProducer_cfi')
process.load('L1Trigger.L1TTrackMatch.l1tTrackVertexAssociationProducer_cfi')
process.load('L1Trigger.L1TTrackMatch.l1tTrackJetsEmulation_cfi')
process.load('L1Trigger.L1TTrackMatch.l1tTrackerEmuHTMiss_cfi')
process.load('L1Trigger.L1TTrackMatch.l1tTrackerEmuEtMiss_cfi')
process.load('L1Trigger.DemonstratorTools.l1tGTTFileWriter_cfi')
                                                            
process.l1tGTTInputProducer.debug = cms.int32(options.debug)

process.l1tTrackSelectionProducer.processSimulatedTracks = cms.bool(False)
process.l1tVertexFinderEmu.VertexReconstruction.VxMinTrackPt = cms.double(0.0)
process.l1tVertexFinderEmu.debug = options.debug
process.l1tTrackVertexAssociationProducer.processSimulatedTracks = cms.bool(False)

process.l1tTrackSelectionProducerForEtMiss.processSimulatedTracks = cms.bool(False)
process.l1tTrackVertexAssociationProducerForEtMiss.processSimulatedTracks = cms.bool(False)
process.l1tTrackerEmuEtMiss.debug = options.debug

process.l1tTrackSelectionProducerForJets.processSimulatedTracks = cms.bool(False)
process.l1tTrackVertexAssociationProducerForJets.processSimulatedTracks = cms.bool(False)
process.l1tTrackerEmuHTMiss.debug = (options.debug > 0)

if options.debug:
    process.MessageLogger.cerr.INFO.limit = cms.untracked.int32(1000000000)
    process.MessageLogger.suppressInfo = cms.untracked.vstring('CondDBESSource', 'PoolDBESSource')
    process.MessageLogger.cerr.CondDBESSource = cms.untracked.PSet(
        limit = cms.untracked.int32(0)
    )

process.l1tGTTFileWriter.format = cms.untracked.string(options.format) #FIXME Put all this into the default GTTFileWriter
process.l1tGTTFileWriter.tracks = cms.untracked.InputTag("l1tTTTracksFromTrackletEmulation", "Level1TTTracks")
process.l1tGTTFileWriter.convertedTracks = cms.untracked.InputTag("l1tGTTInputProducer", "Level1TTTracksConverted")
process.l1tGTTFileWriter.selectedTracks = cms.untracked.InputTag("l1tTrackSelectionProducer", "Level1TTTracksSelectedEmulation")
process.l1tGTTFileWriter.vertices = cms.untracked.InputTag("l1tVertexFinderEmu", "l1tVerticesEmulation")
process.l1tGTTFileWriter.vertexAssociatedTracks = cms.untracked.InputTag("l1tTrackVertexAssociationProducer", "Level1TTTracksSelectedAssociatedEmulation")
process.l1tGTTFileWriter.jets = cms.untracked.InputTag("l1tTrackJetsEmulation","L1TrackJets")
process.l1tGTTFileWriter.htmiss = cms.untracked.InputTag("l1tTrackerEmuHTMiss", "l1tTrackerEmuHTMiss")
process.l1tGTTFileWriter.etmiss = cms.untracked.InputTag("l1tTrackerEmuEtMiss", "l1tTrackerEmuEtMiss")
process.l1tGTTFileWriter.outputCorrelatorFilename = cms.untracked.string("L1GTTOutputToCorrelatorFile")
process.l1tGTTFileWriter.outputGlobalTriggerFilename = cms.untracked.string("L1GTTOutputToGlobalTriggerFile")
process.l1tGTTFileWriter.selectedTracksFilename = cms.untracked.string("L1GTTSelectedTracksFile")
process.l1tGTTFileWriter.vertexAssociatedTracksFilename = cms.untracked.string("L1GTTVertexAssociatedTracksFile")

process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.Timing = cms.Service("Timing", summaryOnly = cms.untracked.bool(True))


process.p = cms.Path(process.l1tGTTFileWriter)
process.p.associate(cms.Task(process.l1tGTTInputProducer, 
                             process.l1tTrackSelectionProducer,
                             process.l1tVertexFinderEmu, 
                             process.l1tTrackVertexAssociationProducer,
                             process.l1tTrackSelectionProducerForJets,
                             process.l1tTrackVertexAssociationProducerForJets,
                             process.l1tTrackJetsEmulation, 
                             process.l1tTrackerEmuHTMiss, 
                             process.l1tTrackSelectionProducerForEtMiss,
                             process.l1tTrackVertexAssociationProducerForEtMiss,
                             process.l1tTrackerEmuEtMiss,
                         )
                )

# Tools for debugging dependencies
#option 1: graph file. process output with ```dot dependency.gv -Tpdf -o dependency.pdf```
process.load("FWCore.Services.DependencyGraph_cfi")
process.DependencyGraph.fileName = 'dependency_mar16.gv'
process.DependencyGraph.showPathDependencies = True
from FWCore.ParameterSet.Utilities import moduleLabelsInSequences
process.DependencyGraph.highlightModules = moduleLabelsInSequences(process.p)
