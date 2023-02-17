import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.ParameterSet.VarParsing as VarParsing
from L1Trigger.VertexFinder.VertexProducer_cff import VertexProducer
from L1Trigger.L1TTrackMatch.L1TrackVertexAssociationProducer_cfi import L1TrackVertexAssociationProducer
from L1Trigger.L1TTrackMatch.L1TrackSelectionProducer_cfi import L1TrackSelectionProducer

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

process.load('Configuration.Geometry.GeometryExtended2026D77Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D77_cff')
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

process.load('L1Trigger.L1TTrackMatch.L1GTTInputProducer_cfi')
process.load("L1Trigger.L1TTrackMatch.L1TrackJetEmulationProducer_cfi")
process.load("L1Trigger.L1TTrackMatch.L1TkHTMissEmulatorProducer_cfi")
process.load("L1Trigger.L1TTrackMatch.L1TrackerEtMissEmulatorProducer_cfi")
process.load('L1Trigger.DemonstratorTools.GTTFileWriter_cff')

process.L1GTTInputProducer.debug = cms.int32(options.debug)

process.L1TrackSelectionProducer = L1TrackSelectionProducer.clone()
process.L1TrackSelectionProducer.processSimulatedTracks = cms.bool(True)
process.L1TrackSelectionProducer.processEmulatedTracks = cms.bool(True)

process.VertexProducerSim = VertexProducer.clone()
process.VertexProducerSim.l1TracksInputTag = cms.InputTag("L1TrackSelectionProducer", "Level1TTTracksSelected")
process.VertexProducerSim.l1VertexCollectionName = cms.string("l1vertices")
process.VertexProducerSim.VertexReconstruction.Algorithm = cms.string("fastHisto")
process.VertexProducerSim.VertexReconstruction.VxMinTrackPt = cms.double(0.0)
process.VertexProducerSim.debug = options.debug

process.VertexProducerEmu = VertexProducer.clone()
process.VertexProducerEmu.l1TracksInputTag = cms.InputTag("L1TrackSelectionProducer", "Level1TTTracksSelectedEmulation")
process.VertexProducerEmu.l1VertexCollectionName = cms.string("l1vertices")
process.VertexProducerEmu.VertexReconstruction.Algorithm = cms.string("fastHistoEmulation")
process.VertexProducerEmu.VertexReconstruction.VxMinTrackPt = cms.double(0.0)
process.VertexProducerEmu.debug = options.debug

process.L1TrackVertexAssociationProducer = L1TrackVertexAssociationProducer.clone()
process.L1TrackVertexAssociationProducer.l1SelectedTracksInputTag = cms.InputTag("L1TrackSelectionProducer", "Level1TTTracksSelected")
process.L1TrackVertexAssociationProducer.l1SelectedTracksEmulationInputTag = cms.InputTag("L1TrackSelectionProducer", "Level1TTTracksSelectedEmulation")
process.L1TrackVertexAssociationProducer.l1VerticesInputTag = cms.InputTag("VertexProducerSim", "l1vertices")
process.L1TrackVertexAssociationProducer.l1VerticesEmulationInputTag = cms.InputTag("VertexProducerEmu", "l1verticesEmulation")
process.L1TrackVertexAssociationProducer.outputCollectionName = cms.string("Level1TTTracksSelectedAssociated")
process.L1TrackVertexAssociationProducer.processSimulatedTracks = cms.bool(True)
process.L1TrackVertexAssociationProducer.processEmulatedTracks = cms.bool(True)

process.L1TrackerEmuEtMiss.L1VertexInputTag = cms.InputTag("VertexProducerEmu", "l1verticesEmulation")
process.L1TrackerEmuEtMiss.L1TrackInputTag = cms.InputTag("L1TrackSelectionProducer", "Level1TTTracksSelectedEmulation")
process.L1TrackerEmuEtMiss.L1TrackAssociatedInputTag = cms.InputTag("L1TrackVertexAssociationProducer", "Level1TTTracksSelectedAssociatedEmulation")
process.L1TrackerEmuEtMiss.debug = options.debug

process.L1TrackJetsEmulation.L1TrackInputTag= cms.InputTag("L1GTTInputProducer", "Level1TTTracksConverted")
process.L1TrackJetsEmulation.VertexInputTag = cms.InputTag("VertexProducerEmu", "l1verticesEmulation")
process.L1TrackerEmuHTMiss.L1TkJetEmulationInputTag = cms.InputTag("L1TrackJetsEmulation", "L1TrackJets")
process.L1TrackerEmuHTMiss.L1MHTCollectionName = cms.string("L1TrackerEmuHTMiss")
process.L1TrackerEmuHTMiss.debug = (options.debug > 0)

if options.debug:
    process.MessageLogger.cerr.INFO.limit = cms.untracked.int32(1000000000)
    process.MessageLogger.suppressInfo = cms.untracked.vstring('CondDBESSource', 'PoolDBESSource')
    process.MessageLogger.cerr.CondDBESSource = cms.untracked.PSet(
        limit = cms.untracked.int32(0)
    )

process.GTTFileWriter.format = cms.untracked.string(options.format)
process.GTTFileWriter.selectedTracks = cms.untracked.InputTag("L1TrackSelectionProducer", "Level1TTTracksSelectedEmulation")
process.GTTFileWriter.vertices = cms.untracked.InputTag("VertexProducerEmu", "l1verticesEmulation")
process.GTTFileWriter.vertexAssociatedTracks = cms.untracked.InputTag("L1TrackVertexAssociationProducer", "Level1TTTracksSelectedAssociatedEmulation")
process.GTTFileWriter.outputCorrelatorFilename = cms.untracked.string("L1GTTOutputToCorrelatorFile")
process.GTTFileWriter.outputGlobalTriggerFilename = cms.untracked.string("L1GTTOutputToGlobalTriggerFile")
process.GTTFileWriter.selectedTracksFilename = cms.untracked.string("L1GTTSelectedTracksFile")
process.GTTFileWriter.vertexAssociatedTracksFilename = cms.untracked.string("L1GTTVertexAssociatedTracksFile")

process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.Timing = cms.Service("Timing", summaryOnly = cms.untracked.bool(True))


process.p = cms.Path(process.GTTFileWriter)
process.p.associate(cms.Task(process.L1GTTInputProducer, 
                             process.L1TrackSelectionProducer,
                             process.VertexProducerEmu, 
                             process.VertexProducerSim, 
                             process.L1TrackVertexAssociationProducer,
                             process.L1TrackJetsEmulation, 
                             process.L1TrackerEmuHTMiss, 
                             process.L1TrackerEmuEtMiss
                         )
                )

# Tools for debugging dependencies
#option 1: graph file. process output with ```dot dependency.gv -Tpdf -o dependency.pdf```
# process.load("FWCore.Services.DependencyGraph_cfi")
# process.DependencyGraph.fileName = 'dependency.gv'
# process.DependencyGraph.showPathDependencies = True
# from FWCore.ParameterSet.Utilities import moduleLabelsInSequences
# process.DependencyGraph.highlightModules = moduleLabelsInSequences(process.p)

#option 2: print consumptions and paths
# process.add_(cms.Service("Tracer", dumpPathsAndConsumes=cms.untracked.bool(True)))

#option 3: print registry
# process.add_(cms.Service("ProductRegistryDumper"))
