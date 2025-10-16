import FWCore.ParameterSet.Config as cms

# configures track finding s/w to use KF emulator instead of KF simulator
def newKFConfig(process):
  process.l1tTTTracksFromTrackletEmulation.Fakefit = True

# configures track finding s/w to behave as track finding f/w
def fwConfig(process):
  newKFConfig(process)
  process.TrackTriggerSetup.Firmware.FreqBE = 240 # Frequency of DTC & KF (determines truncation)
  process.l1tTTTracksFromTrackletEmulation.RemovalType = ""
  process.l1tTTTracksFromTrackletEmulation.DoMultipleMatches = False
  process.l1tTTTracksFromTrackletEmulation.StoreTrackBuilderOutput = True

# configures track finding s/w to behave as a subchain of processing steps
def reducedConfig(process):
  fwConfig(process)
  process.TrackTriggerSetup.KalmanFilter.NumWorker = 1
  process.ChannelAssignment.SeedTypes = cms.vstring( "L1L2" )
  process.ChannelAssignment.SeedTypesSeedLayers = cms.PSet( L1L2 = cms.vint32( 1,  2 ) )
  process.ChannelAssignment.SeedTypesProjectionLayers = cms.PSet( L1L2 = cms.vint32(  3,  4,  5,  6 ) )
  # this are tt::Setup::dtcId in order as in process.l1tTTTracksFromTrackletEmulation.processingModulesFile translated by 
  # reverssing naming logic described in L1FPGATrackProducer
  # TO DO: Eliminate cfg param IRChannelsIn by taking this info from Tracklet wiring map.
  process.ChannelAssignment.IRChannelsIn = cms.vint32( 0, 1, 25, 2, 26, 4, 5, 29, 6, 30, 7, 31, 8, 9, 33 )
  process.l1tTTTracksFromTrackletEmulation.Reduced = True
  process.l1tTTTracksFromTrackletEmulation.memoryModulesFile = 'L1Trigger/TrackFindingTracklet/data/reduced_memorymodules.dat'
  process.l1tTTTracksFromTrackletEmulation.processingModulesFile = 'L1Trigger/TrackFindingTracklet/data/reduced_processingmodules.dat'
  process.l1tTTTracksFromTrackletEmulation.wiresFile = 'L1Trigger/TrackFindingTracklet/data/reduced_wires.dat'

# configures pure tracklet algorithm (as opposed to Hybrid algorithm)
def trackletConfig(process):
    process.l1tTTTracksFromTrackletEmulation.fitPatternFile = cms.FileInPath('L1Trigger/TrackFindingTracklet/data/fitpattern.txt') 

def injectSmartPixelsTrackProducer(process,
                                   smartPixelsEmulatorMode="passthrough",
                                   addAssociation=True,
                                   l1tSmartPixelsTrackProducerLabel="l1tSmartPixelsTrackProducer",
                                   l1tSmartPixelsTrackProducerExtendedLabel="l1tSmartPixelsTrackProducerExtended",
                                   skipModuleTypes=None,
                                   printProcessInfo=False):
  if skipModuleTypes is None:
    skipModuleTypes = ["L1SmartPixelsTrackProducer", "TTTrackAssociator_Phase2TrackerDigi_", "L1FPGATrackProducer", "SimpleL1TTTrackCandidateFlatTableProducer"]
  #print(help(cms.Process))
  def print_helper(item):
    if hasattr(item, "items"):
      for key, value in item.items():
        print("\t>>", key, "== ", value)
    else:
      print("\t==", item)

  if printProcessInfo:
    print("====Process Info====")
    print("source_\n")
    print_helper(process.source_())
    print("schedule_\n")
    print_helper(process.schedule_())
    print("sequences_\n")
    print_helper(process.sequences_())
    print("paths_\n")
    print_helper(process.paths_())
    print("endpaths_\n")
    print_helper(process.endpaths_())
    print("producers_\n")
    print_helper(process.producers_())
    print("filters_\n")
    print_helper(process.filters_())
    print("analyzers_\n")
    print_helper(process.analyzers_())

  for modname, module in process.producers_().items():
    if hasattr(module, '_TypedParameterizable__type') and module._TypedParameterizable__type in skipModuleTypes:
      print("Skipping module type ", module._TypedParameterizable__type, " with module name ", modname)
      continue
    for param_name in module.parameterNames_():
      param = getattr(module, param_name)
      if isinstance(param, cms.InputTag):
        if param.getModuleLabel() == "l1tTTTracksFromTracklet":
          print("l1tTTTracksFromTracklet --> ", modname, param_name, param)
        if param.getModuleLabel() == "l1tTTTracksFromExtendedTrack":
          print("l1tTTTracksFromExtendedTracklet --> ", modname, param_name, param)
        if param.getModuleLabel() == "l1tTTTracksFromTrackletEmulation":
          print("module._TypedParameterizable__type --> ", module._TypedParameterizable__type)
          print("l1tTTTracksFromTrackletEmulation --> ", modname, param_name, param)
          setattr(module, param_name, cms.InputTag(l1tSmartPixelsTrackProducerLabel, "Level1TTTracks"))
          new_param = getattr(module, param_name)
          print("                                 --> ", modname, param_name, new_param)
        if param.getModuleLabel() == "l1tTTTracksFromExtendedTrackletEmulation":
          print("l1tTTTracksFromExtendedTrackletEmulation --> ", modname, param_name, param)
          setattr(module, param_name, cms.InputTag(l1tSmartPixelsTrackProducerExtendedLabel, "Level1TTTracks"))
          new_param = getattr(module, param_name)
          print("                                         --> ", modname, param_name, new_param)


  # globalReplace(self, label, new) -> "Replace the item with label 'label' by object 'new' in the process and all sequences/paths/tasks"
  #for label, module in process.producers_().items():
  #  print(label, module)

  # Ensure L1TrackTrigger is run
  #process.load('Configuration.StandardSequences.L1TrackTrigger_cff')
  #from L1Trigger.TrackTrigger.TrackTrigger_cff import *
  #from SimTracker.TrackTriggerAssociation.TrackTriggerAssociator_cff import *
  #from L1Trigger.TrackerDTC.ProducerED_cff import *
  #from L1Trigger.TrackFindingTracklet.L1HybridEmulationTracks_cff import *

  #L1TrackTrigger=cms.Sequence(TrackTriggerClustersStubs*TrackTriggerAssociatorClustersStubs*TrackerDTCProducer)


  # We Indiana-Jones idol-swap the l1tTTTracksFromTrackletEmulation and our L1SmartPixelsTrackProducer
  # Hahaha no we don't, we get the boulder instead, lets make an auto-scalpel instead...
  #process.load('L1Trigger.TrackFindingTracklet.l1tTTTracksFromTrackletEmulation_cfi')
  #process.l1tTTTracksFromTrackletEmulationOriginal = process.l1tTTTracksFromTrackletEmulation.clone()
  #process.l1tTTTracksFromExtendedTrackletEmulationOriginal = process.l1tTTTracksFromExtendedTrackletEmulation.clone()

  # this should stick to the original tracks, as we use these as input to l1tSmartPixelsTrackProducer/Extended
  #process.load('SimTracker.TrackTriggerAssociation.TTTrackAssociation_cfi')
  #process.TTTrackAssociatorFromPixelDigis.TTTracks = cms.VInputTag(cms.InputTag("l1tTTTracksFromTrackletEmulationOriginal", "Level1TTTracks"))
  #process.TTTrackAssociatorFromPixelDigisExtended = process.TTTrackAssociatorFromPixelDigis.clone(
  #    TTTracks = cms.VInputTag(cms.InputTag("l1tTTTracksFromExtendedTrackletEmulationOriginal", "Level1TTTracks") )
  #)


  if addAssociation:
    process.load('L1Trigger.TrackFindingTracklet.l1tSmartPixelsTrackProducer_cfi')
    #process.l1tTTTracksFromTrackletEmulation = process.l1tSmartPixelsTrackProducer.clone()
    #process.l1tTTTracksFromTrackletEmulation.L1TrackInputTag = cms.InputTag("l1tTTTracksFromTrackletEmulationOriginal", "Level1TTTracks")
    #process.l1tTTTracksFromExtendedTrackletEmulation = process.l1tSmartPixelsTrackProducerExtended.clone()
    #process.l1tTTTracksFromExtendedTrackletEmulation.L1TrackInputTag = cms.InputTag("l1tTTTracksFromExtendedTrackletEmulationOriginal", "Level1TTTracks")
    process.l1tSmartPixelsTrackProducer.smartPixelsEmulatorMode = cms.string(smartPixelsEmulatorMode)
    process.l1tSmartPixelsTrackProducerExtended.smartPixelsEmulatorMode = cms.string(smartPixelsEmulatorMode)
    process.l1tSmartPixelsTrackProducerTask = cms.Task(process.l1tSmartPixelsTrackProducer, process.l1tSmartPixelsTrackProducerExtended)
    #process.l1tSmartPixelsTrackProducerPath = cms.Path(process.l1tSmartPixelsTrackProducerTask)

    print("Associating L1SmartPixelsTrackProducer into the process")
    #process.p.associate(process.l1tSmartPixelsTrackProducerTask)
    for pathname, cmspath in process.paths_().items():
      cmspath.associate(process.l1tSmartPixelsTrackProducerTask)

  return process

def injectSmartPixelsTrackProducer_passthrough(process):
  return injectSmartPixelsTrackProducer(process, smartPixelsEmulatorMode="passthrough", skipModuleTypes=None, printProcessInfo=False)

def injectSmartPixelsTrackProducer_passthroughFloat(process):
  return injectSmartPixelsTrackProducer(process, smartPixelsEmulatorMode="passthroughFloat", skipModuleTypes=None, printProcessInfo=False)

def injectSmartPixelsTrackProducer_passthroughHW(process):
  return injectSmartPixelsTrackProducer(process, smartPixelsEmulatorMode="passthroughHW", skipModuleTypes=None, printProcessInfo=False)

def injectSmartPixelsTrackProducer_trackingParticleTruth(process):
  return injectSmartPixelsTrackProducer(process, smartPixelsEmulatorMode="trackingParticleTruth", skipModuleTypes=None, printProcessInfo=False)

def injectSmartPixelsTrackProducer_correctionlibRegression(process):
  return injectSmartPixelsTrackProducer(process, smartPixelsEmulatorMode="correctionlibRegression", skipModuleTypes=None, printProcessInfo=False)

def referenceSmartPixelsTrackProducer(process):
  return injectSmartPixelsTrackProducer(process, addAssociation=False, skipModuleTypes=None, printProcessInfo=False)

def addSmartPixelsTrackProducerVariants(process):
  regression_configs = ["0000",
                        "0001", "0010", "0100", "1000",
                        "0011", "0101", "0110", "1001", "1010", "1100",
                        "0111", "1011", "1101", "1110",
                        "1111",
                        ]
  
  process.load('L1Trigger.TrackFindingTracklet.l1tSmartPixelsTrackProducer_cfi')
  modules = []
  for sconf in ["passthrough", "passthroughFloat", "passthroughHW", "trackingParticleTruth"]:
    fourparam = f"l1tSmartPixelsTrackProducerW{sconf}"
    fiveparam = f"l1tSmartPixelsTrackProducerExtendedW{sconf}"
    modules.append(fourparam)
    modules.append(fiveparam)

    setattr(process, fourparam, process.l1tSmartPixelsTrackProducer.clone())
    getattr(process, fourparam).smartPixelsEmulatorMode = cms.string(sconf)

    setattr(process, fiveparam, process.l1tSmartPixelsTrackProducerExtended.clone())
    getattr(process, fiveparam).smartPixelsEmulatorMode = cms.string(sconf)

  for rconf in regression_configs:
    fourparam = f"l1tSmartPixelsTrackProducerWcorrectionlibRegression{rconf.replace('0', 'I').replace('1', 'A')}"
    fiveparam = f"l1tSmartPixelsTrackProducerExtendedWcorrectionlibRegression{rconf.replace('0', 'I').replace('1', 'A')}"
    modules.append(fourparam)
    modules.append(fiveparam)

    setattr(process, fourparam, process.l1tSmartPixelsTrackProducer.clone())
    getattr(process, fourparam).smartPixelsEmulatorMode = cms.string("correctionlibRegression")
    getattr(process, fourparam).smartPixelsCorrectionSet = "spixel_smear_all_configs_barrel_CalV1_compound.json"
    getattr(process, fourparam).smartPixelsActiveLayers = rconf

    setattr(process, fiveparam, process.l1tSmartPixelsTrackProducerExtended.clone())
    getattr(process, fiveparam).smartPixelsEmulatorMode = cms.string("correctionlibRegression")
    getattr(process, fiveparam).smartPixelsCorrectionSet = "spixel_smear_all_configs_barrel_CalV1_compound.json"
    getattr(process, fiveparam).smartPixelsActiveLayers = rconf

  return process, modules
