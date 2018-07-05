import FWCore.ParameterSet.Config as cms

process = cms.Process("Test")

process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(
        #'root://cmseos.fnal.gov//store/user/cmsdas/2018/pre_exercises/fourth_set/slimMiniAOD_data_MuEle_1.root',
        #      'file:slimMiniAOD_data_MuEle_1.root',
        'file:/eos/cms/store/user/nmangane/DoubleMuon/crab_CMSDAS_Data_analysis_test0/180704_111446/0000/slimMiniAOD_data_MuEle_1.root',
        'file:/eos/cms/store/user/nmangane/DoubleMuon/crab_CMSDAS_Data_analysis_test0/180704_111446/0000/slimMiniAOD_data_MuEle_2.root',
        'file:/eos/cms/store/user/nmangane/DoubleMuon/crab_CMSDAS_Data_analysis_test0/180704_111446/0000/slimMiniAOD_data_MuEle_3.root',
        'file:/eos/cms/store/user/nmangane/DoubleMuon/crab_CMSDAS_Data_analysis_test0/180704_111446/0000/slimMiniAOD_data_MuEle_4.root',
        'file:/eos/cms/store/user/nmangane/DoubleMuon/crab_CMSDAS_Data_analysis_test0/180704_111446/0000/slimMiniAOD_data_MuEle_5.root',
        'file:/eos/cms/store/user/nmangane/DoubleMuon/crab_CMSDAS_Data_analysis_test0/180704_111446/0000/slimMiniAOD_data_MuEle_6.root',
        'file:/eos/cms/store/user/nmangane/DoubleMuon/crab_CMSDAS_Data_analysis_test0/180704_111446/0000/slimMiniAOD_data_MuEle_7.root',
        'file:/eos/cms/store/user/nmangane/DoubleMuon/crab_CMSDAS_Data_analysis_test0/180704_111446/0000/slimMiniAOD_data_MuEle_8.root',
        'file:/eos/cms/store/user/nmangane/DoubleMuon/crab_CMSDAS_Data_analysis_test0/180704_111446/0000/slimMiniAOD_data_MuEle_9.root',
        'file:/eos/cms/store/user/nmangane/DoubleMuon/crab_CMSDAS_Data_analysis_test0/180704_111446/0000/slimMiniAOD_data_MuEle_10.root',
        'file:/eos/cms/store/user/nmangane/DoubleMuon/crab_CMSDAS_Data_analysis_test0/180704_111446/0000/slimMiniAOD_data_MuEle_11.root',
        'file:/eos/cms/store/user/nmangane/DoubleMuon/crab_CMSDAS_Data_analysis_test0/180704_111446/0000/slimMiniAOD_data_MuEle_12.root',
        'file:/eos/cms/store/user/nmangane/DoubleMuon/crab_CMSDAS_Data_analysis_test0/180704_111446/0000/slimMiniAOD_data_MuEle_13.root',
        'file:/eos/cms/store/user/nmangane/DoubleMuon/crab_CMSDAS_Data_analysis_test0/180704_111446/0000/slimMiniAOD_data_MuEle_14.root',
        'file:/eos/cms/store/user/nmangane/DoubleMuon/crab_CMSDAS_Data_analysis_test0/180704_111446/0000/slimMiniAOD_data_MuEle_15.root',
        'file:/eos/cms/store/user/nmangane/DoubleMuon/crab_CMSDAS_Data_analysis_test0/180704_111446/0000/slimMiniAOD_data_MuEle_16.root',
        'file:/eos/cms/store/user/nmangane/DoubleMuon/crab_CMSDAS_Data_analysis_test0/180704_111446/0000/slimMiniAOD_data_MuEle_17.root',
        'file:/eos/cms/store/user/nmangane/DoubleMuon/crab_CMSDAS_Data_analysis_test0/180704_111446/0000/slimMiniAOD_data_MuEle_18.root',
        'file:/eos/cms/store/user/nmangane/DoubleMuon/crab_CMSDAS_Data_analysis_test0/180704_111446/0000/slimMiniAOD_data_MuEle_19.root',
        'file:/eos/cms/store/user/nmangane/DoubleMuon/crab_CMSDAS_Data_analysis_test0/180704_111446/0000/slimMiniAOD_data_MuEle_20.root',
        'file:/eos/cms/store/user/nmangane/DoubleMuon/crab_CMSDAS_Data_analysis_test0/180704_111446/0000/slimMiniAOD_data_MuEle_21.root',
        'file:/eos/cms/store/user/nmangane/DoubleMuon/crab_CMSDAS_Data_analysis_test0/180704_111446/0000/slimMiniAOD_data_MuEle_22.root',
        'file:/eos/cms/store/user/nmangane/DoubleMuon/crab_CMSDAS_Data_analysis_test0/180704_111446/0000/slimMiniAOD_data_MuEle_23.root',
        'file:/eos/cms/store/user/nmangane/DoubleMuon/crab_CMSDAS_Data_analysis_test0/180704_111446/0000/slimMiniAOD_data_MuEle_24.root',
        'file:/eos/cms/store/user/nmangane/DoubleMuon/crab_CMSDAS_Data_analysis_test0/180704_111446/0000/slimMiniAOD_data_MuEle_25.root',
        'file:/eos/cms/store/user/nmangane/DoubleMuon/crab_CMSDAS_Data_analysis_test0/180704_111446/0000/slimMiniAOD_data_MuEle_26.root',
        'file:/eos/cms/store/user/nmangane/DoubleMuon/crab_CMSDAS_Data_analysis_test0/180704_111446/0000/slimMiniAOD_data_MuEle_27.root',
        'file:/eos/cms/store/user/nmangane/DoubleMuon/crab_CMSDAS_Data_analysis_test0/180704_111446/0000/slimMiniAOD_data_MuEle_28.root',
        'file:/eos/cms/store/user/nmangane/DoubleMuon/crab_CMSDAS_Data_analysis_test0/180704_111446/0000/slimMiniAOD_data_MuEle_29.root',
        'file:/eos/cms/store/user/nmangane/DoubleMuon/crab_CMSDAS_Data_analysis_test0/180704_111446/0000/slimMiniAOD_data_MuEle_30.root',
        'file:/eos/cms/store/user/nmangane/DoubleMuon/crab_CMSDAS_Data_analysis_test0/180704_111446/0000/slimMiniAOD_data_MuEle_31.root',
  )
)

#process.MessageLogger = cms.Service("MessageLogger")
process.MessageLogger = cms.Service("MessageLogger",
    destinations   = cms.untracked.vstring(
        'detailedInfo',
        'critical',
        'cerr'),
    critical       = cms.untracked.PSet(
        threshold = cms.untracked.string('ERROR')
          ),
    detailedInfo   = cms.untracked.PSet(
        threshold = cms.untracked.string('INFO')
        ),
    cerr           = cms.untracked.PSet(
        threshold = cms.untracked.string('WARNING')
        )
)

process.maxEvents = cms.untracked.PSet( 
    input = cms.untracked.int32(-1) 
)

process.analyzeBasicPat = cms.EDAnalyzer("MyZPeakAnalyzer",
  muonSrc = cms.untracked.InputTag("slimmedMuons"),                  
  elecSrc = cms.untracked.InputTag("slimmedElectrons"),
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('myZPeakCRAB.root')
                                   )

process.p = cms.Path(process.analyzeBasicPat)
