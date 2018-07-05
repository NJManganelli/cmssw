import FWCore.ParameterSet.Config as cms

process = cms.Process("FWLitePlots")

process.MuonAnalyzer = cms.PSet(
    ## common input for wrapped analyzers
    fileNames   = cms.vstring(
      #'root://cmseos.fnal.gov//store/user/cmsdas/2018/pre_exercises/fourth_set/slimMiniAOD_data_MuEle_1.root',
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

  ),


    outputFile  = cms.string('myZPeakCRAB_fwlite.root'),## mandatory
    outputEvery  = cms.uint32(15000),
    maxEvents   = cms.int32(-1),                      ## optional
    ##reportAfter = cms.uint32(100),                   ## optional
    ## input specific for this analyzer
    muons = cms.InputTag('slimmedMuons')

)
