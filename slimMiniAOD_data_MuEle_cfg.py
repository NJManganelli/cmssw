## import skeleton process
import FWCore.ParameterSet.Config as cms

process = cms.Process("DAS")

#process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger = cms.Service("MessageLogger",
                     destinations       =  cms.untracked.vstring('messages_slimMiniAOD_data.txt'),
                     debugModules   = cms.untracked.vstring('*'),
                     #messages          = cms.untracked.PSet(
                      #                          threshold =  cms.untracked.vstring('DEBUG')
                       #                             ),
                    #suppressDebug  = cms.untracked.vstring('wordyModule'),
                    suppressInfo       = cms.untracked.vstring('DAS'),
                    #suppressWarning= ms.untracked.vstring('cryWolfModule')
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
      'root://cmseos.fnal.gov//store/user/cmsdas/2017/pre_exercises/DoubleMuon.root'
    )
)


process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('slimMiniAOD_data_MuEle.root'),
    outputCommands = cms.untracked.vstring(['drop *', 'keep *_slimmedMuons__*', 'keep *_slimmedElectrons__*'])
)

process.end = cms.EndPath(process.out)
