import FWCore.ParameterSet.Config as cms

process = cms.Process("COPY")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#        'root://cmseos.fnal.gov//store/user/cmsdas/2017/pre_exercises/DYJetsToLL.root'
        'root://cmseos.fnal.gov//store/user/cmsdas/2017/pre_exercises/DoubleMuon.root'
    )
)

#process.copyAll = cms.OutputModule("PoolOutputModule", fileName = cms.untracked.string("DYJetsToLL_n100.root") )
process.copyAll = cms.OutputModule("PoolOutputModule", fileName = cms.untracked.string("DoubleMuon_n100.root") )

process.out = cms.EndPath(process.copyAll)
