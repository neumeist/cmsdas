import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
#        'file:///eos/uscms/store/user/cmsdas/2016/SHORT_EXERCISES/Muons/dymm.root'
        'file:///eos/uscms/store/user/cmsdas/2016/SHORT_EXERCISES/Muons/qcd.root'
#        "file:///eos/uscms/store/user/cmsdas/2016/SHORT_EXERCISES/Muons/ttjets.root"
    )
)
process.demo = cms.EDAnalyzer('MuonExercise3',
                      muonInputTag_ = cms.InputTag("slimmedMuons"),
                      genInputTag_  = cms.InputTag("packedGenParticles"),
                      vertexInputTag_ = cms.InputTag("offlineSlimmedPrimaryVertices"))

process.TFileService = cms.Service("TFileService", fileName = cms.string("histos_qcd.root"))

process.p = cms.Path(process.demo)
