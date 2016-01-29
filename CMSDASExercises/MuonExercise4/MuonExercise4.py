import FWCore.ParameterSet.Config as cms

process = cms.Process("Muons")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
## signal
        'file:/eos/uscms/store/user/cmsdas/2016/SHORT_EXERCISES/Muons/dymm.root'
## background
#         'file:/eos/uscms/store/user/cmsdas/2016/SHORT_EXERCISES/Muons/ttjets.root'
    )
)

process.TFileService = cms.Service("TFileService",
 fileName = cms.string("histos4.root")
)

process.demo = cms.EDAnalyzer("MuonExercise4",
     vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
     muons = cms.InputTag("slimmedMuons"),
     bits = cms.InputTag("TriggerResults","","HLT"),
     prescales = cms.InputTag("patTrigger"),
     objects = cms.InputTag("selectedPatTrigger"),
     packed = cms.InputTag("packedGenParticles"),
     pruned = cms.InputTag("prunedGenParticles")
)

process.p = cms.Path(process.demo)
