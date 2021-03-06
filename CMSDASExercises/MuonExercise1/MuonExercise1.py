import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('Demo')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit = cms.untracked.int32(-1)
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False),
SkipEvent = cms.untracked.vstring('ProductNotFound'))

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring('file:///eos/uscms/store/user/cmsdas/2015/SHORT_EXERCISES/MUONS/dymm.root'))

process.demo = cms.EDAnalyzer('MuonExercise1')

process.printEventNumber = cms.OutputModule("AsciiOutputModule")

process.TFileService = cms.Service("TFileService",
          fileName = cms.string('histos.root')
)

process.p = cms.Path(process.demo)
