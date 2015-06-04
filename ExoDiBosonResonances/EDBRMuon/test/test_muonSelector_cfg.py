import FWCore.ParameterSet.Config as cms
import sys

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring("file:newPatTuple_ZZ_1500_c1.root")
                            )

process.t = cms.EDFilter("HPTMMuonTest",
                         muons = cms.InputTag("patMuonsWithTrigger")
                         )

process.p = cms.Path(process.t)
                                              
