import FWCore.ParameterSet.Config as cms


selectedEDBRCandFilterMu = cms.EDFilter("CandViewCountFilter",
   src = cms.InputTag('cmgEDBRSelMu'),
   minNumber = cms.uint32(0)
 )


selectedEDBRKinFitCandFilterMu = cms.EDFilter("CandViewCountFilter",
                                      src = cms.InputTag('cmgEDBRSelKinFitMu'),
                                      minNumber = cms.uint32(1)
                                      )

selectedEDBRMergedCandFilterMu = cms.EDFilter("CandViewCountFilter",
   src = cms.InputTag('cmgEDBRMergedSelMu'),
   minNumber = cms.uint32(1)
 )

