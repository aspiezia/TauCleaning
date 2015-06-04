import FWCore.ParameterSet.Config as cms

tunePmuons = cms.EDProducer("HPTMNewTunePProducer",
                            muLabel = cms.InputTag("muons")
                            )

muonTrackError = cms.EDProducer("HPTMTrackErrorProducer",
                                muLabel = cms.InputTag("muons")
                                )
