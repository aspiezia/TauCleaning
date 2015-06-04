import FWCore.ParameterSet.Config as cms

#from ExoDiBosonResonances.PATtupleProduction.PAT_ca8jets_cff import selectedPatJetsCA8CHSpruned

from ExoDiBosonResonances.PATtupleProduction.PAT_ca8jets_cff import *

import ExoDiBosonResonances.ProducerCA8JetsForBoostedTaus.producerca8jetsforboostedtaus_cfi

import ExoDiBosonResonances.ProducerCA8CHSwithQJetsForBoostedTaus.producerca8chswithqjetsforboostedtaus_cfi

selectedPatJetsCA8CHSprunedForBoostedTaus = ExoDiBosonResonances.ProducerCA8JetsForBoostedTaus.producerca8jetsforboostedtaus_cfi.demo.clone(
                                             jets = cms.InputTag("selectedPatJetsCA8CHSpruned"),
)

selectedPatJetsCA8CHSwithQJetsForBoostedTaus = ExoDiBosonResonances.ProducerCA8CHSwithQJetsForBoostedTaus.producerca8chswithqjetsforboostedtaus_cfi.demo.clone(
                                             jets = cms.InputTag("selectedPatJetsCA8CHSwithQjets"),
)

boostedTaus = cms.Sequence(selectedPatJetsCA8CHSprunedForBoostedTaus + selectedPatJetsCA8CHSwithQJetsForBoostedTaus)

