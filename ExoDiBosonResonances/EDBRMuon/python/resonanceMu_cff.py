import FWCore.ParameterSet.Config as cms

from ExoDiBosonResonances.EDBRMuon.factories.cmgdimuondijet_cfi import * 
from ExoDiBosonResonances.EDBRMuon.factories.cmgdimuondijetEDBR_cfi import *
from ExoDiBosonResonances.EDBRMuon.factories.cmgDiMuonSingleJet_cfi import * 
from ExoDiBosonResonances.EDBRMuon.factories.cmgDiMuonSingleJetEDBR_cfi import *
from ExoDiBosonResonances.EDBRMuon.skims.cmgEDBRSel_cff import *
from ExoDiBosonResonances.EDBRMuon.skims.selEventsEDBR_cfi import *
from ExoDiBosonResonances.EDBRMuon.HLTWeights_cff import *

cmgEDBRKinFitWeighted2012AMu = cms.EDProducer("DiMuonDiJetEDBRWeightAdder",
                                        src=cms.InputTag("HLTWeightsKinFitMu"),
                                        weight=cms.InputTag("PUWeights2012A"),
                                        )
cmgEDBRWeighted2012AMu = cms.EDProducer("DiMuonDiJetEDBRWeightAdder",
                                      src=cms.InputTag("HLTWeightsMu"),
                                      weight=cms.InputTag("PUWeights2012A"),
                                      )

cmgEDBRMergedWeighted2012AMu = cms.EDProducer("DiMuonVJetEDBRWeightAdder",
                                            src=cms.InputTag("HLTWeightsMergedMu"),
                                            weight=cms.InputTag("PUWeights2012A")                                      
                                            )

########

cmgEDBRKinFitWeighted2012BMu = cms.EDProducer("DiMuonDiJetEDBRWeightAdder",
                                            src=cms.InputTag("cmgEDBRKinFitWeighted2012AMu"),
                                            weight=cms.InputTag("PUWeights2012B")
                                            )

cmgEDBRWeighted2012BMu = cms.EDProducer("DiMuonDiJetEDBRWeightAdder",
                                      src=cms.InputTag("cmgEDBRWeighted2012AMu"),
                                      weight=cms.InputTag("PUWeights2012B")
                                      )

cmgEDBRMergedWeighted2012BMu = cms.EDProducer("DiMuonVJetEDBRWeightAdder",
                                            src=cms.InputTag("cmgEDBRMergedWeighted2012AMu"),
                                            weight=cms.InputTag("PUWeights2012B"),
                                            )
#####
cmgEDBRKinFitWeightedMu = cms.EDProducer("DiMuonDiJetEDBRWeightAdder",
                                        src=cms.InputTag("cmgEDBRKinFitWeighted2012BMu"),
                                        weight=cms.InputTag("PUWeights"),
                                        )
cmgEDBRWeightedMu = cms.EDProducer("DiMuonDiJetEDBRWeightAdder",
                                  src=cms.InputTag("cmgEDBRWeighted2012BMu"),
                                  weight=cms.InputTag("PUWeights"),
                                  )
cmgEDBRMergedWeightedMu = cms.EDProducer("DiMuonVJetEDBRWeightAdder",
                                  src=cms.InputTag("cmgEDBRMergedWeighted2012BMu"),
                                  weight=cms.InputTag("PUWeights"),
                                  )

#### add extra-kinematic variables to EDBR candidate (only the one with kin-fit)
cmgEDBRExtraMu = cms.EDProducer("DiMuonDiJetEDBRKineAdder",
                                   src=cms.InputTag("cmgEDBRKinFitWeightedMu"),
                                   noKinFitSrc=cms.InputTag("cmgEDBRWeightedMu"), #leave empty if SingleJet EDBR
                                BTagJets=cms.InputTag(""),#leave empty for ignoring
                                BTagCleaningTarget=cms.InputTag("")                                
            
                                   )

cmgEDBRMergedExtraMu = cms.EDProducer("DiMuonSingleJetEDBRKineAdder",
                                   src=cms.InputTag("cmgEDBRMergedWeightedMu"),
                                   noKinFitSrc=cms.InputTag(""), #leave empty if SingleJet EDBR
                                      BTagJets=cms.InputTag(""),#leave empty for ignoring
                                      BTagCleaningTarget=cms.InputTag("")                                
            
                                   )



edbrSequenceMMJJ = cms.Sequence(
    cmgDiMuonDiJet +
    cmgDiMuonDiJetKinFit +
    
    cmgDiMuonDiJetEDBR +
    cmgDiMuonDiJetKinFitEDBR +

    HLTWeightsMu +
    HLTWeightsKinFitMu +

    cmgEDBRKinFitWeighted2012AMu +
    cmgEDBRWeighted2012AMu +
    cmgEDBRKinFitWeighted2012BMu +
    cmgEDBRWeighted2012BMu +
    cmgEDBRKinFitWeightedMu +
    cmgEDBRWeightedMu +

    cmgEDBRExtraMu +
    
    cmgEDBRSelMu +
    cmgEDBRSelKinFitMu 


#    selectedEDBRKinFitCandFilterMu 
#    selectedEDBRMuCandFilter
)

edbrSequenceMergedMMJ = cms.Sequence(
    cmgDiMuonVJet +
    cmgDiMuonVJetEDBR +

    HLTWeightsMergedMu +

    cmgEDBRMergedWeighted2012AMu +
    cmgEDBRMergedWeighted2012BMu +
    cmgEDBRMergedWeightedMu +
    cmgEDBRMergedExtraMu+
    cmgEDBRMergedSelMu 


#    selectedEDBRMergedCandFilterMu 


    )
