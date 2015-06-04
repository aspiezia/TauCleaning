import FWCore.ParameterSet.Config as cms

from ExoDiBosonResonances.EDBRMuon.factories.cmgWmunuDiJet_cfi import * 
from ExoDiBosonResonances.EDBRMuon.factories.cmgWmunuSingleJet_cfi import * 

from ExoDiBosonResonances.EDBRMuon.skims.cmgEDBRSel_Wmu_cff import *
from ExoDiBosonResonances.EDBRMuon.skims.selEventsEDBR_cfi import *
from ExoDiBosonResonances.EDBRMuon.HLTWeights_Wmu_cff import *
from ExoDiBosonResonances.EDBRMuon.BTagWeights_Wmu_cff import *

cmgEDBRKinFitWeighted2012AMu = cms.EDProducer("WmunuDiJetEDBRWeightAdder",
                                        src=cms.InputTag("BTagWeightsKinFitMu"),
                                        weight=cms.InputTag("PUWeights2012A"),
                                        )   
cmgEDBRWeighted2012AMu = cms.EDProducer("WmunuDiJetEDBRWeightAdder",
                                      src=cms.InputTag("BTagWeightsMu"),
                                      weight=cms.InputTag("PUWeights2012A"),
                                      )   

cmgEDBRMergedWeighted2012AMu = cms.EDProducer("WmunuSingleJetEDBRWeightAdder",
                                            src=cms.InputTag("BTagWeightsMergedMu"),
                                            weight=cms.InputTag("PUWeights2012A")    
                                            )   

########

cmgEDBRKinFitWeighted2012BMu = cms.EDProducer("WmunuDiJetEDBRWeightAdder",
                                            src=cms.InputTag("cmgEDBRKinFitWeighted2012AMu"),
                                            weight=cms.InputTag("PUWeights2012B")
                                            )   

cmgEDBRWeighted2012BMu = cms.EDProducer("WmunuDiJetEDBRWeightAdder",
                                      src=cms.InputTag("cmgEDBRWeighted2012AMu"),
                                      weight=cms.InputTag("PUWeights2012B")
                                      )

cmgEDBRMergedWeighted2012BMu = cms.EDProducer("WmunuSingleJetEDBRWeightAdder",
                                            src=cms.InputTag("cmgEDBRMergedWeighted2012AMu"),
                                            weight=cms.InputTag("PUWeights2012B"),
                                            )
#####
cmgEDBRKinFitWeightedMu = cms.EDProducer("WmunuDiJetEDBRWeightAdder",
                                        src=cms.InputTag("cmgEDBRKinFitWeighted2012BMu"),
                                        weight=cms.InputTag("PUWeights"),
                                        )
cmgEDBRWeightedMu = cms.EDProducer("WmunuDiJetEDBRWeightAdder",
                                  src=cms.InputTag("cmgEDBRWeighted2012BMu"),
                                  weight=cms.InputTag("PUWeights"),
                                  )
cmgEDBRMergedWeightedMu = cms.EDProducer("WmunuSingleJetEDBRWeightAdder",
                                  src=cms.InputTag("cmgEDBRMergedWeighted2012BMu"),
                                  weight=cms.InputTag("PUWeights"),
                                  )

#### add extra-kinematic variables to EDBR candidate (only the one with kin-fit)
cmgEDBRExtraMu = cms.EDProducer("WmunuDiJetEDBRKineAdder",
                                src=cms.InputTag("cmgEDBRKinFitWeightedMu"),
                                noKinFitSrc=cms.InputTag("cmgEDBRWeightedMu"), #leave empty if SingleJet EDBR
                                BTagJets=cms.InputTag("jetAK5"),
                                BTagCleaningTarget=cms.InputTag("jetIDJet")                                
                                )

cmgEDBRMergedExtraMu = cms.EDProducer("WmunuSingleJetEDBRKineAdder",
                                      src=cms.InputTag("cmgEDBRMergedWeightedMu"),
                                      noKinFitSrc=cms.InputTag(""), #leave empty if SingleJet EDBR
                                      BTagJets=cms.InputTag("jetAK5"),
                                      BTagCleaningTarget=cms.InputTag("jetIDJet")                                
            
                                   )

edbrSequenceMVJJ = cms.Sequence(
    cmgWmunuDiJet +
    cmgWmunuDiJetKinFit +

    cmgWmunuDiJetEDBR +
    cmgWmunuDiJetKinFitEDBR +

    HLTWeightsMu +
    HLTWeightsKinFitMu +

	BTagWeightsMu +
	BTagWeightsKinFitMu +

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

edbrSequenceMergedMVJ = cms.Sequence(
    cmgWmunuSingleJet +
    cmgWmunuSingleJetEDBR +

    HLTWeightsMergedMu +

	BTagWeightsMergedMu +

    cmgEDBRMergedWeighted2012AMu +
    cmgEDBRMergedWeighted2012BMu +
    cmgEDBRMergedWeightedMu +
    cmgEDBRMergedExtraMu+
    cmgEDBRMergedSelMu


#    selectedEDBRMergedCandFilterMu 


    )
      
