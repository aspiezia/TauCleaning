import FWCore.ParameterSet.Config as cms
#from CMGTools.Common.factories.cmgMuon_cfi import cmgMuon

dimuonVjetFactory = cms.PSet(
       leg1Collection = cms.InputTag("ZmmCand"),
       leg2Collection = cms.InputTag("jetIDMerged"),
       metCollection = cms.InputTag("")
)

#from CMGTools.Common.selections.zmumu_cfi import zmumu
cmgDiMuonVJet = cms.EDFilter(
                   "DiMuonSingleJetPOProducer",
                   cfg = dimuonVjetFactory.clone(),
                   #overlap cut makes sure that Z->J is not a duplicate of Z->ll (no PAT-cleaning performed so far)
                   cuts = cms.PSet(preSel = cms.PSet( noOverlap = cms.string(" deltaR(leg1.leg1.eta,leg1.leg1.phi,leg2.eta,leg2.phi) > 0.8 &&"
                                                                            +"deltaR(leg1.leg2.eta,leg1.leg2.phi,leg2.eta,leg2.phi) > 0.8 "),
                                                      kinematics = cms.PSet( mass = cms.string("mass > 180")
                                                                            ),#XXXX
                                                      mergedJetID = cms.string("leg2.getSelection(\"cuts_looseJetId\") && leg2.getSelection(\"cuts_TOBTECjetsId\")")
                                                       )
                                   
##                     btags      = cms.PSet( btag0= cms.string("leg2.getSelection(\"cuts_btags_btag0\") "),
##                                            btag1= cms.string("leg2.getSelection(\"cuts_btags_btag1\") "),
##                                            btag2= cms.string("leg2.getSelection(\"cuts_btags_btag2\") "),
##                                           ),
                    
                    )
       )
    
