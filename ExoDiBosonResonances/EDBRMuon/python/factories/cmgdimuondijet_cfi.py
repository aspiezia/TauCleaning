import FWCore.ParameterSet.Config as cms
#from CMGTools.Common.factories.cmgMuon_cfi import cmgMuon

dimuondijetFactory = cms.PSet(
       leg1Collection = cms.InputTag("ZmmCand"),
       leg2Collection = cms.InputTag("ZjjCand"),
       metCollection = cms.InputTag("")
)

dimuondijetKinFitFactory = dimuondijetFactory.clone()
dimuondijetKinFitFactory.leg2Collection = cms.InputTag("cmgDiJetKinFit")

#from CMGTools.Common.selections.zmumu_cfi import zmumu
cmgDiMuonDiJet = cms.EDFilter(
    "DiMuonDiJetPOProducer",
    cfg = dimuondijetFactory.clone(),
    cuts = cms.PSet(preSel = cms.PSet( noOverlap = cms.string(" deltaR(leg1.leg1.eta,leg1.leg1.phi,leg2.leg1.sourcePtr.eta,leg2.leg1.sourcePtr.phi) > 0.8 &&"
                                                             +"deltaR(leg1.leg2.eta,leg1.leg2.phi,leg2.leg1.sourcePtr.eta,leg2.leg1.sourcePtr.phi) > 0.8 &&"
                                                             +"deltaR(leg1.leg1.eta,leg1.leg1.phi,leg2.leg2.sourcePtr.eta,leg2.leg2.sourcePtr.phi) > 0.8 &&"
                                                             +"deltaR(leg1.leg2.eta,leg1.leg2.phi,leg2.leg2.sourcePtr.eta,leg2.leg2.sourcePtr.phi) > 0.8 "),
                                       kinematics = cms.PSet( mass = cms.string("mass > 180")),#XXXX
                                       mergedJetID = cms.string("leg2.leg1.getSelection(\"cuts_looseJetId\") && leg2.leg1.getSelection(\"cuts_TOBTECjetsId\") && leg2.leg2.getSelection(\"cuts_looseJetId\") && leg2.leg2.getSelection(\"cuts_TOBTECjetsId\")")

                                        )
                   

##### NO B-TAGGING for the time being                    
###                    btags      = cms.PSet( btag0= cms.string("leg2.getSelection(\"cuts_btags_btag0\") "),
###                                           btag1= cms.string("leg2.getSelection(\"cuts_btags_btag1\") "),
###                                           btag2= cms.string("leg2.getSelection(\"cuts_btags_btag2\") "),
###                    ),
                    )
    )

cmgDiMuonDiJetKinFit=cmgDiMuonDiJet.clone()
cmgDiMuonDiJetKinFit.cfg=dimuondijetKinFitFactory.clone()


