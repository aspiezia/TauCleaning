import FWCore.ParameterSet.Config as cms 

WmunuSignleJetFactory = cms.PSet(
       leg1Collection = cms.InputTag("WmunuCand"),
       leg2Collection = cms.InputTag("jetIDMerged"),
       metCollection = cms.InputTag("")
)

cmgWmunuSingleJet = cms.EDFilter(
                   "WmunuSingleJetPOProducer",
                   cfg = WmunuSignleJetFactory.clone(),
                   #overlap cut makes sure that Z->J is not a duplicate of Z->ll (no PAT-cleaning performed so far)
                   cuts = cms.PSet( preSel = cms.PSet( noOverlap  = cms.string(" deltaR(leg1.leg1.eta,leg1.leg1.phi,leg2.eta,leg2.phi) > 0.8 " ),
                                                                         #   +"deltaR(leg1.leg2.eta,leg1.leg2.phi,leg2.eta,leg2.phi) > 0.8 "),
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


WmunuSingleJetEDBRFactory = cms.PSet(
       inputs = cms.InputTag("cmgWmunuSingleJet"),
       vbftag = cms.InputTag("VBFPairs"),
       overlapcut = cms.string(" deltaR(vbfptr.leg1.eta,vbfptr.leg1.phi,leg2.eta,leg2.phi) > 0.8 &&"
                               +"deltaR(vbfptr.leg2.eta,vbfptr.leg2.phi,leg2.eta,leg2.phi) > 0.8 &&"
                               +"deltaR(vbfptr.leg1.eta,vbfptr.leg1.phi,leg1.leg1.sourcePtr.eta,leg1.leg1.sourcePtr.phi) > 0.8 &&"
                               +"deltaR(vbfptr.leg2.eta,vbfptr.leg2.phi,leg1.leg1.sourcePtr.eta,leg1.leg1.sourcePtr.phi) > 0.8 " )
                         #      +"deltaR(vbfptr.leg1.eta,vbfptr.leg1.phi,leg1.leg2.sourcePtr.eta,leg1.leg2.sourcePtr.phi) > 0.8 &&"
                         #      +"deltaR(vbfptr.leg2.eta,vbfptr.leg2.phi,leg1.leg2.sourcePtr.eta,leg1.leg2.sourcePtr.phi) > 0.8 ")
)

from ExoDiBosonResonances.EDBRCommon.selections.vjetmcmatch_cfi import * 
cmgWmunuSingleJetEDBR = cms.EDFilter(
    "WmunuSingleJetEDBRPOProducer",
    cfg = WmunuSingleJetEDBRFactory.clone(),
    cuts = cms.PSet(  genMatch = genMatchVQQ,
                      genMatchZ = genMatchVV 
                     )   
    
    ) 
