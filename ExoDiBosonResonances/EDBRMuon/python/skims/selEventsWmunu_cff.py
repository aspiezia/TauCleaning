import FWCore.ParameterSet.Config as cms 

WmunuCand = cms.EDFilter(
    "CmgWmunuSelector",
    src = cms.InputTag("cmgWmunuEDBR"),
    cut = cms.string( "getSelection(\"cuts_Wmunu_kinematics\") && getSelection(\"cuts_Wmunu_quality\") && getSelection(\"cuts_Wmunu_validation\") " )
    )   

selectedWmunuCandFilter = cms.EDFilter("CandViewCountFilter",
   src = cms.InputTag('WmunuCand'),
   minNumber = cms.uint32(1)
 )

selectedWmunuSequence = cms.Sequence(WmunuCand+selectedWmunuCandFilter)
