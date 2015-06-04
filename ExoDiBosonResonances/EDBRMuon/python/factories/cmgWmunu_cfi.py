import FWCore.ParameterSet.Config as cms 


WmunuFactory = cms.PSet(
       leg1Collection = cms.InputTag("muonPreselNoIso"),
       leg2Collection = cms.InputTag("MuNeutrinoPresel"),
	   metCollection = cms.InputTag("")
)
from ExoDiBosonResonances.EDBRMuon.selections.Wmunu_cfi import Wmunu_kinematics
from ExoDiBosonResonances.EDBRMuon.selections.Wmunu_cfi import Wmunu_quality
from ExoDiBosonResonances.EDBRMuon.selections.Wmunu_cfi import Wmunu_validation

cmgWmunuEDBR = cms.EDFilter(
    "WmunuPOProducer",
    cfg = WmunuFactory.clone(),
    cuts = cms.PSet(
       Wmunu_kinematics = Wmunu_kinematics.clone(),
       Wmunu_quality = Wmunu_quality.clone(),
       Wmunu_validation =  Wmunu_validation.clone(),
       genP = cms.PSet(genP = cms.string("leg1.sourcePtr.get.hasOverlaps(\"genLeptonsW\")")
       ) 
    )
)
