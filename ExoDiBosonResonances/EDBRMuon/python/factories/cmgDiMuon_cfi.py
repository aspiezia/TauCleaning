import FWCore.ParameterSet.Config as cms
#from ExoDiBosonResonances.EDBRMuon.factories.cmgMuon_cfi import cmgMuon

dimuonFactory = cms.PSet(
       leg1Collection = cms.InputTag("muonPreselNoIso"),
       leg2Collection = cms.InputTag("muonPreselNoIso"),
       metCollection = cms.InputTag("")
)
from ExoDiBosonResonances.EDBRMuon.selections.zmumu_cfi import zmumu, HPTmuGlob

cmgDiMuon = cms.EDFilter(
    "DiMuonPOProducer",
    cfg = dimuonFactory.clone(),
    cuts = cms.PSet(
      zmumu = zmumu.clone(),
      HPTmuGlob = HPTmuGlob.clone(),
      charge = cms.PSet( charge = cms.string("leg1.charge != leg2.charge") ),#XXX
      genP = cms.PSet(genP = cms.string("leg1.sourcePtr.get.hasOverlaps(\"genLeptons\") && leg2.sourcePtr.get.hasOverlaps(\"genLeptons\")"))
      ),
    )


