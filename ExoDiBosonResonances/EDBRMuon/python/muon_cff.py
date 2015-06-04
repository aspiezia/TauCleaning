#from CMGTools.Common.muon_cff import *
from ExoDiBosonResonances.EDBRMuon.factories.cmgMuon_cfi import *


genSelectorZDaughterMu = cms.EDFilter("GenParticleSelector",
    src = cms.InputTag("genParticles"),
    cut = cms.string(' (abs(pdgId)==13 )&& abs(mother.pdgId)==23 ')
)
genSelectorWDaughterMu = cms.EDFilter("GenParticleSelector",
    src = cms.InputTag("genParticles"),
    cut = cms.string(' (abs(pdgId)==13 )&& abs(mother.pdgId)==24 ')
)

selectedPatMuons  = cms.EDProducer("PATMuonCleaner",
                                   src = cms.InputTag("patMuonsWithTrigger"), #selectedPatMuonsAK5
                                   preselection = cms.string(''),
                                   checkOverlaps = cms.PSet( genLeptons = cms.PSet( src = cms.InputTag("genSelectorZDaughterMu"),
                                                                                    algorithm = cms.string("byDeltaR"),
                                                                                    preselection        = cms.string(""),  # don't preselect the muons
                                                                                    deltaR              = cms.double(0.3),
                                                                                    checkRecoComponents = cms.bool(False), # don't check if they share some AOD object ref
                                                                                    pairCut             = cms.string(""),
                                                                                    requireNoOverlaps = cms.bool(False), # overlaps don't cause the electron to be discare
                                                                                    ),
                                                             genLeptonsW = cms.PSet( src = cms.InputTag("genSelectorWDaughterMu"),
                                                                                    algorithm = cms.string("byDeltaR"),
                                                                                    preselection        = cms.string(""),  # don't preselect the muons
                                                                                    deltaR              = cms.double(0.3),
                                                                                    checkRecoComponents = cms.bool(False), # don't check if they share some AOD object ref
                                                                                    pairCut             = cms.string(""),
                                                                                    requireNoOverlaps = cms.bool(False), # overlaps don't cause the electron to be discare
                                                                                    ),
                                                             ),                                      
                                   finalCut = cms.string(''),
                                   )



# cmgMuon.cfg.inputCollection="selectedPatMuons"
# #cmgMuon.cfg.trackType=??? at this moment always glb track used
# #cmgMuon.cfg.muonIDType=??? at this moment TMLastStationLoose ?why??
# cmgMuon.cuts = cms.PSet(
#     vbtfmuon = vbtfmuon.clone(),
#     kinematics = muKinematics.clone(),
#     genLepton = cms.string("sourcePtr().get().hasOverlaps('genLeptons')")
#     )
 
muonSequence = cms.Sequence(genSelectorZDaughterMu + genSelectorWDaughterMu + selectedPatMuons + cmgMuon)


