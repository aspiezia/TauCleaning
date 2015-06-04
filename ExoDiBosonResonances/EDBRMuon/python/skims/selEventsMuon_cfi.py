import FWCore.ParameterSet.Config as cms

muonPreselNoIso = cms.EDFilter(
    "CmgMuonSelector",
    src = cms.InputTag("cmgMuon"),
    cut = cms.string( "getSelection(\"cuts_kinematics\") && getSelection(\"cuts_HPTTKmuon\") " )
    )

muonPreselLoose = cms.EDFilter(
    "CmgMuonSelector",
    src = cms.InputTag("cmgMuon"),
    cut = cms.string( "getSelection(\"cuts_kinematics\") && getSelection(\"cuts_HPTmuonLoose\") " )
    ) 

## isolation cut done by a special module for including PU-corrections
## https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Muon_Isolation
muonPresel = cms.EDProducer(
    "MuonIsoCorrector",
    src=cms.InputTag('muonPreselNoIso'), #electronPreselNoIso
    rho=cms.InputTag('kt6PFJetsCentralNeutral:rho'),
    relisocut=cms.double(0.12),
    EAetaBins=cms.vdouble(0.0, 1.0, 1.5, 2.0, 2.2, 2.3, 2.4, 999.0),
    EAvals=cms.vdouble(0.674, 0.565, 0.442, 0.515, 0.821, 0.660, 0.0),
    verbose=cms.bool(False)
    )

selectedMuonCandFilter = cms.EDFilter("CandViewCountFilter",
   src = cms.InputTag('muonPreselNoIso'),
   minNumber = cms.uint32(0)
 )

selectedMuonLooseCandFilter = cms.EDFilter("CandViewCountFilter",
   src = cms.InputTag('muonPreselLoose'),
   minNumber = cms.uint32(0)
#   maxNumber = cms.uint32(1)
 )

selectedMuonSequence = cms.Sequence(muonPreselNoIso + selectedMuonCandFilter 
                                       +  muonPreselLoose + selectedMuonLooseCandFilter   )

