import FWCore.ParameterSet.Config as cms


cmgEDBRSelMu = cms.EDFilter(
    "CmgDiMuonDiJetEDBRSelector",
    src = cms.InputTag( "cmgEDBRWeightedMu" ), #without extras
#     cut = cms.string( "getSelection(\"cuts_fullSel_noOverlap\") && getSelection(\"cuts_fullSel_kinematics\")" )
   cut = cms.string( "getSelection(\"cuts_preSel\")" )
  )
cmgEDBRSelKinFitMu = cms.EDFilter(
    "CmgDiMuonDiJetEDBRSelector",
    src = cms.InputTag( "cmgEDBRExtraMu" ),
 #   cut = cms.string( "getSelection(\"cuts_fullSel_noOverlap\") && getSelection(\"cuts_fullSel_kinematics\")" )
   cut = cms.string( "getSelection(\"cuts_preSel\")" )
    )


cmgEDBRMergedSelMu = cms.EDFilter(
    "CmgDiMuonSingleJetEDBRSelector",
    src = cms.InputTag( "cmgEDBRMergedExtraMu" ),
#   cut = cms.string( "getSelection(\"cuts_fullSel_noOverlap\") && getSelection(\"cuts_fullSel_kinematics\")" )
   cut = cms.string( "getSelection(\"cuts_preSel\")" )
    )
