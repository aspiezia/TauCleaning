import FWCore.ParameterSet.Config as cms


from ExoDiBosonResonances.EDBRMuon.factories.cmgDiMuon_cfi import *
from ExoDiBosonResonances.EDBRMuon.skims.cmgDiMuonSel_cfi import *
from CMGTools.Common.skims.cmgDiMuonCount_cfi import *
from CMGTools.Common.histograms.cmgDiMuonHistograms_cfi import *

diMuonSequence = cms.Sequence(
    cmgDiMuon  
   + cmgDiMuonSel 
   + cmgDiMuonCount 
###   +cmgDiMuonHistograms
)


###test for Git
