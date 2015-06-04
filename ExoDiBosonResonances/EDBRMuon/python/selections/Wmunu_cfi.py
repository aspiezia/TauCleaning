import FWCore.ParameterSet.Config as cms 

Wmunu_kinematics = cms.PSet(
    #mass = cms.string('mass() >= 77.4 && mass() < 83.4'),
    pt = cms.string('leg1().pt() > 50 && leg2().pt() > 40 '),
	eta = cms.string('leg1().eta()<2.1&&leg1().eta()>-2.1')
)

Wmunu_quality = cms.PSet(
    isGB=cms.string('leg1().getSelection(\"cuts_HPTGBmuon\")'),
    isIsolated=cms.string('leg1.sourcePtr.trackIso/leg1.pt < 0.1')
    )
Wmunu_validation  = cms.PSet(
    ptmatch = cms.string('(leg1().pt()-leg2().getleppt())/leg1().pt()<0.01&&(leg1().pt()-leg2().getleppt())/leg1().pt()>-0.01'),
	etamatch = cms.string('(leg1().eta()-leg2().getlepeta())/leg1().eta()<0.01&&(leg1().eta()-leg2().getlepeta())/leg1().eta()>-0.01')
) 
