import FWCore.ParameterSet.Config as cms

zmumu = cms.PSet(
    mass = cms.string('mass() >= 70 && mass() < 110'),
    leg1_kinematics = cms.PSet(
                               pt = cms.string('(leg1().pt() > 20 && leg2().pt() > 40)||(leg1().pt() > 40 && leg2().pt() > 20)'),
                               ),
    ### Thiago: we have to put the isolation here because we have to modify it if the two muons are close by...
    leg1_isolation = cms.PSet(
        isol = cms.string('(((deltaR(leg1.eta,leg1.phi,leg2.eta,leg2.phi) > 0.3) && (((leg1.sourcePtr.trackIso)/leg1.pt) < 0.1))) || (((deltaR(leg1.eta,leg1.phi,leg2.eta,leg2.phi) < 0.3) && (((leg1.sourcePtr.trackIso - leg2.pt)/leg1.pt) < 0.1)))')
        ),
    leg2_isolation = cms.PSet(
        isol = cms.string('(((deltaR(leg1.eta,leg1.phi,leg2.eta,leg2.phi) > 0.3) && (((leg2.sourcePtr.trackIso)/leg2.pt) < 0.1))) || (((deltaR(leg1.eta,leg1.phi,leg2.eta,leg2.phi) < 0.3) && (((leg2.sourcePtr.trackIso - leg1.pt)/leg2.pt) < 0.1)))')
        ),
    )

HPTmuGlob = cms.PSet(
    isGB=cms.string('leg1().getSelection(\"cuts_HPTGBmuon\") || leg2().getSelection(\"cuts_HPTGBmuon\")')
    )
