import FWCore.ParameterSet.Config as cms

muKinematics = cms.PSet(
    pt = cms.string('pt() > 20'),
    eta = cms.string('abs(eta()) < 2.4'),
    phi = cms.string('abs(phi()) < 3.2')
    )

