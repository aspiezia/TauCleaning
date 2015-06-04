import FWCore.ParameterSet.Config as cms

from CMGTools.Common.selections.vbtfmuon_cfi import vbtfmuon

pfmuonId2012 = cms.PSet(
    isGlobal = cms.string('isGlobalMuon()'),
    isPF = cms.string('isPF()'),
    #Updated for 2012 with new cut from the Muon POG
    #numberOfValidTrackerHits = cms.string('numberOfValidTrackerHits() > 10'),
    trackerLayersWithMeasurement = cms.string('trackerLayersWithMeasurement() > 5'),
    numberOfValidPixelHits = cms.string('numberOfValidPixelHits() > 0'),
    numberOfValidMuonHits = cms.string('numberOfValidMuonHits() > 0'),
    numberOfMatches = cms.string('numberOfMatches() > 1'),
    normalizedChi2 = cms.string('normalizedChi2() < 10'),
 #   dxy = cms.string('abs(dB(pat::Muon::PV3D)) < 0.2'),
    dxy = cms.string('abs(dxy()) < 0.2'),
    dz = cms.string('abs(dz()) < 0.5')
    )


HPTmuonTKId2012 = cms.PSet(
    isTracker = cms.string('isTrackerMuon()'),
    #Updated for 2012 with new cut from the Muon POG
    #numberOfValidTrackerHits = cms.string('numberOfValidTrackerHits() > 10'),
    #    trackerLayersWithMeasurement = cms.string('trackerLayersWithMeasurement() > 8'),
    isGoodFit = cms.string('sourcePtr().userFloat("muonTrackError") < 0.3'),
    trackerLayersWithMeasurement = cms.string('trackerLayersWithMeasurement() > 5'),
    numberOfValidPixelHits = cms.string('numberOfValidPixelHits() > 0'),
    numberOfMatches = cms.string('numberOfMatchedStations() > 1'),
 #   normalizedChi2 = cms.string('normalizedChi2() < 10'),
 #   dxy = cms.string('abs(dB(pat::Muon::PV3D)) < 0.2'),
    dxy = cms.string('abs(dxy()) < 0.2'),
    dz = cms.string('abs(dz()) < 0.5')
    )


HPTmuonGlobalId2012 = cms.PSet(
    isGlobal = cms.string('isGlobalMuon()'),
    #Updated for 2012 with new cut from the Muon POG
    #numberOfValidTrackerHits = cms.string('numberOfValidTrackerHits() > 10'),
    #    trackerLayersWithMeasurement = cms.string('trackerLayersWithMeasurement() > 8'),
    isGoodFit = cms.string('sourcePtr().userFloat("muonTrackError") < 0.3'),
    trackerLayersWithMeasurement = cms.string('trackerLayersWithMeasurement() > 5'),
    numberOfValidPixelHits = cms.string('numberOfValidPixelHits() > 0'),
    numberOfValidMuonHits = cms.string('numberOfValidMuonHits() > 0'),
    numberOfMatches = cms.string('numberOfMatchedStations() > 1'),
 #   normalizedChi2 = cms.string('normalizedChi2() < 10'),
 #   dxy = cms.string('abs(dB(pat::Muon::PV3D)) < 0.2'),
    dxy = cms.string('abs(dxy()) < 0.2'),
    dz = cms.string('abs(dz()) < 0.5')
    )


HPTmuonLooseId = HPTmuonGlobalId2012.clone()
HPTmuonLooseId.isIsolated=cms.string('sourcePtr().trackIso()/pt() < 0.1')

# HPTmuonLooseId = cms.PSet( ###### now it is HPTmuonGlobalId2012 + isolation
#     isGlobal = cms.string('isGlobalMuon()'),
#     #Updated for 2012 with new cut from the Muon POG
#     #numberOfValidTrackerHits = cms.string('numberOfValidTrackerHits() > 10'),
#     trackerLayersWithMeasurement = cms.string('trackerLayersWithMeasurement() > 8'),
#     numberOfValidPixelHits = cms.string('numberOfValidPixelHits() > 0'),
#     numberOfValidMuonHits = cms.string('numberOfValidMuonHits() > 0'),
#     numberOfMatches = cms.string('numberOfMatchedStations() > 1'),
#  #   normalizedChi2 = cms.string('normalizedChi2() < 10'),
#  #   dxy = cms.string('abs(dB(pat::Muon::PV3D)) < 0.2'),
#     dxy = cms.string('abs(dxy()) < 0.2'),
#     dz = cms.string('abs(dz()) < 0.5'),
#     ##isolation
#     isIsolated=cms.string('sourcePtr().trackIso()/pt() < 0.1')
# )


# HPTmuonLooseId = cms.PSet(  ###### now it is HPTmuonTKId2012
#     isTracker = cms.string('isTrackerMuon()'),
#     #Updated for 2012 with new cut from the Muon POG
#     #numberOfValidTrackerHits = cms.string('numberOfValidTrackerHits() > 10'),
#     trackerLayersWithMeasurement = cms.string('trackerLayersWithMeasurement() > 8'),
#     numberOfValidPixelHits = cms.string('numberOfValidPixelHits() > 0'),
#     numberOfMatches = cms.string('numberOfMatchedStations() > 1'),
#  #   normalizedChi2 = cms.string('normalizedChi2() < 10'),
#  #   dxy = cms.string('abs(dB(pat::Muon::PV3D)) < 0.2'),
#     dxy = cms.string('abs(dxy()) < 0.2'),
#     dz = cms.string('abs(dz()) < 0.5')
#     )



