#include "ExoDiBosonResonances/EDBRMuon/interface/HPTMMuonSelector.h"

std::bitset<8>
hptm::MuonSelector::muonBits(const reco::Muon& recoMu, const hptm::MuonSelector::Point& vertex, hptm::MuonIDType idType){
  bool isGlobal = false;
  bool isTracker = false;
  bool muonChamberHit = false;
  bool matchedStations = false;
  bool dBCut = false;
  bool longiCut = false;
  bool pixelHit = false;
  bool trackerLayers = false;

  isGlobal = recoMu.isGlobalMuon();
  isTracker = recoMu.isTrackerMuon();
  matchedStations = (recoMu.numberOfMatchedStations() > 1);

  // The HIGHPT_OLD Muon ID will explicitly ask for the
  // global track, the best track and the inner track separately.
  if(idType == HIGHPT_OLD) {
    const reco::TrackRef& globalTrackRef = recoMu.globalTrack();
    if(globalTrackRef.isNonnull()) {
      muonChamberHit = (globalTrackRef->hitPattern().numberOfValidMuonHits() > 0);
    }
    
    const reco::TrackRef& bestTrackRef = recoMu.muonBestTrack();
    const reco::TrackRef& innerTrackRef = recoMu.innerTrack();
    const reco::TrackRef& trackRef = recoMu.track();

    if(bestTrackRef.isNonnull()) {
      dBCut    = (fabs(bestTrackRef->dxy(vertex)) < 0.2);  
      longiCut = (fabs(bestTrackRef->dz(vertex)) < 0.5);
    }
    else if(innerTrackRef.isNonnull()) {
      dBCut    = (fabs(innerTrackRef->dxy(vertex)) < 0.2);  
      longiCut = (fabs(innerTrackRef->dz(vertex)) < 0.5);
    }

    if(innerTrackRef.isNonnull())
    pixelHit = (innerTrackRef->hitPattern().numberOfValidPixelHits() > 0);
    
    if(trackRef.isNonnull())
      trackerLayers = (trackRef->hitPattern().trackerLayersWithMeasurement() > 8);
  }

  // The TRACKER Muon ID will just ask for the inner track
  // (since there may be no other track).
  if(idType == TRACKER) {
    const reco::TrackRef& innerTrackRef = recoMu.innerTrack();
    if(innerTrackRef.isNonnull()) {
      dBCut         = (fabs(innerTrackRef->dxy(vertex)) < 0.2);
      longiCut      = (fabs(innerTrackRef->dz(vertex)) < 0.5);
      pixelHit      = (innerTrackRef->hitPattern().numberOfValidPixelHits() > 0);
      trackerLayers = (innerTrackRef->hitPattern().trackerLayersWithMeasurement() > 8);
    }
  }
    
  std::bitset<8> result;
  result.reset(); // everything get's set to false
  
  result[0] = isGlobal;
  result[1] = isTracker;
  result[2] = muonChamberHit;
  result[3] = matchedStations;
  result[4] = dBCut;
  result[5] = longiCut;
  result[6] = pixelHit;
  result[7] = trackerLayers;
  
  return result;
}

bool
hptm::MuonSelector::checkMuonID(const reco::Muon& recoMu, const hptm::MuonSelector::Point& vertex, hptm::MuonIDType idType){

  std::bitset<8> result = muonBits(recoMu, vertex, idType);
  
  bool isGlobal = result[0];
  bool isTracker = result[1];
  bool muonChamberHit = result[2];
  bool matchedStations = result[3];
  bool dBCut = result[4];
  bool longiCut = result[5];
  bool pixelHit = result[6];
  bool trackerLayers = result[7];

  bool passed = false;
  if(idType == HIGHPT_OLD)
    passed = (isGlobal and muonChamberHit and matchedStations and dBCut and longiCut and pixelHit and trackerLayers);
  if(idType == TRACKER)
    passed = (isTracker and matchedStations and dBCut and longiCut and pixelHit and trackerLayers);
  
  return passed;
}
