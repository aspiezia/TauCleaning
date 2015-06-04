// -*- C++ -*-
//
// Package:    HPTMNewTunePAnalyzer
// Class:      HPTMNewTunePAnalyzer
// 
/**\class HPTMNewTunePAnalyzer HPTMNewTunePAnalyzer.cc ExoDiBosonResonances/HPTMNewTunePAnalyzer/src/HPTMNewTunePAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Thiago Tomei Fernandez,40 2-B15,+41227671625,
//         Created:  Tue May 21 23:24:55 CEST 2013
// $Id: HPTMNewTunePAnalyzer.cc,v 1.2 2013/05/24 15:31:13 tomei Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
//#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/MuonReco/interface/MuonCocktails.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "TH1F.h"
//
// class declaration
//

/// This is a poor-mans-backporting from CMSSW_5_3_10
namespace muon {
  enum TunePType{defaultTuneP, improvedTuneP};

  /// High-pt muon
  bool isHighPtMuon(const reco::Muon& muon, const reco::Vertex& vtx, TunePType tunePType){
    
    bool muID =  muon.isGlobalMuon() && muon.globalTrack()->hitPattern().numberOfValidMuonHits() >0 && (muon.numberOfMatchedStations() > 1);
    if(!muID) return false;
    
    if(tunePType == improvedTuneP){ 
      // Get the optimized track
      reco::TrackRef cktTrack = (muon::tevOptimized(muon, 200, 17., 40., 0.25)).first;
      bool momQuality = cktTrack->ptError()/cktTrack->pt() < 0.3;
      bool hits = muon.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5 && muon.innerTrack()->hitPattern().numberOfValidPixelHits() > 0;
      bool ip = fabs(cktTrack->dxy(vtx.position())) < 0.2 && fabs(cktTrack->dz(vtx.position())) < 0.5;
      return muID && hits && momQuality && ip;}
    else if(tunePType == defaultTuneP){
      // Get the optimized track
      reco::TrackRef cktTrack = (muon::tevOptimized(muon, 200, 4., 6., -1)).first;
      bool hits = muon.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 8 && muon.innerTrack()->hitPattern().numberOfValidPixelHits() > 0; 
      bool ip = fabs(cktTrack->dxy(vtx.position())) < 0.2 && fabs(cktTrack->dz(vtx.position())) < 0.5;
      return muID && hits && ip;} 
    
    else return false;
  }
  
  bool isTrackerMuon(const reco::Muon& muon, const reco::Vertex& vtx, TunePType tunePType){
    
    bool muID =  muon.isTrackerMuon() && (muon.numberOfMatchedStations() > 1);
    if(!muID) return false;
    
    if(tunePType == improvedTuneP){ 
      // Get the optimized track
      reco::TrackRef cktTrack = (muon::tevOptimized(muon, 200, 17., 40., 0.25)).first;
      bool momQuality = cktTrack->ptError()/cktTrack->pt() < 0.3;
      bool hits = muon.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5 && muon.innerTrack()->hitPattern().numberOfValidPixelHits() > 0;
      bool ip = fabs(cktTrack->dxy(vtx.position())) < 0.2 && fabs(cktTrack->dz(vtx.position())) < 0.5;
      return muID && hits && momQuality && ip;
    }
    else if(tunePType == defaultTuneP){
      // Get the optimized track
      reco::TrackRef cktTrack = (muon::tevOptimized(muon, 200, 4., 6., -1)).first;
      bool hits = muon.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 8 && muon.innerTrack()->hitPattern().numberOfValidPixelHits() > 0; 
      bool ip = fabs(cktTrack->dxy(vtx.position())) < 0.2 && fabs(cktTrack->dz(vtx.position())) < 0.5;
      return muID && hits && ip;
    }
    else return false;
  }
}// end namespace muon

/// Class declaration
class HPTMNewTunePAnalyzer : public edm::EDAnalyzer {
   public:
      explicit HPTMNewTunePAnalyzer(const edm::ParameterSet&);
      ~HPTMNewTunePAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
 

   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

  TH1F* h_deltaPt;
  TH1F* h_deltaEta;
  TH1F* h_deltaPhi;
  TH1F* h_origPt;
  TH1F* h_modiPt;
  TH1F* h_origPtNoCut;
  TH1F* h_modiPtNoCut;
  TH1F* h_deltaCharge;
  TH1F* h_deltaVx;
  TH1F* h_deltaVy;
  TH1F* h_deltaVz;
};

/// Class definition
HPTMNewTunePAnalyzer::HPTMNewTunePAnalyzer(const edm::ParameterSet& iConfig)

{
  edm::Service<TFileService> fs;
  h_deltaPt = fs->make<TH1F>( "deltaPt"  , "", 200,  -1., 1. );
  h_origPt = fs->make<TH1F>( "origPt"  , "", 400,  0., 4000. );
  h_modiPt = fs->make<TH1F>( "modiPt"  , "", 400,  0., 4000. );
  h_origPtNoCut = fs->make<TH1F>( "origPtNoCut"  , "", 400,  0., 4000. );
  h_modiPtNoCut = fs->make<TH1F>( "modiPtNoCut"  , "", 400,  0., 4000. );
  h_deltaEta = fs->make<TH1F>( "deltaEta"  , "", 100,  -0.02, 0.02 );
  h_deltaPhi = fs->make<TH1F>( "deltaPhi"  , "", 100,  -0.02, 0.02 );
  h_deltaCharge = fs->make<TH1F>( "deltaCharge"  , "", 5,  -2.5, 2.5 );
  h_deltaVx = fs->make<TH1F>( "deltaVx"  , "", 100,  -0.02, 0.02 );
  h_deltaVy = fs->make<TH1F>( "deltaVy"  , "", 100,  -0.02, 0.02 );
  h_deltaVz = fs->make<TH1F>( "deltaVz"  , "", 100,  -0.02, 0.02 );
}


HPTMNewTunePAnalyzer::~HPTMNewTunePAnalyzer()
{
}

// ------------ method called for each event  ------------
void
HPTMNewTunePAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

  //load muons
  edm::Handle<edm::View<reco::Muon> > originalMuonHandle;
  iEvent.getByLabel("muons",originalMuonHandle);
  edm::Handle<edm::View<reco::Muon> > modifiedMuonHandle;
  iEvent.getByLabel("tunePmuons",modifiedMuonHandle);

  //load vertices
  edm::Handle<edm::View<reco::Vertex> > verticesHandle;
  iEvent.getByLabel("offlinePrimaryVertices",verticesHandle);
  const edm::View<reco::Vertex>& vertices = *(verticesHandle.product());
  const reco::Vertex& vertex = vertices[0];

  const edm::View<reco::Muon>& originalMuons = *(originalMuonHandle.product());
  const edm::View<reco::Muon>& modifiedMuons = *(modifiedMuonHandle.product());

  /// First test. Assert that both collections have same size.
  if(originalMuons.size() != modifiedMuons.size())
    return;

  size_t numMuons = originalMuons.size();
  
  for (size_t i=0; i!=numMuons; ++i) {
    const reco::Muon& origMu = originalMuons[i];
    const reco::Muon& modiMu = modifiedMuons[i];
    bool isHPT = muon::isHighPtMuon(modiMu, vertex, muon::improvedTuneP);
    bool isTRK = muon::isTrackerMuon(modiMu, vertex, muon::improvedTuneP);
    
    /// Muon in acceptance
    if(modiMu.pt() > 20 and abs(modiMu.eta() < 2.4) and (isHPT or isTRK)) {

      /// Basic pt plot
      h_origPtNoCut->Fill(origMu.pt());
      h_modiPtNoCut->Fill(modiMu.pt());
      
      /// Compare kinematics
      double deltaPt = (origMu.pt() - modiMu.pt());
      double deltaEta = (origMu.eta() - modiMu.eta());
      double dPhi = deltaPhi(origMu.phi(),modiMu.phi());
      int deltaCharge = (origMu.charge() - modiMu.charge());
      int deltaVx = (origMu.vx() - modiMu.vy());
      int deltaVy = (origMu.vy() - modiMu.vx());
      int deltaVz = (origMu.vz() - modiMu.vz());
      h_deltaPt->Fill(deltaPt/origMu.pt());
      h_deltaEta->Fill(deltaEta);
      h_deltaPhi->Fill(dPhi);
      h_deltaCharge->Fill(deltaCharge);
      h_deltaVx->Fill(deltaVx);
      h_deltaVy->Fill(deltaVy);
      h_deltaVz->Fill(deltaVz);

      /// For those muons where the pT changes TOO MUCH... who are they?
      /// What is happening?
      
      if(fabs(deltaPt/origMu.pt()) > 0.10) {
	h_origPt->Fill(origMu.pt());
	h_modiPt->Fill(modiMu.pt());
      }
      /// Even more!
      if(fabs(deltaPt/origMu.pt()) > 0.60) {
	double originalPt  = origMu.pt();
	double originalEta = origMu.eta();
	double originalPhi = origMu.phi();
	double originalSigmaPtOverPt = origMu.track()->ptError()/origMu.pt();
	double modifiedSigmaPtOverPt = modiMu.track()->ptError()/modiMu.pt();
	bool isGlobalMu = origMu.isGlobalMuon();
	int nPixelHits = origMu.innerTrack()->hitPattern().numberOfValidPixelHits();
	int nTrkHits = origMu.track()->hitPattern().trackerLayersWithMeasurement();
	int nMatchedMuonHits = origMu.numberOfMatchedStations();
	double dxy = origMu.muonBestTrack()->dxy(vertex.position());
	double dz =  origMu.muonBestTrack()->dz(vertex.position());
	printf("********************\n"
	       "Original_Pt \t%g\n"
	       "Original_Eta\t%g\n"
	       "Original_Phi\t%g\n"
	       "Delta_Pt    \t%g\n"
	       "Delta_Eta   \t%g\n"
	       "Delta_Phi   \t%g\n"
	       "Original_(SigmaPt/Pt) \t%g\n"
	       "New_(SigmaPt/Pt)      \t%g\n"
	       "IsGlobalMu  \t%i\n"
	       "NPixelHits  \t%i\n"
	       "NTrkHits    \t%i\n"
	       "NMatchedMuonHits \t%i\n"
	       "dxy         \t%g\n"
	       "dz          \t%g\n",
	       originalPt,originalEta,originalPhi,
	       -deltaPt,-deltaEta,dPhi,
	       originalSigmaPtOverPt,modifiedSigmaPtOverPt,
	       isGlobalMu,nPixelHits,nTrkHits,nMatchedMuonHits,dxy,dz);
      }
	
    }// End pt > 20
    
  
  }// Close loop in muons
}


// ------------ method called once each job just before starting event loop  ------------
void 
HPTMNewTunePAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HPTMNewTunePAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
HPTMNewTunePAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
HPTMNewTunePAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
HPTMNewTunePAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
HPTMNewTunePAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HPTMNewTunePAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HPTMNewTunePAnalyzer);
