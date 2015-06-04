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
#include "DataFormats/PatCandidates/interface/Muon.h"
//#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/MuonReco/interface/MuonCocktails.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "TH1F.h"
//
// class declaration
//

/// Class declaration
class HPTMNewTunePPATAnalyzer : public edm::EDAnalyzer {
   public:
      explicit HPTMNewTunePPATAnalyzer(const edm::ParameterSet&);
      ~HPTMNewTunePPATAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
 

   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
};

/// Class definition
HPTMNewTunePPATAnalyzer::HPTMNewTunePPATAnalyzer(const edm::ParameterSet& iConfig)

{
}


HPTMNewTunePPATAnalyzer::~HPTMNewTunePPATAnalyzer()
{
}

// ------------ method called for each event  ------------
void
HPTMNewTunePPATAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

  //load muons
  edm::Handle<reco::MuonCollection> originalMuonHandle;
  iEvent.getByLabel("tunePmuons",originalMuonHandle);
  const reco::MuonCollection& muons = *(originalMuonHandle.product());
  printf("Muons loaded.\n");

  //load valueMap
  edm::Handle<edm::ValueMap<float> > vmHandle;
  iEvent.getByLabel("muonTrackError",vmHandle);
  const edm::ValueMap<float>& vm = *(vmHandle.product());

  
  for (size_t i = 0; i < muons.size(); i++) {
    edm::Ref<reco::MuonCollection> muonRef(originalMuonHandle,i);
    double trackError = vm[muonRef];
    printf("From the value map: muon pt = %g, trackError = %g\n",
	   muonRef->pt(), trackError);
    
  }

  /// With the pat muons
  edm::Handle<edm::View<pat::Muon> > patMuonHandle;
  iEvent.getByLabel("patMuons",patMuonHandle);
  const edm::View<pat::Muon>& patMuons = *(patMuonHandle.product());

  /// First test. Assert that both collections have same size.
 if(muons.size() != patMuons.size())
    return;

 size_t numMuons = muons.size();
 
 for (size_t i=0; i!=numMuons; ++i) {
    const pat::Muon& patMu = patMuons[i];
    double PATtrackError = patMu.userFloat("muonTrackError");
    double PATpt = patMu.pt();
    printf("From the PAT value: muon pt = %g, trackError = %g\n",
	   PATpt, PATtrackError);
 }
 
}


// ------------ method called once each job just before starting event loop  ------------
void 
HPTMNewTunePPATAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HPTMNewTunePPATAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
HPTMNewTunePPATAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
HPTMNewTunePPATAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
HPTMNewTunePPATAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
HPTMNewTunePPATAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HPTMNewTunePPATAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HPTMNewTunePPATAnalyzer);
