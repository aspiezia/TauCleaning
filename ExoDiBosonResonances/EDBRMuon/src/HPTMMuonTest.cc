// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "ExoDiBosonResonances/EDBRMuon/interface/HPTMMuonSelector.h"
//
// class declaration
//

class HPTMMuonTest : public edm::EDFilter {
   public:
      explicit HPTMMuonTest(const edm::ParameterSet&);
      ~HPTMMuonTest();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual bool beginRun(edm::Run&, edm::EventSetup const&);
      virtual bool endRun(edm::Run&, edm::EventSetup const&);
      virtual bool beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual bool endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------

  edm::InputTag muons_;
  hptm::MuonSelector muonSelector_;
  size_t muonsPassedHPT_;
  size_t muonsPassedTRK_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
HPTMMuonTest::HPTMMuonTest(const edm::ParameterSet& iConfig)
{
  muons_ = iConfig.getParameter<edm::InputTag>("muons");
  muonsPassedHPT_ = 0;
  muonsPassedTRK_ = 0;
}


HPTMMuonTest::~HPTMMuonTest()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
HPTMMuonTest::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<edm::View<pat::Muon> > muonHandle;
  iEvent.getByLabel(muons_,muonHandle);
  const edm::View<pat::Muon>& muons = *(muonHandle.product());

  edm::Handle<edm::View<reco::Vertex> > vertexHandle;
  iEvent.getByLabel("offlinePrimaryVertices",vertexHandle);
  const edm::View<reco::Vertex>& vertices = *(vertexHandle.product());

  const reco::Vertex& vertex = vertices[0];

  size_t numMuonsPassed = 0;

  for(size_t muonNr = 0; muonNr != muons.size(); ++muonNr) {
    
    const pat::Muon& recoMu = muons[muonNr];
    
    bool passedHPT = muonSelector_.checkMuonID(recoMu, vertex.position(), hptm::HIGHPT_OLD);
    bool passedTRK = muonSelector_.checkMuonID(recoMu, vertex.position(), hptm::TRACKER);
    if(passedHPT)  muonsPassedHPT_++;
    if(passedTRK)  muonsPassedTRK_++;

    if (passedHPT or passedTRK) numMuonsPassed++;
    
  }
   
  return (numMuonsPassed > 1);
}

// ------------ method called once each job just before starting event loop  ------------
void 
HPTMMuonTest::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HPTMMuonTest::endJob() {
  std::cout << "We got " << int(muonsPassedHPT_) << " high-pt muons and "
	    << int(muonsPassedTRK_) << " tracker muons." << std::endl; 
}

// ------------ method called when starting to processes a run  ------------
bool 
HPTMMuonTest::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
HPTMMuonTest::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
HPTMMuonTest::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
HPTMMuonTest::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HPTMMuonTest::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(HPTMMuonTest);
