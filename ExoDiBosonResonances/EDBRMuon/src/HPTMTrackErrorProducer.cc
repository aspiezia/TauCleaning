// system include files
#include <memory>
#include <iostream>
#include <string>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonCocktails.h"

/// Class declaration
class HPTMTrackErrorProducer : public edm::EDProducer {
public:
  explicit HPTMTrackErrorProducer(const edm::ParameterSet&);
  ~HPTMTrackErrorProducer();

private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data --------------------------
  edm::InputTag muLabel_;
};

/// Class definition
HPTMTrackErrorProducer::HPTMTrackErrorProducer(const edm::ParameterSet& iConfig)
{
  muLabel_=iConfig.getParameter<edm::InputTag>("muLabel");
  produces<edm::ValueMap<float> >().setBranchAlias("TrackErrors");
}


HPTMTrackErrorProducer::~HPTMTrackErrorProducer()
{

}

// ------------ method called to produce the data  ------------
void
HPTMTrackErrorProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //load muons
  edm::Handle<reco::MuonCollection> muHandle;
  iEvent.getByLabel(muLabel_,muHandle);
  const reco::MuonCollection& muons = *(muHandle.product());

  //vector for trackError values
  std::vector<float> values;
  values.reserve(muons.size());
  values.clear();
  
  //loop over muons, get the optimal track again
  for(size_t muonNr=0;muonNr<muons.size();muonNr++){
    reco::Muon * recoMu = muons[muonNr].clone();
    
    if(recoMu->isTrackerMuon()) {
    // call to get the optimal muon track
      reco::TrackRef cktTrack = (muon::tevOptimized(*recoMu, 200, 17., 40., 0.25)).first; 
      double pt = cktTrack->pt();
      double sigmaPt = cktTrack->ptError();
      double trackError = sigmaPt/pt;
      values.push_back(float(trackError));
    }
    else {
      // If it is not a tracker muon, set this value to 9999
      // (this indicates that there is no suitable TeV refit).
      values.push_back(9999);
    }
    
    delete recoMu;

  }

  //printf("++++++++++++++++\n");
  //printf("Added the following track errors:\n");
  //for (size_t ii = 0; ii!= values.size(); ++ii) {
  // printf("Track error = %g\n",values.at(ii));
  //}
  //printf("++++++++++++++++\n");

  std::auto_ptr<edm::ValueMap<float> > out(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler(*out);
  filler.insert(muHandle, values.begin(), values.end());
  filler.fill();

  // put value map into event
  iEvent.put(out);
}

// ------------ method called once each job just before starting event loop  ------------
void
HPTMTrackErrorProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
HPTMTrackErrorProducer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(HPTMTrackErrorProducer);
