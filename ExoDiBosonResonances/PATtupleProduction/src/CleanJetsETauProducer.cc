// -*- C++ -*-
//
// Package:    CleanJetsETauProducer
// Class:      CleanJetsETauProducer
// 
/**\class CleanJetsETauProducer CleanJetsETauProducer.cc ExoDiBosonResonances/EDBRTauAnalyzer/src/CleanJetsETauProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Aniello Spiezia,21 1-007,+41227676459,
//         Created:  Mon Sep  9 13:14:05 CEST 2013
// $Id$
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
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h" 

//new inclusion
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "FWCore/Framework/interface/EDProducer.h" 
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "EGamma/EGammaAnalysisTools/interface/EGammaMvaEleEstimator.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "TrackingTools/IPTools/interface/IPTools.h"


using namespace edm;
using namespace std;
using namespace reco;


//
// class declaration
//

class CleanJetsETauProducer : public edm::EDProducer {
public:
  explicit CleanJetsETauProducer(const edm::ParameterSet&);
  ~CleanJetsETauProducer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

  EGammaMvaEleEstimator* myMVANonTrig;
 

  //HISTOGRAMS
  //TH1D* metPt;

  // ----------member data ---------------------------
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
CleanJetsETauProducer::CleanJetsETauProducer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  //produces<PFJetCollection>( "ak5PFJetsNoMu" );
  produces<PFJetCollection>( );

  // NOTE: It is safer and crab-compliant to get the files locally, i.e in EGamma/EGammaAnalysisTools/data
  // (see the downloard.url file in that directory)
  // Alternatively (for tests), they can be read from AFS:
  std::vector<std::string> myManualWeigths;
  myManualWeigths.push_back("Electrons_BDTG_NonTrigV0_Cat1.weights.xml");
  myManualWeigths.push_back("Electrons_BDTG_NonTrigV0_Cat2.weights.xml");
  myManualWeigths.push_back("Electrons_BDTG_NonTrigV0_Cat3.weights.xml");
  myManualWeigths.push_back("Electrons_BDTG_NonTrigV0_Cat4.weights.xml");
  myManualWeigths.push_back("Electrons_BDTG_NonTrigV0_Cat5.weights.xml");
  myManualWeigths.push_back("Electrons_BDTG_NonTrigV0_Cat6.weights.xml");

  Bool_t manualCat = true;
  myMVANonTrig = new EGammaMvaEleEstimator();
  myMVANonTrig->initialize("BDT",
            EGammaMvaEleEstimator::kNonTrig,
            manualCat, 
            myManualWeigths);

}


CleanJetsETauProducer::~CleanJetsETauProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
CleanJetsETauProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  Handle<reco::PFJetCollection> PFJets;
  iEvent.getByLabel("ak5PFJets", PFJets);
  auto_ptr<PFJetCollection> SetOfJets( new PFJetCollection );

  Handle<reco::GsfElectronCollection> electrons;
  iEvent.getByLabel("gsfElectrons", electrons);

  Handle<reco::VertexCollection> vertices;
  //iEvent.getByLabel("offlinePrimaryVertices", vertices);
  iEvent.getByLabel("primaryVertexFilter", vertices);
  reco::Vertex primaryVertex;
  primaryVertex = vertices->at(0);

  InputTag  reducedEBRecHitCollection(string("reducedEcalRecHitsEB"));
  InputTag  reducedEERecHitCollection(string("reducedEcalRecHitsEE"));
  EcalClusterLazyTools lazyTools(iEvent, iSetup, reducedEBRecHitCollection, reducedEERecHitCollection);
  edm::ESHandle<TransientTrackBuilder> trackBuilder_;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", trackBuilder_);  
  const TransientTrackBuilder thebuilder = *(trackBuilder_.product());
  bool debugMVAclass = false;
  
  for(reco::PFJetCollection::const_iterator iJet = PFJets->begin(); iJet != PFJets->end(); ++iJet){
    
    math::XYZTLorentzVector pfmomentum;  
    std::vector<edm::Ptr<Candidate> > jetConstituents;
    jetConstituents.clear();
    std::vector<reco::PFCandidatePtr> JetPFCands = iJet->getPFConstituents();
    reco::PFJet::Specific specs = iJet->getSpecific();
    
    for(std::vector<edm::Ptr<reco::PFCandidate> >::iterator i = JetPFCands.begin(); i != JetPFCands.end(); ++i){
      reco::PFCandidate pfCand = *i;
      
      if (pfCand.particleId() == 2){ //the PFCandidate is a electron
	reco::GsfElectronRef theRecoElectron = pfCand.gsfElectronRef();
	double myMVANonTrigValue = myMVANonTrig->mvaValue(*theRecoElectron,primaryVertex,thebuilder,lazyTools,debugMVAclass);
	bool selected = false;
	float thiseta = fabs(theRecoElectron->superCluster()->eta());

	//CUTS TAKEN FROM https://twiki.cern.ch/twiki/bin/view/CMS/MultivariateElectronIdentification
	if(theRecoElectron->pt()>7 && theRecoElectron->pt()<10 && thiseta<2.5
	   //&& theRecoElectron->gsfTrack()->trackerExpectedHitsInner().numberOfHits()<=1
	   //&& ip3dSig<4
	   //&& iso<0.4
	   ){
	  if(thiseta<0.8                    && myMVANonTrigValue>0.47 ) selected=true;
	  if(thiseta>0.8   && thiseta<1.479 && myMVANonTrigValue>0.004) selected=true;
	  if(thiseta>1.479 && thiseta<2.5   && myMVANonTrigValue>0.295) selected=true;
	}
	if(theRecoElectron->pt()>10 && thiseta<2.5
	   //&& theRecoElectron->gsfTrack()->trackerExpectedHitsInner().numberOfHits()<=1
	   //&& ip3dSig<4
	   //&& iso<0.4
	   ){
	  if(thiseta<0.8                    && myMVANonTrigValue>-0.34) selected=true;
	  if(thiseta>0.8   && thiseta<1.479 && myMVANonTrigValue>-0.65) selected=true;
	  if(thiseta>1.479 && thiseta<2.5   && myMVANonTrigValue>0.60 ) selected=true;
	}
	/*
	if(theRecoElectron->pt()>8){
	  if(fabs(theRecoElectron->superCluster()->eta())<=1.479 
	     && theRecoElectron->deltaEtaSuperClusterTrackAtVtx()<0.007 
	     && theRecoElectron->deltaPhiSuperClusterTrackAtVtx()<0.8 
	     && theRecoElectron->sigmaIetaIeta()<0.01 && theRecoElectron->hadronicOverEm()<0.15 
	     && fabs(theRecoElectron->gsfTrack()->dxy(primaryVertex.position()))<0.04 
	     && fabs(theRecoElectron->gsfTrack()->dz(primaryVertex.position()))<0.2){
	    selected=true;
	  }
	  if(fabs(theRecoElectron->superCluster()->eta())>1.479 
	     && fabs(theRecoElectron->superCluster()->eta())<2.5 
	     && theRecoElectron->deltaEtaSuperClusterTrackAtVtx()<0.01 
	     && theRecoElectron->deltaPhiSuperClusterTrackAtVtx()<0.7 
	     && theRecoElectron->sigmaIetaIeta()<0.03
	     && fabs(theRecoElectron->gsfTrack()->dxy(primaryVertex.position()))<0.04 
	     && fabs(theRecoElectron->gsfTrack()->dz(primaryVertex.position()))<0.2){
	    selected=true;
	  }
	}
	*/

	if(selected){
	  //cout<<iJet->pt()<<" "<<theRecoElectron->pt()<<endl;
	  specs.mElectronEnergy -= pfCand.p4().e();
	  specs.mElectronMultiplicity -= 1;
	  specs.mChargedEmEnergy -= pfCand.p4().e();
	  specs.mChargedMultiplicity -= 1;
	}
	else{
	  pfmomentum += pfCand.p4();
	  jetConstituents.push_back((*i));
	}
	
      }//the PFCandidate is a electron  
      else {
	pfmomentum += pfCand.p4();
	jetConstituents.push_back((*i));
      }
    }//loop over PFCandidate

    PFJet electronfreePFJet(pfmomentum, specs, jetConstituents);
    SetOfJets->push_back( electronfreePFJet );
  
  }//loop over PFJet

  iEvent.put( SetOfJets );

#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
CleanJetsETauProducer::beginJob()
{
  Service<TFileService> fs;
  //metPt             = fs->make<TH1D>("metPt",            "metPt",            5000, 0, 5000);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
CleanJetsETauProducer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
CleanJetsETauProducer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
CleanJetsETauProducer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
CleanJetsETauProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
CleanJetsETauProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
CleanJetsETauProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


//define this as a plug-in
DEFINE_FWK_MODULE(CleanJetsETauProducer);
