#include "ExoDiBosonResonances/EDBRMuon/src/PATMuonSmearProducer.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "CLHEP/Random/RandGaussQ.h"

PATMuonSmearProducer::PATMuonSmearProducer(const edm::ParameterSet& iConfig):
  inputMuons_(iConfig.getParameter<edm::InputTag>("src")),
  momentumScale_(iConfig.getParameter<double>("momentumScale")),
  highPtRegion_(iConfig.getParameter<double>("highPtRegion")),
  extraMomentumScale_(iConfig.getParameter<double>("extraMomentumScale")),
  muonMass(0.1056583715),//GeV
  factor(1.0),
  isScale_(iConfig.getParameter<bool>("isScale")),
  isPositive_(iConfig.getParameter<bool>("isPositive"))
{
  produces<std::vector<pat::Muon> >();
}

PATMuonSmearProducer::~PATMuonSmearProducer(){
}

void PATMuonSmearProducer::endJob(){}

void PATMuonSmearProducer::beginJob(){
  edm::Service<edm::RandomNumberGenerator> rng;
  if ( ! rng.isAvailable()) {
    throw cms::Exception("Configuration")
      << "PATMuonSmearProducer requires the RandomNumberGeneratorService\n"
      "which is not present in the configuration file.  You must add the service\n"
      "in the configuration file or remove the modules that require it.";
  }
  
  // Get the engine.
  // Note this can be done in any method of the module
  CLHEP::HepRandomEngine& engine = rng->getEngine();

  // engine MUST be a reference here, if a pointer is used the
  // distribution will destroy the engine in its destructor, a major
  // problem because the service owns the engine and will destroy it 
  if(momentumScale_ < 0)
    throw cms::Exception("Configuration")
      << "Use isPositive = cms.bool(False) to make a negative scale.";
  
  // Setting up the gauss distribution
  if(momentumScale_ > 1E-6)
    gaussDistribution_.reset(new CLHEP::RandGaussQ(engine,1.0,momentumScale_));
}

void PATMuonSmearProducer::produce(edm::Event& iEvent, const edm::EventSetup&) {

  //load muons
  edm::Handle<std::vector<pat::Muon> > muHandle;
  iEvent.getByLabel(inputMuons_,muHandle);
  const std::vector<pat::Muon>& muons = *(muHandle.product());

  //smeared muons
  std::auto_ptr<std::vector<pat::Muon> > smearedMuons (new std::vector<pat::Muon>);

  for(size_t ii=0; ii!=muons.size(); ++ii) {

    pat::Muon* clonedMu = (pat::Muon*)muons[ii].clone();

    /// If isScale_, we're shifting the average momentum
    /// else, we're smearing it gaussianly.
    if(isScale_) {
      /// If it's a scale, maybe we should consider the high pT region.
      if(clonedMu->pt() > highPtRegion_) {
	double extraError = ((clonedMu->pt())/1000.0)*extraMomentumScale_;
	double totalError = sqrt(extraError*extraError + momentumScale_*momentumScale_);
	if(isPositive_)
	  factor = 1.0+totalError;
	else
	  factor = 1.0-totalError;
      }      
      else {
	if(isPositive_)
	  factor = 1.0+momentumScale_;
	else
	  factor = 1.0-momentumScale_;
      } // Closes highPtRegion_
    } // Closes isScale_
    else
      factor = gaussDistribution_->fire();
    
    /// First, scale the momentum 
    double pX = factor*(clonedMu->px());
    double pY = factor*(clonedMu->py());
    double pZ = clonedMu->pz(); /// Yes, pZ is not smeared.
    double energy = sqrt((pX*pX)+(pY*pY)+(pZ*pZ)+(muonMass*muonMass));
    
    math::XYZTLorentzVector newP4(pX,pY,pZ,energy);
      
    clonedMu->setP4(newP4); /// This sets the new LorentzVector

    smearedMuons->push_back(*(clonedMu));
  }
  
  iEvent.put(smearedMuons);
}
//define this as a plug-in
DEFINE_FWK_MODULE(PATMuonSmearProducer);
