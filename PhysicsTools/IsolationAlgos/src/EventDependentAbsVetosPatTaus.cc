#include "PhysicsTools/IsolationAlgos/interface/EventDependentAbsVetosPatTaus.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/Common/interface/AssociationMap.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/Tau.h"

bool //OtherCandidatesDeltaRVeto -> BoostedTauDeltaRVeto
reco::isodeposit::BoostedTauDeltaRVeto::veto(double eta, double phi, float value) const 
{
    for (std::vector<Direction>::const_iterator it = items_.begin(), ed = items_.end(); it != ed; ++it) {
        if (::deltaR2(it->eta(), it->phi(), eta, phi) < deltaR2_) return true;
    }
    return false;
}

void
reco::isodeposit::BoostedTauDeltaRVeto::setEvent(const edm::Event &iEvent, const edm::EventSetup &iSetup) {
    items_.clear();
    edm::Handle<edm::View<reco::Candidate> > candidates;
    iEvent.getByLabel(src_, candidates);
    for (edm::View<reco::Candidate>::const_iterator it = candidates->begin(), ed = candidates->end(); it != ed; ++it) {
        items_.push_back(Direction(it->eta(), it->phi()));
    }
}

bool //OtherCandVeto -> BoostedTauVeto
reco::isodeposit::BoostedTauVeto::veto(double eta, double phi, float value) const 
{
    for (std::vector<Direction>::const_iterator it = items_.begin(), ed = items_.end(); it != ed; ++it) {
        veto_->centerOn(it->eta(), it->phi());
        if ( veto_->veto(eta,phi,value) ) return true;
    }
    return false;
}

void
reco::isodeposit::BoostedTauVeto::setEvent(const edm::Event &iEvent, const edm::EventSetup &iSetup) {
    items_.clear();
    edm::Handle<edm::View<reco::Candidate> > candidates;
    iEvent.getByLabel(src_, candidates);
    for (edm::View<reco::Candidate>::const_iterator it = candidates->begin(), ed = candidates->end(); it != ed; ++it) {
        items_.push_back(Direction(it->eta(), it->phi()));
    }
}

bool //OtherJetConstituentsDeltaRVeto -> BoostedTauConstituentsDeltaRVeto
reco::isodeposit::BoostedTauConstituentsDeltaRVeto::veto(double eta, double phi, float value) const 
{
    for (std::vector<Direction>::const_iterator it = items_.begin(), ed = items_.end(); it != ed; ++it) {
        if (::deltaR2(it->eta(), it->phi(), eta, phi) < dR2constituent_) return true;
    }
    return false;
}

void
reco::isodeposit::BoostedTauConstituentsDeltaRVeto::setEvent(const edm::Event& evt, const edm::EventSetup& es) 
{
    //std::cout << "<BoostedTauConstituentsDeltaRVeto::setEvent>:" << std::endl;
    evt_ = &evt;
}
void
reco::isodeposit::BoostedTauConstituentsDeltaRVeto::initialize()
{
    //std::cout << "<BoostedTauConstituentsDeltaRVeto::initialize>:" << std::endl;
    //std::cout << " vetoDir: eta = " << vetoDir_.eta() << ", phi = " << vetoDir_.phi() << std::endl;
    assert(evt_);
    items_.clear();
    edm::Handle<pat::TauCollection> taus;
    evt_->getByLabel(srcTaus_, taus);
    double dR2min = dR2tau_;
    pat::TauRef matchedTau;
    size_t numTaus = taus->size();
    for ( size_t tauIndex = 0; tauIndex < numTaus; ++tauIndex ) {
      pat::TauRef tau(taus, tauIndex);
      double dR2 = ::deltaR2(vetoDir_.eta(), vetoDir_.phi(), tau->eta(), tau->phi());
      //std::cout << "tau #" << tauIndex << ": Pt = " << tau->pt() << ", eta = " << tau->eta() << ", phi = " << tau->phi() << " (dR = " << sqrt(dR2) << ")" << std::endl;
      if ( dR2 < dR2min ) {
	matchedTau = tau;
	dR2min = dR2;
      }
    }
    if ( matchedTau.isNonnull() ) {
      std::vector<PFCandidatePtr> pfTauConstituents = matchedTau->signalPFCands();
      bool TauPassSelection = false;
      if(matchedTau->pt()>20 && 
	 abs(matchedTau->eta())<2.4 && 
	 matchedTau->tauID("decayModeFindingNewDMs")>0.5 && 
	 matchedTau->tauID("againstMuonLoose")>0.5 && 
	 matchedTau->tauID("againstElectronLoose")>0.5 && 
	 matchedTau->tauID("byVLooseIsolationMVA3newDMwoLT")>0.5) TauPassSelection = true;

      if(TauPassSelection){
	int idx = 0;
	for(std::vector<reco::PFCandidatePtr>::const_iterator pfTauConstituent=pfTauConstituents.begin(); pfTauConstituent!=pfTauConstituents.end(); ++pfTauConstituent){
	  items_.push_back(Direction((*pfTauConstituent)->eta(), (*pfTauConstituent)->phi()));
	  //std::cout<<"pfCand #"<<idx<<": Pt = "<<(*pfTauConstituent)->pt()<<", eta = "<<(*pfTauConstituent)->eta()<<", phi = "<<(*pfTauConstituent)->phi()<<std::endl;
	  ++idx;
	}
      }
    }
}

void 
reco::isodeposit::BoostedTauConstituentsDeltaRVeto::centerOn(double eta, double phi) 
{ 
    //std::cout << "<BoostedTauConstituentsDeltaRVeto::centerOn>:" << std::endl;
    //std::cout << " eta = " << eta << std::endl;
    //std::cout << " phi = " << phi << std::endl;
    vetoDir_ = Direction(eta,phi); 
    initialize();
}
