#ifndef PhysicsTools_IsolationAlgos_EventDependentAbsVetosPatTaus_h
#define PhysicsTools_IsolationAlgos_EventDependentAbsVetosPatTaus_h

#include "PhysicsTools/IsolationAlgos/interface/EventDependentAbsVeto.h"
#include "FWCore/Utilities/interface/InputTag.h"

//OtherCandidatesDeltaRVeto      -> BoostedTauDeltaRVeto
//OtherCandVeto                  -> BoostedTauVeto
//OtherJetConstituentsDeltaRVeto -> BoostedTauConstituentsDeltaRVeto
namespace reco {
 namespace isodeposit {
    class BoostedTauDeltaRVeto : public EventDependentAbsVeto {
      public:
          //! Create a veto specifying the input collection of the candidates, and the deltaR
          BoostedTauDeltaRVeto(const edm::InputTag candidates, double deltaR) :
            src_(candidates), deltaR2_(deltaR*deltaR) { }
   
          // Virtual destructor (should always be there) 
          virtual ~BoostedTauDeltaRVeto() {} 

          //! Return "true" if a deposit at specific (eta,phi) with that value must be vetoed in the sum
          //! This is true if the deposit is within the configured deltaR from any item of the source collection
          virtual bool veto(double eta, double phi, float value) const ;

          //! Nothing to do for this
          virtual void centerOn(double eta, double phi) { }

          //! Picks up the directions of the given candidates
          virtual void setEvent(const edm::Event &iEvent, const edm::EventSetup &iSetup) ;

      private:
          edm::InputTag src_;
          float         deltaR2_;
          std::vector<Direction> items_;
    };

    class BoostedTauVeto : public EventDependentAbsVeto {
      public:
          //! Create a veto specifying the input collection of the candidates, and the deltaR
          BoostedTauVeto(const edm::InputTag candidates, AbsVeto *veto) :
            src_(candidates), veto_(veto) { }
   
          // Virtual destructor (should always be there) 
          virtual ~BoostedTauVeto() {} 

          //! Return "true" if a deposit at specific (eta,phi) with that value must be vetoed in the sum
          //! This is true if the deposit is within the stored AbsVeto of any item of the source collection
          virtual bool veto(double eta, double phi, float value) const ;

          //! Nothing to do for this
          virtual void centerOn(double eta, double phi) { }

          //! Picks up the directions of the given candidates
          virtual void setEvent(const edm::Event &iEvent, const edm::EventSetup &iSetup) ;

      private:
          edm::InputTag src_;
          std::vector<Direction> items_;
          std::auto_ptr<AbsVeto> veto_;
    };

    class BoostedTauConstituentsDeltaRVeto : public EventDependentAbsVeto {
      public:
          //! Create a veto specifying the input collection of the taus, the candidates, and the deltaR
          BoostedTauConstituentsDeltaRVeto(Direction dir, const edm::InputTag& taus, double dRtau, const edm::InputTag& pfCandAssocMap, double dRconstituent) 
	    : evt_(0),
	      vetoDir_(dir), 
	      srcTaus_(taus), 
	      dR2tau_(dRtau*dRtau), 
	      srcPFCandAssocMap_(pfCandAssocMap), 
	      dR2constituent_(dRconstituent*dRconstituent) 
	  { 
	    //std::cout << "<BoostedTauConstituentsDeltaRVeto::BoostedTauConstituentsDeltaRVeto>:" << std::endl;
	    //std::cout << " vetoDir: eta = " << vetoDir_.eta() << ", phi = " << vetoDir_.phi() << std::endl;
	    //std::cout << " srcTaus = " << srcTaus_.label() << ":" << srcTaus_.instance() << std::endl;
	    //std::cout << " dRtau = " << sqrt(dR2tau_) << std::endl;
	    //std::cout << " srcPFCandAssocMap = " << srcPFCandAssocMap_.label() << ":" << srcPFCandAssocMap_.instance() << std::endl;
	    //std::cout << " dRconstituent = " << sqrt(dR2constituent_) << std::endl;
	  }
   
          // Virtual destructor (should always be there) 
          virtual ~BoostedTauConstituentsDeltaRVeto() {} 

          //! Return "true" if a deposit at specific (eta,phi) with that value must be vetoed in the sum
          //! This is true if the deposit is within the stored AbsVeto of any item of the source collection
          virtual bool veto(double eta, double phi, float value) const;

          //! Set axis for matching taus
          virtual void centerOn(double eta, double phi);

          //! Picks up the directions of the given candidates
          virtual void setEvent(const edm::Event& evt, const edm::EventSetup& es);

      private:
	  void initialize();

	  const edm::Event* evt_;

	  Direction vetoDir_; 
          edm::InputTag srcTaus_;
	  double dR2tau_;
	  edm::InputTag srcPFCandAssocMap_;
	  double dR2constituent_;
          std::vector<Direction> items_;
    };     
 }
}
#endif
