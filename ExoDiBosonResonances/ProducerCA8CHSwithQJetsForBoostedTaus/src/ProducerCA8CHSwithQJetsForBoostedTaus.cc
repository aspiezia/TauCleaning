// -*- C++ -*-
//
// Package:    ProducerCA8CHSwithQJetsForBoostedTaus
// Class:      ProducerCA8CHSwithQJetsForBoostedTaus
// 
/**\class ProducerCA8CHSwithQJetsForBoostedTaus ProducerCA8CHSwithQJetsForBoostedTaus.cc ExoDiBosonResonances/ProducerCA8CHSwithQJetsForBoostedTaus/src/ProducerCA8CHSwithQJetsForBoostedTaus.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Cesar Bernardes,40 2-B15,+41227671625,
//         Created:  Mon Dec  9 20:34:01 CET 2013
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/JetReco/interface/PFJet.h"

#include <vector>
#include <string>
#include <utility>


//
// class declaration
//

class ProducerCA8CHSwithQJetsForBoostedTaus : public edm::EDProducer {
   public:
      explicit ProducerCA8CHSwithQJetsForBoostedTaus(const edm::ParameterSet&);
      ~ProducerCA8CHSwithQJetsForBoostedTaus();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      static bool sortJetConstituents(const reco::Candidate* cand1, const reco::Candidate* cand2);

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void endRun(edm::Run&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------

      edm::InputTag jets_;
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
ProducerCA8CHSwithQJetsForBoostedTaus::ProducerCA8CHSwithQJetsForBoostedTaus(const edm::ParameterSet& iConfig):
   jets_( iConfig.getParameter<edm::InputTag>( "jets" ) )
{
   //register your products
/* Examples
   produces<ExampleData2>();

   //if do put with a label
   produces<ExampleData2>("label");
 
   //if you want to put into the Run
   produces<ExampleData2,InRun>();
*/
   //now do what ever other initialization is needed
   produces<std::vector<pat::Jet> >();       
}


ProducerCA8CHSwithQJetsForBoostedTaus::~ProducerCA8CHSwithQJetsForBoostedTaus()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
ProducerCA8CHSwithQJetsForBoostedTaus::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   //std::cout << "==================---------NEW EVENT--------+++++++++++++++++++" << std::endl;
   
   // Handle to the jet collection  
   edm::Handle< std::vector<pat::Jet> > Jets;
   iEvent.getByLabel(jets_, Jets);

   //prepare output collection
   std::auto_ptr< std::vector<pat::Jet> > jetColl( new std::vector<pat::Jet> );
   jetColl->reserve(Jets->size()); 

  
   //to store pointers for the charged constituents of the Jet
   std::vector<const reco::Candidate*> VecPtrJetChargedConstituents;
   //to store pointers for the neutral constituents of the Jet
   std::vector<const reco::Candidate*> VecPtrJetNeutralConstituents;
   for(unsigned int i = 0; i < Jets->size();++i){
      //std::cout << "==================---------NEW JET IN THE EVENT--------+++++++++++++++++++" << std::endl;  
      const pat::Jet & jet = (*Jets)[i];

      //std::cout << "Number of daughters of the Jet:  " << jet.numberOfDaughters() << std::endl;

      jetColl->push_back(jet);  

      //Loop in the daughters of the jets to fill the vectors above
      for(unsigned int p=0; p!= jet.numberOfDaughters(); ++p){
         const reco::Candidate* dau = jet.daughter(p);
         const reco::PFJet* JetConstit = (reco::PFJet*) dau;
         //std::cout << "pT of the Constituent:  " << JetConstit->pt() << std::endl;
         if(fabs(JetConstit->charge()) > 0){
               //std::cout << "pT of the Charged Constituents of the Jet : " << JetConstit->pt() << std::endl;
               VecPtrJetChargedConstituents.push_back(JetConstit);
         }else{
               //std::cout << "pT of the Neutral Constituents of the Jet : " << JetConstit->pt() << std::endl;
               VecPtrJetNeutralConstituents.push_back(JetConstit);
         }        
         
      }  
      //do an ordering in pT in the collection of the pointers of charged/neutral constituents
      std::sort(VecPtrJetChargedConstituents.begin(),VecPtrJetChargedConstituents.end(),sortJetConstituents);
      std::sort(VecPtrJetNeutralConstituents.begin(),VecPtrJetNeutralConstituents.end(),sortJetConstituents);
      //Sum of pT for the Charged/Neutral jets constituents
      double SumPT_chConstit = 0;
      double SumPT_neutrConstit = 0;
    
      //Number of charged/neutral constituents of the jets (constituents with pT>0,1,2,3,4,5,10,15,20,25,30,35,40,45,50 GeV/c)
      int N_chConstit_pT0 = 0;
      int N_chConstit_pT1 = 0;
      int N_chConstit_pT2 = 0;
      int N_chConstit_pT3 = 0;
      int N_chConstit_pT4 = 0;
      int N_chConstit_pT5 = 0;
      int N_chConstit_pT10 = 0;
      int N_chConstit_pT15 = 0;
      int N_chConstit_pT20 = 0;
      int N_chConstit_pT25 = 0;
      int N_chConstit_pT30 = 0;
      int N_chConstit_pT35 = 0;
      int N_chConstit_pT40 = 0;
      int N_chConstit_pT45 = 0;
      int N_chConstit_pT50 = 0;

      int N_neutrConstit_pT0 = 0;
      int N_neutrConstit_pT1 = 0;
      int N_neutrConstit_pT2 = 0;
      int N_neutrConstit_pT3 = 0;
      int N_neutrConstit_pT4 = 0;
      int N_neutrConstit_pT5 = 0;
      int N_neutrConstit_pT10 = 0;
      int N_neutrConstit_pT15 = 0;
      int N_neutrConstit_pT20 = 0;
      int N_neutrConstit_pT25 = 0;
      int N_neutrConstit_pT30 = 0;
      int N_neutrConstit_pT35 = 0;
      int N_neutrConstit_pT40 = 0;
      int N_neutrConstit_pT45 = 0;
      int N_neutrConstit_pT50 = 0;
    
      //Loops in the vectors ordered above
      for(std::vector<const reco::Candidate*>::const_iterator JetChargedConstit_ptr = VecPtrJetChargedConstituents.begin();                                                                                          JetChargedConstit_ptr != VecPtrJetChargedConstituents.end(); ++JetChargedConstit_ptr){
         const reco::Candidate* ptr_ChConstit = (*JetChargedConstit_ptr);
         //std::cout << "PT of the Charged Constituents of the Jet (ordering in pT) : " << ptr_ChConstit->pt() << std::endl; 
         SumPT_chConstit += ptr_ChConstit->pt();
          
         N_chConstit_pT0++;
         if(ptr_ChConstit->pt() > 1){
            N_chConstit_pT1++;
            if(ptr_ChConstit->pt() > 2){
               N_chConstit_pT2++;
               if(ptr_ChConstit->pt() > 3){
                  N_chConstit_pT3++; 
                  if(ptr_ChConstit->pt() > 4){
                     N_chConstit_pT4++; 
                     if(ptr_ChConstit->pt() > 5){
                        N_chConstit_pT5++;
                        if(ptr_ChConstit->pt() > 10){
                           N_chConstit_pT10++;
                           if(ptr_ChConstit->pt() > 15){
                              N_chConstit_pT15++;
                              if(ptr_ChConstit->pt() > 20){
                                 N_chConstit_pT20++; 
                                 if(ptr_ChConstit->pt() > 25){
                                    N_chConstit_pT25++; 
                                    if(ptr_ChConstit->pt() > 30){
                                       N_chConstit_pT30++;
                                       if(ptr_ChConstit->pt() > 35){
                                          N_chConstit_pT35++;
                                          if(ptr_ChConstit->pt() > 40){
                                             N_chConstit_pT40++;
                                             if(ptr_ChConstit->pt() > 45){
                                                N_chConstit_pT45++;
                                                if(ptr_ChConstit->pt() > 50){
                                                   N_chConstit_pT50++;
                                                }
                                             }
                                          } 
                                       }
                                    }
                                 } 
                              }
                           } 
                        } 
                     }   
                  }
               }
            }
         }
      }
      for(std::vector<const reco::Candidate*>::const_iterator JetNeutralConstit_ptr = VecPtrJetNeutralConstituents.begin();                                                                                          JetNeutralConstit_ptr != VecPtrJetNeutralConstituents.end(); ++JetNeutralConstit_ptr){
         const reco::Candidate* ptr_NeutrConstit = (*JetNeutralConstit_ptr);
         //std::cout << "PT of the Neutral Constituents of the Jet (ordering in pT) : " << ptr_NeutrConstit->pt() << std::endl; 
         SumPT_neutrConstit += ptr_NeutrConstit->pt();

         N_neutrConstit_pT0++;
         if(ptr_NeutrConstit->pt() > 1){
            N_neutrConstit_pT1++;
            if(ptr_NeutrConstit->pt() > 2){
               N_neutrConstit_pT2++;
               if(ptr_NeutrConstit->pt() > 3){
                  N_neutrConstit_pT3++;
                  if(ptr_NeutrConstit->pt() > 4){
                     N_neutrConstit_pT4++;
                     if(ptr_NeutrConstit->pt() > 5){
                        N_neutrConstit_pT5++;
                        if(ptr_NeutrConstit->pt() > 10){
                           N_neutrConstit_pT10++;
                           if(ptr_NeutrConstit->pt() > 15){
                              N_neutrConstit_pT15++;
                              if(ptr_NeutrConstit->pt() > 20){
                                 N_neutrConstit_pT20++;
                                 if(ptr_NeutrConstit->pt() > 25){
                                    N_neutrConstit_pT25++;
                                    if(ptr_NeutrConstit->pt() > 30){
                                       N_neutrConstit_pT30++;
                                       if(ptr_NeutrConstit->pt() > 35){
                                          N_neutrConstit_pT35++;
                                          if(ptr_NeutrConstit->pt() > 40){
                                             N_neutrConstit_pT40++;
                                             if(ptr_NeutrConstit->pt() > 45){
                                                N_neutrConstit_pT45++;
                                                if(ptr_NeutrConstit->pt() > 50){
                                                   N_neutrConstit_pT50++;
                                                }
                                             }
                                          }
                                       }
                                    }
                                 }
                              }
                           }
                        }
                     }
                  }
               }
            }
         } 
      } 
      //std::cout << "Sum of pT for the charged constituents of the Jets: " << SumPT_chConstit << std::endl;
      jetColl->back().addUserFloat("SumPTchConstit", (float)SumPT_chConstit); 
      //std::cout << "Sum of pT for the neutral constituents of the Jets: " << SumPT_neutrConstit << std::endl;
      jetColl->back().addUserFloat("SumPTneutrConstit", (float)SumPT_neutrConstit);

      jetColl->back().addUserInt("NchConstitPT0", N_chConstit_pT0);
      jetColl->back().addUserInt("NchConstitPT1", N_chConstit_pT1);
      jetColl->back().addUserInt("NchConstitPT2", N_chConstit_pT2);
      jetColl->back().addUserInt("NchConstitPT3", N_chConstit_pT3);
      jetColl->back().addUserInt("NchConstitPT4", N_chConstit_pT4);
      jetColl->back().addUserInt("NchConstitPT5", N_chConstit_pT5);
      jetColl->back().addUserInt("NchConstitPT10", N_chConstit_pT10);
      jetColl->back().addUserInt("NchConstitPT15", N_chConstit_pT15);
      jetColl->back().addUserInt("NchConstitPT20", N_chConstit_pT20);
      jetColl->back().addUserInt("NchConstitPT25", N_chConstit_pT25);
      jetColl->back().addUserInt("NchConstitPT30", N_chConstit_pT30);
      jetColl->back().addUserInt("NchConstitPT35", N_chConstit_pT35);
      jetColl->back().addUserInt("NchConstitPT40", N_chConstit_pT40);
      jetColl->back().addUserInt("NchConstitPT45", N_chConstit_pT45);
      jetColl->back().addUserInt("NchConstitPT50", N_chConstit_pT50);

      jetColl->back().addUserInt("NneutrConstitPT0", N_neutrConstit_pT0);
      jetColl->back().addUserInt("NneutrConstitPT1", N_neutrConstit_pT1);
      jetColl->back().addUserInt("NneutrConstitPT2", N_neutrConstit_pT2);
      jetColl->back().addUserInt("NneutrConstitPT3", N_neutrConstit_pT3);
      jetColl->back().addUserInt("NneutrConstitPT4", N_neutrConstit_pT4);
      jetColl->back().addUserInt("NneutrConstitPT5", N_neutrConstit_pT5);
      jetColl->back().addUserInt("NneutrConstitPT10", N_neutrConstit_pT10);
      jetColl->back().addUserInt("NneutrConstitPT15", N_neutrConstit_pT15);
      jetColl->back().addUserInt("NneutrConstitPT20", N_neutrConstit_pT20);
      jetColl->back().addUserInt("NneutrConstitPT25", N_neutrConstit_pT25);
      jetColl->back().addUserInt("NneutrConstitPT30", N_neutrConstit_pT30);
      jetColl->back().addUserInt("NneutrConstitPT35", N_neutrConstit_pT35);
      jetColl->back().addUserInt("NneutrConstitPT40", N_neutrConstit_pT40);
      jetColl->back().addUserInt("NneutrConstitPT45", N_neutrConstit_pT45);
      jetColl->back().addUserInt("NneutrConstitPT50", N_neutrConstit_pT50);

      //DeltaR between the lead and sublead charged/neutral constituents of the jets
      float deltaR_chConstitLeadXSubLead = -10; 
      if(VecPtrJetChargedConstituents.size() > 1){ 
         const reco::Candidate* LeadChConstit    = VecPtrJetChargedConstituents[0];
         const reco::Candidate* SubLeadChConstit = VecPtrJetChargedConstituents[1];
         deltaR_chConstitLeadXSubLead = (float)deltaR(LeadChConstit->phi(),LeadChConstit->eta(),SubLeadChConstit->phi(),SubLeadChConstit->eta());
      }
      jetColl->back().addUserFloat("dRchConstitLeadXSubLead", deltaR_chConstitLeadXSubLead);
      float deltaR_neutrConstitLeadXSubLead = -10;
      if(VecPtrJetNeutralConstituents.size() > 1){
         const reco::Candidate* LeadNeutrConstit    = VecPtrJetNeutralConstituents[0];
         const reco::Candidate* SubLeadNeutrConstit = VecPtrJetNeutralConstituents[1];
         deltaR_neutrConstitLeadXSubLead = (float) deltaR(LeadNeutrConstit->phi(),LeadNeutrConstit->eta(),SubLeadNeutrConstit->phi(),                                                                                                                            SubLeadNeutrConstit->eta());
      } 
      jetColl->back().addUserFloat("dRneutrConstitLeadXSubLead", deltaR_neutrConstitLeadXSubLead);

      //Pt of the first, second, third and fourth charged/neutral constituent of the Jet
      float pT_LeadChConstit = 0;
      float pT_SubLeadChConstit = 0;
      float pT_ThirdChConstit = 0;
      float pT_FourthChConstit = 0;
      if(VecPtrJetChargedConstituents.size() > 0){
         const reco::Candidate* FirstChConstit    = VecPtrJetChargedConstituents[0];
         pT_LeadChConstit       = (float) FirstChConstit->pt();
         //std::cout << "PT of the first Charged Constituents of the Jet: " << pT_LeadChConstit << std::endl;
         if(VecPtrJetChargedConstituents.size() > 1){
            const reco::Candidate* SecondChConstit = VecPtrJetChargedConstituents[1];  
            pT_SubLeadChConstit    = (float) SecondChConstit->pt();
            //std::cout << "PT of the second Charged Constituents of the Jet: " << pT_SubLeadChConstit << std::endl;
            if(VecPtrJetChargedConstituents.size() > 2){
               const reco::Candidate* ThirdChConstit   = VecPtrJetChargedConstituents[2];
               pT_ThirdChConstit      = (float) ThirdChConstit->pt();
               //std::cout << "PT of the third Charged Constituents of the Jet: " << pT_ThirdChConstit << std::endl;
               if(VecPtrJetChargedConstituents.size() > 3){
                  const reco::Candidate* FourthChConstit  = VecPtrJetChargedConstituents[3];
                  pT_FourthChConstit     = (float) FourthChConstit->pt(); 
                  //std::cout << "PT of the fouth Charged Constituents of the Jet: " << pT_FourthChConstit << std::endl;
               }
            }
         }
      }
      jetColl->back().addUserFloat("PTLeadChConstit",pT_LeadChConstit);
      jetColl->back().addUserFloat("PTSubLeadChConstit",pT_SubLeadChConstit);
      jetColl->back().addUserFloat("PTThirdChConstit",pT_ThirdChConstit);
      jetColl->back().addUserFloat("PTFourthChConstit",pT_FourthChConstit);

      float pT_LeadNeutrConstit = 0;
      float pT_SubLeadNeutrConstit = 0;
      float pT_ThirdNeutrConstit = 0;
      float pT_FourthNeutrConstit = 0;
      if(VecPtrJetNeutralConstituents.size() > 0){
         const reco::Candidate* FirstNeutrConstit    = VecPtrJetNeutralConstituents[0];
         pT_LeadNeutrConstit       = (float) FirstNeutrConstit->pt();
         //std::cout << "PT of the first Neutral Constituents of the Jet: " << pT_LeadNeutrConstit << std::endl;
         if(VecPtrJetNeutralConstituents.size() > 1){
            const reco::Candidate* SecondNeutrConstit = VecPtrJetNeutralConstituents[1];
            pT_SubLeadNeutrConstit    = (float) SecondNeutrConstit->pt();
            //std::cout << "PT of the second Neutral Constituents of the Jet: " << pT_SubLeadNeutrConstit << std::endl;
            if(VecPtrJetNeutralConstituents.size() > 2){
               const reco::Candidate* ThirdNeutrConstit   = VecPtrJetNeutralConstituents[2];
               pT_ThirdNeutrConstit      = (float) ThirdNeutrConstit->pt();
               //std::cout << "PT of the third Neutral Constituents of the Jet: " << pT_ThirdNeutrConstit << std::endl;
               if(VecPtrJetNeutralConstituents.size() > 3){
                  const reco::Candidate* FourthNeutrConstit  = VecPtrJetNeutralConstituents[3];
                  pT_FourthNeutrConstit     = (float) FourthNeutrConstit->pt();
                  //std::cout << "PT of the fouth Neutral Constituents of the Jet: " << pT_FourthNeutrConstit << std::endl;
               }
            }
         }
      }
      jetColl->back().addUserFloat("PTLeadNeutrConstit",pT_LeadNeutrConstit);
      jetColl->back().addUserFloat("PTSubLeadNeutrConstit",pT_SubLeadNeutrConstit);
      jetColl->back().addUserFloat("PTThirdNeutrConstit",pT_ThirdNeutrConstit);
      jetColl->back().addUserFloat("PTFourthNeutrConstit",pT_FourthNeutrConstit);
     
 
      VecPtrJetChargedConstituents.clear();     
      VecPtrJetNeutralConstituents.clear(); 

   

   }

   iEvent.put(jetColl);
 
}

// ------------ method called once each job just before starting event loop  ------------
void 
ProducerCA8CHSwithQJetsForBoostedTaus::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ProducerCA8CHSwithQJetsForBoostedTaus::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
ProducerCA8CHSwithQJetsForBoostedTaus::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
ProducerCA8CHSwithQJetsForBoostedTaus::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
ProducerCA8CHSwithQJetsForBoostedTaus::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
ProducerCA8CHSwithQJetsForBoostedTaus::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ProducerCA8CHSwithQJetsForBoostedTaus::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

bool ProducerCA8CHSwithQJetsForBoostedTaus::sortJetConstituents(const reco::Candidate* cand1, const reco::Candidate* cand2){
  return (cand1->pt() > cand2->pt());
}

//define this as a plug-in
DEFINE_FWK_MODULE(ProducerCA8CHSwithQJetsForBoostedTaus);
