#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include "TH1F.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TLegend.h"

void mergeTree(){
  TFile *fmu = new TFile("../../HLLJJCommon/test/summer11_data_mu_lowmass.root","read");
  TFile *fele = new TFile("../../HLLJJCommon/test/summer11_data_ele_lowmass.root","read");
 TTree* oldtmu = (TTree*) fmu->Get("AngularInfo");
 oldtmu->SetName("AngularInfoMuOld");
 oldtmu->SetBranchStatus("*",1);
 TTree* oldtele = (TTree*) fele->Get("AngularInfo");
 oldtele->SetName("AngularInfoEleOld");
 oldtele->SetBranchStatus("*",1);
 
 TFile *fmuNew = new TFile("tempMu.root","recreate");
 TTree* newtmu = oldtmu->CopyTree("lep==1");
 newtmu->SetName("AngularInfo");
 newtmu->Write();
 TFile *feleNew = new TFile("tempEle.root","recreate");
 TTree* newtele = oldtele->CopyTree("lep==0");
 newtele->SetName("AngularInfo");
 newtele->Write();

 TChain *chain= new TChain("AngularInfo","both leptons, full statistics");
 chain->Add("tempMu.root");
 chain->Add("tempEle.root");
 chain->Merge("../../HLLJJCommon/test/summer11_data_lowmass.root");
 //TFile *fNew = new TFile("../../HLLJJCommon/test/summer11_data_lowmass.root","recreate");
 //chain->Write();
}
