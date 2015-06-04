#include "TPaveText.h"
#include <TFile.h>
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TObject.h"
//#include "TIter.h"
#include "TList.h"
#include <TString.h>
#include <TH1.h>
#include <TF1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TMath.h>
#include <iostream>
#include <sstream>
#include <vector>
TH1D* integralHisto(TH1D* h1, TString option);

//weight for No PU samples
double wH400= (0.025658)*36.02/14396.;
//xsec 2.0267
//BR(H->ZZ) = 0.269
//BR(Z->mumu ) = 0.03366
//BR(Z->jj ) = 0.6991
//2.0267*0.269*2*0.03366*0.6991 = 0.025658

//0.0553 from jhu is for eejj, mumujj
//2.0267*0.269*2*0.067*0.6991 = 0.511

double wZZ= (5.9)*36.02/28560.;
double wZjets= (3048.)*36.02/1.;
double wtt= (157.5)*36.02/3000.;

void superImposeMany(){
  TFile* file = new TFile("dataJuneJuly_lowMass.root","read") ;
  TFile* fileOut = new TFile("plots150.root","recreate") ;
  int rebin=5;
  TIter next(file->GetListOfKeys());
  TH1D *h1; 
  while( h1=(TH1D *)next() ){
    std::cout<<h1->GetName()<<std::endl;
    // if((TString)(h1->GetName()) == (TString)("jetPt1_cuts_W"))
    //{
    TFile* file_data1 = new TFile("dataJuneJuly_lowMass.root","read");
    TFile* file_data2 = new TFile("dataMay_lowMass.root","read");
    TFile* file_H400 = new TFile("powheg150_lowMass.root","read");
    TFile* file_ZZ = new TFile("ZZPythia_lowMass.root","read");
    TFile* file_tt = new TFile("TTJetsMadgraph_lowMass.root","read");
    TFile* file_Zjets = new TFile("DYJetsMadgraph_lowMass.root","read") ;

    TH1D* h_data = (TH1D*)((TH1D*)file_data1->Get(h1->GetName()))->Clone((TString)((TString)(h1->GetName())+"_data"));
    TH1D* h_data2 = (TH1D*)((TH1D*)file_data2->Get(h1->GetName()))->Clone((TString)((TString)(h1->GetName())+"_data"));
    
    TH1D* h_H400 = (TH1D*)((TH1D*)file_H400->Get(h1->GetName()))->Clone((TString)((TString)h1->GetName()+"_H400"));
    TH1D* h_ZZ = (TH1D*)((TH1D*)file_ZZ->Get(h1->GetName()))->Clone((TString)((TString)h1->GetName()+"_ZZ"));
    TH1D* h_tt = (TH1D*)((TH1D*)file_tt->Get(h1->GetName()))->Clone((TString)((TString)h1->GetName()+"_tt"));
    TH1D* h_Zjets = (TH1D*)((TH1D*)file_Zjets->Get(h1->GetName()))->Clone((TString)((TString)h1->GetName()+"_Zjets"));

    TCanvas *c = new TCanvas((TString)((TString)(h1->GetName())+"_dataMC"),(TString)((TString)(h1->GetName())+"_dataMC"));
    h_data->Sumw2();
    h_data2->Sumw2();
    h_data->Add(h_data2);
    h_H400->Sumw2();
    h_ZZ->Sumw2();
    h_Zjets->Sumw2();
    h_tt->Sumw2();
  
   c->SetFillColor(kWhite);
   h_data->SetLineColor(2);
   //h_data->SetTitle("");
   h_data->Rebin(rebin);

   h_tt->SetLineColor(kGreen+1);
   h_tt->Rebin(rebin);
   h_tt->Scale(1.1/1.6);
   //h_tt->Scale(wtt);

    h_Zjets->SetLineColor(4);
    h_Zjets->Rebin(rebin);
    h_Zjets->Scale(1.1/1.6);
    //h_Zjets->Scale(wZjets);
    //h_Zjets->Add(h_tt);
 
    h_ZZ->SetLineColor(6);
    h_ZZ->Rebin(rebin);
    h_ZZ->Scale(1.1/1.6);
    //h_ZZ->Scale(wZZ);
    //h_ZZ->Add(h_Zjets);

    h_H400->Rebin(rebin);
    h_H400->Scale(1.1/1.6);
    //h_H400->Scale(wH400);
    //h_H400->Add(h_ZZ);

    h_data->Draw("EP");
    h_ZZ->Draw("HISTsame");
    h_tt->Add(h_ZZ);
    h_tt->Draw("HISTsame");
    h_Zjets->Add(h_tt);
    h_Zjets->Draw("HISTsame");
    h_H400->Draw("HISTsame");

    TLegend *leg= new TLegend(0.5,0.6,0.7,0.8);
    leg->AddEntry(h_data,"data","LP");
    leg->AddEntry(h_H400,"H150","L");
    leg->AddEntry(h_ZZ,"ZZ","L");
    leg->AddEntry(h_Zjets,"Z+jets","L"); 
    leg->AddEntry(h_tt,"ttbar","L");
    leg->SetFillColor(kWhite);
    leg->Draw("same");
 
    fileOut->cd();
    c->Write();
    }
}
//useful for PU, JES
void compare3versions(TString sample, TString name){
  TFile* file ;
  TFile* fileOut ;
  TFile* file_noPU ;
  TFile* file_PUL1Corr ;
  TFile* file_noL1Corr ;

  if(sample == "H400"){
    fileOut = new TFile("plots_PUeffect_H400_noZjjHiggsCuts.root","recreate") ;
    file_noPU = new TFile("PUstudy_H400_powheg_noPU_noZjjHiggsCuts.root","read");
    file_PUL1Corr = new TFile("PUstudy_H400_powheg_Fastjet_noZjjHiggsCuts.root","read");
    file_noL1Corr = new TFile("PUstudy_H400_powheg_noL1Corr_noZjjHiggsCuts.root","read");
    sample = sample +((TString)"_PUEffect");
  }
  if(sample == "H400Weighted"){
    fileOut = new TFile("plotsStep2_PUeffect_H400Weighted.root","recreate") ;
    file_noPU = new TFile("histosStep2_powheg500_nominal.root","read");
    file_PUL1Corr = new TFile("histosStep2_powheg500_PULessThan5.root","read");
    file_noL1Corr = new TFile("histosStep2_powheg500_PUMoreThan5.root","read");
    sample = sample +((TString)"_PUEffect");
 }
  
  if(sample == "ZJetsMADWeighted"){
    fileOut = new TFile("plots_PUeffect_ZJetsMADWeighted.root","recreate") ;
    file_noPU = new TFile("PUstudy_ZJetsMAD_noPU_powheg.root","read");
    file_PUL1Corr = new TFile("PUstudy_ZJetsMAD_FastjetWeighted_powheg.root","read");
    file_noL1Corr = new TFile("PUstudy_ZJetsMAD_noL1corrWeighted_powheg.root","read");
   sample = sample +((TString)"_PUEffect");
 }

  if(sample == "H400Weighted_JES"){
    fileOut = new TFile("plotsStep1_JESeffect_H400_powheg.root","recreate") ;
    file_noPU = new TFile("histosStep1_H400_powheg_fastjet.root","read");
    file_PUL1Corr = new TFile("histosStep1_H400_powheg_fastjet_JESUP.root","read");
    file_noL1Corr = new TFile("histosStep1_H400_powheg_fastjet_JESDOWN.root","read");
    sample = sample +((TString)"_JESEffect");
 }
  
  if(sample == "dataMay10"){
    fileOut = new TFile("plots_dataMay10_trigger.root","recreate") ;
    file_noPU = new TFile("dataMay10DoubleMu.root","read");
    file_PUL1Corr = new TFile("dataMay10DoubleMu_DoubleMu7.root","read");
    file_noL1Corr = new TFile("dataMay10DoubleMu_DoubleMu7.root","read");
    //file_noPU = new TFile("dataMay10DoubleElectron.root","read");
    //file_PUL1Corr = new TFile("dataMay10DoubleElectron_Ele17_8.root","read");
    //file_noL1Corr = new TFile("dataMay10DoubleElectron_Ele17_8.root","read");
   sample = sample +((TString)"_triggerMu");
 }
  
  if(sample == "dataPromptReco"){
    fileOut = new TFile("plots_dataMay10_trigger.root","recreate") ;
    //file_noPU = new TFile("dataPromptRecoDoubleElectron.root","read");
    //file_PUL1Corr = new TFile("dataPromptRecoDoubleElectron_Ele17_8.root","read");
    //file_noL1Corr = new TFile("dataPromptRecoDoubleElectron_Ele17_8.root","read");
    file_noPU = new TFile("dataPromptRecoDoubleMu.root","read");
    file_PUL1Corr = new TFile("dataPromptRecoDoubleMu_DoubleMu7_Mu13_8.root","read");
    file_noL1Corr = new TFile("dataPromptRecoDoubleMu_DoubleMu7_Mu13_8.root","read");
   sample = sample +((TString)"_triggerMu");
 }
  
  int rebin=1;
  TIter next(file_noPU->GetListOfKeys());
  TH1D *h1; 
  while( h1=(TH1D *)next() ){
    std::cout<<h1->GetName()<<std::endl;
    if(((TString)h1->GetName())== name || name =="") {

      rebin=1;
    TH1D* h_noPU = (TH1D*)((TH1D*)file_noPU->Get(h1->GetName()))->Clone((TString)((TString)h1->GetName()+"_noPU"));
    TH1D* h_PUL1Corr = (TH1D*)((TH1D*)file_PUL1Corr->Get(h1->GetName()))->Clone((TString)((TString)h1->GetName()+"_PUL1Corr"));
    TH1D* h_noL1Corr = (TH1D*)((TH1D*)file_noL1Corr->Get(h1->GetName()))->Clone((TString)((TString)h1->GetName()+"_noL1Corr"));
    
    TString name = h1->GetName()+((TString)"_");

    TCanvas *c = new TCanvas((TString)(name+sample),(TString)(name+sample));
    h_noPU->Sumw2();
    h_PUL1Corr->Sumw2();
    h_noL1Corr->Sumw2();
    h_noPU->Rebin(rebin);
    h_PUL1Corr->Rebin(rebin);
    h_noL1Corr->Rebin(rebin);

    //h_noPU->Scale(1./h_noPU->Integral());
    //h_PUL1Corr->Scale(1./h_PUL1Corr->Integral());
    //h_noL1Corr->Scale(1./h_noL1Corr->Integral());

    h_noPU->SetLineWidth(2);
    h_PUL1Corr->SetLineColor(2);
    h_noL1Corr->SetLineColor(4);
  
    TPad* spad2 = new TPad("spad2","The second subpad",.0,.3,1.,1.);
    spad2->cd();
    h_noPU->GetYaxis()->SetTitle("Events / 10 GeV");
    h_noPU->Draw("EHIST");
    h_noPU->SetTitle("");
    h_PUL1Corr->Draw("EHISTsame");
    //h_noL1Corr->Draw("EHISTsame");

    spad2->SetFillColor(kWhite);

    TLegend *leg= new TLegend(0.5,0.6,0.7,0.8);
    if(sample.Contains("_PUEffect")){
      /*leg->AddEntry(h_noPU,"no PU","L");
      leg->AddEntry(h_PUL1Corr,"PU corrected","L");
      leg->AddEntry(h_noL1Corr,"with PU","L");*/
       leg->AddEntry(h_noPU,"nominal","L");
      leg->AddEntry(h_PUL1Corr,"less than 5","L");
      leg->AddEntry(h_noL1Corr,"more than 5","L");
    }
    if(sample.Contains("_JESEffect")){
      leg->AddEntry(h_noPU,"central value","L");
      leg->AddEntry(h_PUL1Corr,"JES up","L");
      leg->AddEntry(h_noL1Corr,"JES down","L");
    }
    if(sample.Contains("_trigger")){
      leg->AddEntry(h_noPU,"no trigger requirement","L");
      //leg->AddEntry(h_PUL1Corr,"HLT_DoubleMu7_v*","L");
      leg->AddEntry(h_PUL1Corr,"HLT_DoubleMu7_v* OR HLT_Mu13_Mu8_v*","L");
      //leg->AddEntry(h_PUL1Corr,"HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v*","L");
    }
    leg->SetFillColor(kWhite);
    leg->Draw("same");

    TPad* spad1 = new TPad("spad1","The first subpad",.0,.0,1.,.3);
    spad1->cd();
    TH1D* h_PUL1Corr_ratio=(TH1D*) (h_PUL1Corr->Clone((TString)h1->GetName()+"L1corrRatio"));
    TH1D* h_noL1Corr_ratio=(TH1D*) (h_noL1Corr->Clone((TString)h1->GetName()+"PURatio"));
    TH1D* h_noPU_ratio=(TH1D*) (h_noPU->Clone((TString)h1->GetName()+"noPURatio"));
    //h_PUL1Corr_ratio->Rebin(2);
    //h_noL1Corr_ratio->Rebin(2);
    //h_noPU_ratio->Rebin(2);
    h_PUL1Corr_ratio->Sumw2();
    h_noL1Corr_ratio->Sumw2();
    h_noPU_ratio->Sumw2();
    
    h_PUL1Corr_ratio->Divide(h_noPU_ratio);
    h_noL1Corr_ratio->Divide(h_noPU_ratio);
    if(sample.Contains("_PUEffect"))
      h_PUL1Corr_ratio->GetYaxis()->SetTitle("ratio to noPU");
    if(sample.Contains("_JESEffect"))
      h_PUL1Corr_ratio->GetYaxis()->SetTitle("ratio to central value");
     h_PUL1Corr_ratio->GetYaxis()->SetTitleOffset(0.5);
    h_PUL1Corr_ratio->GetYaxis()->SetTitleSize(0.08);
    h_PUL1Corr_ratio->Draw("EHIST");
    h_PUL1Corr_ratio->SetTitle("");
    h_PUL1Corr_ratio->GetXaxis()->SetLabelSize(0.08);
    h_PUL1Corr_ratio->GetYaxis()->SetLabelSize(0.08);
    //h_noL1Corr_ratio->Draw("EHISTsame");
    h_noPU_ratio->Divide(h_noPU_ratio);
    h_noPU_ratio->Draw("HISTsame");

    spad1->SetFillColor(kWhite);

    c->cd();
    spad1->Draw("same");
    spad2->Draw("same");
    fileOut->cd();
    c->Write();

    
    }
  }
}


void comparejhu_powheg(){
  TFile* file ;
  TFile* fileOut ;
  TFile* file_noPU ;
  TFile* file_noPU_jhu ;

    file = new TFile("h300ana_jhu_noPU.root","read") ;
    fileOut = new TFile("plots_jhuVSpowheg_H300.root","recreate") ;
    file_noPU = new TFile("h300ana_powheg_lowPU_L1fj.root","read");
    file_noPU_jhu = new TFile("h300ana_jhu_noPU.root","read");

  int rebin=2;
  TIter next(file->GetListOfKeys());
  TH1D *h1; 
  while( h1=(TH1D *)next() ){
    std::cout<<h1->GetName()<<std::endl;
    //    if((TString)(h1->GetName()) == (TString)("jetPt1_cuts_W"))
    if(((TString)(h1->GetName())).Contains((TString)("gen")))
    {
    
    TH1D* h_noPU = (TH1D*)((TH1D*)file_noPU->Get(h1->GetName()))->Clone((TString)((TString)h1->GetName()+"_noPU"));
    TH1D* h_noPU_jhu = (TH1D*)((TH1D*)file_noPU_jhu->Get(h1->GetName()))->Clone((TString)((TString)h1->GetName()+"_noPU_jhu"));
    
    TString name = h1->GetName()+((TString)"_jhuVSpowheg");

    TCanvas *c = new TCanvas(name,name);
    c->SetFillColor(kWhite);
    h_noPU->Sumw2();
    h_noPU_jhu->Sumw2();

    h_noPU->Rebin(rebin);
    h_noPU_jhu->Rebin(rebin);

    h_noPU->Scale(1./h_noPU->Integral());
    h_noPU_jhu->Scale(1./h_noPU_jhu->Integral());

    h_noPU->SetLineColor(2);
 
    TPad* spad2 = new TPad("spad2","The second subpad",.0,.3,1.,1.);
    spad2->cd();
    spad2->SetFillColor(kWhite);

    h_noPU->Draw("EHIST");
    h_noPU_jhu->Draw("EHISTsame");


    TLegend *leg= new TLegend(0.5,0.6,0.7,0.8);
    //leg->AddEntry(h_noPU,"powheg (low PU)","L");
    //leg->AddEntry(h_noPU_jhu,"JHU generator (noPU)","L");
    leg->AddEntry(h_noPU,"powheg","L");
    leg->AddEntry(h_noPU_jhu,"JHU generator","L");
    leg->SetFillColor(kWhite);
    leg->Draw("same");
 
     TPad* spad1 = new TPad("spad1","The first subpad",.0,.0,1.,.3);
    spad1->cd();
    TH1D* h_noPU_jhu_ratio=(TH1D*) (h_noPU_jhu->Clone((TString)h1->GetName()+"jhuRatio"));
    TH1D* h_noPU_ratio=(TH1D*) (h_noPU->Clone((TString)h1->GetName()+"noPURatio"));
    //h_noPU_jhu_ratio->Rebin(2);
    //h_noL1Corr_ratio->Rebin(2);
    //h_noPU_ratio->Rebin(2);
    h_noPU_jhu_ratio->Sumw2();
    h_noPU_ratio->Sumw2();
    
    h_noPU_jhu_ratio->Divide(h_noPU_ratio);
    h_noPU_jhu_ratio->GetYaxis()->SetTitle("ratio to powheg");
    h_noPU_jhu_ratio->GetYaxis()->SetTitleOffset(0.5);
    h_noPU_jhu_ratio->GetYaxis()->SetTitleSize(0.08);
    h_noPU_jhu_ratio->Draw("EHIST");
    h_noPU_jhu_ratio->SetTitle("");
    h_noPU_jhu_ratio->GetXaxis()->SetLabelSize(0.08);
    h_noPU_jhu_ratio->GetYaxis()->SetLabelSize(0.08);

    h_noPU_ratio->Divide(h_noPU_ratio);
    h_noPU_ratio->Draw("HISTsame");

    spad1->SetFillColor(kWhite);

    c->cd();
    spad1->Draw("");
    spad2->Draw("same");
  fileOut->cd();
    c->Write();
   }
  }
}

void comparejhu_powheg_integral(TString histoName, TString option){
  TFile* file ;
  TFile* fileOut ;
  TFile* file_noPU ;
  TFile* file_noPU_jhu ;

    file = new TFile("h300ana_jhu_noPU.root","read") ;
    fileOut = new TFile("plots_jhuVSpowheg_H300.root","recreate") ;
    file_noPU = new TFile("h300ana_powheg_lowPU_L1fj.root","read");
    file_noPU_jhu = new TFile("h300ana_jhu_noPU.root","read");

  int rebin=1;
  TIter next(file->GetListOfKeys());
  TH1D *h1; 
    TH1D* h_noPU = (TH1D*)((TH1D*)file_noPU->Get(histoName))->Clone((TString)((TString)histoName+"_noPU"));
    TH1D* h_noPU_jhu = (TH1D*)((TH1D*)file_noPU_jhu->Get(histoName))->Clone((TString)((TString)histoName+"_noPU_jhu"));
    
    TString name = histoName+((TString)"_jhuVSpowheg_integral");

    TCanvas *c = new TCanvas(name,name);
    c->SetFillColor(kWhite);
    h_noPU->Sumw2();
    h_noPU_jhu->Sumw2();

    h_noPU->Rebin(rebin);
    h_noPU_jhu->Rebin(rebin);

    h_noPU->Scale(1./h_noPU->Integral(0,h_noPU->GetXaxis()->GetNbins()+1));
    h_noPU_jhu->Scale(1./h_noPU_jhu->Integral(0,h_noPU_jhu->GetXaxis()->GetNbins()+1));

    h_noPU->SetLineColor(2);
    h_noPU->SetTitle("");
    h_noPU_jhu->SetTitle("");

    TH1D* h_noPU_integral = integralHisto(h_noPU,option);
    TH1D* h_noPU_jhu_integral = integralHisto(h_noPU_jhu,option);
    
    h_noPU_integral->Draw("EHIST");
    h_noPU_jhu_integral->Draw("EHISTsame");


    TLegend *leg= new TLegend(0.5,0.6,0.7,0.8);
    //leg->AddEntry(h_noPU,"powheg (low PU)","L");
    //leg->AddEntry(h_noPU_jhu,"JHU generator (noPU)","L");
    leg->AddEntry(h_noPU_integral,"powheg","L");
    leg->AddEntry(h_noPU_jhu_integral,"JHU generator","L");
    leg->SetFillColor(kWhite);
    leg->Draw("same");
 
   
    fileOut->cd();
    c->Write();
}



TH1D* integralHisto(TH1D* h1, TString option){
  TH1D* hIntegro=(TH1D*) (h1->Clone((TString)h1->GetName()+"_integro"));
  for(int i=0; i<h1->GetXaxis()->GetNbins()+2; i++){
    hIntegro->SetBinContent(i,0);
  }
 if(option == "fromLeft")
    for(int i=0; i<=h1->GetXaxis()->GetNbins(); i++){
      hIntegro->SetBinContent(i+1,h1->Integral(0,i+1));
    }
  else if(option == "fromRight")
    for(int i=h1->GetXaxis()->GetNbins(); i>-1; i--){
      cout<<"bin "<<i<<" integral from "<<i<<" to "<<h1->GetXaxis()->GetNbins()+1<<" = "
	  <<h1->Integral(i,h1->GetXaxis()->GetNbins()+1)<<"  total = "<<h1->Integral()<<endl;
      hIntegro->SetBinContent(i,h1->Integral(i,(h1->GetXaxis()->GetNbins()+1)));
    }
  else if(option == "fromBoth")
    for(int i=0; i<(h1->GetXaxis()->GetNbins()/2); i++){
      double integral = h1->Integral(0,i+1);
      integral = integral+h1->Integral(h1->GetXaxis()->GetNbins()-i,h1->GetXaxis()->GetNbins()+1);
      hIntegro->SetBinContent(i+1,integral);
    }
  else if(option == "fromInside"){
    int A= h1->FindBin(300);
    for(int i=0; i<=h1->GetXaxis()->GetNbins(); i++){
      if(A-i>-1 && A+i<h1->GetXaxis()->GetNbins()+1){
	hIntegro->SetBinContent(A+i,h1->Integral(A-i,A+i));
	hIntegro->SetBinContent(A-i,h1->Integral(A-i,A+i));
      }
      else if(A-i<0 && A+i<h1->GetXaxis()->GetNbins()+2){
	hIntegro->SetBinContent(A+i,h1->Integral(0,A+i));
	hIntegro->SetBinContent(A-i,h1->Integral(0,A+i));
      }
      if(A-i>-1 && A+i>h1->GetXaxis()->GetNbins()+1){
	hIntegro->SetBinContent(A+i,h1->Integral(A-i,h1->GetXaxis()->GetNbins()+1));
	hIntegro->SetBinContent(A-i,h1->Integral(A-i,h1->GetXaxis()->GetNbins()+1));
      }
    }
  }
  return hIntegro;
}

void extractPUWeight(){

    TFile* file_MC= new TFile("/afs/cern.ch/user/s/sbologne/scratch0/CMSSW/CMGTools/CMSSW_4_2_3/src/HiggsAna/HLLJJCommon/data/PU_MC_Summer11.root","read") ;
    TH1D* h_MC = (TH1D*)((TH1D*)file_MC->Get("PU"))->Clone("MC_clone");
    TFile* file_data= new TFile("/afs/cern.ch/user/s/sbologne/scratch0/CMSSW/CMGTools/CMSSW_4_2_3/src/HiggsAna/HLLJJCommon/data/Pileup_2011_to_172802_LP_LumiScale.root","read") ;
    TH1D* h_data = (TH1D*)((TH1D*)file_data->Get("pileup"))->Clone("data_clone");
    TH1D* h_ratio = (TH1D*)((TH1D*)file_data->Get("pileup"))->Clone("ratio_clone");

    TCanvas *c = new TCanvas("PUDist","PUDist");
    c->SetFillColor(kWhite);
    h_MC->Scale(1./h_MC->Integral());
    h_MC->Draw();
    h_data->SetLineColor(2);
    h_data->Scale(1./h_data->Integral());
    h_data->Draw("same");
    TLegend *leg= new TLegend(0.5,0.6,0.7,0.8);
    leg->AddEntry(h_MC,"MC","L");
    leg->AddEntry(h_data,"data","L");
    leg->SetFillColor(kWhite);
    leg->Draw("same");
  

    h_ratio->Scale(1./h_ratio->Integral());
    h_ratio->Divide(h_MC);
    TCanvas *c2 = new TCanvas("PUratio","PUratio");
    h_ratio->Draw();

    double totWeight=0;
    for(int i=1;i<52;i++){
      cout<<h_data->GetBinCenter(i)<<"   ";
      cout<<h_data->GetBinContent(i)/h_MC->GetBinContent(i)<<endl;
      totWeight=totWeight+h_data->GetBinContent(i)/h_MC->GetBinContent(i);
    }
    cout<<"sum of weight "<<totWeight<<endl;
}

void computeEffErr(double num, double den){
  double eff = num/den;
  double err = sqrt(eff*(1-eff)/den);
  cout<<"eff "<<eff<<" +/- "<<err<<endl;
}

void computeAccept(TString name, double startingEvents){				
  TFile *file = new TFile(name,"read");
  TH1D* h = (TH1D*)((TH1D*)file->Get("hfinalHiggsNumber"))->Clone("clone");
  double eventsAfterPreselectrion = h->Integral(2,2);
  computeEffErr(eventsAfterPreselectrion,startingEvents);
}

void computeAccept(){
  /*computeAccept("histosStep2_powheg500_nominal2.root",109982/3.);
  computeAccept("histosStep2_powheg500_JESUP2.root",109982/3.);
  computeAccept("histosStep2_powheg500_JESDOWN2.root",109982/3.);*/
  computeAccept("histosStep2_powheg500_nominal2.root",109982/3.);
  computeAccept("histosStep2_powheg500_JESUP2.root",109982/3.);
  computeAccept("histosStep2_powheg500_JESDOWN2.root",109982/3.);
}

vector<double> extractNEvents(TString fileName){
  cout<<fileName<<endl;
  vector<double> DYresults;
  TFile* file = new TFile(fileName,"read") ;
  TH1D* hMass = (TH1D*)((TH1D*)file->Get("hMassBestCand"));
  TH1D* hMassMu = (TH1D*)((TH1D*)file->Get("hMassBestCand_Mu"));
  TH1D* hMassEle = (TH1D*)((TH1D*)file->Get("hMassBestCand_Ele"));
  TH1D* hMassMu0btag = (TH1D*)((TH1D*)file->Get("hMassBestCand_Mu0btag"));
  TH1D* hMassEle0btag = (TH1D*)((TH1D*)file->Get("hMassBestCand_Ele0btag"));
  TH1D* hMassMu1btag = (TH1D*)((TH1D*)file->Get("hMassBestCand_Mu1btag"));
  TH1D* hMassEle1btag = (TH1D*)((TH1D*)file->Get("hMassBestCand_Ele1btag"));
  TH1D* hMassMu2btag = (TH1D*)((TH1D*)file->Get("hMassBestCand_Mu2btag"));
  TH1D* hMassEle2btag = (TH1D*)((TH1D*)file->Get("hMassBestCand_Ele2btag"));
  cout<<hMassMu0btag<<" "<<hMassEle0btag<<endl;
  /*cout<<"MUON "<<endl;
    cout<<"total "<<hMassMu->Integral()<<":  0 btag "<<hBCatMu->Integral(1,1)
    <<"    1 btag "<<hBCatMu->Integral(2,2)<<"    2 btag "<<hBCatMu->Integral(3,3)<<endl;
    cout<<"ELE"<<endl;
    cout<<"total "<<hMassEle->Integral()<<":  0 btag "<<hBCatEle->Integral(1,1)
    <<"    1 btag "<<hBCatEle->Integral(2,2)<<"    2 btag "<<hBCatEle->Integral(3,3)<<endl;*/
  double mu0btag = hMassMu0btag->Integral(140,160)/1.6 * 1.1;
  double mu1btag = hMassMu1btag->Integral(140,160)/1.6 * 1.1;
  double mu2btag = hMassMu2btag->Integral(140,160)/1.6 * 1.1;
  //double mu = hMassMu->Integral(14,16)/1.6 * 1.1;
  double mu = mu0btag+mu1btag+mu2btag;
  double ele0btag = hMassEle0btag->Integral(140,160)/1.6 * 1.1;
  double ele1btag = hMassEle1btag->Integral(140,160)/1.6 * 1.1;
  double ele2btag = hMassEle2btag->Integral(140,160)/1.6 * 1.1;
  //double ele = hMassMu->Integral(14,16)/1.6 * 1.1;
  double ele = ele0btag+ele1btag+ele2btag;
  double tot = hMass->Integral(140,160)/1.6 * 1.1;
  DYresults.push_back(mu0btag);
  DYresults.push_back(mu1btag);
  DYresults.push_back(mu2btag);
  DYresults.push_back(mu);
  DYresults.push_back(ele0btag);
  DYresults.push_back(ele1btag);
  DYresults.push_back(ele2btag);
  DYresults.push_back(ele);
  DYresults.push_back(tot);
  return DYresults;
}

void printNevents(){
  vector<double> DYresults = extractNEvents("DYJetsMadgraph_lowMass.root");
  vector<double> TTresults = extractNEvents("TTJetsMadgraph_lowMass.root");
  vector<double> ZZresults = extractNEvents("ZZPythia_lowMass.root");
  vector<double> powheg150results = extractNEvents("powheg150_lowMass.root");
  vector<double> dataJuneJulyresults = extractNEvents("dataJuneJuly_lowMass.root");
  vector<double> dataMayresults = extractNEvents("dataMay_lowMass.root");
  
  cout<<"            Signal  |  Z+jets |  TT+jets | ZZ | Background | Tot MC |  DATA"<<endl;
  cout<<"mu 0 btag   "<<powheg150results[0]<<"  |  "<<DYresults[0]<<"  |  "<<TTresults[0]
      <<"  |  "<<ZZresults[0]<<"  |  "<<DYresults[0]+TTresults[0]+ZZresults[0]
      <<"  |  "<<DYresults[0]+TTresults[0]+ZZresults[0]+powheg150results[0]
      <<" | "<<dataJuneJulyresults[0]+dataMayresults[0]<<endl;
  cout<<"mu 1 btag   "<<powheg150results[1]<<"  |  "<<DYresults[1]<<"  |  "<<TTresults[1]
      <<"  |  "<<ZZresults[1]<<"  |  "<<DYresults[1]+TTresults[1]+ZZresults[1]
      <<"  |  "<<DYresults[1]+TTresults[1]+ZZresults[1]+powheg150results[1]
      <<" | "<<dataJuneJulyresults[1]+dataMayresults[1]<<endl;

  cout<<"mu 2 btag   "<<powheg150results[2]<<"  |  "<<DYresults[2]<<"  |  "<<TTresults[2]
      <<"  |  "<<ZZresults[2]<<"  |  "<<DYresults[2]+TTresults[2]+ZZresults[2]
      <<"  |  "<<DYresults[2]+TTresults[2]+ZZresults[2]+powheg150results[2]
      <<" | "<<dataJuneJulyresults[2]+dataMayresults[2]<<endl;

  cout<<"mu          "<<powheg150results[3]<<"  |  "<<DYresults[3]<<"  |  "<<TTresults[3]
      <<"  |  "<<ZZresults[3]<<"  |  "<<DYresults[3]+TTresults[3]+ZZresults[3]
      <<"  |  "<<DYresults[3]+TTresults[3]+ZZresults[3]+powheg150results[3]
      <<" | "<<dataJuneJulyresults[3]+dataMayresults[3]<<endl;
  cout<<"ele 0 btag   "<<powheg150results[4]<<"  |  "<<DYresults[4]<<"  |  "<<TTresults[4]
      <<"  |  "<<ZZresults[4]<<"  |  "<<DYresults[4]+TTresults[4]+ZZresults[4]
      <<"  |  "<<DYresults[4]+TTresults[4]+ZZresults[4]+powheg150results[4]
      <<" | "<<dataJuneJulyresults[4]+dataMayresults[4]<<endl;
 cout<<"ele 1 btag   "<<powheg150results[5]<<"  |  "<<DYresults[5]<<"  |  "<<TTresults[5]
      <<"  |  "<<ZZresults[5]<<"  |  "<<DYresults[5]+TTresults[5]+ZZresults[5]
      <<"  |  "<<DYresults[5]+TTresults[5]+ZZresults[5]+powheg150results[5]
      <<" | "<<dataJuneJulyresults[5]+dataMayresults[5]<<endl;

  cout<<"ele 2 btag   "<<powheg150results[6]<<"  |  "<<DYresults[6]<<"  |  "<<TTresults[6]
      <<"  |  "<<ZZresults[6]<<"  |  "<<DYresults[6]+TTresults[6]+ZZresults[6]
      <<"  |  "<<DYresults[6]+TTresults[6]+ZZresults[6]+powheg150results[6]
      <<" | "<<dataJuneJulyresults[6]+dataMayresults[6]<<endl;

  cout<<"ele          "<<powheg150results[7]<<"  |  "<<DYresults[7]<<"  |  "<<TTresults[7]
      <<"  |  "<<ZZresults[7]<<"  |  "<<DYresults[7]+TTresults[7]+ZZresults[7]
      <<"  |  "<<DYresults[7]+TTresults[7]+ZZresults[7]+powheg150results[7]
      <<" | "<<dataJuneJulyresults[7]+dataMayresults[7]<<endl;
  cout<<"tot          "<<powheg150results[8]<<"  |  "<<DYresults[8]<<"  |  "<<TTresults[8]
      <<"  |  "<<ZZresults[8]<<"  |  "<<DYresults[8]+TTresults[8]+ZZresults[8]
      <<"  |  "<<DYresults[8]+TTresults[8]+ZZresults[8]+powheg150results[8]
      <<" | "<<dataJuneJulyresults[8]+dataMayresults[8]<<endl;
}

