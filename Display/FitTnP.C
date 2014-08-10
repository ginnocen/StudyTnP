#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TH1D.h"
#define NUM_BX 9000


void FitTnP(){

  const int nMuPtBin = 7;
  double mumulow=2.;
  double mumuhigh=4.;
  
  TH1D* hTrigPtPass[nMuPtBin];
  TH1D* hTrigPtFail[nMuPtBin];
  TH1D* hTrigPtAll[nMuPtBin];
  
  TH1D* hTrkPtPass[nMuPtBin];
  TH1D* hTrkPtFail[nMuPtBin];
  TH1D* hTrkPtAll[nMuPtBin];
  
  TH1D* hMuonIdPtPass[nMuPtBin];
  TH1D* hMuonIdPtFail[nMuPtBin];
  TH1D* hMuonIdPtAll[nMuPtBin];
  
  TF1* tf1mumu_trg_all[nMuPtBin];
  TF1* tf1mumu_trg_pass[nMuPtBin];
  TF1* tf1mumu_trk_all[nMuPtBin];

  for(int i = 0; i < nMuPtBin; i++){
  }


  
  TFile*finput=new TFile("../Code/Results/foutput.root","read");
  finput->cd();
  for(int i = 0; i < nMuPtBin; i++){
    hTrigPtPass[i]=(TH1D*)finput->Get(Form("hTrigPtPass%d",i));
    hTrigPtFail[i]=(TH1D*)finput->Get(Form("hTrigPtFail%d",i));
    hTrigPtAll[i]=(TH1D*)finput->Get(Form("hTrigPtAll%d",i));
    
  }
  
  TF1 *signal = new TF1("Gaus(x,[0],[1])","gaus");  
  signal->SetParameter(0,3.1);
  hTrigPtPass[0]->Fit("signal", "L q", "",2.8,3.2);

  
  TCanvas*ctest=new TCanvas("ctest","ctest",650,600);
  ctest->cd();
  hTrigPtPass[0]->Draw();
  signal->Draw("same");
  ctest->SaveAs("ctest.eps");

}


