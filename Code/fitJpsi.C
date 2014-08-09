#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "utilities.h"
#include "TH1D.h"
#define NUM_BX 9000


void fitJpsi(){

  const int nMuPtBin = 7;
  const double MuPtBin[nMuPtBin+1] = {0.0,1.5,3.0,4.5,6.0,9.0,20.0,30.0};
  
  int nBin = 100;
  double mumulow=2.;
  double mumuhigh=4.;
  
  bool IsMuonInAcceptance(Float_t,Float_t,Float_t);
  
  TString infname_mc="/data/bmeson/MC_Jpsi/jpsi-Kp.root";
  TFile *inf_mc = new TFile(infname_mc.Data());
  
  TTree *ntuple = (TTree*) inf_mc->Get("ntJpsi");
  TTree *nt_mcgen = (TTree*)inf_mc->Get("ntGen");
  ntuple->AddFriend(nt_mcgen);
  
  TH1D* hTrigPtPass[nMuPtBin];
  TH1D* hTrigPtFail[nMuPtBin];
  TH1D* hTrigPtAll[nMuPtBin];
  
  TH1D* hTrkPtPass[nMuPtBin];
  TH1D* hTrkPtFail[nMuPtBin];
  TH1D* hTrkPtAll[nMuPtBin];
  
  TH1D* hMuonIdPtPass[nMuPtBin];
  TH1D* hMuonIdPtFail[nMuPtBin];
  TH1D* hMuonIdPtAll[nMuPtBin];

  for(int i = 0; i < nMuPtBin; i++){
    hTrigPtPass[i] = new TH1D(Form("hTrigPtPass%d",i),"",nBin,mumulow,mumuhigh);
    hTrigPtFail[i] = new TH1D(Form("hTrigPtFail%d",i),"",nBin,mumulow,mumuhigh);
    hTrigPtAll[i] = new TH1D(Form("hTrigPtAll%d",i),"",nBin,mumulow,mumuhigh);
    hTrkPtPass[i] = new TH1D(Form("hTrkPtPass%d",i),"",nBin,mumulow,mumuhigh);
    hTrkPtFail[i] = new TH1D(Form("hTrkPtFail%d",i),"",nBin,mumulow,mumuhigh);
    hTrkPtAll[i] = new TH1D(Form("hTrkPtAll%d",i),"",nBin,mumulow,mumuhigh);
    hMuonIdPtPass[i] = new TH1D(Form("hMuonIdPtPass%d",i),"",nBin,mumulow,mumuhigh);
    hMuonIdPtFail[i] = new TH1D(Form("hMuonIdPtFail%d",i),"",nBin,mumulow,mumuhigh);
    hMuonIdPtAll[i] = new TH1D(Form("hMuonIdPtAll%d",i),"",nBin,mumulow,mumuhigh);
  }
  
  int Run,size,Event;
  int isTriggered1[NUM_BX];
  int isTriggered2[NUM_BX];
  Float_t pt[NUM_BX];
  Float_t mass[NUM_BX];
  Float_t genpt[NUM_BX];
  Float_t eta[NUM_BX];
  Float_t y[NUM_BX];
  Float_t phi[NUM_BX];
  int isTracker1[NUM_BX];
  int isTracker2[NUM_BX];
  Float_t pt1[NUM_BX];
  Float_t pt2[NUM_BX];
  Float_t eta1[NUM_BX];
  Float_t eta2[NUM_BX];
  Float_t phi1[NUM_BX];
  Float_t phi2[NUM_BX];
  int id1[NUM_BX];
  int id2[NUM_BX];
  Bool_t outerTrackisNonnull1[NUM_BX];
  Bool_t outerTrackisNonnull2[NUM_BX];
  int gen[NUM_BX];

  ntuple->SetBranchAddress("Run",&Run);
  ntuple->SetBranchAddress("Event",&Event);
  ntuple->SetBranchAddress("size",&size);
  ntuple->SetBranchAddress("mass",mass);
  ntuple->SetBranchAddress("pt",pt);
  ntuple->SetBranchAddress("genpt",genpt);
  ntuple->SetBranchAddress("eta",eta);
  ntuple->SetBranchAddress("y",y);
  ntuple->SetBranchAddress("phi",phi);
  ntuple->SetBranchAddress("isTracker1",isTracker1);
  ntuple->SetBranchAddress("isTracker2",isTracker2);
  ntuple->SetBranchAddress("pt1",pt1);
  ntuple->SetBranchAddress("pt2",pt2);
  ntuple->SetBranchAddress("eta1",eta1);
  ntuple->SetBranchAddress("eta2",eta2);
  ntuple->SetBranchAddress("phi1",phi1);
  ntuple->SetBranchAddress("phi2",phi2);
  ntuple->SetBranchAddress("id1",id1);
  ntuple->SetBranchAddress("id2",id2);
  ntuple->SetBranchAddress("outerTrackisNonnull1",outerTrackisNonnull1);
  ntuple->SetBranchAddress("outerTrackisNonnull2",outerTrackisNonnull2);
  ntuple->SetBranchAddress("gen",gen);
  ntuple->SetBranchAddress("isTriggered1",isTriggered1);
  ntuple->SetBranchAddress("isTriggered2",isTriggered2);

  
  Int_t entries = (Int_t)ntuple->GetEntries();
  
  for (int i=0; i<entries; i++) {
    ntuple->GetEntry(i);
    for(int j=0;j<size;j++){
      if (IsMuonInAcceptance(pt1[j],pt1[j],eta1[j])&&IsMuonInAcceptance(pt2[j],pt2[j],eta2[j])){
        for(int m = 0; m < nMuPtBin; m++){
          if((pt1[j]>MuPtBin[m]&&pt1[j]>MuPtBin[m+1])&&(pt2[j]>MuPtBin[m]&&pt2[j]>MuPtBin[m+1])){
            
            hTrigPtAll[m]->Fill(mass[j]);
            if(isTriggered1[j]&&isTriggered2[j]) hTrigPtPass[m]->Fill(mass[j]);
		    else hTrigPtFail[m]->Fill(mass[j]);
          
          }//if pt in proper range
        }//loop over pt of muons
      }//muons in acceptance
    }//loop over candidates  
  }// loop over events
  
  TFile*foutput=new TFile("foutput.root","recreate");
  foutput->cd();
  
  for(int i = 0; i < nMuPtBin; i++){
    hTrigPtPass[i]->Write();
    hTrigPtFail[i]->Write();
    hTrigPtAll[i]->Write();
    hTrkPtPass[i]->Write();
    hTrkPtFail[i]->Write();
    hTrkPtAll[i]->Write();
    hMuonIdPtPass[i]->Write();
    hMuonIdPtFail[i]->Write();
    hMuonIdPtAll[i]->Write();

  }
  
  foutput->Close();
  delete foutput;
  
  TCanvas*canvasPass=new TCanvas("canvasPass","canvasPass",650,600);
  canvasPass->Divide(4,2);
  canvasPass->cd(1);
  hTrigPtPass[0]->Draw();
  canvasPass->cd(2);
  hTrigPtPass[1]->Draw();
  canvasPass->cd(3);
  hTrigPtPass[2]->Draw();
  canvasPass->cd(4);
  hTrigPtPass[3]->Draw();
  canvasPass->cd(5);
  hTrigPtPass[4]->Draw();
  canvasPass->cd(6);
  hTrigPtPass[5]->Draw();
  canvasPass->cd(7);
  hTrigPtPass[6]->Draw();
  canvasPass->SaveAs("canvasPass.eps");

  TCanvas*canvasFail=new TCanvas("canvasFail","canvasFail",650,600);
  canvasFail->Divide(4,2);
  canvasFail->cd(1);
  hTrigPtFail[0]->Draw();
  canvasFail->cd(2);
  hTrigPtFail[1]->Draw();
  canvasFail->cd(3);
  hTrigPtFail[2]->Draw();
  canvasFail->cd(4);
  hTrigPtFail[3]->Draw();
  canvasFail->cd(5);
  hTrigPtFail[4]->Draw();
  canvasFail->cd(6);
  hTrigPtFail[5]->Draw();
  canvasFail->cd(7);
  hTrigPtFail[6]->Draw();
  canvasFail->SaveAs("canvasFail.eps");

  TCanvas*canvasAll=new TCanvas("canvasAll","canvasAll",650,600);
  canvasAll->Divide(4,2);
  canvasAll->cd(1);
  hTrigPtAll[0]->Draw();
  canvasAll->cd(2);
  hTrigPtAll[1]->Draw();
  canvasAll->cd(3);
  hTrigPtAll[2]->Draw();
  canvasAll->cd(4);
  hTrigPtAll[3]->Draw();
  canvasAll->cd(5);
  hTrigPtAll[4]->Draw();
  canvasAll->cd(6);
  hTrigPtAll[5]->Draw();
  canvasAll->cd(7);
  hTrigPtAll[6]->Draw();
  canvasAll->SaveAs("canvasAll.eps");

}

bool IsMuonInAcceptance(Float_t pt,Float_t p,Float_t eta){
  bool isselected=false;
  isselected=(abs(eta)<1.3&&pt>3.3)||(abs(eta)>1.3&&abs(eta)<2.2&&p>2.9)||(abs(eta)>2.2&&abs(eta)<2.4&&pt>0.8);
  return isselected;
}
