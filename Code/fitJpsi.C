#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "utilities.h"
#include "TH1D.h"
#define NUM_BX 9000


void fitJpsi(){

  const int nMuPtBin = 7;
  const int nMuEtaBin = 5;
  const double MuPtBin[nMuPtBin+1] = {0.0,1.5,3.0,4.5,6.0,9.0,20.0,30.0};
  const double MuEtaBin[nMuEtaBin+1] = {-2.4,-1.5,-0.5,0.5,1.5,2.4};
  
  int nBin = 100;
  double mumulow=2.6;
  double mumuhigh=3.5;
  double mumuTrklow=2.0;
  double mumuTrkhigh=5.0;
  
  bool IsMuonInAcceptance(Float_t,Float_t,Float_t);
  bool IsTag(Float_t, Float_t, Float_t, Int_t, bool, bool);

  //TString infname_mc="/data/bmeson/MC_Jpsi/jpsi-Kp.root"; 
  TString infname_mc="/data/ginnocen/TnPinputsMC/nt_BoostedMC_20140707_BdJpsiKstar_pPb_TnP.root";
  TFile *inf_mc = new TFile(infname_mc.Data());
  
  TTree *ntuple = (TTree*) inf_mc->Get("ntJpsi");
  TTree *nt_mcgen = (TTree*)inf_mc->Get("ntGen");
  ntuple->AddFriend(nt_mcgen);
  
  //Pt probes
  //trigger efficiency
  TH1D* hTrigPtPass[nMuPtBin];
  TH1D* hTrigPtFail[nMuPtBin];
  TH1D* hTrigPtAll[nMuPtBin];
  
  //tracking efficiency
  TH1D* hTrkPtPass[nMuPtBin];
  TH1D* hTrkPtFail[nMuPtBin];
  TH1D* hTrkPtAll[nMuPtBin];
  
  //Muon ID efficiency
  TH1D* hMuIdPtPass[nMuPtBin];
  TH1D* hMuIdPtFail[nMuPtBin];
  TH1D* hMuIdPtAll[nMuPtBin];

  //Eta probes
  //trigger efficiency
  TH1D* hTrigEtaPass[nMuEtaBin];
  TH1D* hTrigEtaFail[nMuEtaBin];
  TH1D* hTrigEtaAll[nMuEtaBin];
  
  //tracking efficiency
  TH1D* hTrkEtaPass[nMuEtaBin];
  TH1D* hTrkEtaFail[nMuEtaBin];
  TH1D* hTrkEtaAll[nMuEtaBin];
  
  //Muon ID efficinecy
  TH1D* hMuIdEtaPass[nMuEtaBin];
  TH1D* hMuIdEtaFail[nMuEtaBin];
  TH1D* hMuIdEtaAll[nMuEtaBin];

  for(int i = 0; i < nMuPtBin; i++){
    hTrigPtPass[i] = new TH1D(Form("hTrigPtPass%d",i),"",nBin,mumulow,mumuhigh);
    hTrigPtFail[i] = new TH1D(Form("hTrigPtFail%d",i),"",nBin,mumulow,mumuhigh);
    hTrigPtAll[i] = new TH1D(Form("hTrigPtAll%d",i),"",nBin,mumulow,mumuhigh);
    hTrkPtPass[i] = new TH1D(Form("hTrkPtPass%d",i),"",nBin,mumuTrklow,mumuTrkhigh);
    hTrkPtFail[i] = new TH1D(Form("hTrkPtFail%d",i),"",nBin,mumuTrklow,mumuTrkhigh);
    hTrkPtAll[i] = new TH1D(Form("hTrkPtAll%d",i),"",nBin,mumuTrklow,mumuTrkhigh);
    hMuIdPtPass[i] = new TH1D(Form("hMuIdPtPass%d",i),"",nBin,mumulow,mumuhigh);
    hMuIdPtFail[i] = new TH1D(Form("hMuIdPtFail%d",i),"",nBin,mumulow,mumuhigh);
    hMuIdPtAll[i] = new TH1D(Form("hMuIdPtAll%d",i),"",nBin,mumulow,mumuhigh);
  }
  
  for(int i = 0; i < nMuEtaBin; i++){
    hTrigEtaPass[i] = new TH1D(Form("hTrigEtaPass%d",i),"",nBin,mumulow,mumuhigh);
    hTrigEtaFail[i] = new TH1D(Form("hTrigEtaFail%d",i),"",nBin,mumulow,mumuhigh);
    hTrigEtaAll[i] = new TH1D(Form("hTrigEtaAll%d",i),"",nBin,mumulow,mumuhigh);
    hTrkEtaPass[i] = new TH1D(Form("hTrkEtaPass%d",i),"",nBin,mumuTrklow,mumuTrkhigh);
    hTrkEtaFail[i] = new TH1D(Form("hTrkEtaFail%d",i),"",nBin,mumuTrklow,mumuTrkhigh);
    hTrkEtaAll[i] = new TH1D(Form("hTrkEtaAll%d",i),"",nBin,mumuTrklow,mumuTrkhigh);
    hMuIdEtaPass[i] = new TH1D(Form("hMuIdEtaPass%d",i),"",nBin,mumulow,mumuhigh);
    hMuIdEtaFail[i] = new TH1D(Form("hMuIdEtaFail%d",i),"",nBin,mumulow,mumuhigh);
    hMuIdEtaAll[i] = new TH1D(Form("hMuIdEtaAll%d",i),"",nBin,mumulow,mumuhigh);
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
  Float_t p1[NUM_BX];
  Float_t p2[NUM_BX];
  Float_t phi1[NUM_BX];
  Float_t phi2[NUM_BX];
  int id1[NUM_BX];
  int id2[NUM_BX];
  Bool_t outerTrackisNonnull1[NUM_BX];
  Bool_t outerTrackisNonnull2[NUM_BX];
  int gen[NUM_BX];
  Bool_t isTrackerMuArbitrated1[NUM_BX];
  Bool_t isTrackerMuArbitrated2[NUM_BX];
  Bool_t isTMOneStationTight1[NUM_BX];
  Bool_t isTMOneStationTight2[NUM_BX];
  Bool_t calomuon1[NUM_BX];
  Bool_t calomuon2[NUM_BX];

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
  ntuple->SetBranchAddress("p1",p1);
  ntuple->SetBranchAddress("p2",p2);
  ntuple->SetBranchAddress("phi1",phi1);
  ntuple->SetBranchAddress("phi2",phi2);
  ntuple->SetBranchAddress("id1",id1);
  ntuple->SetBranchAddress("id2",id2);
  ntuple->SetBranchAddress("outerTrackisNonnull1",outerTrackisNonnull1);
  ntuple->SetBranchAddress("outerTrackisNonnull2",outerTrackisNonnull2);
  ntuple->SetBranchAddress("isTrackerMuArbitrated1",isTrackerMuArbitrated1);
  ntuple->SetBranchAddress("isTrackerMuArbitrated2",isTrackerMuArbitrated2);
  ntuple->SetBranchAddress("isTMOneStationTight1",isTMOneStationTight1);
  ntuple->SetBranchAddress("isTMOneStationTight2",isTMOneStationTight2);
  ntuple->SetBranchAddress("gen",gen);
  ntuple->SetBranchAddress("isTriggered1",isTriggered1);
  ntuple->SetBranchAddress("isTriggered2",isTriggered2);
  ntuple->SetBranchAddress("calomuon1",calomuon1);
  ntuple->SetBranchAddress("calomuon2",calomuon2);
  
  Int_t entries = (Int_t)ntuple->GetEntries();
  
  for (int i=0; i<entries; i++)
  {
    ntuple->GetEntry(i);

    for(int j=0;j<size;j++)
    {
      Bool_t track_cut1 = false;
      Bool_t track_cut2 = false;
      Bool_t glb_cut1 = false;
      Bool_t glb_cut2 = false;

      if(*id1==1) track_cut1 = true;
      if(*id2==1) track_cut2 = true;
      if(isTrackerMuArbitrated1&&isTMOneStationTight1) glb_cut1 = true;
      if(isTrackerMuArbitrated2&&isTMOneStationTight2) glb_cut2 = true;

      if(IsTag(pt1[j],p1[j],eta1[j],isTriggered1[j],track_cut1,glb_cut1) || IsTag(pt2[j],p2[j],eta2[j],isTriggered2[j],track_cut2,glb_cut2))
      {
        if(IsTag(pt1[j],p1[j],eta1[j],isTriggered1[j],track_cut1,glb_cut1))
        {
          for(int m = 0; m < nMuPtBin; m++)
          {
            if(pt2[j] > MuPtBin[m] && pt2[j] < MuPtBin[m+1])
            {
              //tracking efficiency
              if(outerTrackisNonnull2[j])
              {
                hTrkPtAll[m]->Fill(mass[j]);
                if(track_cut2&&isTracker2[j]) hTrkPtPass[m]->Fill(mass[j]);
                else hTrkPtFail[m]->Fill(mass[j]);
              }//tracking efficinecy over

              //muon ID efficiency
              if(calomuon2[j]&&track_cut2)
              {
                hMuIdPtAll[m]->Fill(mass[j]);
                if(track_cut2&&isTracker2[j]&&IsMuonInAcceptance(pt2[j],p2[j],eta2[j])) hMuIdPtPass[m]->Fill(mass[j]);
                else hMuIdPtFail[m]->Fill(mass[j]);
              }//muon ID efficiency over

              //trigger efficiency
              if(track_cut2&&glb_cut2&&isTracker2[j]&&IsMuonInAcceptance(pt2[j],p2[j],eta2[j]))
              {
                hTrigPtAll[m]->Fill(mass[j]);
                if(isTriggered2[j]) hTrigPtPass[m]->Fill(mass[j]);
                else hTrigPtFail[m]->Fill(mass[j]);
              }//trigger efficiency over
            }
          }//loop over pt

          for(int m = 0; m < nMuEtaBin; m++)
          {
            if(eta2[j] > MuEtaBin[m] && eta2[j] < MuEtaBin[m+1])
            {
              //tracking efficiency
              if(outerTrackisNonnull1[j])
              {
                hTrkEtaAll[m]->Fill(mass[j]);
                if(track_cut2&&isTracker2[j]) hTrkEtaPass[m]->Fill(mass[j]);
                else hTrkEtaFail[m]->Fill(mass[j]);
              }//tracking efficinecy over

              //muon ID efficiency
              if(calomuon2[j]&&track_cut2)
              {
                hMuIdEtaAll[m]->Fill(mass[j]);
                if(track_cut2&&isTracker2[j]&&IsMuonInAcceptance(pt2[j],p2[j],eta2[j])) hMuIdEtaPass[m]->Fill(mass[j]);
                else hMuIdEtaFail[m]->Fill(mass[j]);
              }//muon ID efficiency over

              //trigger efficiency
              if(track_cut2&&glb_cut2&&isTracker2[j]&&IsMuonInAcceptance(pt2[j],p2[j],eta2[j]))
              {
                hTrigEtaAll[m]->Fill(mass[j]);
                if(isTriggered2[j]) hTrigEtaPass[m]->Fill(mass[j]);
                else hTrigEtaFail[m]->Fill(mass[j]);
              }//trigger efficiency over
            }
          }//loop over eta
        }//loop over tag is 1
        else
        {
          for(int m = 0; m < nMuPtBin; m++)
          {
            if(pt1[j] > MuPtBin[m] && pt1[j] < MuPtBin[m+1])
            {
              //tracking efficiency
              if(outerTrackisNonnull1[j])
              {
                hTrkPtAll[m]->Fill(mass[j]);
                if(track_cut1&&isTracker1[j]) hTrkPtPass[m]->Fill(mass[j]);
                else hTrkPtFail[m]->Fill(mass[j]);
              }//tracking efficinecy over

              //muon ID efficiency
              if(calomuon1[j]&&track_cut1)
              {
                hMuIdPtAll[m]->Fill(mass[j]);
                if(track_cut1&&isTracker1[j]&&IsMuonInAcceptance(pt1[j],p1[j],eta1[j])) hMuIdPtPass[m]->Fill(mass[j]);
                else hMuIdPtFail[m]->Fill(mass[j]);
              }//muon ID efficiency over

              //trigger efficiency
              if(track_cut1&&glb_cut1&&isTracker1[j]&&IsMuonInAcceptance(pt1[j],p1[j],eta1[j]))
              {
                hTrigPtAll[m]->Fill(mass[j]);
                if(isTriggered1[j]) hTrigPtPass[m]->Fill(mass[j]);
                else hTrigPtFail[m]->Fill(mass[j]);
              }//trigger efficiency over
            }
          }//loop over pt

          for(int m = 0; m < nMuEtaBin; m++)
          {
            if(eta1[j] > MuEtaBin[m] && eta1[j] < MuEtaBin[m+1])
            {
              //tracking efficiency
              if(outerTrackisNonnull1[j])
              {
                hTrkEtaAll[m]->Fill(mass[j]);
                if(track_cut1&&isTracker1[j]) hTrkEtaPass[m]->Fill(mass[j]);
                else hTrkEtaFail[m]->Fill(mass[j]);
              }//tracking efficinecy over

              //muon ID efficiency
              if(calomuon1[j]&&track_cut1)
              {
                hMuIdEtaAll[m]->Fill(mass[j]);
                if(track_cut1&&isTracker1[j]&&IsMuonInAcceptance(pt1[j],p1[j],eta1[j])) hMuIdEtaPass[m]->Fill(mass[j]);
                else hMuIdEtaFail[m]->Fill(mass[j]);
              }//muon ID efficiency over

              //trigger efficiency
              if(track_cut1&&glb_cut1&&isTracker1[j]&&IsMuonInAcceptance(pt1[j],p1[j],eta1[j]))
              {
                hTrigEtaAll[m]->Fill(mass[j]);
                if(isTriggered1[j]) hTrigEtaPass[m]->Fill(mass[j]);
                else hTrigEtaFail[m]->Fill(mass[j]);
              }//trigger efficiency over
            }
          }//loop over eta
        }//loop over tag is 2
      }
    }//loop over candidates  
  }// loop over events
  
  TFile*foutput=new TFile("Results/foutput.root","recreate");
  foutput->cd();
  
  for(int i = 0; i < nMuPtBin; i++){
    hTrigPtPass[i]->Write();
    hTrigPtFail[i]->Write();
    hTrigPtAll[i]->Write();
    hTrkPtPass[i]->Write();
    hTrkPtFail[i]->Write();
    hTrkPtAll[i]->Write();
    hMuIdPtPass[i]->Write();
    hMuIdPtFail[i]->Write();
    hMuIdPtAll[i]->Write();

    hTrigEtaPass[i]->Write();
    hTrigEtaFail[i]->Write();
    hTrigEtaAll[i]->Write();
    hTrkEtaPass[i]->Write();
    hTrkEtaFail[i]->Write();
    hTrkEtaAll[i]->Write();
    hMuIdEtaPass[i]->Write();
    hMuIdEtaFail[i]->Write();
    hMuIdEtaAll[i]->Write();
  }
  
  foutput->Close();
  delete foutput;

}

bool IsMuonInAcceptance(Float_t pt,Float_t p,Float_t eta){
  bool isselected=false;
  isselected=(abs(eta)<1.3&&pt>3.3)||(abs(eta)>1.3&&abs(eta)<2.2&&p>2.9)||(abs(eta)>2.2&&abs(eta)<2.4&&pt>0.8);
  return isselected;
}
bool IsTag(Float_t pt, Float_t p, Float_t eta, Int_t trigger, bool track_cut, bool glb_cut)
{
  bool isTag=false;
  isTag=(IsMuonInAcceptance(pt, p, eta) && trigger && track_cut && glb_cut);
  return isTag;
}
