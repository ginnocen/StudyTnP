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
  
  Bool_t IsMuonInAcceptance(Float_t,Float_t,Float_t);
  Bool_t IsTag(Bool_t, Int_t, Bool_t, Bool_t);

  //TString infname_mc="/data/bmeson/Data_Jpsi/jpsi.root"; 
  TString infname_mc="/data/ginnocen/TnPinputsMC/nt_BoostedMC_20140806_HIJINGemb_BuJpsiK_TuneZ2star_5TeV.root";
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
  Int_t isCalo1[NUM_BX];
  Int_t isCalo2[NUM_BX];
  

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
  ntuple->SetBranchAddress("isCalo1",isCalo1);
  ntuple->SetBranchAddress("isCalo2",isCalo2);

  Bool_t qualitycut1;
  Bool_t qualitycut2;
  Bool_t glb_cut1;
  Bool_t glb_cut2;
  Bool_t isTag1;
  Bool_t isTag2;
  Bool_t isacceptance1;
  Bool_t isacceptance2;

  Int_t entries = (Int_t)ntuple->GetEntries();
  
  for (int i=0; i<entries; i++){
    ntuple->GetEntry(i);

    for(int j=0;j<size;j++){
    
      qualitycut1 = false;
      qualitycut2 = false;
      glb_cut1 = false;
      glb_cut2 = false;
      isTag1=false;
      isTag2=false;
      isacceptance1=false;
      isacceptance2=false;


      if(id1[j]) qualitycut1 = true;
      if(id2[j]) qualitycut2 = true;
      if(isTrackerMuArbitrated1[j]&&isTMOneStationTight1[j]) glb_cut1 = true;
      if(isTrackerMuArbitrated2[j]&&isTMOneStationTight2[j]) glb_cut2 = true;
      isacceptance1=IsMuonInAcceptance(pt1[j],p1[j],eta1[j]);
      isacceptance2=IsMuonInAcceptance(pt2[j],p2[j],eta2[j]);
      
      isTag1=IsTag(isacceptance1,isTriggered1[j],qualitycut1,glb_cut1);
      isTag2=IsTag(isacceptance2,isTriggered2[j],qualitycut2,glb_cut2);

      if(isTag1||isTag2){
        if(isTag1){
          
          for(int m = 0; m < nMuPtBin; m++){
            if(pt2[j] > MuPtBin[m] && pt2[j] < MuPtBin[m+1]){
              //tracking efficiency
              if(outerTrackisNonnull2[j]){
                hTrkPtAll[m]->Fill(mass[j]);
                if(qualitycut2&&isTracker2[j]) hTrkPtPass[m]->Fill(mass[j]);
                else hTrkPtFail[m]->Fill(mass[j]);
              }//tracking efficinecy over

              //muon ID efficiency
              if(qualitycut2&&isCalo2[j]){
                hMuIdPtAll[m]->Fill(mass[j]);
                if(isTracker2[j]&&isacceptance2) hMuIdPtPass[m]->Fill(mass[j]);
                else hMuIdPtFail[m]->Fill(mass[j]);
              }//muon ID efficiency over

              //trigger efficiency
              if(qualitycut2&&glb_cut2&&isTracker2[j]&&isacceptance2){
                hTrigPtAll[m]->Fill(mass[j]);
                if(isTriggered2[j]) hTrigPtPass[m]->Fill(mass[j]);
                else hTrigPtFail[m]->Fill(mass[j]);
              }//trigger efficiency over
            }//if proper pt bin
          }//loop over pt bins
        }//if isTag1
        else{
          
          for(int m = 0; m < nMuPtBin; m++){
            if(pt1[j] > MuPtBin[m] && pt1[j] < MuPtBin[m+1]){
              //tracking efficiency
              if(outerTrackisNonnull1[j]){
                hTrkPtAll[m]->Fill(mass[j]);
                if(qualitycut1&&isTracker1[j]) hTrkPtPass[m]->Fill(mass[j]);
                else hTrkPtFail[m]->Fill(mass[j]);
              }//tracking efficinecy over

              //muon ID efficiency
              if(qualitycut1&&isCalo1[j]){
                hMuIdPtAll[m]->Fill(mass[j]);
                if(isTracker1[j]&&isacceptance1) hMuIdPtPass[m]->Fill(mass[j]);
                else hMuIdPtFail[m]->Fill(mass[j]);
              }//muon ID efficiency over

              //trigger efficiency
              if(qualitycut1&&glb_cut1&&isTracker1[j]&&isacceptance1){
                hTrigPtAll[m]->Fill(mass[j]);
                if(isTriggered1[j]) hTrigPtPass[m]->Fill(mass[j]);
                else hTrigPtFail[m]->Fill(mass[j]);
              }//trigger efficiency over
            }//if proper pt bin
          }//loop over pt bins
        }//if isTag2
      }//if isTag1||isTag2
      
      if(isTag1||isTag2){
        if(isTag1){
          
          for(int m = 0; m < nMuEtaBin; m++){
            if(eta2[j] > MuEtaBin[m] && eta2[j] < MuEtaBin[m+1]){
              //tracking efficiency
              if(outerTrackisNonnull2[j]){
                hTrkEtaAll[m]->Fill(mass[j]);
                if(qualitycut2&&isTracker2[j]) hTrkEtaPass[m]->Fill(mass[j]);
                else hTrkEtaFail[m]->Fill(mass[j]);
              }//tracking efficinecy over

              //muon ID efficiency
              if(qualitycut2&&isCalo2[j]){
                hMuIdEtaAll[m]->Fill(mass[j]);
                if(isTracker2[j]&&isacceptance2) hMuIdEtaPass[m]->Fill(mass[j]);
                else hMuIdEtaFail[m]->Fill(mass[j]);
              }//muon ID efficiency over

              //trigger efficiency
              if(qualitycut2&&glb_cut2&&isTracker2[j]&&isacceptance2){
                hTrigEtaAll[m]->Fill(mass[j]);
                if(isTriggered2[j]) hTrigEtaPass[m]->Fill(mass[j]);
                else hTrigEtaFail[m]->Fill(mass[j]);
              }//trigger efficiency over
            }//if proper eta bin
          }//loop over eta bins
        }//if isTag1
        else{
          
          for(int m = 0; m < nMuEtaBin; m++){
            if(eta1[j] > MuEtaBin[m] && eta1[j] < MuEtaBin[m+1]){
              //tracking efficiency
              if(outerTrackisNonnull1[j]){
                hTrkEtaAll[m]->Fill(mass[j]);
                if(qualitycut1&&isTracker1[j]) hTrkEtaPass[m]->Fill(mass[j]);
                else hTrkEtaFail[m]->Fill(mass[j]);
              }//tracking efficinecy over

              //muon ID efficiency
              if(qualitycut1&&isCalo1[j]){
                hMuIdEtaAll[m]->Fill(mass[j]);
                if(isTracker1[j]&&isacceptance1) hMuIdEtaPass[m]->Fill(mass[j]);
                else hMuIdEtaFail[m]->Fill(mass[j]);
              }//muon ID efficiency over

              //trigger efficiency
              if(qualitycut1&&glb_cut1&&isTracker1[j]&&isacceptance1){
                hTrigEtaAll[m]->Fill(mass[j]);
                if(isTriggered1[j]) hTrigEtaPass[m]->Fill(mass[j]);
                else hTrigEtaFail[m]->Fill(mass[j]);
              }//trigger efficiency over
            }//if proper eta bin
          }//loop over eta bins
        }//if isTag2
      }//if isTag1||isTag2

      
      
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
  }

  for(int i = 0; i < nMuEtaBin; i++){

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

Bool_t IsMuonInAcceptance(Float_t pt,Float_t p,Float_t eta){
  Bool_t isselected=false;
  isselected=(abs(eta)<1.3&&pt>3.3)||(abs(eta)>1.3&&abs(eta)<2.2&&p>2.9)||(abs(eta)>2.2&&abs(eta)<2.4&&pt>0.8);
  return isselected;
}
Bool_t IsTag(Bool_t isacceptance, Int_t istrigger, Bool_t qualitycut, Bool_t glb_cut)
{
  Bool_t isTag=false;
  isTag=(isacceptance&& istrigger && qualitycut&& glb_cut);
  return isTag;
}
