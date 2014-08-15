#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TH1D.h"

#include "format.h"
#include "utilities.h"

#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooDataSet.h>
#include <RooPlot.h>
#include <RooGaussian.h>
#include <RooExponential.h>
#include <RooChevychev.h>
#include <RooAddPdf.h>
#include <RooDecay.h>
#include <RooGaussModel.h>
#include <RooProdPdf.h>
#include <RooAddModel.h>
#include <RooCBShape.h>
#include <RooCategory.h>
#include <RooSimultaneous.h>

#define NUM_BX 9000


void FitTnP(){

  using namespace std;
  using namespace RooFit;

  const int nMuPtBin = 7;
  const int nMuEtaBin = 5;
  double mumulow=2.6;
  double mumuhigh=3.5.;
  double mumutrklow=2.;
  double mumutrkhigh=5.;
  
  const double MuPtBin[nMuPtBin+1] = {0.0, 1.5, 3.0, 4.5, 6.0, 9.0, 20.0, 30.0};
  const double MuEtaBin[nMuEtaBin+1] = {-2.4, -1.5, -0.5, 0.5, 1.5, 2.4};

  RooCategory sampleTrig("sampleTrig", "sampleTrig");
  sampleTrig.defineType("Pass");
  sampleTrig.defineType("Fail");
  sampleTrig.defineType("All");
  RooCategory sampleMuId("sampleMuId", "sampleMuId");
  sampleMuId.defineType("Pass");
  sampleMuId.defineType("Fail");
  sampleMuId.defineType("All");
  RooCategory sampleTrk("sampleTrk", "sampleTrk");
  sampleTrk.defineType("Pass");
  sampleTrk.defineType("Fail");
  sampleTrk.defineType("All");

  //pt histogram
  TH1D* hTrigPtPass[nMuPtBin];
  TH1D* hTrigPtFail[nMuPtBin];
  TH1D* hTrigPtAll[nMuPtBin];
  
  TH1D* hTrkPtPass[nMuPtBin];
  TH1D* hTrkPtFail[nMuPtBin];
  TH1D* hTrkPtAll[nMuPtBin];
  
  TH1D* hMuIdPtPass[nMuPtBin];
  TH1D* hMuIdPtFail[nMuPtBin];
  TH1D* hMuIdPtAll[nMuPtBin];
  
  //eta histogram
  TH1D* hTrigEtaPass[nMuPtBin];
  TH1D* hTrigEtaFail[nMuPtBin];
  TH1D* hTrigEtaAll[nMuPtBin];
  
  TH1D* hTrkEtaPass[nMuPtBin];
  TH1D* hTrkEtaFail[nMuPtBin];
  TH1D* hTrkEtaAll[nMuPtBin];
  
  TH1D* hMuIdEtaPass[nMuPtBin];
  TH1D* hMuIdEtaFail[nMuPtBin];
  TH1D* hMuIdEtaAll[nMuPtBin];

  TFile*finput=new TFile("../Code/Results/foutput.root","read");
  finput->cd();

  //get pt histogram
  for(int i = 0; i < nMuPtBin; i++)
  {
    hTrigPtPass[i]=(TH1D*)finput->Get(Form("hTrigPtPass%d",i));
    hTrigPtFail[i]=(TH1D*)finput->Get(Form("hTrigPtFail%d",i));
    hTrigPtAll[i]=(TH1D*)finput->Get(Form("hTrigPtAll%d",i));
    hTrkPtPass[i]=(TH1D*)finput->Get(Form("hTrkPtPass%d",i));
    hTrkPtFail[i]=(TH1D*)finput->Get(Form("hTrkPtFail%d",i));
    hTrkPtAll[i]=(TH1D*)finput->Get(Form("hTrkPtAll%d",i));
    hMuIdPtPass[i]=(TH1D*)finput->Get(Form("hMuIdPtPass%d",i));
    hMuIdPtFail[i]=(TH1D*)finput->Get(Form("hMuIdPtFail%d",i));
    hMuIdPtAll[i]=(TH1D*)finput->Get(Form("hMuIdPtAll%d",i));
  }

  //get eta histogram
  for(int i = 0; i < nMuEtaBin; i++)
  {
    hTrigPass[i]=(TH1D*)finput->Get(Form("hTrigEtaPass%d",i));
    hTrigFail[i]=(TH1D*)finput->Get(Form("hTrigEtaFail%d",i));
    hTrigAll[i]=(TH1D*)finput->Get(Form("hTrigEtaAll%d",i));
    hTrkPass[i]=(TH1D*)finput->Get(Form("hTrkEtaPass%d",i));
    hTrkFail[i]=(TH1D*)finput->Get(Form("hTrkEtaFail%d",i));
    hTrkAll[i]=(TH1D*)finput->Get(Form("hTrkEtaAll%d",i));
    hMuIdPass[i]=(TH1D*)finput->Get(Form("hMuIdEtaPass%d",i));
    hMuIdFail[i]=(TH1D*)finput->Get(Form("hMuIdEtaFail%d",i));
    hMuIdAll[i]=(TH1D*)finput->Get(Form("hMuIdEtaAll%d",i));
  }

  //make CB function
  RooRealVar mean_CB("mean_CB", "mean_CB", 3.1, 3.0, 3.2);
  RooRealVar sigma_CB("sigma_CB", "sigma_CB", 0.05);
  isgma_CB.setConstant(kFALSE);
  RooRealVar alpha_CB("alpha_CB", "alpha_CB", 2.0, 1.0, 5.0);
  RooRealVar n_CB("n_CB", "n_CB", 1, 0.5, 100.0);

  RooCBShape signal_CB_Pass("signal_CB_Pass", "signal_CB_Pass", mass, mean_CB, sigma_CB, alpha_CB, n_CB);
  RooCBShape signal_CB_Fail("signal_CB_Fail", "signal_CB_Fail", mass, mean_CB, sigma_CB, alpha_CB, n_CB);
  RooCBShape signal_CB_All("signal_CB_All", "signal_CB_sAll", mass, mean_CB, sigma_CB, alpha_CB, n_CB);

  //make 2gauss function
  RooRealVar mean_gaus("mean_gaus", "mean_gaus", 3.1, 3.0, 3.2);
  RooRealVar sigma1_gaus("sigma1_gaus", "sigma1_gaus", 0.15, 0.05, 0.25);
  RooRealVar sigma2_gaus("sigma2_gaus", "sigma2_gaus", 0.02, 0.01, 0.1);
  RooGaussian gauss1("gauss1", "gauss1", mass, mean_gaus, sigma1_gaus);
  RooGaussian gauss2("gauss2", "gauss2", mass, mean_gaus, sigma2_gaus);
  RooRealVar mfrac_gaus("mfrac_gaus", "mfrac_gaus", 0.2, 0.0, 1.0);

  RooAddPdf signal_gauss_Pass("signal_gauss_Pass", "signal_gauss_Pass", RooArgList(gauss1, gauss2), RooArgList(mfrac_gaus));
  RooAddPdf signal_gauss_Fail("signal_gauss_Fail", "signal_gauss_Fail", RooArgList(gauss1, gauss2), RooArgList(mfrac_gaus));
  RooAddPdf signal_gauss_All("signal_gauss_All", "signal_gauss_All", RooArgList(gauss1, gauss2), RooArgList(mfrac_gaus));

  //make Chebyschev function
  RooRealVar cheb_p1Pass("cheb_paPass", "cheb_p1Pass", 0., -1., +1.);
  RooRealVar cheb_p2Pass("cheb_p2Pass", "cheb_p2Pass", 0., -1., +1.);
  RooChebychev background_cheb_Pass("background_cheb_Pass", "background_cheb_Pass", mass, RooArgList(cheb_p1Pass, cheb_p2Pass));
  RooRealVar cheb_p1Fail("cheb_p1Fail", "cheb_p1Fail", 0., -1., +1.);
  RooRealVar cheb_p2Fail("cheb_p2Fail", "cheb_p2Fail", 0., -1., +1.);
  RooChebychev background_cheb_Fail("background_cheb_Fail", "background_cheb_Fail", mass, RooArgList(cheb_p1Fail, cheb_p2Fail));
  RooRealVar cheb_p1All("cheb_p1All", "cheb_p1All", 0., -1., +1.);
  RooRealVar cheb_p2Fail("cheb_p2All", "cheb_p2All", 0., -1., +1.);
  RooChebychev background_cheb_All("background_cheb_All", "background_cheb_All", mass, RooArgList(cheb_p1All, cheb_p2All));

  //make exponential function
  RooRealVar lp("lp", "lp", 0.0, -5., 5.);
  RooRealVar lf("lf", "lf", 0.0, -5., 5.);
  RooRealVar la("la", "la", 0.0, -5., 5.);
  RooExponential background_exp_Pass("background_exp_Pass", "background_exp_Pass", mass, lp);
  RooExponential background_exp_Fail("background_exp_Fail", "background_exp_Fail", mass, lf);
  RooExponential background_exp_All("background_exp_All", "background_exp_All", mass, la);

  //define PDF
  RooAddPdf model_CB_cheb_Pass("model_CB_cheb_Pass", "model_CB_cheb_Pass", RooArgList(signal_CB_Pass, background_cheb_Pass), RooArgList(nsigPass, nbkgPass));
  RooAddPdf model_CB_cheb_Fail("model_CB_cheb_Fail", "model_CB_cheb_Fail", RooArgList(signal_CB_Fail, background_cheb_Fail), RooArgList(nsigFail, nbkgFail));
  RooAddPdf model_CB_cheb_All("model_CB_cheb_All", "model_CB_cheb_All", RooArgList(signal_CB_All, background_cheb_All), RooArgList(nsigAll, nbkgAll));

  RooAddPdf model_CB_exp_Pass("model_CB_exp_Pass", "model_CB_exp_Pass", RooArgList(signal_CB_Pass, background_exp_Pass), RooArgList(nsigPass, nbkgPass));
  RooAddPdf model_CB_exp_Fail("model_CB_exp_Fail", "model_CB_exp_Fail", RooArgList(signal_CB_Fail, background_exp_Fail), RooArgList(nsigFail, nbkgFail));
  RooAddPdf model_CB_exp_All("model_CB_exp_All", "model_CB_exp_All", RooArgList(signal_CB_All, background_exp_All), RooArgList(nsigAll, nbkgAll));

  RooAddPdf model_gaus_cheb_Pass("model_gaus_cheb_Pass", "model_gaus_cheb_Pass", RooArgList(signal_gaus_Pass, background_cheb_Pass), RooArgList(nsigPass, nbkgPass));
  RooAddPdf model_gaus_cheb_Fail("model_gaus_cheb_Fail", "model_gaus_cheb_Fail", RooArgList(signal_gaus_Fail, background_cheb_Fail), RooArgList(nsigFail, nbkgFail));
  RooAddPdf model_gaus_cheb_All("model_gaus_cheb_All", "model_gaus_cheb_All", RooArgList(signal_gaus_All, background_cheb_All), RooArgList(nsigAll, nbkgAll));

  //fitting pt
  for(int i = 0; i < nMuPtBin; i++)
  {
    //set axis variable
    RooRealVar x("x", "mass", mumulow, mumuhigh, "GeV");
    RooRealVar xtrk("x", "mass", mumutrklow, mumutrkhigh, "GeV");

    //set histogram
    RooDataHist dataTrigPass("dataTrigPass", "dataTrigPass", x, hTrigPtPass[i]);
    RooDataHist dataTrigFail("dataTrigFail", "dataTrigFail", x, hTrigPtFail[i]);
    RooDataHist dataTrigAll("dataTrigAll", "dataTrigAll", x, hTrigPtAll[i]);

    RooDataHist dataMuIdPass("dataMuIdPass", "dataMuIdPass", x, hPtMuIdPass[i]);
    RooDataHist dataMuIdFail("dataMuIdFail", "dataMuIdFail", x, hPtMuIdFail[i]);
    RooDataHist dataMuIdAll("dataMuIdAll", "dataMuIdAll", x, hPtMuIdAll[i]);

    RooDataHist dataTrkPass("dataTrkPass", "dataTrkPass", xtrk, hTrkPtPass[i]);
    RooDataHist dataTrkFail("dataTrkFail", "dataTrkFail", xtrk, hTrkPtFail[i]);
    RooDataHist dataTrkAll("dataTrkAll", "dataTrkAll", xtrk, hTrkPtAll[i]);

    //fit
    RooSimultaneous simPdf("simPdf", "simPdf", sample);
  }
  
  //fitting eta
  for(int i = 0; i < nMuEtaBin; i++)
  {
    //set axis variable
    RooRealVar x("x", "mass", mumulow, mumuhigh, "GeV");
    RooRealVar xtrk("x", "mass", mumutrklow, mumutrkhigh, "GeV");

    //set histogram
    RooDataHist dataTrigPass("dataTrigPass", "dataTrigPass", x, hTrigEtaPass[i]);
    RooDataHist dataTrigFail("dataTrigFail", "dataTrigFail", x, hTrigEtaFail[i]);
    RooDataHist dataTrigAll("dataTrigAll", "dataTrigAll", x, hTrigEtaAll[i]);

    RooDataHist dataMuIdPass("dataMuIdPass", "dataMuIdPass", x, hEtaMuIdPass[i]);
    RooDataHist dataMuIdFail("dataMuIdFail", "dataMuIdFail", x, hEtaMuIdFail[i]);
    RooDataHist dataMuIdAll("dataMuIdAll", "dataMuIdAll", x, hEtaMuIdAll[i]);

    RooDataHist dataTrkPass("dataTrkPass", "dataTrkPass", x, hTrkEtaPass[i]);
    RooDataHist dataTrkFail("dataTrkFail", "dataTrkFail", x, hTrkEtaFail[i]);
    RooDataHist dataTrkAll("dataTrkAll", "dataTrkAll", x, hTrkEtaAll[i]);
  }




  RooSimultaneous simPdf("simPdf", "simPdf", sample);


  simPdf.fitTo(combData, Extended(true));

  RooPlot *frame1 = mass.frame();
  combData.plotOn(frame1, Binning(100), Cut("sample==sample::Pass"));
  simPdf.plotOn(frame1, Slice(sample, "Pass"), ProjWData(sample, combData));
  simPdf.PlotOn(frame1, Slice(sample, "Pass"), Components("signal_CB_Pass"), ProjWData(sample, combData), LineStyle(kDashed), LineColor(kGreen));
  simPdf.PlotOn(frame1, Slice(sample, "Pass"), Components("background_cheb_Pass"), ProjWData(sample, combData), LineStyle(kDashed), LineColor(kGreen));

  RooPlot *frame2 = mass.frame();
  combData.plotOn(frame2, Binning(100), Cut("sample==sample::Fail"));
  simPdf.plotOn(frame2, Slice(sample, "Fail"), ProjWData(sample, combData));
  simPdf.plotOn(frame2, Slice(sample, "Fail"), Components("signal_CB_Fail"), ProjWData(sample, combData), LineStyle(kDashed), LineColor(kGreen));

  TF1 *signal = new TF1("Gaus(x,[0],[1])","gaus");  
  signal->SetParameter(0,3.1);
  hTrigPtPass[0]->Fit("signal", "L q", "",2.8,3.2);

  
  TCanvas*ctest=new TCanvas("ctest","ctest",650,600);
  ctest->cd();
  hTrigPtPass[0]->Draw();
  signal->Draw("same");
  ctest->SaveAs("ctest.eps");

}


