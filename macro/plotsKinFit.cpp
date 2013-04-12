#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH2D.h"
#include "TLegend.h"
#include <iostream>

void plotKineFit() {

  TFile file("finalizedTrees_Radion/Radion_Radion_M-300_madgraph_default_CSV.root");
  TH2D *h2_mggjj_vs_mjj_1btag        = (TH2D*)file.Get("mggjj_vs_mjj_1btag");
  TH2D *h2_mggjj_vs_mjj_kinfit_1btag = (TH2D*)file.Get("mggjj_vs_mjj_kinfit_1btag"); 
  TH1F *h1_mggjj_kinfit_1btag        = (TH1F*)file.Get("mggjj_kinfit_1btag");
  TH1F *h1_mggjj_nokinfit_1btag      = (TH1F*)file.Get("mggjj_nokinfit_1btag");

  TH2D *h2_mggjj_vs_mjj_2btag        = (TH2D*)file.Get("mggjj_vs_mjj_2btag");
  TH2D *h2_mggjj_vs_mjj_kinfit_2btag = (TH2D*)file.Get("mggjj_vs_mjj_kinfit_2btag"); 
  TH1F *h1_mggjj_kinfit_2btag        = (TH1F*)file.Get("mggjj_kinfit_2btag");
  TH1F *h1_mggjj_nokinfit_2btag      = (TH1F*)file.Get("mggjj_nokinfit_2btag");

  h2_mggjj_vs_mjj_1btag->SetTitle("1b-jet category, no kine fit");
  h2_mggjj_vs_mjj_1btag->GetXaxis()->SetTitle("mjj [GeV]");
  h2_mggjj_vs_mjj_1btag->GetYaxis()->SetTitle("mggjj [GeV]");

  h2_mggjj_vs_mjj_kinfit_1btag->SetTitle("1b-jet category, with kine fit");
  h2_mggjj_vs_mjj_kinfit_1btag->GetXaxis()->SetTitle("mjj [GeV]");
  h2_mggjj_vs_mjj_kinfit_1btag->GetYaxis()->SetTitle("mggjj [GeV]");

  h1_mggjj_kinfit_1btag->SetTitle("1b-jet category");
  h1_mggjj_kinfit_1btag->GetXaxis()->SetTitle("mggjj [GeV]");
  h1_mggjj_kinfit_1btag->SetLineColor(2);

  h1_mggjj_nokinfit_1btag->SetTitle("1b-jet category");
  h1_mggjj_nokinfit_1btag->GetXaxis()->SetTitle("mggjj [GeV]");
  h1_mggjj_nokinfit_1btag->SetLineColor(4);

  h2_mggjj_vs_mjj_2btag->SetTitle("2b-jet category, no kine fit");
  h2_mggjj_vs_mjj_2btag->GetXaxis()->SetTitle("mjj [GeV]");
  h2_mggjj_vs_mjj_2btag->GetYaxis()->SetTitle("mggjj [GeV]");

  h2_mggjj_vs_mjj_kinfit_2btag->SetTitle("2b-jet category, with kine fit");
  h2_mggjj_vs_mjj_kinfit_2btag->GetXaxis()->SetTitle("mjj [GeV]");
  h2_mggjj_vs_mjj_kinfit_2btag->GetYaxis()->SetTitle("mggjj [GeV]");

  h1_mggjj_kinfit_2btag->SetTitle("2b-jet category");
  h1_mggjj_kinfit_2btag->GetXaxis()->SetTitle("mggjj [GeV]");
  h1_mggjj_kinfit_2btag->SetLineColor(2);

  h1_mggjj_nokinfit_2btag->SetTitle("2b-jet category");
  h1_mggjj_nokinfit_2btag->GetXaxis()->SetTitle("mggjj [GeV]");
  h1_mggjj_nokinfit_2btag->SetLineColor(4);

  gStyle->SetOptStat(0);

  TLegend *leg;
  leg = new TLegend(0.6,0.6,0.85,0.85);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.05);
  leg->SetFillColor(0);
  leg->AddEntry(h1_mggjj_kinfit_2btag,   "with kin fit", "l");
  leg->AddEntry(h1_mggjj_nokinfit_2btag, "no kin fit",   "l");

  TCanvas myC1("c1", "c1", 1);
  h2_mggjj_vs_mjj_kinfit_1btag->Draw("colz");
  myC1.SaveAs("corrK_1.png");

  TCanvas myC2("c2", "c2", 1);
  h2_mggjj_vs_mjj_1btag->Draw("colz");
  myC2.SaveAs("corrNoK_1.png");

  TCanvas myC11("c11", "c11", 1);
  h2_mggjj_vs_mjj_kinfit_2btag->Draw("colz");
  myC11.SaveAs("corrK_2.png");

  TCanvas myC22("c22", "c22", 1);
  h2_mggjj_vs_mjj_2btag->Draw("colz");
  myC22.SaveAs("corrNoK_2.png");

  TCanvas myC31("c31", "c31", 1);
  h1_mggjj_kinfit_1btag ->Draw("hist");  
  h1_mggjj_nokinfit_1btag ->Draw("histsame");  
  leg->Draw();
  myC31.SaveAs("comp_1.png");

  TCanvas myC32("c32", "c32", 1);
  h1_mggjj_kinfit_2btag ->Draw("hist");  
  h1_mggjj_nokinfit_2btag ->Draw("histsame");  
  leg->Draw();
  myC32.SaveAs("comp_2.png");
}




