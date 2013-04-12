#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH2D.h"
#include "TLegend.h"
#include <iostream>

void doPlotUL() {

  float eff[5] = { 0.5, 0.6, 0.7, 0.8, 0.9 };

  float ul_1b_3var[5] = { 98.5612, 90.5495, 84.6464, 85.6297, 85.4228 };
  float ul_1b_4var[5] = { 88.8435, 83.8169, 82.3924, 83.135,  83.8256 };
  float ul_2b_3var[5] = { 25.2893, 23.6267, 22.1692, 21.5376, 23.3353 };
  float ul_2b_4var[5] = { 24.4318, 23.3642, 21.7822, 23.2495, 23.9495 } ;

  float err[5] = { 0., 0., 0., 0., 0. };

  TGraphErrors *grUL_1b_3var = new TGraphErrors(5, eff, ul_1b_3var, err, err);
  TGraphErrors *grUL_1b_4var = new TGraphErrors(5, eff, ul_1b_4var, err, err);
  TGraphErrors *grUL_2b_3var = new TGraphErrors(5, eff, ul_2b_3var, err, err);
  TGraphErrors *grUL_2b_4var = new TGraphErrors(5, eff, ul_2b_4var, err, err);

  grUL_1b_3var -> SetMarkerColor(2);
  grUL_1b_3var -> SetMarkerStyle(20);
  grUL_1b_4var -> SetMarkerColor(4);
  grUL_1b_4var -> SetMarkerStyle(21);

  grUL_2b_3var -> SetMarkerColor(2);
  grUL_2b_3var -> SetMarkerStyle(20);
  grUL_2b_4var -> SetMarkerColor(4);
  grUL_2b_4var -> SetMarkerStyle(21);

  TH2F *myH_1b = new TH2F("myH_1b","myH_1b",100, 0.45, 0.95, 100, 80., 100.);
  TH2F *myH_2b = new TH2F("myH_2b","myH_2b",100, 0.45, 0.95, 100, 15., 35.);

  myH_1b->SetTitle("1b-jet category");
  myH_1b->GetXaxis()->SetTitle("signal eff");
  myH_1b->GetYaxis()->SetTitle("UL/xsec");

  myH_2b->SetTitle("2b-jet category");
  myH_2b->GetXaxis()->SetTitle("signal eff");
  myH_2b->GetYaxis()->SetTitle("UL/xsec");

  gStyle->SetOptStat(0);

  TLegend *leg;
  leg = new TLegend(0.6,0.6,0.85,0.85);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.05);
  leg->SetFillColor(0);
  leg->AddEntry(grUL_1b_3var, "3 variables", "lmp");
  leg->AddEntry(grUL_1b_4var, "4 variables", "lmp");

  TCanvas myC1("c_1b", "c_1b", 1);
  myH_1b->Draw();
  grUL_1b_3var->Draw("Psame");
  grUL_1b_4var->Draw("Psame");
  leg->Draw();
  myC1.SaveAs("UL_1b.png");

  TCanvas myC2("c_2b", "c_2b", 1);
  myH_2b->Draw();
  grUL_2b_3var->Draw("Psame");
  grUL_2b_4var->Draw("Psame");
  leg->Draw();
  myC2.SaveAs("UL_2b.png");
  

}




