{
  gStyle->SetOptStat(0);
  
  TFile* fileConR=TFile::Open("../finalizedTrees_m300_fitToMggbb_perOptim_conRegression/Radion_Radion_M-300_regr_default_CSV.root");
  TTree* treeConR=(TTree*)fileConR->Get("myTrees");

  TFile* fileNoR=TFile::Open("../finalizedTrees_m300_fitToMggbb_perOptim_noRegressionForCompar/Radion_Radion_M-300_regr_default_CSV.root");
  TTree* treeNoR=(TTree*)fileNoR->Get("myTrees");

  // ======================================================
  treeConR -> Draw("mggjj>>h1_mggjj_1b_conR(30,100.,700.)",  "btagCategory==1");
  treeConR -> Draw("mggjj>>h1_mggjj_2b_conR(30,100.,700.)",  "btagCategory==2");
  h1_mggjj_1b_conR -> Sumw2();
  h1_mggjj_2b_conR -> Sumw2();
  h1_mggjj_1b_conR -> SetLineWidth(2);   
  h1_mggjj_2b_conR -> SetLineWidth(2);   
  h1_mggjj_1b_conR -> SetLineColor(2);
  h1_mggjj_2b_conR -> SetLineColor(2);
  h1_mggjj_1b_conR -> SetFillColor(0);
  h1_mggjj_2b_conR -> SetFillColor(0);
  h1_mggjj_1b_conR -> GetXaxis()->SetTitle("m_{#gamma #gamma jj} [GeV]");
  h1_mggjj_2b_conR -> GetXaxis()->SetTitle("m_{#gamma #gamma jj} [GeV]");

  treeNoR -> Draw("mggjj>>h1_mggjj_1b_noR(30,100.,700.)",  "btagCategory==1");
  treeNoR -> Draw("mggjj>>h1_mggjj_2b_noR(30,100.,700.)",  "btagCategory==2");
  h1_mggjj_1b_noR -> Sumw2();
  h1_mggjj_2b_noR -> Sumw2();
  h1_mggjj_1b_noR -> SetLineWidth(2);   
  h1_mggjj_2b_noR -> SetLineWidth(2);   
  h1_mggjj_1b_noR -> SetLineColor(4);
  h1_mggjj_2b_noR -> SetLineColor(4);
  h1_mggjj_1b_noR -> SetFillColor(0);
  h1_mggjj_2b_noR -> SetFillColor(0);
  h1_mggjj_1b_noR -> GetXaxis()->SetTitle("m_{#gamma #gamma jj} [GeV]");
  h1_mggjj_2b_noR -> GetXaxis()->SetTitle("m_{#gamma #gamma jj} [GeV]");

  treeConR -> Draw("mjj>>h1_mjj_1b_conR(30,50.,300.)",  "btagCategory==1");
  treeConR -> Draw("mjj>>h1_mjj_2b_conR(30,50.,300.)",  "btagCategory==2");
  h1_mjj_1b_conR -> Sumw2();
  h1_mjj_2b_conR -> Sumw2();
  h1_mjj_1b_conR -> SetLineWidth(2);   
  h1_mjj_2b_conR -> SetLineWidth(2);   
  h1_mjj_1b_conR -> SetLineColor(2);
  h1_mjj_2b_conR -> SetLineColor(2);
  h1_mjj_1b_conR -> SetFillColor(0);
  h1_mjj_2b_conR -> SetFillColor(0);
  h1_mjj_1b_conR -> GetXaxis()->SetTitle("m_{jj} [GeV]");
  h1_mjj_2b_conR -> GetXaxis()->SetTitle("m_{jj} [GeV]");

  treeNoR -> Draw("mjj>>h1_mjj_1b_noR(30,50.,300.)",  "btagCategory==1");
  treeNoR -> Draw("mjj>>h1_mjj_2b_noR(30,50.,300.)",  "btagCategory==2");
  h1_mjj_1b_noR -> Sumw2();
  h1_mjj_2b_noR -> Sumw2();
  h1_mjj_1b_noR -> SetLineWidth(2);   
  h1_mjj_2b_noR -> SetLineWidth(2);   
  h1_mjj_1b_noR -> SetLineColor(4);
  h1_mjj_2b_noR -> SetLineColor(4);
  h1_mjj_1b_noR -> SetFillColor(0);
  h1_mjj_2b_noR -> SetFillColor(0);
  h1_mjj_1b_noR -> GetXaxis()->SetTitle("m_{jj} [GeV]");
  h1_mjj_2b_noR -> GetXaxis()->SetTitle("m_{jj} [GeV]");



  // ======================================================
  TLegend *leg;
  leg = new TLegend(0.5,0.5,0.75,0.75);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.05);
  leg->SetFillColor(0);
  leg->AddEntry(h1_mggjj_2b_conR, "with regression", "l");
  leg->AddEntry(h1_mggjj_2b_noR,  "no regression",   "l");

  TCanvas c1a("c1a","m=300, 1btag", 1);
  h1_mjj_1b_noR ->Draw("hist"); 
  h1_mjj_1b_conR->Draw("samehist");
  leg->Draw();
  c1a->SaveAs("mjj_1btag.png");
  c1a->SaveAs("mjj_1btag.pdf");

  TCanvas c1b("c1b","m=300, 2btag", 1);
  h1_mjj_2b_conR->Draw("hist");
  h1_mjj_2b_noR ->Draw("samehist"); 
  leg->Draw();
  c1b->SaveAs("mjj_2btag.png");
  c1b->SaveAs("mjj_2btag.pdf");

  TCanvas c2a("c2a","m=300, 1btag", 1);
  h1_mggjj_1b_conR->Draw("hist"); 
  h1_mggjj_1b_noR ->Draw("samehist");
  leg->Draw();
  c2a->SaveAs("mggjj_1btag.png");
  c2a->SaveAs("mggjj_1btag.pdf");

  TCanvas c2b("c2b","m=300, 2btag", 1);
  h1_mggjj_2b_conR->Draw("hist"); 
  h1_mggjj_2b_noR ->Draw("samehist");
  leg->Draw();
  c2b->SaveAs("mggjj_2btag.png");
  c2b->SaveAs("mggjj_2btag.pdf");
}
