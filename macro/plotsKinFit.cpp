{
  gStyle->SetOptStat(0);
  
  TFile* file300=TFile::Open("../finalizedTrees_Radion_preselFitToMggbb/Radion_Radion_M-300_madgraph_default_CSV.root");
  TTree* tree300=(TTree*)file300->Get("myTrees");

  TFile* file500=TFile::Open("../finalizedTrees_Radion_preselFitToMggbb/Radion_Radion_M-500_madgraph_default_CSV.root");
  TTree* tree500=(TTree*)file500->Get("myTrees");

  TFile* file700=TFile::Open("../finalizedTrees_Radion_preselFitToMggbb/Radion_Radion_M-700_madgraph_default_CSV.root");
  TTree* tree700=(TTree*)file700->Get("myTrees");

  TFile* fileVH=TFile::Open("../finalizedTrees_Radion_preselFitToMggbb/Radion_WH_HToGG_M-125_8TeV_default_CSV.root");
  TTree* treeVH=(TTree*)fileVH->Get("myTrees");

  // for comparison w/wo regression
  TFile* file300_noReg_noKF=TFile::Open("../finalizedTrees_Radion_noRegression_noKinFit/Radion_Radion_M-300_regr_default_CSV.root");
  TTree* tree300_noReg_noKF=(TTree*)file300_noReg_noKF->Get("myTrees");

  TFile* file300_siReg_noKF=TFile::Open("../finalizedTrees_Radion_conRegression_noKinFit/Radion_Radion_M-300_regr_default_CSV.root");
  TTree* tree300_siReg_noKF=(TTree*)file300_siReg_noKF->Get("myTrees");

  TFile* file300_noReg_siKF=TFile::Open("../finalizedTrees_Radion_noRegression_conKinFit/Radion_Radion_M-300_regr_default_CSV.root");
  TTree* tree300_noReg_siKF=(TTree*)file300_noReg_siKF->Get("myTrees");

  TFile* file300_siReg_siKF=TFile::Open("../finalizedTrees_Radion_conRegression_conKinFit/Radion_Radion_M-300_regr_default_CSV.root");
  TTree* tree300_siReg_siKF=(TTree*)file300_siReg_siKF->Get("myTrees");


  // ======================================================
  // 1dim histos for signal and background 
  tree300 -> Draw("mggjj>>h1_1b_300_noKinFit(30,100.,700.)",  "(massggnewvtx>115 && massggnewvtx<135 && btagCategory==1)");
  tree300 -> Draw("mggjj_kin>>h1_1b_300_kinFit(30,100.,700.)","(massggnewvtx>115 && massggnewvtx<135 && btagCategory==1)");
  h1_1b_300_noKinFit -> Sumw2();
  h1_1b_300_kinFit   -> Sumw2();

  tree500 -> Draw("mggjj>>h1_1b_500_noKinFit(30,200.,800.)",  "(massggnewvtx>115 && massggnewvtx<135 && btagCategory==1)");
  tree500 -> Draw("mggjj_kin>>h1_1b_500_kinFit(30,200.,800.)","(massggnewvtx>115 && massggnewvtx<135 && btagCategory==1)");
  h1_1b_500_noKinFit -> Sumw2();
  h1_1b_500_kinFit   -> Sumw2();

  tree700 -> Draw("mggjj>>h1_1b_700_noKinFit(30,400.,1000.)",  "(massggnewvtx>115 && massggnewvtx<135 && btagCategory==1)");
  tree700 -> Draw("mggjj_kin>>h1_1b_700_kinFit(30,400.,1000.)","(massggnewvtx>115 && massggnewvtx<135 && btagCategory==1)");
  h1_1b_700_noKinFit -> Sumw2();
  h1_1b_700_kinFit   -> Sumw2();

  treeVH -> Draw("mggjj>>h1_1b_VH_noKinFit(30,100.,700.)",  "(massggnewvtx>115 && massggnewvtx<135 && btagCategory==1)");
  treeVH -> Draw("mggjj_kin>>h1_1b_VH_kinFit(30,100.,700.)","(massggnewvtx>115 && massggnewvtx<135 && btagCategory==1)");
  h1_1b_VH_noKinFit -> Sumw2();
  h1_1b_VH_kinFit   -> Sumw2();

  // for comparison w/wo regression
  tree300_noReg_noKF -> Draw("mggjj>>h1_1b_300_noReg_noKF(60,100.,700.)",  "(massggnewvtx>115 && massggnewvtx<135 && btagCategory==1)");
  tree300_siReg_noKF -> Draw("mggjj>>h1_1b_300_siReg_noKF(60,100.,700.)",  "(massggnewvtx>115 && massggnewvtx<135 && btagCategory==1)");
  tree300_noReg_siKF -> Draw("mggjj>>h1_1b_300_noReg_siKF(60,100.,700.)",  "(massggnewvtx>115 && massggnewvtx<135 && btagCategory==1)");
  tree300_siReg_siKF -> Draw("mggjj>>h1_1b_300_siReg_siKF(60,100.,700.)",  "(massggnewvtx>115 && massggnewvtx<135 && btagCategory==1)");
  h1_1b_300_noReg_noKF -> Sumw2();
  h1_1b_300_siReg_noKF -> Sumw2();
  h1_1b_300_noReg_siKF -> Sumw2();
  h1_1b_300_siReg_siKF -> Sumw2();

  // ======================================================
  tree300 -> Draw("mggjj>>h1_2b_300_noKinFit(45,100.,700.)",  "(massggnewvtx>115 && massggnewvtx<135 && btagCategory==2)");
  tree300 -> Draw("mggjj_kin>>h1_2b_300_kinFit(45,100.,700.)","(massggnewvtx>115 && massggnewvtx<135 && btagCategory==2)");
  h1_2b_300_noKinFit -> Sumw2();
  h1_2b_300_kinFit   -> Sumw2();

  tree500 -> Draw("mggjj>>h1_2b_500_noKinFit(45,200.,800.)",  "(massggnewvtx>115 && massggnewvtx<135 && btagCategory==2)");
  tree500 -> Draw("mggjj_kin>>h1_2b_500_kinFit(45,200.,800.)","(massggnewvtx>115 && massggnewvtx<135 && btagCategory==2)");
  h1_2b_500_noKinFit -> Sumw2();
  h1_2b_500_kinFit   -> Sumw2();

  tree700 -> Draw("mggjj>>h1_2b_700_noKinFit(45,300.,900.)",  "(massggnewvtx>115 && massggnewvtx<135 && btagCategory==2)");
  tree700 -> Draw("mggjj_kin>>h1_2b_700_kinFit(45,300.,900.)","(massggnewvtx>115 && massggnewvtx<135 && btagCategory==2)");
  h1_2b_700_noKinFit -> Sumw2();
  h1_2b_700_kinFit   -> Sumw2();

  treeVH -> Draw("mggjj>>h1_2b_VH_noKinFit(45,100.,700.)",  "(massggnewvtx>115 && massggnewvtx<135 && btagCategory==2)");
  treeVH -> Draw("mggjj_kin>>h1_2b_VH_kinFit(45,100.,700.)","(massggnewvtx>115 && massggnewvtx<135 && btagCategory==2)");
  h1_2b_VH_noKinFit -> Sumw2();
  h1_2b_VH_kinFit   -> Sumw2();

  // for comparison w/wo regression
  tree300_noReg_noKF -> Draw("mggjj>>h1_2b_300_noReg_noKF(60,100.,700.)",  "(massggnewvtx>115 && massggnewvtx<135 && btagCategory==2)");
  tree300_siReg_noKF -> Draw("mggjj>>h1_2b_300_siReg_noKF(60,100.,700.)",  "(massggnewvtx>115 && massggnewvtx<135 && btagCategory==2)");
  tree300_noReg_siKF -> Draw("mggjj>>h1_2b_300_noReg_siKF(60,100.,700.)",  "(massggnewvtx>115 && massggnewvtx<135 && btagCategory==2)");
  tree300_siReg_siKF -> Draw("mggjj>>h1_2b_300_siReg_siKF(60,100.,700.)",  "(massggnewvtx>115 && massggnewvtx<135 && btagCategory==2)");
  h1_2b_300_noReg_noKF -> Sumw2();
  h1_2b_300_siReg_noKF -> Sumw2();
  h1_2b_300_noReg_siKF -> Sumw2();
  h1_2b_300_siReg_siKF -> Sumw2();

  // ======================================================
  // 2dim plot
  tree300->Draw("mjj:mggjj>>h2D_noKinFitCorr_1b_300(100,100.,700.,100,0.,300.)","(massggnewvtx>115 && massggnewvtx<135 && btagCategory==1)","colz");
  tree300->Draw("mjj:mggjj_kin>>h2D_kinFitCorr_1b_300(100,100.,700.,100,0.,300.)","(massggnewvtx>115 && massggnewvtx<135 && btagCategory==1)","colz");
  tree300->Draw("mjj:mggjj>>h2D_noKinFitCorr_2b_300(100,100.,700.,100,0.,300.)","(massggnewvtx>115 && massggnewvtx<135 && btagCategory==2)","colz");
  tree300->Draw("mjj:mggjj_kin>>h2D_kinFitCorr_2b_300(100,100.,700.,100,0.,300.)","(massggnewvtx>115 && massggnewvtx<135 && btagCategory==2)","colz");

  // =======================================================================
  // cosmetics
  h1_1b_300_noKinFit -> GetXaxis()->SetTitle("m_{#gamma #gamma jj} [GeV]");
  h1_1b_500_noKinFit -> GetXaxis()->SetTitle("m_{#gamma #gamma jj} [GeV]");
  h1_1b_700_noKinFit -> GetXaxis()->SetTitle("m_{#gamma #gamma jj} [GeV]");
  h1_1b_VH_noKinFit  -> GetXaxis()->SetTitle("m_{#gamma #gamma jj} [GeV]");

  h1_2b_300_noKinFit -> GetXaxis()->SetTitle("m_{#gamma #gamma jj} [GeV]");
  h1_2b_500_noKinFit -> GetXaxis()->SetTitle("m_{#gamma #gamma jj} [GeV]");
  h1_2b_700_noKinFit -> GetXaxis()->SetTitle("m_{#gamma #gamma jj} [GeV]");
  h1_2b_VH_noKinFit  -> GetXaxis()->SetTitle("m_{#gamma #gamma jj} [GeV]");

  h1_1b_300_kinFit -> GetXaxis()->SetTitle("m_{#gamma #gamma jj} [GeV]");
  h1_1b_500_kinFit -> GetXaxis()->SetTitle("m_{#gamma #gamma jj} [GeV]");
  h1_1b_700_kinFit -> GetXaxis()->SetTitle("m_{#gamma #gamma jj} [GeV]");
  h1_1b_VH_kinFit  -> GetXaxis()->SetTitle("m_{#gamma #gamma jj} [GeV]");

  h1_2b_300_kinFit -> GetXaxis()->SetTitle("m_{#gamma #gamma jj} [GeV]");
  h1_2b_500_kinFit -> GetXaxis()->SetTitle("m_{#gamma #gamma jj} [GeV]");
  h1_2b_700_kinFit -> GetXaxis()->SetTitle("m_{#gamma #gamma jj} [GeV]");
  h1_2b_VH_kinFit  -> GetXaxis()->SetTitle("m_{#gamma #gamma jj} [GeV]");

  h1_1b_300_noReg_noKF -> GetXaxis()->SetTitle("m_{#gamma #gamma jj} [GeV]");
  h1_1b_300_siReg_noKF -> GetXaxis()->SetTitle("m_{#gamma #gamma jj} [GeV]");
  h1_1b_300_noReg_siKF -> GetXaxis()->SetTitle("m_{#gamma #gamma jj} [GeV]");
  h1_1b_300_siReg_siKF -> GetXaxis()->SetTitle("m_{#gamma #gamma jj} [GeV]");

  h1_2b_300_noReg_noKF -> GetXaxis()->SetTitle("m_{#gamma #gamma jj} [GeV]");
  h1_2b_300_siReg_noKF -> GetXaxis()->SetTitle("m_{#gamma #gamma jj} [GeV]");
  h1_2b_300_noReg_siKF -> GetXaxis()->SetTitle("m_{#gamma #gamma jj} [GeV]");
  h1_2b_300_siReg_siKF -> GetXaxis()->SetTitle("m_{#gamma #gamma jj} [GeV]");

  h1_1b_300_noKinFit -> SetFillColor(0);
  h1_1b_500_noKinFit -> SetFillColor(0);
  h1_1b_700_noKinFit -> SetFillColor(0);
  h1_1b_VH_noKinFit  -> SetFillColor(0);
  h1_1b_300_noKinFit -> SetLineColor(1);
  h1_1b_500_noKinFit -> SetLineColor(1);
  h1_1b_700_noKinFit -> SetLineColor(1);
  h1_1b_VH_noKinFit  -> SetLineColor(1);
  h1_1b_300_noKinFit -> SetLineWidth(2);
  h1_1b_500_noKinFit -> SetLineWidth(2);
  h1_1b_700_noKinFit -> SetLineWidth(2);
  h1_1b_VH_noKinFit  -> SetLineWidth(2);

  h1_2b_300_noKinFit -> SetFillColor(0);
  h1_2b_500_noKinFit -> SetFillColor(0);
  h1_2b_700_noKinFit -> SetFillColor(0);
  h1_2b_VH_noKinFit  -> SetFillColor(0);
  h1_2b_300_noKinFit -> SetLineColor(1);
  h1_2b_500_noKinFit -> SetLineColor(1);
  h1_2b_700_noKinFit -> SetLineColor(1);
  h1_2b_VH_noKinFit  -> SetLineColor(1);
  h1_2b_300_noKinFit -> SetLineWidth(2);
  h1_2b_500_noKinFit -> SetLineWidth(2);
  h1_2b_700_noKinFit -> SetLineWidth(2);
  h1_2b_VH_noKinFit  -> SetLineWidth(2);

  h1_1b_300_kinFit -> SetFillColor(0);
  h1_1b_500_kinFit -> SetFillColor(0);
  h1_1b_700_kinFit -> SetFillColor(0);
  h1_1b_VH_kinFit  -> SetFillColor(0);
  h1_1b_300_kinFit -> SetLineColor(2);
  h1_1b_500_kinFit -> SetLineColor(2);
  h1_1b_700_kinFit -> SetLineColor(2);
  h1_1b_VH_kinFit  -> SetLineColor(2);
  h1_1b_300_kinFit -> SetLineWidth(2);
  h1_1b_500_kinFit -> SetLineWidth(2);
  h1_1b_700_kinFit -> SetLineWidth(2);
  h1_1b_VH_kinFit  -> SetLineWidth(2);

  h1_2b_300_kinFit -> SetFillColor(0);
  h1_2b_500_kinFit -> SetFillColor(0);
  h1_2b_700_kinFit -> SetFillColor(0);
  h1_2b_VH_kinFit  -> SetFillColor(0);
  h1_2b_300_kinFit -> SetLineColor(2);
  h1_2b_500_kinFit -> SetLineColor(2);
  h1_2b_700_kinFit -> SetLineColor(2);
  h1_2b_VH_kinFit  -> SetLineColor(2);
  h1_2b_300_kinFit -> SetLineWidth(2);
  h1_2b_500_kinFit -> SetLineWidth(2);
  h1_2b_700_kinFit -> SetLineWidth(2);
  h1_2b_VH_kinFit  -> SetLineWidth(2);

  h1_1b_300_noReg_noKF -> SetFillColor(0);
  h1_1b_300_siReg_noKF -> SetFillColor(0);
  h1_1b_300_noReg_siKF -> SetFillColor(0);
  h1_1b_300_siReg_siKF -> SetFillColor(0);
  h1_1b_300_noReg_noKF -> SetLineColor(1);
  h1_1b_300_siReg_noKF -> SetLineColor(2);
  h1_1b_300_noReg_siKF -> SetLineColor(3);
  h1_1b_300_siReg_siKF -> SetLineColor(4);
  h1_1b_300_noReg_noKF -> SetLineWidth(2);
  h1_1b_300_siReg_noKF -> SetLineWidth(2);
  h1_1b_300_noReg_siKF -> SetLineWidth(2);
  h1_1b_300_siReg_siKF -> SetLineWidth(2);

  h1_2b_300_noReg_noKF -> SetFillColor(0);
  h1_2b_300_siReg_noKF -> SetFillColor(0);
  h1_2b_300_noReg_siKF -> SetFillColor(0);
  h1_2b_300_siReg_siKF -> SetFillColor(0);
  h1_2b_300_noReg_noKF -> SetLineColor(1);
  h1_2b_300_siReg_noKF -> SetLineColor(2);
  h1_2b_300_noReg_siKF -> SetLineColor(3);
  h1_2b_300_siReg_siKF -> SetLineColor(4);
  h1_2b_300_noReg_noKF -> SetLineWidth(2);
  h1_2b_300_siReg_noKF -> SetLineWidth(2);
  h1_2b_300_noReg_siKF -> SetLineWidth(2);
  h1_2b_300_siReg_siKF -> SetLineWidth(2);

  // ======================================================
  h2D_noKinFitCorr_1b_300 -> GetXaxis()->SetTitle("m_{jj} [GeV]");
  h2D_noKinFitCorr_1b_300 -> GetYaxis()->SetTitle("m_{#gamma #gamma jj} [GeV]");

  h2D_noKinFitCorr_2b_300 -> GetXaxis()->SetTitle("m_{jj} [GeV]");
  h2D_noKinFitCorr_2b_300 -> GetYaxis()->SetTitle("m_{#gamma #gamma jj} [GeV]");

  h2D_kinFitCorr_1b_300 -> GetXaxis()->SetTitle("m_{jj} [GeV]");
  h2D_kinFitCorr_1b_300 -> GetYaxis()->SetTitle("m_{#gamma #gamma jj} [GeV]");

  h2D_kinFitCorr_2b_300 -> GetXaxis()->SetTitle("m_{jj} [GeV]");
  h2D_kinFitCorr_2b_300 -> GetYaxis()->SetTitle("m_{#gamma #gamma jj} [GeV]");

  // ======================================================
  TCanvas c1a("c1a","m=300, 1btag", 1);
  h1_1b_300_kinFit->Draw("hist"); 
  h1_1b_300_noKinFit->Draw("samehist");
  c1a->SaveAs("signal1btag_300.png");

  TCanvas c1b("c1b","m=300, 2btag", 1);
  h1_2b_300_kinFit->Draw("hist");
  h1_2b_300_noKinFit->Draw("samehist");
  c1b->SaveAs("signal2btag_300.png");

  TCanvas c2a("c2a","m=500, 1btag", 1);
  h1_1b_500_kinFit->Draw("hist");
  h1_1b_500_noKinFit->Draw("samehist");
  c2a->SaveAs("signal1btag_500.png");

  TCanvas c2b("c2b","m=500, 2btag", 1);
  h1_2b_500_kinFit->Draw("hist");
  h1_2b_500_noKinFit->Draw("samehist");
  c2b->SaveAs("signal2btag_500.png");

  TCanvas c3a("c3a","m=700, 1btag", 1);
  h1_1b_700_kinFit->Draw("hist");
  h1_1b_700_noKinFit->Draw("samehist");
  c3a->SaveAs("signal1btag_700.png");

  TCanvas c3b("c3b","m=700, 2btag", 1);
  h1_2b_700_kinFit->Draw("hist");
  h1_2b_700_noKinFit->Draw("samehist");
  c3b->SaveAs("signal2btag_700.png");

  TCanvas c4a("c4a","VH, 1btag", 1);
  h1_1b_VH_noKinFit->Draw("hist");
  h1_1b_VH_kinFit->Draw("samehist");
  c4a->SaveAs("signal1btag_VH.png");

  TCanvas c4b("c4b","VH, 2btag", 1);
  h1_2b_VH_kinFit->Draw("hist");
  h1_2b_VH_noKinFit->Draw("samehist");
  c4b->SaveAs("signal2btag_VH.png");

  TCanvas c5a("c5a","m=300, 1btag", 1);
  h2D_noKinFitCorr_1b_300 -> Draw("colz");
  c5a->SaveAs("corrNoKinFit1btag_300.png");

  TCanvas c5b("c5b","m=300, 1btag", 1);
  h2D_kinFitCorr_1b_300 -> Draw("colz");
  c5b->SaveAs("corrKinFit1btag_300.png");

  TCanvas c5c("c5c","m=300, 2btag", 1);
  h2D_noKinFitCorr_2b_300 -> Draw("colz");
  c5c->SaveAs("corrNoKinFit2btag_300.png");

  TCanvas c5d("c5d","m=300, 2btag", 1);
  h2D_kinFitCorr_2b_300 -> Draw("colz");
  c5d->SaveAs("corrKinFit2btag_300.png");

  TCanvas c6a("c6a","1btag", 1);
  h1_1b_300_kinFit->SetLineColor(2);
  h1_1b_300_kinFit->DrawNormalized("hist"); 
  h1_1b_VH_kinFit->SetLineColor(4);
  h1_1b_VH_kinFit->DrawNormalized("samehist");
  c6a->SaveAs("signalVsVH_1btag_kinFit.png");

  TCanvas c6b("c6b","2btag", 1);
  h1_2b_300_kinFit->SetLineColor(2);
  h1_2b_300_kinFit->DrawNormalized("hist"); 
  h1_2b_VH_kinFit->SetLineColor(4);
  h1_2b_VH_kinFit->DrawNormalized("samehist");
  c6b->SaveAs("signalVsVH_2btag_kinFit.png");

  TCanvas c6c("c6c","1btag", 1);
  h1_1b_300_noKinFit->SetLineColor(2);
  h1_1b_300_noKinFit->DrawNormalized("hist"); 
  h1_1b_VH_noKinFit->SetLineColor(4);
  h1_1b_VH_noKinFit->DrawNormalized("samehist");
  c6c->SaveAs("signalVsVH_1btag_noKinFit.png");

  TCanvas c6d("c6d","2btag", 1);
  h1_2b_300_noKinFit->SetLineColor(2);
  h1_2b_300_noKinFit->DrawNormalized("hist"); 
  h1_2b_VH_noKinFit->SetLineColor(4);
  h1_2b_VH_noKinFit->DrawNormalized("samehist");
  c6d->SaveAs("signalVsVH_2btag_noKinFit.png");


  // w/wo kinfit and regression
  TLegend *leg;
  leg = new TLegend(0.5,0.5,0.75,0.75);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.05);
  leg->SetFillColor(0);
  leg->AddEntry(h1_1b_300_siReg_siKF, "regression+kinfit",   "l");
  leg->AddEntry(h1_1b_300_noReg_noKF, "no regression, no kinfit",   "l");
  leg->AddEntry(h1_1b_300_noReg_siKF, "kinfit only",    "l");
  leg->AddEntry(h1_1b_300_siReg_noKF, "regression only","l");

  TCanvas c7a("c7a","1btag", 1);
  h1_1b_300_siReg_siKF->DrawNormalized("hist");
  h1_1b_300_siReg_noKF->DrawNormalized("samehist");
  h1_1b_300_noReg_siKF->DrawNormalized("samehist");
  h1_1b_300_noReg_noKF->DrawNormalized("samehist");
  leg->Draw();
  c7a->SaveAs("CompWithRegression_1btag.png");

  TCanvas c7b("c7b","2btag", 1);
  h1_2b_300_siReg_siKF->DrawNormalized("hist");
  h1_2b_300_siReg_noKF->DrawNormalized("samehist");
  h1_2b_300_noReg_siKF->DrawNormalized("samehist");
  h1_2b_300_noReg_noKF->DrawNormalized("samehist");
  leg->Draw();
  c7b->SaveAs("CompWithRegression_2btag.png");
}
