{
  gStyle->SetOptStat(0);
  
  TFile* dataFile=TFile::Open("../finalizedTrees_Radion_presel/Radion_Data2012_default_CSV.root");
  TTree* tree_data=(TTree*)dataFile->Get("myTrees");

  TFile* csFile=TFile::Open("../finalizedTrees_Radion_presel_CS/Radion_Data2012_default_CSV.root");
  TTree* tree_cs=(TTree*)csFile->Get("myTrees");

  TFile* csFileWeightG=TFile::Open("treesFromCS_presel_withWeights/csWithWeightFromGammas.root");
  TTree* tree_cswG=(TTree*)csFileWeightG->Get("myTrees_withWeight");

  TFile* csFileWeightJ=TFile::Open("treesFromCS_presel_withWeights/csWithWeightFromJets.root");
  TTree* tree_cswJ=(TTree*)csFileWeightJ->Get("myTrees_withWeight");


  // ------------------------------------------------------------------------------  
  // pre-weight, to understand what to rescale
  tree_data -> Draw("ptPhot1>>h1_ptGamma2_data(35,20.,160.)", "weight*(massggnewvtx>100 && massggnewvtx<180 && (massggnewvtx>135 || massggnewvtx<115))");
  tree_cs   -> Draw("ptPhot1>>h1_ptGamma2_cs(35,20.,160.)",   "weight*(massggnewvtx>100 && massggnewvtx<180)");
  tree_cswG -> Draw("ptPhot1>>h1_ptGamma2_csG(35,20.,160.)",  "pt_scaled_2D_weight_data*(massggnewvtx>100 && massggnewvtx<180)");
  tree_cswJ -> Draw("ptPhot1>>h1_ptGamma2_csJ(35,20.,160.)",  "pt_scaled_2D_weight_data*(massggnewvtx>100 && massggnewvtx<180)");
  //
  tree_data -> Draw("ptPhot2>>h1_ptGamma1_data(35,20.,160.)", "weight*(massggnewvtx>100 && massggnewvtx<180 && (massggnewvtx>135 || massggnewvtx<115))");
  tree_cs   -> Draw("ptPhot2>>h1_ptGamma1_cs(35,20.,160.)",   "weight*(massggnewvtx>100 && massggnewvtx<180)");
  tree_cswG -> Draw("ptPhot2>>h1_ptGamma1_csG(35,20.,160.)",  "pt_scaled_2D_weight_data*(massggnewvtx>100 && massggnewvtx<180)");
  tree_cswJ -> Draw("ptPhot2>>h1_ptGamma1_csJ(35,20.,160.)",  "pt_scaled_2D_weight_data*(massggnewvtx>100 && massggnewvtx<180)");
  //
  tree_data -> Draw("ptCorrJet1>>h1_ptJet2_data(20,20.,160.)", "weight*(massggnewvtx>100 && massggnewvtx<180 && (massggnewvtx>135 || massggnewvtx<115))");
  tree_cs   -> Draw("ptCorrJet1>>h1_ptJet2_cs(20,20.,160.)",   "weight*(massggnewvtx>100 && massggnewvtx<180)");
  tree_cswG -> Draw("ptCorrJet1>>h1_ptJet2_csG(20,20.,160.)",  "pt_scaled_2D_weight_data*(massggnewvtx>100 && massggnewvtx<180)");
  tree_cswJ -> Draw("ptCorrJet1>>h1_ptJet2_csJ(20,20.,160.)",  "pt_scaled_2D_weight_data*(massggnewvtx>100 && massggnewvtx<180)");
  //
  tree_data -> Draw("ptCorrJet2>>h1_ptJet1_data(35,20.,160.)", "weight*(massggnewvtx>100 && massggnewvtx<180 && (massggnewvtx>135 || massggnewvtx<115))");
  tree_cs   -> Draw("ptCorrJet2>>h1_ptJet1_cs(35,20.,160.)",   "weight*(massggnewvtx>100 && massggnewvtx<180)");
  tree_cswG -> Draw("ptCorrJet2>>h1_ptJet1_csG(35,20.,160.)",  "pt_scaled_2D_weight_data*(massggnewvtx>100 && massggnewvtx<180)");
  tree_cswJ -> Draw("ptCorrJet2>>h1_ptJet1_csJ(35,20.,160.)",  "pt_scaled_2D_weight_data*(massggnewvtx>100 && massggnewvtx<180)");
  // 
  tree_data -> Draw("HT_jet>>h1_HT_jet_data(25,0.,400.)", "weight*(massggnewvtx>100 && massggnewvtx<180 && (massggnewvtx>135 || massggnewvtx<115))");
  tree_cs   -> Draw("HT_jet>>h1_HT_jet_cs(25,0.,400.)",   "weight*(massggnewvtx>100 && massggnewvtx<180)");
  tree_cswG -> Draw("HT_jet>>h1_HT_jet_csG(25,0.,400.)",  "pt_scaled_2D_weight_data*(massggnewvtx>100 && massggnewvtx<180)");
  tree_cswJ -> Draw("HT_jet>>h1_HT_jet_csJ(25,0.,400.)",  "pt_scaled_2D_weight_data*(massggnewvtx>100 && massggnewvtx<180)");
  //
  tree_data -> Draw("deltaphiggjj>>h1_deltaphiggjj_data(20,0.,3.14)", "weight*(massggnewvtx>100 && massggnewvtx<180 && (massggnewvtx>135 || massggnewvtx<115))");
  tree_cs   -> Draw("deltaphiggjj>>h1_deltaphiggjj_cs(20,0.,3.14)",   "weight*(massggnewvtx>100 && massggnewvtx<180)");
  tree_cswG -> Draw("deltaphiggjj>>h1_deltaphiggjj_csG(20,0.,3.14)",  "weight*pt_scaled_2D_weight_data*(massggnewvtx>100 && massggnewvtx<180)");
  tree_cswJ -> Draw("deltaphiggjj>>h1_deltaphiggjj_csJ(20,0.,3.14)",  "weight*pt_scaled_2D_weight_data*(massggnewvtx>100 && massggnewvtx<180)");
  //
  tree_data -> Draw("btagCategory>>h1_btag_data(4,0.,4.)", "weight*(massggnewvtx>100 && massggnewvtx<180 && (massggnewvtx>135 || massggnewvtx<115))");
  tree_cs   -> Draw("btagCategory>>h1_btag_cs(4,0.,4.)",   "weight*(massggnewvtx>100 && massggnewvtx<180)");
  tree_cswG -> Draw("btagCategory>>h1_btag_csG(4,0.,4.)",  "weight*pt_scaled_2D_weight_data*(massggnewvtx>100 && massggnewvtx<180)");
  tree_cswJ -> Draw("btagCategory>>h1_btag_csJ(4,0.,4.)",  "weight*pt_scaled_2D_weight_data*(massggnewvtx>100 && massggnewvtx<180)");

  h1_ptGamma1_data     -> Sumw2();  h1_ptGamma1_cs     -> Sumw2();  h1_ptGamma1_csG     -> Sumw2();  h1_ptGamma1_csJ     -> Sumw2();
  h1_ptGamma2_data     -> Sumw2();  h1_ptGamma2_cs     -> Sumw2();  h1_ptGamma2_csG     -> Sumw2();  h1_ptGamma2_csJ     -> Sumw2();
  h1_ptJet1_data       -> Sumw2();  h1_ptJet1_cs       -> Sumw2();  h1_ptJet1_csG       -> Sumw2();  h1_ptJet1_csJ       -> Sumw2();
  h1_ptJet2_data       -> Sumw2();  h1_ptJet2_cs       -> Sumw2();  h1_ptJet2_csG       -> Sumw2();  h1_ptJet2_csJ       -> Sumw2();
  h1_HT_jet_data       -> Sumw2();  h1_HT_jet_cs       -> Sumw2();  h1_HT_jet_csG       -> Sumw2();  h1_HT_jet_csJ       -> Sumw2();
  h1_deltaphiggjj_data -> Sumw2();  h1_deltaphiggjj_cs -> Sumw2();  h1_deltaphiggjj_csG -> Sumw2();  h1_deltaphiggjj_csJ -> Sumw2();
  h1_btag_data         -> Sumw2();  h1_btag_cs         -> Sumw2();  h1_btag_csG         -> Sumw2();  h1_btag_csJ         -> Sumw2();

  h1_ptGamma1_data -> GetXaxis()->SetTitle("p_{T} (gamma1) [GeV]");
  h1_ptGamma2_data -> GetXaxis()->SetTitle("p_{T} (gamma2) [GeV]");
  h1_ptGamma1_cs   -> GetXaxis()->SetTitle("p_{T} (gamma1) [GeV]");
  h1_ptGamma2_cs   -> GetXaxis()->SetTitle("p_{T} (gamma2) [GeV]");
  h1_ptGamma1_csG  -> GetXaxis()->SetTitle("p_{T} (gamma1) [GeV]");
  h1_ptGamma2_csG  -> GetXaxis()->SetTitle("p_{T} (gamma2) [GeV]");
  h1_ptGamma1_csJ  -> GetXaxis()->SetTitle("p_{T} (gamma1) [GeV]");
  h1_ptGamma2_csJ  -> GetXaxis()->SetTitle("p_{T} (gamma2) [GeV]");

  h1_ptJet1_data -> GetXaxis()->SetTitle("p_{T} (jet1) [GeV]");
  h1_ptJet2_data -> GetXaxis()->SetTitle("p_{T} (jet2) [GeV]");
  h1_ptJet1_cs   -> GetXaxis()->SetTitle("p_{T} (jet1) [GeV]");
  h1_ptJet2_cs   -> GetXaxis()->SetTitle("p_{T} (jet2) [GeV]");
  h1_ptJet1_csG  -> GetXaxis()->SetTitle("p_{T} (jet1) [GeV]");
  h1_ptJet2_csG  -> GetXaxis()->SetTitle("p_{T} (jet2) [GeV]");
  h1_ptJet1_csJ  -> GetXaxis()->SetTitle("p_{T} (jet1) [GeV]");
  h1_ptJet2_csJ  -> GetXaxis()->SetTitle("p_{T} (jet2) [GeV]");

  h1_HT_jet_data -> GetXaxis()->SetTitle("p_{T} (jet1) + p_{T} (jet2) [GeV]");
  h1_HT_jet_cs   -> GetXaxis()->SetTitle("p_{T} (jet1) + p_{T} (jet2) [GeV]");
  h1_HT_jet_csG  -> GetXaxis()->SetTitle("p_{T} (jet1) + p_{T} (jet2) [GeV]");
  h1_HT_jet_csJ  -> GetXaxis()->SetTitle("p_{T} (jet1) + p_{T} (jet2) [GeV]");

  h1_deltaphiggjj_data -> GetXaxis()->SetTitle("#Delta #phi (H,H)");
  h1_deltaphiggjj_cs   -> GetXaxis()->SetTitle("#Delta #phi (H,H)");
  h1_deltaphiggjj_csG  -> GetXaxis()->SetTitle("#Delta #phi (H,H)");
  h1_deltaphiggjj_csJ  -> GetXaxis()->SetTitle("#Delta #phi (H,H)");

  h1_btag_data -> GetXaxis()->SetTitle("b-tag category");
  h1_btag_cs   -> GetXaxis()->SetTitle("b-tag category");
  h1_btag_csG  -> GetXaxis()->SetTitle("b-tag category");
  h1_btag_csJ  -> GetXaxis()->SetTitle("b-tag category");

  h1_ptGamma1_data -> SetMarkerStyle(20); 
  h1_ptGamma2_data -> SetMarkerStyle(20); 
  h1_ptJet1_data   -> SetMarkerStyle(20); 
  h1_ptJet2_data   -> SetMarkerStyle(20); 

  h1_ptGamma1_data -> SetMarkerSize(0.8); 
  h1_ptGamma2_data -> SetMarkerSize(0.8); 
  h1_ptJet1_data   -> SetMarkerSize(0.8); 
  h1_ptJet2_data   -> SetMarkerSize(0.8); 

  h1_HT_jet_data       -> SetMarkerStyle(20); 
  h1_HT_jet_data       -> SetMarkerSize(0.8); 
  h1_deltaphiggjj_data -> SetMarkerStyle(20); 
  h1_deltaphiggjj_data -> SetMarkerSize(0.8); 
  h1_btag_data         -> SetMarkerStyle(20); 
  h1_btag_data         -> SetMarkerSize(0.8); 

  h1_ptGamma1_cs     -> SetFillColor(kRed-9);  h1_ptGamma1_cs     -> SetLineColor(kRed+3); 
  h1_ptGamma2_cs     -> SetFillColor(kRed-9);  h1_ptGamma2_cs     -> SetLineColor(kRed+3); 
  h1_ptJet1_cs       -> SetFillColor(kRed-9);  h1_ptJet1_cs       -> SetLineColor(kRed+3); 
  h1_ptJet2_cs       -> SetFillColor(kRed-9);  h1_ptJet2_cs       -> SetLineColor(kRed+3); 
  h1_HT_jet_cs       -> SetFillColor(kRed-9);  h1_HT_jet_cs       -> SetLineColor(kRed+3); 
  h1_deltaphiggjj_cs -> SetFillColor(kRed-9);  h1_deltaphiggjj_cs -> SetLineColor(kRed+3); 
  h1_btag_cs         -> SetFillColor(kRed-9);  h1_btag_cs         -> SetLineColor(kRed+3); 

  h1_ptGamma1_csG     -> SetFillColor(kAzure-9);  h1_ptGamma1_csG     -> SetLineColor(kAzure+3); 
  h1_ptGamma2_csG     -> SetFillColor(kAzure-9);  h1_ptGamma2_csG     -> SetLineColor(kAzure+3); 
  h1_ptJet1_csG       -> SetFillColor(kAzure-9);  h1_ptJet1_csG       -> SetLineColor(kAzure+3); 
  h1_ptJet2_csG       -> SetFillColor(kAzure-9);  h1_ptJet2_csG       -> SetLineColor(kAzure+3); 
  h1_HT_jet_csG       -> SetFillColor(kAzure-9);  h1_HT_jet_csG       -> SetLineColor(kAzure+3); 
  h1_deltaphiggjj_csG -> SetFillColor(kAzure-9);  h1_deltaphiggjj_csG -> SetLineColor(kAzure+3); 
  h1_btag_csG         -> SetFillColor(kAzure-9);  h1_btag_csG         -> SetLineColor(kAzure+3); 

  h1_ptGamma1_csJ -> SetFillColor(kGreen-9);  h1_ptGamma1_csJ -> SetLineColor(kGreen+3); 
  h1_ptGamma2_csJ -> SetFillColor(kGreen-9);  h1_ptGamma2_csJ -> SetLineColor(kGreen+3); 
  h1_ptJet1_csJ   -> SetFillColor(kGreen-9);  h1_ptJet1_csJ   -> SetLineColor(kGreen+3); 
  h1_ptJet2_csJ   -> SetFillColor(kGreen-9);  h1_ptJet2_csJ   -> SetLineColor(kGreen+3);   
  h1_HT_jet_csJ   -> SetFillColor(kGreen-9);  h1_HT_jet_csJ   -> SetLineColor(kGreen+3);   
  h1_deltaphiggjj_csJ -> SetFillColor(kGreen-9);  h1_deltaphiggjj_csJ -> SetLineColor(kGreen+3);   
  h1_btag_csJ     -> SetFillColor(kGreen-9);  h1_btag_csJ     -> SetLineColor(kGreen+3); 

  // ------------------------------------------------------------------------------  
  // mgg  
  tree_data -> Draw("massggnewvtx>>h1_mgg_data(20,100.,180.)", "weight*(massggnewvtx>100 && massggnewvtx<180 && (mjj>150 || mjj<100))");
  tree_cs   -> Draw("massggnewvtx>>h1_mgg_cs(20,100.,180.)",   "weight*(massggnewvtx>100 && massggnewvtx<180)");
  tree_cswG -> Draw("massggnewvtx>>h1_mgg_csw1(20,100.,180.)", "weight*pt_eta_scaled_2D_weight_data*(massggnewvtx>100 && massggnewvtx<180)");
  tree_cswG -> Draw("massggnewvtx>>h1_mgg_csw2(20,100.,180.)", "weight*pt_scaled_2D_weight_data*(massggnewvtx>100 && massggnewvtx<180)");
  tree_cswG -> Draw("massggnewvtx>>h1_mgg_csw3(20,100.,180.)", "weight*eta_scaled_2D_weight_data*(massggnewvtx>100 && massggnewvtx<180)");
  h1_mgg_data -> Sumw2();
  h1_mgg_cs   -> Sumw2();
  h1_mgg_csw1 -> Sumw2();
  h1_mgg_csw2 -> Sumw2();
  h1_mgg_csw3 -> Sumw2();
  
  h1_mgg_data -> GetXaxis()->SetTitle("m_{#gamma #gamma} [GeV]");
  h1_mgg_cs   -> GetXaxis()->SetTitle("m_{#gamma #gamma} [GeV]");
  h1_mgg_csw1 -> GetXaxis()->SetTitle("m_{#gamma #gamma} [GeV]");
  h1_mgg_csw2 -> GetXaxis()->SetTitle("m_{#gamma #gamma} [GeV]");
  h1_mgg_csw3 -> GetXaxis()->SetTitle("m_{#gamma #gamma} [GeV]");
  h1_mgg_data -> SetMarkerStyle(20); 
  h1_mgg_data -> SetMarkerSize(0.8); 
  h1_mgg_cs   -> SetFillColor(kAzure-9); h1_mgg_cs   -> SetLineColor(kAzure+3); 
  h1_mgg_csw1 -> SetFillColor(kAzure-9); h1_mgg_csw1 -> SetLineColor(kAzure+3); 
  h1_mgg_csw2 -> SetFillColor(kAzure-9); h1_mgg_csw2 -> SetLineColor(kAzure+3); 
  h1_mgg_csw3 -> SetFillColor(kAzure-9); h1_mgg_csw3 -> SetLineColor(kAzure+3); 


  // ------------------------------------------------------------------------------
  // mjj  
  tree_data -> Draw("mjj>>h1_mjj_data(30,0.,300.)", "weight*(massggnewvtx>100 && massggnewvtx<180 && (massggnewvtx>135 || massggnewvtx<115))");
  tree_cs   -> Draw("mjj>>h1_mjj_cs(30,0.,300.)",   "weight*(massggnewvtx>100 && massggnewvtx<180)");
  tree_cswG -> Draw("mjj>>h1_mjj_csw1(30,0.,300.)", "weight*pt_eta_scaled_2D_weight_data*(massggnewvtx>100 && massggnewvtx<180)");
  tree_cswG -> Draw("mjj>>h1_mjj_csw2(30,0.,300.)", "weight*pt_scaled_2D_weight_data*(massggnewvtx>100 && massggnewvtx<180)");
  tree_cswG -> Draw("mjj>>h1_mjj_csw3(30,0.,300.)", "weight*eta_scaled_2D_weight_data*(massggnewvtx>100 && massggnewvtx<180)");
  h1_mjj_data -> Sumw2();
  h1_mjj_cs   -> Sumw2();
  h1_mjj_csw1 -> Sumw2();
  h1_mjj_csw2 -> Sumw2();
  h1_mjj_csw3 -> Sumw2();

  h1_mjj_data -> GetXaxis()->SetTitle("m_{#gamma #gamma} [GeV]");
  h1_mjj_cs   -> GetXaxis()->SetTitle("m_{#gamma #gamma} [GeV]");
  h1_mjj_csw1 -> GetXaxis()->SetTitle("m_{#gamma #gamma} [GeV]");
  h1_mjj_csw2 -> GetXaxis()->SetTitle("m_{#gamma #gamma} [GeV]");
  h1_mjj_csw3 -> GetXaxis()->SetTitle("m_{#gamma #gamma} [GeV]");
  h1_mjj_data -> SetMarkerStyle(20); 
  h1_mjj_data -> SetMarkerSize(0.8); 
  h1_mjj_cs   -> SetFillColor(kAzure-9); h1_mjj_cs   -> SetLineColor(kAzure+3); 
  h1_mjj_csw1 -> SetFillColor(kAzure-9); h1_mjj_csw1 -> SetLineColor(kAzure+3); 
  h1_mjj_csw2 -> SetFillColor(kAzure-9); h1_mjj_csw2 -> SetLineColor(kAzure+3); 
  h1_mjj_csw3 -> SetFillColor(kAzure-9); h1_mjj_csw3 -> SetLineColor(kAzure+3); 


  // ------------------------------------------------------------------------------
  // mggjj  
  tree_data -> Draw("mggjj>>h1_mggjj_data(30,100.,600.)", "weight*(massggnewvtx>100 && massggnewvtx<180 && (massggnewvtx>135 || massggnewvtx<115))");
  tree_cs   -> Draw("mggjj>>h1_mggjj_cs(30,100.,600.)",   "weight*(massggnewvtx>100 && massggnewvtx<180)");
  tree_cswG -> Draw("mggjj>>h1_mggjj_csw1(30,100.,600.)", "weight*pt_eta_scaled_2D_weight_data*(massggnewvtx>100 && massggnewvtx<180)");
  tree_cswG -> Draw("mggjj>>h1_mggjj_csw2(30,100.,600.)", "weight*pt_scaled_2D_weight_data*(massggnewvtx>100 && massggnewvtx<180)");
  tree_cswG -> Draw("mggjj>>h1_mggjj_csw3(30,100.,600.)", "weight*eta_scaled_2D_weight_data*(massggnewvtx>100 && massggnewvtx<180)");
  h1_mggjj_data -> Sumw2();
  h1_mggjj_cs   -> Sumw2();
  h1_mggjj_csw1 -> Sumw2();
  h1_mggjj_csw2 -> Sumw2();
  h1_mggjj_csw3 -> Sumw2();

  h1_mggjj_data -> GetXaxis()->SetTitle("m_{jj #gamma #gamma} [GeV]");
  h1_mggjj_cs   -> GetXaxis()->SetTitle("m_{jj #gamma #gamma} [GeV]");
  h1_mggjj_csw1 -> GetXaxis()->SetTitle("m_{jj #gamma #gamma} [GeV]");
  h1_mggjj_csw2 -> GetXaxis()->SetTitle("m_{jj #gamma #gamma} [GeV]");
  h1_mggjj_csw3 -> GetXaxis()->SetTitle("m_{jj #gamma #gamma} [GeV]");
  h1_mggjj_data -> SetMarkerStyle(20); 
  h1_mggjj_data -> SetMarkerSize(0.8); 
  h1_mggjj_cs   -> SetFillColor(kAzure-9);   h1_mggjj_cs   -> SetLineColor(kAzure+3); 
  h1_mggjj_csw1 -> SetFillColor(kAzure-9);   h1_mggjj_csw1 -> SetLineColor(kAzure+3); 
  h1_mggjj_csw2 -> SetFillColor(kAzure-9);   h1_mggjj_csw2 -> SetLineColor(kAzure+3); 
  h1_mggjj_csw3 -> SetFillColor(kAzure-9);   h1_mggjj_csw3 -> SetLineColor(kAzure+3); 


  // ------------------------------------------------------------------------------
  // mjj per cat with the chosen rescaling 
  tree_data -> Draw("mjj>>h1_mjj_data_1(20,0.,300.)", "weight*(massggnewvtx>100 && massggnewvtx<180 && (massggnewvtx>135 || massggnewvtx<115) && btagCategory==1)");
  tree_data -> Draw("mjj>>h1_mjj_data_2(20,0.,300.)", "weight*(massggnewvtx>100 && massggnewvtx<180 && (massggnewvtx>135 || massggnewvtx<115) && btagCategory==2)");
  tree_cswG -> Draw("mjj>>h1_mjj_cs_1(20,0.,300.)",   "weight*pt_scaled_2D_weight_data*(massggnewvtx>100 && massggnewvtx<180 && btagCategory==1)");
  tree_cswG -> Draw("mjj>>h1_mjj_cs_2(20,0.,300.)",   "weight*pt_scaled_2D_weight_data*(massggnewvtx>100 && massggnewvtx<180 && btagCategory==2)");

  h1_mjj_data_1 -> Sumw2();
  h1_mjj_data_2 -> Sumw2();
  h1_mjj_cs_1   -> Sumw2();
  h1_mjj_cs_2   -> Sumw2();

  h1_mjj_data_1 -> GetXaxis()->SetTitle("m_{#gamma #gamma} [GeV]");
  h1_mjj_data_2 -> GetXaxis()->SetTitle("m_{#gamma #gamma} [GeV]");
  h1_mjj_cs_1   -> GetXaxis()->SetTitle("m_{#gamma #gamma} [GeV]");
  h1_mjj_cs_2   -> GetXaxis()->SetTitle("m_{#gamma #gamma} [GeV]");

  h1_mjj_data_1 -> SetMarkerStyle(20); 
  h1_mjj_data_2 -> SetMarkerStyle(20); 
  h1_mjj_data_1 -> SetMarkerSize(0.8); 
  h1_mjj_data_2 -> SetMarkerSize(0.8); 
  h1_mjj_cs_1   -> SetFillColor(kAzure-9); h1_mjj_cs_1   -> SetLineColor(kAzure+3); 
  h1_mjj_cs_2   -> SetFillColor(kAzure-9); h1_mjj_cs_2   -> SetLineColor(kAzure+3); 


  // ------------------------------------------------------------------------------
  // mggjj per cat with the chosen rescaling 
  tree_data -> Draw("mggjj>>h1_mggjj_data_1(20,100.,600.)", "weight*(massggnewvtx>100 && massggnewvtx<180 && (massggnewvtx>135 || massggnewvtx<115) && btagCategory==1)");
  tree_data -> Draw("mggjj>>h1_mggjj_data_2(20,100.,600.)", "weight*(massggnewvtx>100 && massggnewvtx<180 && (massggnewvtx>135 || massggnewvtx<115) && btagCategory==2)");
  tree_cswG -> Draw("mggjj>>h1_mggjj_cs_1(20,100.,600.)",   "weight*pt_scaled_2D_weight_data*(massggnewvtx>100 && massggnewvtx<180 && btagCategory==1)");
  tree_cswG -> Draw("mggjj>>h1_mggjj_cs_2(20,100.,600.)",   "weight*pt_scaled_2D_weight_data*(massggnewvtx>100 && massggnewvtx<180 && btagCategory==2)");

  h1_mggjj_data_1 -> Sumw2();
  h1_mggjj_data_2 -> Sumw2();
  h1_mggjj_cs_1   -> Sumw2();
  h1_mggjj_cs_2   -> Sumw2();

  h1_mggjj_data_1 -> GetXaxis()->SetTitle("m_{jj #gamma #gamma} [GeV]");
  h1_mggjj_data_2 -> GetXaxis()->SetTitle("m_{jj #gamma #gamma} [GeV]");
  h1_mggjj_cs_1   -> GetXaxis()->SetTitle("m_{jj #gamma #gamma} [GeV]");
  h1_mggjj_cs_2   -> GetXaxis()->SetTitle("m_{jj #gamma #gamma} [GeV]");

  h1_mggjj_data_1 -> SetMarkerStyle(20); 
  h1_mggjj_data_2 -> SetMarkerStyle(20); 
  h1_mggjj_data_1 -> SetMarkerSize(0.8); 
  h1_mggjj_data_2 -> SetMarkerSize(0.8); 
  h1_mggjj_cs_1   -> SetFillColor(kAzure-9); h1_mggjj_cs_1   -> SetLineColor(kAzure+3); 
  h1_mggjj_cs_2   -> SetFillColor(kAzure-9); h1_mggjj_cs_2   -> SetLineColor(kAzure+3); 


  // -----------------------------------------------
  TLegend *leg;
  leg = new TLegend(0.6,0.6,0.85,0.85);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.05);
  leg->SetFillColor(0);
  leg->AddEntry(h1_ptGamma1_cs,   "no reweight", "f");
  leg->AddEntry(h1_ptGamma1_csG,  "pT(g1) x pT(g2) weight", "f");
  leg->AddEntry(h1_ptGamma1_csJ,  "pT(j1) x pT(j2) weight", "f");

  TLegend *leg2;
  leg2 = new TLegend(0.5,0.5,0.75,0.75);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.05);
  leg2->SetFillColor(0);
  leg2->AddEntry(h1_mgg_cs, "no reweight", "f");

  TLegend *leg3;
  leg3 = new TLegend(0.5,0.5,0.75,0.75);
  leg3->SetFillStyle(0);
  leg3->SetBorderSize(0);
  leg3->SetTextSize(0.05);
  leg3->SetFillColor(0);
  leg3->AddEntry(h1_mgg_csw1,  "gammas - pT x#eta weight", "f");

  TLegend *leg4;
  leg4 = new TLegend(0.5,0.5,0.75,0.75);
  leg4->SetFillStyle(0);
  leg4->SetBorderSize(0);
  leg4->SetTextSize(0.05);
  leg4->SetFillColor(0);
  leg4->AddEntry(h1_mgg_csw2,  "pT(g1) x pT(g2) weight", "f");

  TLegend *leg5;
  leg5 = new TLegend(0.5,0.5,0.75,0.75);
  leg5->SetFillStyle(0);
  leg5->SetBorderSize(0);
  leg5->SetTextSize(0.05);
  leg5->SetFillColor(0);
  leg5->AddEntry(h1_mgg_csw3,  "#eta (g1) x #eta (g2) weight", "f");  

  TLegend *leg6;
  leg6 = new TLegend(0.6,0.6,0.85,0.85);
  leg6->SetFillStyle(0);
  leg6->SetBorderSize(0);
  leg6->SetTextSize(0.05);
  leg6->SetFillColor(0);
  leg6->AddEntry(h1_mjj_cs_1,  "1btag", "f");

  TLegend *leg7;
  leg7 = new TLegend(0.6,0.6,0.85,0.85);
  leg7->SetFillStyle(0);
  leg7->SetBorderSize(0);
  leg7->SetTextSize(0.05);
  leg7->SetFillColor(0);
  leg7->AddEntry(h1_mjj_cs_2,  "2btag", "f");

  // plots
  TCanvas c1("c1","gamma1 pT",1);  
  c1->Divide(3,1);
  c1->cd(1);  h1_ptGamma1_cs  -> DrawNormalized("hist"); h1_ptGamma1_data -> DrawNormalized("same pE");
  c1->cd(2);  h1_ptGamma1_csG -> DrawNormalized("hist"); h1_ptGamma1_data -> DrawNormalized("same pE");
  c1->cd(3);  h1_ptGamma1_csJ -> DrawNormalized("hist"); h1_ptGamma1_data -> DrawNormalized("same pE");
  leg->Draw();
  c1->SaveAs("gamma1Pt.png");
  c1->SaveAs("gamma1Pt.pdf");
  c1->SaveAs("gamma1Pt.eps");

  TCanvas c1a("c1a","gamma2 pT",1);  
  c1a->Divide(3,1);
  c1a->cd(1);  h1_ptGamma2_cs  -> DrawNormalized("hist"); h1_ptGamma2_data -> DrawNormalized("same pE");
  c1a->cd(2);  h1_ptGamma2_csG -> DrawNormalized("hist"); h1_ptGamma2_data -> DrawNormalized("same pE");
  c1a->cd(3);  h1_ptGamma2_csJ -> DrawNormalized("hist"); h1_ptGamma2_data -> DrawNormalized("same pE");
  leg->Draw();
  c1a->SaveAs("gamma2Pt.png");
  c1a->SaveAs("gamma2Pt.pdf");
  c1a->SaveAs("gamma2Pt.eps");

  TCanvas c2("c2","jet1 pT",1);  
  c2->Divide(3,1);
  c2->cd(1);  h1_ptJet1_cs  -> DrawNormalized("hist"); h1_ptJet1_data -> DrawNormalized("same pE");
  c2->cd(2);  h1_ptJet1_csG -> DrawNormalized("hist"); h1_ptJet1_data -> DrawNormalized("same pE");
  c2->cd(3);  h1_ptJet1_csJ -> DrawNormalized("hist"); h1_ptJet1_data -> DrawNormalized("same pE");
  leg->Draw();
  c2->SaveAs("jet1Pt.png");
  c2->SaveAs("jet1Pt.pdf");
  c2->SaveAs("jet1Pt.eps");

  TCanvas c2a("c2a","jet2 pT",1);  
  c2a->Divide(3,1);
  c2a->cd(1);  h1_ptJet2_cs  -> DrawNormalized("hist"); h1_ptJet2_data -> DrawNormalized("same pE");
  c2a->cd(2);  h1_ptJet2_csG -> DrawNormalized("hist"); h1_ptJet2_data -> DrawNormalized("same pE");
  c2a->cd(3);  h1_ptJet2_csJ -> DrawNormalized("hist"); h1_ptJet2_data -> DrawNormalized("same pE");
  leg->Draw();
  c2a->SaveAs("jet2Pt.png");
  c2a->SaveAs("jet2Pt.pdf");
  c2a->SaveAs("jet2Pt.eps");

  TCanvas c3("c3","HT_jet",1);
  c3->Divide(3,1);
  c3->cd(1); h1_HT_jet_cs->DrawNormalized("hist");   h1_HT_jet_data->DrawNormalized("same pE");
  c3->cd(2); h1_HT_jet_csG->DrawNormalized("hist");  h1_HT_jet_data->DrawNormalized("same pE");
  c3->cd(3); h1_HT_jet_csJ->DrawNormalized("hist");  h1_HT_jet_data->DrawNormalized("same pE");
  leg->Draw();
  c3->SaveAs("HT_jet.png");
  c3->SaveAs("HT_jet.pdf");
  c3->SaveAs("HT_jet.eps");

  TCanvas c4("c4","deltaphiggjj",1);
  c4->Divide(3,1);
  c4->cd(1); h1_deltaphiggjj_cs->DrawNormalized("hist");  h1_deltaphiggjj_data->DrawNormalized("same pE");
  c4->cd(2); h1_deltaphiggjj_csG->DrawNormalized("hist"); h1_deltaphiggjj_data->DrawNormalized("same pE");
  c4->cd(3); h1_deltaphiggjj_csJ->DrawNormalized("hist"); h1_deltaphiggjj_data->DrawNormalized("same pE");
  leg->Draw();
  c4->SaveAs("deltaphiggjj.png");
  c4->SaveAs("deltaphiggjj.pdf");
  c4->SaveAs("deltaphiggjj.eps");

  TCanvas c4b("c4b","btag",1);
  c4b->Divide(3,1);
  c4b->cd(1); h1_btag_cs->DrawNormalized("hist");  h1_btag_data->DrawNormalized("same pE");
  c4b->cd(2); h1_btag_csG->DrawNormalized("hist"); h1_btag_data->DrawNormalized("same pE");
  c4b->cd(3); h1_btag_csJ->DrawNormalized("hist"); h1_btag_data->DrawNormalized("same pE");
  leg->Draw();
  c4b->SaveAs("btag.png");
  c4b->SaveAs("btag.pdf");
  c4b->SaveAs("btag.eps");

  TCanvas c5("c5","mgg",1);
  c5->Divide(2,2);
  c5->cd(1); h1_mgg_cs->DrawNormalized("hist");   h1_mgg_data->DrawNormalized("same pE"); leg2->Draw();
  c5->cd(2); h1_mgg_csw1->DrawNormalized("hist"); h1_mgg_data->DrawNormalized("same pE"); leg3->Draw();
  c5->cd(3); h1_mgg_csw2->DrawNormalized("hist"); h1_mgg_data->DrawNormalized("same pE"); leg4->Draw();
  c5->cd(4); h1_mgg_csw3->DrawNormalized("hist"); h1_mgg_data->DrawNormalized("same pE"); leg5->Draw();
  c5->SaveAs("mgg.png");
  c5->SaveAs("mgg.pdf");
  c5->SaveAs("mgg.eps");

  TCanvas c6("c6","mjj",1);
  c6->Divide(2,2);
  c6->cd(1); h1_mjj_cs->DrawNormalized("hist");   h1_mjj_data->DrawNormalized("same pE"); leg2->Draw();
  c6->cd(2); h1_mjj_csw1->DrawNormalized("hist"); h1_mjj_data->DrawNormalized("same pE"); leg3->Draw();
  c6->cd(3); h1_mjj_csw2->DrawNormalized("hist"); h1_mjj_data->DrawNormalized("same pE"); leg4->Draw();
  c6->cd(4); h1_mjj_csw3->DrawNormalized("hist"); h1_mjj_data->DrawNormalized("same pE"); leg5->Draw();
  c6->SaveAs("mjj.png");
  c6->SaveAs("mjj.pdf");
  c6->SaveAs("mjj.eps");

  TCanvas c6a("c6a","mjj",1);
  c6a->Divide(2,1);
  c6a->cd(1); h1_mjj_cs_1->DrawNormalized("hist");   h1_mjj_data_1->DrawNormalized("same pE"); leg6->Draw();
  c6a->cd(2); h1_mjj_cs_2->DrawNormalized("hist");   h1_mjj_data_2->DrawNormalized("same pE"); leg7->Draw();
  c6a->SaveAs("mjj_perCat.png");
  c6a->SaveAs("mjj_perCat.pdf");
  c6a->SaveAs("mjj_perCat.eps");

  TCanvas c7("c7","mggjj",1);
  c7->Divide(2,2);
  c7->cd(1); h1_mggjj_cs->DrawNormalized("hist");   h1_mggjj_data->DrawNormalized("same pE"); leg2->Draw();
  c7->cd(2); h1_mggjj_csw1->DrawNormalized("hist"); h1_mggjj_data->DrawNormalized("same pE"); leg3->Draw();
  c7->cd(3); h1_mggjj_csw2->DrawNormalized("hist"); h1_mggjj_data->DrawNormalized("same pE"); leg4->Draw();
  c7->cd(4); h1_mggjj_csw3->DrawNormalized("hist"); h1_mggjj_data->DrawNormalized("same pE"); leg5->Draw();
  c7->SaveAs("mggjj.png");
  c7->SaveAs("mggjj.pdf");
  c7->SaveAs("mggjj.eps");

  TCanvas c7a("c7a","mggjj",1);
  c7a->Divide(2,1);
  c7a->cd(1); h1_mggjj_cs_1->DrawNormalized("hist");   h1_mggjj_data_1->DrawNormalized("same pE"); leg6->Draw();
  c7a->cd(2); h1_mggjj_cs_2->DrawNormalized("hist");   h1_mggjj_data_2->DrawNormalized("same pE"); leg7->Draw();
  c7a->SaveAs("mggjj_perCat.png");
  c7a->SaveAs("mggjj_perCat.pdf");
  c7a->SaveAs("mggjj_perCat.eps");
}
