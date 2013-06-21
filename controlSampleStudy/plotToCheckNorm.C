{
  gStyle->SetOptStat(0);
  
  TFile* dataFile=TFile::Open("../finalizedTrees_Radion_presel/Radion_Data2012_default_CSV.root");
  TTree* tree_data=(TTree*)dataFile->Get("myTrees");

  TFile* csFileWeightG=TFile::Open("treesFromCS_presel_withWeights/csWithWeightFromGammas.root");
  TTree* tree_cswG=(TTree*)csFileWeightG->Get("myTrees_withWeight");

  // ------------------------------------------------------------------------------  
  cout << endl;
  tree_data -> Draw("massggnewvtx>>h1_check_data(20,100.,180.)","(massggnewvtx>100 && massggnewvtx<180)");
  tree_cswG -> Draw("massggnewvtx>>h1_check_cs(20,100.,180.)",  "(pt_scaled_2D_weight_data)*(massggnewvtx>100 && massggnewvtx<180)");
  cout << "dati: " << h1_check_data -> Integral() << endl;
  cout << "CS: "   << h1_check_cs   -> Integral() << endl;
  cout << "weight CS with " << (float)(h1_check_data->Integral())/(float)(h1_check_cs->Integral())  << endl;
  cout << endl;

  // ------------------------------------------------------------------------------  
  // mgg  
  tree_data -> Draw("massggnewvtx>>h1_mgg_data(20,100.,180.)","(massggnewvtx>100 && massggnewvtx<180)");
  tree_cswG -> Draw("massggnewvtx>>h1_mgg_csR(20,100.,180.)", "(0.70*pt_scaled_2D_weight_data)*(massggnewvtx>100 && massggnewvtx<180)");
  h1_mgg_data -> Sumw2();
  h1_mgg_csR  -> Sumw2();

  h1_mgg_data -> GetXaxis()->SetTitle("m_{#gamma #gamma}");
  h1_mgg_csR  -> GetXaxis()->SetTitle("m_{#gamma #gamma}");
  h1_mgg_data -> SetMarkerStyle(20); 
  h1_mgg_data -> SetMarkerSize(0.8); 
  h1_mgg_csR  -> SetFillColor(kAzure-9); 
  h1_mgg_csR  -> SetLineColor(kAzure+3); 


  // ------------------------------------------------------------------------------
  // mjj  
  tree_data -> Draw("mjj>>h1_mjj_data(30,0.,300.)", "(massggnewvtx>100 && massggnewvtx<180)");
  tree_cswG -> Draw("mjj>>h1_mjj_cs(30,0.,300.)",   "(0.70*pt_scaled_2D_weight_data)*(massggnewvtx>100 && massggnewvtx<180)");
  h1_mjj_data -> Sumw2();
  h1_mjj_cs   -> Sumw2();

  h1_mjj_data -> GetXaxis()->SetTitle("m_{#gamma #gamma}");
  h1_mjj_cs   -> GetXaxis()->SetTitle("m_{#gamma #gamma}");
  h1_mjj_data -> SetMarkerStyle(20); 
  h1_mjj_data -> SetMarkerSize(0.8); 
  h1_mjj_cs   -> SetFillColor(kAzure-9); 
  h1_mjj_cs   -> SetLineColor(kAzure+3); 

  // ------------------------------------------------------------------------------
  // mggjj  
  tree_data -> Draw("mggjj>>h1_mggjj_data(30,100.,600.)", "(massggnewvtx>100 && massggnewvtx<180)");
  tree_cswG -> Draw("mggjj>>h1_mggjj_cs(30,100.,600.)",   "(0.70*pt_scaled_2D_weight_data)*(massggnewvtx>100 && massggnewvtx<180)");
  h1_mggjj_data -> Sumw2();
  h1_mggjj_cs   -> Sumw2();

  h1_mggjj_data -> GetXaxis()->SetTitle("m_{jj #gamma #gamma}");
  h1_mggjj_cs   -> GetXaxis()->SetTitle("m_{jj #gamma #gamma}");
  h1_mggjj_data -> SetMarkerStyle(20); 
  h1_mggjj_data -> SetMarkerSize(0.8); 
  h1_mggjj_cs   -> SetFillColor(kAzure-9);   
  h1_mggjj_cs   -> SetLineColor(kAzure+3); 

  // ------------------------------------------------------------------------------  
  TCanvas c1("c1","mgg",1);
  h1_mgg_csR->Draw("hist");   
  h1_mgg_data->Draw("same pE"); 
  c1->SaveAs("mgg.png");
  c1->SaveAs("mgg.root");

  TCanvas c2("c2","mjj",1);
  h1_mjj_cs->Draw("hist");   
  h1_mjj_data->Draw("same pE"); 
  c2->SaveAs("mjj.png");
  c2->SaveAs("mjj.root");

  TCanvas c3("c3","mggjj",1);
  h1_mggjj_cs->Draw("hist");   
  h1_mggjj_data->Draw("same pE"); 
  c3->SaveAs("mggjj.png");
  c3->SaveAs("mggjj.root");
}
