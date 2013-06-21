{
  // weights_from_data
  TFile* dataFile=TFile::Open("../finalizedTrees_Radion_presel/Radion_Data2012_default_CSV.root");
  TTree* tree_data=(TTree*)dataFile->Get("myTrees");

  TFile* csFile=TFile::Open("../finalizedTrees_Radion_presel_CS/Radion_Data2012_default_CSV.root");
  TTree* tree_cs=(TTree*)csFile->Get("myTrees");
  
  // pt gamma 
  tree_data->Draw("ptPhot1:ptPhot2>>h2D_pt_data(35,20.,160.,35,20.,160.)","weight*(massggnewvtx>100 && massggnewvtx<180)*(massggnewvtx<115 || massggnewvtx>135)","colz");
  tree_cs->Draw("ptPhot1:ptPhot2>>h2D_pt_cs(35,20.,160.,35,20.,160.)","weight","colz");
  float integral_2D_pt_data = h2D_pt_data->Integral();
  h2D_pt_data->Scale(1./integral_2D_pt_data);
  float integral_2D_pt_cs = h2D_pt_cs->Integral();
  h2D_pt_cs->Scale(1./integral_2D_pt_cs);
  TH2F* cs_norm_pt = h2D_pt_cs;
  h2D_pt_data->Divide(cs_norm_pt);
  h2D_pt_data->GetXaxis()->SetTitle("p_{T} Sublead Photon");
  h2D_pt_data->GetYaxis()->SetTitle("p_{T} Lead Photon");
  h2D_pt_data->SaveAs("scales_2D_pt_data_4GeVbinning.root");
  
  // eta gamma
  tree_data->Draw("etaPhot1:etaPhot2>>h2D_eta_data(30,-3,3,30,-3,3)","weight*(massggnewvtx>100 && massggnewvtx<180)*(massggnewvtx<115 || massggnewvtx>135)","colz");
  tree_cs->Draw("etaPhot1:etaPhot2>>h2D_eta_cs(30,-3,3,30,-3,3)","weight","colz");
  float integral_2D_eta_data = h2D_eta_data->Integral();
  h2D_eta_data->Scale(1./integral_2D_eta_data);
  float integral_2D_eta_cs = h2D_eta_cs->Integral();
  h2D_eta_cs->Scale(1./integral_2D_eta_cs);
  TH2F* cs_norm_eta = h2D_eta_cs;
  h2D_eta_data->Divide(cs_norm_eta);
  h2D_eta_data->GetXaxis()->SetTitle("#eta Sublead Photon");
  h2D_eta_data->GetYaxis()->SetTitle("#eta Lead Photon");
  h2D_eta_data->SaveAs("scales_2D_eta_data_01binning.root");

  // pt jets 
  tree_data->Draw("ptCorrJet1:ptCorrJet2>>h2D_ptJ_data(35,20.,160.,35,20.,160.)","weight*(massggnewvtx>100 && massggnewvtx<180)*(massggnewvtx<115 || massggnewvtx>135)","colz");
  tree_cs->Draw("ptCorrJet1:ptCorrJet2>>h2D_ptJ_cs(35,20.,160.,35,20.,160.)","weight","colz");
  float integral_2D_ptJ_data = h2D_ptJ_data->Integral();
  h2D_ptJ_data->Scale(1./integral_2D_ptJ_data);
  float integral_2D_ptJ_cs = h2D_ptJ_cs->Integral();
  h2D_ptJ_cs->Scale(1./integral_2D_ptJ_cs);
  TH2F* cs_norm_ptJ = h2D_ptJ_cs;
  h2D_ptJ_data->Divide(cs_norm_ptJ);
  h2D_ptJ_data->GetXaxis()->SetTitle("p_{T} Sublead Jet");
  h2D_ptJ_data->GetYaxis()->SetTitle("p_{T} Lead Jet");
  h2D_ptJ_data->SaveAs("scales_2D_ptJ_data_4GeVbinning.root");

  // eta jets
  tree_data->Draw("etaJet1:etaJet2>>h2D_etaJ_data(30,-3,3,30,-3,3)","weight*(massggnewvtx>100 && massggnewvtx<180)*(massggnewvtx<115 || massggnewvtx>135)","colz");
  tree_cs->Draw("etaJet1:etaJet2>>h2D_etaJ_cs(30,-3,3,30,-3,3)","weight","colz");
  float integral_2D_etaJ_data = h2D_etaJ_data->Integral();
  h2D_etaJ_data->Scale(1./integral_2D_etaJ_data);
  float integral_2D_etaJ_cs = h2D_etaJ_cs->Integral();
  h2D_etaJ_cs->Scale(1./integral_2D_etaJ_cs);
  TH2F* cs_norm_etaJ = h2D_etaJ_cs;
  h2D_etaJ_data->Divide(cs_norm_etaJ);
  h2D_etaJ_data->GetXaxis()->SetTitle("#eta Sublead Jet");
  h2D_etaJ_data->GetYaxis()->SetTitle("#eta Lead Jet");
  h2D_etaJ_data->SaveAs("scales_2D_etaJ_data_01binning.root");
}
