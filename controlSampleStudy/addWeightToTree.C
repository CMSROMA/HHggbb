#define addWeightToTree_cxx
#include "addWeightToTree.h"
#include <TH2.h>
#include "TFile.h"
#include "TTree.h"
#include <TChain.h>
#include <iostream>
#include <fstream>
#include <math.h>

void addWeightToTree::Loop()
{
  if (fChain == 0) return;
  
  // pt and eta weight files
  std::string ptweight2DFileName_data  = "weights_withRegression_withKinFit_noPreselCut/scales_2D_pt_data_4GeVbinning.root";
  std::string etaweight2DFileName_data = "weights_withRegression_withKinFit_noPreselCut/scales_2D_eta_data_01binning.root";
  // std::string ptweight2DFileName_data  = "weights_withRegression_withKinFit_noPreselCut/scales_2D_ptJ_data_4GeVbinning.root";
  // std::string etaweight2DFileName_data = "weights_withRegression_withKinFit_noPreselCut/scales_2D_etaJ_data_01binning.root";
  
  TFile* ptweight2DFile_data  = TFile::Open(ptweight2DFileName_data.c_str());
  TFile* etaweight2DFile_data = TFile::Open(etaweight2DFileName_data.c_str()); 
  
  TH2F* h2_ptweight_data  = (TH2F*)ptweight2DFile_data->Get("h2D_pt_data");
  TH2F* h2_etaweight_data = (TH2F*)etaweight2DFile_data->Get("h2D_eta_data");
  // TH2F* h2_ptweight_data  = (TH2F*)ptweight2DFile_data->Get("h2D_ptJ_data");
  // TH2F* h2_etaweight_data = (TH2F*)etaweight2DFile_data->Get("h2D_etaJ_data");

  // output tree, with modified weight
  TFile* outFile=TFile::Open("csWithWeight.root","recreate");
  TTree* treeWithWeights= new TTree();
  treeWithWeights->SetName("myTrees_withWeight");
  createBranches(treeWithWeights);
  
  Long64_t nentries = fChain->GetEntriesFast();
  
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    
    if(jentry%500==0)cout<<jentry<<endl;
    
    // reweighting for cs
    double ptweight2D_data=h2_ptweight_data->GetBinContent(h2_ptweight_data->GetXaxis()->FindBin(ptPhot2),h2_ptweight_data->GetYaxis()->FindBin(ptPhot1));
    double etaweight2D_data=h2_etaweight_data->GetBinContent(h2_etaweight_data->GetXaxis()->FindBin(etaPhot2),h2_etaweight_data->GetYaxis()->FindBin(etaPhot1));
    // double ptweight2D_data=h2_ptweight_data->GetBinContent(h2_ptweight_data->GetXaxis()->FindBin(ptCorrJet2),h2_ptweight_data->GetYaxis()->FindBin(ptCorrJet1));
    // double etaweight2D_data=h2_etaweight_data->GetBinContent(h2_etaweight_data->GetXaxis()->FindBin(etaJet2),h2_etaweight_data->GetYaxis()->FindBin(etaJet1));
    
    if (ptweight2D_data*etaweight2D_data<5.){   //skip few big weights in pt due to low stat
      pt_scaled_2D_weight_data_t     = ptweight2D_data;
      eta_scaled_2D_weight_data_t    = etaweight2D_data;
      pt_eta_scaled_2D_weight_data_t = ptweight2D_data*etaweight2D_data;
    } else {
      pt_scaled_2D_weight_data_t     = 1.;
      eta_scaled_2D_weight_data_t    = 1.;
      pt_eta_scaled_2D_weight_data_t = 1.;
    }
    
    treeWithWeights->Fill();
  }
  
  treeWithWeights->Write();
  outFile->Write();
  outFile->Close();
}
