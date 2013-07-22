//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Jun  5 16:43:42 2013 by ROOT version 5.32/00
// from TTree myTrees/
// found on file: ../finalizedTrees_Radion_presel_CS/Radion_Data2012_default_CSV.root
//////////////////////////////////////////////////////////

#ifndef addWeightToTree_h
#define addWeightToTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>

using namespace std;

// Fixed size dimensions of array or collections stored in the TTree if any.

class addWeightToTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // output stuffs
   TString optFileName;
   TString outFileName;

   // new weights
   Float_t pt_scaled_2D_weight_data_t;
   Float_t eta_scaled_2D_weight_data_t;
   Float_t pt_eta_scaled_2D_weight_data_t;

   // Declaration of leaf types
   Int_t           run;
   Int_t           event;
   Float_t         massggnewvtx;
   Float_t         ptPhot1;
   Float_t         runptPhot1;
   Float_t         absCosThetaStar;
   Float_t         ptPhot2;
   Float_t         etaPhot1;
   Float_t         etaPhot2;
   Int_t           cicPhot1;
   Int_t           cicPhot2;
   Float_t         r9Phot1;
   Float_t         r9Phot2;
   Float_t         ptgg;
   Float_t         etagg;
   Float_t         absetagg;
   Int_t           njets;
   Float_t         runPtCorrJet1;
   Float_t         ptCorrJet1;
   Float_t         ptCorrJet2;
   Float_t         etaJet1;
   Float_t         etaJet2;
   Float_t         deltaphijj;
   Float_t         deltaetajj;
   Float_t         mjj;
   Float_t         ptjj;
   Float_t         minDeltaR_gb;
   Float_t         etajj;
   Float_t         mggjj;
   Float_t         deltaphiggjj;
   Float_t         deltaetaggjj;
   Float_t         deltaRggjj;
   Int_t           btagCategory;
   Int_t           theCategory;
   Int_t           theGammaCategory;
   Int_t           nbjets_loose;
   Int_t           nbjets_medium;
   Int_t           nbjets_tight;
   Float_t         chiSquareProbH;
   Int_t           nvtx;
   Int_t           vertex;
   Int_t           isBtagJet1;
   Int_t           isBtagJet2;
   Float_t         mjj_kin;
   Float_t         mggjj_kin;
   Float_t         HT_jet;
   Float_t         HT_gjet;
   Float_t         weight;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_massggnewvtx_t;   //!
   TBranch        *b_ptphot1_t;   //!
   TBranch        *b_runptPhot1_t;   //!
   TBranch        *b_absCosThetaStar_t;   //!
   TBranch        *b_ptphot2_t;   //!
   TBranch        *b_etaphot1_t;   //!
   TBranch        *b_etaphot2_t;   //!
   TBranch        *b_cicphot1_t;   //!
   TBranch        *b_cicphot2_t;   //!
   TBranch        *b_r9phot1_t;   //!
   TBranch        *b_r9phot2_t;   //!
   TBranch        *b_ptgg_t;   //!
   TBranch        *b_etagg_t;   //!
   TBranch        *b_absetagg_t;   //!
   TBranch        *b_njets_t;   //!
   TBranch        *b_runptcorrJet1_t;   //!
   TBranch        *b_ptcorrJet1_t;   //!
   TBranch        *b_ptcorrJet2_t;   //!
   TBranch        *b_etajet1_t;   //!
   TBranch        *b_etajet2_t;   //!
   TBranch        *b_deltaphijj_t;   //!
   TBranch        *b_deltaetajj_t;   //!
   TBranch        *b_invmassjet_t;   //!
   TBranch        *b_ptjj_t;   //!
   TBranch        *b_minDeltaR_gb_t;   //!
   TBranch        *b_etajj_t;   //!
   TBranch        *b_massggjj_t;   //!
   TBranch        *b_deltaphijjgg_t;   //!
   TBranch        *b_deltaetajjgg_t;   //!
   TBranch        *b_deltaRjjgg_t;   //!
   TBranch        *b_btagCategory_t;   //!
   TBranch        *b_theCategory_t;   //!
   TBranch        *b_theGammaCategory_t;   //!
   TBranch        *b_nbjets_loose_t;   //!
   TBranch        *b_nbjets_medium_t;   //!
   TBranch        *b_nbjets_tight_t;   //!
   TBranch        *b_chiSquareProbH_t;   //!
   TBranch        *b_nvtx_t;   //!
   TBranch        *b_theVertex_t;   //!
   TBranch        *b_isj1btagged_t;   //!
   TBranch        *b_isj2btagged_t;   //!
   TBranch        *b_mjj_kin_t;   //!
   TBranch        *b_mggjj_kin_t;   //!
   TBranch        *b_HT_jet_t;   //!
   TBranch        *b_HT_gjet_t;   //!
   TBranch        *b_weight_t;   //!

   addWeightToTree(TTree *tree=0);
   virtual ~addWeightToTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   void createBranches(TTree* treeWithWeights);
};

#endif

#ifdef addWeightToTree_cxx
addWeightToTree::addWeightToTree(TTree *tree) : fChain(0) 
{
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../finalizedTrees_m300_noCutPerOptim_conRegression_conKinFit_CS/Radion_DataABCD_regr_default_CSV.root");
    if (!f || !f->IsOpen()) {
      f = new TFile("../finalizedTrees_m300_noCutPerOptim_conRegression_conKinFit_CS/Radion_DataABCD_regr_default_CSV.root");
    }
    f->GetObject("myTrees",tree);
    
   }
   Init(tree);
}

addWeightToTree::~addWeightToTree()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t addWeightToTree::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}

Long64_t addWeightToTree::LoadTree(Long64_t entry)
{
  // Set the environment to read one entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (fChain->GetTreeNumber() != fCurrent) {
    fCurrent = fChain->GetTreeNumber();
    Notify();
  }
  return centry;
}

void addWeightToTree::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).
  
  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);
  
  fChain->SetBranchAddress("run", &run, &b_run);
  fChain->SetBranchAddress("event", &event, &b_event);
  fChain->SetBranchAddress("massggnewvtx", &massggnewvtx, &b_massggnewvtx_t);
  fChain->SetBranchAddress("ptPhot1", &ptPhot1, &b_ptphot1_t);
  fChain->SetBranchAddress("runptPhot1", &runptPhot1, &b_runptPhot1_t);
  fChain->SetBranchAddress("absCosThetaStar", &absCosThetaStar, &b_absCosThetaStar_t);
  fChain->SetBranchAddress("ptPhot2", &ptPhot2, &b_ptphot2_t);
  fChain->SetBranchAddress("etaPhot1", &etaPhot1, &b_etaphot1_t);
  fChain->SetBranchAddress("etaPhot2", &etaPhot2, &b_etaphot2_t);
  fChain->SetBranchAddress("cicPhot1", &cicPhot1, &b_cicphot1_t);
  fChain->SetBranchAddress("cicPhot2", &cicPhot2, &b_cicphot2_t);
  fChain->SetBranchAddress("r9Phot1", &r9Phot1, &b_r9phot1_t);
  fChain->SetBranchAddress("r9Phot2", &r9Phot2, &b_r9phot2_t);
  fChain->SetBranchAddress("ptgg", &ptgg, &b_ptgg_t);
  fChain->SetBranchAddress("etagg", &etagg, &b_etagg_t);
  fChain->SetBranchAddress("absetagg", &absetagg, &b_absetagg_t);
  fChain->SetBranchAddress("njets", &njets, &b_njets_t);
  fChain->SetBranchAddress("runPtCorrJet1", &runPtCorrJet1, &b_runptcorrJet1_t);
  fChain->SetBranchAddress("ptCorrJet1", &ptCorrJet1, &b_ptcorrJet1_t);
  fChain->SetBranchAddress("ptCorrJet2", &ptCorrJet2, &b_ptcorrJet2_t);
  fChain->SetBranchAddress("etaJet1", &etaJet1, &b_etajet1_t);
  fChain->SetBranchAddress("etaJet2", &etaJet2, &b_etajet2_t);
  fChain->SetBranchAddress("deltaphijj", &deltaphijj, &b_deltaphijj_t);
  fChain->SetBranchAddress("deltaetajj", &deltaetajj, &b_deltaetajj_t);
  fChain->SetBranchAddress("mjj", &mjj, &b_invmassjet_t);
  fChain->SetBranchAddress("ptjj", &ptjj, &b_ptjj_t);
  fChain->SetBranchAddress("minDeltaR_gb", &minDeltaR_gb, &b_minDeltaR_gb_t);
  fChain->SetBranchAddress("etajj", &etajj, &b_etajj_t);
  fChain->SetBranchAddress("mggjj", &mggjj, &b_massggjj_t);
  fChain->SetBranchAddress("deltaphiggjj", &deltaphiggjj, &b_deltaphijjgg_t);
  fChain->SetBranchAddress("deltaetaggjj", &deltaetaggjj, &b_deltaetajjgg_t);
  fChain->SetBranchAddress("deltaRggjj", &deltaRggjj, &b_deltaRjjgg_t);
  fChain->SetBranchAddress("btagCategory", &btagCategory, &b_btagCategory_t);
  fChain->SetBranchAddress("theCategory", &theCategory, &b_theCategory_t);
  fChain->SetBranchAddress("theGammaCategory", &theGammaCategory, &b_theGammaCategory_t);
  fChain->SetBranchAddress("nbjets_loose", &nbjets_loose, &b_nbjets_loose_t);
  fChain->SetBranchAddress("nbjets_medium", &nbjets_medium, &b_nbjets_medium_t);
  fChain->SetBranchAddress("nbjets_tight", &nbjets_tight, &b_nbjets_tight_t);
  fChain->SetBranchAddress("chiSquareProbH", &chiSquareProbH, &b_chiSquareProbH_t);
  fChain->SetBranchAddress("nvtx", &nvtx, &b_nvtx_t);
  fChain->SetBranchAddress("vertex", &vertex, &b_theVertex_t);
  fChain->SetBranchAddress("isBtagJet1", &isBtagJet1, &b_isj1btagged_t);
  fChain->SetBranchAddress("isBtagJet2", &isBtagJet2, &b_isj2btagged_t);
  fChain->SetBranchAddress("mjj_kin", &mjj_kin, &b_mjj_kin_t);
  fChain->SetBranchAddress("mggjj_kin", &mggjj_kin, &b_mggjj_kin_t);
  fChain->SetBranchAddress("HT_jet", &HT_jet, &b_HT_jet_t);
  fChain->SetBranchAddress("HT_gjet", &HT_gjet, &b_HT_gjet_t);
  fChain->SetBranchAddress("weight", &weight, &b_weight_t);
  Notify();
}

Bool_t addWeightToTree::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.
  
  return kTRUE;
}

void addWeightToTree::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}

Int_t addWeightToTree::Cut(Long64_t entry)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}

void addWeightToTree::createBranches(TTree* treeWithWeights){

  treeWithWeights->Branch("massggnewvtx", &massggnewvtx, "massggnewvtx/F");
  treeWithWeights->Branch("runptPhot1", &runptPhot1, "runptPhot1/F");
  treeWithWeights->Branch("absCosThetaStar", &absCosThetaStar, "absCosThetaStar/F");
  treeWithWeights->Branch("ptPhot1", &ptPhot1, "ptPhot1/F");
  treeWithWeights->Branch("ptPhot2", &ptPhot2, "ptphot2/F");
  treeWithWeights->Branch("ptgg", &ptgg, "ptgg/F");
  treeWithWeights->Branch("ptjj", &ptjj, "ptjj/F");
  treeWithWeights->Branch("etaPhot1", &etaPhot1, "etaphot1/F");
  treeWithWeights->Branch("etaPhot2", &etaPhot2, "etaphot2/F");
  treeWithWeights->Branch("r9Phot1", &r9Phot1, "r9phot1/F");
  treeWithWeights->Branch("r9Phot2", &r9Phot2, "r9phot2/F");
  treeWithWeights->Branch("runPtCorrJet1", &runPtCorrJet1, "runptcorrJet1/F");
  treeWithWeights->Branch("ptCorrJet1", &ptCorrJet1, "ptCorrJet1/F");
  treeWithWeights->Branch("ptCorrJet2", &ptCorrJet2, "ptCorrJet2/F");
  treeWithWeights->Branch("etaJet1", &etaJet1, "etaJet1/F");
  treeWithWeights->Branch("etaJet2", &etaJet2, "etaJet2/F");
  treeWithWeights->Branch("mjj", &mjj, "mjj/F");
  treeWithWeights->Branch("mggjj", &mggjj, "mggjj/F");
  treeWithWeights->Branch("deltaphiggjj", &deltaphiggjj, "deltaphijjgg/F");
  treeWithWeights->Branch("btagCategory", &btagCategory, "btagCategory/I");
  treeWithWeights->Branch("HT_jet", &HT_jet, "HT_jet/F");
  treeWithWeights->Branch("minDeltaR_gb", &minDeltaR_gb, "minDeltaR_gb/F");
  treeWithWeights->Branch("weight", &weight, "weight/F");
  treeWithWeights->Branch("pt_scaled_2D_weight_data", &pt_scaled_2D_weight_data_t, "pt_scaled_2D_weight_data/F");
  treeWithWeights->Branch("eta_scaled_2D_weight_data", &eta_scaled_2D_weight_data_t, "eta_scaled_2D_weight_data/F");
  treeWithWeights->Branch("pt_eta_scaled_2D_weight_data", &pt_eta_scaled_2D_weight_data_t, "pt_eta_scaled_2D_weight_data/F");
}

#endif // #ifdef addWeightToTree_cxx
