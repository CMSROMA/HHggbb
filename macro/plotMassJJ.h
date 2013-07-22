//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Apr 17 15:30:38 2013 by ROOT version 5.32/00
// from TTree myTrees/
// found on file: finalizedTrees_Radion/Radion_Radion_M-300_madgraph_default_CSV.root
//////////////////////////////////////////////////////////

#ifndef plotMassJJ_h
#define plotMassJJ_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class plotMassJJ {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Float_t         mjj_gen;
   Float_t         mjj_pt;
   Float_t         mjj_maxptjj;
   Float_t         mjj_maxptmjj;
   Float_t         mjj_mgg;
   Float_t         mjj_btag_pt;
   Float_t         mjj_btag_ptCSV;
   Float_t         mjj_btag_ptjj;
   Float_t         mjj_btag_ptmjj;
   Float_t         mjj_btag_mgg;
   Int_t           btagCategory;
   Int_t           ngoodJets;
   Float_t         dMmin;
   Float_t         dMmin_btag;

   // List of branches
   TBranch        *b_invMassJJ_gen;   //!
   TBranch        *b_invMassJJ_pt;   //!
   TBranch        *b_invMassJJ_maxptjj;   //!
   TBranch        *b_invMassJJ_maxptmjj;   //!
   TBranch        *b_invMassJJ_mgg;   //!
   TBranch        *b_invMassJJ_btag_pt;   //!
   TBranch        *b_invMassJJ_btag_ptCSV;   //!
   TBranch        *b_invMassJJ_btag_ptjj;   //!
   TBranch        *b_invMassJJ_btag_ptmjj;   //!
   TBranch        *b_invMassJJ_btag_mgg;   //!
   TBranch        *b_btagCategory_t;   //!
   TBranch        *b_ngoodJets_t;  //!
   TBranch        *b_dMmin_t;
   TBranch        *b_dMmin_btag_t;

   plotMassJJ(TTree *tree=0);
   virtual ~plotMassJJ();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef plotMassJJ_cxx
plotMassJJ::plotMassJJ(TTree *tree) : fChain(0) 
{
  if (tree == 0) {                                                            
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../trees_for_jetsChoice/treeMjj_sign300_mediumBtag_conOlivierRegression/Radion_Radion_M-300_regr_default_CSV.root");
    if (!f || !f->IsOpen()) {
      f = new TFile("../trees_for_jetsChoice/treeMjj_sign300_mediumBtag_conOlivierRegression/Radion_Radion_M-300_regr_default_CSV.root");
    }
    f->GetObject("myTrees",tree);
    
  }
  Init(tree);
}

plotMassJJ::~plotMassJJ()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t plotMassJJ::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}

Long64_t plotMassJJ::LoadTree(Long64_t entry)
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

void plotMassJJ::Init(TTree *tree)
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
  
  fChain->SetBranchAddress("mjj_gen", &mjj_gen, &b_invMassJJ_gen);
  fChain->SetBranchAddress("mjj_pt", &mjj_pt, &b_invMassJJ_pt);
  fChain->SetBranchAddress("mjj_maxptjj", &mjj_maxptjj, &b_invMassJJ_maxptjj);
  fChain->SetBranchAddress("mjj_maxptmjj", &mjj_maxptmjj, &b_invMassJJ_maxptmjj);
  fChain->SetBranchAddress("mjj_mgg", &mjj_mgg, &b_invMassJJ_mgg);
  fChain->SetBranchAddress("mjj_btag_pt", &mjj_btag_pt, &b_invMassJJ_btag_pt);
  fChain->SetBranchAddress("mjj_btag_ptCSV", &mjj_btag_ptCSV, &b_invMassJJ_btag_ptCSV);
  fChain->SetBranchAddress("mjj_btag_ptjj", &mjj_btag_ptjj, &b_invMassJJ_btag_ptjj);
  fChain->SetBranchAddress("mjj_btag_ptmjj", &mjj_btag_ptmjj, &b_invMassJJ_btag_ptmjj);
  fChain->SetBranchAddress("mjj_btag_mgg", &mjj_btag_mgg, &b_invMassJJ_btag_mgg);
  fChain->SetBranchAddress("btagCategory", &btagCategory, &b_btagCategory_t);
  fChain->SetBranchAddress("ngoodJets", &ngoodJets, &b_ngoodJets_t);
  fChain->SetBranchAddress("dMmin", &dMmin, &b_dMmin_t);
  fChain->SetBranchAddress("dMmin_btag", &dMmin_btag, &b_dMmin_btag_t);
  Notify();
}

Bool_t plotMassJJ::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.
  
  return kTRUE;
}

void plotMassJJ::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}

Int_t plotMassJJ::Cut(Long64_t entry)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}
#endif // #ifdef plotMassJJ_cxx
