//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Apr 19 10:15:51 2013 by ROOT version 5.32/00
// from TTree Radion_m300_8TeV/Radion_m300_8TeV
// found on file: Radion.root
//////////////////////////////////////////////////////////

#ifndef fillPlot2012_radion_commonNtp_h
#define fillPlot2012_radion_commonNtp_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TRandom3.h>

#include "RedNtpFinalizer_commonNtp.h"

using namespace std;

class fillPlot2012_radion_commonNtp : public RedNtpFinalizer_commonNtp {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           itype;
   Float_t         run;
   Float_t         lumis;
   Float_t         weight;
   Float_t         evweight;
   Float_t         pu_weight;
   Float_t         pu_n;
   Float_t         nvtx;
   Float_t         rho;
   Int_t           category;
   Float_t         ph1_e;
   Float_t         ph2_e;
   Float_t         ph1_pt;
   Float_t         ph2_pt;
   Float_t         ph1_phi;
   Float_t         ph2_phi;
   Float_t         ph1_eta;
   Float_t         ph2_eta;
   Float_t         ph1_r9;
   Float_t         ph2_r9;
   Int_t           ph1_isPrompt;
   Int_t           ph2_isPrompt;
   Float_t         ph1_SCEta;
   Float_t         ph2_SCEta;
   Float_t         ph1_SCPhi;
   Float_t         ph2_SCPhi;
   Float_t         ph1_hoe;
   Float_t         ph2_hoe;
   Float_t         ph1_sieie;
   Float_t         ph2_sieie;
   Float_t         ph1_pfchargedisogood03;
   Float_t         ph2_pfchargedisogood03;
   Float_t         ph1_pfchargedisobad04;
   Float_t         ph2_pfchargedisobad04;
   Float_t         ph1_etawidth;
   Float_t         ph2_etawidth;
   Float_t         ph1_phiwidth;
   Float_t         ph2_phiwidth;
   Float_t         ph1_eseffssqrt;
   Float_t         ph2_eseffssqrt;
   Float_t         ph1_pfchargedisobad03;
   Float_t         ph2_pfchargedisobad03;
   Float_t         ph1_sieip;
   Float_t         ph2_sieip;
   Float_t         ph1_sipip;
   Float_t         ph2_sipip;
   Float_t         ph1_ecaliso;
   Float_t         ph2_ecaliso;
   Float_t         ph1_ecalisobad;
   Float_t         ph2_ecalisobad;
   Float_t         ph1_badvtx_Et;
   Float_t         ph2_badvtx_Et;
   Float_t         ph1_isconv;
   Float_t         ph2_isconv;
   Int_t           ph1_ciclevel;
   Int_t           ph2_ciclevel;
   Float_t         ph1_sigmaEoE;
   Float_t         ph2_sigmaEoE;
   Float_t         ph1_ptoM;
   Float_t         ph2_ptoM;
   Int_t           ph1_isEB;
   Int_t           ph2_isEB;
   Float_t         ph1_s4ratio;
   Float_t         ph2_s4ratio;
   Float_t         ph1_e3x3;
   Float_t         ph2_e3x3;
   Float_t         ph1_e5x5;
   Float_t         ph2_e5x5;
   Float_t         PhotonsMass;
   Float_t         dipho_E;
   Float_t         dipho_pt;
   Float_t         dipho_eta;
   Float_t         dipho_phi;
   Float_t         dipho_cosThetaStar_CS;
   Float_t         dipho_tanhYStar;
   Float_t         dipho_Y;
   Int_t           vtx_ind;
   Float_t         vtx_x;
   Float_t         vtx_y;
   Float_t         vtx_z;
   Float_t         vtx_mva;
   Float_t         vtx_mva_2;
   Float_t         vtx_mva_3;
   Float_t         vtx_ptbal;
   Float_t         vtx_ptasym;
   Float_t         vtx_logsumpt2;
   Float_t         vtx_pulltoconv;
   Float_t         vtx_prob;
   Int_t           njets_passing_kLooseID;
   Float_t         j1_e;
   Float_t         j1_pt;
   Float_t         j1_phi;
   Float_t         j1_eta;
   Float_t         j1_beta;
   Float_t         j1_betaStar;
   Float_t         j1_betaStarClassic;
   Float_t         j1_dR2Mean;
   Float_t         j1_csvBtag;
   Float_t         j1_csvMvaBtag;
   Float_t         j1_jetProbBtag;
   Float_t         j1_tcheBtag;
   Float_t         j2_e;
   Float_t         j2_pt;
   Float_t         j2_phi;
   Float_t         j2_eta;
   Float_t         j2_beta;
   Float_t         j2_betaStar;
   Float_t         j2_betaStarClassic;
   Float_t         j2_dR2Mean;
   Float_t         j2_csvBtag;
   Float_t         j2_csvMvaBtag;
   Float_t         j2_jetProbBtag;
   Float_t         j2_tcheBtag;
   Float_t         j3_e;
   Float_t         j3_pt;
   Float_t         j3_phi;
   Float_t         j3_eta;
   Float_t         j3_beta;
   Float_t         j3_betaStar;
   Float_t         j3_betaStarClassic;
   Float_t         j3_dR2Mean;
   Float_t         j3_csvBtag;
   Float_t         j3_csvMvaBtag;
   Float_t         j3_jetProbBtag;
   Float_t         j3_tcheBtag;
   Float_t         j4_e;
   Float_t         j4_pt;
   Float_t         j4_phi;
   Float_t         j4_eta;
   Float_t         j4_beta;
   Float_t         j4_betaStar;
   Float_t         j4_betaStarClassic;
   Float_t         j4_dR2Mean;
   Float_t         j4_csvBtag;
   Float_t         j4_csvMvaBtag;
   Float_t         j4_jetProbBtag;
   Float_t         j4_tcheBtag;
   Float_t         JetsMass;
   Float_t         dijet_E;
   Float_t         dijet_Pt;
   Float_t         dijet_Eta;
   Float_t         dijet_Phi;
   Float_t         RadMass;
   Float_t         radion_E;
   Float_t         radion_Pt;
   Float_t         radion_Eta;
   Float_t         radion_Phi;

   // List of branches
   TBranch        *b_itype;   //!
   TBranch        *b_run;   //!
   TBranch        *b_lumis;   //!
   TBranch        *b_weight;   //!
   TBranch        *b_evweight;   //!
   TBranch        *b_pu_weight;   //!
   TBranch        *b_pu_n;   //!
   TBranch        *b_nvtx;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_category;   //!
   TBranch        *b_ph1_e;   //!
   TBranch        *b_ph2_e;   //!
   TBranch        *b_ph1_pt;   //!
   TBranch        *b_ph2_pt;   //!
   TBranch        *b_ph1_phi;   //!
   TBranch        *b_ph2_phi;   //!
   TBranch        *b_ph1_eta;   //!
   TBranch        *b_ph2_eta;   //!
   TBranch        *b_ph1_r9;   //!
   TBranch        *b_ph2_r9;   //!
   TBranch        *b_ph1_isPrompt;   //!
   TBranch        *b_ph2_isPrompt;   //!
   TBranch        *b_ph1_SCEta;   //!
   TBranch        *b_ph2_SCEta;   //!
   TBranch        *b_ph1_SCPhi;   //!
   TBranch        *b_ph2_SCPhi;   //!
   TBranch        *b_ph1_hoe;   //!
   TBranch        *b_ph2_hoe;   //!
   TBranch        *b_ph1_sieie;   //!
   TBranch        *b_ph2_sieie;   //!
   TBranch        *b_ph1_pfchargedisogood03;   //!
   TBranch        *b_ph2_pfchargedisogood03;   //!
   TBranch        *b_ph1_pfchargedisobad04;   //!
   TBranch        *b_ph2_pfchargedisobad04;   //!
   TBranch        *b_ph1_etawidth;   //!
   TBranch        *b_ph2_etawidth;   //!
   TBranch        *b_ph1_phiwidth;   //!
   TBranch        *b_ph2_phiwidth;   //!
   TBranch        *b_ph1_eseffssqrt;   //!
   TBranch        *b_ph2_eseffssqrt;   //!
   TBranch        *b_ph1_pfchargedisobad03;   //!
   TBranch        *b_ph2_pfchargedisobad03;   //!
   TBranch        *b_ph1_sieip;   //!
   TBranch        *b_ph2_sieip;   //!
   TBranch        *b_ph1_sipip;   //!
   TBranch        *b_ph2_sipip;   //!
   TBranch        *b_ph1_ecaliso;   //!
   TBranch        *b_ph2_ecaliso;   //!
   TBranch        *b_ph1_ecalisobad;   //!
   TBranch        *b_ph2_ecalisobad;   //!
   TBranch        *b_ph1_badvtx_Et;   //!
   TBranch        *b_ph2_badvtx_Et;   //!
   TBranch        *b_ph1_isconv;   //!
   TBranch        *b_ph2_isconv;   //!
   TBranch        *b_ph1_ciclevel;   //!
   TBranch        *b_ph2_ciclevel;   //!
   TBranch        *b_ph1_sigmaEoE;   //!
   TBranch        *b_ph2_sigmaEoE;   //!
   TBranch        *b_ph1_ptoM;   //!
   TBranch        *b_ph2_ptoM;   //!
   TBranch        *b_ph1_isEB;   //!
   TBranch        *b_ph2_isEB;   //!
   TBranch        *b_ph1_s4ratio;   //!
   TBranch        *b_ph2_s4ratio;   //!
   TBranch        *b_ph1_e3x3;   //!
   TBranch        *b_ph2_e3x3;   //!
   TBranch        *b_ph1_e5x5;   //!
   TBranch        *b_ph2_e5x5;   //!
   TBranch        *b_PhotonsMass;   //!
   TBranch        *b_dipho_E;   //!
   TBranch        *b_dipho_pt;   //!
   TBranch        *b_dipho_eta;   //!
   TBranch        *b_dipho_phi;   //!
   TBranch        *b_dipho_cosThetaStar_CS;   //!
   TBranch        *b_dipho_tanhYStar;   //!
   TBranch        *b_dipho_Y;   //!
   TBranch        *b_vtx_ind;   //!
   TBranch        *b_vtx_x;   //!
   TBranch        *b_vtx_y;   //!
   TBranch        *b_vtx_z;   //!
   TBranch        *b_vtx_mva;   //!
   TBranch        *b_vtx_mva_2;   //!
   TBranch        *b_vtx_mva_3;   //!
   TBranch        *b_vtx_ptbal;   //!
   TBranch        *b_vtx_ptasym;   //!
   TBranch        *b_vtx_logsumpt2;   //!
   TBranch        *b_vtx_pulltoconv;   //!
   TBranch        *b_vtx_prob;   //!
   TBranch        *b_njets_passing_kLooseID;   //!
   TBranch        *b_j1_e;   //!
   TBranch        *b_j1_pt;   //!
   TBranch        *b_j1_phi;   //!
   TBranch        *b_j1_eta;   //!
   TBranch        *b_j1_beta;   //!
   TBranch        *b_j1_betaStar;   //!
   TBranch        *b_j1_betaStarClassic;   //!
   TBranch        *b_j1_dR2Mean;   //!
   TBranch        *b_j1_csvBtag;   //!
   TBranch        *b_j1_csvMvaBtag;   //!
   TBranch        *b_j1_jetProbBtag;   //!
   TBranch        *b_j1_tcheBtag;   //!
   TBranch        *b_j2_e;   //!
   TBranch        *b_j2_pt;   //!
   TBranch        *b_j2_phi;   //!
   TBranch        *b_j2_eta;   //!
   TBranch        *b_j2_beta;   //!
   TBranch        *b_j2_betaStar;   //!
   TBranch        *b_j2_betaStarClassic;   //!
   TBranch        *b_j2_dR2Mean;   //!
   TBranch        *b_j2_csvBtag;   //!
   TBranch        *b_j2_csvMvaBtag;   //!
   TBranch        *b_j2_jetProbBtag;   //!
   TBranch        *b_j2_tcheBtag;   //!
   TBranch        *b_j3_e;   //!
   TBranch        *b_j3_pt;   //!
   TBranch        *b_j3_phi;   //!
   TBranch        *b_j3_eta;   //!
   TBranch        *b_j3_beta;   //!
   TBranch        *b_j3_betaStar;   //!
   TBranch        *b_j3_betaStarClassic;   //!
   TBranch        *b_j3_dR2Mean;   //!
   TBranch        *b_j3_csvBtag;   //!
   TBranch        *b_j3_csvMvaBtag;   //!
   TBranch        *b_j3_jetProbBtag;   //!
   TBranch        *b_j3_tcheBtag;   //!
   TBranch        *b_j4_e;   //!
   TBranch        *b_j4_pt;   //!
   TBranch        *b_j4_phi;   //!
   TBranch        *b_j4_eta;   //!
   TBranch        *b_j4_beta;   //!
   TBranch        *b_j4_betaStar;   //!
   TBranch        *b_j4_betaStarClassic;   //!
   TBranch        *b_j4_dR2Mean;   //!
   TBranch        *b_j4_csvBtag;   //!
   TBranch        *b_j4_csvMvaBtag;   //!
   TBranch        *b_j4_jetProbBtag;   //!
   TBranch        *b_j4_tcheBtag;   //!
   TBranch        *b_JetsMass;   //!
   TBranch        *b_dijet_E;   //!
   TBranch        *b_dijet_Pt;   //!
   TBranch        *b_dijet_Eta;   //!
   TBranch        *b_dijet_Phi;   //!
   TBranch        *b_RadMass;   //!
   TBranch        *b_radion_E;   //!
   TBranch        *b_radion_Pt;   //!
   TBranch        *b_radion_Eta;   //!
   TBranch        *b_radion_Phi;   //!

   fillPlot2012_radion_commonNtp( const std::string& dataset, const std::string& selectionType, const std::string& bTaggerType="JP" );
   virtual ~fillPlot2012_radion_commonNtp();
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init();
   void setSelectionType( const std::string& selectionType );
   virtual void     finalize();
   virtual int myHighestPtJet(vector<int> jets);
   virtual bool passCutBasedJetId(int jet);
   double delta_phi(double phi1, double phi2);

   std::string selectionType_;
   std::string bTaggerType_;

   bool dopureeventWeight_;
   int cicselection;

   float ptphot1cut;
   float ptphot2cut;

   float pthiggsmincut;
   float pthiggsmaxcut;

   float ptjetacccut;
   float etajetacccut;

   float ptdiphotcut_1bj;
   float ptdiphotcut_2bj;

   float etadiphotcut_1bj;
   float etadiphotcut_2bj;

   float ptphot1cut_1bj;
   float ptphot1cut_2bj;

   float ptjet1cut_1bj;
   float ptjet1cut_2bj;

   float costhetascut_1bj;
   float costhetascut_2bj;

   float mjjmincut;
   float mjjmaxcut;

   // vectors for jets
   float ptcorrjet[4], ecorrjet[4];
   float etajet[4], phijet[4];              
   float btagjprobjet[4], btagcsvjet[4];
};

#endif

