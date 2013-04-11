//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Feb 22 09:12:08 2013 by ROOT version 5.32/00
// from TTree AnaTree/Reduced tree for final analysis
// found on file: redntp_DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1_47.root
//////////////////////////////////////////////////////////

#ifndef fillPlot2012_radion_h
#define fillPlot2012_radion_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TRandom3.h>

#include "RedNtpFinalizer.h"

using namespace std;

class fillPlot2012_radion : public RedNtpFinalizer {
  
 public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain
  
  // Declaration of leaf types
  Int_t           run;
  Int_t           event;
  Int_t           lumi;
  Bool_t          H_event;
  Bool_t          V_event;
  Bool_t          WH_event;
  Bool_t          ZH_event;
  Bool_t          Zbb_event;
  Bool_t          Vqq_event;
  Float_t         rhoPF;
  Float_t         massgg;
  Float_t         ptgg;
  Float_t         ptggnewvtx;
  Float_t         phigg;
  Float_t         etagg;
  Float_t         massggnewvtx;
  Float_t         ptphot1;
  Float_t         ptphot2;
  Float_t         deltaRToTrackphot1;
  Float_t         deltaRToTrackphot2;
  Float_t         timephot1;
  Float_t         timephot2;
  Float_t         etaphot1;
  Float_t         etaphot2;
  Float_t         phiphot1;
  Float_t         phiphot2;
  Float_t         etascphot1;
  Float_t         etascphot2;
  Float_t         phiscphot1;
  Float_t         phiscphot2;
  Float_t         E1phot1;
  Float_t         E1phot2;
  Float_t         E9phot1;
  Float_t         E9phot2;
  Float_t         energyErrphot1;
  Float_t         energyErrphot2;
  Float_t         energySmearingphot1;
  Float_t         energySmearingphot2;
  Float_t         r9phot1;
  Float_t         r9phot2;
  Int_t           isemEGphot1;
  Int_t           isemEGphot2;
  Int_t           promptGamma;
  Int_t           LOGamma;
  Int_t           ISRGamma;
  Int_t           FSRGamma;
  Float_t         idmvaphot1;
  Float_t         idmvaphot2;
  Int_t           idcicphot1;
  Int_t           idcicphot2;
  Int_t           idcicnoelvetophot1;
  Int_t           idcicnoelvetophot2;
  Int_t           idcicpfphot1;
  Int_t           idcicpfphot2;
  Int_t           idcicpfnoelvetophot1;
  Int_t           idcicpfnoelvetophot2;
  Int_t           idelephot1;
  Int_t           idelephot2;
  Int_t           pid_isEMphot1;
  Int_t           pid_isEMphot2;
  Int_t           pid_haspixelseedphot1;
  Int_t           pid_haspixelseedphot2;
  Float_t         pid_jurECALphot1;
  Float_t         pid_jurECALphot2;
  Float_t         pid_twrHCALphot1;
  Float_t         pid_twrHCALphot2;
  Float_t         pid_HoverEphot1;
  Float_t         pid_HoverEphot2;
  Float_t         pid_hlwTrackphot1;
  Float_t         pid_hlwTrackphot2;
  Float_t         pid_etawidphot1;
  Float_t         pid_etawidphot2;
  Float_t         pid_hlwTrackNoDzphot1;
  Float_t         pid_hlwTrackNoDzphot2;
  Int_t           pid_hasMatchedConvphot1;
  Int_t           pid_hasMatchedConvphot2;
  Int_t           pid_hasMatchedPromptElephot1;
  Int_t           pid_hasMatchedPromptElephot2;
  Float_t         pid_sminphot1;
  Float_t         pid_sminphot2;
  Float_t         pid_smajphot1;
  Float_t         pid_smajphot2;
  Int_t           pid_ntrkphot1;
  Int_t           pid_ntrkphot2;
  Float_t         pid_ptisophot1;
  Float_t         pid_ptisophot2;
  Int_t           pid_ntrkcsphot1;
  Int_t           pid_ntrkcsphot2;
  Float_t         pid_ptisocsphot1;
  Float_t         pid_ptisocsphot2;
  Float_t         pid_ecalisophot1;
  Float_t         pid_ecalisophot2;
  Float_t         pid_hcalisophot1;
  Float_t         pid_hcalisophot2;
  Int_t           njets;
  Float_t         ecorrjet[10];   //[njets]
  Float_t         ptjet[10];   //[njets]
  Float_t         ptcorrjet[10];   //[njets]
  Float_t         etajet[10];   //[njets]
  Float_t         phijet[10];   //[njets]
  Float_t         betajet[10];   //[njets]
  Float_t         betastarjet[10];   //[njets]
  Float_t         btagvtxjet[10];   //[njets]
  Float_t         btagcsvjet[10];   //[njets]
  Float_t         btagjprobjet[10];   //[njets]
  Float_t         btagcsvmvajet[10];   //[njets]
  Float_t         btagbjprobjet[10];   //[njets]
  Float_t         btagvtxhighpurjet[10];   //[njets]
  Float_t         btagtchejet[10];   //[njets]
  Float_t         btagtchpjet[10];   //[njets]
  Float_t         ptDjet[10];   //[njets]
  Float_t         ptD_QCjet[10];   //[njets]
  Float_t         axis2_QCjet[10];   //[njets]
  Float_t         rmsjet[10];   //[njets]
  Int_t           ntrkjet[10];   //[njets]
  Int_t           nneutjet[10];   //[njets]
  Int_t           nChg_QCjet[10];   //[njets]
  Int_t           nNeutral_ptCutjet[10];   //[njets]
  Float_t         jetIdSimple_mvajet[10];   //[njets]
  Float_t         jetIdFull_mvajet[10];   //[njets]
  Float_t         jetId_dR2Meanjet[10];   //[njets]
  Float_t         jetId_betaStarClassicjet[10];   //[njets]
  Float_t         jetId_frac01jet[10];   //[njets]
  Float_t         jetId_frac02jet[10];   //[njets]
  Float_t         jetId_frac03jet[10];   //[njets]
  Float_t         jetId_frac04jet[10];   //[njets]
  Float_t         jetId_frac05jet[10];   //[njets]
  Float_t         jetId_betajet[10];   //[njets]
  Float_t         jetId_betaStarjet[10];   //[njets]
  Float_t         jetIdCutBased_wpjet[10];   //[njets]
  Float_t         jetIdSimple_wpjet[10];   //[njets]
  Float_t         jetIdFull_wpjet[10];   //[njets]
  Int_t           assjet[10];   //[njets]
  Int_t           partPdgIDjet[10];   //[njets]
  Int_t           partMomPdgIDjet[10];   //[njets]
  Float_t         nvtx;
  Int_t           vtxId;
  Float_t         vtxPos_x;
  Float_t         vtxPos_y;
  Float_t         vtxPos_z;
  Float_t         vtxIdMVA;
  Float_t         vtxIdEvtProb;
  Float_t         diPhotMVA;
  Float_t         diPhotMVA_vtx0;
  Float_t         diPhotMVA_vtxPair;
  Int_t           preselPairId;
  Float_t         tmva_dipho_MIT_dmom;
  Float_t         tmva_dipho_MIT_dmom_wrong_vtx;
  Float_t         tmva_dipho_MIT_vtxprob;
  Float_t         tmva_dipho_MIT_ptom1;
  Float_t         tmva_dipho_MIT_ptom2;
  Float_t         tmva_dipho_MIT_eta1;
  Float_t         tmva_dipho_MIT_eta2;
  Float_t         tmva_dipho_MIT_dphi;
  Float_t         tmva_dipho_MIT_ph1mva;
  Float_t         tmva_dipho_MIT_ph2mva;
  Float_t         ePUMet;
  Float_t         ePUMet2;
  Float_t         ePUMet3;
  Float_t         ePUMet4;
  Float_t         ePUMet5;
  Float_t         ecorrPUMet5;
  Float_t         phiPUMet;
  Float_t         phiPUMet2;
  Float_t         phiPUMet3;
  Float_t         phiPUMet4;
  Float_t         phiPUMet5;
  Float_t         phiCorrPUMet5;
  Float_t         phot1Metx;
  Float_t         phot2Metx;
  Float_t         leptonsMetx;
  Float_t         part_in_jetsMetx;
  Float_t         chg_vtx_unclMetx;
  Float_t         chg_novtx_unclMetx;
  Float_t         neutrals_unclMetx;
  Float_t         part_fail_puidMetx;
  Float_t         phot1Mety;
  Float_t         phot2Mety;
  Float_t         leptonsMety;
  Float_t         part_in_jetsMety;
  Float_t         chg_vtx_unclMety;
  Float_t         chg_novtx_unclMety;
  Float_t         neutrals_unclMety;
  Float_t         part_fail_puidMety;
  Float_t         scaling;
  Float_t         sMet;
  Float_t         eMet;
  Float_t         phiMet;
  Float_t         signifMet;
  Float_t         eSmearedMet;
  Float_t         phiSmearedMet;
  Float_t         eShiftedMet;
  Float_t         phiShiftedMet;
  Float_t         eShiftedScaledMet;
  Float_t         phiShiftedScaledMet;
  Float_t         eSmearedShiftedMet;
  Float_t         phiSmearedShiftedMet;
  Float_t         eShiftedScaledMetPUcorr;
  Float_t         phiShiftedScaledMetPUcorr;
  Float_t         eSmearedShiftedMePUcorrt;
  Float_t         phiSmearedShiftedMetPUcorr;
  Float_t         sCorrMet;
  Float_t         eCorrMet;
  Float_t         phiCorrMet;
  Float_t         signifCorrMet;
  Float_t         smuCorrMet;
  Float_t         emuCorrMet;
  Float_t         phimuCorrMet;
  Float_t         signifmuCorrMet;
  Float_t         sNoHFMet;
  Float_t         eNoHFMet;
  Float_t         phiNoHFMet;
  Float_t         signifNoHFMet;
  Float_t         stcMet;
  Float_t         etcMet;
  Float_t         phitcMet;
  Float_t         signiftcMet;
  Float_t         sglobalPfMet;
  Float_t         eglobalPfMet;
  Float_t         phiglobalPfMet;
  Float_t         signifglobalPfMet;
  Float_t         scentralPfMet;
  Float_t         ecentralPfMet;
  Float_t         phicentralPfMet;
  Float_t         signifcentralPfMet;
  Float_t         eassocPfMet;
  Float_t         phiassocPfMet;
  Float_t         signifassocPfMet;
  Float_t         eassocOtherVtxPfMet;
  Float_t         phiassocOtherVtxPfMet;
  Float_t         signifassocOtherVtxPfMet;
  Float_t         etrkPfMet;
  Float_t         phitrkPfMet;
  Float_t         signiftrkPfMet;
  Float_t         ecleanPfMet;
  Float_t         phicleanPfMet;
  Float_t         signifcleanPfMet;
  Float_t         ecleanedSaclayPfMet;
  Float_t         phicleanedSaclayPfMet;
  Float_t         signifcleanedSaclayPfMet;
  Float_t         eminTypeICleanSaclayPfMet;
  Float_t         phiminTypeICleanSaclayPfMet;
  Float_t         signifminTypeICleanSaclayPfMet;
  Float_t         globalPfSums;
  Float_t         spfMet;
  Float_t         epfMet;
  Float_t         phipfMet;
  Float_t         signifpfMet;
  Float_t         spfMetType1;
  Float_t         epfMetType1;
  Float_t         phipfMetType1;
  Float_t         signifpfMetType1;
  Float_t         sMetGen;
  Float_t         eMetGen;
  Float_t         phiMetGen;
  Float_t         signifMetGen;
  Float_t         sMetGen2;
  Float_t         eMetGen2;
  Float_t         phiMetGen2;
  Int_t           npu;
  Int_t           NtotEvents;
  Float_t         xsection;
  Float_t         EquivLumi;
  Int_t           SampleID;
  Float_t         pu_weight;
  Float_t         pt_weight;
  Int_t           gen_custom_processId;
  Float_t         gen_pt_gamma1;
  Float_t         gen_pt_gamma2;
  Float_t         gen_eta_gamma1;
  Float_t         gen_eta_gamma2;
  Float_t         gen_phi_gamma1;
  Float_t         gen_phi_gamma2;
  Float_t         gen_mass_diphoton;
  Float_t         gen_pt_diphoton;
  Float_t         gen_eta_diphoton;
  Float_t         gen_phi_diphoton;
  Float_t         gen_mass_dijet;
  Float_t         gen_pt_dijet;
  Float_t         gen_eta_dijet;
  Float_t         gen_phi_dijet;
  Float_t         gen_zeppenfeld;
  Float_t         gen_pt_lep1;
  Float_t         gen_pt_lep2;
  Float_t         gen_eta_lep1;
  Float_t         gen_eta_lep2;
  Float_t         gen_phi_lep1;
  Float_t         gen_phi_lep2;
  Int_t           gen_pid_lep1;
  Int_t           gen_pid_lep2;
  Float_t         ptele1;
  Float_t         ptele2;
  Float_t         etaele1;
  Float_t         etaele2;
  Float_t         phiele1;
  Float_t         phiele2;
  Float_t         eneele1;
  Float_t         eneele2;
  Float_t         sIeIeele1;
  Float_t         sIeIeele2;
  Float_t         dphiele1;
  Float_t         dphiele2;
  Float_t         detaele1;
  Float_t         detaele2;
  Float_t         hoeele1;
  Float_t         hoeele2;
  Int_t           mhitsele1;
  Int_t           mhitsele2;
  Float_t         d0ele1;
  Float_t         d0ele2;
  Float_t         dzele1;
  Float_t         dzele2;
  Float_t         invMassele1g1;
  Float_t         invMassele1g2;
  Float_t         invMassele2g1;
  Float_t         invMassele2g2;
  Float_t         oEmoPele1;
  Float_t         oEmoPele2;
  Float_t         mvanotrigele1;
  Float_t         mvanotrigele2;
  Float_t         mvatrigele1;
  Float_t         mvatrigele2;
  Int_t           matchconvele1;
  Int_t           matchconvele2;
  Float_t         chHadIso03ele1;
  Float_t         chHadIso03ele2;
  Float_t         nHadIso03ele1;
  Float_t         nHadIso03ele2;
  Float_t         photIso03ele1;
  Float_t         photIso03ele2;
  Float_t         pteleloose1;
  Float_t         pteleloose2;
  Float_t         etaeleloose1;
  Float_t         etaeleloose2;
  Float_t         phieleloose1;
  Float_t         phieleloose2;
  Float_t         eneeleloose1;
  Float_t         eneeleloose2;
  Float_t         sIeIeeleloose1;
  Float_t         sIeIeeleloose2;
  Float_t         dphieleloose1;
  Float_t         dphieleloose2;
  Float_t         detaeleloose1;
  Float_t         detaeleloose2;
  Float_t         hoeeleloose1;
  Float_t         hoeeleloose2;
  Int_t           mhitseleloose1;
  Int_t           mhitseleloose2;
  Float_t         d0eleloose1;
  Float_t         d0eleloose2;
  Float_t         dzeleloose1;
  Float_t         dzeleloose2;
  Float_t         invMasseleloose1g1;
  Float_t         invMasseleloose1g2;
  Float_t         invMasseleloose2g1;
  Float_t         invMasseleloose2g2;
  Float_t         oEmoPeleloose1;
  Float_t         oEmoPeleloose2;
  Float_t         mvanotrigeleloose1;
  Float_t         mvanotrigeleloose2;
  Float_t         mvatrigeleloose1;
  Float_t         mvatrigeleloose2;
  Int_t           matchconveleloose1;
  Int_t           matchconveleloose2;
  Float_t         chHadIso03eleloose1;
  Float_t         chHadIso03eleloose2;
  Float_t         nHadIso03eleloose1;
  Float_t         nHadIso03eleloose2;
  Float_t         photIso03eleloose1;
  Float_t         photIso03eleloose2;
  Float_t         ptelenontr801;
  Float_t         ptelenontr802;
  Float_t         etaelenontr801;
  Float_t         etaelenontr802;
  Float_t         phielenontr801;
  Float_t         phielenontr802;
  Float_t         eneelenontr801;
  Float_t         eneelenontr802;
  Float_t         ptelenontr901;
  Float_t         ptelenontr902;
  Float_t         etaelenontr901;
  Float_t         etaelenontr902;
  Float_t         phielenontr901;
  Float_t         phielenontr902;
  Float_t         eneelenontr901;
  Float_t         eneelenontr902;
  Float_t         ptmu1;
  Float_t         ptmu2;
  Float_t         etamu1;
  Float_t         etamu2;
  Float_t         phimu1;
  Float_t         phimu2;
  Float_t         enemu1;
  Float_t         enemu2;
  Int_t           pixhitsmu1;
  Int_t           pixhitsmu2;
  Int_t           trkhitsmu1;
  Int_t           trkhitsmu2;
  Int_t           hitsmu1;
  Int_t           hitsmu2;
  Float_t         chi2mu1;
  Float_t         chi2mu2;
  Int_t           matchmu1;
  Int_t           matchmu2;
  Float_t         d0mu1;
  Float_t         d0mu2;
  Float_t         dzmu1;
  Float_t         dzmu2;
  Float_t         chHadmu1;
  Float_t         chHadmu2;
  Float_t         nHadmu1;
  Float_t         nHadmu2;
  Float_t         photmu1;
  Float_t         photmu2;
  Float_t         puptmu1;
  Float_t         puptmu2;
  Float_t         ptmuloose1;
  Float_t         ptmuloose2;
  Float_t         etamuloose1;
  Float_t         etamuloose2;
  Float_t         phimuloose1;
  Float_t         phimuloose2;
  Float_t         enemuloose1;
  Float_t         enemuloose2;
  Int_t           pixhitsmuloose1;
  Int_t           pixhitsmuloose2;
  Int_t           trkhitsmuloose1;
  Int_t           trkhitsmuloose2;
  Int_t           hitsmuloose1;
  Int_t           hitsmuloose2;
  Float_t         chi2muloose1;
  Float_t         chi2muloose2;
  Int_t           matchmuloose1;
  Int_t           matchmuloose2;
  Float_t         d0muloose1;
  Float_t         d0muloose2;
  Float_t         dzmuloose1;
  Float_t         dzmuloose2;
  Float_t         ptmuvloose1;
  Float_t         ptmuvloose2;
  Float_t         etamuvloose1;
  Float_t         etamuvloose2;
  Float_t         phimuvloose1;
  Float_t         phimuvloose2;
  Float_t         enemuvloose1;
  Float_t         enemuvloose2;
  Int_t           pixhitsmuvloose1;
  Int_t           pixhitsmuvloose2;
  Int_t           trkhitsmuvloose1;
  Int_t           trkhitsmuvloose2;
  Int_t           hitsmuvloose1;
  Int_t           hitsmuvloose2;
  Float_t         chi2muvloose1;
  Float_t         chi2muvloose2;
  Int_t           matchmuvloose1;
  Int_t           matchmuvloose2;
  Float_t         d0muvloose1;
  Float_t         d0muvloose2;
  Float_t         dzmuvloose1;
  Float_t         dzmuvloose2;
  Int_t           hasPassedSinglePhot;
  Int_t           hasPassedDoublePhot;
  Float_t         chHadmuloose1;
  Float_t         chHadmuloose2;
  Float_t         nHadmuloose1;
  Float_t         nHadmuloose2;
  Float_t         photmuloose1;
  Float_t         photmuloose2;
  Float_t         puptmuloose1;
  Float_t         puptmuloose2;
  Float_t         chHadmuvloose1;
  Float_t         chHadmuvloose2;
  Float_t         nHadmuvloose1;
  Float_t         nHadmuvloose2;
  Float_t         photmuvloose1;
  Float_t         photmuvloose2;
  Float_t         puptmuvloose1;
  Float_t         puptmuvloose2;

  // List of branches
  TBranch        *b_run;   //!
  TBranch        *b_event;   //!
  TBranch        *b_lumi;   //!
  TBranch        *b_H_event;   //!
  TBranch        *b_V_event;   //!
  TBranch        *b_WH_event;   //!
  TBranch        *b_ZH_event;   //!
  TBranch        *b_Zbb_event;   //!
  TBranch        *b_Vqq_event;   //!
  TBranch        *b_rhoPF;   //!
  TBranch        *b_massgg;   //!
  TBranch        *b_ptgg;   //!
  TBranch        *b_ptggnewvtx;   //!
  TBranch        *b_phigg;   //!
  TBranch        *b_etagg;   //!
  TBranch        *b_massggnewvtx;   //!
  TBranch        *b_ptphot1;   //!
  TBranch        *b_ptphot2;   //!
  TBranch        *b_deltaRToTrackphot1;   //!
  TBranch        *b_deltaRToTrackphot2;   //!
  TBranch        *b_timephot1;   //!
  TBranch        *b_timephot2;   //!
  TBranch        *b_etaphot1;   //!
  TBranch        *b_etaphot2;   //!
  TBranch        *b_phiphot1;   //!
  TBranch        *b_phiphot2;   //!
  TBranch        *b_etascphot1;   //!
  TBranch        *b_etascphot2;   //!
  TBranch        *b_phiscphot1;   //!
  TBranch        *b_phiscphot2;   //!
  TBranch        *b_E1phot1;   //!
  TBranch        *b_E1phot2;   //!
  TBranch        *b_E9phot1;   //!
  TBranch        *b_E9phot2;   //!
  TBranch        *b_energyErrphot1;   //!
  TBranch        *b_energyErrphot2;   //!
  TBranch        *b_energySmearingphot1;   //!
  TBranch        *b_energySmearingphot2;   //!
  TBranch        *b_r9phot1;   //!
  TBranch        *b_r9phot2;   //!
  TBranch        *b_isemEGphot1;   //!
  TBranch        *b_isemEGphot2;   //!
  TBranch        *b_promptGamma;   //!
  TBranch        *b_LOGamma;   //!
  TBranch        *b_ISRGamma;   //!
  TBranch        *b_FSRGamma;   //!
  TBranch        *b_idmvaphot1;   //!
  TBranch        *b_idmvaphot2;   //!
  TBranch        *b_idcicphot1;   //!
  TBranch        *b_idcicphot2;   //!
  TBranch        *b_idcicnoelvetophot1;   //!
  TBranch        *b_idcicnoelvetophot2;   //!
  TBranch        *b_idcicpfphot1;   //!
  TBranch        *b_idcicpfphot2;   //!
  TBranch        *b_idcicpfnoelvetophot1;   //!
  TBranch        *b_idcicpfnoelvetophot2;   //!
  TBranch        *b_idelephot1;   //!
  TBranch        *b_idelephot2;   //!
  TBranch        *b_pid_isEMphot1;   //!
  TBranch        *b_pid_isEMphot2;   //!
  TBranch        *b_pid_haspixelseedphot1;   //!
  TBranch        *b_pid_haspixelseedphot2;   //!
  TBranch        *b_pid_jurECALphot1;   //!
  TBranch        *b_pid_jurECALphot2;   //!
  TBranch        *b_pid_twrHCALphot1;   //!
  TBranch        *b_pid_twrHCALphot2;   //!
  TBranch        *b_pid_HoverEphot1;   //!
  TBranch        *b_pid_HoverEphot2;   //!
  TBranch        *b_pid_hlwTrackphot1;   //!
  TBranch        *b_pid_hlwTrackphot2;   //!
  TBranch        *b_pid_etawidphot1;   //!
  TBranch        *b_pid_etawidphot2;   //!
  TBranch        *b_pid_hlwTrackNoDzphot1;   //!
  TBranch        *b_pid_hlwTrackNoDzphot2;   //!
  TBranch        *b_pid_hasMatchedConvphot1;   //!
  TBranch        *b_pid_hasMatchedConvphot2;   //!
  TBranch        *b_pid_hasMatchedPromptElephot1;   //!
  TBranch        *b_pid_hasMatchedPromptElephot2;   //!
  TBranch        *b_pid_sminphot1;   //!
  TBranch        *b_pid_sminphot2;   //!
  TBranch        *b_pid_smajphot1;   //!
  TBranch        *b_pid_smajphot2;   //!
  TBranch        *b_pid_ntrkphot1;   //!
  TBranch        *b_pid_ntrkphot2;   //!
  TBranch        *b_pid_ptisophot1;   //!
  TBranch        *b_pid_ptisophot2;   //!
  TBranch        *b_pid_ntrkcsphot1;   //!
  TBranch        *b_pid_ntrkcsphot2;   //!
  TBranch        *b_pid_ptisocsphot1;   //!
  TBranch        *b_pid_ptisocsphot2;   //!
  TBranch        *b_pid_ecalisophot1;   //!
  TBranch        *b_pid_ecalisophot2;   //!
  TBranch        *b_pid_hcalisophot1;   //!
  TBranch        *b_pid_hcalisophot2;   //!
  TBranch        *b_njets;   //!
  TBranch        *b_ecorrjet;   //!
  TBranch        *b_ptjet;   //!
  TBranch        *b_ptcorrjet;   //!
  TBranch        *b_etajet;   //!
  TBranch        *b_phijet;   //!
  TBranch        *b_betajet;   //!
  TBranch        *b_betastarjet;   //!
  TBranch        *b_btagvtxjet;   //!
  TBranch        *b_btagcsvjet;   //!
  TBranch        *b_btagjprobjet;   //!
  TBranch        *b_btagcsvmvajet;   //!
  TBranch        *b_btagbjprobjet;   //!
  TBranch        *b_btagvtxhighpurjet;   //!
  TBranch        *b_btagtchejet;   //!
  TBranch        *b_btagtchpjet;   //!
  TBranch        *b_ptDjet;   //!
  TBranch        *b_ptD_QCjet;   //!
  TBranch        *b_axis2_QCjet;   //!
  TBranch        *b_rmsjet;   //!
  TBranch        *b_ntrkjet;   //!
  TBranch        *b_nneutjet;   //!
  TBranch        *b_nChg_QCjet;   //!
  TBranch        *b_nNeutral_ptCutjet;   //!
  TBranch        *b_jetIdSimple_mvajet;   //!
  TBranch        *b_jetIdFull_mvajet;   //!
  TBranch        *b_jetId_dR2Meanjet;   //!
  TBranch        *b_jetId_betaStarClassicjet;   //!
  TBranch        *b_jetId_frac01jet;   //!
  TBranch        *b_jetId_frac02jet;   //!
  TBranch        *b_jetId_frac03jet;   //!
  TBranch        *b_jetId_frac04jet;   //!
  TBranch        *b_jetId_frac05jet;   //!
  TBranch        *b_jetId_betajet;   //!
  TBranch        *b_jetId_betaStarjet;   //!
  TBranch        *b_jetIdCutBased_wpjet;   //!
  TBranch        *b_jetIdSimple_wpjet;   //!
  TBranch        *b_jetIdFull_wpjet;   //!
  TBranch        *b_assjet;   //!
  TBranch        *b_partPdgIDjet;   //!
  TBranch        *b_partMomPdgIDjet;   //!
  TBranch        *b_nvtx;   //!
  TBranch        *b_vtxId;   //!
  TBranch        *b_vtxPos_x;   //!
  TBranch        *b_vtxPos_y;   //!
  TBranch        *b_vtxPos_z;   //!
  TBranch        *b_vtxIdMVA;   //!
  TBranch        *b_vtxIdEvtProb;   //!
  TBranch        *b_diPhotMVA;   //!
  TBranch        *b_diPhotMVA_vtx0;   //!
  TBranch        *b_diPhotMVA_vtxPair;   //!
  TBranch        *b_preselPairId;   //!
  TBranch        *b_tmva_dipho_MIT_dmom;   //!
  TBranch        *b_tmva_dipho_MIT_dmom_wrong_vtx;   //!
  TBranch        *b_tmva_dipho_MIT_vtxprob;   //!
  TBranch        *b_tmva_dipho_MIT_ptom1;   //!
  TBranch        *b_tmva_dipho_MIT_ptom2;   //!
  TBranch        *b_tmva_dipho_MIT_eta1;   //!
  TBranch        *b_tmva_dipho_MIT_eta2;   //!
  TBranch        *b_tmva_dipho_MIT_dphi;   //!
  TBranch        *b_tmva_dipho_MIT_ph1mva;   //!
  TBranch        *b_tmva_dipho_MIT_ph2mva;   //!
  TBranch        *b_ePUMet;   //!
  TBranch        *b_ePUMet2;   //!
  TBranch        *b_ePUMet3;   //!
  TBranch        *b_ePUMet4;   //!
  TBranch        *b_ePUMet5;   //!
  TBranch        *b_ecorrPUMet5;   //!
  TBranch        *b_phiPUMet;   //!
  TBranch        *b_phiPUMet2;   //!
  TBranch        *b_phiPUMet3;   //!
  TBranch        *b_phiPUMet4;   //!
  TBranch        *b_phiPUMet5;   //!
  TBranch        *b_phiCorrPUMet5;   //!
  TBranch        *b_phot1Metx;   //!
  TBranch        *b_phot2Metx;   //!
  TBranch        *b_leptonsMetx;   //!
  TBranch        *b_part_in_jetsMetx;   //!
  TBranch        *b_chg_vtx_unclMetx;   //!
  TBranch        *b_chg_novtx_unclMetx;   //!
  TBranch        *b_neutrals_unclMetx;   //!
  TBranch        *b_part_fail_puidMetx;   //!
  TBranch        *b_phot1Mety;   //!
  TBranch        *b_phot2Mety;   //!
  TBranch        *b_leptonsMety;   //!
  TBranch        *b_part_in_jetsMety;   //!
  TBranch        *b_chg_vtx_unclMety;   //!
  TBranch        *b_chg_novtx_unclMety;   //!
  TBranch        *b_neutrals_unclMety;   //!
  TBranch        *b_part_fail_puidMety;   //!
  TBranch        *b_scaling;   //!
  TBranch        *b_sMet;   //!
  TBranch        *b_eMet;   //!
  TBranch        *b_phiMet;   //!
  TBranch        *b_signifMet;   //!
  TBranch        *b_eSmearedMet;   //!
  TBranch        *b_phiSmearedMet;   //!
  TBranch        *b_eShiftedMet;   //!
  TBranch        *b_phiShiftedMet;   //!
  TBranch        *b_eShiftedScaledMet;   //!
  TBranch        *b_phiShiftedScaledMet;   //!
  TBranch        *b_eSmearedShiftedMet;   //!
  TBranch        *b_phiSmearedShiftedMet;   //!
  TBranch        *b_eShiftedScaledMetPUcorr;   //!
  TBranch        *b_phiShiftedScaledMetPUcorr;   //!
  TBranch        *b_eSmearedShiftedMetPUcorr;   //!
  TBranch        *b_phiSmearedShiftedMetPUcorr;   //!
  TBranch        *b_sCorrMet;   //!
  TBranch        *b_eCorrMet;   //!
  TBranch        *b_phiCorrMet;   //!
  TBranch        *b_signifCorrMet;   //!
  TBranch        *b_smuCorrMet;   //!
  TBranch        *b_emuCorrMet;   //!
  TBranch        *b_phimuCorrMet;   //!
  TBranch        *b_signifmuCorrMet;   //!
  TBranch        *b_sNoHFMet;   //!
  TBranch        *b_eNoHFMet;   //!
  TBranch        *b_phiNoHFMet;   //!
  TBranch        *b_signifNoHFMet;   //!
  TBranch        *b_stcMet;   //!
  TBranch        *b_etcMet;   //!
  TBranch        *b_phitcMet;   //!
  TBranch        *b_signiftcMet;   //!
  TBranch        *b_sglobalPfMet;   //!
  TBranch        *b_eglobalPfMet;   //!
  TBranch        *b_phiglobalPfMet;   //!
  TBranch        *b_signifglobalPfMet;   //!
  TBranch        *b_scentralPfMet;   //!
  TBranch        *b_ecentralPfMet;   //!
  TBranch        *b_phicentralPfMet;   //!
  TBranch        *b_signifcentralPfMet;   //!
  TBranch        *b_eassocPfMet;   //!
  TBranch        *b_phiassocPfMet;   //!
  TBranch        *b_signifassocPfMet;   //!
  TBranch        *b_eassocOtherVtxPfMet;   //!
  TBranch        *b_phiassocOtherVtxPfMet;   //!
  TBranch        *b_signifassocOtherVtxPfMet;   //!
  TBranch        *b_etrkPfMet;   //!
  TBranch        *b_phitrkPfMet;   //!
  TBranch        *b_signiftrkPfMet;   //!
  TBranch        *b_ecleanPfMet;   //!
  TBranch        *b_phicleanPfMet;   //!
  TBranch        *b_signifcleanPfMet;   //!
  TBranch        *b_ecleanedSaclayPfMet;   //!
  TBranch        *b_phicleanedSaclayPfMet;   //!
  TBranch        *b_signifcleanedSaclayPfMet;   //!
  TBranch        *b_eminTypeICleanSaclayPfMet;   //!
  TBranch        *b_phiminTypeICleanSaclayPfMet;   //!
  TBranch        *b_signifminTypeICleanSaclayPfMet;   //!
  TBranch        *b_globalPfSums;   //!
  TBranch        *b_spfMet;   //!
  TBranch        *b_epfMet;   //!
  TBranch        *b_phipfMet;   //!
  TBranch        *b_signifpfMet;   //!
  TBranch        *b_spfMetType1;   //!
  TBranch        *b_epfMetType1;   //!
  TBranch        *b_phipfMetType1;   //!
  TBranch        *b_signifpfMetType1;   //!
  TBranch        *b_sMetGen;   //!
  TBranch        *b_eMetGen;   //!
  TBranch        *b_phiMetGen;   //!
  TBranch        *b_signifMetGen;   //!
  TBranch        *b_sMetGen2;   //!
  TBranch        *b_eMetGen2;   //!
  TBranch        *b_phiMetGen2;   //!
  TBranch        *b_npu;   //!
  TBranch        *b_NtotEvents;   //!
  TBranch        *b_xsection;   //!
  TBranch        *b_EquivLumi;   //!
  TBranch        *b_SampleID;   //!
  TBranch        *b_pu_weight;   //!
  TBranch        *b_pt_weight;   //!
  TBranch        *b_gen_custom_processId;   //!
  TBranch        *b_gen_pt_gamma1;   //!
  TBranch        *b_gen_pt_gamma2;   //!
  TBranch        *b_gen_eta_gamma1;   //!
  TBranch        *b_gen_eta_gamma2;   //!
  TBranch        *b_gen_phi_gamma1;   //!
  TBranch        *b_gen_phi_gamma2;   //!
  TBranch        *b_gen_mass_diphoton;   //!
  TBranch        *b_gen_pt_diphoton;   //!
  TBranch        *b_gen_eta_diphoton;   //!
  TBranch        *b_gen_phi_diphoton;   //!
  TBranch        *b_gen_mass_dijet;   //!
  TBranch        *b_gen_pt_dijet;   //!
  TBranch        *b_gen_eta_dijet;   //!
  TBranch        *b_gen_phi_dijet;   //!
  TBranch        *b_gen_zeppenfeld;   //!
  TBranch        *b_gen_pt_lep1;   //!
  TBranch        *b_gen_pt_lep2;   //!
  TBranch        *b_gen_eta_lep1;   //!
  TBranch        *b_gen_eta_lep2;   //!
  TBranch        *b_gen_phi_lep1;   //!
  TBranch        *b_gen_phi_lep2;   //!
  TBranch        *b_gen_pid_lep1;   //!
  TBranch        *b_gen_pid_lep2;   //!
  TBranch        *b_ptele1;   //!
  TBranch        *b_ptele2;   //!
  TBranch        *b_etaele1;   //!
  TBranch        *b_etaele2;   //!
  TBranch        *b_phiele1;   //!
  TBranch        *b_phiele2;   //!
  TBranch        *b_eneele1;   //!
  TBranch        *b_eneele2;   //!
  TBranch        *b_sIeIeele1;   //!
  TBranch        *b_sIeIeele2;   //!
  TBranch        *b_dphiele1;   //!
  TBranch        *b_dphiele2;   //!
  TBranch        *b_detaele1;   //!
  TBranch        *b_detaele2;   //!
  TBranch        *b_hoeele1;   //!
  TBranch        *b_hoeele2;   //!
  TBranch        *b_mhitsele1;   //!
  TBranch        *b_mhitsele2;   //!
  TBranch        *b_d0ele1;   //!
  TBranch        *b_d0ele2;   //!
  TBranch        *b_dzele1;   //!
  TBranch        *b_dzele2;   //!
  TBranch        *b_invMassele1g1;   //!
  TBranch        *b_invMassele1g2;   //!
  TBranch        *b_invMassele2g1;   //!
  TBranch        *b_invMassele2g2;   //!
  TBranch        *b_oEmoPele1;   //!
  TBranch        *b_oEmoPele2;   //!
  TBranch        *b_mvanotrigele1;   //!
  TBranch        *b_mvanotrigele2;   //!
  TBranch        *b_mvatrigele1;   //!
  TBranch        *b_mvatrigele2;   //!
  TBranch        *b_matchconvele1;   //!
  TBranch        *b_matchconvele2;   //!
  TBranch        *b_chHadIso03ele1;   //!
  TBranch        *b_chHadIso03ele2;   //!
  TBranch        *b_nHadIso03ele1;   //!
  TBranch        *b_nHadIso03ele2;   //!
  TBranch        *b_photIso03ele1;   //!
  TBranch        *b_photIso03ele2;   //!
  TBranch        *b_pteleloose1;   //!
  TBranch        *b_pteleloose2;   //!
  TBranch        *b_etaeleloose1;   //!
  TBranch        *b_etaeleloose2;   //!
  TBranch        *b_phieleloose1;   //!
  TBranch        *b_phieleloose2;   //!
  TBranch        *b_eneeleloose1;   //!
  TBranch        *b_eneeleloose2;   //!
  TBranch        *b_sIeIeeleloose1;   //!
  TBranch        *b_sIeIeeleloose2;   //!
  TBranch        *b_dphieleloose1;   //!
  TBranch        *b_dphieleloose2;   //!
  TBranch        *b_detaeleloose1;   //!
  TBranch        *b_detaeleloose2;   //!
  TBranch        *b_hoeeleloose1;   //!
  TBranch        *b_hoeeleloose2;   //!
  TBranch        *b_mhitseleloose1;   //!
  TBranch        *b_mhitseleloose2;   //!
  TBranch        *b_d0eleloose1;   //!
  TBranch        *b_d0eleloose2;   //!
  TBranch        *b_dzeleloose1;   //!
  TBranch        *b_dzeleloose2;   //!
  TBranch        *b_invMasseleloose1g1;   //!
  TBranch        *b_invMasseleloose1g2;   //!
  TBranch        *b_invMasseleloose2g1;   //!
  TBranch        *b_invMasseleloose2g2;   //!
  TBranch        *b_oEmoPeleloose1;   //!
  TBranch        *b_oEmoPeleloose2;   //!
  TBranch        *b_mvanotrigeleloose1;   //!
  TBranch        *b_mvanotrigeleloose2;   //!
  TBranch        *b_mvatrigeleloose1;   //!
  TBranch        *b_mvatrigeleloose2;   //!
  TBranch        *b_matchconveleloose1;   //!
  TBranch        *b_matchconveleloose2;   //!
  TBranch        *b_chHadIso03eleloose1;   //!
  TBranch        *b_chHadIso03eleloose2;   //!
  TBranch        *b_nHadIso03eleloose1;   //!
  TBranch        *b_nHadIso03eleloose2;   //!
  TBranch        *b_photIso03eleloose1;   //!
  TBranch        *b_photIso03eleloose2;   //!
  TBranch        *b_ptelenontr801;   //!
  TBranch        *b_ptelenontr802;   //!
  TBranch        *b_etaelenontr801;   //!
  TBranch        *b_etaelenontr802;   //!
  TBranch        *b_phielenontr801;   //!
  TBranch        *b_phielenontr802;   //!
  TBranch        *b_eneelenontr801;   //!
  TBranch        *b_eneelenontr802;   //!
  TBranch        *b_ptelenontr901;   //!
  TBranch        *b_ptelenontr902;   //!
  TBranch        *b_etaelenontr901;   //!
  TBranch        *b_etaelenontr902;   //!
  TBranch        *b_phielenontr901;   //!
  TBranch        *b_phielenontr902;   //!
  TBranch        *b_eneelenontr901;   //!
  TBranch        *b_eneelenontr902;   //!
  TBranch        *b_ptmu1;   //!
  TBranch        *b_ptmu2;   //!
  TBranch        *b_etamu1;   //!
  TBranch        *b_etamu2;   //!
  TBranch        *b_phimu1;   //!
  TBranch        *b_phimu2;   //!
  TBranch        *b_enemu1;   //!
  TBranch        *b_enemu2;   //!
  TBranch        *b_pixhitsmu1;   //!
  TBranch        *b_pixhitsmu2;   //!
  TBranch        *b_trkhitsmu1;   //!
  TBranch        *b_trkhitsmu2;   //!
  TBranch        *b_hitsmu1;   //!
  TBranch        *b_hitsmu2;   //!
  TBranch        *b_chi2mu1;   //!
  TBranch        *b_chi2mu2;   //!
  TBranch        *b_matchmu1;   //!
  TBranch        *b_matchmu2;   //!
  TBranch        *b_d0mu1;   //!
  TBranch        *b_d0mu2;   //!
  TBranch        *b_dzmu1;   //!
  TBranch        *b_dzmu2;   //!
  TBranch        *b_chHadmu1;   //!
  TBranch        *b_chHadmu2;   //!
  TBranch        *b_nHadmu1;   //!
  TBranch        *b_nHadmu2;   //!
  TBranch        *b_photmu1;   //!
  TBranch        *b_photmu2;   //!
  TBranch        *b_puptmu1;   //!
  TBranch        *b_puptmu2;   //!
  TBranch        *b_ptmuloose1;   //!
  TBranch        *b_ptmuloose2;   //!
  TBranch        *b_etamuloose1;   //!
  TBranch        *b_etamuloose2;   //!
  TBranch        *b_phimuloose1;   //!
  TBranch        *b_phimuloose2;   //!
  TBranch        *b_enemuloose1;   //!
  TBranch        *b_enemuloose2;   //!
  TBranch        *b_pixhitsmuloose1;   //!
  TBranch        *b_pixhitsmuloose2;   //!
  TBranch        *b_trkhitsmuloose1;   //!
  TBranch        *b_trkhitsmuloose2;   //!
  TBranch        *b_hitsmuloose1;   //!
  TBranch        *b_hitsmuloose2;   //!
  TBranch        *b_chi2muloose1;   //!
  TBranch        *b_chi2muloose2;   //!
  TBranch        *b_matchmuloose1;   //!
  TBranch        *b_matchmuloose2;   //!
  TBranch        *b_d0muloose1;   //!
  TBranch        *b_d0muloose2;   //!
  TBranch        *b_dzmuloose1;   //!
  TBranch        *b_dzmuloose2;   //!
  TBranch        *b_ptmuvloose1;   //!
  TBranch        *b_ptmuvloose2;   //!
  TBranch        *b_etamuvloose1;   //!
  TBranch        *b_etamuvloose2;   //!
  TBranch        *b_phimuvloose1;   //!
  TBranch        *b_phimuvloose2;   //!
  TBranch        *b_enemuvloose1;   //!
  TBranch        *b_enemuvloose2;   //!
  TBranch        *b_pixhitsmuvloose1;   //!
  TBranch        *b_pixhitsmuvloose2;   //!
  TBranch        *b_trkhitsmuvloose1;   //!
  TBranch        *b_trkhitsmuvloose2;   //!
  TBranch        *b_hitsmuvloose1;   //!
  TBranch        *b_hitsmuvloose2;   //!
  TBranch        *b_chi2muvloose1;   //!
  TBranch        *b_chi2muvloose2;   //!
  TBranch        *b_matchmuvloose1;   //!
  TBranch        *b_matchmuvloose2;   //!
  TBranch        *b_d0muvloose1;   //!
  TBranch        *b_d0muvloose2;   //!
  TBranch        *b_dzmuvloose1;   //!
  TBranch        *b_dzmuvloose2;   //!
  TBranch        *b_hasPassedSinglePhot;   //!
  TBranch        *b_hasPassedDoublePhot;   //!
  TBranch        *b_chHadmuloose1;   //!
  TBranch        *b_chHadmuloose2;   //!
  TBranch        *b_nHadmuloose1;   //!
  TBranch        *b_nHadmuloose2;   //!
  TBranch        *b_photmuloose1;   //!
  TBranch        *b_photmuloose2;   //!
  TBranch        *b_puptmuloose1;   //!
  TBranch        *b_puptmuloose2;   //!
  TBranch        *b_chHadmuvloose1;   //!
  TBranch        *b_chHadmuvloose2;   //!
  TBranch        *b_nHadmuvloose1;   //!
  TBranch        *b_nHadmuvloose2;   //!
  TBranch        *b_photmuvloose1;   //!
  TBranch        *b_photmuvloose2;   //!
  TBranch        *b_puptmuvloose1;   //!
  TBranch        *b_puptmuvloose2;   //!
  
  fillPlot2012_radion( const std::string& dataset, const std::string& selectionType, const std::string& bTaggerType="JP" );
  virtual ~fillPlot2012_radion();
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init();
  void setSelectionType( const std::string& selectionType );
  virtual void     finalize();
  virtual std::pair<int,int> myTwoHighestPtJet(vector<int> jets);
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
};

#endif

