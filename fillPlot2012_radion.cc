#include "fillPlot2012_radion.h"

#include "../KinematicFit/DiJetKinFitter.h"

#include "../HelicityLikelihoodDiscriminant/HelicityLikelihoodDiscriminant.h"

using namespace std;

fillPlot2012_radion::fillPlot2012_radion( const std::string& dataset, const std::string& selectionType, const std::string& bTaggerType ) : RedNtpFinalizer( "Radion", dataset ) {
  
  bTaggerType_ = bTaggerType;
  
  setSelectionType(selectionType);
}

fillPlot2012_radion::~fillPlot2012_radion() {
  
  outFile_->Close();
  
  if (!tree_) return;
  delete tree_->GetCurrentFile();
}

double fillPlot2012_radion::delta_phi(double phi1, double phi2) {

  double dphi = TMath::Abs(phi1 - phi2);
  return (dphi <= TMath::Pi())? dphi : TMath::TwoPi() - dphi;
}

void fillPlot2012_radion::finalize() {
  
  this->Init();
  
  std::string fullFlags = selectionType_ + "_";
  fullFlags+=bTaggerType_;
  this->set_flags(fullFlags); 
  this->createOutputFile();

  outFile_->cd();


  // ------------------------------------------------------
  // histograms for pu id studies

  TH1F *myDenomEtaAss  = new TH1F("myDenomEtaAss",  "myDenomEtaAss",  41,-5.,5.);
  myDenomEtaAss->Sumw2();
  TH1F *myDenomEtaNAss = new TH1F("myDenomEtaNAss", "myDenomEtaNAss", 41,-5.,5.);
  myDenomEtaNAss->Sumw2();
  TH1F *myNumEtaAss    = new TH1F("myNumEtaAss",    "myNumEtaAss",    41,-5.,5.);
  myNumEtaAss->Sumw2();
  TH1F *myNumEtaNAss   = new TH1F("myNumEtaNAss",   "myNumEtaNAss",   41,-5.,5.);
  myNumEtaNAss->Sumw2();
  //
  TH1F *myDenomPtAss  = new TH1F("myDenomPtAss",  "myDenomPtAss",  25,25.,100.);
  myDenomPtAss->Sumw2();
  TH1F *myDenomPtNAss = new TH1F("myDenomPtNAss", "myDenomPtNAss", 25,25.,100.);
  myDenomPtNAss->Sumw2();
  TH1F *myNumPtAss    = new TH1F("myNumPtAss",    "myNumPtAss",    25,25.,100.);
  myNumPtAss->Sumw2();
  TH1F *myNumPtNAss   = new TH1F("myNumPtNAss",   "myNumPtNAss",   25,25.,100.);
  myNumPtNAss->Sumw2();
  // 
  TH1F *myAssEta    = new TH1F("myAssEta",    "myAssEta",    41,-5.,5.);
  myAssEta->Sumw2();
  TH1F *myNAssEta   = new TH1F("myNAssEta",   "myNAssEta",   41,-5.,5.);
  myNAssEta->Sumw2();
  TH1F *myAssRMS    = new TH1F("myAssRMS",    "myAssRMS",    50,0.,0.2);
  myAssRMS->Sumw2();
  TH1F *myNAssRMS   = new TH1F("myNAssRMS",   "myNAssRMS",   50,0.,0.2);
  myNAssRMS->Sumw2();
  TH1F *myAssBetaS  = new TH1F("myAssBetaS",  "myAssBetaS",  50,0.,1.);
  myAssBetaS->Sumw2();
  TH1F *myNAssBetaS = new TH1F("myNAssBetaS", "myNAssBetaS", 50,0.,1.);
  myNAssBetaS->Sumw2();
  //
  TH1F *myAssEtalt25    = new TH1F("myAssEtalt25",    "myAssEtalt25",    41,-5.,5.);
  myAssEtalt25->Sumw2();
  TH1F *myNAssEtalt25   = new TH1F("myNAssEtalt25",   "myNAssEtalt25",   41,-5.,5.);
  myNAssEtalt25->Sumw2();
  TH1F *myAssRMSlt25    = new TH1F("myAssRMSlt25",    "myAssRMSlt25",    50,0.,0.2);
  myAssRMSlt25->Sumw2();
  TH1F *myNAssRMSlt25   = new TH1F("myNAssRMSlt25",   "myNAssRMSlt25",   50,0.,0.2);
  myNAssRMSlt25->Sumw2();
  TH1F *myAssBetaSlt25  = new TH1F("myAssBetaSlt25",  "myAssBetaSlt25",  50,0.,1.);
  myAssBetaSlt25->Sumw2();
  TH1F *myNAssBetaSlt25 = new TH1F("myNAssBetaSlt25", "myNAssBetaSlt25", 50,0.,1.);
  myNAssBetaSlt25->Sumw2();
  //
  TH1F *myAssEtalt25_highRho    = new TH1F("myAssEtalt25_highRho",    "myAssEtalt25_highRho",    41,-5.,5.);
  myAssEtalt25_highRho->Sumw2();
  TH1F *myNAssEtalt25_highRho   = new TH1F("myNAssEtalt25_highRho",   "myNAssEtalt25_highRho",   41,-5.,5.);
  myNAssEtalt25_highRho->Sumw2();
  TH1F *myAssRMSlt25_highRho    = new TH1F("myAssRMSlt25_highRho",    "myAssRMSlt25_highRho",    50,0.,0.2);
  myAssRMSlt25_highRho->Sumw2();
  TH1F *myNAssRMSlt25_highRho   = new TH1F("myNAssRMSlt25_highRho",   "myNAssRMSlt25_highRho",   50,0.,0.2);
  myNAssRMSlt25_highRho->Sumw2();
  TH1F *myAssBetaSlt25_highRho  = new TH1F("myAssBetaSlt25_highRho",  "myAssBetaSlt25_highRho",  50,0.,1.);
  myAssBetaSlt25_highRho->Sumw2();
  TH1F *myNAssBetaSlt25_highRho = new TH1F("myNAssBetaSlt25_highRho", "myNAssBetaSlt25_highRho", 50,0.,1.);
  myNAssBetaSlt25_highRho->Sumw2();
  // 
  TH2F *myAssBetaSvsNvtxlt25 = new TH2F("myAssBetaSvsNvtxlt25","myAssBetaSvsNvtxlt25",20,0.,20.,50,0.,1.);
  myAssBetaSvsNvtxlt25->Sumw2();


  // ------------------------------------------------------
  // histograms for kinematic optimization

  TH1D*  h1_njets = new TH1D("njets", "", 11, -0.5, 10.5);
  h1_njets->Sumw2();
  TH1D*  h1_nbjets_loose = new TH1D("nbjets_loose", "", 11, -0.5, 10.5);
  h1_nbjets_loose->Sumw2();
  TH1D*  h1_nbjets_medium = new TH1D("nbjets_medium", "", 11, -0.5, 10.5);
  h1_nbjets_medium->Sumw2();
  TH1D*  h1_nbjets_tight = new TH1D("nbjets_tight", "", 11, -0.5, 10.5);
  h1_nbjets_tight->Sumw2();

  TH1D* h1_ptphot0 = new TH1D("ptphot0", "", 100, 0., 200.);
  h1_ptphot0->Sumw2();
  TH1D* h1_ptphot1 = new TH1D("ptphot1", "", 100, 0., 200.);
  h1_ptphot1->Sumw2();

  TH1D* h1_runptphot0 = new TH1D("runptphot0", "", 100, 0., 200.);
  h1_runptphot0->Sumw2();

  TH1D* h1_mgg_preselG = new TH1D("mgg_preselG", "", 80, 100., 180.);
  h1_mgg_preselG->Sumw2();
  TH1D* h1_mgg_preselJ = new TH1D("mgg_preselJ", "", 80, 100., 180.);
  h1_mgg_preselJ->Sumw2();
  TH1D* h1_mjj_preselJ = new TH1D("mjj_preselJ", "", 200, 0., 500.);
  h1_mjj_preselJ->Sumw2();

  TH1D* h1_ptjet0 = new TH1D("ptjet0", "", 60, 0., 400.);
  h1_ptjet0->Sumw2();
  TH1D* h1_runptjet0 = new TH1D("runptjet0", "", 60, 0., 400.);
  h1_runptjet0->Sumw2();
  TH1D* h1_ptjet1 = new TH1D("ptjet1", "", 30, 0., 200.);
  h1_ptjet1->Sumw2();
  TH1D* h1_etajet0 = new TH1D("etajet0", "", 30, -3., 3.);
  h1_etajet0->Sumw2();
  TH1D* h1_etajet1 = new TH1D("etajet1", "", 30, -3., 3.);
  h1_etajet1->Sumw2();

  TH1D* h1_kinfit_chiSquareProbH = new TH1D("kinfit_chiSquareProbH", "", 1000, 0., 1.0001);
  h1_kinfit_chiSquareProbH->Sumw2();

  TH1D* h1_mjj_0btag = new TH1D("mjj_0btag", "", 200, 0., 500.);
  h1_mjj_0btag->Sumw2();
  TH1D* h1_mjj_1btag = new TH1D("mjj_1btag", "", 200, 0., 500.);
  h1_mjj_1btag->Sumw2();
  TH1D* h1_mjj_2btag = new TH1D("mjj_2btag", "", 200, 0., 500.);
  h1_mjj_2btag->Sumw2();

  TH1D* h1_mggjj = new TH1D("mggjj", "", 100, 0., 1000.);
  h1_mggjj->Sumw2();
  TH1D* h1_mggjj_0btag = new TH1D("mggjj_0btag", "", 100, 0., 1000.);
  h1_mggjj_0btag->Sumw2();
  TH1D* h1_mggjj_1btag = new TH1D("mggjj_1btag", "", 100, 0., 1000.);
  h1_mggjj_1btag->Sumw2();
  TH1D* h1_mggjj_2btag = new TH1D("mggjj_2btag", "", 100, 0., 1000.);
  h1_mggjj_2btag->Sumw2();

  TH1D* h1_mgg_0btag = new TH1D("mgg_0btag", "", 80, 100., 180.);
  h1_mgg_0btag->Sumw2();
  TH1D* h1_mgg_1btag = new TH1D("mgg_1btag", "", 80, 100., 180.);
  h1_mgg_1btag->Sumw2();
  TH1D* h1_mgg_2btag = new TH1D("mgg_2btag", "", 80, 100., 180.);
  h1_mgg_2btag->Sumw2();

  TH1D* h1_mgg_0btag_ebeb = new TH1D("mgg_0btag_ebeb", "", 80, 100., 180.);
  h1_mgg_0btag_ebeb->Sumw2();
  TH1D* h1_mgg_1btag_ebeb = new TH1D("mgg_1btag_ebeb", "", 80, 100., 180.);
  h1_mgg_1btag_ebeb->Sumw2();
  TH1D* h1_mgg_2btag_ebeb = new TH1D("mgg_2btag_ebeb", "", 80, 100., 180.);
  h1_mgg_2btag_ebeb->Sumw2();

  TH1D* h1_mgg_0btag_nebeb = new TH1D("mgg_0btag_nebeb", "", 80, 100., 180.);
  h1_mgg_0btag_nebeb->Sumw2();
  TH1D* h1_mgg_1btag_nebeb = new TH1D("mgg_1btag_nebeb", "", 80, 100., 180.);
  h1_mgg_1btag_nebeb->Sumw2();
  TH1D* h1_mgg_2btag_nebeb = new TH1D("mgg_2btag_nebeb", "", 80, 100., 180.);
  h1_mgg_2btag_nebeb->Sumw2();

  TH1D*  h1_ptDiphot = new TH1D("ptDiphot", "", 100, 0., 500.);
  h1_ptDiphot->Sumw2();
  TH1D*  h1_etaDiphot = new TH1D("etaDiphot", "", 40, -10., 10.);
  h1_etaDiphot->Sumw2();
  TH1D* h1_deltaEtaDiphot = new TH1D("deltaEtaDiphot", "", 40, -5., 5.);
  h1_deltaEtaDiphot->Sumw2();

  TH1D*  h1_deltaPhi = new TH1D("deltaPhi", "", 100, 0., 3.1416);
  h1_deltaPhi->Sumw2();
  TH1D*  h1_deltaEta = new TH1D("deltaEta", "", 50, -5., 5.);
  h1_deltaEta->Sumw2();
  TH1D*  h1_ptDijet = new TH1D("ptDijet", "", 100, 0., 500.);
  h1_ptDijet->Sumw2();
  TH1D*  h1_etaDijet = new TH1D("etaDijet", "", 40, -10., 10.);
  h1_etaDijet->Sumw2();
  TH1D*  h1_ptRatio = new TH1D("ptRatio", "", 100, 0., 3.);
  h1_ptRatio->Sumw2();
  TH1D*  h1_ptDifference = new TH1D("ptDifference", "", 100, -200., 200.);
  h1_ptDifference->Sumw2();
  TH1D* h1_zeppen = new TH1D("zeppen", "", 100, -8., 8.);
  h1_zeppen->Sumw2();

  TH1D* h1_deltaPhiJets = new TH1D("deltaPhiJets", "", 50, -5., 5.);
  h1_deltaPhiJets->Sumw2();
  TH1D* h1_deltaEtaJets = new TH1D("deltaEtaJets", "", 40, -5., 5.);
  h1_deltaEtaJets->Sumw2();
  TH1D* h1_deltaFabsEtaJets = new TH1D("deltaFabsEtaJets", "", 50, -5., 5.);
  h1_deltaFabsEtaJets->Sumw2();

  TH1D*  h1_deltaPhi_kinfit = new TH1D("deltaPhi_kinfit", "", 100, 0., 3.1416);
  h1_deltaPhi_kinfit->Sumw2();
  TH1D*  h1_ptDijet_kinfit = new TH1D("ptDijet_kinfit", "", 100, 0., 500.);
  h1_ptDijet_kinfit->Sumw2();
  TH1D*  h1_ptRatio_kinfit = new TH1D("ptRatio_kinfit", "", 100, 0., 3.);
  h1_ptRatio_kinfit->Sumw2();
  TH1D*  h1_ptDifference_kinfit = new TH1D("ptDifference_kinfit", "", 100, -200., 200.);
  h1_ptDifference_kinfit->Sumw2();
  TH1D* h1_zeppen_kinfit = new TH1D("zeppen_kinfit", "", 100, -8., 8.);
  h1_zeppen_kinfit->Sumw2();

  TH1D* h1_deltaEtaJets_kinfit = new TH1D("deltaEtaJets_kinfit", "", 50, -5., 5.);
  h1_deltaEtaJets_kinfit->Sumw2();
  TH1D* h1_deltaFabsEtaJets_kinfit = new TH1D("deltaFabsEtaJets_kinfit", "", 50, -5., 5.);
  h1_deltaFabsEtaJets_kinfit->Sumw2();

  TH1D* h1_cosTheta1 = new TH1D("cosTheta1", "", 50, -1.0001, 1.0001);
  h1_cosTheta1->Sumw2();
  TH1D* h1_cosTheta2 = new TH1D("cosTheta2", "", 25, 0., 1.0001);
  h1_cosTheta2->Sumw2();
  TH1D* h1_cosThetaStar = new TH1D("cosThetaStar", "", 50, -1.0001, 1.0001);
  h1_cosThetaStar->Sumw2();
  TH1D* h1_helphi = new TH1D("helphi", "", 100, 0., 3.1416);
  h1_helphi->Sumw2();
  TH1D* h1_helphi1 = new TH1D("helphi1", "", 100, 0., 3.1416);
  h1_helphi1->Sumw2();

  TH1D* h1_cosThetaStar_jets = new TH1D("cosThetaStar_jets", "", 50, -1.0001, 1.0001);
  h1_cosThetaStar_jets->Sumw2();

  TH1D* h1_helicityAngle_V = new TH1D("helicityAngle_V", "", 100, -1.0001, 1.0001);
  h1_helicityAngle_V->Sumw2();

  TH1D* h1_mVstar = new TH1D("mVstar", "", 500, 0., 1000.);
  h1_mVstar->Sumw2();
  TH1D* h1_ptVstar = new TH1D("ptVstar", "", 500, 0., 500.);
  h1_ptVstar->Sumw2();
  TH1D* h1_etaVstar = new TH1D("etaVstar", "", 100, -5., 5.);
  h1_etaVstar->Sumw2();
  TH1D* h1_phiVstar = new TH1D("phiVstar", "", 100, 0., 3.1416);
  h1_phiVstar->Sumw2();

  TH1D* h1_mVstar_kinfit = new TH1D("mVstar_kinfit", "", 500, 0., 1000.);
  h1_mVstar_kinfit->Sumw2();
  TH1D* h1_ptVstar_kinfit = new TH1D("ptVstar_kinfit", "", 500, 0., 500.);
  h1_ptVstar_kinfit->Sumw2();
  TH1D* h1_etaVstar_kinfit = new TH1D("etaVstar_kinfit", "", 100, -5., 5.);
  h1_etaVstar_kinfit->Sumw2();
  TH1D* h1_phiVstar_kinfit = new TH1D("phiVstar_kinfit", "", 100, 0., 3.1416);
  h1_phiVstar_kinfit->Sumw2();


  // kin fit study
  TH2D* h2_mggjj_vs_mjj_1btag = new TH2D("mggjj_vs_mjj_1btag", "", 100, 0., 300., 100, 100., 600.);
  h2_mggjj_vs_mjj_1btag->Sumw2();
  TH2D* h2_mggjj_vs_mjj_2btag = new TH2D("mggjj_vs_mjj_2btag", "", 100, 0., 300., 100, 100., 600.);
  h2_mggjj_vs_mjj_2btag->Sumw2();

  TH2D* h2_mggjj_vs_mjj_kinfit_1btag = new TH2D("mggjj_vs_mjj_kinfit_1btag", "", 100, 0., 300., 100, 100., 600.);
  h2_mggjj_vs_mjj_kinfit_1btag->Sumw2();
  TH2D* h2_mggjj_vs_mjj_kinfit_2btag = new TH2D("mggjj_vs_mjj_kinfit_2btag", "", 100, 0., 300., 100, 100., 600.);
  h2_mggjj_vs_mjj_kinfit_2btag->Sumw2();

  TH1D* h1_mggjj_kinfit_1btag = new TH1D("mggjj_kinfit_1btag", "", 100, 100., 600.);
  h1_mggjj_kinfit_1btag->Sumw2();
  TH1D* h1_mggjj_kinfit_2btag = new TH1D("mggjj_kinfit_2btag", "", 100, 100., 600.);
  h1_mggjj_kinfit_2btag->Sumw2();

  TH1D* h1_mggjj_nokinfit_1btag = new TH1D("mggjj_nokinfit_1btag", "", 100, 100., 600.);
  h1_mggjj_nokinfit_1btag->Sumw2();
  TH1D* h1_mggjj_nokinfit_2btag = new TH1D("mggjj_nokinfit_2btag", "", 100, 100., 600.);
  h1_mggjj_nokinfit_2btag->Sumw2();



  //-----------------------------------------------------------------------------
  // for the tree with selected events
  float massggnewvtx_t;
  float ptphot1_t, runptphot1_t, ptphot2_t;
  float etaphot1_t, etaphot2_t;
  float ptgg_t, etagg_t, absetagg_t;
  int njets_t;
  float ptcorrjet1_t, runptcorrjet1_t, ptcorrjet2_t;
  float etajet1_t, etajet2_t;
  float deltaphijj_t, deltaetajj_t;
  float invmassjet_t, ptjj_t, etajj_t, massggjj_t, deltaphijjgg_t, deltaetajjgg_t;
  int btagCategory_t, nbjets_loose_t, nbjets_medium_t, nbjets_tight_t;  
  int theCategory_t;
  float zeppen_t; 
  float chiSquareProbH_t, absCosThetaStar_t;
  int nvtx_t;
  int theVertex_t;
  float rhoPF_t;
  int isj1btagged_t, isj2btagged_t;
  float weight;

  TTree* myTrees = new TTree();
  myTrees->SetName("myTrees");
  myTrees->Branch( "run", &run, "run/I" );
  myTrees->Branch( "lumi", &lumi, "lumi/I" );
  myTrees->Branch( "event", &event, "event/I" );
  myTrees->Branch( "massggnewvtx", &massggnewvtx_t, "massggnewvtx_t/F" );
  myTrees->Branch( "ptPhot1", &ptphot1_t, "ptphot1_t/F" );
  myTrees->Branch( "runptPhot1", &runptphot1_t, "runptPhot1_t/F" );
  myTrees->Branch( "ptPhot2", &ptphot2_t, "ptphot2_t/F" );
  myTrees->Branch( "etaPhot1", &etaphot1_t, "etaphot1_t/F" );
  myTrees->Branch( "etaPhot2", &etaphot2_t, "etaphot2_t/F" );
  myTrees->Branch( "ptgg", &ptgg_t, "ptgg_t/F" );
  myTrees->Branch( "etagg", &etagg_t, "etagg_t/F" );
  myTrees->Branch( "absetagg", &absetagg_t, "absetagg_t/F" );
  myTrees->Branch( "njets", &njets_t, "njets_t/I" );
  myTrees->Branch( "runPtCorrJet1", &runptcorrjet1_t, "runptcorrJet1_t/F" );
  myTrees->Branch( "ptCorrJet1", &ptcorrjet1_t, "ptcorrJet1_t/F" );
  myTrees->Branch( "ptCorrJet2", &ptcorrjet2_t, "ptcorrJet2_t/F" );
  myTrees->Branch( "etaJet1", &etajet1_t, "etajet1_t/F" );
  myTrees->Branch( "etaJet2", &etajet2_t, "etajet2_t/F" );
  myTrees->Branch( "deltaphijj", &deltaphijj_t, "deltaphijj_t/F");
  myTrees->Branch( "deltaetajj", &deltaetajj_t, "deltaetajj_t/F");
  myTrees->Branch( "mjj", &invmassjet_t, "invmassjet_t/F" );
  myTrees->Branch( "ptjj", &ptjj_t, "ptjj_t/F" );
  myTrees->Branch( "etajj", &etajj_t, "etajj_t/F" );
  myTrees->Branch( "mggjj", &massggjj_t, "massggjj_t/F" );
  myTrees->Branch( "deltaphiggjj", &deltaphijjgg_t, "deltaphijjgg_t/F");
  myTrees->Branch( "deltaetaggjj", &deltaetajjgg_t, "deltaetajjgg_t/F");
  myTrees->Branch( "btagCategory", &btagCategory_t, "btagCategory_t/I" );
  myTrees->Branch( "theCategory",  &theCategory_t,  "theCategory_t/I" );
  myTrees->Branch( "nbjets_loose",  &nbjets_loose_t,  "nbjets_loose_t/I" );
  myTrees->Branch( "nbjets_medium", &nbjets_medium_t, "nbjets_medium_t/I" );
  myTrees->Branch( "nbjets_tight",  &nbjets_tight_t,  "nbjets_tight_t/I" );
  myTrees->Branch( "zeppen", &zeppen_t, "zeppen_t/F" );
  myTrees->Branch( "chiSquareProbH", &chiSquareProbH_t, "chiSquareProbH_t/F" );
  myTrees->Branch( "absCosThetaStar", &absCosThetaStar_t, "absCosThetaStar_t/F" );
  myTrees->Branch( "nvtx", &nvtx_t, "nvtx_t/I" );
  myTrees->Branch( "vertex", &theVertex_t, "theVertex_t/I" );
  myTrees->Branch( "rhoPF", &rhoPF_t, "rhoPF_t/F" );
  myTrees->Branch( "isBtagJet1", &isj1btagged_t, "isj1btagged_t/I" );
  myTrees->Branch( "isBtagJet2", &isj2btagged_t, "isj2btagged_t/I" );
  myTrees->Branch( "weight", &weight, "weight/F" );

  // for central limits
  // TTree* TCVARS = new TTree("TCVARS","two photon two jet selection");
  // TCVARS->SetName("TCVARS");
  // TCVARS->Branch( "mgg",          &massggnewvtx_t, "massggnewvtx_t/F" );
  // TCVARS->Branch( "mjj",          &invmassjet_t,   "invmassjet_t/F" );
  // TCVARS->Branch( "mtot",         &massggjj_t,     "massggjj_t/F" );
  // TCVARS->Branch( "cut_based_ct", &theCategory_t,  "theCategory_t/I" );
  // TCVARS->Branch( "evWeight",     &weight,         "weight/F" );


  // ------------------------------------------------------
  // for the kinematic fits, assuming the two jets come from a 125 GeV Higgs 
  float Hmass = 125.;
  DiJetKinFitter* fitter_jetsH = new DiJetKinFitter( "fitter_jetsH", "fitter_jets", Hmass );

  // for the helicity angles study
  HelicityLikelihoodDiscriminant *helicityDiscriminator = new HelicityLikelihoodDiscriminant();
  int seed = 110;
  TRandom3* coin = new TRandom3(seed);

  // analysis
  if (tree_ == 0) return;

  Long64_t nentries = tree_->GetEntries();
  Long64_t nbytes = 0, nb = 0;

  if( xSection_!=1. && xSection_!=0. ) {
    std::cout << std::endl;
    std::cout << "-> Cross-Section: "      << xSection_             << " pb" << std::endl;
    std::cout << "-> # Generated Events: " << nGenEvents_           << std::endl;
    std::cout << "-> Event Weight: "       << xSection_/nGenEvents_ << std::endl << std::endl;
  } 

  for (Long64_t jentry=0; jentry<nentries;jentry++) {

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = tree_->GetEntry(jentry);   nbytes += nb;
    
    // event weight for MC events: xsec and PU
    bool isMC = ( run<5 );
    if( isMC ) {
      if( nGenEvents_<=0 ) {
	std::cout << std::endl << "-> There must be a problem: nGenEvents is not positive!! Exiting" << std::endl;
	exit(99);
      }
      weight = xSection_ / nGenEvents_ ;

      // pu reWeighting
      if(dopureeventWeight_) weight *= pu_weight;

      // chiara 
      // weight *= 11700.;
    }

    // remove duplicate events - chiara, per il momento non rimuovo niente, ma va fatto sia qcd che gjet che wg che zg

    // ---------------------------------------------------
    // gamma-gamma analysis

    // photons acceptance
    if((TMath::Abs(etascphot1)>1.4442&&TMath::Abs(etascphot1)<1.566)||(TMath::Abs(etascphot2)>1.4442&&TMath::Abs(etascphot2)<1.566)
       || TMath::Abs(etascphot1)>2.5 || TMath::Abs(etascphot2)>2.5) continue;  // acceptance
    
    // photon id
    bool idphot1(0), idphot2(0), pxlphot1(1), pxlphot2(1);
    idphot1 = (idcicpfphot1 >= cicselection);
    idphot2 = (idcicpfphot2 >= cicselection);

    if(cicselection>0) {
      if(!(idphot1)) continue;
      if(!(idphot2)) continue;
    } else{
      if(!(idphot1 && pxlphot1)) continue;
      if(!(idphot2 && pxlphot2)) continue;
    }

    // extra cuts on photons: splitting events per photon class if needed
    int myEB=2;
    if( (fabs(etascphot1)>1.4442 || fabs(etascphot2)>1.4442) ) myEB=0;
    if( (fabs(etascphot1)<1.4442 && fabs(etascphot2)<1.4442) ) myEB=1;

    // the two selected photons for the analysis
    TLorentzVector t4phot1, t4phot2;  
    t4phot1.SetPtEtaPhiM(ptphot1,etaphot1,phiphot1,0.);
    t4phot2.SetPtEtaPhiM(ptphot2,etaphot2,phiphot2,0.);
    TLorentzVector t4diPhot;
    t4diPhot.SetPtEtaPhiM( ptggnewvtx, etagg, phigg, massggnewvtx );
    TVector3 t3diPhot;
    t3diPhot.SetPtEtaPhi( ptggnewvtx, etagg, phigg );

    // further cuts on photons -----------------------

    // invariant mass cut on photons
    if (t4diPhot.M()<100 || t4diPhot.M()>180) continue;
    
    // control plots to check the preselection
    h1_ptphot0->Fill(ptphot1, weight);
    h1_ptphot1->Fill(ptphot2, weight);
    h1_runptphot0->Fill(ptphot1*120./massggnewvtx, weight);

    // photons pt cuts 
    if(ptphot1<ptphot1cut* massggnewvtx/120.) continue;         // pt first photon
    if(ptphot2<ptphot2cut)                    continue;         // pt second photon

    // control plots after preselection on photons (pT, acceptance and ID)
    h1_mgg_preselG->Fill( massggnewvtx, weight );


    // --------------------------------------------------------------------------------
    // jets
    
    /*  
    // two highest pT jets, no btagging, no jetID. This is just for PU jetID studies
    vector<int> v_allJets;
    for (int ij=0; ij<njets; ij++) { 
    if ( ptcorrjet[ij] < ptjetacccut ) continue;   
    if ( fabs(etajet[ij])>4.7 )        continue;   
    v_allJets.push_back(ij);
    }

    if (v_allJets.size()<2) continue;   
    std::pair<int,int> allJets_highPt_noId = myTwoHighestPtJet(v_allJets);
    
    // PU jetID studies
    int jet1noId      = allJets_highPt_noId.first;
    int jet2noId      = allJets_highPt_noId.second;    
    float etaJ1noId   = etajet[jet1noId];
    float etaJ2noId   = etajet[jet2noId];
    float rmsJ1noId   = rmsjet[jet1noId]; 
    float rmsJ2noId   = rmsjet[jet2noId]; 
    float betaSJ1noId = betastarjet[jet1noId];
    float betaSJ2noId = betastarjet[jet2noId];

    // control plots
    if (assjet[jet1noId])  { 
    myAssEta   -> Fill(etaJ1noId,weight); 
    myAssRMS   -> Fill(rmsJ1noId,weight);
    myAssBetaS -> Fill(betaSJ1noId,weight); 
    if (fabs(etaJ1noId)<2.5) {
    myAssEtalt25   -> Fill(etaJ1noId,weight); 
    myAssRMSlt25   -> Fill(rmsJ1noId,weight);
    myAssBetaSlt25 -> Fill(betaSJ1noId,weight); 
    myAssBetaSvsNvtxlt25 -> Fill(nvtx,betaSJ1noId,weight);
    if (rhoPF>15) {
    myAssEtalt25_highRho   -> Fill(etaJ1noId,weight); 
    myAssRMSlt25_highRho   -> Fill(rmsJ1noId,weight);
    myAssBetaSlt25_highRho -> Fill(betaSJ1noId,weight); 
    }
    }
    }
    if (assjet[jet2])  { 
    myAssEta   -> Fill(etaJ2noId,weight); 
    myAssRMS   -> Fill(rmsJ2noId,weight);
    myAssBetaS -> Fill(betaSJ2noId,weight); 
    if (fabs(etaJ2noId)<2.5) {
    myAssEtalt25   -> Fill(etaJ2noId,weight); 
    myAssRMSlt25   -> Fill(rmsJ2noId,weight);
    myAssBetaSlt25 -> Fill(betaSJ2noId,weight); 
    myAssBetaSvsNvtxlt25 -> Fill(nvtx,betaSJ2noId,weight);
    if (rhoPF>15) {
    myAssEtalt25_highRho   -> Fill(etaJ2noId,weight); 
    myAssRMSlt25_highRho   -> Fill(rmsJ2noId,weight);
    myAssBetaSlt25_highRho -> Fill(betaSJ2noId,weight); 
    }
    }
    }
    if (!assjet[jet1]) { 
    myNAssEta   -> Fill(etaJ1noId,weight); 
    myNAssRMS   -> Fill(rmsJ1noId,weight);
    myNAssBetaS -> Fill(betaSJ1noId,weight); 
    if (fabs(etaJ1noId)<2.5) {
    myNAssEtalt25   -> Fill(etaJ1noId,weight); 
    myNAssRMSlt25   -> Fill(rmsJ1noId,weight);
    myNAssBetaSlt25 -> Fill(betaSJ1noId,weight); 
    if (rhoPF>15) {
    myNAssEtalt25_highRho   -> Fill(etaJ1noId,weight); 
    myNAssRMSlt25_highRho   -> Fill(rmsJ1noId,weight);
    myNAssBetaSlt25_highRho -> Fill(betaSJ1noId,weight); 
    }
    }
    }
    if (!assjet[jet2]) { 
    myNAssEta   -> Fill(etaJ2noId,weight);
    myNAssRMS   -> Fill(rmsJ2noId,weight);
    myNAssBetaS -> Fill(betaSJ2noId,weight); 
    if (fabs(etaJ2noId)<2.5) {
    myNAssEtalt25   -> Fill(etaJ2noId,weight); 
    myNAssRMSlt25   -> Fill(rmsJ2noId,weight);
    myNAssBetaSlt25 -> Fill(betaSJ2noId,weight); 
    if (rhoPF>15) {
    myNAssEtalt25_highRho   -> Fill(etaJ2noId,weight); 
    myNAssRMSlt25_highRho   -> Fill(rmsJ2noId,weight);
    myNAssBetaSlt25_highRho -> Fill(betaSJ2noId,weight); 
    }
    }
    } 
    
    // control efficiency plots
    if (assjet[jet1]) { 
    myDenomEtaAss -> Fill(etaJ1noId,weight); 
    myDenomPtAss  -> Fill(ptJ1noId,weight); 
    if (passCutBasedJetId(jet1noId)) { 
    myNumEtaAss -> Fill(etaJ1noId,weight); 
    myNumPtAss  -> Fill(ptJ1noId,weight); 
    }
    }
    if (!assjet[jet1]) { 
    myDenomEtaNAss -> Fill(etaJ1noId,weight); 
    myDenomPtNAss  -> Fill(ptJ1noId,weight); 
    if (passCutBasedJetId(jet1noId)) { 
    myNumEtaNAss -> Fill(etaJ1noId,weight); 
    myNumPtNAss  -> Fill(ptJ1noId,weight); 
    }
    }
    if (assjet[jet2]) { 
    myDenomEtaAss -> Fill(etaJ2noId,weight); 
    myDenomPtAss  -> Fill(ptJ2noId,weight); 
    if (passCutBasedJetId(jet2noId)) { 
    myNumEtaAss -> Fill(etaJ2noId,weight); 
    myNumPtAss  -> Fill(ptJ2noId,weight); 
    }
    }
    if (!assjet[jet2]) { 
    myDenomEtaNAss -> Fill(etaJ2noId,weight); 
    myDenomPtNAss  -> Fill(ptJ2noId,weight); 
    if (passCutBasedJetId(jet2noId)) { 
    myNumEtaNAss -> Fill(etaJ2noId,weight); 
    myNumPtNAss  -> Fill(ptJ2noId,weight); 
    }
    }
    */


    // --------------------------------------------------------
    // kinematic analysis 
    
    // jets, no btagging, passing cut based jetID  
    // now restricting to |eta|<2.5
    vector<int> v_puIdJets;
    for (int ij=0; ij<njets; ij++) { 
      if ( ptcorrjet[ij] < ptjetacccut )   continue;   
      if ( fabs(etajet[ij])>etajetacccut ) continue;   
      if ( !passCutBasedJetId(ij) )        continue;
      v_puIdJets.push_back(ij);
    }

    // jets passing btagging ( + eta/pT/PUid cuts)
    vector<int> v_looseJP,  v_mediumJP,  v_tightJP; 
    vector<int> v_looseCSV, v_mediumCSV, v_tightCSV;
    for (int ij=0; ij<int(v_puIdJets.size());ij++) {   
      int index = v_puIdJets[ij];
      if (btagjprobjet[index]>0.275) v_looseJP.push_back(index);
      if (btagjprobjet[index]>0.545) v_mediumJP.push_back(index);
      if (btagjprobjet[index]>0.790) v_tightJP.push_back(index);
      if (btagcsvjet[index]>0.244)   v_looseCSV.push_back(index);
      if (btagcsvjet[index]>0.679)   v_mediumCSV.push_back(index);
      if (btagcsvjet[index]>0.898)   v_tightCSV.push_back(index);
    }

    // control plots (before any cuts on jets)
    h1_njets->Fill( v_puIdJets.size(), weight );
    if (bTaggerType_=="JP") {
      h1_nbjets_loose ->Fill( v_looseJP.size(),  weight );
      h1_nbjets_medium->Fill( v_mediumJP.size(), weight );
      h1_nbjets_tight ->Fill( v_tightJP.size(),  weight );
    } else if (bTaggerType_=="CSV") {
      h1_nbjets_loose ->Fill( v_looseCSV.size(),  weight );
      h1_nbjets_medium->Fill( v_mediumCSV.size(), weight );
      h1_nbjets_tight ->Fill( v_tightCSV.size(),  weight );
    } else {
      cout << "this btag algo does not exists" << endl;
    }

    // at least 2 preselected jets
    if (v_puIdJets.size()<2) continue;   

    // choice of analysis jets ---------------------

    // highest pT btagged jet
    int jet1btag   = -1;
    bool isJustOne = false;
    if (bTaggerType_=="CSV") {
      if (v_mediumCSV.size()>=1) jet1btag = myHighestPtJet(v_mediumCSV);
      if (v_mediumCSV.size()==1) isJustOne = true;
    } else if (bTaggerType_=="JP") {
      if (v_mediumJP.size()>=1) jet1btag = myHighestPtJet(v_mediumJP);       
      if (v_mediumJP.size()==1) isJustOne = true;
    } else {
      cout << "this btag algo does not exists" << endl;
    }

    // at least 1 btagged jet
    if ( jet1btag<0 ) continue;  

    // choose the two jets with highest pT(jj), giving preference to btagged ones
    int jet1 = -1;
    int jet2 = -1;
    int isj1btagged = 0;  
    int isj2btagged = 0;  
    float maxPtBtag = -999.;

    if( isJustOne ) {  // =1 btagged jets
      for (int jet=0; jet<(v_puIdJets.size()); jet++) {
	int index = v_puIdJets[jet];
	if (index==jet1btag) continue;
	TLorentzVector t4jet, t4bjet;
	t4jet.SetPtEtaPhiE(ptcorrjet[index],etajet[index],phijet[index],ecorrjet[index]);
	t4bjet.SetPtEtaPhiE(ptcorrjet[jet1btag],etajet[jet1btag],phijet[jet1btag],ecorrjet[jet1btag]); 
	
	TLorentzVector t4_bnotb = t4jet + t4bjet;
	if ( t4_bnotb.Pt() > maxPtBtag) {
	  maxPtBtag = t4_bnotb.Pt();
	  jet1 = jet1btag;
	  jet2 = index;
	  isj1btagged = true;
	  isj2btagged = false;
	}
      }
      
    } else {  // >1 btagged jets
      
      if (bTaggerType_=="CSV") {
	for (int jetA=0; jetA<(v_mediumCSV.size()-1); jetA++) {
	  for (int jetB=jetA+1; jetB<v_mediumCSV.size(); jetB++) {
	    TLorentzVector t4jetA, t4jetB;
	    int indexA = v_mediumCSV[jetA];
	    int indexB = v_mediumCSV[jetB];
	    t4jetA.SetPtEtaPhiE(ptcorrjet[indexA],etajet[indexA],phijet[indexA],ecorrjet[indexA]);
	    t4jetB.SetPtEtaPhiE(ptcorrjet[indexB],etajet[indexB],phijet[indexB],ecorrjet[indexB]);
	    TLorentzVector t4AB = t4jetA + t4jetB;
	    if ( t4AB.Pt() > maxPtBtag) {
	      maxPtBtag = t4AB.Pt();
	      jet1 = indexA;
	      jet2 = indexB;
	      isj1btagged = true;
	      isj2btagged = true;
	    }
	  }
	}
      } else if (bTaggerType_=="JP") {
	for (int jetA=0; jetA<(v_mediumJP.size()-1); jetA++) {
	  for (int jetB=jetA+1; jetB<v_mediumJP.size(); jetB++) {
	    TLorentzVector t4jetA, t4jetB;
	    int indexA = v_mediumJP[jetA];
	    int indexB = v_mediumJP[jetB];
	    t4jetA.SetPtEtaPhiE(ptcorrjet[indexA],etajet[indexA],phijet[indexA],ecorrjet[indexA]);
	    t4jetB.SetPtEtaPhiE(ptcorrjet[indexB],etajet[indexB],phijet[indexB],ecorrjet[indexB]);
	    TLorentzVector t4AB = t4jetA + t4jetB;
	    if ( t4AB.Pt() > maxPtBtag) {
	      maxPtBtag = t4AB.Pt();
	      jet1 = indexA;
	      jet2 = indexB;
	      isj1btagged = true;
	      isj2btagged = true;
	    }
	  }
	}
      } else {
	cout << "this btag algo does not exists" << endl;	  
      }
    }

    if (jet1<0 || jet2<0) cout << "problem! at least 1 jet not correctly selected" << endl;

    // swap jet1 / jet2 if needed to order in pT
    if (ptcorrjet[jet1]<ptcorrjet[jet2]) {
      int myjet        = jet1;
      bool ismyjetbtag = isj1btagged;
      jet1 = jet2;
      jet2 = myjet;
      isj1btagged = isj2btagged;
      isj2btagged = ismyjetbtag;
    }

    // selected analysis jets 
    TLorentzVector t4jet1, t4jet2;
    t4jet1.SetPtEtaPhiE(ptcorrjet[jet1],etajet[jet1],phijet[jet1],ecorrjet[jet1]);
    t4jet2.SetPtEtaPhiE(ptcorrjet[jet2],etajet[jet2],phijet[jet2],ecorrjet[jet2]);
    TLorentzVector t4diJet = t4jet1 + t4jet2;
    
    // invariant mass cut on jets 
    float invMassJJ = t4diJet.M();
    if( t4diJet.M()<0. || t4diJet.M()>300. )  continue;

    // invariant mass plots after the two preselection on jets and photons
    h1_mgg_preselJ->Fill( massggnewvtx, weight );    
    h1_mjj_preselJ->Fill( t4diJet.M(), weight );

    


    // -------------------------------------------------------------------------
    // jet related kinematic variables to further study the jet selection
    float ptJ1       = t4jet1.Pt();
    float ptJ2       = t4jet2.Pt();
    float etaJ1      = t4jet1.Eta();
    float etaJ2      = t4jet2.Eta();
    float deltaEtaJJ = t4jet1.Eta() - t4jet2.Eta();
    float deltaPhiJJ = t4jet1.DeltaPhi(t4jet2);
    float ptJJ       = t4diJet.Pt();
    float etaJJ      = t4diJet.Eta();

    // jet/photon related kinematic variables
    float zeppen        = t4diPhot.Eta() - 0.5*( etaJ1 + etaJ2); 
    float deltaPhi_ggjj = t4diPhot.DeltaPhi(t4diJet);
    float deltaEta_ggjj = t4diPhot.Eta() - t4diJet.Eta();

    // further diphoton variables
    float deltaEtaGG = t4phot1.Eta() - t4phot2.Eta();

    // radion 4vector
    TLorentzVector t4Radion = t4diJet + t4diPhot;
    float radMass = t4Radion.M();

    // control plots on basic kinematic variables after the two jets and two photons are preselected

    // jets
    h1_ptjet0->Fill( ptJ1, weight );
    h1_ptjet1->Fill( ptJ2, weight );
    h1_runptjet0->Fill(ptJ1*120./t4diJet.M(), weight);

    h1_etajet0->Fill( etaJ1, weight );
    h1_etajet1->Fill( etaJ2, weight );

    h1_ptDijet          -> Fill( ptJJ,  weight );
    h1_etaDijet         -> Fill( etaJJ,  weight );
    h1_deltaEtaJets     -> Fill( deltaEtaJJ, weight );
    h1_deltaFabsEtaJets -> Fill( fabs(deltaEtaJJ), weight );
    h1_deltaPhiJets     -> Fill( deltaPhiJJ, weight );

    // jets - gammas
    h1_zeppen       -> Fill( zeppen, weight);
    h1_deltaPhi     -> Fill( deltaPhi_ggjj, weight );
    h1_deltaEta     -> Fill( deltaEta_ggjj, weight );  
    h1_ptRatio      -> Fill( t4diJet.Pt()/t4diPhot.Pt(), weight );
    h1_ptDifference -> Fill( t4diJet.Pt()-t4diPhot.Pt(), weight );

    // gammas
    h1_ptDiphot       -> Fill( t4diPhot.Pt(), weight );
    h1_etaDiphot      -> Fill( t4diPhot.Eta(), weight );
    h1_deltaEtaDiphot -> Fill( deltaEtaGG, weight );

    
    // ------------------------------------------------------------
    // more refined analyses

    // perform two kinfits
    std::pair<TLorentzVector,TLorentzVector> jets_kinfitH = fitter_jetsH->fit(t4jet1, t4jet2);
    float chiSquareProbH = TMath::Prob(fitter_jetsH->getS(), fitter_jetsH->getNDF());

    // helicity angles
    HelicityLikelihoodDiscriminant::HelicityAngles hangles;
    if( coin->Uniform(1.)<0.5 ) hangles = helicityDiscriminator->computeHelicityAngles(t4phot1, t4phot2, t4jet1, t4jet2);
    else                        hangles = helicityDiscriminator->computeHelicityAngles(t4phot1, t4phot2, t4jet1, t4jet2);

    float cosThetaStar = hangles.helCosThetaStar;
    h1_cosThetaStar -> Fill( hangles.helCosThetaStar, weight );
    h1_cosTheta1    -> Fill( hangles.helCosTheta1, weight );
    h1_cosTheta2    -> Fill( hangles.helCosTheta2, weight );
    h1_helphi       -> Fill( hangles.helPhi, weight );
    h1_helphi1      -> Fill( hangles.helPhi1, weight );

    // boost stuff in the radion frame                                                                                    
    TLorentzVector Vstar_Vstar(t4Radion);
    Vstar_Vstar.Boost(-t4Radion.BoostVector());
    TLorentzVector V_Vstar(t4diJet);
    V_Vstar.Boost(-t4Radion.BoostVector());

    // boost stuff in the V_Vstar frame:                                                                             
    TLorentzVector Vstar_V(Vstar_Vstar);
    Vstar_V.Boost(-V_Vstar.BoostVector());

    TLorentzVector jet1_V, jet2_V;
    // randomize:                                                                                                    
    if( coin->Uniform(1.)<0.5) {
      jet1_V = t4jet1;
      jet2_V = t4jet2;
    } else {
      jet1_V = t4jet2;
      jet2_V = t4jet1;
    }
    jet1_V.Boost(-t4diJet.BoostVector());
    jet2_V.Boost(-t4diJet.BoostVector());

    TVector3 v3_jet1_V  = jet1_V.Vect();
    TVector3 v3_jet2_V  = jet2_V.Vect();
    TVector3 v3_Vstar_V = Vstar_V.Vect();
    TVector3 v3_V_Vstar = V_Vstar.Vect();

    float cosThetaStar_jets = cos( v3_jet1_V.Angle(v3_V_Vstar) );
    float helicityAngle_V   = sin( v3_jet1_V.Angle(v3_Vstar_V) );

    // control plots
    h1_cosThetaStar_jets -> Fill( cosThetaStar_jets, weight );
    h1_helicityAngle_V   -> Fill( helicityAngle_V, weight );
    h1_mVstar            -> Fill( t4Radion.M(), weight );
    h1_ptVstar           -> Fill( t4Radion.Pt(), weight );
    h1_etaVstar          -> Fill( t4Radion.Eta(), weight );
    h1_phiVstar          -> Fill( t4Radion.Phi(), weight );

    
    // refitted jets
    TLorentzVector jet1_kinfit, jet2_kinfit;
    jet1_kinfit = jets_kinfitH.first;
    jet2_kinfit = jets_kinfitH.second;
    TLorentzVector dijet_kinfit = jet1_kinfit + jet2_kinfit;

    // control plots after kin fit
    h1_kinfit_chiSquareProbH->Fill( chiSquareProbH, weight );

    float zeppen_kinfit = t3diPhot.Eta() - 0.5*( jet1_kinfit.Eta() + jet2_kinfit.Eta() );
    h1_zeppen_kinfit->Fill( zeppen_kinfit, weight );
 
    float deltaphi_kinfit = fabs(dijet_kinfit.DeltaPhi(t4diPhot));
    h1_deltaPhi_kinfit         -> Fill( deltaphi_kinfit, weight );
    h1_ptDijet_kinfit          -> Fill( dijet_kinfit.Pt(), weight );
    h1_ptRatio_kinfit          -> Fill( dijet_kinfit.Pt()/t4diPhot.Pt(), weight );
    h1_ptDifference_kinfit     -> Fill( dijet_kinfit.Pt()-t4diPhot.Pt(), weight );
    h1_deltaEtaJets_kinfit     -> Fill( jet1_kinfit.Eta()-jet2_kinfit.Eta(), weight );
    h1_deltaFabsEtaJets_kinfit -> Fill( fabs(jet1_kinfit.Eta())-fabs(jet2_kinfit.Eta()), weight );

    TLorentzVector Vstar_kinfit = dijet_kinfit + t4diPhot;
    h1_mVstar_kinfit   -> Fill( Vstar_kinfit.M(),   weight );
    h1_ptVstar_kinfit  -> Fill( Vstar_kinfit.Pt(),  weight );
    h1_etaVstar_kinfit -> Fill( Vstar_kinfit.Eta(), weight );
    h1_phiVstar_kinfit -> Fill( Vstar_kinfit.Phi(), weight );


    // -------------------------------------------------------------
    // counting the number of bjets to categorize the events
    int btagCategory = -1;
    if ( isj1btagged==1 && isj2btagged==1 ) btagCategory = 2; 
    if ( (isj1btagged==1 && isj2btagged==0) || (isj1btagged==0 && isj2btagged==1) ) btagCategory = 1; 

    // for test: only EBEB events
    // if (!myEB) continue;

    // checking if photons are in EB or in EE to categorize the events    
    int myR9=-1;
    if(r9phot1>.94 && r9phot2>.94) myR9 = 1;
    if(r9phot1<.94 || r9phot2<.94) myR9 = 0;

    // total category combining photons / jets
    int theCategory = -1;
    if (myR9==1 && btagCategory==2) theCategory = 0;
    if (myR9==1 && btagCategory==1) theCategory = 1;
    if (myR9==0 && btagCategory==2) theCategory = 2;
    if (myR9==0 && btagCategory==1) theCategory = 3;

    // to compare with an only photon-based analysis 
    // if (myR9==1 && myEB)  theCategory = 0;
    // if (myR9==0 && myEB)  theCategory = 1;
    // if (myR9==1 && !myEB) theCategory = 2;
    // if (myR9==0 && !myEB) theCategory = 3;

    // -------------------------------------------------------------
    // here we apply further cuts according to the btag category
    /*
    if (btagCategory<1 || btagCategory>2) continue;

    if (btagCategory==1) {
      if ( ptdiphotcut_1bj>0  && t4diPhot.Pt()<ptdiphotcut_1bj ) continue;
      if ( etadiphotcut_1bj>0 && fabs(etagg)>etadiphotcut_1bj  ) continue;
      if ( ptphot1cut_1bj>0   && ((120./massggnewvtx)*ptphot1)<ptphot1cut_1bj ) continue;
      if ( ptjet1cut_1bj>0    && ptJ1<ptjet1cut_1bj ) continue;
      if ( costhetascut_1bj>0 && fabs(cosThetaStar)>costhetascut_1bj ) continue;

    } else if (btagCategory==2) {
      if ( ptdiphotcut_2bj>0  && t4diPhot.Pt()<ptdiphotcut_2bj ) continue;
      if ( etadiphotcut_2bj>0 && fabs(etagg)>etadiphotcut_2bj  ) continue;
      if ( ptphot1cut_2bj>0   && ((120./massggnewvtx)*ptphot1)<ptphot1cut_2bj ) continue;
      if ( ptjet1cut_2bj>0    && ptJ1<ptjet1cut_2bj ) continue;
      if ( costhetascut_2bj>0 && fabs(cosThetaStar)>costhetascut_2bj ) continue;
    }
    */

    // -------------------------------------------------------------

    // invariant mass of jets by btag category
    if( btagCategory==0 ) {
      h1_mjj_0btag->Fill( invMassJJ, weight );
    } else if( btagCategory==1 ) {
      h1_mjj_1btag->Fill( invMassJJ, weight );
    } else {
      h1_mjj_2btag->Fill( invMassJJ, weight );
    }

    // invariant mass of photons by btag category
    if( btagCategory==0 ) {
      h1_mgg_0btag -> Fill( massggnewvtx, weight );
      if( myEB ) h1_mgg_0btag_ebeb ->Fill( massggnewvtx, weight );
      else       h1_mgg_0btag_nebeb->Fill( massggnewvtx, weight );
    } else if( btagCategory==1 ) {
      h1_mgg_1btag -> Fill( massggnewvtx, weight );
      if( myEB ) h1_mgg_1btag_ebeb ->Fill( massggnewvtx, weight );
      else       h1_mgg_1btag_nebeb->Fill( massggnewvtx, weight );
    } else if( btagCategory==2 ) {
      h1_mgg_2btag -> Fill( massggnewvtx, weight );
      if( myEB ) h1_mgg_2btag_ebeb ->Fill( massggnewvtx, weight );
      else       h1_mgg_2btag_nebeb->Fill( massggnewvtx, weight );
    }

    // invariant mass of radion by btag category 
    h1_mggjj->Fill( radMass, weight );
    if( btagCategory==0 ) {
      h1_mggjj_0btag->Fill( radMass, weight );
    } else if( btagCategory==1 ) {
      h1_mggjj_1btag->Fill( radMass, weight );
    } else {
      h1_mggjj_2btag->Fill( radMass, weight );
    }

    // ---------------------------------------
    // for a comparison of masses with / wo kin fit
    if( btagCategory==1 ) {
      h2_mggjj_vs_mjj_1btag        -> Fill( invMassJJ, radMass, weight );
      h2_mggjj_vs_mjj_kinfit_1btag -> Fill( invMassJJ, Vstar_kinfit.M(), weight );
      h1_mggjj_kinfit_1btag        -> Fill( Vstar_kinfit.M(), weight );
      h1_mggjj_nokinfit_1btag      -> Fill( radMass, weight );
    } else if (btagCategory==2 ) {
      h2_mggjj_vs_mjj_2btag        -> Fill( invMassJJ, radMass, weight );
      h2_mggjj_vs_mjj_kinfit_2btag -> Fill( invMassJJ, Vstar_kinfit.M(), weight );
      h1_mggjj_kinfit_2btag        -> Fill( Vstar_kinfit.M(), weight );
      h1_mggjj_nokinfit_2btag      -> Fill( radMass, weight );
    }

    // filling the tree for selected events 
    massggnewvtx_t = massggnewvtx;
    ptphot1_t  = ptphot1;
    runptphot1_t = (120./massggnewvtx)*ptphot1;
    ptphot2_t  = ptphot2;
    etaphot1_t = etaphot1;
    etaphot2_t = etaphot2;
    ptgg_t  = ptgg;
    etagg_t = etagg; 
    absetagg_t = fabs(etagg); 
    njets_t = v_puIdJets.size();
    ptcorrjet1_t = ptJ1;
    ptcorrjet2_t = ptJ2;
    runptcorrjet1_t = (120./invMassJJ)*ptJ1;
    etajet1_t = etaJ1;
    etajet2_t = etaJ2;
    deltaphijj_t = deltaPhiJJ;
    deltaetajj_t = deltaEtaJJ;
    invmassjet_t = invMassJJ;
    ptjj_t = ptJJ;
    etajj_t = etaJJ;
    massggjj_t      = radMass;
    deltaphijjgg_t  = deltaPhi_ggjj;
    deltaetajjgg_t  = deltaEta_ggjj;
    btagCategory_t  = btagCategory;
    theCategory_t   = theCategory;
    if (bTaggerType_=="JP") {
      nbjets_loose_t  = v_looseJP.size();
      nbjets_medium_t = v_mediumJP.size();
      nbjets_tight_t  = v_tightJP.size();
    } else if (bTaggerType_=="CSV") {
      nbjets_loose_t  = v_looseCSV.size();
      nbjets_medium_t = v_mediumCSV.size();
      nbjets_tight_t  = v_tightCSV.size();
    } else {
      cout << "this btag algo does not exist" << endl;
    }
    zeppen_t = zeppen;
    chiSquareProbH_t = chiSquareProbH;
    absCosThetaStar_t  = fabs(cosThetaStar);
    nvtx_t      = nvtx;
    theVertex_t = vtxId;
    rhoPF_t  = rhoPF;
    isj1btagged_t = isj1btagged;
    isj2btagged_t = isj2btagged;
    // chiara
    // if( !isMC ) { weight = 1.; }

    myTrees->Fill();

    // TCVARS->Fill();

  } // loop over entries 


  outFile_->cd();
  // TCVARS->Write();  

  myTrees->Write();

  h1_njets->Write();
  h1_nbjets_loose->Write();
  h1_nbjets_medium->Write();
  h1_nbjets_tight->Write();
  
  h1_ptphot0->Write();
  h1_ptphot1->Write();
  h1_runptphot0->Write();

  h1_mgg_preselG->Write();
  h1_mgg_preselJ->Write();
  h1_mjj_preselJ->Write();

  h1_ptjet0->Write();
  h1_ptjet1->Write();
  h1_runptjet0->Write();

  h1_etajet0->Write();
  h1_etajet1->Write();
  
  h1_kinfit_chiSquareProbH->Write();
  
  h1_mjj_0btag->Write();
  h1_mjj_1btag->Write();
  h1_mjj_2btag->Write();
  
  h1_mggjj->Write();
  h1_mggjj_0btag->Write();
  h1_mggjj_1btag->Write();
  h1_mggjj_2btag->Write();

  h1_mgg_0btag->Write();
  h1_mgg_1btag->Write();
  h1_mgg_2btag->Write();

  h1_mgg_0btag_ebeb->Write();
  h1_mgg_1btag_ebeb->Write();
  h1_mgg_2btag_ebeb->Write();

  h1_mgg_0btag_nebeb->Write();
  h1_mgg_1btag_nebeb->Write();
  h1_mgg_2btag_nebeb->Write();

  h1_ptDiphot->Write();
  h1_etaDiphot->Write();
  h1_deltaEtaDiphot->Write();
  
  h1_deltaPhi->Write();
  h1_deltaEta->Write();
  h1_ptDijet->Write();
  h1_etaDijet->Write();
  h1_ptRatio->Write();
  h1_ptDifference->Write();
  
  h1_deltaPhiJets->Write();
  h1_deltaEtaJets->Write();
  h1_deltaFabsEtaJets->Write();
  h1_zeppen->Write();
  
  h1_deltaPhi_kinfit->Write();
  h1_ptDijet_kinfit->Write();
  h1_ptRatio_kinfit->Write();
  h1_ptDifference_kinfit->Write();
  
  h1_deltaEtaJets_kinfit->Write();
  h1_deltaFabsEtaJets_kinfit->Write();
  h1_zeppen_kinfit->Write();
  
  h1_cosTheta1->Write();
  h1_cosTheta2->Write();
  h1_cosThetaStar->Write();
  h1_helphi->Write();
  h1_helphi1->Write();
  
  h1_cosThetaStar_jets->Write();
  h1_helicityAngle_V->Write();
  
  h1_mVstar->Write();
  h1_ptVstar->Write();
  h1_etaVstar->Write();
  h1_phiVstar->Write();
  
  h1_mVstar_kinfit->Write();
  h1_ptVstar_kinfit->Write();
  h1_etaVstar_kinfit->Write();
  h1_phiVstar_kinfit->Write();

  h2_mggjj_vs_mjj_1btag        -> Write();
  h2_mggjj_vs_mjj_kinfit_1btag -> Write();
  h1_mggjj_kinfit_1btag        -> Write();
  h1_mggjj_nokinfit_1btag      -> Write();

  h2_mggjj_vs_mjj_2btag        -> Write();
  h2_mggjj_vs_mjj_kinfit_2btag -> Write();
  h1_mggjj_kinfit_2btag        -> Write();
  h1_mggjj_nokinfit_2btag      -> Write();

  // eff vs # eta
  TH1F *myEffEtaAss = (TH1F*)myNumEtaAss->Clone("myEffEtaAss");
  myEffEtaAss   -> Divide(myDenomEtaAss);
  myNumEtaAss   -> Write();
  myDenomEtaAss -> Write();
  myEffEtaAss   -> Write();
  
  TH1F *myEffEtaNAss = (TH1F*)myNumEtaNAss->Clone("myEffEtaNAss");
  myEffEtaNAss   -> Divide(myDenomEtaNAss);
  myNumEtaNAss   -> Write();
  myDenomEtaNAss -> Write();
  myEffEtaNAss   -> Write();
  
  // eff vs # pt
  TH1F *myEffPtAss = (TH1F*)myNumPtAss->Clone("myEffPtAss");
  myEffPtAss   -> Divide(myDenomPtAss);
  myNumPtAss   -> Write();
  myDenomPtAss -> Write();
  myEffPtAss   -> Write();
  
  TH1F *myEffPtNAss = (TH1F*)myNumPtNAss->Clone("myEffPtNAss");
  myEffPtNAss   -> Divide(myDenomPtNAss);
  myNumPtNAss   -> Write();
  myDenomPtNAss -> Write();
  myEffPtNAss   -> Write();
  
  myNAssEta   -> Write();
  myAssEta    -> Write();
  myNAssRMS   -> Write();
  myAssRMS    -> Write();
  myNAssBetaS -> Write();
  myAssBetaS  -> Write();
  
  myNAssEtalt25   -> Write();
  myAssEtalt25    -> Write();
  myNAssRMSlt25   -> Write();
  myAssRMSlt25    -> Write();
  myNAssBetaSlt25 -> Write();
  myAssBetaSlt25  -> Write();
  
  myNAssEtalt25_highRho   -> Write();
  myAssEtalt25_highRho    -> Write();
  myNAssRMSlt25_highRho   -> Write();
  myAssRMSlt25_highRho    -> Write();
  myNAssBetaSlt25_highRho -> Write();
  myAssBetaSlt25_highRho  -> Write();

  myAssBetaSvsNvtxlt25 -> Write();

} // finalize

Int_t fillPlot2012_radion::GetEntry(Long64_t entry) {

  if (!tree_) return 0;
  return tree_->GetEntry(entry);
}

Long64_t fillPlot2012_radion::LoadTree(Long64_t entry) {

  if (!tree_) return -5;
  Long64_t centry = tree_->LoadTree(entry);
  if (centry < 0) return centry;
  if (tree_->GetTreeNumber() != fCurrent) {
    fCurrent = tree_->GetTreeNumber();
  }
  return centry;
}

void fillPlot2012_radion::Init() {

  // The Init() function is called when the selector needs to initialize                                                                      
  // a new tree or chain. Typically here the branch addresses and branch                                                                      
  // pointers of the tree will be set.                                                                                                        
  // It is normally not necessary to make changes to the generated                                                                            
  // code, but the routine can be extended by the user if needed.                                                                             
  // Init() will be called many times when running on PROOF                                                                                   
  // (once per file to be processed).                                                                                                         

  // Set branch addresses and branch pointers                                                                          
  fCurrent = -1;
  tree_->SetMakeClass(1);

  tree_->SetBranchAddress("run", &run, &b_run);
  tree_->SetBranchAddress("event", &event, &b_event);
  tree_->SetBranchAddress("lumi", &lumi, &b_lumi);
  tree_->SetBranchAddress("H_event", &H_event, &b_H_event);
  tree_->SetBranchAddress("V_event", &V_event, &b_V_event);
  tree_->SetBranchAddress("WH_event", &WH_event, &b_WH_event);
  tree_->SetBranchAddress("ZH_event", &ZH_event, &b_ZH_event);
  tree_->SetBranchAddress("Zbb_event", &Zbb_event, &b_Zbb_event);
  tree_->SetBranchAddress("Vqq_event", &Vqq_event, &b_Vqq_event);
  tree_->SetBranchAddress("rhoPF", &rhoPF, &b_rhoPF);
  tree_->SetBranchAddress("massgg", &massgg, &b_massgg);
  tree_->SetBranchAddress("ptgg", &ptgg, &b_ptgg);
  tree_->SetBranchAddress("ptggnewvtx", &ptggnewvtx, &b_ptggnewvtx);
  tree_->SetBranchAddress("phigg", &phigg, &b_phigg);
  tree_->SetBranchAddress("etagg", &etagg, &b_etagg);
  tree_->SetBranchAddress("massggnewvtx", &massggnewvtx, &b_massggnewvtx);
  tree_->SetBranchAddress("ptphot1", &ptphot1, &b_ptphot1);
  tree_->SetBranchAddress("ptphot2", &ptphot2, &b_ptphot2);
  tree_->SetBranchAddress("deltaRToTrackphot1", &deltaRToTrackphot1, &b_deltaRToTrackphot1);
  tree_->SetBranchAddress("deltaRToTrackphot2", &deltaRToTrackphot2, &b_deltaRToTrackphot2);
  tree_->SetBranchAddress("timephot1", &timephot1, &b_timephot1);
  tree_->SetBranchAddress("timephot2", &timephot2, &b_timephot2);
  tree_->SetBranchAddress("etaphot1", &etaphot1, &b_etaphot1);
  tree_->SetBranchAddress("etaphot2", &etaphot2, &b_etaphot2);
  tree_->SetBranchAddress("phiphot1", &phiphot1, &b_phiphot1);
  tree_->SetBranchAddress("phiphot2", &phiphot2, &b_phiphot2);
  tree_->SetBranchAddress("etascphot1", &etascphot1, &b_etascphot1);
  tree_->SetBranchAddress("etascphot2", &etascphot2, &b_etascphot2);
  tree_->SetBranchAddress("phiscphot1", &phiscphot1, &b_phiscphot1);
  tree_->SetBranchAddress("phiscphot2", &phiscphot2, &b_phiscphot2);
  tree_->SetBranchAddress("E1phot1", &E1phot1, &b_E1phot1);
  tree_->SetBranchAddress("E1phot2", &E1phot2, &b_E1phot2);
  tree_->SetBranchAddress("E9phot1", &E9phot1, &b_E9phot1);
  tree_->SetBranchAddress("E9phot2", &E9phot2, &b_E9phot2);
  tree_->SetBranchAddress("energyErrphot1", &energyErrphot1, &b_energyErrphot1);
  tree_->SetBranchAddress("energyErrphot2", &energyErrphot2, &b_energyErrphot2);
  tree_->SetBranchAddress("energySmearingphot1", &energySmearingphot1, &b_energySmearingphot1);
  tree_->SetBranchAddress("energySmearingphot2", &energySmearingphot2, &b_energySmearingphot2);
  tree_->SetBranchAddress("r9phot1", &r9phot1, &b_r9phot1);
  tree_->SetBranchAddress("r9phot2", &r9phot2, &b_r9phot2);
  tree_->SetBranchAddress("isemEGphot1", &isemEGphot1, &b_isemEGphot1);
  tree_->SetBranchAddress("isemEGphot2", &isemEGphot2, &b_isemEGphot2);
  tree_->SetBranchAddress("promptGamma", &promptGamma, &b_promptGamma);
  tree_->SetBranchAddress("LOGamma", &LOGamma, &b_LOGamma);
  tree_->SetBranchAddress("ISRGamma", &ISRGamma, &b_ISRGamma);
  tree_->SetBranchAddress("FSRGamma", &FSRGamma, &b_FSRGamma);
  tree_->SetBranchAddress("idmvaphot1", &idmvaphot1, &b_idmvaphot1);
  tree_->SetBranchAddress("idmvaphot2", &idmvaphot2, &b_idmvaphot2);
  tree_->SetBranchAddress("idcicphot1", &idcicphot1, &b_idcicphot1);
  tree_->SetBranchAddress("idcicphot2", &idcicphot2, &b_idcicphot2);
  tree_->SetBranchAddress("idcicnoelvetophot1", &idcicnoelvetophot1, &b_idcicnoelvetophot1);
  tree_->SetBranchAddress("idcicnoelvetophot2", &idcicnoelvetophot2, &b_idcicnoelvetophot2);
  tree_->SetBranchAddress("idcicpfphot1", &idcicpfphot1, &b_idcicpfphot1);
  tree_->SetBranchAddress("idcicpfphot2", &idcicpfphot2, &b_idcicpfphot2);
  tree_->SetBranchAddress("idcicpfnoelvetophot1", &idcicpfnoelvetophot1, &b_idcicpfnoelvetophot1);
  tree_->SetBranchAddress("idcicpfnoelvetophot2", &idcicpfnoelvetophot2, &b_idcicpfnoelvetophot2);
  tree_->SetBranchAddress("idelephot1", &idelephot1, &b_idelephot1);
  tree_->SetBranchAddress("idelephot2", &idelephot2, &b_idelephot2);
  tree_->SetBranchAddress("pid_isEMphot1", &pid_isEMphot1, &b_pid_isEMphot1);
  tree_->SetBranchAddress("pid_isEMphot2", &pid_isEMphot2, &b_pid_isEMphot2);
  tree_->SetBranchAddress("pid_haspixelseedphot1", &pid_haspixelseedphot1, &b_pid_haspixelseedphot1);
  tree_->SetBranchAddress("pid_haspixelseedphot2", &pid_haspixelseedphot2, &b_pid_haspixelseedphot2);
  tree_->SetBranchAddress("pid_jurECALphot1", &pid_jurECALphot1, &b_pid_jurECALphot1);
  tree_->SetBranchAddress("pid_jurECALphot2", &pid_jurECALphot2, &b_pid_jurECALphot2);
  tree_->SetBranchAddress("pid_twrHCALphot1", &pid_twrHCALphot1, &b_pid_twrHCALphot1);
  tree_->SetBranchAddress("pid_twrHCALphot2", &pid_twrHCALphot2, &b_pid_twrHCALphot2);
  tree_->SetBranchAddress("pid_HoverEphot1", &pid_HoverEphot1, &b_pid_HoverEphot1);
  tree_->SetBranchAddress("pid_HoverEphot2", &pid_HoverEphot2, &b_pid_HoverEphot2);
  tree_->SetBranchAddress("pid_hlwTrackphot1", &pid_hlwTrackphot1, &b_pid_hlwTrackphot1);
  tree_->SetBranchAddress("pid_hlwTrackphot2", &pid_hlwTrackphot2, &b_pid_hlwTrackphot2);
  tree_->SetBranchAddress("pid_etawidphot1", &pid_etawidphot1, &b_pid_etawidphot1);
  tree_->SetBranchAddress("pid_etawidphot2", &pid_etawidphot2, &b_pid_etawidphot2);
  tree_->SetBranchAddress("pid_hlwTrackNoDzphot1", &pid_hlwTrackNoDzphot1, &b_pid_hlwTrackNoDzphot1);
  tree_->SetBranchAddress("pid_hlwTrackNoDzphot2", &pid_hlwTrackNoDzphot2, &b_pid_hlwTrackNoDzphot2);
  tree_->SetBranchAddress("pid_hasMatchedConvphot1", &pid_hasMatchedConvphot1, &b_pid_hasMatchedConvphot1);
  tree_->SetBranchAddress("pid_hasMatchedConvphot2", &pid_hasMatchedConvphot2, &b_pid_hasMatchedConvphot2);
  tree_->SetBranchAddress("pid_hasMatchedPromptElephot1", &pid_hasMatchedPromptElephot1, &b_pid_hasMatchedPromptElephot1);
  tree_->SetBranchAddress("pid_hasMatchedPromptElephot2", &pid_hasMatchedPromptElephot2, &b_pid_hasMatchedPromptElephot2);
  tree_->SetBranchAddress("pid_sminphot1", &pid_sminphot1, &b_pid_sminphot1);
  tree_->SetBranchAddress("pid_sminphot2", &pid_sminphot2, &b_pid_sminphot2);
  tree_->SetBranchAddress("pid_smajphot1", &pid_smajphot1, &b_pid_smajphot1);
  tree_->SetBranchAddress("pid_smajphot2", &pid_smajphot2, &b_pid_smajphot2);
  tree_->SetBranchAddress("pid_ntrkphot1", &pid_ntrkphot1, &b_pid_ntrkphot1);
  tree_->SetBranchAddress("pid_ntrkphot2", &pid_ntrkphot2, &b_pid_ntrkphot2);
  tree_->SetBranchAddress("pid_ptisophot1", &pid_ptisophot1, &b_pid_ptisophot1);
  tree_->SetBranchAddress("pid_ptisophot2", &pid_ptisophot2, &b_pid_ptisophot2);
  tree_->SetBranchAddress("pid_ntrkcsphot1", &pid_ntrkcsphot1, &b_pid_ntrkcsphot1);
  tree_->SetBranchAddress("pid_ntrkcsphot2", &pid_ntrkcsphot2, &b_pid_ntrkcsphot2);
  tree_->SetBranchAddress("pid_ptisocsphot1", &pid_ptisocsphot1, &b_pid_ptisocsphot1);
  tree_->SetBranchAddress("pid_ptisocsphot2", &pid_ptisocsphot2, &b_pid_ptisocsphot2);
  tree_->SetBranchAddress("pid_ecalisophot1", &pid_ecalisophot1, &b_pid_ecalisophot1);
  tree_->SetBranchAddress("pid_ecalisophot2", &pid_ecalisophot2, &b_pid_ecalisophot2);
  tree_->SetBranchAddress("pid_hcalisophot1", &pid_hcalisophot1, &b_pid_hcalisophot1);
  tree_->SetBranchAddress("pid_hcalisophot2", &pid_hcalisophot2, &b_pid_hcalisophot2);
  tree_->SetBranchAddress("njets", &njets, &b_njets);
  tree_->SetBranchAddress("ecorrjet", ecorrjet, &b_ecorrjet);
  tree_->SetBranchAddress("ptjet", ptjet, &b_ptjet);
  tree_->SetBranchAddress("ptcorrjet", ptcorrjet, &b_ptcorrjet);
  tree_->SetBranchAddress("etajet", etajet, &b_etajet);
  tree_->SetBranchAddress("phijet", phijet, &b_phijet);
  tree_->SetBranchAddress("betajet", betajet, &b_betajet);
  tree_->SetBranchAddress("betastarjet", betastarjet, &b_betastarjet);
  tree_->SetBranchAddress("btagvtxjet", btagvtxjet, &b_btagvtxjet);
  tree_->SetBranchAddress("btagcsvjet", btagcsvjet, &b_btagcsvjet);
  tree_->SetBranchAddress("btagjprobjet", btagjprobjet, &b_btagjprobjet);
  tree_->SetBranchAddress("btagcsvmvajet", btagcsvmvajet, &b_btagcsvmvajet);
  tree_->SetBranchAddress("btagbjprobjet", btagbjprobjet, &b_btagbjprobjet);
  tree_->SetBranchAddress("btagvtxhighpurjet", btagvtxhighpurjet, &b_btagvtxhighpurjet);
  tree_->SetBranchAddress("btagtchejet", btagtchejet, &b_btagtchejet);
  tree_->SetBranchAddress("btagtchpjet", btagtchpjet, &b_btagtchpjet);
  tree_->SetBranchAddress("ptDjet", ptDjet, &b_ptDjet);
  tree_->SetBranchAddress("ptD_QCjet", ptD_QCjet, &b_ptD_QCjet);
  tree_->SetBranchAddress("axis2_QCjet", axis2_QCjet, &b_axis2_QCjet);
  tree_->SetBranchAddress("rmsjet", rmsjet, &b_rmsjet);
  tree_->SetBranchAddress("ntrkjet", ntrkjet, &b_ntrkjet);
  tree_->SetBranchAddress("nneutjet", nneutjet, &b_nneutjet);
  tree_->SetBranchAddress("nChg_QCjet", nChg_QCjet, &b_nChg_QCjet);
  tree_->SetBranchAddress("nNeutral_ptCutjet", nNeutral_ptCutjet, &b_nNeutral_ptCutjet);
  tree_->SetBranchAddress("jetIdSimple_mvajet", jetIdSimple_mvajet, &b_jetIdSimple_mvajet);
  tree_->SetBranchAddress("jetIdFull_mvajet", jetIdFull_mvajet, &b_jetIdFull_mvajet);
  tree_->SetBranchAddress("jetId_dR2Meanjet", jetId_dR2Meanjet, &b_jetId_dR2Meanjet);
  tree_->SetBranchAddress("jetId_betaStarClassicjet", jetId_betaStarClassicjet, &b_jetId_betaStarClassicjet);
  tree_->SetBranchAddress("jetId_frac01jet", jetId_frac01jet, &b_jetId_frac01jet);
  tree_->SetBranchAddress("jetId_frac02jet", jetId_frac02jet, &b_jetId_frac02jet);
  tree_->SetBranchAddress("jetId_frac03jet", jetId_frac03jet, &b_jetId_frac03jet);
  tree_->SetBranchAddress("jetId_frac04jet", jetId_frac04jet, &b_jetId_frac04jet);
  tree_->SetBranchAddress("jetId_frac05jet", jetId_frac05jet, &b_jetId_frac05jet);
  tree_->SetBranchAddress("jetId_betajet", jetId_betajet, &b_jetId_betajet);
  tree_->SetBranchAddress("jetId_betaStarjet", jetId_betaStarjet, &b_jetId_betaStarjet);
  tree_->SetBranchAddress("jetIdCutBased_wpjet", jetIdCutBased_wpjet, &b_jetIdCutBased_wpjet);
  tree_->SetBranchAddress("jetIdSimple_wpjet", jetIdSimple_wpjet, &b_jetIdSimple_wpjet);
  tree_->SetBranchAddress("jetIdFull_wpjet", jetIdFull_wpjet, &b_jetIdFull_wpjet);
  tree_->SetBranchAddress("assjet", assjet, &b_assjet);
  tree_->SetBranchAddress("partPdgIDjet", partPdgIDjet, &b_partPdgIDjet);
  tree_->SetBranchAddress("partMomPdgIDjet", partMomPdgIDjet, &b_partMomPdgIDjet);
  tree_->SetBranchAddress("nvtx", &nvtx, &b_nvtx);
  tree_->SetBranchAddress("vtxId", &vtxId, &b_vtxId);
  tree_->SetBranchAddress("vtxPos_x", &vtxPos_x, &b_vtxPos_x);
  tree_->SetBranchAddress("vtxPos_y", &vtxPos_y, &b_vtxPos_y);
  tree_->SetBranchAddress("vtxPos_z", &vtxPos_z, &b_vtxPos_z);
  tree_->SetBranchAddress("vtxIdMVA", &vtxIdMVA, &b_vtxIdMVA);
  tree_->SetBranchAddress("vtxIdEvtProb", &vtxIdEvtProb, &b_vtxIdEvtProb);
  tree_->SetBranchAddress("diPhotMVA", &diPhotMVA, &b_diPhotMVA);
  tree_->SetBranchAddress("diPhotMVA_vtx0", &diPhotMVA_vtx0, &b_diPhotMVA_vtx0);
  tree_->SetBranchAddress("diPhotMVA_vtxPair", &diPhotMVA_vtxPair, &b_diPhotMVA_vtxPair);
  tree_->SetBranchAddress("preselPairId", &preselPairId, &b_preselPairId);
  tree_->SetBranchAddress("tmva_dipho_MIT_dmom", &tmva_dipho_MIT_dmom, &b_tmva_dipho_MIT_dmom);
  tree_->SetBranchAddress("tmva_dipho_MIT_dmom_wrong_vtx", &tmva_dipho_MIT_dmom_wrong_vtx, &b_tmva_dipho_MIT_dmom_wrong_vtx);
  tree_->SetBranchAddress("tmva_dipho_MIT_vtxprob", &tmva_dipho_MIT_vtxprob, &b_tmva_dipho_MIT_vtxprob);
  tree_->SetBranchAddress("tmva_dipho_MIT_ptom1", &tmva_dipho_MIT_ptom1, &b_tmva_dipho_MIT_ptom1);
  tree_->SetBranchAddress("tmva_dipho_MIT_ptom2", &tmva_dipho_MIT_ptom2, &b_tmva_dipho_MIT_ptom2);
  tree_->SetBranchAddress("tmva_dipho_MIT_eta1", &tmva_dipho_MIT_eta1, &b_tmva_dipho_MIT_eta1);
  tree_->SetBranchAddress("tmva_dipho_MIT_eta2", &tmva_dipho_MIT_eta2, &b_tmva_dipho_MIT_eta2);
  tree_->SetBranchAddress("tmva_dipho_MIT_dphi", &tmva_dipho_MIT_dphi, &b_tmva_dipho_MIT_dphi);
  tree_->SetBranchAddress("tmva_dipho_MIT_ph1mva", &tmva_dipho_MIT_ph1mva, &b_tmva_dipho_MIT_ph1mva);
  tree_->SetBranchAddress("tmva_dipho_MIT_ph2mva", &tmva_dipho_MIT_ph2mva, &b_tmva_dipho_MIT_ph2mva);
  tree_->SetBranchAddress("ePUMet", &ePUMet, &b_ePUMet);
  tree_->SetBranchAddress("ePUMet2", &ePUMet2, &b_ePUMet2);
  tree_->SetBranchAddress("ePUMet3", &ePUMet3, &b_ePUMet3);
  tree_->SetBranchAddress("ePUMet4", &ePUMet4, &b_ePUMet4);
  tree_->SetBranchAddress("ePUMet5", &ePUMet5, &b_ePUMet5);
  tree_->SetBranchAddress("ecorrPUMet5", &ecorrPUMet5, &b_ecorrPUMet5);
  tree_->SetBranchAddress("phiPUMet", &phiPUMet, &b_phiPUMet);
  tree_->SetBranchAddress("phiPUMet2", &phiPUMet2, &b_phiPUMet2);
  tree_->SetBranchAddress("phiPUMet3", &phiPUMet3, &b_phiPUMet3);
  tree_->SetBranchAddress("phiPUMet4", &phiPUMet4, &b_phiPUMet4);
  tree_->SetBranchAddress("phiPUMet5", &phiPUMet5, &b_phiPUMet5);
  tree_->SetBranchAddress("phiCorrPUMet5", &phiCorrPUMet5, &b_phiCorrPUMet5);
  tree_->SetBranchAddress("phot1Metx", &phot1Metx, &b_phot1Metx);
  tree_->SetBranchAddress("phot2Metx", &phot2Metx, &b_phot2Metx);
  tree_->SetBranchAddress("leptonsMetx", &leptonsMetx, &b_leptonsMetx);
  tree_->SetBranchAddress("part_in_jetsMetx", &part_in_jetsMetx, &b_part_in_jetsMetx);
  tree_->SetBranchAddress("chg_vtx_unclMetx", &chg_vtx_unclMetx, &b_chg_vtx_unclMetx);
  tree_->SetBranchAddress("chg_novtx_unclMetx", &chg_novtx_unclMetx, &b_chg_novtx_unclMetx);
  tree_->SetBranchAddress("neutrals_unclMetx", &neutrals_unclMetx, &b_neutrals_unclMetx);
  tree_->SetBranchAddress("part_fail_puidMetx", &part_fail_puidMetx, &b_part_fail_puidMetx);
  tree_->SetBranchAddress("phot1Mety", &phot1Mety, &b_phot1Mety);
  tree_->SetBranchAddress("phot2Mety", &phot2Mety, &b_phot2Mety);
  tree_->SetBranchAddress("leptonsMety", &leptonsMety, &b_leptonsMety);
  tree_->SetBranchAddress("part_in_jetsMety", &part_in_jetsMety, &b_part_in_jetsMety);
  tree_->SetBranchAddress("chg_vtx_unclMety", &chg_vtx_unclMety, &b_chg_vtx_unclMety);
  tree_->SetBranchAddress("chg_novtx_unclMety", &chg_novtx_unclMety, &b_chg_novtx_unclMety);
  tree_->SetBranchAddress("neutrals_unclMety", &neutrals_unclMety, &b_neutrals_unclMety);
  tree_->SetBranchAddress("part_fail_puidMety", &part_fail_puidMety, &b_part_fail_puidMety);
  tree_->SetBranchAddress("scaling", &scaling, &b_scaling);
  tree_->SetBranchAddress("sMet", &sMet, &b_sMet);
  tree_->SetBranchAddress("eMet", &eMet, &b_eMet);
  tree_->SetBranchAddress("phiMet", &phiMet, &b_phiMet);
  tree_->SetBranchAddress("signifMet", &signifMet, &b_signifMet);
  tree_->SetBranchAddress("eSmearedMet", &eSmearedMet, &b_eSmearedMet);
  tree_->SetBranchAddress("phiSmearedMet", &phiSmearedMet, &b_phiSmearedMet);
  tree_->SetBranchAddress("eShiftedMet", &eShiftedMet, &b_eShiftedMet);
  tree_->SetBranchAddress("phiShiftedMet", &phiShiftedMet, &b_phiShiftedMet);
  tree_->SetBranchAddress("eShiftedScaledMet", &eShiftedScaledMet, &b_eShiftedScaledMet);
  tree_->SetBranchAddress("phiShiftedScaledMet", &phiShiftedScaledMet, &b_phiShiftedScaledMet);
  tree_->SetBranchAddress("eSmearedShiftedMet", &eSmearedShiftedMet, &b_eSmearedShiftedMet);
  tree_->SetBranchAddress("phiSmearedShiftedMet", &phiSmearedShiftedMet, &b_phiSmearedShiftedMet);
  tree_->SetBranchAddress("eShiftedScaledMetPUcorr", &eShiftedScaledMetPUcorr, &b_eShiftedScaledMetPUcorr);
  tree_->SetBranchAddress("phiShiftedScaledMetPUcorr", &phiShiftedScaledMetPUcorr, &b_phiShiftedScaledMetPUcorr);
  tree_->SetBranchAddress("eSmearedShiftedMePUcorrt", &eSmearedShiftedMePUcorrt, &b_eSmearedShiftedMetPUcorr);
  tree_->SetBranchAddress("phiSmearedShiftedMetPUcorr", &phiSmearedShiftedMetPUcorr, &b_phiSmearedShiftedMetPUcorr);
  tree_->SetBranchAddress("sCorrMet", &sCorrMet, &b_sCorrMet);
  tree_->SetBranchAddress("eCorrMet", &eCorrMet, &b_eCorrMet);
  tree_->SetBranchAddress("phiCorrMet", &phiCorrMet, &b_phiCorrMet);
  tree_->SetBranchAddress("signifCorrMet", &signifCorrMet, &b_signifCorrMet);
  tree_->SetBranchAddress("smuCorrMet", &smuCorrMet, &b_smuCorrMet);
  tree_->SetBranchAddress("emuCorrMet", &emuCorrMet, &b_emuCorrMet);
  tree_->SetBranchAddress("phimuCorrMet", &phimuCorrMet, &b_phimuCorrMet);
  tree_->SetBranchAddress("signifmuCorrMet", &signifmuCorrMet, &b_signifmuCorrMet);
  tree_->SetBranchAddress("sNoHFMet", &sNoHFMet, &b_sNoHFMet);
  tree_->SetBranchAddress("eNoHFMet", &eNoHFMet, &b_eNoHFMet);
  tree_->SetBranchAddress("phiNoHFMet", &phiNoHFMet, &b_phiNoHFMet);
  tree_->SetBranchAddress("signifNoHFMet", &signifNoHFMet, &b_signifNoHFMet);
  tree_->SetBranchAddress("stcMet", &stcMet, &b_stcMet);
  tree_->SetBranchAddress("etcMet", &etcMet, &b_etcMet);
  tree_->SetBranchAddress("phitcMet", &phitcMet, &b_phitcMet);
  tree_->SetBranchAddress("signiftcMet", &signiftcMet, &b_signiftcMet);
  tree_->SetBranchAddress("sglobalPfMet", &sglobalPfMet, &b_sglobalPfMet);
  tree_->SetBranchAddress("eglobalPfMet", &eglobalPfMet, &b_eglobalPfMet);
  tree_->SetBranchAddress("phiglobalPfMet", &phiglobalPfMet, &b_phiglobalPfMet);
  tree_->SetBranchAddress("signifglobalPfMet", &signifglobalPfMet, &b_signifglobalPfMet);
  tree_->SetBranchAddress("scentralPfMet", &scentralPfMet, &b_scentralPfMet);
  tree_->SetBranchAddress("ecentralPfMet", &ecentralPfMet, &b_ecentralPfMet);
  tree_->SetBranchAddress("phicentralPfMet", &phicentralPfMet, &b_phicentralPfMet);
  tree_->SetBranchAddress("signifcentralPfMet", &signifcentralPfMet, &b_signifcentralPfMet);
  tree_->SetBranchAddress("eassocPfMet", &eassocPfMet, &b_eassocPfMet);
  tree_->SetBranchAddress("phiassocPfMet", &phiassocPfMet, &b_phiassocPfMet);
  tree_->SetBranchAddress("signifassocPfMet", &signifassocPfMet, &b_signifassocPfMet);
  tree_->SetBranchAddress("eassocOtherVtxPfMet", &eassocOtherVtxPfMet, &b_eassocOtherVtxPfMet);
  tree_->SetBranchAddress("phiassocOtherVtxPfMet", &phiassocOtherVtxPfMet, &b_phiassocOtherVtxPfMet);
  tree_->SetBranchAddress("signifassocOtherVtxPfMet", &signifassocOtherVtxPfMet, &b_signifassocOtherVtxPfMet);
  tree_->SetBranchAddress("etrkPfMet", &etrkPfMet, &b_etrkPfMet);
  tree_->SetBranchAddress("phitrkPfMet", &phitrkPfMet, &b_phitrkPfMet);
  tree_->SetBranchAddress("signiftrkPfMet", &signiftrkPfMet, &b_signiftrkPfMet);
  tree_->SetBranchAddress("ecleanPfMet", &ecleanPfMet, &b_ecleanPfMet);
  tree_->SetBranchAddress("phicleanPfMet", &phicleanPfMet, &b_phicleanPfMet);
  tree_->SetBranchAddress("signifcleanPfMet", &signifcleanPfMet, &b_signifcleanPfMet);
  tree_->SetBranchAddress("ecleanedSaclayPfMet", &ecleanedSaclayPfMet, &b_ecleanedSaclayPfMet);
  tree_->SetBranchAddress("phicleanedSaclayPfMet", &phicleanedSaclayPfMet, &b_phicleanedSaclayPfMet);
  tree_->SetBranchAddress("signifcleanedSaclayPfMet", &signifcleanedSaclayPfMet, &b_signifcleanedSaclayPfMet);
  tree_->SetBranchAddress("eminTypeICleanSaclayPfMet", &eminTypeICleanSaclayPfMet, &b_eminTypeICleanSaclayPfMet);
  tree_->SetBranchAddress("phiminTypeICleanSaclayPfMet", &phiminTypeICleanSaclayPfMet, &b_phiminTypeICleanSaclayPfMet);
  tree_->SetBranchAddress("signifminTypeICleanSaclayPfMet", &signifminTypeICleanSaclayPfMet, &b_signifminTypeICleanSaclayPfMet);
  tree_->SetBranchAddress("globalPfSums", &globalPfSums, &b_globalPfSums);
  tree_->SetBranchAddress("spfMet", &spfMet, &b_spfMet);
  tree_->SetBranchAddress("epfMet", &epfMet, &b_epfMet);
  tree_->SetBranchAddress("phipfMet", &phipfMet, &b_phipfMet);
  tree_->SetBranchAddress("signifpfMet", &signifpfMet, &b_signifpfMet);
  tree_->SetBranchAddress("spfMetType1", &spfMetType1, &b_spfMetType1);
  tree_->SetBranchAddress("epfMetType1", &epfMetType1, &b_epfMetType1);
  tree_->SetBranchAddress("phipfMetType1", &phipfMetType1, &b_phipfMetType1);
  tree_->SetBranchAddress("signifpfMetType1", &signifpfMetType1, &b_signifpfMetType1);
  tree_->SetBranchAddress("sMetGen", &sMetGen, &b_sMetGen);
  tree_->SetBranchAddress("eMetGen", &eMetGen, &b_eMetGen);
  tree_->SetBranchAddress("phiMetGen", &phiMetGen, &b_phiMetGen);
  tree_->SetBranchAddress("signifMetGen", &signifMetGen, &b_signifMetGen);
  tree_->SetBranchAddress("sMetGen2", &sMetGen2, &b_sMetGen2);
  tree_->SetBranchAddress("eMetGen2", &eMetGen2, &b_eMetGen2);
  tree_->SetBranchAddress("phiMetGen2", &phiMetGen2, &b_phiMetGen2);
  tree_->SetBranchAddress("npu", &npu, &b_npu);
  tree_->SetBranchAddress("NtotEvents", &NtotEvents, &b_NtotEvents);
  tree_->SetBranchAddress("xsection", &xsection, &b_xsection);
  tree_->SetBranchAddress("EquivLumi", &EquivLumi, &b_EquivLumi);
  tree_->SetBranchAddress("SampleID", &SampleID, &b_SampleID);
  tree_->SetBranchAddress("pu_weight", &pu_weight, &b_pu_weight);
  tree_->SetBranchAddress("pt_weight", &pt_weight, &b_pt_weight);
  tree_->SetBranchAddress("gen_custom_processId", &gen_custom_processId, &b_gen_custom_processId);
  tree_->SetBranchAddress("gen_pt_gamma1", &gen_pt_gamma1, &b_gen_pt_gamma1);
  tree_->SetBranchAddress("gen_pt_gamma2", &gen_pt_gamma2, &b_gen_pt_gamma2);
  tree_->SetBranchAddress("gen_eta_gamma1", &gen_eta_gamma1, &b_gen_eta_gamma1);
  tree_->SetBranchAddress("gen_eta_gamma2", &gen_eta_gamma2, &b_gen_eta_gamma2);
  tree_->SetBranchAddress("gen_phi_gamma1", &gen_phi_gamma1, &b_gen_phi_gamma1);
  tree_->SetBranchAddress("gen_phi_gamma2", &gen_phi_gamma2, &b_gen_phi_gamma2);
  tree_->SetBranchAddress("gen_mass_diphoton", &gen_mass_diphoton, &b_gen_mass_diphoton);
  tree_->SetBranchAddress("gen_pt_diphoton", &gen_pt_diphoton, &b_gen_pt_diphoton);
  tree_->SetBranchAddress("gen_eta_diphoton", &gen_eta_diphoton, &b_gen_eta_diphoton);
  tree_->SetBranchAddress("gen_phi_diphoton", &gen_phi_diphoton, &b_gen_phi_diphoton);
  tree_->SetBranchAddress("gen_mass_dijet", &gen_mass_dijet, &b_gen_mass_dijet);
  tree_->SetBranchAddress("gen_pt_dijet", &gen_pt_dijet, &b_gen_pt_dijet);
  tree_->SetBranchAddress("gen_eta_dijet", &gen_eta_dijet, &b_gen_eta_dijet);
  tree_->SetBranchAddress("gen_phi_dijet", &gen_phi_dijet, &b_gen_phi_dijet);
  tree_->SetBranchAddress("gen_zeppenfeld", &gen_zeppenfeld, &b_gen_zeppenfeld);
  tree_->SetBranchAddress("gen_pt_lep1", &gen_pt_lep1, &b_gen_pt_lep1);
  tree_->SetBranchAddress("gen_pt_lep2", &gen_pt_lep2, &b_gen_pt_lep2);
  tree_->SetBranchAddress("gen_eta_lep1", &gen_eta_lep1, &b_gen_eta_lep1);
  tree_->SetBranchAddress("gen_eta_lep2", &gen_eta_lep2, &b_gen_eta_lep2);
  tree_->SetBranchAddress("gen_phi_lep1", &gen_phi_lep1, &b_gen_phi_lep1);
  tree_->SetBranchAddress("gen_phi_lep2", &gen_phi_lep2, &b_gen_phi_lep2);
  tree_->SetBranchAddress("gen_pid_lep1", &gen_pid_lep1, &b_gen_pid_lep1);
  tree_->SetBranchAddress("gen_pid_lep2", &gen_pid_lep2, &b_gen_pid_lep2);
  tree_->SetBranchAddress("ptele1", &ptele1, &b_ptele1);
  tree_->SetBranchAddress("ptele2", &ptele2, &b_ptele2);
  tree_->SetBranchAddress("etaele1", &etaele1, &b_etaele1);
  tree_->SetBranchAddress("etaele2", &etaele2, &b_etaele2);
  tree_->SetBranchAddress("phiele1", &phiele1, &b_phiele1);
  tree_->SetBranchAddress("phiele2", &phiele2, &b_phiele2);
  tree_->SetBranchAddress("eneele1", &eneele1, &b_eneele1);
  tree_->SetBranchAddress("eneele2", &eneele2, &b_eneele2);
  tree_->SetBranchAddress("sIeIeele1", &sIeIeele1, &b_sIeIeele1);
  tree_->SetBranchAddress("sIeIeele2", &sIeIeele2, &b_sIeIeele2);
  tree_->SetBranchAddress("dphiele1", &dphiele1, &b_dphiele1);
  tree_->SetBranchAddress("dphiele2", &dphiele2, &b_dphiele2);
  tree_->SetBranchAddress("detaele1", &detaele1, &b_detaele1);
  tree_->SetBranchAddress("detaele2", &detaele2, &b_detaele2);
  tree_->SetBranchAddress("hoeele1", &hoeele1, &b_hoeele1);
  tree_->SetBranchAddress("hoeele2", &hoeele2, &b_hoeele2);
  tree_->SetBranchAddress("mhitsele1", &mhitsele1, &b_mhitsele1);
  tree_->SetBranchAddress("mhitsele2", &mhitsele2, &b_mhitsele2);
  tree_->SetBranchAddress("d0ele1", &d0ele1, &b_d0ele1);
  tree_->SetBranchAddress("d0ele2", &d0ele2, &b_d0ele2);
  tree_->SetBranchAddress("dzele1", &dzele1, &b_dzele1);
  tree_->SetBranchAddress("dzele2", &dzele2, &b_dzele2);
  tree_->SetBranchAddress("invMassele1g1", &invMassele1g1, &b_invMassele1g1);
  tree_->SetBranchAddress("invMassele1g2", &invMassele1g2, &b_invMassele1g2);
  tree_->SetBranchAddress("invMassele2g1", &invMassele2g1, &b_invMassele2g1);
  tree_->SetBranchAddress("invMassele2g2", &invMassele2g2, &b_invMassele2g2);
  tree_->SetBranchAddress("oEmoPele1", &oEmoPele1, &b_oEmoPele1);
  tree_->SetBranchAddress("oEmoPele2", &oEmoPele2, &b_oEmoPele2);
  tree_->SetBranchAddress("mvanotrigele1", &mvanotrigele1, &b_mvanotrigele1);
  tree_->SetBranchAddress("mvanotrigele2", &mvanotrigele2, &b_mvanotrigele2);
  tree_->SetBranchAddress("mvatrigele1", &mvatrigele1, &b_mvatrigele1);
  tree_->SetBranchAddress("mvatrigele2", &mvatrigele2, &b_mvatrigele2);
  tree_->SetBranchAddress("matchconvele1", &matchconvele1, &b_matchconvele1);
  tree_->SetBranchAddress("matchconvele2", &matchconvele2, &b_matchconvele2);
  tree_->SetBranchAddress("chHadIso03ele1", &chHadIso03ele1, &b_chHadIso03ele1);
  tree_->SetBranchAddress("chHadIso03ele2", &chHadIso03ele2, &b_chHadIso03ele2);
  tree_->SetBranchAddress("nHadIso03ele1", &nHadIso03ele1, &b_nHadIso03ele1);
  tree_->SetBranchAddress("nHadIso03ele2", &nHadIso03ele2, &b_nHadIso03ele2);
  tree_->SetBranchAddress("photIso03ele1", &photIso03ele1, &b_photIso03ele1);
  tree_->SetBranchAddress("photIso03ele2", &photIso03ele2, &b_photIso03ele2);
  tree_->SetBranchAddress("pteleloose1", &pteleloose1, &b_pteleloose1);
  tree_->SetBranchAddress("pteleloose2", &pteleloose2, &b_pteleloose2);
  tree_->SetBranchAddress("etaeleloose1", &etaeleloose1, &b_etaeleloose1);
  tree_->SetBranchAddress("etaeleloose2", &etaeleloose2, &b_etaeleloose2);
  tree_->SetBranchAddress("phieleloose1", &phieleloose1, &b_phieleloose1);
  tree_->SetBranchAddress("phieleloose2", &phieleloose2, &b_phieleloose2);
  tree_->SetBranchAddress("eneeleloose1", &eneeleloose1, &b_eneeleloose1);
  tree_->SetBranchAddress("eneeleloose2", &eneeleloose2, &b_eneeleloose2);
  tree_->SetBranchAddress("sIeIeeleloose1", &sIeIeeleloose1, &b_sIeIeeleloose1);
  tree_->SetBranchAddress("sIeIeeleloose2", &sIeIeeleloose2, &b_sIeIeeleloose2);
  tree_->SetBranchAddress("dphieleloose1", &dphieleloose1, &b_dphieleloose1);
  tree_->SetBranchAddress("dphieleloose2", &dphieleloose2, &b_dphieleloose2);
  tree_->SetBranchAddress("detaeleloose1", &detaeleloose1, &b_detaeleloose1);
  tree_->SetBranchAddress("detaeleloose2", &detaeleloose2, &b_detaeleloose2);
  tree_->SetBranchAddress("hoeeleloose1", &hoeeleloose1, &b_hoeeleloose1);
  tree_->SetBranchAddress("hoeeleloose2", &hoeeleloose2, &b_hoeeleloose2);
  tree_->SetBranchAddress("mhitseleloose1", &mhitseleloose1, &b_mhitseleloose1);
  tree_->SetBranchAddress("mhitseleloose2", &mhitseleloose2, &b_mhitseleloose2);
  tree_->SetBranchAddress("d0eleloose1", &d0eleloose1, &b_d0eleloose1);
  tree_->SetBranchAddress("d0eleloose2", &d0eleloose2, &b_d0eleloose2);
  tree_->SetBranchAddress("dzeleloose1", &dzeleloose1, &b_dzeleloose1);
  tree_->SetBranchAddress("dzeleloose2", &dzeleloose2, &b_dzeleloose2);
  tree_->SetBranchAddress("invMasseleloose1g1", &invMasseleloose1g1, &b_invMasseleloose1g1);
  tree_->SetBranchAddress("invMasseleloose1g2", &invMasseleloose1g2, &b_invMasseleloose1g2);
  tree_->SetBranchAddress("invMasseleloose2g1", &invMasseleloose2g1, &b_invMasseleloose2g1);
  tree_->SetBranchAddress("invMasseleloose2g2", &invMasseleloose2g2, &b_invMasseleloose2g2);
  tree_->SetBranchAddress("oEmoPeleloose1", &oEmoPeleloose1, &b_oEmoPeleloose1);
  tree_->SetBranchAddress("oEmoPeleloose2", &oEmoPeleloose2, &b_oEmoPeleloose2);
  tree_->SetBranchAddress("mvanotrigeleloose1", &mvanotrigeleloose1, &b_mvanotrigeleloose1);
  tree_->SetBranchAddress("mvanotrigeleloose2", &mvanotrigeleloose2, &b_mvanotrigeleloose2);
  tree_->SetBranchAddress("mvatrigeleloose1", &mvatrigeleloose1, &b_mvatrigeleloose1);
  tree_->SetBranchAddress("mvatrigeleloose2", &mvatrigeleloose2, &b_mvatrigeleloose2);
  tree_->SetBranchAddress("matchconveleloose1", &matchconveleloose1, &b_matchconveleloose1);
  tree_->SetBranchAddress("matchconveleloose2", &matchconveleloose2, &b_matchconveleloose2);
  tree_->SetBranchAddress("chHadIso03eleloose1", &chHadIso03eleloose1, &b_chHadIso03eleloose1);
  tree_->SetBranchAddress("chHadIso03eleloose2", &chHadIso03eleloose2, &b_chHadIso03eleloose2);
  tree_->SetBranchAddress("nHadIso03eleloose1", &nHadIso03eleloose1, &b_nHadIso03eleloose1);
  tree_->SetBranchAddress("nHadIso03eleloose2", &nHadIso03eleloose2, &b_nHadIso03eleloose2);
  tree_->SetBranchAddress("photIso03eleloose1", &photIso03eleloose1, &b_photIso03eleloose1);
  tree_->SetBranchAddress("photIso03eleloose2", &photIso03eleloose2, &b_photIso03eleloose2);
  tree_->SetBranchAddress("ptelenontr801", &ptelenontr801, &b_ptelenontr801);
  tree_->SetBranchAddress("ptelenontr802", &ptelenontr802, &b_ptelenontr802);
  tree_->SetBranchAddress("etaelenontr801", &etaelenontr801, &b_etaelenontr801);
  tree_->SetBranchAddress("etaelenontr802", &etaelenontr802, &b_etaelenontr802);
  tree_->SetBranchAddress("phielenontr801", &phielenontr801, &b_phielenontr801);
  tree_->SetBranchAddress("phielenontr802", &phielenontr802, &b_phielenontr802);
  tree_->SetBranchAddress("eneelenontr801", &eneelenontr801, &b_eneelenontr801);
  tree_->SetBranchAddress("eneelenontr802", &eneelenontr802, &b_eneelenontr802);
  tree_->SetBranchAddress("ptelenontr901", &ptelenontr901, &b_ptelenontr901);
  tree_->SetBranchAddress("ptelenontr902", &ptelenontr902, &b_ptelenontr902);
  tree_->SetBranchAddress("etaelenontr901", &etaelenontr901, &b_etaelenontr901);
  tree_->SetBranchAddress("etaelenontr902", &etaelenontr902, &b_etaelenontr902);
  tree_->SetBranchAddress("phielenontr901", &phielenontr901, &b_phielenontr901);
  tree_->SetBranchAddress("phielenontr902", &phielenontr902, &b_phielenontr902);
  tree_->SetBranchAddress("eneelenontr901", &eneelenontr901, &b_eneelenontr901);
  tree_->SetBranchAddress("eneelenontr902", &eneelenontr902, &b_eneelenontr902);
  tree_->SetBranchAddress("ptmu1", &ptmu1, &b_ptmu1);
  tree_->SetBranchAddress("ptmu2", &ptmu2, &b_ptmu2);
  tree_->SetBranchAddress("etamu1", &etamu1, &b_etamu1);
  tree_->SetBranchAddress("etamu2", &etamu2, &b_etamu2);
  tree_->SetBranchAddress("phimu1", &phimu1, &b_phimu1);
  tree_->SetBranchAddress("phimu2", &phimu2, &b_phimu2);
  tree_->SetBranchAddress("enemu1", &enemu1, &b_enemu1);
  tree_->SetBranchAddress("enemu2", &enemu2, &b_enemu2);
  tree_->SetBranchAddress("pixhitsmu1", &pixhitsmu1, &b_pixhitsmu1);
  tree_->SetBranchAddress("pixhitsmu2", &pixhitsmu2, &b_pixhitsmu2);
  tree_->SetBranchAddress("trkhitsmu1", &trkhitsmu1, &b_trkhitsmu1);
  tree_->SetBranchAddress("trkhitsmu2", &trkhitsmu2, &b_trkhitsmu2);
  tree_->SetBranchAddress("hitsmu1", &hitsmu1, &b_hitsmu1);
  tree_->SetBranchAddress("hitsmu2", &hitsmu2, &b_hitsmu2);
  tree_->SetBranchAddress("chi2mu1", &chi2mu1, &b_chi2mu1);
  tree_->SetBranchAddress("chi2mu2", &chi2mu2, &b_chi2mu2);
  tree_->SetBranchAddress("matchmu1", &matchmu1, &b_matchmu1);
  tree_->SetBranchAddress("matchmu2", &matchmu2, &b_matchmu2);
  tree_->SetBranchAddress("d0mu1", &d0mu1, &b_d0mu1);
  tree_->SetBranchAddress("d0mu2", &d0mu2, &b_d0mu2);
  tree_->SetBranchAddress("dzmu1", &dzmu1, &b_dzmu1);
  tree_->SetBranchAddress("dzmu2", &dzmu2, &b_dzmu2);
  tree_->SetBranchAddress("chHadmu1", &chHadmu1, &b_chHadmu1);
  tree_->SetBranchAddress("chHadmu2", &chHadmu2, &b_chHadmu2);
  tree_->SetBranchAddress("nHadmu1", &nHadmu1, &b_nHadmu1);
  tree_->SetBranchAddress("nHadmu2", &nHadmu2, &b_nHadmu2);
  tree_->SetBranchAddress("photmu1", &photmu1, &b_photmu1);
  tree_->SetBranchAddress("photmu2", &photmu2, &b_photmu2);
  tree_->SetBranchAddress("puptmu1", &puptmu1, &b_puptmu1);
  tree_->SetBranchAddress("puptmu2", &puptmu2, &b_puptmu2);
  tree_->SetBranchAddress("ptmuloose1", &ptmuloose1, &b_ptmuloose1);
  tree_->SetBranchAddress("ptmuloose2", &ptmuloose2, &b_ptmuloose2);
  tree_->SetBranchAddress("etamuloose1", &etamuloose1, &b_etamuloose1);
  tree_->SetBranchAddress("etamuloose2", &etamuloose2, &b_etamuloose2);
  tree_->SetBranchAddress("phimuloose1", &phimuloose1, &b_phimuloose1);
  tree_->SetBranchAddress("phimuloose2", &phimuloose2, &b_phimuloose2);
  tree_->SetBranchAddress("enemuloose1", &enemuloose1, &b_enemuloose1);
  tree_->SetBranchAddress("enemuloose2", &enemuloose2, &b_enemuloose2);
  tree_->SetBranchAddress("pixhitsmuloose1", &pixhitsmuloose1, &b_pixhitsmuloose1);
  tree_->SetBranchAddress("pixhitsmuloose2", &pixhitsmuloose2, &b_pixhitsmuloose2);
  tree_->SetBranchAddress("trkhitsmuloose1", &trkhitsmuloose1, &b_trkhitsmuloose1);
  tree_->SetBranchAddress("trkhitsmuloose2", &trkhitsmuloose2, &b_trkhitsmuloose2);
  tree_->SetBranchAddress("hitsmuloose1", &hitsmuloose1, &b_hitsmuloose1);
  tree_->SetBranchAddress("hitsmuloose2", &hitsmuloose2, &b_hitsmuloose2);
  tree_->SetBranchAddress("chi2muloose1", &chi2muloose1, &b_chi2muloose1);
  tree_->SetBranchAddress("chi2muloose2", &chi2muloose2, &b_chi2muloose2);
  tree_->SetBranchAddress("matchmuloose1", &matchmuloose1, &b_matchmuloose1);
  tree_->SetBranchAddress("matchmuloose2", &matchmuloose2, &b_matchmuloose2);
  tree_->SetBranchAddress("d0muloose1", &d0muloose1, &b_d0muloose1);
  tree_->SetBranchAddress("d0muloose2", &d0muloose2, &b_d0muloose2);
  tree_->SetBranchAddress("dzmuloose1", &dzmuloose1, &b_dzmuloose1);
  tree_->SetBranchAddress("dzmuloose2", &dzmuloose2, &b_dzmuloose2);
  tree_->SetBranchAddress("ptmuvloose1", &ptmuvloose1, &b_ptmuvloose1);
  tree_->SetBranchAddress("ptmuvloose2", &ptmuvloose2, &b_ptmuvloose2);
  tree_->SetBranchAddress("etamuvloose1", &etamuvloose1, &b_etamuvloose1);
  tree_->SetBranchAddress("etamuvloose2", &etamuvloose2, &b_etamuvloose2);
  tree_->SetBranchAddress("phimuvloose1", &phimuvloose1, &b_phimuvloose1);
  tree_->SetBranchAddress("phimuvloose2", &phimuvloose2, &b_phimuvloose2);
  tree_->SetBranchAddress("enemuvloose1", &enemuvloose1, &b_enemuvloose1);
  tree_->SetBranchAddress("enemuvloose2", &enemuvloose2, &b_enemuvloose2);
  tree_->SetBranchAddress("pixhitsmuvloose1", &pixhitsmuvloose1, &b_pixhitsmuvloose1);
  tree_->SetBranchAddress("pixhitsmuvloose2", &pixhitsmuvloose2, &b_pixhitsmuvloose2);
  tree_->SetBranchAddress("trkhitsmuvloose1", &trkhitsmuvloose1, &b_trkhitsmuvloose1);
  tree_->SetBranchAddress("trkhitsmuvloose2", &trkhitsmuvloose2, &b_trkhitsmuvloose2);
  tree_->SetBranchAddress("hitsmuvloose1", &hitsmuvloose1, &b_hitsmuvloose1);
  tree_->SetBranchAddress("hitsmuvloose2", &hitsmuvloose2, &b_hitsmuvloose2);
  tree_->SetBranchAddress("chi2muvloose1", &chi2muvloose1, &b_chi2muvloose1);
  tree_->SetBranchAddress("chi2muvloose2", &chi2muvloose2, &b_chi2muvloose2);
  tree_->SetBranchAddress("matchmuvloose1", &matchmuvloose1, &b_matchmuvloose1);
  tree_->SetBranchAddress("matchmuvloose2", &matchmuvloose2, &b_matchmuvloose2);
  tree_->SetBranchAddress("d0muvloose1", &d0muvloose1, &b_d0muvloose1);
  tree_->SetBranchAddress("d0muvloose2", &d0muvloose2, &b_d0muvloose2);
  tree_->SetBranchAddress("dzmuvloose1", &dzmuvloose1, &b_dzmuvloose1);
  tree_->SetBranchAddress("dzmuvloose2", &dzmuvloose2, &b_dzmuvloose2);
  tree_->SetBranchAddress("hasPassedSinglePhot", &hasPassedSinglePhot, &b_hasPassedSinglePhot);
  tree_->SetBranchAddress("hasPassedDoublePhot", &hasPassedDoublePhot, &b_hasPassedDoublePhot);
  tree_->SetBranchAddress("chHadmuloose1", &chHadmuloose1, &b_chHadmuloose1);
  tree_->SetBranchAddress("chHadmuloose2", &chHadmuloose2, &b_chHadmuloose2);
  tree_->SetBranchAddress("nHadmuloose1", &nHadmuloose1, &b_nHadmuloose1);
  tree_->SetBranchAddress("nHadmuloose2", &nHadmuloose2, &b_nHadmuloose2);
  tree_->SetBranchAddress("photmuloose1", &photmuloose1, &b_photmuloose1);
  tree_->SetBranchAddress("photmuloose2", &photmuloose2, &b_photmuloose2);
  tree_->SetBranchAddress("puptmuloose1", &puptmuloose1, &b_puptmuloose1);
  tree_->SetBranchAddress("puptmuloose2", &puptmuloose2, &b_puptmuloose2);
  tree_->SetBranchAddress("chHadmuvloose1", &chHadmuvloose1, &b_chHadmuvloose1);
  tree_->SetBranchAddress("chHadmuvloose2", &chHadmuvloose2, &b_chHadmuvloose2);
  tree_->SetBranchAddress("nHadmuvloose1", &nHadmuvloose1, &b_nHadmuvloose1);
  tree_->SetBranchAddress("nHadmuvloose2", &nHadmuvloose2, &b_nHadmuvloose2);
  tree_->SetBranchAddress("photmuvloose1", &photmuvloose1, &b_photmuvloose1);
  tree_->SetBranchAddress("photmuvloose2", &photmuvloose2, &b_photmuvloose2);
  tree_->SetBranchAddress("puptmuvloose1", &puptmuvloose1, &b_puptmuvloose1);
  tree_->SetBranchAddress("puptmuvloose2", &puptmuvloose2, &b_puptmuvloose2);
}

void fillPlot2012_radion::setSelectionType( const std::string& selectionType ) {

  selectionType_ = selectionType;
  
  // default values                                                                                                     
  dopureeventWeight_ = true;
  cicselection = 4;

  ptphot1cut = 40.;
  ptphot2cut = 25.;

  ptjetacccut  = 25.;
  etajetacccut = 2.5;  

  if( selectionType=="default" ) {

    ptdiphotcut_1bj = 0;
    ptdiphotcut_2bj = 0;

    etadiphotcut_1bj = 2.1;
    etadiphotcut_2bj = 2.0;  
    
    ptphot1cut_1bj = 52.;
    ptphot1cut_2bj = 57.;

    ptjet1cut_1bj = 41.;
    ptjet1cut_2bj = 48.;

    costhetascut_1bj = 0.9;
    costhetascut_2bj = 1.0;

  } else {
    
    std::cout << std::endl << std::endl << "Selection '" << selectionType << "' currently not implemented. Exiting." <<	std::endl;
    exit(12345);
  }
}

std::pair<int,int> fillPlot2012_radion::myTwoHighestPtJet(std::vector<int> jets) {

  float secondJpt = -999.;
  float firstJpt  = -998.;
  int secondJ     = -999;
  int firstJ      = -999;
  
  for (int ij=0; ij<int(jets.size());ij++) {
    int jetIndex = jets[ij]; 

    if (ptcorrjet[jetIndex]>firstJpt) {
      secondJpt = firstJpt;
      firstJpt  = ptcorrjet[jetIndex];
      secondJ   = firstJ;
      firstJ    = jetIndex;
    } else if (ptcorrjet[ij]>secondJpt) {
      secondJpt = ptcorrjet[jetIndex];
      secondJ   = jetIndex;
    }
  }

  return make_pair(firstJ,secondJ);
}

int fillPlot2012_radion::myHighestPtJet(std::vector<int> jets) {

  float firstJpt = -998.;
  int firstJ     = -999;
  
  for (int ij=0; ij<int(jets.size());ij++) {
    int jetIndex = jets[ij]; 
    if (ptcorrjet[jetIndex]>firstJpt) {
      firstJpt  = ptcorrjet[jetIndex];
      firstJ    = jetIndex;
    } 
  }

  return firstJ;
}

bool fillPlot2012_radion::passCutBasedJetId(int jet) {

  bool isGood = true;

  if ( fabs(etajet[jet]) < 2.5 ) {
    if ( betastarjet[jet] > 0.2 * log( nvtx - 0.64 ) ) isGood = false;
    if ( rmsjet[jet] > 0.06)                           isGood = false;
  } else if (fabs(etajet[jet]) < 3.){
    if ( rmsjet[jet] > 0.05)  isGood =false;
  } else {
    if ( rmsjet[jet] > 0.055) isGood =false;
  }
  
  return isGood;
}
