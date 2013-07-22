// C++ includes
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <math.h>
#include <TMath.h>
#include <TF1.h>
#include <TH1.h>
#include <TFile.h>
#include <TTree.h>

#include "cl95cms.C"

using namespace std;

double UL(double expsig, double expbkg, double xsec){

  float effS = expsig / (19600.*xsec);
  float ul = CLA( 19600., 0., effS, 0., expbkg, 0. );                                                                                            
  return ul/xsec;
 }

pair<float,float> runOptimization(TTree *theTree, TString cut0, TString cut1) {

  TH1F *mjj0 = new TH1F("mjj0", "", 100, 0.0, 10000.0);
  TH1F *mjj1 = new TH1F("mjj1", "", 100, 0.0, 10000.0);
  
  theTree->Project("mjj0","mjj", cut0);
  theTree->Project("mjj1","mjj", cut1);

  float afterCut   = mjj1->Integral();
  float efficiency = mjj1->Integral()/mjj0->Integral();

  cout << "cut before: " << cut0 << endl;
  cout << "cut after:  " << cut1 << endl;
  cout << "eff = " << efficiency << ", mjj->Integral() = " << mjj1->Integral() << endl;

  return std::make_pair<float,float>(afterCut, efficiency);  
}

int main(int argc, char* argv[]) {

  // chiara
  bool isCS = false;
  int theMass = 300;

  // loading files
  TFile *fileSig, *fileBkg;
  TTree *treeSig, *treeBkg;
  if (!isCS) fileBkg = TFile::Open("finalizedTrees_Radion_presel/Radion_Data2012_default_CSV.root");
  // if (!isCS) fileBkg = TFile::Open("finalizedTrees_Radion_presel/Radion_HToGG_M-125_8TeV-pythia6_default_CSV.root");
  // if (!isCS) fileBkg = TFile::Open("finalizedTrees_Radion_presel/Radion_GluGluToHToGG_M-125_8TeV_default_CSV.root"); 
  if ( isCS) fileBkg = TFile::Open("controlSampleStudy/treesFromCS_presel_withWeights_fitToMggbb/csWithWeightFromGammas.root");
  fileSig = TFile::Open("finalizedTrees_Radion_presel/Radion_Radion_M-300_madgraph_default_CSV.root");
  // fileSig = TFile::Open("finalizedTrees_Radion_presel/Radion_Radion_M-500_madgraph_default_CSV.root");                                  
  // fileSig = TFile::Open("finalizedTrees_Radion_presel/Radion_Radion_M-700_madgraph_default_CSV.root");                                  
  // fileSig = TFile::Open("finalizedTrees_Radion_presel/Radion_Radion_M-1000_madgraph_default_CSV.root");

  if( fileSig && fileBkg) {
    fileBkg->cd();
    if (!isCS) treeBkg = (TTree*)fileBkg->Get("myTrees");
    if (isCS)  treeBkg = (TTree*)fileBkg->Get("myTrees_withWeight");
    fileSig->cd();
    treeSig = (TTree*)fileSig->Get("myTrees");
  } else {
    cout << "File " << fileSig << " or " << fileBkg << " not existing !" << endl;
    return 0;
  }
  if(!treeSig || !treeBkg) {
    cout << "AnaTree not existing in background and/or signal files!" << endl;
    return 0;
  }

  cout << "All files OK!" << endl;
  cout << "signal: "     << treeSig->GetEntries() << " entries" << endl;
  cout << "background: " << treeBkg->GetEntries() << " entries" << endl;
  cout << endl;


  // to run cuts optimization 
  /*
  // per costheta*
  TH1F *signEffMap_1btag = new TH1F("signEffMap_1btag", "",14, 0.3, 1.);
  TH1F *backEffMap_1btag = new TH1F("backEffMap_1btag", "",14, 0.3, 1.);
  TH1F *ulMap_1btag      = new TH1F("ulMap_1btag",      "",14, 0.3, 1.);
  TH1F *signEffMap_2btag = new TH1F("signEffMap_2btag", "",14, 0.3, 1.);
  TH1F *backEffMap_2btag = new TH1F("backEffMap_2btag", "",14, 0.3, 1.);
  TH1F *ulMap_2btag      = new TH1F("ulMap_2btag",      "",14, 0.3, 1.);
  */
  /* 
  // per pT(gg) e pT(jj)    
  TH1F *signEffMap_1btag = new TH1F("signEffMap_1btag", "",14, 20., 90.);
  TH1F *backEffMap_1btag = new TH1F("backEffMap_1btag", "",14, 20., 90.);
  TH1F *ulMap_1btag      = new TH1F("ulMap_1btag",      "",14, 20., 90.);
  TH1F *signEffMap_2btag = new TH1F("signEffMap_2btag", "",14, 20., 90.);
  TH1F *backEffMap_2btag = new TH1F("backEffMap_2btag", "",14, 20., 90.);
  TH1F *ulMap_2btag      = new TH1F("ulMap_2btag",      "",14, 20., 90.);
  */
  // per minDeltaR
  TH1F *signEffMap_1btag = new TH1F("signEffMap_1btag", "",14, 0.2, 2.4);
  TH1F *backEffMap_1btag = new TH1F("backEffMap_1btag", "",14, 0.2, 2.4);
  TH1F *ulMap_1btag      = new TH1F("ulMap_1btag",      "",14, 0.2, 2.4);
  TH1F *signEffMap_2btag = new TH1F("signEffMap_2btag", "",14, 0.2, 2.4);
  TH1F *backEffMap_2btag = new TH1F("backEffMap_2btag", "",14, 0.2, 2.4);
  TH1F *ulMap_2btag      = new TH1F("ulMap_2btag",      "",14, 0.2, 2.4);  

  ulMap_1btag ->GetXaxis() -> SetTitle("inf cut [GeV]");
  ulMap_2btag ->GetXaxis() -> SetTitle("inf cut [GeV]");
  
  for (int ii=0; ii<15; ii++ ){

    float thePtDigammaCut       = 0.;     // 0.;
    float thePtDijetCut         = 0.;     // 0.;
    float theAbsCosThetaStarCut = 0.9;    // 0.9;
    float theDeltaRbgCut        = 0.;     // 5.;
    
    // float thisDijetPtCut = 20. + ii*5;
    float thisDijetPtCut = thePtDijetCut;
    stringstream ssJetPtCut (stringstream::in | stringstream::out);
    ssJetPtCut << thisDijetPtCut;
    
    // float thisDigammaPtCut = 20. + ii*5;
    float thisDigammaPtCut = thePtDigammaCut;
    stringstream ssGammaPtCut (stringstream::in | stringstream::out);
    ssGammaPtCut << thisDigammaPtCut;

    // float thisAbsCosThetaStarCut = 0.3 + ii*0.05;
    float thisAbsCosThetaStarCut = theAbsCosThetaStarCut;
    stringstream ssAbsCosThetaStarCut (stringstream::in | stringstream::out);
    ssAbsCosThetaStarCut << thisAbsCosThetaStarCut;

    float thisDeltaRbgCut = 0.2 + ii*0.16;
    // float thisDeltaRbgCut = theDeltaRbgCut;
    stringstream ssDeltaRbgCut (stringstream::in | stringstream::out);
    ssDeltaRbgCut << thisDeltaRbgCut;


    // cuts
    vector<TString> cutOpt0,  cutOpt1;
    vector<TString> cutOpt0b, cutOpt1b;
    
    // cutOpt0.push_back( TString("1.*weight*(btagCategory==1 && mjj>90. && mjj<170. && mggjj>270. && mggjj<330. && ptjj>") + ssJetPtCut.str() + TString(" && minDeltaR_gb>") + ssDeltaRbgCut.str() + TString(" && absCosThetaStar<") + ssAbsCosThetaStarCut.str() + TString(")"));
    // cutOpt0.push_back( TString("1.*weight*(btagCategory==1 && mjj>90. && mjj<170. && mggjj>270. && mggjj<330. && ptgg>") + ssGammaPtCut.str() + TString(" && minDeltaR_gb>") + ssDeltaRbgCut.str() + TString(" && absCosThetaStar<") + ssAbsCosThetaStarCut.str() + TString(")") );
    // cutOpt0.push_back( TString("1.*weight*(btagCategory==1 && mjj>90. && mjj<170. && mggjj>270. && mggjj<330. && ptgg>") + ssGammaPtCut.str() + TString(" && minDeltaR_gb>") + ssDeltaRbgCut.str() + TString(" && ptjj>") + ssJetPtCut.str() + TString(")") );
    cutOpt0.push_back( TString("1.*weight*(btagCategory==1 && mjj>90. && mjj<170. && mggjj>270. && mggjj<330. && ptgg>") + ssGammaPtCut.str() + TString(" && ptjj>") + ssJetPtCut.str() + TString(" && absCosThetaStar<") + ssAbsCosThetaStarCut.str() + TString(")") );
    

    // cutOpt0.push_back( TString("1.*weight*(btagCategory==2 && mjj>95. && mjj<150. && mggjj>270. && mggjj<330. && ptjj>") + ssJetPtCut.str() + TString(" && minDeltaR_gb>") + ssDeltaRbgCut.str() + TString(" && absCosThetaStar<") + ssAbsCosThetaStarCut.str() + TString(")"));
    // cutOpt0.push_back( TString("1.*weight*(btagCategory==2 && mjj>95. && mjj<150. && mggjj>270. && mggjj<330. && ptgg>") + ssGammaPtCut.str() + TString(" && minDeltaR_gb>") + ssDeltaRbgCut.str() + TString(" && absCosThetaStar<") + ssAbsCosThetaStarCut.str() + TString(")") );
    // cutOpt0.push_back( TString("1.*weight*(btagCategory==2 && mjj>95. && mjj<150. && mggjj>270. && mggjj<330. && ptgg>") + ssGammaPtCut.str() + TString(" && minDeltaR_gb>") + ssDeltaRbgCut.str() + TString(" && ptjj>") + ssJetPtCut.str() + TString(")") );
    cutOpt0.push_back( TString("1.*weight*(btagCategory==2 && mjj>95. && mjj<150. && mggjj>270. && mggjj<330. && ptgg>") + ssGammaPtCut.str() + TString(" && ptjj>") + ssJetPtCut.str() + TString(" && absCosThetaStar<") + ssAbsCosThetaStarCut.str() + TString(")") );

    // cutOpt0b.push_back( TString("0.59*pt_scaled_2D_weight_data*(btagCategory==1 && mjj>90. && mjj<170. && mggjj>270. && mggjj<330. && ptjj>") + ssJetPtCut.str() + TString(" && minDeltaR_gb>") + ssDeltaRbgCut.str() + TString(" && absCosThetaStar<") + ssAbsCosThetaStarCut.str() + TString(")"));
    // cutOpt0b.push_back( TString("0.59*pt_scaled_2D_weight_data*(btagCategory==1 && mjj>90. && mjj<170. && mggjj>270. && mggjj<330. && ptgg>") + ssGammaPtCut.str() + TString(" && minDeltaR_gb>") + ssDeltaRbgCut.str() + TString(" && absCosThetaStar<") + ssAbsCosThetaStarCut.str() + TString(")") );
    // cutOpt0b.push_back( TString("1*0.59*pt_scaled_2D_weight_data*(btagCategory==1 && mjj>90. && mjj<170. && mggjj>270. && mggjj<330. && ptgg>") + ssGammaPtCut.str() + TString(" && minDeltaR_gb>") + ssDeltaRbgCut.str() + TString(" && ptjj>") + ssJetPtCut.str() + TString(")") );
    cutOpt0b.push_back( TString("1*0.59*pt_scaled_2D_weight_data*(btagCategory==1 && mjj>90. && mjj<170. && mggjj>270. && mggjj<330. && ptgg>") + ssGammaPtCut.str() + TString(" && ptjj>") + ssJetPtCut.str() + TString(" && absCosThetaStar<") + ssAbsCosThetaStarCut.str() + TString(")") );    

    // cutOpt0b.push_back( TString("0.59*pt_scaled_2D_weight_data*(btagCategory==2 && mjj>95. && mjj<150. && mggjj>270. && mggjj<330. && ptjj>") + ssJetPtCut.str() + TString(" && minDeltaR_gb>") + ssDeltaRbgCut.str() + TString(" && absCosThetaStar<") + ssAbsCosThetaStarCut.str() + TString(")"));
    // cutOpt0b.push_back( TString("0.59*pt_scaled_2D_weight_data*(btagCategory==2 && mjj>95. && mjj<150. && mggjj>270. && mggjj<330. && ptgg>") + ssGammaPtCut.str() + TString(" && minDeltaR_gb>") + ssDeltaRbgCut.str() + TString(" && absCosThetaStar<") + ssAbsCosThetaStarCut.str() + TString(")") );
    // cutOpt0b.push_back( TString("1*0.59*pt_scaled_2D_weight_data*(btagCategory==2 && mjj>95. && mjj<150. && mggjj>270. && mggjj<330. && ptgg>") + ssGammaPtCut.str() + TString(" && minDeltaR_gb>") + ssDeltaRbgCut.str() + TString(" && ptjj>") + ssJetPtCut.str() + TString(")") );
    cutOpt0b.push_back( TString("1*0.59*pt_scaled_2D_weight_data*(btagCategory==2 && mjj>95. && mjj<150. && mggjj>270. && mggjj<330. && ptgg>") + ssGammaPtCut.str() + TString(" && ptjj>") + ssJetPtCut.str() + TString(" && absCosThetaStar<") + ssAbsCosThetaStarCut.str() + TString(")") );


    cutOpt1.push_back( TString("1.*weight*(btagCategory==1 && mjj>90. && mjj<170. && mggjj>270. && mggjj<330. && ptgg>") + ssGammaPtCut.str() + TString(" && minDeltaR_gb>") + ssDeltaRbgCut.str() + TString(" && ptjj>") + ssJetPtCut.str() + TString(" && absCosThetaStar<") + ssAbsCosThetaStarCut.str() + TString(")") );
    cutOpt1.push_back( TString("1.*weight*(btagCategory==2 && mjj>95. && mjj<150. && mggjj>270. && mggjj<330. && ptgg>") + ssGammaPtCut.str() + TString(" && minDeltaR_gb>") + ssDeltaRbgCut.str() + TString(" && ptjj>") + ssJetPtCut.str() + TString(" && absCosThetaStar<") + ssAbsCosThetaStarCut.str() + TString(")") );
    cutOpt1b.push_back( TString("1*0.59*pt_scaled_2D_weight_data*(btagCategory==1 && mjj>90. && mjj<170. && mggjj>270. && mggjj<330. && ptgg>") + ssGammaPtCut.str() + TString(" && minDeltaR_gb>") + ssDeltaRbgCut.str() + TString(" && ptjj>") + ssJetPtCut.str() + TString(" && absCosThetaStar<") + ssAbsCosThetaStarCut.str() + TString(")") );
    cutOpt1b.push_back( TString("1*0.59*pt_scaled_2D_weight_data*(btagCategory==2 && mjj>95. && mjj<150. && mggjj>270. && mggjj<330. && ptgg>") + ssGammaPtCut.str() + TString(" && minDeltaR_gb>") + ssDeltaRbgCut.str() + TString(" && ptjj>") + ssJetPtCut.str() + TString(" && absCosThetaStar<") + ssAbsCosThetaStarCut.str() + TString(")") );
    
    vector<TString> cutOptAll0,  cutOptAll1;
    vector<TString> cutOptAll0b, cutOptAll1b;
    for(int i=0;i<(int)cutOpt0.size();++i)  cutOptAll0.push_back(cutOpt0[i]);
    for(int i=0;i<(int)cutOpt1.size();++i)  cutOptAll1.push_back(cutOpt1[i]);
    for(int i=0;i<(int)cutOpt0b.size();++i) cutOptAll0b.push_back(cutOpt0b[i]);
    for(int i=0;i<(int)cutOpt1b.size();++i) cutOptAll1b.push_back(cutOpt1b[i]);

    cout << endl;

    // signal
    cout << "signal" << endl;
    pair<float,float> thisSpair_1btag = runOptimization(treeSig, cutOptAll0[0], cutOptAll1[0]);
    pair<float,float> thisSpair_2btag = runOptimization(treeSig, cutOptAll0[1], cutOptAll1[1]);
    float theS_1btag    = thisSpair_1btag.first;
    float theS_2btag    = thisSpair_2btag.first;
    float theSeff_1btag = thisSpair_1btag.second;
    float theSeff_2btag = thisSpair_2btag.second;

    // background
    cout << "background" << endl;
    pair<float,float> thisBpair_1btag, thisBpair_2btag;
    if (!isCS) {
      thisBpair_1btag = runOptimization(treeBkg, cutOptAll0[0], cutOptAll1[0]);
      thisBpair_2btag = runOptimization(treeBkg, cutOptAll0[1], cutOptAll1[1]);
    } else {
      thisBpair_1btag = runOptimization(treeBkg, cutOptAll0b[0], cutOptAll1b[0]);
      thisBpair_2btag = runOptimization(treeBkg, cutOptAll0b[1], cutOptAll1b[1]);
    }

    float theB_1btag    = thisBpair_1btag.first;
    float theB_2btag    = thisBpair_2btag.first;
    // to deal with low stat
    if (theB_1btag<1) { cout << "increasing back to 1" << endl; theB_1btag=1; }
    if (theB_2btag<1) { cout << "increasing back to 1" << endl; theB_2btag=1; }
    float theBeff_1btag = thisBpair_1btag.second;
    float theBeff_2btag = thisBpair_2btag.second;
    
    float theUL_1btag;
    float theUL_2btag;
    if (theMass==300) {
      theUL_1btag = UL(theS_1btag, theB_1btag, 0.000271);
      theUL_2btag = UL(theS_2btag, theB_2btag, 0.000271);
    } else if (theMass==500) {
      theUL_1btag = UL(theS_1btag, theB_1btag, 0.0000471);
      theUL_2btag = UL(theS_2btag, theB_2btag, 0.0000471);
    } else if (theMass==700) {
      theUL_1btag = UL(theS_1btag, theB_1btag, 0.0000158);
      theUL_2btag = UL(theS_2btag, theB_2btag, 0.0000158);
    } else if (theMass==1000) {
      theUL_1btag = UL(theS_1btag, theB_1btag, 0.00000263);
      theUL_2btag = UL(theS_2btag, theB_2btag, 0.00000263);
    } else {
      cout << "wrong mass. I don't know the xsec" << endl;
    }

    signEffMap_1btag -> SetBinContent(ii,theSeff_1btag);
    backEffMap_1btag -> SetBinContent(ii,theBeff_1btag);
    ulMap_1btag      -> SetBinContent(ii,theUL_1btag);
    
    signEffMap_2btag -> SetBinContent(ii,theSeff_2btag);
    backEffMap_2btag -> SetBinContent(ii,theBeff_2btag);
    ulMap_2btag      -> SetBinContent(ii,theUL_2btag);
  }



  TFile myFile("myFile.root","RECREATE");
  myFile.cd();
  signEffMap_1btag -> Write();
  signEffMap_2btag -> Write();
  backEffMap_1btag -> Write();
  backEffMap_2btag -> Write();
  ulMap_1btag      -> Write();
  ulMap_2btag      -> Write();
  myFile.Close();
}

