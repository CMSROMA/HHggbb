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

#include <TH2.h>
#include <TFile.h>
#include <TTree.h>

#include "cl95cms.C"

using namespace std;

double UL(double expsig, double expbkg){

  double prob=0.;
  TF1 *mypoisson = new TF1("mypoisson","TMath::Poisson(x,[0])",0,10000);
  double ul(0);
  
  for(int i=0; i<100000; i++){
    double exp = i*(0.0002*expsig) + expbkg;
    mypoisson->SetParameter(0, exp);
    prob = mypoisson->Integral(0,expbkg-0.5);
    if(prob<0.05){
      ul = exp-expbkg;
      break;
    }
  }
  return ul/expsig;
}

// chiara: cross-sections are for mH=300
double UL2(double expsig, double expbkg){

  float effS = expsig / (19600.*0.000271);
  float ul = CLA( 19600., 0., effS, 0., expbkg, 0. );                                                                                            
  return ul/0.000271;
 }

double probability_disc(double obs, double exp){

  double prob=0.;
  TF1 *mypoisson = new TF1("mypoisson","TMath::Poisson(x,[0])",0,10000);
  mypoisson->SetParameter(0, exp);
  if (exp<100)
    prob= mypoisson->Integral(obs-0.5,999);
  else
    prob= mypoisson->Integral(obs-0.5,9999);
  return TMath::NormQuantile(1-prob);
}

pair<float,float> runOptimization(TTree *theTree, TString cut0, TString cut1) {

  TH1F *mjj0 = new TH1F("mjj0", "", 100, 0.0, 10000.0);
  TH1F *mjj1 = new TH1F("mjj1", "", 100, 0.0, 10000.0);
  
  theTree->Project("mjj0","mjj", cut0);
  theTree->Project("mjj1","mjj", cut1);

  float afterCut   = mjj1->Integral();
  float efficiency = mjj1->Integral()/mjj0->Integral();

  cout << "cut = " << cut1 << ", eff = " << efficiency << ", mjj->Integral() = " << mjj1->Integral() << endl;

  return std::make_pair<float,float>(afterCut, efficiency);  
}

int main(int argc, char* argv[]) {

  TFile *fileSig, *fileBkg;
  TTree *treeSig, *treeBkg;
  fileBkg = TFile::Open("finalizedTrees_Radion_presel/Radion_Data2012_default_CSV.root");
  fileSig = TFile::Open("finalizedTrees_Radion_presel/Radion_Radion_M-300_madgraph_default_CSV.root");

  if( fileSig && fileBkg) {
    fileBkg->cd();
    treeBkg = (TTree*)fileBkg->Get("myTrees");
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
  cout << endl;

  // to run cuts optimization 
  /*
  // chiara: for mjj
  TH2F *signEffMap_1btag = new TH2F ("signEffMap_1btag", "",10, 30., 80., 10, 130., 180.); 
  TH2F *backEffMap_1btag = new TH2F ("backEffMap_1btag", "",10, 30., 80., 10, 130., 180.); 
  TH2F *ulMap_1btag      = new TH2F ("ulMap_1btag",      "",10, 30., 80., 10, 130., 180.); 
  TH2F *signifMap_1btag  = new TH2F ("signifMap_1btag",  "",10, 30., 80., 10, 130., 180.); 
  TH2F *SsqrtBMap_1btag  = new TH2F ("SsqrtBMap_1btag",  "",10, 30., 80., 10, 130., 180.); 
  TH2F *SoverBMap_1btag  = new TH2F ("SoverBMap_1btag",  "",10, 30., 80., 10, 130., 180.); 

  TH2F *signEffMap_2btag = new TH2F ("signEffMap_2btag", "",10, 50., 100., 10, 115., 165.); 
  TH2F *backEffMap_2btag = new TH2F ("backEffMap_2btag", "",10, 50., 100., 10, 115., 165.); 
  TH2F *ulMap_2btag      = new TH2F ("ulMap_2btag",      "",10, 50., 100., 10, 115., 165.); 
  TH2F *signifMap_2btag  = new TH2F ("signifMap_2btag",  "",10, 50., 100., 10, 115., 165.); 
  TH2F *SsqrtBMap_2btag  = new TH2F ("SsqrtBMap_2btag",  "",10, 50., 100., 10, 115., 165.); 
  TH2F *SoverBMap_2btag  = new TH2F ("SoverBMap_2btag",  "",10, 50., 100., 10, 115., 165.); 
  */

  // chiara: for mggjj 
  TH2F *signEffMap_1btag = new TH2F ("signEffMap_1btag", "",10, 230., 280., 10, 330., 380.); 
  TH2F *backEffMap_1btag = new TH2F ("backEffMap_1btag", "",10, 230., 280., 10, 330., 380.); 
  TH2F *ulMap_1btag      = new TH2F ("ulMap_1btag",      "",10, 230., 280., 10, 330., 380.); 
  TH2F *signifMap_1btag  = new TH2F ("signifMap_1btag",  "",10, 230., 280., 10, 330., 380.); 
  TH2F *SsqrtBMap_1btag  = new TH2F ("SsqrtBMap_1btag",  "",10, 230., 280., 10, 330., 380.); 
  TH2F *SoverBMap_1btag  = new TH2F ("SoverBMap_1btag",  "",10, 230., 280., 10, 330., 380.); 

  TH2F *signEffMap_2btag = new TH2F ("signEffMap_2btag", "",10, 250., 300., 10, 315., 365.); 
  TH2F *backEffMap_2btag = new TH2F ("backEffMap_2btag", "",10, 250., 300., 10, 315., 365.); 
  TH2F *ulMap_2btag      = new TH2F ("ulMap_2btag",      "",10, 250., 300., 10, 315., 365.); 
  TH2F *signifMap_2btag  = new TH2F ("signifMap_2btag",  "",10, 250., 300., 10, 315., 365.); 
  TH2F *SsqrtBMap_2btag  = new TH2F ("SsqrtBMap_2btag",  "",10, 250., 300., 10, 315., 365.); 
  TH2F *SoverBMap_2btag  = new TH2F ("SoverBMap_2btag",  "",10, 250., 300., 10, 315., 365.); 

  for (int mjjInf=0; mjjInf<10; mjjInf++) {
    for (int mjjSup=0; mjjSup<10; mjjSup++) {
      
      // chiara
      float mjjInfCutInf_1 = 230. + mjjInf*5.;
      float mjjInfCutSup_1 = 230. + (mjjInf+1.)*5.;
      float mjjSupCutInf_1 = 330. + mjjSup*5.;
      float mjjSupCutSup_1 = 330. + (mjjSup+1.)*5.;

      float mjjInfCutInf_2 = 250. + mjjInf*5.;
      float mjjInfCutSup_2 = 250. + (mjjInf+1.)*5.;
      float mjjSupCutInf_2 = 315. + mjjSup*5.;
      float mjjSupCutSup_2 = 315. + (mjjSup+1.)*5.;

      float binMjjInf_1 = (mjjInfCutInf_1+mjjInfCutSup_1)/2.;  
      float binMjjSup_1 = (mjjSupCutInf_1+mjjSupCutSup_1)/2.;  

      cout << "1btag: mjjInf = " << mjjInf << ", mjjInfCutInf_1 = " << mjjInfCutInf_1 << ", mjjInfCutSup_1 = " << mjjInfCutSup_1 << ", binMjjInf_1 = " << binMjjInf_1 << endl;  
      cout << "1btag: mjjSup = " << mjjSup << ", mjjSupCutInf_1 = " << mjjSupCutInf_1 << ", mjjSupCutSup_1 = " << mjjSupCutSup_1 << ", binMjjSup_1 = " << binMjjSup_1 << endl; 

      float binMjjInf_2 = (mjjInfCutInf_2+mjjInfCutSup_2)/2.;  
      float binMjjSup_2 = (mjjSupCutInf_2+mjjSupCutSup_2)/2.;  

      cout << endl;
      cout << "2btag: mjjInf = " << mjjInf << ", mjjInfCutInf_2 = " << mjjInfCutInf_2 << ", mjjInfCutSup_2 = " << mjjInfCutSup_2 << ", binMjjInf_2 = " << binMjjInf_2 << endl;  
      cout << "2btag: mjjSup = " << mjjSup << ", mjjSupCutInf_2 = " << mjjSupCutInf_2 << ", mjjSupCutSup_2 = " << mjjSupCutSup_2 << ", binMjjSup_2 = " << binMjjSup_2 << endl; 

      stringstream ssMjjInfCutInf_1 (stringstream::in | stringstream::out);
      ssMjjInfCutInf_1 << mjjInfCutInf_1;
      stringstream ssMjjSupCutInf_1 (stringstream::in | stringstream::out);
      ssMjjSupCutInf_1 << mjjSupCutInf_1;

      stringstream ssMjjInfCutInf_2 (stringstream::in | stringstream::out);
      ssMjjInfCutInf_2 << mjjInfCutInf_2;
      stringstream ssMjjSupCutInf_2 (stringstream::in | stringstream::out);
      ssMjjSupCutInf_2 << mjjSupCutInf_2;

      /*
      // chiara: for mjj
      vector<TString> cutOpt0;
      cutOpt0.push_back( TString("1.*weight*(btagCategory==1)") );
      cutOpt0.push_back( TString("1.*weight*(btagCategory==2)") );

      vector<TString> cutOpt1;
      cutOpt1.push_back( TString("1.*weight*(btagCategory==1 && mjj>") + ssMjjInfCutInf_1.str() + TString(" && mjj<") + ssMjjSupCutInf_1.str() +TString(")"));
      cutOpt1.push_back( TString("1.*weight*(btagCategory==2 && mjj>") + ssMjjInfCutInf_2.str() + TString(" && mjj<") + ssMjjSupCutInf_2.str() +TString(")"));
      */

      // chiara: for mggjj
      vector<TString> cutOpt0;
      cutOpt0.push_back( TString("1.*weight*(btagCategory==1 && mjj>70 && mjj<165)") );
      cutOpt0.push_back( TString("1.*weight*(btagCategory==2 && mjj>95 && mjj<140)") );

      vector<TString> cutOpt1;
      cutOpt1.push_back( TString("1.*weight*(btagCategory==1 && mjj>70 && mjj<165 && mggjj_kin>") + ssMjjInfCutInf_1.str() + TString(" && mggjj_kin<") + ssMjjSupCutInf_1.str() +TString(")"));
      cutOpt1.push_back( TString("1.*weight*(btagCategory==2 && mjj>95 && mjj<140 && mggjj_kin>") + ssMjjInfCutInf_2.str() + TString(" && mggjj_kin<") + ssMjjSupCutInf_2.str() +TString(")"));


      vector<TString> cutOptAll0;
      for(int i=0;i<(int)cutOpt0.size();++i) cutOptAll0.push_back(cutOpt0[i]);

      vector<TString> cutOptAll1;
      for(int i=0;i<(int)cutOpt1.size();++i) cutOptAll1.push_back(cutOpt1[i]);

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
      pair<float,float> thisBpair_1btag = runOptimization(treeBkg, cutOptAll0[0], cutOptAll1[0]);
      pair<float,float> thisBpair_2btag = runOptimization(treeBkg, cutOptAll0[1], cutOptAll1[1]);
      float theB_1btag    = thisBpair_1btag.first;
      float theB_2btag    = thisBpair_2btag.first;
      float theBeff_1btag = thisBpair_1btag.second;
      float theBeff_2btag = thisBpair_2btag.second;

      float theUL2_1btag = UL2(theS_1btag, theB_1btag);      
      float theUL2_2btag = UL2(theS_2btag, theB_2btag);      
      cout << "UL2: " << theUL2_1btag << " " << theUL2_2btag << endl;  
      float theSignif_1btag = probability_disc(theS_1btag+theB_1btag, theB_1btag);   // 1 - expected p-value
      float theSignif_2btag = probability_disc(theS_2btag+theB_2btag, theB_2btag);   // 1 - expected p-value
      float theSsqrtB_1btag = theS_1btag/sqrt(theB_1btag);
      float theSsqrtB_2btag = theS_2btag/sqrt(theB_2btag);
      float theSoverB_1btag = theS_1btag/theB_1btag;
      float theSoverB_2btag = theS_2btag/theB_2btag;

      signEffMap_1btag -> Fill(binMjjInf_1,binMjjSup_1,theSeff_1btag);
      backEffMap_1btag -> Fill(binMjjInf_1,binMjjSup_1,theBeff_1btag);
      signifMap_1btag  -> Fill(binMjjInf_1,binMjjSup_1,theSignif_1btag);
      ulMap_1btag      -> Fill(binMjjInf_1,binMjjSup_1,theUL2_1btag);
      SsqrtBMap_1btag  -> Fill(binMjjInf_1,binMjjSup_1,theSsqrtB_1btag);
      SoverBMap_1btag  -> Fill(binMjjInf_1,binMjjSup_1,theSoverB_1btag);

      signEffMap_2btag -> Fill(binMjjInf_2,binMjjSup_2,theSeff_2btag);
      backEffMap_2btag -> Fill(binMjjInf_2,binMjjSup_2,theBeff_2btag);
      signifMap_2btag  -> Fill(binMjjInf_2,binMjjSup_2,theSignif_2btag);
      ulMap_2btag      -> Fill(binMjjInf_2,binMjjSup_2,theUL2_2btag);
      SsqrtBMap_2btag  -> Fill(binMjjInf_2,binMjjSup_2,theSsqrtB_2btag);
      SoverBMap_2btag  -> Fill(binMjjInf_2,binMjjSup_2,theSoverB_2btag);
    }
  }


  TFile myFile("myFile.root","RECREATE");
  myFile.cd();
  signEffMap_1btag -> Write();
  signEffMap_2btag -> Write();
  backEffMap_1btag -> Write();
  backEffMap_2btag -> Write();
  signifMap_1btag  -> Write();
  signifMap_2btag  -> Write();
  ulMap_1btag      -> Write();
  ulMap_2btag      -> Write();
  SsqrtBMap_1btag  -> Write();
  SsqrtBMap_2btag  -> Write();
  SoverBMap_1btag  -> Write();
  SoverBMap_2btag  -> Write();
  myFile.Close();
}

