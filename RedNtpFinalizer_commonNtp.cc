#include <iostream>
#include "TString.h"
#include "RedNtpFinalizer_commonNtp.h"

RedNtpFinalizer_commonNtp::RedNtpFinalizer_commonNtp( const std::string& analyzerType, const std::string& dataset, const std::string& flags) {

  DEBUG_ = false;

  // chiara: aggiungi gli altri quando ci saranno
  TString dataset_tstr(dataset);

  if      ( dataset_tstr.Contains("Radion_M-300_madgraph") )    tree_ = new TChain("Radion_m300_8TeV_nm");   
  else if ( dataset_tstr.Contains("Radion_M-500_madgraph") )    tree_ = new TChain("Radion_m500_8TeV_nm");   
  else if ( dataset_tstr.Contains("Radion_M-700_madgraph") )    tree_ = new TChain("Radion_m700_8TeV_nm");   
  else if ( dataset_tstr.Contains("Radion_M-1000_madgraph") )   tree_ = new TChain("Radion_m1000_8TeV_nm");   

  else if ( dataset_tstr.Contains("GluGluToHToGG_M-125_8TeV") ) tree_ = new TChain("ggh_m125_8TeV");   
  else if ( dataset_tstr.Contains("VBF_HToGG_M-125_8TeV")  )    tree_ = new TChain("vbf_m125_8TeV");   
  else if ( dataset_tstr.Contains("WH_ZH_HToGG_M-125_8TeV")  )  tree_ = new TChain("wzh_m125_8TeV");   
  else if ( dataset_tstr.Contains("WH_HToGG_M-125_8TeV")  )     tree_ = new TChain("wzh_m125_8TeV");   
  else if ( dataset_tstr.Contains("ZH_HToGG_M-125_8TeV")  )     tree_ = new TChain("wzh_m125_8TeV");   
  else if ( dataset_tstr.Contains("TTH_HToGG_M-125_8TeV")  )    tree_ = new TChain("tth_m125_8TeV");   

  else if( dataset_tstr.Contains("DiPhotonJets_8TeV-madgraph") ) tree_ = new TChain("diphojet_8TeV");

  else if( dataset_tstr.Contains("DiPhotonBox_Pt-250ToInf_8TeV-pythia") ) tree_ = new TChain("dipho_Box_250_8TeV");
  else if( dataset_tstr.Contains("DiPhotonBox_Pt-25To250_8TeV-pythia") )  tree_ = new TChain("dipho_Box_25_8TeV");

  else if( dataset_tstr.Contains("GJet_Pt-20to40_doubleEMEnriched_8TeV-pf") ) tree_ = new TChain("gjet_20_8TeV_pf");
  else if( dataset_tstr.Contains("GJet_Pt-20to40_doubleEMEnriched_8TeV-pp") ) tree_ = new TChain("gjet_20_8TeV_pp");
  else if( dataset_tstr.Contains("GJet_Pt40_doubleEMEnriched_8TeV-pf") )      tree_ = new TChain("gjet_40_8TeV_pf");  
  else if( dataset_tstr.Contains("GJet_Pt40_doubleEMEnriched_8TeV-pp") )      tree_ = new TChain("gjet_40_8TeV_pp");  

  else if( dataset_tstr.Contains("QCD_Pt-30to40_doubleEMEnriched_8TeV-ff") )  tree_ = new TChain("qcd_30_8TeV_ff");
  else if( dataset_tstr.Contains("QCD_Pt-30to40_doubleEMEnriched_8TeV-pf") )  tree_ = new TChain("qcd_30_8TeV_pf");
  else if( dataset_tstr.Contains("QCD_Pt-40_doubleEMEnriched_8TeV-ff") )      tree_ = new TChain("qcd_40_8TeV_ff");
  else if( dataset_tstr.Contains("QCD_Pt-40_doubleEMEnriched_8TeV-pf") )      tree_ = new TChain("qcd_40_8TeV_pf");

  else if( dataset_tstr.Contains("ZZ_TuneZ2star_8TeV") ) tree_ = new TChain("ZZJetsTo2L2Q");                   

  else if( dataset_tstr.Contains("TTbarGG") )                         tree_ = new TChain("ttgg_dR02");
  else if( dataset_tstr.Contains("TTJets_TuneZ2star_8TeV-madgraph") ) tree_ = new TChain("TTJets");
  else if( dataset_tstr.Contains("TTGJets_8TeV-madgraph") )           tree_ = new TChain("TTGJets");

  else if( dataset_tstr.Contains("DYJetsToLL") ) tree_ = new TChain("DYJetsToLL"); 

  else if( dataset_tstr.Contains("WGToLNuG") ) tree_ = new TChain("WGToLNuG"); 

  else if( dataset_tstr.Contains("ZGG") ) tree_ = new TChain("Zgg_dR02");

  else if( dataset_tstr.Contains("Data2012") ) tree_ = new TChain("Data");

  // to study the regression
  else if ( dataset_tstr.Contains("Radion_M-300_regr") ) tree_ = new TChain("Radion_m300_8TeV_nm");   

  else cout << "wrong dataset name" << endl;

  analyzerType_ = analyzerType;
  redNtpDir_ = "";
  dataset_ = dataset;
  flags_ = flags;
  outFile_ = 0;
  outputDir_ = "";
} 

RedNtpFinalizer_commonNtp::~RedNtpFinalizer_commonNtp() {

  std::cout << std::endl << "-> Histograms saved in: " << outFile_->GetName() << std::endl;
  
  if( tree_!=0 ) {
    delete tree_;
    tree_=0;
  }
} 

void RedNtpFinalizer_commonNtp::createOutputFile( const std::string& additionalFlags ) {

   std::string outfileName;

   if( DEBUG_ ) outfileName = "prova_"+dataset_;
   else {
    if(dataset_!="") outfileName =  analyzerType_ + "_" + dataset_;
    else outfileName = analyzerType_;
   }
   
   if( flags_!="" )
     outfileName = outfileName + "_" + flags_;
   if( additionalFlags!="" )
     outfileName = outfileName + "_" + additionalFlags;
   
   outfileName = outfileName + ".root";

   if( outputDir_!="" ) {
     outfileName = outputDir_ + "/" + outfileName;
     std::string mkdir_command = "mkdir -p " + outputDir_;
     system( mkdir_command.c_str() );
   }

   outFile_ = TFile::Open(outfileName.c_str(), "RECREATE");
   
   outFile_->cd();
}

void RedNtpFinalizer_commonNtp::addFile(const std::string& dataset, const std::string& selection) {

  // std::string infileName;
  TString dataset_tstr(dataset);

  // chiara: aggiungi gli altri quando ci saranno
  std::string treeName;
  if      ( dataset_tstr.Contains("Radion_M-300_madgraph") )  treeName = redNtpDir_ + "/Graviton_Radion-nm.root/Radion_m300_8TeV_nm";   
  else if ( dataset_tstr.Contains("Radion_M-500_madgraph") )  treeName = redNtpDir_ + "/Graviton_Radion-nm.root/Radion_m500_8TeV_nm";   
  else if ( dataset_tstr.Contains("Radion_M-700_madgraph") )  treeName = redNtpDir_ + "/Graviton_Radion-nm.root/Radion_m700_8TeV_nm";   
  else if ( dataset_tstr.Contains("Radion_M-1000_madgraph") ) treeName = redNtpDir_ + "/Graviton_Radion-nm.root/Radion_m1000_8TeV_nm";   

  else if ( dataset_tstr.Contains("GluGluToHToGG_M-125_8TeV")  ) treeName = redNtpDir_ + "/SMHiggs_m125.root/ggh_m125_8TeV";   
  else if ( dataset_tstr.Contains("VBF_HToGG_M-125_8TeV")  )     treeName = redNtpDir_ + "/SMHiggs_m125.root/vbf_m125_8TeV";   
  else if ( dataset_tstr.Contains("WH_ZH_HToGG_M-125_8TeV")  )   treeName = redNtpDir_ + "/SMHiggs_m125.root/wzh_m125_8TeV";   
  else if ( dataset_tstr.Contains("WH_HToGG_M-125_8TeV")  )      treeName = redNtpDir_ + "/SMHiggs_m125_WH.root/wzh_m125_8TeV";   
  else if ( dataset_tstr.Contains("ZH_HToGG_M-125_8TeV")  )      treeName = redNtpDir_ + "/SMHiggs_m125_ZH.root/wzh_m125_8TeV";   
  else if ( dataset_tstr.Contains("TTH_HToGG_M-125_8TeV")  )     treeName = redNtpDir_ + "/SMHiggs_m125.root/tth_m125_8TeV";   

  else if( dataset_tstr.Contains("DiPhotonJets_8TeV-madgraph") ) treeName = redNtpDir_ + "/MCBackgrounds_DiPhotonJets-madgraph.root/diphojet_8TeV";

  else if( dataset_tstr.Contains("DiPhotonBox_Pt-250ToInf_8TeV-pythia") ) treeName = redNtpDir_ + "/MCBackgrounds.root/dipho_Box_250_8TeV";
  else if( dataset_tstr.Contains("DiPhotonBox_Pt-25To250_8TeV-pythia") )  treeName = redNtpDir_ + "/MCBackgrounds.root/dipho_Box_25_8TeV"; 

  else if( dataset_tstr.Contains("GJet_Pt-20to40_doubleEMEnriched_8TeV-pf") ) treeName = redNtpDir_ + "/MCBackgrounds.root/gjet_20_8TeV_pf";
  else if( dataset_tstr.Contains("GJet_Pt-20to40_doubleEMEnriched_8TeV-pp") ) treeName = redNtpDir_ + "/MCBackgrounds.root/gjet_20_8TeV_pp";
  else if( dataset_tstr.Contains("GJet_Pt40_doubleEMEnriched_8TeV-pf") )      treeName = redNtpDir_ + "/MCBackgrounds.root/gjet_40_8TeV_pf";
  else if( dataset_tstr.Contains("GJet_Pt40_doubleEMEnriched_8TeV-pp") )      treeName = redNtpDir_ + "/MCBackgrounds.root/gjet_40_8TeV_pp";

  else if( dataset_tstr.Contains("QCD_Pt-30to40_doubleEMEnriched_8TeV-ff") ) treeName = redNtpDir_ + "/MCBackgrounds.root/qcd_30_8TeV_ff";
  else if( dataset_tstr.Contains("QCD_Pt-30to40_doubleEMEnriched_8TeV-pf") ) treeName = redNtpDir_ + "/MCBackgrounds.root/qcd_30_8TeV_pf";  
  else if( dataset_tstr.Contains("QCD_Pt-40_doubleEMEnriched_8TeV-ff") )     treeName = redNtpDir_ + "/MCBackgrounds.root/qcd_40_8TeV_ff";
  else if( dataset_tstr.Contains("QCD_Pt-40_doubleEMEnriched_8TeV-pf") )     treeName = redNtpDir_ + "/MCBackgrounds.root/qcd_40_8TeV_pf";

  else if( dataset_tstr.Contains("ZZ_TuneZ2star_8TeV") ) treeName = redNtpDir_ + "/MCBackgrounds.root/ZZJetsTo2L2Q";                   

  else if( dataset_tstr.Contains("TTbarGG") )                         treeName = redNtpDir_ + "/MCBackgrounds_ttgg_Zgg.root/ttgg_dR02";
  else if( dataset_tstr.Contains("TTJets_TuneZ2star_8TeV-madgraph") ) treeName = redNtpDir_ + "/MCBackgrounds.root/TTJets";
  else if( dataset_tstr.Contains("TTGJets_8TeV-madgraph") )           treeName = redNtpDir_ + "/MCBackgrounds.root/TTGJets";

  else if( dataset_tstr.Contains("DYJetsToLL") ) treeName = redNtpDir_ + "/MCBackgrounds_DY_WG.root/DYJetsToLL"; 

  else if( dataset_tstr.Contains("WGToLNuG") ) treeName = redNtpDir_ + "/MCBackgrounds_DY_WG.root/WGToLNuG"; 

  else if( dataset_tstr.Contains("ZGG") ) treeName = redNtpDir_ + "/MCBackgrounds_ttgg_Zgg.root/Zgg_dR02";

  else if( dataset_tstr.Contains("Data2012") ) treeName = redNtpDir_ + "/Data_full2012.root/Data";
  
  // to study regression
  else if( dataset_tstr.Contains("Radion_M-300_regr") ) treeName = redNtpDir_ + "/Radion_m300_withRegression.root/Radion_m300_8TeV_nm";   

  else cout << "wrong dataset name" << endl;
  
  tree_->Add(treeName.c_str());
  std::cout << "-> Added " << treeName << ". Tree now has " << tree_->GetEntries() << " entries." << std::endl;
}

