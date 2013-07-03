#include <TMath.h>
#include <iostream>
#include <fillPlot2012_radion_commonNtp.h>

struct RedntpDirStruct {
  
  std::string maindatadir;
  std::string mainmcdir;
  std::string datadir;
  std::string mcdir;
};

void finalize_oneDataset_commonNtp(const std::string& redntpVersion,  const std::string& dataset, const std::string& selectionType, const std::string& bTaggerType, std::vector<std::string> *datasets );

void do_haddCommand_commonNtp( const std::string& redntpVersion, const std::string& dataset, std::vector<std::string> *datasets, const std::string& selectionType, const std::string& bTaggerType );

std::string getSingleFileName_commonNtp( const std::string& redntpVersion, const std::string& dataset, const std::string& selectionType, const std::string& bTaggerType );

RedntpDirStruct get_dirs_commonNtp( const std::string& prodVersion );

int main( int argc, char* argv[] ) {

  if( argc!=4 && argc!=5 ) {
    std::cout << "USAGE: ./finalize2012_radion_commonNtp [redntpVersion] [dataset] [selectionType] [bTaggerType=JP or CSV]" <<std::endl;
    return 13;
  }

  std::string redntpVersion(argv[1]);
  std::string dataset(argv[2]);
  std::string selectionType(argv[3]);

  std::string bTaggerType="JP";
  if( argc>4 ) {
    std::string bTaggerType_str(argv[4]);
    bTaggerType = bTaggerType_str;
  }
  cout << "Running analysis for bTaggerType = " << bTaggerType << endl;
  
  std::vector<std::string> *datasets = new std::vector<std::string>;

  // chiara: qui andranno messi gli altri a mano a mano che arrivano 
  if( dataset=="Data2012" ) {  
    finalize_oneDataset_commonNtp(redntpVersion, "Data2012", selectionType, bTaggerType, datasets); 
    
  } else if( dataset=="DiPhotonBox_8TeV-pythia6" ) {
    // finalize_oneDataset_commonNtp(redntpVersion, "DiPhotonBox_Pt-10To25_8TeV-pythia6",   selectionType, bTaggerType, datasets);
    finalize_oneDataset_commonNtp(redntpVersion, "DiPhotonBox_Pt-250ToInf_8TeV-pythia", selectionType, bTaggerType, datasets);
    finalize_oneDataset_commonNtp(redntpVersion, "DiPhotonBox_Pt-25To250_8TeV-pythia",  selectionType, bTaggerType, datasets);

  } else if( dataset=="DiPhoton_8TeV" ) {

    // finalize_oneDataset_commonNtp(redntpVersion, "DiPhotonBox_Pt-10To25_8TeV-pythia6",   selectionType, bTaggerType, datasets);
    finalize_oneDataset_commonNtp(redntpVersion, "DiPhotonBox_Pt-250ToInf_8TeV-pythia6",  selectionType,  bTaggerType, datasets);
    finalize_oneDataset_commonNtp(redntpVersion, "DiPhotonBox_Pt-25To250_8TeV-pythia6",   selectionType,  bTaggerType, datasets);
    finalize_oneDataset_commonNtp(redntpVersion, "DiPhotonJets_8TeV-madgraph", selectionType,  bTaggerType, datasets);

  } else if( dataset=="DiPhoton_madgraph" ) {

    finalize_oneDataset_commonNtp(redntpVersion, "DiPhotonJets_8TeV-madgraph", selectionType,  bTaggerType, datasets);

  } else if( dataset=="V_8TeV") {
       
    finalize_oneDataset_commonNtp(redntpVersion, "DYJetsToLL",	selectionType, bTaggerType, datasets); 
    // finalize_oneDataset_commonNtp(redntpVersion, "WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1",      selectionType, bTaggerType, datasets);

  } else if( dataset=="GJet_doubleEMEnriched_TuneZ2star_8TeV-pythia6" ) {

    finalize_oneDataset_commonNtp(redntpVersion, "GJet_Pt-20to40_doubleEMEnriched_8TeV-pf", selectionType, bTaggerType, datasets);
    finalize_oneDataset_commonNtp(redntpVersion, "GJet_Pt-20to40_doubleEMEnriched_8TeV-pp", selectionType, bTaggerType, datasets);
    finalize_oneDataset_commonNtp(redntpVersion, "GJet_Pt40_doubleEMEnriched_8TeV-pf",      selectionType, bTaggerType, datasets);
    finalize_oneDataset_commonNtp(redntpVersion, "GJet_Pt40_doubleEMEnriched_8TeV-pp",      selectionType, bTaggerType, datasets);

  } else if( dataset=="QCD_doubleEMEnriched_TuneZ2star_8TeV-pythia6") {

    finalize_oneDataset_commonNtp(redntpVersion, "QCD_Pt-30to40_doubleEMEnriched_8TeV-ff", selectionType, bTaggerType, datasets);
    finalize_oneDataset_commonNtp(redntpVersion, "QCD_Pt-30to40_doubleEMEnriched_8TeV-pf", selectionType, bTaggerType, datasets);
    finalize_oneDataset_commonNtp(redntpVersion, "QCD_Pt-40_doubleEMEnriched_8TeV-ff", selectionType, bTaggerType, datasets);
    finalize_oneDataset_commonNtp(redntpVersion, "QCD_Pt-40_doubleEMEnriched_8TeV-pf", selectionType, bTaggerType, datasets);

  } else if( dataset=="TT_8TeV" ) {
       
    finalize_oneDataset_commonNtp(redntpVersion, "TTbarGG", selectionType, bTaggerType, datasets);
    finalize_oneDataset_commonNtp(redntpVersion, "TTGJets_8TeV-madgraph", selectionType, bTaggerType, datasets);
    finalize_oneDataset_commonNtp(redntpVersion, "TTJets_TuneZ2star_8TeV-madgraph", selectionType, bTaggerType, datasets); 

  } else if( dataset=="VV_8TeV" ) {

    finalize_oneDataset_commonNtp(redntpVersion, "ZZ_TuneZ2star_8TeV", selectionType, bTaggerType, datasets);
    // finalize_oneDataset_commonNtp(redntpVersion, "WW_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1", selectionType, bTaggerType, datasets); 
    // finalize_oneDataset_commonNtp(redntpVersion, "WZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1", selectionType, bTaggerType, datasets);
    finalize_oneDataset_commonNtp(redntpVersion, "WGToLNuG_8TeV", selectionType, bTaggerType, datasets);
    // finalize_oneDataset_commonNtp(redntpVersion, "ZGToLLG_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1",             selectionType, bTaggerType, datasets);

  } else if( dataset=="VGG_8TeV" ) {
    
    // finalize_oneDataset_commonNtp(redntpVersion, "WmGG", selectionType, bTaggerType, datasets); 
    // finalize_oneDataset_commonNtp(redntpVersion, "WpGG", selectionType, bTaggerType, datasets);
    finalize_oneDataset_commonNtp(redntpVersion, "ZGG",  selectionType, bTaggerType, datasets);  

  } else if( dataset=="HToGG_M-125_8TeV-pythia6" ) {

    finalize_oneDataset_commonNtp(redntpVersion, "WH_ZH_HToGG_M-125_8TeV",    selectionType, bTaggerType, datasets);
    finalize_oneDataset_commonNtp(redntpVersion, "GluGluToHToGG_M-125_8TeV",  selectionType, bTaggerType, datasets);  
    finalize_oneDataset_commonNtp(redntpVersion, "VBF_HToGG_M-125_8TeV",      selectionType, bTaggerType, datasets);  
    finalize_oneDataset_commonNtp(redntpVersion, "TTH_HToGG_M-125_8TeV",      selectionType, bTaggerType, datasets);

  } else if( dataset=="WH_ZH_HToGG_M-125_8TeV" ) {

    finalize_oneDataset_commonNtp(redntpVersion, "WH_ZH_HToGG_M-125_8TeV",    selectionType, bTaggerType, datasets);

  } else if( dataset=="WH_HToGG_M-125_8TeV" ) {

    finalize_oneDataset_commonNtp(redntpVersion, "WH_HToGG_M-125_8TeV",    selectionType, bTaggerType, datasets);

  } else if( dataset=="ZH_HToGG_M-125_8TeV" ) {

    finalize_oneDataset_commonNtp(redntpVersion, "ZH_HToGG_M-125_8TeV",    selectionType, bTaggerType, datasets);

  } else if( dataset=="VBF_HToGG_M-125_8TeV" ) {

    finalize_oneDataset_commonNtp(redntpVersion, "VBF_HToGG_M-125_8TeV",    selectionType, bTaggerType, datasets);

  } else if( dataset=="GluGluToHToGG_M-125_8TeV" ) {

    finalize_oneDataset_commonNtp(redntpVersion, "GluGluToHToGG_M-125_8TeV",  selectionType, bTaggerType, datasets);

  } else if( dataset=="TTH_HToGG_M-125_8TeV" ) {

    finalize_oneDataset_commonNtp(redntpVersion, "TTH_HToGG_M-125_8TeV",  selectionType, bTaggerType, datasets);

  } else if( dataset=="Radion_M-300_madgraph" ) {  

    finalize_oneDataset_commonNtp(redntpVersion, "Radion_M-300_madgraph", selectionType, bTaggerType, datasets);    

  } else if( dataset=="Radion_M-500_madgraph" ) {  

    finalize_oneDataset_commonNtp(redntpVersion, "Radion_M-500_madgraph", selectionType, bTaggerType, datasets);    

  } else if( dataset=="Radion_M-700_madgraph" ) {  

    finalize_oneDataset_commonNtp(redntpVersion, "Radion_M-700_madgraph", selectionType, bTaggerType, datasets);    

  } else if( dataset=="Radion_M-1000_madgraph" ) {  

    finalize_oneDataset_commonNtp(redntpVersion, "Radion_M-1000_madgraph", selectionType, bTaggerType, datasets);    

  } else if( dataset=="Radion_M-300_regr" ) {  

    finalize_oneDataset_commonNtp(redntpVersion, "Radion_M-300_regr", selectionType, bTaggerType, datasets);    

  } else {
    
    cout << "this dataset does not exists" << endl;
  }
  
  do_haddCommand_commonNtp(redntpVersion, dataset, datasets, selectionType, bTaggerType );

  return 0;
}

void finalize_oneDataset_commonNtp( const std::string& redntpProdVersion, const std::string& dataset, const std::string& selectionType, const std::string& bTaggerType, std::vector<std::string> *datasets ) {

  TString dataset_tstr(dataset);
  
  std::cout << std::endl << std::endl << std::endl << "####     Finalizing " << dataset << std::endl;
  std::cout << "####  Selection: " << selectionType << std::endl;
  std::cout << "####  b-Tagger: "  << bTaggerType << std::endl << std::endl;

  RedntpDirStruct dirs = get_dirs_commonNtp( redntpProdVersion );

  fillPlot2012_radion_commonNtp* rf = new fillPlot2012_radion_commonNtp( dataset, selectionType, bTaggerType );    

  bool isData = ( dataset_tstr.Contains("2011") || dataset_tstr.Contains("2012") ); 

  std::string redNtpDir;
  if( isData ) {
    redNtpDir = dirs.maindatadir;
    redNtpDir = redNtpDir + "/" + dirs.datadir;
    cout << "Running on redNtpDir = " << redNtpDir << endl;
  } else {
    redNtpDir = dirs.mainmcdir;
    redNtpDir = redNtpDir + "/" + dirs.mcdir;
    cout << "Running on redNtpDir = " << redNtpDir << endl;
  }

  rf->set_redNtpDir(redNtpDir);   
  rf->set_outputDir("finalizedTrees_"+redntpProdVersion);    
  rf->addFile(dataset);    
  rf->finalize();
  delete rf;

  datasets->push_back(dataset);
}

void do_haddCommand_commonNtp( const std::string& redntpVersion, const std::string& dataset, std::vector<std::string> *datasets, const std::string& selectionType, const std::string& bTaggerType ) {

  if( datasets->size()<=1 ) return;

  std::string hadd_command = "hadd -f " + getSingleFileName_commonNtp( redntpVersion, dataset, selectionType, bTaggerType );

  std::string suffix = "";
  for( unsigned i=0; i<datasets->size(); ++i )
    suffix = suffix + " " + getSingleFileName_commonNtp( redntpVersion, datasets->at(i), selectionType, bTaggerType );

  hadd_command += suffix;
  std::string rm_command = "rm " + suffix;
  system(hadd_command.c_str());
  if( dataset!="HToGG_M-125_8TeV" ) //keep also separate higgs processes to quantify VH                                                                          
    system(rm_command.c_str());
}

std::string getSingleFileName_commonNtp( const std::string& redntpVersion, const std::string& dataset, const std::string& selectionType, const std::string& bTaggerType ) {

  std::string fileName = "finalizedTrees_" + redntpVersion + "/Radion_" + dataset + "_" + selectionType + "_" + bTaggerType + ".root";
  return fileName;
}

RedntpDirStruct get_dirs_commonNtp( const std::string& prodVersion ) {

  RedntpDirStruct returnStruct;

  if( prodVersion=="Radion" ) {

    // returnStruct.maindatadir = "/xrootdfs/cms/local/crovelli/Radion/common";  
    // returnStruct.mainmcdir   = "/xrootdfs/cms/local/crovelli/Radion/common";
    // returnStruct.datadir     = "radion_tree_v03";
    // returnStruct.mcdir       = "radion_tree_v03";

    returnStruct.maindatadir = "/afs/cern.ch/work/c/crovelli/Radion/CMSSW_5_3_6/src/HHggbb";
    returnStruct.mainmcdir   = "/afs/cern.ch/work/c/crovelli/Radion/CMSSW_5_3_6/src/HHggbb";
    returnStruct.datadir     = "dataCommonNtp";
    returnStruct.mcdir       = "dataCommonNtp";

  } else {

    std::cout << "-> Unknown prodVersion: '" << prodVersion << "'! Exiting." << std::endl;
    exit(11);
  }

  return returnStruct;
}
