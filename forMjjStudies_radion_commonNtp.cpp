#include <TMath.h>
#include <iostream>
#include <bestJets_commonNtp.h>

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
    std::cout << "USAGE: ./forMjjStudies_radion_commonNtp [redntpVersion] [dataset] [selectionType] [bTaggerType=JP or CSV]" <<std::endl;
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

  if( dataset=="DiPhoton_madgraph" ) {
    finalize_oneDataset_commonNtp(redntpVersion, "DiPhotonJets_8TeV-madgraph", selectionType,  bTaggerType, datasets);
    
  } else if( dataset=="Radion_M-300_madgraph" ) {  
    finalize_oneDataset_commonNtp(redntpVersion, "Radion_M-300_madgraph", selectionType, bTaggerType, datasets);    
    
  } else if( dataset=="Radion_M-500_madgraph" ) {  
    finalize_oneDataset_commonNtp(redntpVersion, "Radion_M-500_madgraph", selectionType, bTaggerType, datasets);    

  } else if( dataset=="Radion_M-700_madgraph" ) {  
    finalize_oneDataset_commonNtp(redntpVersion, "Radion_M-700_madgraph", selectionType, bTaggerType, datasets);    

  } else if( dataset=="Radion_M-1000_madgraph" ) {  
    finalize_oneDataset_commonNtp(redntpVersion, "Radion_M-1000_madgraph", selectionType, bTaggerType, datasets);    

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

  bestJets_commonNtp* rf = new bestJets_commonNtp( dataset, selectionType, bTaggerType );

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
  if( dataset!="HToGG_M-125_8TeV-pythia6" ) //keep also separate higgs processes to quantify VH                                                                          
    system(rm_command.c_str());
}

std::string getSingleFileName_commonNtp( const std::string& redntpVersion, const std::string& dataset, const std::string& selectionType, const std::string& bTaggerType ) {

  std::string fileName = "finalizedTrees_" + redntpVersion + "/Radion_" + dataset + "_" + selectionType + "_" + bTaggerType + ".root";
  return fileName;
}

RedntpDirStruct get_dirs_commonNtp( const std::string& prodVersion ) {

  RedntpDirStruct returnStruct;

  if( prodVersion=="Radion" ) {

    returnStruct.maindatadir = "/xrootdfs/cms/local/crovelli/Radion/common";  
    returnStruct.mainmcdir   = "/xrootdfs/cms/local/crovelli/Radion/common";
    returnStruct.datadir     = "radion_tree_v03";
    returnStruct.mcdir       = "radion_tree_v03";
  } else {

    std::cout << "-> Unknown prodVersion: '" << prodVersion << "'! Exiting." << std::endl;
    exit(11);
  }

  return returnStruct;
}
