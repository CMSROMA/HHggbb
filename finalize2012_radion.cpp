#include <TMath.h>
#include <iostream>
#include <fillPlot2012_radion.h>

struct RedntpDirStruct {

  std::string maindatadir;
  std::string mainmcdir;
  std::string datadir;
  std::string mcdir;
};

void finalize_oneDataset(const std::string& redntpVersion,  const std::string& dataset, const std::string& selectionType, const std::string& bTaggerType, std::vector<std::string> *datasets );

void do_haddCommand( const std::string& redntpVersion, const std::string& dataset, std::vector<std::string> *datasets, const std::string& selectionType, const std::string& bTaggerType );

std::string getSingleFileName( const std::string& redntpVersion, const std::string& dataset, const std::string& selectionType, const std::string& bTaggerType );

RedntpDirStruct get_dirs( const std::string& prodVersion );

int main( int argc, char* argv[] ) {

  if( argc!=4 && argc!=5 ) {
    std::cout << "USAGE: ./finalize2012_radion [redntpVersion] [dataset] [selectionType] [bTaggerType=JP or CSV]" <<std::endl;
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

  if( dataset=="DATA_Run2012_FULL" ) {  
    finalize_oneDataset(redntpVersion, "full2012_minusAug24_minusDec11_minusD", selectionType, bTaggerType, datasets); 

  } else if( dataset=="DiPhotonBox_8TeV-pythia6" ) {
    finalize_oneDataset(redntpVersion, "DiPhotonBox_Pt-10To25_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1",   selectionType, bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "DiPhotonBox_Pt-250ToInf_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1", selectionType, bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "DiPhotonBox_Pt-25To250_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1",  selectionType, bTaggerType, datasets);

  } else if( dataset=="DiPhoton_8TeV-pythia6" ) {

    finalize_oneDataset(redntpVersion, "DiPhotonBox_Pt-10To25_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1",    selectionType,  bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "DiPhotonBox_Pt-250ToInf_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1",  selectionType,  bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "DiPhotonBox_Pt-25To250_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1",   selectionType,  bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "DiPhotonJets_8TeV-madgraph-tarball-v2_Summer12_DR53X-PU_S10_START53_V7A-v1", selectionType,  bTaggerType, datasets);

  } else if( dataset=="V_8TeV") {
                                         
    finalize_oneDataset(redntpVersion, "DYJetsToLL_M-10To50filter_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1",  selectionType, bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1",	selectionType, bTaggerType, datasets); 
    finalize_oneDataset(redntpVersion, "WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1",      selectionType, bTaggerType, datasets);

  } else if( dataset=="GJet_doubleEMEnriched_TuneZ2star_8TeV-pythia6" ) {

    finalize_oneDataset(redntpVersion, "GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1", selectionType, bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "GJet_Pt40_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1",      selectionType, bTaggerType, datasets);

  } else if( dataset=="QCD_doubleEMEnriched_TuneZ2star_8TeV-pythia6") {

    finalize_oneDataset(redntpVersion, "QCD_Pt-40_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1",     selectionType, bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "QCD_Pt-30to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1", selectionType, bTaggerType, datasets);

  } else if( dataset=="TT_8TeV" ) {
                                        
    finalize_oneDataset(redntpVersion, "TTbarGG_0Jet_S1" , selectionType, bTaggerType, datasets); //chiara:52vx3
    finalize_oneDataset(redntpVersion, "TTGJets_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1", selectionType, bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "TTJets_TuneZ2star_8TeV-madgraph-tauola_Summer12-PU_S7_START52_V9-v1", selectionType, bTaggerType, datasets); //chiara:52xv3

  } else if( dataset=="VV_8TeV" ) {

    finalize_oneDataset(redntpVersion, "WW_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1", selectionType, bTaggerType, datasets); 
    finalize_oneDataset(redntpVersion, "ZZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1", selectionType, bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "WZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1", selectionType, bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "WGToLNuG_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1", selectionType, bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "ZGToLLG_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1",             selectionType, bTaggerType, datasets);

  } else if( dataset=="VGG_8TeV" ) {
    
    finalize_oneDataset(redntpVersion, "WmGG", selectionType, bTaggerType, datasets);  //chiara:52xv3     
    finalize_oneDataset(redntpVersion, "WpGG", selectionType, bTaggerType, datasets);  //chiara:52xv3     
    finalize_oneDataset(redntpVersion, "ZGG",  selectionType, bTaggerType, datasets);  //chiara:52xv3     

  } else if( dataset=="HToGG_M-125_8TeV-pythia6" ) {

    finalize_oneDataset(redntpVersion, "WH_ZH_HToGG_M-125_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1",   selectionType, bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "GluGluToHToGG_M-125_8TeV-powheg-pythia6_Summer12-PU_S7_START52_V9-v1",  selectionType, bTaggerType, datasets);  //chiara:52xv3
    finalize_oneDataset(redntpVersion, "VBF_HToGG_M-125_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1", selectionType, bTaggerType, datasets);  //chiara:53xv2
    finalize_oneDataset(redntpVersion, "TTH_HToGG_M-125_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1",     selectionType, bTaggerType, datasets);

  } else if( dataset=="Radion_M-300_madgraph" ) {  

    finalize_oneDataset(redntpVersion, "Radion_M-300_madgraph", selectionType, bTaggerType, datasets);    

  } else if( dataset=="Radion_M-500_madgraph" ) {  

    finalize_oneDataset(redntpVersion, "Radion_M-500_madgraph", selectionType, bTaggerType, datasets);    

  } else if( dataset=="Radion_M-700_madgraph" ) {  

    finalize_oneDataset(redntpVersion, "Radion_M-700_madgraph", selectionType, bTaggerType, datasets);    

  } else if( dataset=="Radion_M-1000_madgraph" ) {  

    finalize_oneDataset(redntpVersion, "Radion_M-1000_madgraph", selectionType, bTaggerType, datasets);    

  }
  
  do_haddCommand(redntpVersion, dataset, datasets, selectionType, bTaggerType );

  return 0;
}

void finalize_oneDataset( const std::string& redntpProdVersion, const std::string& dataset, const std::string& selectionType, const std::string& bTaggerType, std::vector<std::string> *datasets ) {

  TString dataset_tstr(dataset);
  
  std::cout << std::endl << std::endl << std::endl << "####     Finalizing " << dataset << std::endl;
  std::cout << "####  Selection: " << selectionType << std::endl;
  std::cout << "####  b-Tagger: "  << bTaggerType << std::endl << std::endl;

  RedntpDirStruct dirs = get_dirs( redntpProdVersion );

  fillPlot2012_radion* rf = new fillPlot2012_radion( dataset, selectionType, bTaggerType );

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


void do_haddCommand( const std::string& redntpVersion, const std::string& dataset, std::vector<std::string> *datasets, const std::string& selectionType, const std::string& bTaggerType ) {

  if( datasets->size()<=1 ) return;

  std::string hadd_command = "hadd -f " + getSingleFileName( redntpVersion, dataset, selectionType, bTaggerType );

  std::string suffix = "";
  for( unsigned i=0; i<datasets->size(); ++i )
    suffix = suffix + " " + getSingleFileName( redntpVersion, datasets->at(i), selectionType, bTaggerType );

  hadd_command += suffix;
  std::string rm_command = "rm " + suffix;
  system(hadd_command.c_str());
  if( dataset!="HToGG_M-125_8TeV-pythia6" ) //keep also separate higgs processes to quantify VH                                                                          
    system(rm_command.c_str());
}

std::string getSingleFileName( const std::string& redntpVersion, const std::string& dataset, const std::string& selectionType, const std::string& bTaggerType ) {

  std::string fileName = "finalizedTrees_" + redntpVersion + "/Radion_" + dataset + "_" + selectionType + "_" + bTaggerType + ".root";
  return fileName;
}

RedntpDirStruct get_dirs( const std::string& prodVersion ) {

  RedntpDirStruct returnStruct;

  if( prodVersion=="Radion" ) {

    returnStruct.maindatadir = "/xrootdfs/cms/local/crovelli/Radion/reduced";  
    returnStruct.mainmcdir   = "/xrootdfs/cms/local/crovelli/Radion/reduced";
    returnStruct.datadir     = "redntp.53xv2_data.cicpfloose.scales-Lisbon-Hgg.moriond_dataset/merged";  
    returnStruct.mcdir       = "redntp.52xv3_53xv2_53xv3.cicpfloose.scales-Lisbon-Hgg.radion_v3/merged";

  } else {

    std::cout << "-> Unknown prodVersion: '" << prodVersion << "'! Exiting." << std::endl;
    exit(11);
  }

  return returnStruct;
}
