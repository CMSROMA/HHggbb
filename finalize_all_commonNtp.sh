#     usage: source finalize_all_commonNtp.sh [redntpProdVersion] [selectionType] [bTaggerType=\"CSV\"]"
# eg: usage: source finalize_all_commonNtp.sh Radion default CSV

 cp finalize2012_radion_commonNtp finalize2012_radion_commonNtp_tmp

 ./finalize2012_radion_commonNtp_tmp  $1 Data2012 $2 $3

 #./finalize2012_radion_commonNtp_tmp  $1 Radion_M-300_madgraph $2 $3    
 #./finalize2012_radion_commonNtp_tmp  $1 Radion_M-500_madgraph $2 $3    
 #./finalize2012_radion_commonNtp_tmp  $1 Radion_M-700_madgraph $2 $3    
 #./finalize2012_radion_commonNtp_tmp  $1 Radion_M-1000_madgraph $2 $3    

 #./finalize2012_radion_commonNtp_tmp  $1 HToGG_M-125_8TeV-pythia6 $2 $3
 #./finalize2012_radion_commonNtp_tmp  $1 WH_ZH_HToGG_M-125_8TeV $2 $3
 #./finalize2012_radion_commonNtp_tmp  $1 GluGluToHToGG_M-125_8TeV $2 $3
 #./finalize2012_radion_commonNtp_tmp  $1 VBF_HToGG_M-125_8TeV $2 $3
 #./finalize2012_radion_commonNtp_tmp  $1 TTH_HToGG_M-125_8TeV $2 $3

 #./finalize2012_radion_commonNtp_tmp  $1 DiPhoton_8TeV $2 $3

 #./finalize2012_radion_commonNtp_tmp  $1 DiPhoton_madgraph $2 $3

 ##./finalize2012_radion_commonNtp_tmp  $1 V_8TeV $2 $3
 #./finalize2012_radion_commonNtp_tmp  $1 GJet_doubleEMEnriched_TuneZ2star_8TeV-pythia6 $2 $3
 #./finalize2012_radion_commonNtp_tmp  $1 QCD_doubleEMEnriched_TuneZ2star_8TeV-pythia6 $2 $3
 ##./finalize2012_radion_commonNtp_tmp  $1 TT_8TeV $2 $3
 ##./finalize2012_radion_commonNtp_tmp  $1 VV_8TeV $2 $3
 ##./finalize2012_radion_commonNtp_tmp  $1 VGG_8TeV $2 $3
 rm finalize2012_radion_commonNtp_tmp
