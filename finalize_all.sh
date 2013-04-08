#     usage: source finalize_all.sh [redntpProdVersion] [selectionType] [bTaggerType=\"CSV\"]"
# eg: usage: source finalize_all Radion default CSV

 cp finalize2012_radion finalize2012_radion_tmp
 ./finalize2012_radion_tmp  $1 DATA_Run2012_FULL $2 $3
 ./finalize2012_radion_tmp  $1 Radion_M-300_madgraph $2 $3    
 ./finalize2012_radion_tmp  $1 HToGG_M-125_8TeV-pythia6 $2 $3
 ./finalize2012_radion_tmp  $1 DiPhoton_8TeV-pythia6 $2 $3
 ./finalize2012_radion_tmp  $1 V_8TeV $2 $3
 ./finalize2012_radion_tmp  $1 GJet_doubleEMEnriched_TuneZ2star_8TeV-pythia6 $2 $3
 ./finalize2012_radion_tmp  $1 QCD_doubleEMEnriched_TuneZ2star_8TeV-pythia6 $2 $3
 ./finalize2012_radion_tmp  $1 TT_8TeV $2 $3
 ./finalize2012_radion_tmp  $1 VV_8TeV $2 $3
 ./finalize2012_radion_tmp  $1 VGG_8TeV $2 $3
 rm finalize2012_radion_tmp
