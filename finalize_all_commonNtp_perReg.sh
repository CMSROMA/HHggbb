#     usage: source finalize_all_commonNtp_preReg.sh [redntpProdVersion] [selectionType] [bTaggerType=\"CSV\"]"
# eg: usage: source finalize_all_commonNtp_perReg.sh Radion_V05 default CSV
            
 cp finalize2012_radion_commonNtp finalize2012_radion_commonNtp_tmp

 ./finalize2012_radion_commonNtp_tmp  $1 DataABCD_regr $2 $3

 ./finalize2012_radion_commonNtp_tmp  $1 Radion_M-300_regr $2 $3    
 ./finalize2012_radion_commonNtp_tmp  $1 Radion_M-500_regr $2 $3    
 ./finalize2012_radion_commonNtp_tmp  $1 Radion_M-700_regr $2 $3    
 ./finalize2012_radion_commonNtp_tmp  $1 Radion_M-1000_regr $2 $3    

 rm finalize2012_radion_commonNtp_tmp
