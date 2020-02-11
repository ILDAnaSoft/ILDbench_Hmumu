#!bin/bash

for process in {hmumu,higgs_ffh,2f_Z_bhabhag,2f_Z_hadronic,2f_Z_leptonic,4f_singleW_leptonic,4f_singleW_semileptonic,4f_singleZee_leptonic,4f_singleZee_semileptonic,4f_singleZnunu_leptonic,4f_singleZnunu_semileptonic,4f_singleZsingleWMix_leptonic,4f_WW_hadronic,4f_WW_leptonic,4f_WW_semileptonic,4f_ZZ_hadronic,4f_ZZ_leptonic,4f_ZZ_semileptonic,4f_ZZWWMix_hadronic,4f_ZZWWMix_leptonic,5f,aa_4f,6f_eeWW,6f_llWW,6f_vvWW,6f_xxWW,6f_xxxxZ,6f_yyyyZ,6f_ttbar}
do
    echo "check for ${process}"
    ls ./tmp/${process}/*jobstatus.log | wc -l
    grep "END OF PROCESSING" ./tmp/${process}/*output.log | wc -l
done
