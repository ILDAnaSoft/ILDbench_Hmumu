#!/bin/bash

echo "START hadd..."

echo "combining signal process"
hadd hmumu.root hmumu/hmumu-???.root

echo "combining ffh background"
hadd higgs_ffh.root higgs_ffh/higgs_ffh-???.root

echo "combining 2f backgrounds"
hadd 2f_Z_bhabhag.root 2f_Z_bhabhag/2f_Z_bhabhag-???.root
hadd 2f_Z_hadronic.root 2f_Z_hadronic/2f_Z_hadronic-???.root
hadd 2f_Z_leptonic.root 2f_Z_leptonic/2f_Z_leptonic-???.root
echo "combine into single 2f file"
hadd 2f_all.root 2f_Z_bhabhag.root 2f_Z_hadronic.root 2f_Z_leptonic.root
echo "deleting intermediate files"
rm 2f_Z_bhabhag.root 2f_Z_hadronic.root 2f_Z_leptonic.root

echo "combining 4f backgrounds"
hadd 4f_singleW_leptonic.root 4f_singleW_leptonic/4f_singleW_leptonic-???.root
hadd 4f_singleW_semileptonic.root 4f_singleW_semileptonic/4f_singleW_semileptonic-???.root
hadd 4f_singleZee_leptonic.root 4f_singleZee_leptonic/4f_singleZee_leptonic-???.root
hadd 4f_singleZee_semileptonic.root 4f_singleZee_semileptonic/4f_singleZee_semileptonic-???.root
hadd 4f_singleZnunu_leptonic.root 4f_singleZnunu_leptonic/4f_singleZnunu_leptonic-???.root
hadd 4f_singleZnunu_semileptonic.root 4f_singleZnunu_semileptonic/4f_singleZnunu_semileptonic-???.root
hadd 4f_singleZsingleWMix_leptonic.root 4f_singleZsingleWMix_leptonic/4f_singleZsingleWMix_leptonic-???.root
hadd 4f_WW_leptonic.root 4f_WW_leptonic/4f_WW_leptonic-???.root
hadd 4f_WW_semileptonic.root 4f_WW_semileptonic/4f_WW_semileptonic-???.root
hadd 4f_WW_hadronic.root 4f_WW_hadronic/4f_WW_hadronic-???.root
hadd 4f_ZZ_leptonic.root 4f_ZZ_leptonic/4f_ZZ_leptonic-???.root
hadd 4f_ZZ_semileptonic.root 4f_ZZ_semileptonic/4f_ZZ_semileptonic-???.root
hadd 4f_ZZ_hadronic.root 4f_ZZ_hadronic/4f_ZZ_hadronic-???.root
hadd 4f_ZZWWMix_leptonic.root 4f_ZZWWMix_leptonic/4f_ZZWWMix_leptonic-???.root
hadd 4f_ZZWWMix_hadronic.root 4f_ZZWWMix_hadronic/4f_ZZWWMix_hadronic-???.root
echo "combine into single 4f file and delete intermediate files"
hadd 4f_all.root 4f_singleW_leptonic.root 4f_singleW_semileptonic.root 4f_singleZee_leptonic.root 4f_singleZee_semileptonic.root 4f_singleZnunu_leptonic.root 4f_singleZnunu_semileptonic.root 4f_singleZsingleWMix_leptonic.root 4f_WW_leptonic.root 4f_WW_semileptonic.root 4f_WW_hadronic.root 4f_ZZ_leptonic.root 4f_ZZ_semileptonic.root 4f_ZZ_hadronic.root 4f_ZZWWMix_leptonic.root 4f_ZZWWMix_hadronic.root
echo "delete intermediate files"
rm 4f_singleW_leptonic.root 4f_singleW_semileptonic.root 4f_singleZee_leptonic.root 4f_singleZee_semileptonic.root 4f_singleZnunu_leptonic.root 4f_singleZnunu_semileptonic.root 4f_singleZsingleWMix_leptonic.root 4f_WW_leptonic.root 4f_WW_semileptonic.root 4f_WW_hadronic.root 4f_ZZ_leptonic.root 4f_ZZ_semileptonic.root 4f_ZZ_hadronic.root 4f_ZZWWMix_leptonic.root 4f_ZZWWMix_hadronic.root

echo "combining files except aa_2f and 1f_3f"
hadd alldata1.root hmumu.root higgs_ffh.root 2f_all.root 4f_all.root

echo "delete intermediate files"
rm hmumu.root higgs_ffh.root 2f_all.root 4f_all.root

echo "combining 1f_3f background"
hadd 1f_3f.root 1f_3f/1f_3f-???.root

echo "combine aa_2f background"
hadd aa_2f.root aa_2f/aa_2f-???.root

echo "combine 2 processes"
hadd alldata2.root 1f_3f.root aa_2f.root

echo "delete intermediate files"
rm 1f_3f.root aa_2f.root

echo "finished!"
