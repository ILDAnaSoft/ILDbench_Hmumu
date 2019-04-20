#!/bin/bash

echo "START hadd..."

echo "combining signal process"
hadd hmumu.root hmumu/hmumu-???.root

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
echo "combine into single 4f file and delete intermediate files"
hadd 4f_all.root 4f_singleW_leptonic.root 4f_singleW_semileptonic.root 4f_singleZee_leptonic.root 4f_singleZee_semileptonic.root 4f_singleZnunu_leptonic.root 4f_singleZnunu_semileptonic.root 4f_singleZsingleWMix_leptonic.root 4f_WW_leptonic.root 4f_WW_semileptonic.root 4f_WW_hadronic.root 4f_ZZ_leptonic.root 4f_ZZ_semileptonic.root 4f_ZZ_hadronic.root 4f_ZZWWMix_leptonic.root
echo "delete intermediate files"
rm 4f_singleW_leptonic.root 4f_singleW_semileptonic.root 4f_singleZee_leptonic.root 4f_singleZee_semileptonic.root 4f_singleZnunu_leptonic.root 4f_singleZnunu_semileptonic.root 4f_singleZsingleWMix_leptonic.root 4f_WW_leptonic.root 4f_WW_semileptonic.root 4f_WW_hadronic.root 4f_ZZ_leptonic.root 4f_ZZ_semileptonic.root 4f_ZZ_hadronic.root 4f_ZZWWMix_leptonic.root

echo "combining 5f background"
hadd 5f.root 5f/5f-???.root

echo "combining part of 6f_ttbar"
hadd 6f_ttbar1.root 6f_ttbar/6f_ttbar-9??.root

echo "combine everything"
hadd alldata1.root hmumu.root 4f_all.root 5f.root 6f_ttbar1.root

echo "delete intermediate files"
rm hmumu.root higgs_ffh.root 2f_all.root 4f_all.root 5f.root 6f_ttbar1.root

echo "finished!"
