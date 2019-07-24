#!/bin/bash

source /afs/desy.de/user/s/skawada/sonas_work/stage_nnh500/analysis/init_ilcsoft.sh
cd /afs/desy.de/user/s/skawada/sonas_work/stage_nnh500/analysis/output

echo "combining part of 6f_ttbar"
hadd 6f_ttbar4.root 6f_ttbar/6f_ttbar-1[4-5]??.root 6f_ttbar/6f_ttbar-8??.root

echo "combining aa_4f"
hadd aa_4f.root aa_4f/*.root

echo "combine everything"
hadd alldata4.root 6f_ttbar4.root aa_4f.root

echo "delete intermediate files"
rm 6f_ttbar4.root aa_4f.root

echo "finished!"
