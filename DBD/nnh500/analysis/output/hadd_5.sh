#!/bin/bash

source /afs/desy.de/user/s/skawada/sonas_work/stage_nnh500/analysis/init_ilcsoft.sh
cd /afs/desy.de/user/s/skawada/sonas_work/stage_nnh500/analysis/output

echo "combining part of 6f_ttbar"
hadd 6f_ttbar5.root 6f_ttbar/6f_ttbar-[2-4]??.root

echo "combining higgs_ffhf"
hadd higgs_ffh.root higgs_ffh/*.root

echo "combine everything"
hadd alldata5.root 6f_ttbar5.root higgs_ffh.root

echo "delete intermediate files"
rm 6f_ttbar5.root higgs_ffh.root

echo "finished!"
