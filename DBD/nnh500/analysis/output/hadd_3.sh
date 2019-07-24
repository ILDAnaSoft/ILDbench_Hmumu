#!/bin/bash

source /afs/desy.de/user/s/skawada/sonas_work/stage_nnh500/analysis/init_ilcsoft.sh
cd /afs/desy.de/user/s/skawada/sonas_work/stage_nnh500/analysis/output

echo "combining part of 6f_ttbar"
hadd 6f_ttbar3.root 6f_ttbar/6f_ttbar-1[2-3]??.root

echo "combining 6f_vvWW, xxWW"
hadd 6f_2.root 6f_vvWW/*.root 6f_xxWW/*.root

echo "combine everything"
hadd alldata3.root 6f_ttbar3.root 6f_2.root

echo "delete intermediate files"
rm 6f_ttbar3.root 6f_2.root

echo "finished!"
