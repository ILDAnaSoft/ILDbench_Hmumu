#!/bin/bash

source /afs/desy.de/user/s/skawada/sonas_work/stage_qqh250/analysis/init_ilcsoft.sh
cd /afs/desy.de/user/s/skawada/sonas_work/stage_qqh250/analysis/output

echo "combining 1f_3f background"
hadd 1f_3f.root 1f_3f/1f_3f-???.root

echo "combine aa_2f background"
hadd aa_2f.root aa_2f/aa_2f-???.root

echo "combine 2 processes"
hadd alldata2.root 1f_3f.root aa_2f.root

echo "delete intermediate files"
rm 1f_3f.root aa_2f.root

echo "finished!"
