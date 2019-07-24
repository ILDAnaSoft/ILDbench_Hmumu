#!/bin/bash

source /afs/desy.de/user/s/skawada/sonas_work/stage_qqh500/analysis/init_ilcsoft.sh
cd /afs/desy.de/user/s/skawada/sonas_work/stage_qqh500/analysis/output

echo "combining 2f backgrounds"
hadd 2f_Z_bhabhag.root 2f_Z_bhabhag/2f_Z_bhabhag-???.root
hadd 2f_Z_hadronic.root 2f_Z_hadronic/2f_Z_hadronic-???.root
hadd 2f_Z_leptonic.root 2f_Z_leptonic/2f_Z_leptonic-???.root
echo "combine into single 2f file"
hadd 2f_all.root 2f_Z_bhabhag.root 2f_Z_hadronic.root 2f_Z_leptonic.root
echo "deleting intermediate files"
rm 2f_Z_bhabhag.root 2f_Z_hadronic.root 2f_Z_leptonic.root

echo "combining part of 6f_ttbar"
hadd 6f_ttbar7.root 6f_ttbar/6f_ttbar-[0-1]??.root

echo "combine everything"
hadd alldata7.root 6f_ttbar7.root 2f_all.root

echo "delete intermediate files"
rm 2f_all.root 6f_ttbar7.root

echo "finished!"
