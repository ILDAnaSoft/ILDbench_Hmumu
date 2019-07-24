#!/bin/bash

source /afs/desy.de/user/s/skawada/sonas_work/stage_qqh500/analysis/init_ilcsoft.sh
cd /afs/desy.de/user/s/skawada/sonas_work/stage_qqh500/analysis/output

echo "combining all files"
hadd alldata_select.root alldata1_select.root alldata2_select.root alldata3_select.root alldata4_select.root alldata5_select.root alldata6_select.root alldata7_select.root

sleep 5

echo "delete intermediate files"
rm alldata1_select.root alldata2_select.root alldata3_select.root alldata4_select.root alldata5_select.root alldata6_select.root alldata7_select.root

sleep 5

echo "proc-ing..."
root -l -b <<EOF
.L proc.C+
process()
.q
EOF

sleep 5

echo "copy to LEFT and RIGHT"
cp alldata_select.root LEFT/.
cp alldata_select.root RIGHT/.

echo "finished!"