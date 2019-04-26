#!/bin/bash

source /afs/desy.de/user/s/skawada/sonas_work/stage_nnh250/analysis/init_ilcsoft.sh
cd /afs/desy.de/user/s/skawada/sonas_work/stage_nnh250/analysis/output

echo "combining 2 files"
hadd alldata_select.root alldata1_select.root alldata2_select.root

echo "delete intermediate files"
rm alldata1_select.root alldata2_select.root

echo "proc-ing..."
root -l -b <<EOF
.L proc.C+
process()
.q
EOF

echo "copy to LEFT and RIGHT"
cp alldata_select.root LEFT/.
cp alldata_select.root RIGHT/.

echo "finished!"
