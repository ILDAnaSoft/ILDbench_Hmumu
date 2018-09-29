#!/bin/bash

source /afs/desy.de/user/s/skawada/sonas_work/stage_nnh500/analysis/init_ilcsoft.sh
cd /afs/desy.de/user/s/skawada/sonas_work/stage_nnh500/analysis/output/RIGHT

rm *~

root -l -b alldata_select.root <<EOF
.L AddWeight.C+
MakePolInfo("procid500.txt")
AddWeight("alldata_select.root",hev)
.q
EOF

sleep 5

root -l -b <<EOF
.L SkimCut.C+
SkimCut("alldata_select.root","alldata_wo.root","dataTree","!((processid>=106515&&processid<=106526)&&abs(higgsdecay1)==13)")
.q
EOF

sleep 5

root -l -b <<EOF
.L SkimCut.C+
SkimCut("alldata_wo.root","alldata_cut1.root","dataTree","n_muminus==1&&n_muplus==1")
.q
EOF
