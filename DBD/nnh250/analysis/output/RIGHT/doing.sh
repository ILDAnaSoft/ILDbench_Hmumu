#!/bin/bash

source /afs/desy.de/user/s/skawada/sonas_work/stage_nnh250/analysis/init_ilcsoft.sh
cd /afs/desy.de/user/s/skawada/sonas_work/stage_nnh250/analysis/output/RIGHT

rm *~ *.log

root -l -b alldata_select.root <<EOF
.L AddWeight.C+
MakePolInfo("procid250.txt")
AddWeight("alldata_select.root",hev)
.q
EOF

sleep 5

root -l -b <<EOF
.L SkimCut.C+
SkimCut("alldata_select.root","alldata_wo.root","dataTree","!((processid>=106475&&processid<=106486)&&abs(higgsdecay1)==13)")
.q
EOF

sleep 5

root -l -b <<EOF
.L SkimCut.C+
SkimCut("alldata_wo.root","alldata_cut1.root","dataTree","n_muminus==1&&n_muplus==1")
.q
EOF
