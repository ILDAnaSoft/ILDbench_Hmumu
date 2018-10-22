#!/bin/bash

source /afs/desy.de/user/s/skawada/sonas_work/bench_nnh500/analysis/init_ilcsoft.sh
cd /afs/desy.de/user/s/skawada/sonas_work/bench_nnh500/analysis/output/RIGHT

rm *~

echo "start add weight"
sleep 5
root -l -b alldata_select.root <<EOF
.L AddWeight.C+
AddWeight("alldata_select.root",hev)
.q
EOF
echo "finished"

sleep 5

echo "start remove h->mumu events from general qqh/llh/nnh samples"
sleep 5
root -l -b <<EOF
.L SkimCut.C+
SkimCut("alldata_select.root","alldata_wo.root","dataTree","!((processid>=106515&&processid<=106526)&&abs(higgsdecay1)==13)")
.q
EOF
echo "finished"

sleep 5

echo "start select events with # mu+- = 1"
sleep 5
root -l -b <<EOF
.L SkimCut.C+
SkimCut("alldata_wo.root","alldata_cut1.root","dataTree","n_muminus==1&&n_muplus==1")
.q
EOF

echo "all finished"
