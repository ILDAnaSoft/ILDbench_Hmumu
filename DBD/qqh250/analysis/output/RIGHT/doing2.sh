#!/bin/bash

source /afs/desy.de/user/s/skawada/sonas_work/stage_qqh250/analysis/init_ilcsoft.sh
cd /afs/desy.de/user/s/skawada/sonas_work/stage_qqh250/analysis/output/RIGHT

root -l -b <<EOF
.L driverwithcut.C+
driver("drivertest.dat","sig_drivertest.dat","alldata_wo.root","dataTree",100,0,100,"ntrks")
.q
EOF
