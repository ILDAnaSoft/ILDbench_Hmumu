#!/bin/bash

cd /afs/desy.de/user/s/skawada/sonas_work/stage_qqh500/analysis/output

rm *.root *~ LEFT/*.root RIGHT/*.root
var=0
if [ $var -eq 1 ]; then
rm LEFT/*.log LEFT/BDTG/*.root LEFT/BDTG/test/*.root LEFT/BDTG/test/Analysis_result.dat LEFT/BDTG/test/*/* LEFT/BDTG/125fix5/*.root
rm RIGHT/*.log RIGHT/BDTG/*.root RIGHT/BDTG/test/*.root RIGHT/BDTG/test/Analysis_result.dat RIGHT/BDTG/test/*/* RIGHT/BDTG/125fix5/*.root
fi