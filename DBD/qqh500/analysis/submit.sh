#!/bin/bash

echo "submission for hmumu"
./submit-hmumu.sh
echo "done."

sleep 1
echo "Next run in 5 seconds."
sleep 5

echo "submission for higgs_ffh"
./submit-higgs_ffh.sh
echo "done."

sleep 1
echo "Next run in 20 seconds."
sleep 20

echo "submission for 2f_Z_bhabhag"
./submit-2f_Z_bhabhag.sh
echo "done."

sleep 1
echo "Next run in 20 seconds."
sleep 20

echo "submission for 2f_Z_hadronic"
./submit-2f_Z_hadronic.sh
echo "done."

sleep 1
echo "Next run in 60 seconds."
sleep 20
echo "Next run in 40 seconds."
sleep 20
echo "Next run in 20 seconds."
sleep 20

echo "submission for 2f_Z_leptonic"
./submit-2f_Z_leptonic.sh
echo "done."

sleep 1
echo "Next run in 20 seconds."
sleep 20

echo "submission for 4f_singleW_leptonic"
./submit-4f_singleW_leptonic.sh
echo "done."

sleep 1
echo "Next run in 20 seconds."
sleep 20

echo "submission for 4f_singleW_semileptonic"
./submit-4f_singleW_semileptonic.sh
echo "done."

sleep 1
echo "Next run in 20 seconds."
sleep 20

echo "submission for 4f_singleZee_leptonic"
./submit-4f_singleZee_leptonic.sh
echo "done."

sleep 1
echo "Next run in 60 seconds."
sleep 20
echo "Next run in 40 seconds."
sleep 20
echo "Next run in 20 seconds."
sleep 20

echo "submission for 4f_singleZee_semileptonic"
./submit-4f_singleZee_semileptonic.sh
echo "done."

sleep 1
echo "Next run in 20 seconds."
sleep 20

echo "submission for 4f_singleZnunu_leptonic"
./submit-4f_singleZnunu_leptonic.sh
echo "done."

sleep 1
echo "Next run in 10 seconds."
sleep 10

echo "submission for 4f_singleZnunu_semileptonic"
./submit-4f_singleZnunu_semileptonic.sh
echo "done."

sleep 1
echo "Next run in 10 seconds."
sleep 10

echo "submission for 4f_singleZsingleWMix_leptonic"
./submit-4f_singleZsingleWMix_leptonic.sh
echo "done."

sleep 1
echo "Next run in 20 seconds."
sleep 20

echo "submission for 4f_WW_leptonic"
./submit-4f_WW_leptonic.sh
echo "done."

sleep 1
echo "Next run in 5 seconds."
sleep 5

echo "submission for 4f_WW_semileptonic"
./submit-4f_WW_semileptonic.sh
echo "done."

sleep 1
echo "Next run in 20 seconds."
sleep 20

echo "submission for 4f_WW_hadronic"
./submit-4f_WW_hadronic.sh
echo "done."

sleep 1
echo "Next run in 20 seconds."
sleep 20

echo "submission for 4f_ZZ_leptonic"
./submit-4f_ZZ_leptonic.sh
echo "done."

sleep 1
echo "Next run in 5 seconds."
sleep 5

echo "submission for 4f_ZZ_semileptonic"
./submit-4f_ZZ_semileptonic.sh
echo "done."

sleep 1
echo "Next run in 10 seconds."
sleep 10

echo "submission for 4f_ZZ_hadronic"
./submit-4f_ZZ_hadronic.sh
echo "done."

sleep 1
echo "Next run in 5 seconds."
sleep 5

echo "submission for 4f_ZZWWMix_leptonic"
./submit-4f_ZZWWMix_leptonic.sh
echo "done."

sleep 1
echo "Next run in 10 seconds."
sleep 10

echo "submission for 5f"
./submit-5f.sh
echo "done."

sleep 1
echo "Next run in 60 seconds."
sleep 20
echo "Next run in 40 seconds."
sleep 20
echo "Next run in 20 seconds."
sleep 20

echo "submission for aa_4f"
./submit-aa_4f.sh
echo "done."

echo "everything done."