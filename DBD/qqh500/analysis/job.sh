#!bin/sh

rm -rf tmp/ *~
#./runMarlin alllist init_ilcsoft.sh nnh_250GeV.xml list_all_250GeV.txt 1 1506

./runMarlin hmumu /afs/desy.de/user/s/skawada/sonas_work/stage_qqh500/analysis/init_ilcsoft.sh hmumu.xml list_hmumu_full.txt 1 4
chmod 755 condor_submit.sh
mv condor_submit.sh submit-hmumu.sh
echo "preparation done for hmumu."

./runMarlin higgs_ffh /afs/desy.de/user/s/skawada/sonas_work/stage_qqh500/analysis/init_ilcsoft.sh hmumu.xml list_higgs.txt 1 34
chmod 755 condor_submit.sh
mv condor_submit.sh submit-higgs_ffh.sh
echo "preparation done for higgs_ffh."

./runMarlin 2f_Z_bhabhag /afs/desy.de/user/s/skawada/sonas_work/stage_qqh500/analysis/init_ilcsoft.sh hmumu.xml list_2f_Z_bhabhag.txt 1 28
chmod 755 condor_submit.sh
mv condor_submit.sh submit-2f_Z_bhabhag.sh
echo "preparation done for 2f_Z_bhabhag."

./runMarlin 2f_Z_hadronic /afs/desy.de/user/s/skawada/sonas_work/stage_qqh500/analysis/init_ilcsoft.sh hmumu.xml list_2f_Z_hadronic.txt 1 120
chmod 755 condor_submit.sh
mv condor_submit.sh submit-2f_Z_hadronic.sh
echo "preparation done for 2f_Z_hadronic."

./runMarlin 2f_Z_leptonic /afs/desy.de/user/s/skawada/sonas_work/stage_qqh500/analysis/init_ilcsoft.sh hmumu.xml list_2f_Z_leptonic.txt 1 26
chmod 755 condor_submit.sh
mv condor_submit.sh submit-2f_Z_leptonic.sh
echo "preparation done for 2f_Z_leptonic."

./runMarlin 4f_singleW_leptonic /afs/desy.de/user/s/skawada/sonas_work/stage_qqh500/analysis/init_ilcsoft.sh hmumu.xml list_4f_singleW_leptonic.txt 1 34
chmod 755 condor_submit.sh
mv condor_submit.sh submit-4f_singleW_leptonic.sh
echo "preparation done for 4f_singleW_leptonic."

./runMarlin 4f_singleW_semileptonic /afs/desy.de/user/s/skawada/sonas_work/stage_qqh500/analysis/init_ilcsoft.sh hmumu.xml list_4f_singleW_semileptonic.txt 1 28
chmod 755 condor_submit.sh
mv condor_submit.sh submit-4f_singleW_semileptonic.sh
echo "preparation done for 4f_singleW_semileptonic."

./runMarlin 4f_singleZee_leptonic /afs/desy.de/user/s/skawada/sonas_work/stage_qqh500/analysis/init_ilcsoft.sh hmumu.xml list_4f_singleZee_leptonic.txt 1 129
chmod 755 condor_submit.sh
mv condor_submit.sh submit-4f_singleZee_leptonic.sh
echo "preparation done for 4f_singleZee_leptonic."

./runMarlin 4f_singleZee_semileptonic /afs/desy.de/user/s/skawada/sonas_work/stage_qqh500/analysis/init_ilcsoft.sh hmumu.xml list_4f_singleZee_semileptonic.txt 1 20
chmod 755 condor_submit.sh
mv condor_submit.sh submit-4f_singleZee_semileptonic.sh
echo "preparation done for 4f_singleZee_semileptonic."

./runMarlin 4f_singleZnunu_leptonic /afs/desy.de/user/s/skawada/sonas_work/stage_qqh500/analysis/init_ilcsoft.sh hmumu.xml list_4f_singleZnunu_leptonic.txt 1 8
chmod 755 condor_submit.sh
mv condor_submit.sh submit-4f_singleZnunu_leptonic.sh
echo "preparation done for 4f_singleZnunu_leptonic."

./runMarlin 4f_singleZnunu_semileptonic /afs/desy.de/user/s/skawada/sonas_work/stage_qqh500/analysis/init_ilcsoft.sh hmumu.xml list_4f_singleZnunu_semileptonic.txt 1 6
chmod 755 condor_submit.sh
mv condor_submit.sh submit-4f_singleZnunu_semileptonic.sh
echo "preparation done for 4f_singleZnunu_semileptonic."

./runMarlin 4f_singleZsingleWMix_leptonic /afs/desy.de/user/s/skawada/sonas_work/stage_qqh500/analysis/init_ilcsoft.sh hmumu.xml list_4f_singleZsingleWMix_leptonic.txt 1 20
chmod 755 condor_submit.sh
mv condor_submit.sh submit-4f_singleZsingleWMix_leptonic.sh
echo "preparation done for 4f_singleZsingleWMix_leptonic."

./runMarlin 4f_WW_leptonic /afs/desy.de/user/s/skawada/sonas_work/stage_qqh500/analysis/init_ilcsoft.sh hmumu.xml list_4f_WW_leptonic.txt 1 5
chmod 755 condor_submit.sh
mv condor_submit.sh submit-4f_WW_leptonic.sh
echo "preparation done for 4f_WW_leptonic."

./runMarlin 4f_WW_semileptonic /afs/desy.de/user/s/skawada/sonas_work/stage_qqh500/analysis/init_ilcsoft.sh hmumu.xml list_4f_WW_semileptonic.txt 1 26
chmod 755 condor_submit.sh
mv condor_submit.sh submit-4f_WW_semileptonic.sh
echo "preparation done for 4f_WW_semileptonic."

./runMarlin 4f_WW_hadronic /afs/desy.de/user/s/skawada/sonas_work/stage_qqh500/analysis/init_ilcsoft.sh hmumu.xml list_4f_WW_hadronic.txt 1 21
chmod 755 condor_submit.sh
mv condor_submit.sh submit-4f_WW_hadronic.sh
echo "preparation done for 4f_WW_hadronic."

./runMarlin 4f_ZZ_leptonic /afs/desy.de/user/s/skawada/sonas_work/stage_qqh500/analysis/init_ilcsoft.sh hmumu.xml list_4f_ZZ_leptonic.txt 1 6
chmod 755 condor_submit.sh
mv condor_submit.sh submit-4f_ZZ_leptonic.sh
echo "preparation done for 4f_ZZ_leptonic."

./runMarlin 4f_ZZ_semileptonic /afs/desy.de/user/s/skawada/sonas_work/stage_qqh500/analysis/init_ilcsoft.sh hmumu.xml list_4f_ZZ_semileptonic.txt 1 18
chmod 755 condor_submit.sh
mv condor_submit.sh submit-4f_ZZ_semileptonic.sh
echo "preparation done for 4f_ZZ_semileptonic."

./runMarlin 4f_ZZ_hadronic /afs/desy.de/user/s/skawada/sonas_work/stage_qqh500/analysis/init_ilcsoft.sh hmumu.xml list_4f_ZZ_hadronic.txt 1 7
chmod 755 condor_submit.sh
mv condor_submit.sh submit-4f_ZZ_hadronic.sh
echo "preparation done for 4f_ZZ_hadronic."

./runMarlin 4f_ZZWWMix_leptonic /afs/desy.de/user/s/skawada/sonas_work/stage_qqh500/analysis/init_ilcsoft.sh hmumu.xml list_4f_ZZWWMix_leptonic.txt 1 13
chmod 755 condor_submit.sh
mv condor_submit.sh submit-4f_ZZWWMix_leptonic.sh
echo "preparation done for 4f_ZZWWMix_leptonic."

./runMarlin 4f_ZZWWMix_hadronic /afs/desy.de/user/s/skawada/sonas_work/stage_qqh500/analysis/init_ilcsoft.sh hmumu.xml list_4f_ZZWWMix_hadronic.txt 1 18
chmod 755 condor_submit.sh
mv condor_submit.sh submit-4f_ZZWWMix_hadronic.sh
echo "preparation done for 4f_ZZWWMix_hadronic."

./runMarlin 5f /afs/desy.de/user/s/skawada/sonas_work/stage_qqh500/analysis/init_ilcsoft.sh hmumu.xml list_5f.txt 1 400
chmod 755 condor_submit.sh
mv condor_submit.sh submit-5f.sh
echo "preparation done for 5f."

./runMarlin aa_4f /afs/desy.de/user/s/skawada/sonas_work/stage_qqh500/analysis/init_ilcsoft.sh hmumu.xml list_aa_4f.txt 1 160
chmod 755 condor_submit.sh
mv condor_submit.sh submit-aa_4f.sh
echo "preparation done for aa_4f."

./runMarlin 6f_eeWW /afs/desy.de/user/s/skawada/sonas_work/stage_qqh500/analysis/init_ilcsoft.sh hmumu.xml list_6f_eeWW.txt 1 48
chmod 755 condor_submit.sh
mv condor_submit.sh submit-6f_eeWW.sh
echo "preparation done for 6f_eeWW."

./runMarlin 6f_llWW /afs/desy.de/user/s/skawada/sonas_work/stage_qqh500/analysis/init_ilcsoft.sh hmumu.xml list_6f_llWW.txt 1 24
chmod 755 condor_submit.sh
mv condor_submit.sh submit-6f_llWW.sh
echo "preparation done for 6f_llWW."

./runMarlin 6f_vvWW /afs/desy.de/user/s/skawada/sonas_work/stage_qqh500/analysis/init_ilcsoft.sh hmumu.xml list_6f_vvWW.txt 1 32
chmod 755 condor_submit.sh
mv condor_submit.sh submit-6f_vvWW.sh
echo "preparation done for 6f_vvWW."

./runMarlin 6f_xxWW /afs/desy.de/user/s/skawada/sonas_work/stage_qqh500/analysis/init_ilcsoft.sh hmumu.xml list_6f_xxWW.txt 1 35
chmod 755 condor_submit.sh
mv condor_submit.sh submit-6f_xxWW.sh
echo "preparation done for 6f_xxWW."

./runMarlin 6f_xxxxZ /afs/desy.de/user/s/skawada/sonas_work/stage_qqh500/analysis/init_ilcsoft.sh hmumu.xml list_6f_xxxxZ.txt 1 14
chmod 755 condor_submit.sh
mv condor_submit.sh submit-6f_xxxxZ.sh
echo "preparation done for 6f_xxxxZ."

./runMarlin 6f_yyyyZ /afs/desy.de/user/s/skawada/sonas_work/stage_qqh500/analysis/init_ilcsoft.sh hmumu.xml list_6f_yyyyZ.txt 1 40
chmod 755 condor_submit.sh
mv condor_submit.sh submit-6f_yyyyZ.sh
echo "preparation done for 6f_yyyyZ."

./runMarlin 6f_ttbar /afs/desy.de/user/s/skawada/sonas_work/stage_qqh500/analysis/init_ilcsoft.sh hmumu.xml list_6f_ttbar.txt 1 1550
chmod 755 condor_submit.sh
mv condor_submit.sh submit-6f_ttbar.sh
echo "preparation done for 6f_ttbar."

echo "All preparation done."
