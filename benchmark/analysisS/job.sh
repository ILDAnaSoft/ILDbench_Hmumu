#!bin/sh

rm -rf tmp/ *~

./runMarlin hmumu /afs/desy.de/user/s/skawada/sonas_work/bench_nnh500/analysisS/init_ilcsoft.sh hmumu.xml list_hmumu.txt 1 4
chmod 755 condor_submit.sh
mv condor_submit.sh submit-hmumu.sh
echo "preparation done for hmumu."

./runMarlin higgs_ffh /afs/desy.de/user/s/skawada/sonas_work/bench_nnh500/analysisS/init_ilcsoft.sh hmumu.xml list_higgs.txt 1 60
chmod 755 condor_submit.sh
mv condor_submit.sh submit-higgs_ffh.sh
echo "preparation done for higgs_ffh."

./runMarlin 2f_Z_bhabhag /afs/desy.de/user/s/skawada/sonas_work/bench_nnh500/analysisS/init_ilcsoft.sh hmumu.xml list_2f_Z_bhabhag.txt 1 40
chmod 755 condor_submit.sh
mv condor_submit.sh submit-2f_Z_bhabhag.sh
echo "preparation done for 2f_Z_bhabhag."

./runMarlin 2f_Z_leptonic /afs/desy.de/user/s/skawada/sonas_work/bench_nnh500/analysisS/init_ilcsoft.sh hmumu.xml list_2f_Z_leptonic.txt 1 24
chmod 755 condor_submit.sh
mv condor_submit.sh submit-2f_Z_leptonic.sh
echo "preparation done for 2f_Z_leptonic."

./runMarlin 4f_singleW_leptonic /afs/desy.de/user/s/skawada/sonas_work/bench_nnh500/analysisS/init_ilcsoft.sh hmumu.xml list_4f_singleW_leptonic.txt 1 56
chmod 755 condor_submit.sh
mv condor_submit.sh submit-4f_singleW_leptonic.sh
echo "preparation done for 4f_singleW_leptonic."

./runMarlin 4f_singleW_semileptonic /afs/desy.de/user/s/skawada/sonas_work/bench_nnh500/analysisS/init_ilcsoft.sh hmumu.xml list_4f_singleW_semileptonic.txt 1 39
chmod 755 condor_submit.sh
mv condor_submit.sh submit-4f_singleW_semileptonic.sh
echo "preparation done for 4f_singleW_semileptonic."

./runMarlin 4f_singleZee_leptonic /afs/desy.de/user/s/skawada/sonas_work/bench_nnh500/analysisS/init_ilcsoft.sh hmumu.xml list_4f_singleZee_leptonic.txt 1 283
chmod 755 condor_submit.sh
mv condor_submit.sh submit-4f_singleZee_leptonic.sh
echo "preparation done for 4f_singleZee_leptonic."

./runMarlin 4f_singleZee_semileptonic /afs/desy.de/user/s/skawada/sonas_work/bench_nnh500/analysisS/init_ilcsoft.sh hmumu.xml list_4f_singleZee_semileptonic.txt 1 36
chmod 755 condor_submit.sh
mv condor_submit.sh submit-4f_singleZee_semileptonic.sh
echo "preparation done for 4f_singleZee_semileptonic."

./runMarlin 4f_singleZnunu_leptonic /afs/desy.de/user/s/skawada/sonas_work/bench_nnh500/analysisS/init_ilcsoft.sh hmumu.xml list_4f_singleZnunu_leptonic.txt 1 5
chmod 755 condor_submit.sh
mv condor_submit.sh submit-4f_singleZnunu_leptonic.sh
echo "preparation done for 4f_singleZnunu_leptonic."

./runMarlin 4f_singleZnunu_semileptonic /afs/desy.de/user/s/skawada/sonas_work/bench_nnh500/analysisS/init_ilcsoft.sh hmumu.xml list_4f_singleZnunu_semileptonic.txt 1 11
chmod 755 condor_submit.sh
mv condor_submit.sh submit-4f_singleZnunu_semileptonic.sh
echo "preparation done for 4f_singleZnunu_semileptonic."

./runMarlin 4f_singleZsingleWMix_leptonic /afs/desy.de/user/s/skawada/sonas_work/bench_nnh500/analysisS/init_ilcsoft.sh hmumu.xml list_4f_singleZsingleWMix_leptonic.txt 1 36
chmod 755 condor_submit.sh
mv condor_submit.sh submit-4f_singleZsingleWMix_leptonic.sh
echo "preparation done for 4f_singleZsingleWMix_leptonic."

./runMarlin 4f_WW_leptonic /afs/desy.de/user/s/skawada/sonas_work/bench_nnh500/analysisS/init_ilcsoft.sh hmumu.xml list_4f_WW_leptonic.txt 1 7
chmod 755 condor_submit.sh
mv condor_submit.sh submit-4f_WW_leptonic.sh
echo "preparation done for 4f_WW_leptonic."

./runMarlin 4f_WW_semileptonic /afs/desy.de/user/s/skawada/sonas_work/bench_nnh500/analysisS/init_ilcsoft.sh hmumu.xml list_4f_WW_semileptonic.txt 1 43
chmod 755 condor_submit.sh
mv condor_submit.sh submit-4f_WW_semileptonic.sh
echo "preparation done for 4f_WW_semileptonic."

./runMarlin 4f_ZZ_leptonic /afs/desy.de/user/s/skawada/sonas_work/bench_nnh500/analysisS/init_ilcsoft.sh hmumu.xml list_4f_ZZ_leptonic.txt 1 2
chmod 755 condor_submit.sh
mv condor_submit.sh submit-4f_ZZ_leptonic.sh
echo "preparation done for 4f_ZZ_leptonic."

./runMarlin 4f_ZZ_semileptonic /afs/desy.de/user/s/skawada/sonas_work/bench_nnh500/analysisS/init_ilcsoft.sh hmumu.xml list_4f_ZZ_semileptonic.txt 1 41
chmod 755 condor_submit.sh
mv condor_submit.sh submit-4f_ZZ_semileptonic.sh
echo "preparation done for 4f_ZZ_semileptonic."

./runMarlin 4f_ZZWWMix_leptonic /afs/desy.de/user/s/skawada/sonas_work/bench_nnh500/analysisS/init_ilcsoft.sh hmumu.xml list_4f_ZZWWMix_leptonic.txt 1 15
chmod 755 condor_submit.sh
mv condor_submit.sh submit-4f_ZZWWMix_leptonic.sh
echo "preparation done for 4f_ZZWWMix_leptonic."

./runMarlin aa_4f /afs/desy.de/user/s/skawada/sonas_work/bench_nnh500/analysisS/init_ilcsoft.sh hmumu.xml list_aa_4f.txt 1 83
chmod 755 condor_submit.sh
mv condor_submit.sh submit-aa_4f.sh
echo "preparation done for aa_4f."

echo "All preparation done."
