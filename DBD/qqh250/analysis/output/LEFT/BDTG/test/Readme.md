# Caution
This part contains DESY-specific part.

# Hou To Use
0. Assume you have the file named alldata_precuts.root, also this file has 2 trees named datatrain and datatest.
1. Make following directories: `C/`, `dat/`, `exec/`, `log/`, `root/`, and `scripts/`.
2. Edit `HiggsToMuMutrain.C.orig` and `result.C`. Make sure that the order of input variables should be the same between these 2 files, otherwise it does not work. You can test it will work by using `. test_job.sh`.
3. `condor_submit job.submit` (DESY-specific), or create your own job submission script because it takes typically 8-12 hours with in total 720 combination cases.
4. Do `Analysis_job.sh`.
5. Do `. submit_all.sh` (DESY-specific), or create your own job submission script because the total number of jobs will be 720 when you use these scripts as it was (for each job typically takes a few ten minutes).
6. Do `cat dat/* > Analysis_result.dat`.
7. Do `root -l -b -q Maximum.C`. Finally you will obtain the best significance and the combination of TMVA internal variables.

# Explanation
This place is used for TMVA analysis, especially to determine the TMVA internal variables.

## Each file
- `job.sh`: setup the TMVA internal variable, train and test, then change the combination of internal variable, then train and test. In my analysis, I used typical numbers for these TMVA internal variables, and iys total is 720 combinations. You can modify tochange the combination variables, but of course be careful when you increase the number of combinations, it takes even much longer time.
- `Analysis_job.sh`: prepare scripts to evaluate the significance
- `Maximum.C`: pickup best significance case, and show corresponding internal variables
