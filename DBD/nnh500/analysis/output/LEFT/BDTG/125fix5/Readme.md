# Caution
This part contains DESY-specific part.
Sorry for ugly name.

# How To Use
0. Assume you have the root file named alldata_after.root.
1. Do `doing.sh > doing.output`.
2. Check `doing.output`, like chi2/ndf value, and remove strange cases like too high chi2/ndf. Also need to remove strange cases in `data.txt`.
3. Do `. autoscan.sh` (DESY-specific), or create your own job submission script. If you set higher number of pseudo-experiments, it will take very long time. In my case (50000 times, hard-coded), typically it takes 6-8 hours. In maximum (my case), this will submit 78 jobs, thus you need 468-624 CPU hours typically.

You can do similar things to smeared cases.
1. Create your necessary directories like `1_3` and `1_3_fix`.
2. Do `momres_doing.sh > momres_doing.output`.
3. Check `momres_doing.output`, and remove strange cases from each data.
4. Do `add.sh`.
5. Do `. automomres.sh` (DESY-specific), or create your own job submission script. In maximum (my case), this will submit 1014 jobs, thus you need 6084-8112 CPU hours typically.

# Explanation
This is the final part of analysis. This place is used for toy MC to extract final precision.

## Each file
- `param_tmp.C`: do modeling for muon pair invariant mass, Crystal Ball + Gaussian for signal, pol1 for background
- `momres_param_tmp.C`: similar as `param_tmp.C`, for smeared momentum resolution case
- `doing.sh`: apply BDTG score cut, do modeling, and combine all fitting data into single file as `data.txt`
- `momres_doing.sh`: similar as `doing.sh`, for smeared momentum resolution case
- `scan/`: both yield of signal and background are free parameters
- `tmp/`: similar as `scan/`, for smeared momentum resolution case
- `scan_fix/`: only yield of signal is free parameter
- `tmp_fix/`: similar as `scan_fix/`, for smeared momentum resolution case 
