# Caution
This part contains DESY-specific part.
Sorry for ugly name.

# How To Use
0. Assume you have the root file named alldata_after.root.
1. Do `doing.sh > doing.output`.
2. Check `doing.output`, like chi2/ndf value, and remove strange cases like too high chi2/ndf. Also need to remove strange cases in `data.txt`.
3. Do `. autoscan.sh` (DESY-specific), or create your own job submission script. If you set higher number of pseudo-experiments, it will take very long time. In my case (50000 times), typically need 6-8 hours.

# Explanation
This is the final part of analysis. This place is used for toy MC to extract final precision.

## Each file
- `param_tmp.C`: do modeling for muon pair invariant mass, Crystal Ball + Gaussian for signal, pol1 for background
- 'doing.sh': apply BDTG score cut, do modeling, and combine all fitting data into single file as `data.txt`
- `scan/`: both yield of signal and background are free parameters
- `scan_fix/`: only yield of signal is free parameter
