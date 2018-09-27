# How To Use
1. preparation
- condor_submit doing.submit (DESY-specific)
- Or do `. doing.sh`(Not recommend to do this in worker node at usual working time because it takes long time. Do it in lunchtime/midnight or create your own job submission script and submit it is much better for your colleague.)
2. cut-based analysis
```
root -l alldata_cut1.root
.L stackwithcut.C+
MakeProcinfo("procid250.txt")
MakeAllWithWeight("file1","file2",dataTree,nbin,nlow,nhigh,"variable","cut_condition")
```
Then you should see the histogram of variable with cut_condition applied (if cut_condition == 1, then no cuts are applied).
Now you can investigate your favorite cuts.
```
root -l
.L driverwithcut.C+
driver("file1","file2","inputrootfile","dataTree",nbin,nlow,nhigh,"variable")
```
The `driverwithcut.C` will provide the result when you apply the cut sequentially.

# Explanation
This place is used for further preparation and cut-based analysis.

## Each file
1. preparation
- `doing.sh`: 3 steps will be performed
  - Add event weight to each event, this weight is calculated to adjust 0.9 ab-1 with left-handed beam polarization. This step uses `AddWeight.C`.
  - Remove h->mu+mu- event in DBD samples to avoid confusion, because we have dedicated ffh_mumu samples. This uses alldata_select.root, and produces alldata_wo.root.
  - Select events which only have one mu+ and one mu-. This used alldata_wo.root, and produces alldata_cut1.root.
2. cut-based analysis
- `stackwithcut.C`: create histogram of variable, using different color for different processes, hard-coded many things to specify process
- `simple.C`: similar to `stackwithcut.C` but all backgrounds are treated inclusively
- `driverwithcut.C`: can use for cut optimization, create entire cut table