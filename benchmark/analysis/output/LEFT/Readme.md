# Caution
This part contains DESY-specific part.

# How To Use
1. preparation
- `condor_submit doing.submit` (DESY-specific)
- Or do `. doing.sh`(Not recommend to do this in worker node at usual working time because it takes long time. Do it in lunchtime/midnight or create your own job submission script and submit it is much better for your colleague.)

2. cut-based analysis
```
root -l alldata_cut1.root
.L stackwithcut.C+
MakeAllWithWeight("file1","file2",dataTree,nbin,nlow,nhigh,"variable","cut_condition")
```
Then you should see the histogram of variable with cut_condition applied (if cut_condition == 1, then no cuts are applied).
Now you can investigate your favorite cuts.

3. create your precuts-applied root file
```
root -l -b
.L SkimCut.C+
SkimCut("alldata_cut1.root","alldata_precuts.root","dataTree","YOUR_PRESELECTION_CUTS")
.q
```
When you finish your cut-based analysis for preselection, it is recommended to create your own preselected file with using `SkimCut.C`.
The final product alldata_precuts.root is used for further analysis.

# Explanation
This place is used for further preparation and cut-based analysis.
Mainly this part is used to define preselection before TMVA analysis.

## Each file
1. preparation
- `doing.sh`: 3 steps will be performed
  - Add event weight to each event, this weight is calculated to adjust 1.6 ab-1 with left-handed beam polarization. This step uses `AddWeight.C`.
  - Remove h->mu+mu- event in ffh samples to avoid confusion, because we have dedicated ffh_mumu samples. This uses alldata_select.root, and produces alldata_wo.root.
  - Select events which only have one mu+ and one mu-. This uses alldata_wo.root, and produces alldata_cut1.root.
2. cut-based analysis
- `stackwithcut.C`: create histogram of variable, using different color for different processes, hard-coded many things to specify process
- `simple.C`: similar to `stackwithcut.C` but all backgrounds are treated inclusively
3. tools
- `norm.C`: create histogram of variable, using different color for different processes, hard-coded many things to specify process, all histograms are normalized to 1
- `norm_sig.C`: similar to `norm.C`, but only show normalized signal histogram, useful for comparison between IDR-L and IDR-S
