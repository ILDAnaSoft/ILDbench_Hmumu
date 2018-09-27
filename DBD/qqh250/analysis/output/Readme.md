# Caution
This part contains many DESY-specific part.
You need to be careful.

# How To Use
- condor_submit doing.submit (DESY-specific)
- Or, do `. doing.sh` (Not recommend to do this in worker node at usual working time because it takes very long time. Do it in midnight or create your own job submission script and submit it is much better for your colleague.)

## Step by step case
```
. hadd.sh
root -l -b
.L SkimVar.C+
SkimVar("alldata1.root","alldata1_select.root","dataTree","var1:var2:var3:...")
SkimVar("alldata2.root","alldata2_select.root","dataTree","var1:var2:var3:...")
.q
. hadd2.sh
root -l -b
.L proc.C+
process()
.q
cp alldata_select.root LEFT/.
cp alldata_select.root RIGHT/.
```

# Explanation
This place is used for prepation of analysis: combine all small root files to single large root file, extract necessary variables, add process ID histogram, and copy the processed root file at `LEFT/.` and `RIGHT/.`.
This is my favorite style, you can do it more clever way I think.

## Each file
- `hadd.sh`,`hadd2.sh`: for combining each root file into single large file (When the total file size is too huge, then ROOT will split this huge file. This is the reason why I have 2 scripts for merging.)
- `SkimVar.C`: for extracting necessary variables from root file, somehow self-documented
- `proc.C`: creating one histogram of process ID, will be used for further prepartion
- `selectvar.txt`: it would be useful to summarize your necessary variables, and easy to use `SkimVar.C`, just copy-and-paste it