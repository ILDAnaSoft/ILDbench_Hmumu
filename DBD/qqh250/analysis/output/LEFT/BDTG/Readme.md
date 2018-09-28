# Caution
This part contains DESY-specific part.

# How To Use
0. Assume you have a root file after the preselection named as alldata_precuts.root.
1. Do `root -l -b -q treeseparate.C`.
2. Do `root -l -b -q count.C`, and store these numbers in somewhere.


# Explanation
This place is used fot TMVA analysis.


## Each file
- `treeseparate.C`: divide 1 data tree into 2 data trees. In my case, dataTree is divided into datatrain and datatest.
- `count.C`: count number of MC events and weighted events for 2 trees
