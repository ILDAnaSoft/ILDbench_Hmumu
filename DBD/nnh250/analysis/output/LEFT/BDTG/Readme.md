# How To Use
0. Assume you have a root file after the preselection named as alldata_precuts.root.
1. Do `root -l -b -q treeseparate.C`.
2. Do `root -l -b -q count.C`, and store these numbers in somewhere.
3. Check the distribution of variable using `simple.C` or even `stackwithcut.C`. Be careful that if you are looking the events in datatest tree, then all event weights should be multiplied by 2 because of the splitting in procedure 1. This factor 2 is already considered and hard-coded in `simple.C` and `stackwithcut.C`.
4. Edit `HiggsToMuMutrain.C`, especially the part of `factory->AddVariable("SOMETHING","F")`.
5. Do `root -l -b -q HiggsToMuMutrain.C`. Usually it takes a few minutes if your root file is small enough.
6. Do following things. Then the control panel will be appeared.
```
root -l
.L TMVAGui.C
TMVAGui("HiggsToMuMu.root")
```
7. Go to `test/` directory. This is for detemining TMVA internal parameters. Details can be found there.  
8. Return to this directory. Now you should have your input variables and TMVA internal parameters. Adjust them in `HiggsToMuMutrain.C`, and do `root -l -b -q HiggsToMuMutrain.C`.  
9. Do following things.
```
cp alldata_precuts.root alldata_after.root
root -l
.L result.C
result("alldata_after.root")
.q
```
Now you have 2 root files: one is just the result of precuts and other is also contain the result of TMVA analysis.

# Explanation
This place is used fot TMVA analysis, especially for determination of input variables to TMVA.
Also final TMVA analysis is performed at here after the study at `test/` directory.
In all cases, I only used BDTG technique.
If you want to use other TMVA techniques, you need to modify proper places to perform TMVA analysis.

## Each file
- `treeseparate.C`: divide 1 data tree into 2 data trees. In my case, dataTree is divided into datatrain and datatest.
- `count.C`: count number of MC events and weighted events for 2 trees
- `HiggsToMuMutrain.C`: controling for TMVA method, input variables, internal TMVA parameters like nCuts, MaxDepth, and so on
- `result.C`: to add TMVA result to your root file. Be careful that the order of input variables should be the same as you set at `HiggsToMuMutrain.C`, otherwise it does not work properly.
