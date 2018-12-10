# Caution
These codes will only work with ilcsoft v02-00-02.
For other versions maybe does not work, or tons of warning message will appear.

The `VertexInfo` functionality has some memory-related problem which I don't know how to solve it.
At the end of the job, `Marlin` will give us the crush information, while expected numbers are correctly filled in NTuple.
For analyzers point of view it is not a problem, but this problem should be fixed at some point.


# How To Use
0. Assume you are in the plain terminal. Initialize your ilcsoft by doing `source /YOUR_ILCSOFT_PATH/init_ilcsoft.sh`.
1. Do following things.  
```
mkdir build
cd build  
cmake -C $ILCSOFT/ILCSoft.cmake ..  
make install  
cd ../  
export MARLIN_DLL=./lib/libmymarlin.so
```  
2. Edit xml file (check your input file(s) (`LCIOInputFiles`) are correct), then run it using `Marlin hmumu.xml`
3. Now you should have output.root.

Alternatively (and my favorite style), you can copy-and-paste init_ilcsoft.sh from somewhere, and add next sentence at the end of `init_ilcsoft.sh`.  
```
export MARLIN_DLL=/YOUR_WORKING_DIRECTORY/mymarlin/lib/libmymarlin.so
```
Then, when you restart your work, just type following command to initialize your ilcsoft, together with including your own library.
```
source /YOUR_WORKING_DIRECTORY/mymarlin/init_ilcsoft.sh
```

# Explanation
These codes are used for nnh500 analysis of IDR benchmark.
Basically you only need to edit following files:
```
src/HiggsToMuMuProcessor.cc  
include/HiggsToMuMuProcessor.h  
hmumu.xml
```

## Each file
- `include/HiggsToMuMuProcessor.h`: defining necessary variables
- `src/HiggsToMuMuProcessor.cc`: extracting/calculating necessary variables, and store in NTuple
- `VertexInfo`: extracting primary vertex position from 2 tracks; in this analysis 2 tracks correspond to 2 muons which is selected in `IsolatedLeptonTagging`
- `hmumu.xml`: control input files, processors which you want to use

## Analysis flow
0. `InitDD4hep`: necessary for `VertexInfo` functionality
1. `Add4MomCovMatrixCharged`: calculate convariance matrix in momenta space from track parameter space
2. `ThrustReconstruction`,`Sphere`: calculate event-shape variables
3. `IsolatedLeptonTagging`: to select h->mu+mu- candidate from `PandoraPFOs`
4. `HiggsToMuMuProcessor`: extract/calculate necessary variables, and store in NTuple
