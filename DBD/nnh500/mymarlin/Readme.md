# Caution
These codes will only work up to ilcsoft v01-17-09 or around (I did everything with v01-17-09).
For higher version maybe does not work, or tons of warning message will appear.

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

Alternatively (and my favorite style), you can copy-and-paste init_ilcsoft.sh from somewhere, and add next sentence at the end of init_ilcsoft.sh.  
```
export MARLIN_DLL=/YOUR_WORKING_DIRECTORY/mymarlin/lib/libmymarlin.so
```
Then, when you restart your work, just type following command to initialize your ilcsoft, together with including your own library.
```
source /YOUR_WORKING_DIRECTORY/mymarlin/init_ilcsoft.sh
```

# Explanation
These codes are used for nnh500 analysis.
Basically you only need to edit following files:
```
src/HiggsToMuMuProcessor.cc  
include/HiggsToMuMuProcessor.h  
hmumu.xml
```
Other files are used to calculate evet-shape variables.

## Each file
- `include/HiggsToMuMuProcessor.h`: defining necessary variables
- `src/HiggsToMuMuProcessor.cc`: extracting/calculating necessary variables, and store in NTuple
- `hmumu.xml`: control input files, processors which you want to use
- `procid500.txt`: external file contains process name/ID, beam polarization, and cross section

## Analysis flow
1. `Add4MomCovMatrixCharged`: calculate convariance matrix in momenta space from track parameter space
2. `ThrustReconstruction`, `Sphericity`: calculate event-shape variables
3. `IsolatedLeptonTagging`: to select h->mu+mu- candidate from `PandoraPFOs`
4. `HiggsToMuMuProcessor`: extract/calculate necessary variables, and store in NTuple
