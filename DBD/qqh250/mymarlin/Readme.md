# Caution
These analysis processors will only work up to ilcsoft v01-17-09 or around (I did everything with v01-17-09).
For higher version maybe does not work, because the FORTRAN-related things are not supported anymore.
This is due to `SatoruJetFinder` which uses some FORTRAN code.
If you want to use other jet clustering package, modify it by yourself.

# How To Use
0. Assume you are in the plain terminal. Initialize your ilcsoft.
```
source /YOUR_ILCSOFT_PATH/init_ilcsoft.sh
```
1. Do following things.  
```
mkdir build
cd build
cmake -C $ILCSOFT/ILCSoft.cmake ..
make install
cd ../
export MARLIN_DLL=./lib/libmymarlin.so
```  
2. Edit xml file (check your input file(s) (`LCIOInputFiles`) are correct), then run it using next command:
```
Marlin hmumu.xml
```
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
These analysis processors are used for qqh250 analysis.
Basically you only need to edit following files:
```
src/HiggsToMuMuProcessor.cc
include/HiggsToMuMuProcessor.h
src/ISRFinder.cc
include/ISRFinder.h
hmumu.xml
```

## Each file
- `include/HiggsToMuMuProcessor.h`, `include/ISRFinder.h`: defining necessary variables
- `src/ISRFinder.cc`: very simple ISR finder
- `src/HiggsToMuMuProcessor.cc`: extracting/calculating necessary variables, and store in NTuple
- `hmumu.xml`: control input files, processors which you want to use

## Analysis flow
1. `Add4MomCovMatrixCharged`: calculate convariance matrix in momenta space from track parameter space
2. `ThrustReconstruction`, `Sphericity`: calculate event-shape variables
3. `IsolatedLeptonTagging`: to select h->mu+mu- candidate from `PandoraPFOs`
4. `ISRFinder`: to remove ISR photons from `PFOsWithoutIsoleps`
5. `IsolatedLeptonTagging`: count number of isolated leptons in `PFOsWithoutISR`, and use it for vetoing in further analysis
6. `SatoruJetFinder`: reconstruct 2-jets from `PFOsWithoutISR` using Durham algorithm
7. `HiggsToMuMuProcessor`: extracte/calculate necessary variables, and store in NTuple
