### ILDbench_Hmumu
Main analyzer: Shin-ichi Kawada (DESY)
shin-ichi.kawada@desy.de
(wrote on 2018/Sep./24)

### Introduction
These codes are used for h->mu+mu- branching ratio measurement at the ILD.
Since the structure of directory looks very different compare to other usual usage case,
I will try to explain how to use these.
Use these codes are your own risk.

### Structure
You can see 2 different directories; DBD/ and benchmark/ .
You can easily imagine that the files in DBD/ directory are used for DBD-style samples.
Several different types of codes are used, so I will explain it at there.

The files in bechmark/ directory are used for IDR analysis.






### Installation

Explain here:

- what are the package dependencies (iLCSoft, others ?)
- how to compile your package. Should normally be something like:

```shell
source /path/to/ilcsoft/init_ilcsoft.sh
mkdir build
cd build
cmake -C $ILCSOFT/ILCSoft.cmake ..
make install
```

### How to run the analysis

Explain here:

- where to find data needed for your analysis or how to produce them
- how to run you analysis: 
   - Marlin processors to run ?
   - ROOT macros to run ?
   - Shell scripts ?
   - Run the analysis on grid if you provide scripts for that

Example:

```shell
export MARLIN_DLL=./lib/libILDbench_Hmumu.so
Marlin ./scripts/ExampleProcessor.xml
```

If you want to provide a lot of details on your analysis, use the doc/Readme.md and point to it from this Readme.md file:

More documentation available here in [doc/Readme.md](doc/Readme.md) !

