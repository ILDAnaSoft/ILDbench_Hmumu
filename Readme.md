# ILDbench_Hmumu
Main analyzer: Shin-ichi Kawada (DESY)  
shin-ichi.kawada[atmk]desy.de  
(updated on 2018/Oct./1, many things under construction, but for DBD part is finished)

# Introduction
These codes are used for h->mu+mu- branching ratio measurement at the ILD.
Since the structure of directory looks very different compare to other usual usage case, I will try to explain how to use these at each place.
Use these codes are your own risk.

# Structure
You can see 2 different directories; `DBD/` and `benchmark/`.
You can easily imagine that the files in `DBD/` directory are used for DBD-style samples.
Several different types of codes are used, so I will explain it at there.

The files in `bechmark/` directory are used for IDR analysis.
This is WIP.

# Analysis Flow
Basically in all cases, the same analysis flow is used.
In short:
1. initialize your ilcsoft
2. edit codes in `src/` and `include/` directories
3. compile (if error occurs, you can enjoy bugfix)
4. edit steering file, then run it
5. technical things; merge root files, extract necessary variables, add weights, and so on (this is my favorite style, you can do it more clever way)
6. enjoy cut-based analysis before TMVA, this is preselection part
7. enjoy TMVA analysis
8. enjoy toy MC analysis

From procedure 5, I always created directory to perform next analysis. e.g,: the directory structure is getting like this: 5/6/7/8/ .
In each channel, I always created `LEFT/` and `RIGHT/` directories.
LEFT means the beam polarization of P(e+,e-) = (+0.3,-0.8), and RIGHT means P(e+,e-) = (-0.3,+0.8).

Some codes/macros contain DESY-specific and ilcsoft-version-specific contents.
You have to adjust to your own environment.
