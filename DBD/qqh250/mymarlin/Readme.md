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
2. Edit xml file (check your input (`LCIOInputFiles`) file is correct), then run it using next command:
```
Marlin hmumu.xml
```
3. Now you should have output.root.

Alternatively (and my favorite style), you can copy-and-paste init_ilcsoft.sh from somewhere, and add next sentence at the end of init_ilcsoft.sh.  
```
export MARLIN_DLL=/YOUR_WORKING_DIRECTORY/mymarlin/lib/libmymarlin.so
```
Then, when you restart your work, just do following things to initialize your ilcsoft, together with your own library.
```
source /YOUR_WORKING_DIRECTORY/mymarlin/init_ilcsoft.sh
```

# Explanation

