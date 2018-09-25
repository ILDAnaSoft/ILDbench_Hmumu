# How To Use
0. Assume you are in the plain terminal. Initialize your ilcsoft.  
'''bash:
source /YOUR_ILCSOFT_PATH/init_ilcsoft.sh
'''  
1. Do following things.  
'''bash:
mkdir build
cd build
cmake -C $ILCSOFT/ILCSoft.cmake ..
make install
cd ../
export MARLIN_DLL=./lib/libmymarlin.so
'''  
2. Edit xml file, then run it using next command:.
'''bash:
Marlin hmumu.xml
'''

Alternatively, you can copy-and-paste init_ilcsoft.sh from somewhere, and add next sentence at the end of file.  
'''bash:
export MARLIN_DLL=/YOUR_WORKING_DIRECTORY/mymarlin/lib/libmymarlin.so
'''

# Explanation

