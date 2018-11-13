#!/bin/bash

echo "Job started:" > time.txt
date >> time.txt

#for nCuts AAAAA
#for Shrinkage BBBBB
#for MaxDepth CCCCC
#for NTrees DDDDD
#for MinNodeSize EEEEE%

shrinkage=0.1
depth=3
tree=1000
node=5
for cut in {10,20,30}
do
    cp HiggsToMuMutrain.C.orig HiggsToMuMutrain.C
    sed -i "s/AAAAA/$cut/g" HiggsToMuMutrain.C
    sed -i "s/BBBBB/$shrinkage/g" HiggsToMuMutrain.C
    sed -i "s/CCCCC/$depth/g" HiggsToMuMutrain.C
    sed -i "s/DDDDD/$tree/g" HiggsToMuMutrain.C
    sed -i "s/EEEEE/$node/g" HiggsToMuMutrain.C
    root -l -b -q HiggsToMuMutrain.C
    filename="Result${cut}_${shrinkage}_${depth}_${tree}_${node}.root"
    cp alldata_precuts.root $filename
    root -l -b << EOF
.L result.C
result("$filename")
EOF
done

echo "Job ended:" >> time.txt
date >> time.txt
