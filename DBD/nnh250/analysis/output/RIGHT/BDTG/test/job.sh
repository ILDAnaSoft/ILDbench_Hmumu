#!/bin/bash

echo "Job started:" > time.txt
date >> time.txt

source /afs/desy.de/user/s/skawada/sonas_work/stage_nnh250/analysis/init_ilcsoft.sh

#for nCuts AAAAA
#for Shrinkage BBBBB
#for MaxDepth CCCCC
#for NTrees DDDDD
#for MinNodeSize EEEEE%
for cut in {10,20,30}
do
    for shrinkage in {0.1,0.15,0.2,0.25,0.3}
    do
	for depth in {2,3,4,5}
	do
	    for tree in {300,500,800,1000}
	    do
		for node in {3,5,10}
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
	    done
	done
    done
echo "cut ${cut} ended:" >> time.txt
date >> time.txt
done

echo "Job ended:" >> time.txt
date >> time.txt
