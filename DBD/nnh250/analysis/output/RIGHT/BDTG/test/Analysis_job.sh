#!/bin/bash

HERE="/afs/desy.de/user/s/skawada/sonas_work/stage_nnh250/analysis/output/RIGHT/BDTG/test"
ALL="submit_all.sh"
echo "#!/bin/bash" > $ALL

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
                    filename="Result${cut}_${shrinkage}_${depth}_${tree}_${node}.root"
                    output="Analysis${cut}_${shrinkage}_${depth}_${tree}_${node}.dat"
		    script="Analysis${cut}_${shrinkage}_${depth}_${tree}_${node}.C"
		    submit="submit${cut}_${shrinkage}_${depth}_${tree}_${node}.submit"
		    exec="exec${cut}_${shrinkage}_${depth}_${tree}_${node}.sh"
		    cp Analysis.C.orig $script
		    sed -i "s/AAAAA/$filename/g" $script
		    sed -i "s/BBBBB/$output/g" $script

		    echo "condor_submit ${HERE}/${submit}" >> $ALL

		    echo "executable = ${exec}" > $submit
		    echo "universe = vanilla" >> $submit
		    echo "notify_user = shin-ichi.kawada@desy.de" >> $submit
		    echo "JobNotification = Error" >> $submit
                    echo "output = ${HERE}/log/output${cut}_${shrinkage}_${depth}_${tree}_${node}.log" >> $submit
                    echo "log = ${HERE}/log/log${cut}_${shrinkage}_${depth}_${tree}_${node}.log" >> $submit
                    echo "error = ${HERE}/log/error${cut}_${shrinkage}_${depth}_${tree}_${node}.log" >> $submit
                    echo "+RequestRuntime = 3600*3" >> $submit
                    echo "RequestMemory = 2G" >> $submit
                    echo "RequestDisk = 2G" >> $submit
                    echo "Requirements = OpSysAndVer == \"SL6\"" >> $submit
		    echo "queue" >> $submit

		    echo "source /afs/desy.de/user/s/skawada/sonas_work/stage_nnh250/analysis/init_ilcsoft.sh" >> $exec
		    echo "cd ${HERE}" >> $exec
		    echo "root -l -b <<EOF" >> $exec
		    echo ".L $script" >> $exec
		    echo "Analysis()" >> $exec
		    echo "EOF" >> $exec
		    echo "echo \"$cut $shrinkage $depth $tree $node\" >> $output" >> $exec
		    echo "mv ${filename} root/." >> $exec
		    echo "mv ${output} dat/." >> $exec
		    echo "mv ${script} C/." >> $exec
		    echo "mv ${submit} script/." >> $exec
		    echo "mv ${exec} exec/." >> $exec

		    chmod 755 $exec
		done
	    done
	done
    done
done