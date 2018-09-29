#!/bin/bash

for momres in {1_3,5_4,3_4,2_4,1_4,5_5,3_5,2_5,1_5,5_6,3_6,2_6,1_6}
do
    rm ./${momres}/data2.txt ./${momres}/truedata.txt
    paste -d " " ./${momres}/data.txt ./scan/data.txt > ./${momres}/data2.txt

    while read line
    do
	a=($line)
	echo "${a[0]} ${a[1]} ${a[2]} ${a[3]} ${a[4]} ${a[5]} ${a[6]} ${a[7]} ${a[16]}" >> ./${momres}/truedata.txt
    done < ./${momres}/data2.txt
done
