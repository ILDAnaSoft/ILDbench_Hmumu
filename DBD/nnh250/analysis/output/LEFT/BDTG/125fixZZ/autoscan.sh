#!/bin/bash

cp data.txt ./scan/.
cp data.txt ./scan_fix/.

cd scan/
rm *~ *log exec*submit submit*sh toyMC*C result*txt
./auto.sh data.txt
sleep 2
echo "start submitting..."
. submit_all.sh
echo "submission finished for scan/"

cd ../scan_fix/
rm *~ *log exec*submit submit*sh toyMC*C result*txt
./auto.sh data.txt
sleep 2
echo "start submitting..."
. submit_all.sh
echo "submission finished for scan_fix/"

cd ../
echo "DONE"
