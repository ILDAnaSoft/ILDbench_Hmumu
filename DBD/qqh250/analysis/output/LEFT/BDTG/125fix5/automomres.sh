#!/bin/bash

for momres in {1_3,5_4,3_4,2_4,1_4,5_5,3_5,2_5,1_5,5_6,3_6,2_6,1_6}
do
    cd /afs/desy.de/user/s/skawada/sonas_work/stage_qqh250/analysis/output/LEFT/BDTG/125fix5/${momres}
    cp truedata.txt ../${momres}_fix/.
    rm *~ *orig *log exec*submit submit*sh toyMC*C
    echo "prepare scripts for ${momres}..."
    cp ../tmp/auto.sh .
    cp ../tmp/rootlogon.C .
    cp ../tmp/*orig .
    sed -i "s/MOMRESMOMRES/${momres}/g" exec.submit.orig
    sed -i "s/MOMRESMOMRES/${momres}/g" submit.sh.orig
    sed -i "s/MOMRESMOMRES/${momres}/g" toyMC.C.orig
    ./auto.sh truedata.txt
    sleep 2
    echo "start submitting..."
    . submit_all.sh
    echo "submission finished for ${momres}"
    echo "sleeping..."
    sleep 10

    cd ../${momres}_fix
    rm *~ *orig *log exec*submit submit*sh toyMC*C
    echo "prepare scripts for ${momres}_fix..."
    cp ../tmp_fix/auto.sh .
    cp ../tmp_fix/rootlogon.C .
    cp ../tmp_fix/*orig .
    sed -i "s/MOMRESMOMRES/${momres}/g" exec.submit.orig
    sed -i "s/MOMRESMOMRES/${momres}/g" submit.sh.orig
    sed -i "s/MOMRESMOMRES/${momres}/g" toyMC.C.orig
    ./auto.sh truedata.txt
    sleep 2
    echo "start submitting..."
    . submit_all.sh
    echo "submission finished for ${momres}_fix"
    echo "sleeping..."
    sleep 10
done

cd ../
echo "here is..."
pwd
echo "DONE"
