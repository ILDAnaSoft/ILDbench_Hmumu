#!/bin/bash

#nbins=200
nbinsbkg=50
for cut in {-0.95,-0.9,-0.85,-0.8,-0.75,-0.7,-0.65,-0.6,-0.55,-0.5,-0.45,-0.4,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95}
do
    echo "BDTGcut is ${cut}"
    cp param_tmp.C param.C
    sed -i "s/BDTGBDTG/${cut}/g" param.C
    sed -i "s/NBINSBKGNBINSBKG/${nbinsbkg}/g" param.C
    root -l -b <<EOF
.L param.C
param()
.q
EOF
    echo "BDTGcut ${cut} finished"
sleep 1
done

echo "combining..."
cat output_-0.95.txt output_-0.9.txt output_-0.85.txt output_-0.8.txt output_-0.75.txt output_-0.7.txt output_-0.65.txt output_-0.6.txt output_-0.55.txt output_-0.5.txt output_-0.45.txt output_-0.4.txt output_-0.35.txt output_-0.3.txt output_-0.25.txt output_-0.2.txt output_-0.15.txt output_-0.1.txt output_-0.05.txt output_0.txt output_0.05.txt output_0.1.txt output_0.15.txt output_0.2.txt output_0.25.txt output_0.3.txt output_0.35.txt output_0.4.txt output_0.45.txt output_0.5.txt output_0.55.txt output_0.6.txt output_0.65.txt output_0.7.txt output_0.75.txt output_0.8.txt output_0.85.txt output_0.9.txt output_0.95.txt > data.txt

sleep 1

rm output_-0.95.txt output_-0.9.txt output_-0.85.txt output_-0.8.txt output_-0.75.txt output_-0.7.txt output_-0.65.txt output_-0.6.txt output_-0.55.txt output_-0.5.txt output_-0.45.txt output_-0.4.txt output_-0.35.txt output_-0.3.txt output_-0.25.txt output_-0.2.txt output_-0.15.txt output_-0.1.txt output_-0.05.txt output_0.txt output_0.05.txt output_0.1.txt output_0.15.txt output_0.2.txt output_0.25.txt output_0.3.txt output_0.35.txt output_0.4.txt output_0.45.txt output_0.5.txt output_0.55.txt output_0.6.txt output_0.65.txt output_0.7.txt output_0.75.txt output_0.8.txt output_0.85.txt output_0.9.txt output_0.95.txt
echo "All finished!"
