#!/bin/bash

for file in $(find /pnfs/desy.de/ilc/prod/ilc/mc-opt-3/ild/dst-merged/500-TDR_ws -name "*.slcio")
do
    timeout 60 lcio_event_counter $file > /dev/null
    if [ $? == 124 ] ; then
	echo "fucking timed out file: $file"
    fi
done
