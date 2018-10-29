#!/bin/bash

echo "check Higgs files"
for file in $(find /pnfs/desy.de/ilc/prod/ilc/mc-opt-3/ild/dst-merged/500-TDR_ws/higgs_ffh/ILD_s5_o1_v02/ -name "*.slcio")
do
    timeout 60 lcio_event_counter $file > /dev/null
    if [ $? == 124 ] ; then
	echo "fucking timed out file: $file"
    fi
done

echo "check 2f_Z_bhabhg"
for file in $(find /pnfs/desy.de/ilc/prod/ilc/mc-opt-3/ild/dst-merged/500-TDR_ws/2f_Z_bhabhag/ILD_s5_o1_v02/ -name "*.slcio")
do
    timeout 60 lcio_event_counter $file > /dev/null
    if [ $? == 124 ] ; then
	echo "fucking timed out file: $file"
    fi
done

echo "check 2f_Z_leptonic"
for file in $(find /pnfs/desy.de/ilc/prod/ilc/mc-opt-3/ild/dst-merged/500-TDR_ws/2f_Z_leptonic/ILD_s5_o1_v02/ -name "*.slcio")
do
    timeout 60 lcio_event_counter $file > /dev/null
    if [ $? == 124 ] ; then
	echo "fucking timed out file: $file"
    fi
done

echo "check aa_4f"
for file in $(find /pnfs/desy.de/ilc/prod/ilc/mc-opt-3/ild/dst-merged/500-TDR_ws/aa_4f/ILD_s5_o1_v02/ -name "*.slcio")
do
    timeout 60 lcio_event_counter $file > /dev/null
    if [ $? == 124 ] ; then
	echo "fucking timed out file: $file"
    fi
done

echo "check 4f_singleW_leptonic"
for file in $(find /pnfs/desy.de/ilc/prod/ilc/mc-opt-3/ild/dst-merged/500-TDR_ws/4f_singleW_leptonic/ILD_s5_o1_v02/ -name "*.slcio")
do
    timeout 60 lcio_event_counter $file > /dev/null
    if [ $? == 124 ] ; then
	echo "fucking timed out file: $file"
    fi
done

echo "check 4f_singleW_semileptonic"
for file in $(find /pnfs/desy.de/ilc/prod/ilc/mc-opt-3/ild/dst-merged/500-TDR_ws/4f_singleW_semileptonic/ILD_s5_o1_v02/ -name "*.slcio")
do
    timeout 60 lcio_event_counter $file > /dev/null
    if [ $? == 124 ] ; then
	echo "fucking timed out file: $file"
    fi
done

echo "check 4f_singleZee_leptonic"
for file in $(find /pnfs/desy.de/ilc/prod/ilc/mc-opt-3/ild/dst-merged/500-TDR_ws/4f_singleZee_leptonic/ILD_s5_o1_v02/ -name "*.slcio")
do
    timeout 60 lcio_event_counter $file > /dev/null
    if [ $? == 124 ] ; then
	echo "fucking timed out file: $file"
    fi
done

echo "check 4f_singleZee_semileptonic"
for file in $(find /pnfs/desy.de/ilc/prod/ilc/mc-opt-3/ild/dst-merged/500-TDR_ws/4f_singleZee_semileptonic/ILD_s5_o1_v02/ -name "*.slcio")
do
    timeout 60 lcio_event_counter $file > /dev/null
    if [ $? == 124 ] ; then
	echo "fucking timed out file: $file"
    fi
done

echo "check 4f_singleZnunu_leptonic"
for file in $(find /pnfs/desy.de/ilc/prod/ilc/mc-opt-3/ild/dst-merged/500-TDR_ws/4f_singleZnunu_leptonic/ILD_s5_o1_v02/ -name "*.slcio")
do
    timeout 60 lcio_event_counter $file > /dev/null
    if [ $? == 124 ] ; then
	echo "fucking timed out file: $file"
    fi
done

echo "check 4f_singleZnunu_semileptonic"
for file in $(find /pnfs/desy.de/ilc/prod/ilc/mc-opt-3/ild/dst-merged/500-TDR_ws/4f_singleZnunu_semileptonic/ILD_s5_o1_v02/ -name "*.slcio")
do
    timeout 60 lcio_event_counter $file > /dev/null
    if [ $? == 124 ] ; then
	echo "fucking timed out file: $file"
    fi
done

echo "check 4f_singleZsingleWMix_leptonic"
for file in $(find /pnfs/desy.de/ilc/prod/ilc/mc-opt-3/ild/dst-merged/500-TDR_ws/4f_singleZsingleWMix_leptonic/ILD_s5_o1_v02/ -name "*.slcio")
do
    timeout 60 lcio_event_counter $file > /dev/null
    if [ $? == 124 ] ; then
	echo "fucking timed out file: $file"
    fi
done

echo "check 4f_WW_leptonic"
for file in $(find /pnfs/desy.de/ilc/prod/ilc/mc-opt-3/ild/dst-merged/500-TDR_ws/4f_WW_leptonic/ILD_s5_o1_v02/ -name "*.slcio")
do
    timeout 60 lcio_event_counter $file > /dev/null
    if [ $? == 124 ] ; then
	echo "fucking timed out file: $file"
    fi
done

echo "check 4f_WW_semileptonic"
for file in $(find /pnfs/desy.de/ilc/prod/ilc/mc-opt-3/ild/dst-merged/500-TDR_ws/4f_WW_semileptonic/ILD_s5_o1_v02/ -name "*.slcio")
do
    timeout 60 lcio_event_counter $file > /dev/null
    if [ $? == 124 ] ; then
	echo "fucking timed out file: $file"
    fi
done

echo "check 4f_ZZ_leptonic"
for file in $(find /pnfs/desy.de/ilc/prod/ilc/mc-opt-3/ild/dst-merged/500-TDR_ws/4f_ZZ_leptonic/ILD_s5_o1_v02/ -name "*.slcio")
do
    timeout 60 lcio_event_counter $file > /dev/null
    if [ $? == 124 ] ; then
	echo "fucking timed out file: $file"
    fi
done

echo "check 4f_ZZ_semileptonic"
for file in $(find /pnfs/desy.de/ilc/prod/ilc/mc-opt-3/ild/dst-merged/500-TDR_ws/4f_ZZ_semileptonic/ILD_s5_o1_v02/ -name "*.slcio")
do
    timeout 60 lcio_event_counter $file > /dev/null
    if [ $? == 124 ] ; then
	echo "fucking timed out file: $file"
    fi
done

echo "check 4f_ZZWWMix_leptonic"
for file in $(find /pnfs/desy.de/ilc/prod/ilc/mc-opt-3/ild/dst-merged/500-TDR_ws/4f_ZZWWMix_leptonic/ILD_s5_o1_v02/ -name "*.slcio")
do
    timeout 60 lcio_event_counter $file > /dev/null
    if [ $? == 124 ] ; then
	echo "fucking timed out file: $file"
    fi
done

echo "check finished"
