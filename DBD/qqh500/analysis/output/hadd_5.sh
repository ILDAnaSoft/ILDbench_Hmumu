#!/bin/bash

echo "START hadd..."

echo "combining part of 6f_ttbar"
hadd 6f_ttbar5.root 6f_ttbar/6f_ttbar-[2-4]??.root

echo "combining higgs_ffhf"
hadd higgs_ffh.root higgs_ffh/*.root

echo "combine everything"
hadd alldata5.root 6f_ttbar5.root higgs_ffh.root

echo "delete intermediate files"
rm 6f_ttbar5.root higgs_ffh.root

echo "finished!"
