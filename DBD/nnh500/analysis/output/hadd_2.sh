#!/bin/bash

echo "START hadd..."

echo "combining part of 6f_ttbar"
hadd 6f_ttbar2.root 6f_ttbar/6f_ttbar-1[0-1]??.root

echo "combining 6f_eeWW, llWW, xxxxZ"
hadd 6f_1.root 6f_eeWW/*.root 6f_llWW/*.root 6f_xxxxZ/*.root

echo "combine everything"
hadd alldata2.root 6f_ttbar2.root 6f_1.root

echo "delete intermediate files"
rm 6f_ttbar2.root 6f_1.root

echo "finished!"
