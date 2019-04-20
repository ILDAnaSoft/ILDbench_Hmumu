#!/bin/bash

echo "START hadd..."

echo "combining part of 6f_ttbar"
hadd 6f_ttbar6.root 6f_ttbar/6f_ttbar-[5-7]??.root

echo "combining 6f_yyyyZ"
hadd yyyyZ.root 6f_yyyyZ/*.root

echo "combine everything"
hadd alldata6.root 6f_ttbar6.root yyyyZ.root

echo "delete intermediate files"
rm 6f_ttbar6.root yyyyZ.root

echo "finished!"
