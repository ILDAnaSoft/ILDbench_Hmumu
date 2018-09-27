#!/bin/bash

echo "START hadd..."

echo "combining 2 files"
hadd alldata_select.root alldata1_select.root alldata2_select.root

echo "delete intermediate files"
rm alldata1_select.root alldata2_select.root alldata1.root alldata2.root

echo "finished!"
