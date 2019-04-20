#!/bin/bash

echo "START hadd..."

echo "combining 1f_3f background"
hadd 1f_3f.root 1f_3f/1f_3f-???.root

echo "combine aa_2f background"
hadd aa_2f.root aa_2f/aa_2f-???.root

echo "combine 2 processes"
hadd alldata2.root 1f_3f.root aa_2f.root

echo "delete intermediate files"
rm 1f_3f.root aa_2f.root

echo "finished!"
