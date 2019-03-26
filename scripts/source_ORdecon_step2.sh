#!/bin/bash

# Select read counts per umi for OR genes only
for i in *.umi.distributions.txt; do grep Olfr $i > ORs_$i; done
# Correct umi assignments when the same umi or a umi that is one hamming distance away maps to more than one OR
for i in ORs_*; do python cleanLowEndUmis_ORdeconvolution.py -i $i -o clean_$i -n 2; done
# 

