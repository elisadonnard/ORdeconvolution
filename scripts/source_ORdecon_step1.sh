#!/bin/bash

# convert OR reads to bed
for i in or_reads*; do bedtools bamtobed -split -i $i > bed_$i; done
# print read name for reads that map to more than one OR
bedtools intersect -wo -a $1 -b bed_* | awk -F "\t" '{print $4,$17}' | sort | uniq | awk '{print $2}' | sort | uniq -d > multimapped_and_gap_readnames
# get OR name for reads that map to more than one OR
bedtools intersect -wo -a $1 -b bed_* | fgrep -w -f multimapped_and_gap_readnames | awk '{print $4,$17}' | sort | uniq > multimapped_and_gap_readnames_ORnames
# parse the OR pairs and print out pairs that have at least two reads supporting
Rscript parse_OR_pairs.R
