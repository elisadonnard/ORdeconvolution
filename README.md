# ORdeconvolution
Cleaning up Olfactory Receptor expression from Olfactory Sensory Neurons single cell data

## Getting Started

You will need the umi.distributions.txt file resulting from the [quantification from ESAT](https://github.com/garber-lab/inDrop_Processing)

### Prerequisites

R;
bedtools 2.27.1 or higher;
Reference files for OR gene annotation (used here gencode VM18 for mm10)

### Step 1: Identify similar ORs from aligned reads

Starting from the reads that map to OR genes, identify common OR genes that can lead to quantification errors due to similarity, based on multimapped reads and gapped reads that span more than one OR

```
#Input files necessary:
#OR_reads.bam for each sample
#gencodeVM18_3UTR_ext_ORs.bed reference file with OR positions in genome including a 1000bp 3'UTR extension

./source_ORdecon_step1.sh gencodeVM18_3UTR_ext_ORs.bed
```
the output will be a file named empirical_or_pairs_filtered

### Step 2: Reassign identical UMIs that mapped to 2 different ORs to the top expressed OR

The input files are the .umi.distribution files resulting from ESAT

```
./source_ORdecon_step2.sh
``` 

### Step 3: Correct the count files

```
Rscript source_ORdecon_step3.R clean_OR_sampleNNN.umi.distributions.txt 0.75 
```



