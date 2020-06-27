#!/bin/bash

# This script creates a bloom filter from rRNA (FPR 0.0075), and runs trimFilter 
# as in the example trimFilter_Sreport/ but using the bloom filter with a score
# of 0.4 instead of the tree.
#
# This script helps automatid testing

cd examples
fastqPath=fa_fq_files/
fastq=EColi_rRNA.fq.gz
fasta=$fastqPath/rRNA_modified.fa
sequence=$fastqPath/rRNA_modified.fa

output=testdir
mkdir -p $output
bf=$output/rRNA

Qreport -i $fastqPath/$fastq -l 50 -o $output/$fastq

# trim with tree
trimFilter -l 50 --ifq "$fastqPath/$fastq" --ifa "$fasta":0.9:50 --method TREE --trimQ ENDSFRAC --trimN STRIP -o $output/tree
Qreport -i $output/tree_good.fq.gz -l 50 -o $output/tree_good


# trim with bloom filter
makeBloom -f $sequence -o $bf -p 0.0075
trimFilter -f $fastqPath/$fastq -l 50 --method BLOOM --idx $bf.bf:0.4 --trimQ ENDSFRAC --trimN ENDS -o $output/bloom
Qreport -i $output/bloom_good.fq.gz -l 50 -o $output/bloom_good
