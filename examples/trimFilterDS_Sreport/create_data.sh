#!/bin/bash


EColi_gen="../fa_fq_files/EColi_genome.fa"
rRNA_seq="../fa_fq_files/rRNA_modified.fa"
NrEC=10000  # Number of reads from Ecoli genome
NrRNA=5000 # Number of reads from rRNA sequence 
Lr=70       # Read length
Lf=100      # fragment length
sdf=10      # fragment length standard deviation
bwa=1       # output per-read-end (bwa) only
pre="EColi_rRNA" # prefix of the output file
suf1=".read1.fq" #suffix output1
suf2=".read2.fq" #

# USE Claudios tool 
# Generate EColi reads
../fa_fq_files/simReads.pl ${EColi_gen} 0.01 0 70 0 > ${name}$suf1
../fa_fq_files/simReads.pl ${EColi_gen} 0.01 0 70 0 > ${name}$suf2

# Generates rRNA reads
../fa_fq_files/simReads.pl ${rRNA_seq} 0.005 0 70 0 >> ${name}$suf1
../fa_fq_files/simReads.pl ${rRNA_seq} 0.005 0 70 0 >> $name$suf2


# Check --trimQ 
echo "@read ENDS trimQ:3:67, FRAC discard, ENDSFRAC trimQ:3:67

"
