#!/bin/bash


EColi_gen="../fa_fq_files/EColi_genome.fa"
rRNA_seq="../fa_fq_files/rRNA_modified.fa"
NrEC=10000  # Number of reads from Ecoli genome
NrRNA=5000 # Number of reads from rRNA sequence 
Lr=70       # Read length
Lf=100      # fragment length
sdf=10      # fragment length standard deviation
bwa=1       # output per-read-end (bwa) only
name_EC="EColi"  # prefix of EColi reads output
name_rRNA="rRNA" # prefix of rRNA reads output
name_both="EColi_rRNA" # prefix for merged data output
suffix1=".bwa.read1.fastq.gz"
suffix2=".bwa.read2.fastq.gz"

# USE Claudios tool 
# Generate EColi reads
../fa_fq_files/simReads.pl ${EColi_gen} 0.01 0 70 0 > ${name_EC}.read1.fq
../fa_fq_files/simReads.pl ${EColi_gen} 0.01 0 70 0 > ${name_EC}.read2.fq

# Generates rRNA reads
../fa_fq_files/simReads.pl ${rRNA_seq} 0.005 0 70 0 > ${name_rRNA}.read1.fq
../fa_fq_files/simReads.pl ${rRNA_seq} 0.005 0 70 0 > ${name_rRNA}.read2.fq

# Concatenate the two 
cat ${name_EC}.read1.fq ${name_rRNA}.read1.fq | gzip > ${name_both}.read1.fq.gz
cat ${name_EC}.read2.fq ${name_rRNA}.read2.fq | gzip > ${name_both}.read2.fq.gz


## Generate EColi reads
#dwgsim -o $bwa -1 $Lr -2 $Lr -N $NrEC -s $sdf -d $Lf  ${EColi_gen} ${name_EC}
#
## Generates rRNA reads
#dwgsim -o $bwa -1 $Lr -2 $Lr -N $NrRNA -s $sdf -d $Lf -r 0 ${rRNA_seq} ${name_rRNA}
#
## Concatenate the two 
#zcat ${name_EC}$suffix1 ${name_rRNA}$suffix1 | gzip > ${name_both}.read1.fq.gz
#zcat ${name_EC}$suffix2 ${name_rRNA}$suffix2 | gzip > ${name_both}.read2.fq.gz
