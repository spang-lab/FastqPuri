# This script calls the executable and saves the output 
fastq1=EColi_rRNA.read1.fq.gz
fastq2=EColi_rRNA.read2.fq.gz
fasta=../fa_fq_files/rRNA_modified.fa
trimFilterDS=../../bin/trimFilterDS
trimFilter=../../bin/trimFilter


# Run trimFilterDS using a tree
$trimFilterDS -l 70 --ifq "${fastq1}:${fastq2}" --ifa "$fasta":0.10:25 \
 --method TREE --trimQ ENDSFRAC --trimN STRIP -o treeDS -q 5 --adapter \
 Ad_read1.fa:Ad_read2.fa:2:40

# Run trimFilter using a tree
#$trimFilter -l 70 --ifq "${fastq1}" --ifa "$fasta":0.10:25 \
# --method TREE --trimQ ENDSFRAC --trimN STRIP -o tree -q 5 
