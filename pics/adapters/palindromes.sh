#!/bin/bash 

# We describe the cases in palindrome mode

# Definition of colors
n="\e[0m"
b="\e[1m"
cr1="\e[33;1m" cr2="\e[31;1m"
ca1="\e[34;1m" ca2="\e[32;1m"
cr1r2="\e[43;1m"


#Adapters
a1fw="CAGTTGG";a2fw="AAGGCTC"; 
a1bw="CCAACTG";a2bw="GAGCCTT";
a2rv=$ca1`echo -n "$a2fw"`$n; 
a1fw=$ca1$a1fw$n; a1bw=$ca1$a1bw$n; 
a2fw=$ca2$a2fw$n; a2bw=$ca2$a2bw$n; 

printf "Adapter 1: $a1fw\tReverse complement: $a1bw\n"
printf "Adapter 2: $a2fw\tReverse complement: $a2fw\n"

# Case 1: no overlap, no contamination
printf "${b}CASE1:$n No overlap, no contamination.\n"
r1fw="AAAATATAAAGGGGTCTG"; r1bw="CAGACCCCTTTATATTTT"
r2fw="CCCGATTTTACCCCGAGG"; r2bw="CCTCGGGGTAAAATCGGG"
gapfw="CACGTACGAC"; gapbw="GTCGTACGTG"
r1fw=$cr1$r1fw$n; r1bw=$cr1$r1bw$n; 
r2fw=$cr2$r2fw$n; r2bw=$cr2$r2bw$n;
gapfw=$b$gapfw$n
gapbw=$b$gapbw$n
empty="$b------$n"
printf " * Full fragments in both senses (with adapters) \n"
printf "   1.  $a1fw$r1fw$gapfw$r2bw$a2bw \n"
printf "   2.  $a2fw$r2fw$gapbw$r1bw$a1bw \n"
printf " * Fastq files\n"
printf "   1. $r1fw\n" 
printf "   2. $r2fw\n" 
printf " * Sequence alignment: \n"
printf "    $a1fw$r1fw$empty\n"
printf "\t\t\t     $empty$r2bw$a2bw\n"
printf " * Output:  nothing done\n"

# Case 2: overlap, no contamination
printf "${b}CASE2:$n No overlap, no contamination.\n"
r1fw="AAAATATAAAGGGGTCTG"; r1bw="CAGACCCCTTTATATTTT"
r2fw="CCCGATTTTACCCCGAGG"; r2bw="CCTCGGGGTAAAATCGGG"
ovfw="CACGTACGAC"; ovbw="GTCGTACGTG"
r2rv=$cr2`echo -n "$r2fw" | rev`$n
ovrv=$cr1r2`echo -n "$ovbw" | rev`$n
r1fw=$cr1$r1fw$n; r1bw=$cr1$r1bw$n; 
r2fw=$cr2$r2fw$n; r2bw=$cr2$r2bw$n;     
ovfw=$cr1r2$ovfw$n; ovbw=$cr1r2$ovbw$n; 
printf " * Full fragments in both senses (with adapters) \n"
printf "   1.  $a1fw$r1fw$ovfw$r2bw$a2bw \n"
printf "   2.  $a2fw$r2fw$ovbw$r1bw$a1bw \n"
printf " * Fastq files\n"
printf "   1. $r1fw$ovfw\n" 
printf "   2. $r2fw$ovbw\n" 
printf " * Sequence alignment: \n"
printf "    $a1fw$r1fw$ovfw\n"
printf "\t\t\t     $ovrv$r2rv$a2rv\n"
printf " * Output:  nothing done\n"


#to test: 
#for i in `seq 30 47; seq 90 96`; do for j in `seq 1 30`;  do printf "\e[${i};${j}m $i $j; \e[0m"; done ; echo ""; done

