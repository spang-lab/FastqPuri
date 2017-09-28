#!/bin/bash

# Set variables
trimFilter="../../../bin/trimFilter"
adevL="../../fa_fq_files/adapter_even_long.fa"
adevS="../../fa_fq_files/adapter_even_short.fa"
adodL="../../fa_fq_files/adapter_odd_long.fa"
adodS="../../fa_fq_files/adapter_odd_short.fa"
fqev_adevL="../../fa_fq_files/human_even_wad_even_long.fq"
fqev_adevS="../../fa_fq_files/human_even_wad_even_short.fq"
fqev_adodL="../../fa_fq_files/human_even_wad_odd_long.fq"
fqev_adodS="../../fa_fq_files/human_even_wad_odd_short.fq"
fqod_adevL="../../fa_fq_files/human_odd_wad_even_long.fq"
fqod_adevS="../../fa_fq_files/human_odd_wad_even_short.fq"
fqod_adodL="../../fa_fq_files/human_odd_wad_odd_long.fq"
fqod_adodS="../../fa_fq_files/human_odd_wad_odd_short.fq"


# Run even fq with long even adapter
$trimFilter --ifq ${fqev_adevL}  --length 150 --adapter ${adevL}:2:5.3 \
            --output evev_long
$trimFilter --ifq ${fqod_adevL}  --length 151 --adapter ${adevL}:2:5.3 \
            --output odev_long
$trimFilter --ifq ${fqod_adodL}  --length 149 --adapter ${adodL}:2:5.3 \
            --output odod_long

