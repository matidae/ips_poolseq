#!/usr/bin/env bash

cat regions_1M | parallel -j40 "region={}; chunk=\$(echo \$region | sed 's/[:-]/_/g'); bcftools mpileup -Ou -I -a FORMAT/AD --max-depth 5000 -r \$region -f ../../data/reference/Ips_typograpgus_LG16corrected.final.fasta ../../data/03
_dedup/WFIN_Ea_2015.dedup.sort.bam  ../../data/03_dedup/WFIN_Eb_2015.dedup.sort.bam ../../data/03_dedup/WFIN_Ea_2017.dedup.sort.bam  ../../data/03_dedup/WFIN_Eb_2017.dedup.sort.bam ../../data/03_dedup/WFIN_Ea_2020.dedup.sort.bam  ../../dat
a/03_dedup/WFIN_Eb_2020.dedup.sort.bam ../../data/03_dedup/WFIN_Ea_2022.dedup.sort.bam  ../../data/03_dedup/WFIN_Eb_2022.dedup.sort.bam ../../data/03_dedup/WFIN_Ea_2024.dedup.sort.bam  ../../data/03_dedup/WFIN_Eb_2024.dedup.sort.bam  | bcf
tools call -vmO v -o ../../results/04_varcalls/WFIN_E_2015.dedup.{#}.\$chunk.vcf"