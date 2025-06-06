#!/usr/bin/env bash

dir="../../data/03_dedup"
out="../../results/04_varcalls"

time head -n 10  regions_10k | parallel -j10 "region={}; chunk=\$(echo \$region | sed 's/[:-]/_/g'); bcftools mpileup -Ou -I -a FORMAT/AD --max-depth 5000 -r \$region -f ../../data/reference/Ips_typograpgus_LG16corrected.final.fasta ../../data/03_dedup/WFIN_Ea_2015.dedup.sort.bam  ../../data/03_dedup/WFIN_Eb_2015.dedup.sort.bam  | bcftools call -vmO v -o ../../results/04_varcalls/WFIN_E_2015.dedup.parallel.\$chunk.vcf" 1>log1_1 2>log2_1

time bcftools mpileup --threads 10 -Ou -I -a FORMAT/AD --max-depth 5000 -r LG1:1-100000 -f ../../data/reference/Ips_typograpgus_LG16corrected.final.fasta $dir/WFIN_Ea_2015.dedup.sort.bam  $dir/WFIN_Eb_2015.dedup.sort.bam  | bcftools --threads 10 call -vmO v -o $dir/WFIN_E_2015.dedup.threads.vcf 1>log_1_2 2>log2_2
