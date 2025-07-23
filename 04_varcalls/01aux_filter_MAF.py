#!/usr/bin/env python3

#-------------------------------------------------------------------------------
# Filter the VCF file by : min(p) <= (1 - min_MAF) and max(p) >= min_MAF
#
# Input:
#   - ips_merged.vcf.gz: vcf file for all samples
# Output:
#   - ips_merged.m05.vcf.gz: vcf file filtered by min_MAF for all samples
#-------------------------------------------------------------------------------

import gzip

work_dir = "../../results/04_varcalls/"
min_MAF = 0.05

# Input file
vcf_filter_in = f"{work_dir}/ips.biallelic_q20_m20.vcf.gz"
# Output file
vcf_m05_out = f"{work_dir}/ips.biallelic_q20_m20.maf05.vcf.gz"

def calc_ref_freq(ad_field):
    ref, alt = map(int, ad_field.split(','))
    total = ref + alt
    if total == 0:
        return None
    return ref / total

def process_vcf(vcf_filter_in, vcf_m05_out, min_MAF):
    with gzip.open(vcf_filter_in, 'rt') as vcf_filter_fh, gzip.open(vcf_m05_out, "wt") as vcf_m05_fh:
        for line in vcf_filter_fh:
            if line.startswith("#"):
                    vcf_m05_fh.write(line)
                    continue                                  

            fields = line.strip().split('\t')            
            samples = fields[9:]

            p = []
            for sample in samples:
                ad_field = sample.split(":")[2]                
                freq = calc_ref_freq(ad_field)
                if freq is not None:
                    p.append(freq)

            if not p:
                continue

            if min(p) <= (1 - min_MAF) and max(p) >= min_MAF:
                vcf_m05_fh.write(line)

def main():
    process_vcf(vcf_filter_in, vcf_m05_out, min_MAF)

if __name__ == "__main__":
    main()

