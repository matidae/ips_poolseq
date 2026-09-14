#!/usr/bin/env python3

#-------------------------------------------------------------------------------
# Filter the VCF file by : min(p) <= (1 - min_MAF) and max(p) >= min_MAF
#
# Input:
#   - $VARCALL_RESULTS/ips.biallelic_q30_m30.vcf.gz
# Output:
#   - $VARCALL_RESULTS/ips.biallelic_q30_m30.maf05.vcf.gz
#-------------------------------------------------------------------------------

import gzip
import sys
sys.path.append("./utils")
from utils import MIN_MAF, load_config, log

cfg = load_config()
work_dir = cfg["VARCALL_RESULTS"]

# Input file
vcf_filter_in = f"{work_dir}/ips.biallelic_q30_m30.vcf.gz"
# Output file
vcf_m05_out = f"{work_dir}/ips.biallelic_q30_m30.maf05.vcf.gz"

def calc_ref_freq(ad_field):
    ref, alt = map(int, ad_field.split(','))
    total = ref + alt
    if total == 0:
        return None
    return ref / total

def process_vcf(vcf_filter_in, vcf_m05_out, MIN_MAF):
    kept = 0
    total = 0
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
            total += 1
            if min(p) <= (1 - MIN_MAF) and max(p) >= MIN_MAF:
                vcf_m05_fh.write(line)
                kept += 1
    log(f"SNPs passing MAF filter: {kept}/{total} ({100*kept/total:.1f}%)")

def main():
    log(f"=== MAF filtering start (MIN_MAF={MIN_MAF}) ===")
    log(f"input:  {vcf_filter_in}")
    log(f"output: {vcf_m05_out}")
    process_vcf(vcf_filter_in, vcf_m05_out, MIN_MAF)
    log("=== MAF filtering complete ===")

if __name__ == "__main__":
    main()

