#!/usr/bin/env python3

#----------------------------------------------------------------------
# Compares fluctuating SNPs across 4 populations in 3 ways:
#   - Exact SNP overlap (thinned) - match on CHROM_POS
#   - Gene-level overlap (thinned) - match on gene name 
#   - Exact SNP overlap (pre-thinning)- match on CHROM_POS
#
# Outputs: 3 TSV tables with pairwise overlaps
#  - overlap_thinned_snp.tsv
#  - overlap_thinned_gene.tsv
#  - overlap_prethinned_snp.tsv
#----------------------------------------------------------------------

import os

prefixes = ["SFIN_E", "SFIN_L", "WFIN_E", "WFIN_L"]
thinned_dir = "../results/08_models/s2_thinning"
filter_dir = "../results/08_models/s1_filter"
out_dir = "../results/08_models/overlap"

# Loaders
# Return set of (CHROM, POS) from thinned file.
def load_thinned_snps(pre):   
    path = f"{thinned_dir}/tests.{pre}.FDR.fluctuating.thinned.tsv"
    snps = set()
    with open(path) as f:
        for line in f:
            cols = line.strip().split("\t")
            snps.add((cols[1], cols[2]))
    return snps

#Return set of gene names (col 0) from thinned file.
def load_thinned_genes(pre):    
    path = f"{thinned_dir}/tests.{pre}.FDR.fluctuating.thinned.tsv"
    genes = set()
    with open(path) as f:
        for line in f:
            cols = line.strip().split("\t")
            genes.add(cols[0])
    return genes

#Return set of (CHROM, POS) from pre-thinning FDR file (skip header).
def load_prethinned_snps(pre):    
    path = f"{filter_dir}/tests.{pre}.FDR.fluctuating.tsv"
    snps = set()
    with open(path) as f:
        f.readline()  # skip header
        for line in f:
            cols = line.strip().split("\t")
            snps.add((cols[0], cols[1]))
    return snps

# Build matrix: diagonal = population size, lower triangle = shared (% of smaller pop).
def build_table(sets, label): 
    # header row
    rows = [[""] + prefixes]

    for i, p1 in enumerate(prefixes):
        row = [p1]
        for j, p2 in enumerate(prefixes):
            if j > i:
                row.append("") # upper triangle empty
            elif j == i:
                row.append(str(len(sets[p1]))) # diagonal = size
            else:
                shared = len(sets[p1] & sets[p2])
                smaller = min(len(sets[p1]), len(sets[p2]))
                pct = round(100 * shared / smaller, 1)
                row.append(f"{shared} | {pct}%") # lower triangle
        rows.append(row)
    return rows

def write_tsv(rows, path):
    if not rows:
        return
    with open(path, "w") as f:
        for row in rows:
            f.write("\t".join(row) + "\n")

def main():
    os.makedirs(out_dir, exist_ok=True)

    # - Exact SNP overlap (thinned) 
    thinned_snps = {pre: load_thinned_snps(pre) for pre in prefixes}    
    rows1 = build_table(thinned_snps, "thinned_snp")
    write_tsv(rows1, f"{out_dir}/overlap_thinned_snp.tsv")

    # - Gene-level overlap (thinned)
    thinned_genes = {pre: load_thinned_genes(pre) for pre in prefixes}    
    rows2 = build_table(thinned_genes, "thinned_gene")
    write_tsv(rows2, f"{out_dir}/overlap_thinned_gene.tsv")

    # - Exact SNP overlap (pre-thinning)
    prethinned_snps = {pre: load_prethinned_snps(pre) for pre in prefixes}    
    rows3 = build_table(prethinned_snps, "prethinned_snp")
    write_tsv(rows3, f"{out_dir}/overlap_prethinned_snp.tsv")

if __name__ == "__main__":
    main()