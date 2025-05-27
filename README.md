### Scripts for *Ips typographus* population genomics data analysis
Bioinformatics pipeline for population genomics analysis using Illumina Pool-seq data from DNA sequencing of *Ips typographus* (spruce bark beetle).

#### Directory structure
- `ips_poolseq/` – Scripts for data analisys
	- `01_proc_reads/` – Run `fastp`, generate summary tables and HTML reports
	- `02_mappings/` – Align with `bwa-mem2` and sort with `samtools`