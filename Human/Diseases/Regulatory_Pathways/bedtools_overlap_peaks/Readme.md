
# SNP–Peak Overlap Pipeline

This pipeline performs an overlap between genome-wide significant SNPs and cell-type-specific chromatin accessibility peaks predicted using the CRESTED model.

It is designed to identify which disease-associated SNPs fall within regulatory regions specific to a given BG cell type.

---
## Input Format — Significant SNPs in .bed format

To perform the overlap, you need a BED file containing significant SNPs (PVAL < 5e-8), formatted as follows: (Ensure that the SNP coordinates in the BED file are based on the same reference genome as the peak coordinates):

| Column | Description                                                  |
|--------|--------------------------------------------------------------|
| chr    | Chromosome name, prefixed with `chr` (e.g., `chr1`)         |
| start  | SNP position minus 1 (**0-based start**, as per BED format) |
| end    | SNP position (**1-based end**)                               |
| name   | SNP identifier (e.g., `rsID`)

---
## Overlap with CRESTED-predicted peaks

Once the .bed file is correctly formatted, run the following command to intersect it with predicted chromatin accessibility peaks for a specific cell type (e.g., Lamp5):

```bash
bedtools intersect \
  -a path_to_disease_significant_snps/significant_snps.bed \
  -b path_to_cell_type_specific_peaks/Lamp5_prioritized_peaks.bed  \
  -wa -wb > outcome_dir_path/snps_in_Lamp5_predicted_peaks.bed
```
---
The output contains one row per SNP–peak overlap, with columns indicating the SNP coordinates and ID (snp_chr, snp_start, snp_end_POS, rsid), followed by the overlapping peak region (peak_chr, peak_start, peak_end).
