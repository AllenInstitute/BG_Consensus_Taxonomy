## CRESTED Analyses: prediction of the impact of SNPs mutations on chromatin accessibility

### Data Requirements

To run the CRESTED-based SNP impact pipeline, the following data are required:

- **Genome FASTA**: A reference genome file indexed for sequence retrieval.
- **Pre-trained CRESTED model**: The model used to generate regulatory predictions.
- **SNP–peak annotation DataFrame**: A table containing genome-wide significant SNPs along with peak and allele annotations. The following columns are required:

  | Column               | Description                                           |
  |----------------------|-------------------------------------------------------|
  | `snp_chr`            | Chromosome of the SNP                                |
  | `snp_end_POS_hg38`   | SNP end position (hg38 coordinates)                  |
  | `rsid`               | SNP identifier                         |
  | `peak_chr`           | Chromosome of the overlapping regulatory peak        |
  | `peak_start`         | Start coordinate of the peak                         |
  | `peak_end`           | End coordinate of the peak                           |
  | `effect_allele`      | Allele associated with the phenotype (effect allele) |
  | `other_allele`       | Alternative allele                                   |
  | `PVAL`               | GWAS p-value for the SNP                             |

For the analysis of the preprint manuscript, SNP–peak annotation DataFrames are available here : `/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/BasalGanglia/external/Disease_related/CRESTED_results/snps_peaks_BG/`

### Run the Analysis of SNP Effects Across Cell Types

#### Summary View of SNPs effects across cell types
The notebook `RECAP_effects_all_SNPs.ipynb` takes the prepared input data and computes the predicted regulatory effect of each SNP on chromatin accessibility across cell types.

The core function `plot_variant_specific_scores` compares the chromatin accessibility predictions of the reference (non-mutated) and mutated versions of each peak sequence using a CRESTED model.

The notebook identifies the top 20 SNPs with the largest absolute effect (the greatest predicted change in accessibility), and visualizes them using a dot plot:
- **Color** encodes the direction of the effect (positive or negative)
- **Dot size** reflects the magnitude of the effect across cell types


#### Zoom-In on a Specific Variant

The notebook `CRESTED_prediction_score_CT_specificity_ISM_effect_allele.ipynb` allows focused analysis of a single variant in a given phenotype and cell type.

It takes as input:
- `pheno`: the phenotype of interest  
- `cell_type`: the cell type of interest  
- `snp`: the rsID of the variant  

The notebook performs the following:
- Computes contribution scores (attribution maps) for the reference and mutated sequences in the specified cell type.
- Visualizes how the mutation alters the predicted chromatin accessibility.
- Plots the distribution of predicted effects for all possible point mutations across the peak sequence (in silico mutagenesis).

This enables interpretation of SNP impact at single-base resolution in a cell-type-specific context.



### Accessing Contribution Score Matrices for the background across cell types:

To explore the effect of a given variant on all cell types (contribution scores background), the full matrix of contribution scores is available at the following locations:

- **Reference sequence**:  
  `/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/BasalGanglia/external/Disease_related/CRESTED_results/contribution_reference_scores_{snp}_{cell_type}_{pheno}.h5`

- **Mutated sequence**:  
  `/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/BasalGanglia/external/Disease_related/CRESTED_results/contribution_mutated_scores_{snp}_{cell_type}_{pheno}.h5`

> Replace `{snp_id}` with the corresponding SNP identifier (e.g., `rs1883123`).

Each `.h5` file contains contribution scores across all cell types, with dimensions `[cell_type x position x base]`.
