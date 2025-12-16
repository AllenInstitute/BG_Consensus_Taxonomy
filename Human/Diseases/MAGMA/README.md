## Magma Analyses : Running the enrichment

### Data requirements (Inspired by Duncan et al.) :

**1. GWAS summary statistics**  
The summary statistics for the different phenotypes can be found in that directory: `/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/BasalGanglia/external/Disease_related/raw_GWAS_data/`

**2. Updated Basal Ganglia Dataset** 
`/allen/programs/celltypes/workgroups/rnaseqanalysis/hct_ux3_cellxgene/anndata_080/BasalGanglia/BICAN-releases/final/Human_HMBA_basalganglia_AIT_pre-print.h5ad`

**3. MAGMA auxiliary files of the same genome build and ancestry as GWAS summary statistics**  
MAGMA auxiliary files of the same genome build and ancestry as GWAS summary statistics to be used.  
Details about the genome build and ancestry for each GWAS phenotype are provided here:[Recap GWAS Data](https://docs.google.com/spreadsheets/d/11oP8UJAEofMAzROBXA0sFL90YaKX98-JvVZQIWrHim4/edit?usp=sharing)

- Gene locations, build 37 and 38 files: NCBI37.3.gene.loc and NCBI38.3.gene.loc (can be downloaded from [here](https://cncr.nl/research/magma/))
- The extended MHC was removed from the auxiliary gene location files (extendedMHCexcluded), for the Step 2 of the "MAGMA input preparation".
For build 37, the file NCBI37.3.gene.loc.extendedMHCexcluded comes from [Bryois et al](https://github.com/jbryois/scRNA_disease/blob/master/Code_Paper/Data/NCBI/NCBI37.3.gene.loc.extendedMHCexcluded).

- Reference data per ancestry: Can be downloaded from [here](https://cncr.nl/research/magma/). As an example for the the European Ancestry, download the folder g1000_eur which contains:  
  g1000_eur.bed  
  g1000_eur.bim  
  g1000_eur.fam  
  g1000_eur.synonyms


### Get MAGMA inputs:

**For the BG dataset:**  
1/ [Get the ln(1+x)-transformed cluster-by-gene matrix](https://github.com/AllenInstitute/HMBA_Genomics/blob/main/BasalGanglia/Human/Diseases/MAGMA/Preparation_input_MAGMA/create_L2-log_dataset.py)  
2/ [Preprocess the matrix and calculate specificity](https://github.com/AllenInstitute/HMBA_Genomics/blob/main/BasalGanglia/Human/Diseases/MAGMA/Preparation_input_MAGMA/get_Siletti_continuous_input.R)  

Set the `base_dir_BG_dataset` variable in the `GWAS_config.env` to the folder where your processed Basal Ganglia data is stored.


**For the GWAS summary statistics**  
1/ Check that GWAS summary statistics files have the following columns: rsID, chromosome, base pair location, p-value, and n
(n = Nb of cases + Nb of controls) for the phenotype and ancestry of interest  
Include a header row on this file. This is the file with the suffix _no_heading.  
2/ Create a SNP location file (snploc_{GWAS_file_name})  
This file should contain three columns of the GWAS summary statistics in the following order: rsID, chromosome, and base pair position
Remove any header rows from this file  

Set the `folder` variable in the `GWAS_config.env` to the folder where your processed GWAS summary data is stored (no_heading file and snp_loc file):
`/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/BasalGanglia/external/Disease_related/MAGMA_results/{phenotype}_final_taxonomy/{ancestry}/`

### Save MAGMA outputs:
All the MAGMA intermediate files and outputs are saved in the same directory : `/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/BasalGanglia/external/Disease_related/MAGMA_results/{phenotype}_final_taxonomy/{ancestry}/`

MAGMA outputs: 
- `{file}.step1.genes.annot` Contains the list of SNPs mapped to each gene based on physical proximity 
- `{file}.step2.genes.out` Provides gene-level association results. Each row corresponds to a gene, with summary statistics including the p-value (PVALUE), number of SNPs, and mean chi-squared score.
- `final_taxonomy_BG_l2_conti-spe_{phenotype}.gsa.out`Contains the cluster-level association results. Each row corresponds to a cluster, with summary statistics, number of genes, PVALUE
