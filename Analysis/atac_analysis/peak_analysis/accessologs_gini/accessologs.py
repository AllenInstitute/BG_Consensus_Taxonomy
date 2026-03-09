import pandas as pd
import numpy as np
import anndata as ad
from tqdm import tqdm
import h5py
import os
import csv
import re

##
base_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/basal-ganglia/data/"
analysis_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis/atac/accessologs/"

## -------------------------------------------
## Load ortholog peak ids
work_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/"
anno_dir = os.path.join(work_dir, "analysis", "annotations")
liftover = pd.read_csv(os.path.join(anno_dir, "human_ref_liftover_HALPER_minMatch_0-5.tsv"), sep="\t")
liftover["region"] = liftover["human_ID"].copy()
liftover.set_index("region", inplace=True)

##
ortholog_peaks = liftover.loc[liftover.ortholog == True,:]

## Macaque chrom info 
macaque_chrom_alias = pd.read_csv("/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/genomes/macaque/ncbi/rheMac10.chromAlias.txt", sep="\t", header=None, comment="#")

## -------------------------------------------
## Identify species files, could also be a list from glob
species_h5ads = {
    "human" : os.path.join(base_dir, "human", "crested_adata", "human_basalganglia_hmba_pre-print_crested_500bp.h5ad"),
    "macaque" : os.path.join(base_dir, "macaque", "crested_adata", "macaque_basalganglia_hmba_pre-print_crested_500bp.h5ad"),
    "marmoset": os.path.join(base_dir, "marmoset", "crested_adata", "marmoset_basalganglia_hmba_pre-print_crested_500bp.h5ad"),
}

## -------------------------------------------
## Load accessibility per species and compute accessologs
expr = {}; genes={}; celltypes = {}; adatas = {}
for species in species_h5ads:
    ## Load species data
    print(f"Loading {species}")
    adata = ad.read_h5ad(species_h5ads[species])
    ## Apply peak lookup to adata.var_names
    if species != "human":
        ## Gather just species ortholog peak names
        peak_lookup_species = ortholog_peaks[["human_ID", f"{species}_ID"]].copy()
        peak_lookup_species = peak_lookup_species.loc[peak_lookup_species[f"{species}_ID"].notna()]
        peak_lookup_species = dict(zip(peak_lookup_species[f"{species}_ID"], peak_lookup_species["human_ID"]))
        ## For non-human species, map the gene names to human orthologs
        adata.var["region"] = adata.var_names
        adata.var["human_peak_name"] = adata.var["region"].map(peak_lookup_species)
    else:
        ## For human species, use the existing names
        adata.var["human_peak_name"] = adata.var_names
    ## Update Macaque using chrom alias
    # if species == "macaque":
    #     chrom_alias_dict = dict(zip(macaque_chrom_alias[0], macaque_chrom_alias[1]))
    #     adata.var["chr"] = adata.var["chr"].replace(chrom_alias_dict)
    ## Filter to only ortholog peaks
    adata = adata[:,adata.var["human_peak_name"].notna()].copy()
    ## Set index to human peak names
    adata.var.set_index("human_peak_name", inplace=True, drop=False)
    ## Print shape
    print(f"Species {species} has {adata.n_vars} ortholog peaks and {adata.n_obs} cell types")
    ## Store data
    expr[species] = pd.DataFrame(adata.X, 
                                        index=adata.obs_names, 
                                        columns=adata.var["human_peak_name"]).T
    genes[species] = list(adata.var["human_peak_name"])
    celltypes[species] = list(adata.obs_names)
    adatas[species] = adata

## Save data
for species, data in adatas.items():
    print(f"Saving expression data for {species} with shape {data.shape}")
    data.write_h5ad(os.path.join(analysis_dir, f"{species}_crested_adata_ortholog_accessibility.h5ad"))

## Determine common genes
common_genes = [set(gene) for gene in genes.values()]
common_genes = list(set.intersection(*common_genes))

## Determine cell types
common_celltypes = [set(ct) for ct in celltypes.values()]
common_celltypes = list(set.union(*common_celltypes))

## For each species, ensure all cell types exist. If not then set to 0
for species_name in expr.keys():
  missing_ct = list(set(common_celltypes) - set(celltypes[species_name]))
  print(f"Missing cell types in {species_name}: {missing_ct}")
  if len(missing_ct) > 0:
    for ct in missing_ct:
        expr[species_name].loc[:,ct] = 0
  expr[species_name] = expr[species_name].loc[:,common_celltypes].T
  expr[species_name] = expr[species_name].loc[:, ~expr[species_name].columns.duplicated()]

## -------------------------------------------
## Expressolog data computation
accessolog = []
gene_universe = common_genes
species = list(expr.keys())

## Compute accessologs
for i in range(len(species) - 1):
    sp1 = species[i]
    print(sp1)
    for j in range(i + 1, len(species)):
        sp2 = species[j]
        expr_cor = np.empty((len(gene_universe), 3))
        for idx, gene in enumerate(tqdm(gene_universe, desc="Processing ortholog peaks")):
            ## 
            if gene not in expr[sp1].columns:
                sp1_expr = np.zeros([len(common_celltypes),1])
            else:
                sp1_expr = expr[sp1].loc[:,gene].to_numpy().reshape(-1, 1)
            ##
            if gene not in expr[sp2].columns:
                sp2_expr = np.zeros([len(common_celltypes),1])
            else:
                sp2_expr = expr[sp2].loc[:, gene].to_numpy().reshape(-1, 1) ## Revist, why is peak duplicated
            ##
            sp1_expr_stat = np.std(sp1_expr)
            sp2_expr_stat = np.std(sp2_expr)
            cor1 = np.corrcoef(sp1_expr[:, 0], sp2_expr[:, 0])[0, 1]
            expr_cor[idx, :] = [sp1_expr_stat, sp2_expr_stat, cor1]
        ##
        corm = pd.DataFrame({
            'region': gene_universe,
            'species1': sp1,
            'species2': sp2,
            'sp1_expr_stat': expr_cor[:, 0],
            'sp2_expr_stat': expr_cor[:, 1],
            'expr_cor': expr_cor[:, 2],
            'sp_order': 1
        })
        ## Symmetric comparison
        corm2 = pd.DataFrame({
            'region': gene_universe,
            'species1': sp2,
            'species2': sp1,
            'sp1_expr_stat': expr_cor[:, 1],
            'sp2_expr_stat': expr_cor[:, 0],
            'expr_cor': expr_cor[:, 2],
            'sp_order': 2
        })
        ##
        accessolog.append(pd.concat([corm, corm2], ignore_index=True))

## Combine all results into a single DataFrame
accessologs = pd.concat(accessolog, ignore_index=True)

## Compute rank and auroc
accessologs['median_cor'] = accessologs.groupby('region')['expr_cor'].transform('median')
accessologs = accessologs.sort_values(by='median_cor', ascending=True).reset_index(drop=True)
accessologs['rank'] = accessologs['median_cor'].rank(method='dense', ascending=True)
accessologs['auroc'] = accessologs['rank'] / accessologs['region'].nunique()

## -------------------------------------------
## Data saviing and organization

## Get species order from metadata
species_order = pd.Series(["human", "macaque", "marmoset"]).tolist()

## Convert species1 and species2 to categorical data types with the given order
accessologs['species1'] = pd.Categorical(accessologs['species1'], categories=species_order, ordered=True)
accessologs['species2'] = pd.Categorical(accessologs['species2'], categories=species_order, ordered=True)

## Write accessologs data.frame to csv
accessologs.to_csv(os.path.join(analysis_dir, "accessologs_primate_orthologs.csv"), index=False)



