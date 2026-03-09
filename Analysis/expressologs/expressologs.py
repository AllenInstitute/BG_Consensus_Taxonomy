import pandas as pd
import numpy as np
import h5py
import os
import csv
import re
from tqdm import tqdm
import glob
import anndata as ad

wkdir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis/expressologs/"
data_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis/mean_expr/"

## -------------------------------------------
## mean expression across species

for anno_term in ["Neighborhood", "Class", "Subclass", "Group"]:
    ## Identify files for the given annotation term
    expr_files = glob.glob(os.path.join(data_dir, anno_term.lower(), f"*_{anno_term.lower()}_mean_expr_orthologs.h5ad"))
    ## Load expression data from hdf
    expr = {}; genes={}; celltypes = {}
    for file in expr_files:
        ##
        adata = ad.read_h5ad(file)
        ##
        species_name = re.search(r'/([^/]+)_'+anno_term.lower()+'_mean_expr_orthologs\.h5ad$', file).group(1)
        expr[species_name] = pd.DataFrame(adata.layers["mean"], 
                                            index=adata.obs_names, 
                                            columns=adata.var_names)
        genes[species_name] = list(adata.var_names)
        celltypes[species_name] = list(adata.obs_names)
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
                expr[species_name].loc[ct,:] = 0
        expr[species_name] = expr[species_name].loc[common_celltypes,:]
    ## -------------------------------------------
    ## Expressolog data computation
    expressolog = []
    gene_universe = common_genes
    species = list(expr.keys())
    ## Compute expressologs
    for i in range(len(species) - 1):
        sp1 = species[i]
        print(sp1)
        for j in range(i + 1, len(species)):
            print(f"Comparing {sp1} with {species[j]}")
            sp2 = species[j]
            expr_cor = np.empty((len(gene_universe), 3))
            for idx, gene in tqdm(enumerate(gene_universe)):
                ## 
                if gene not in expr[sp1].columns:
                    sp1_expr = np.zeros([len(common_celltypes),1])
                else:
                    sp1_expr = expr[sp1].loc[:,gene].to_numpy().reshape(-1, 1)
                ## 
                if gene not in expr[sp2].columns:
                    sp2_expr = np.zeros([len(common_celltypes),1])
                else:
                    sp2_expr = expr[sp2].loc[:, gene].to_numpy().reshape(-1, 1)
                ##
                sp1_expr_stat = np.std(sp1_expr)
                sp2_expr_stat = np.std(sp2_expr)
                cor1 = np.corrcoef(sp1_expr[:, 0], sp2_expr[:, 0])[0, 1]
                expr_cor[idx, :] = [sp1_expr_stat, sp2_expr_stat, cor1]
            ##
            corm = pd.DataFrame({
                'gene': gene_universe,
                'species1': sp1,
                'species2': sp2,
                'sp1_expr_stat': expr_cor[:, 0],
                'sp2_expr_stat': expr_cor[:, 1],
                'expr_cor': expr_cor[:, 2],
                'sp_order': 1
            })
            ## Symmetric comparison
            corm2 = pd.DataFrame({
                'gene': gene_universe,
                'species1': sp2,
                'species2': sp1,
                'sp1_expr_stat': expr_cor[:, 1],
                'sp2_expr_stat': expr_cor[:, 0],
                'expr_cor': expr_cor[:, 2],
                'sp_order': 2
            })
            ##
            expressolog.append(pd.concat([corm, corm2], ignore_index=True))
    ## Combine all results into a single DataFrame
    expressologs = pd.concat(expressolog, ignore_index=True)
    ## -------------------------------------------
    ## Data saviing and organization
    ## Get species order from metadata
    species_order = pd.Series(["human", "macaque", "marmoset"]).tolist()
    ## Convert species1 and species2 to categorical data types with the given order
    expressologs['species1'] = pd.Categorical(expressologs['species1'], categories=species_order, ordered=True)
    expressologs['species2'] = pd.Categorical(expressologs['species2'], categories=species_order, ordered=True)
    ## Ensure save path exits
    if os.path.exists(os.path.join(wkdir, anno_term.lower())) is False:
        os.mkdir(os.path.join(wkdir, anno_term.lower()))
    ## Write expressologs data.frame to csv
    expressologs.to_csv(os.path.join(wkdir, anno_term.lower(), "expressologs_"+anno_term.lower()+".csv"), index=False)
