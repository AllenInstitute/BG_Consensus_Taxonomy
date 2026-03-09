import xarray as xr
import pandas as pd

de_results = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/Analysis/Markers/log_reg_markers.netcdf"
# de_results = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/Analysis/DE/deseq2_spc.netcdf"

ds = xr.open_dataset(de_results)
# ds['lfc'].sel({'contrast':'Human-Macaque'}).to_dataframe().reset_index()

# ## Pull lfc and pval for Human-Macaque creating a data.frame with celltype, species1, species2, lfc, pval
# human_macaque_df = ds['lfc'].sel({'contrast':'Human-Macaque'}).to_dataframe().reset_index()
# human_macaque_pval_df = ds['padj'].sel({'contrast':'Human-Macaque'}).to_dataframe().reset_index()
# human_macaque_df = human_macaque_df.merge(human_macaque_pval_df, on=['gene', 'Group', 'contrast'], suffixes=('_lfc', '_padj'))

# ## Split contrast into species1 and species2
# human_macaque_df[['species1', 'species2']] = human_macaque_df['contrast'].str.split('-', expand=True)
# human_macaque_df = human_macaque_df[['Group', 'gene', 'species1', 'species2', 'lfc', 'padj']]

# ## Create a dataframe with significant DE genes
# sig_de_df = human_macaque_df[
#     (human_macaque_df['padj'] < 0.01) & 
#     (human_macaque_df['lfc'] > 4)
# ].copy()

##