import scanpy as sc
import anndata as ad
from matplotlib import rc_context

adata_GPi_shell = adata[adata.obs.loc[adata.obs.Group.isin(["GPi Shell"]),:].index]

markers = {"GPi Shell": ["SLC17A6", "TBR1", "SST"],
           "GPi Core": ["PVALB"],
           }

with rc_context(rc={'figure.dpi': 600}):    
    sc.pl.dotplot(adata_GPi_shell, markers, "organism", save="_GPi_shell_dotplot.png") 


adata_GPi_core = adata[adata.obs.loc[adata.obs.Group.isin(["GPi Core"]),:].index]

markers = {"GPi Shell": ["SLC17A6", "TBR1", "SST"],
           "GPi Core": ["PVALB"],
           }

with rc_context(rc={'figure.dpi': 600}):    
    sc.pl.dotplot(adata_GPi_core, markers, "organism", save="_GPi_core_dotplot.png") 