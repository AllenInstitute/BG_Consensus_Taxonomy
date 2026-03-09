import os, sys
import anndata as ad
import scanpy as sc

species = "xspecies"
work_dir = f"/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/{species}/BasalGanglia"

## Gather integrated object with ortholog genes
adata = ad.read_h5ad(os.path.join(work_dir, "HMBA_Human_Macaque_Marmoset_Mouse_snRNA-seq_BG_alignment_v6.h5ad"))

## Subset to TAC3
TAC3_idx = (adata.obs["AIT21.subclass"].isin(["055 STR Lhx8 Gaba"])) | (adata.obs.Group.isin(["STR-BF TAC3-PLPP4-LHX8 GABA", "STR TAC3-PLPP4 GABA"]))
adata = adata[TAC3_idx]

## Basics
adata.X = adata.layers["UMIs"].copy()
sc.pp.normalize_total(adata, target_sum=1e6)
sc.pp.log1p(adata)

## Save
adata.write_h5ad("/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis/TAC3/TAC3_data_merged_ortholog.h5ad")