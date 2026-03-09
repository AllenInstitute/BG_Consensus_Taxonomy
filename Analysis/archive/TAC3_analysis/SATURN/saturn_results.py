import os
import scanpy as sc
import pandas as pd
import pickle
import matplotlib.pyplot as plt

##
work_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis/saturn/TAC3"
cxg_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/hct_ux3_cellxgene/anndata_080/BasalGanglia/xspecies"

##
adata = sc.read_h5ad(f"{work_dir}/v2/saturn_results/test256_data_HMBA_Human_for_saturn_HMBA_Macaque_for_saturn_HMBA_Marmoset_for_saturn_Mouse_Broad_for_saturn_org_saturn_seed_0.h5ad")

## Macrogene weights
# with open(os.path.join(work_dir, "v2/saturn_results/test256_data_HMBA_Human_for_saturn_HMBA_Macaque_for_saturn_HMBA_Marmoset_for_saturn_Mouse_Broad_for_saturn_org_saturn_seed_0_genes_to_macrogenes.pkl"), "rb") as f:
#     macrogene_weights = pickle.load(f)

# macrogene_adata = sc.AnnData(adata.obsm["macrogenes"])
# macrogene_adata.obs = adata.obs

# sc.tl.rank_genes_groups(macrogene_adata, groupby="labels", method="wilcoxon")
# sc.pl.rank_genes_groups_dotplot(macrogene_adata, n_genes=5, figsize=(15,10))

# plt.savefig(os.path.join(work_dir, "v2/saturn_results/TAC3_saturn_macrogene_dotplot_human_macaque_marmoset_mouse_broad.png"), dpi=300)

# def get_scores(macrogene):
#     '''
#     Given the index of a macrogene, return the scores by gene for that centroid
#     '''
#     scores = {}
#     for (gene), score in macrogene_weights.items():
#         scores[gene] = score[int(macrogene)]
#     return scores

# ## Macrogenes space is a joint space composed of genes inferred to be functionally related based on the similarity of their protein embeddings.
# for macrogene in ["40"]:
#     print(f"Macrogene {macrogene}")
#     pd.DataFrame(get_scores(macrogene).items(), columns=["gene", "weight"])\
#             .sort_values("weight", ascending=False)\
#             .head(25)

# pd.DataFrame(macrogene_weights["mouse_Tac2"].tolist(), columns=["weight"]).sort_values("weight", ascending=False).head(25)

sc.pp.neighbors(adata, use_rep="X")
sc.tl.umap(adata)

## Write to cxg
adata.write(os.path.join(cxg_dir, "TAC3_saturn_results_human_macaque_marmoset_mouse_broad.h5ad"))
