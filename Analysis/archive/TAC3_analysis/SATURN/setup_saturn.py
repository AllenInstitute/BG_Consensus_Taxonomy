import os, glob
import pandas as pd
import anndata as ad

##
data_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis/samap"
emb_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis/saturn/proteome/embeddings"
analysis_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis/saturn/TAC3"

## --------------------------------------------
## Create dataframe for SATURN training
## --------------------------------------------
df = pd.DataFrame(columns=["path", "species", "embedding_path"])
df["species"] = ["human", "macaque", "marmoset", "mouse"]
df["path"] = [os.path.join(analysis_dir, "HMBA_Human_for_saturn.h5ad"),
                os.path.join(analysis_dir, "HMBA_Macaque_for_saturn.h5ad"),
                os.path.join(analysis_dir, "HMBA_Marmoset_for_saturn.h5ad"),
                os.path.join(analysis_dir, "Mouse_Broad_for_saturn.h5ad")]
df["embedding_path"] = [os.path.join(emb_dir, "Homo_sapiens.GRCh38.gene_symbol_to_embedding_ESM1b.pt"),
                        os.path.join(emb_dir, "Macaca_mulatta.Mmul_10.gene_symbol_to_embedding_ESM1b.pt"),
                        os.path.join(emb_dir, "Callithrix_jacchus.mCalJac1.pat.X.gene_symbol_to_embedding_ESM1b.pt"),
                        os.path.join(emb_dir, "Mus_musculus.GRCm39.gene_symbol_to_embedding_ESM1b.pt")]

df.to_csv(os.path.join(analysis_dir, "TAC3_run.csv"), index=False)

## --------------------------------------------
## Create cell type dataframe for SATURN training
## --------------------------------------------
## Preprocess data for SAMap

h5ad_files = glob.glob(f"{data_dir}/*for_samap.h5ad")
adatas = []
for filename in h5ad_files:
    adata = ad.read_h5ad(filename)
    if "Group" in adata.obs.columns:
        adata.obs["cell_type"] = adata.obs.Group
    else:
        adata.obs["cell_type"] = adata.obs["AIT21.cluster"]
    adata.write_h5ad(os.path.join(analysis_dir, os.path.basename(filename).replace("for_samap", "for_saturn")))

##
combined_orig = ad.concat(adatas, label="species", keys=["hu", "mq", "ms", "mm"], index_unique=None, join="outer", merge="unique")
combined_obs = combined_orig.obs

python ./SATURN/train-saturn.py --in_data=./TAC3/TAC3_run.csv \
                              --in_label_col=cell_type --ref_label_col=cell_type \
                              --num_macrogenes=2000     --hv_genes=8000          \
                              --work_dir=./TAC3/v2 \
                              --device_num=1 \
                              --epochs=200