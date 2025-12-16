import anndata as ad
import pandas as pd
import scipy
import h5py
import os

work_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/changkyul/IntegrativeClustering/WholeBrain/QCed/WMB2/Consensus/ANNData"
hmba_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/data/mouse/"

## ABC atlas metadata
abc_meta = pd.read_csv(os.path.join(hmba_dir, "BG_AIT21_cl.list_scRNAseq.csv"))

## Potential GPi shell info
gpi_shell = pd.read_csv(os.path.join(hmba_dir, "GPi_candidates.Macosko.sample_id.csv"))

## Broad data
files = [
    'WMB_HY-EA-Glut-GABA.Macosko.h5ad',
    'WMB_MB-HB-CB-GABA.Macosko.h5ad',
    'WMB_MB-HB-Glut-Sero-Dopa.Macosko.h5ad',
    'WMB_NN-IMN-GC.Macosko.h5ad',
    'WMB_Pallium-Glut.Macosko.h5ad',
    'WMB_Subpallium-GABA.Macosko.h5ad',
    'WMB_TH-EPI-Glut.Macosko.h5ad',
]

## Load all the obs
obs = []
for file in files:
    with h5py.File(os.path.join(work_dir, file)) as f:
        tmp_obs = ad.io.read_elem(f['obs'])
        tmp_obs["dataset"] = file.split('.')[0]  # Add dataset name
    obs.append(tmp_obs)

##
obs_combined = pd.concat(obs, axis=0)

## Broad metadata
broad_meta = pd.read_csv(os.path.join(hmba_dir, "broad_mouse_meta.csv"), index_col=0)
broad_meta_bg = broad_meta.loc[broad_meta["AIT21.cl"].isin(abc_meta["x"]) | broad_meta.sample_id.isin(gpi_shell.sample_id.tolist()),:]
broad_meta_bg.set_index("sample_id", drop=True, inplace=True)

## Broad metadata to obs
broad_meta_bg = broad_meta_bg.merge(obs_combined, left_index=True, right_index=True, how="left")
broad_meta_bg = broad_meta_bg.loc[~broad_meta_bg.index.duplicated(keep='first'), :]

## Load in all anndata files and filter to broad_meta_bg index on obs_names
for file in files:
    print(f"Processing file: {file}")
    adata = ad.read_h5ad(os.path.join(work_dir, file))
    adata = adata[adata.obs_names.isin(broad_meta_bg.index), :]
    adata.obs = adata.obs.merge(broad_meta_bg, left_index=True, right_index=True, how="left")
    adata.write_h5ad(os.path.join(hmba_dir, f"filtered_{file}"))

## Now load the subset files and combine them
adatas = []
for file in files:
    print(f"Loading filtered file: {file}")
    adata = ad.read_h5ad(os.path.join(hmba_dir, f"filtered_{file}"))
    adatas.append(adata)

## Combine!
combined_adata = ad.concat(adatas, label="dataset", keys=[f.split('.')[0] for f in files], merge="same")
combined_adata.X = scipy.sparse.csr_matrix(combined_adata.X)
combined_adata.write_h5ad(os.path.join(hmba_dir, "Broad_Mouse_BasalGanglia_v2.h5ad"))

## Remove filtered cells
work_dir = f"/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/xspecies/BasalGanglia"
with h5py.File(os.path.join(work_dir, "HMBA_Human_Macaque_Marmoset_Mouse_snRNA-seq_BG_alignment_v3.h5ad")) as f:
    trimmed_obs = ad.io.read_elem(f['obs'])

##
final_adata = combined_adata[combined_adata.obs_names.isin(trimmed_obs.index) | combined_adata.obs_names.isin(gpi_shell.sample_id.tolist()), :]
final_adata.write_h5ad(os.path.join(hmba_dir, "Broad_Mouse_BasalGanglia_pre-print.h5ad"))
