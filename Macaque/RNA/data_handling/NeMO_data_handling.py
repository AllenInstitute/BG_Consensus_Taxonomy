import anndata as ad
import pandas as pd
import numpy as np
from tqdm import tqdm
import glob
import os

##
species = "Macaque"

##
data_dir = f"/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/{species}"
work_dir = f"/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/{species}/BasalGanglia"

##
libraries = glob.glob(os.path.join(data_dir, "NeMO/GEX/**/*.h5ad"), recursive=True)

##
adatas = []
with open(f"errorLog{species}_NeMO_data_merge.txt", 'w') as file:
    for library in tqdm(libraries):
        try:
            tmp = ad.read_h5ad(library)
            mean_expression = np.array(tmp.X.mean(axis=0)).flatten()
            gene_names = tmp.var_names 
            mean_expression_df = pd.DataFrame({'Gene': gene_names, 
                                                    'Mean Expression': mean_expression})
            print(mean_expression_df.loc[mean_expression_df.Gene == "KCNIP1",:])
            tmp = tmp[tmp.obs["star_IsCell"] == True]
            adatas.append(tmp)
        except Exception as e:
            print("Library: ", os.path.basename(library))
            print("An error occurred:", e)
            file.write("An error occurred for library: " + os.path.basename(library) + "\n")
            file.write(str(e) + "\n")
            file.flush()

##
adata = ad.concat(adatas, merge="same", uns_merge="same")

#######################################
## Adjust meta.data

## Rename
adata.obs.rename(columns={"CellID": "cell_barcode"}, inplace=True)
del adata.obs["cell_names"] ## This is 1 of 2 spots that RNA barcode is recorded
del adata.obs["barcoded_cell_sample"]

## Pull barcoded_cell_sample and alignment_job_id
adata.obs["barcoded_cell_sample_label"] = adata.obs["input_id"].str.split("_Job", expand=True)[0]
adata.obs["alignment_job_id"] = "Job" + adata.obs["input_id"].str.replace("_gex", "").str.split("_Job", expand=True)[1]
adata.obs["alignment_job_database"] = "NeMO"

if species is 'Human':
    adata.obs["barcoded_cell_sample_label"] = adata.obs["barcoded_cell_sample_label"].str.upper()

## Add cell_label and update index to match
adata.obs["cell_label"] = adata.obs["cell_barcode"].astype(str) + "-" + adata.obs["barcoded_cell_sample_label"].astype(str)
adata.obs.index = adata.obs["cell_label"].copy()

#######################################
## Base file
adata.write(os.path.join(work_dir, f"{species}_basalganglia_NeMO.h5ad"))

#######################################
## Add metadata from Kim's tracker and adjust to AIT standards now
tracker = pd.read_csv("/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BICAN_tracker_NIMP-NeMO_240909.csv")
tracker = tracker.loc[:,["library_avg_size_bp", "library_concentration_nm", "library_input_ng",
                        "library_prep_pass_fail", "library_prep_set", "library_quantification_fmol",
                        "library_quantification_ng", "barcoded_cell_sample_local_name", "barcoded_cell_sample_tag_local_name", 
                        "barcoded_cell_sample_number_of_expected_cells", "barcoded_cell_sample_technique",
                        "barcoded_cell_input_quantity_count", "tissue_sample_local_name", "tissue_structure",	
                        "species", "project_identifier", "donor_local_name", "age_at_death_value", "age_at_death_unit",
                        "biological_sex", "Species common", "ROI abbrv", "NeMO Bucket"]]
tracker.rename(columns={"donor_local_name": "donor_id",
                        "Species common": "organism",
                        "species": "organism_ontology_term_id",
                        "age_at_death_value": "donor_age",
                        "ROI abbrv": "anatomical_region",
                        "tissue_structure": "brain_region_ontology_term_id",
                        "biological_sex": "self_reported_sex",
                        "barcoded_cell_sample_technique": "assay",
                        "NeMO Bucket": "nemo_bucket"
                        }, 
                inplace=True)

tracker.drop_duplicates(subset=["barcoded_cell_sample_local_name"], inplace=True)
tracker = tracker.loc[tracker["organism"] == species]
tracker = tracker.loc[tracker["barcoded_cell_sample_local_name"].isin(adata.obs["barcoded_cell_sample_label"].unique())]

## Merge metadata
metadata = adata.obs.merge(tracker, left_on="barcoded_cell_sample_label", right_on="barcoded_cell_sample_local_name")
metadata.index = adata.obs["cell_label"].copy()

##
np.all(adata.obs.index == metadata.index)

## Update metadata
adata.obs = metadata.copy()

#######################################
## Save NeMO_BICAN file
adata.write(os.path.join(work_dir, f"{species}_basalganglia_NeMO_BICAN.h5ad"))
