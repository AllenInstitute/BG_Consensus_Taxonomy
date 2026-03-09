import anndata as ad
import os
import gc
import nsforest as ns
from nsforest import preprocessing as pp
from nsforest import nsforesting

species_list = ["Macaque", "Marmoset", "Human"]
anno = "cluster_id"

base_input_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/hct_ux3_cellxgene/anndata_080/BasalGanglia/BICAN-releases/final/downsample/"  # same as downsample output
base_result_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/Human/BasalGanglia/analysis/nsforest/"

for species in species_list:
    print(f"Running NS-Forest on {species}...")

    # Load downsampled data
    input_file = os.path.join(base_input_dir, f"{species}_downsampled_{anno}_10k.h5ad")
    adata = ad.read_h5ad(input_file)

    # Compute medians and scores
    adata = pp.prep_medians(adata, anno, use_mean=False, positive_genes_only=True)
    adata = pp.prep_binary_scores(adata, anno, medians_header="medians_")

    # Output path
    output_path = os.path.join(base_result_dir, species, anno)
    os.makedirs(output_path, exist_ok=True)

    # Run NS-Forest
    results = nsforesting.NSForest(
        adata,
        cluster_header=anno,
        save_supplementary=True,
        output_folder=output_path,
        outputfilename_prefix=f"{species}_BG_{anno}_NSForest"
    )

    del adata
    del results
    gc.collect()

    print(f"Completed NS-Forest for {species}")

