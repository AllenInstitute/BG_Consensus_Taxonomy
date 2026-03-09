import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial import cKDTree
from scipy import sparse

def compute_ABC_scores(
    peaks_df,
    tss_df,
    peak_matrix,        # n_peaks x n_celltypes matrix (numpy or sparse)
    celltypes=None,          # list of celltype names
    method='exp',            # 'exp' or 'power'
    scale=10000.0,           # for exponential decay
    power=1.0,               # for power law
    max_distance=500000,     # maximum bp distance
    normalize_per_gene=True  # normalize per gene per celltype
):
    """
    Compute ABC scores for each peak-gene-celltype triplet.

    Returns
    -------
    abc_df : DataFrame
        columns: ['gene_id', 'chrom', 'peak_start', 'peak_end',
                  'distance', 'celltype', 'activity', 'contact', 'ABC']
    """
    # Convert accessibility matrix to dense if needed
    if sparse.issparse(peak_matrix):
        peak_matrix = peak_matrix.toarray()
    n_peaks, n_celltypes = peak_matrix.shape
    if celltypes is None:
        celltypes = [f"celltype_{i}" for i in range(n_celltypes)]
    # Prepare coordinates
    peak_chroms = peaks_df['chrom'].astype(str).values
    peak_mids = ((peaks_df['start'].values + peaks_df['end'].values) // 2).astype(np.int64)
    tss_chroms = tss_df['chrom'].astype(str).values
    tss_pos = tss_df['tss'].astype(np.int64).values
    gene_ids = tss_df.index.to_numpy()
    ##
    rows = []
    # Precompute distance-based contact weights
    for chrom in tqdm(sorted(set(peak_chroms) & set(tss_chroms)), desc="Chromosomes"):
        p_idx = np.where(peak_chroms == chrom)[0]
        g_idx = np.where(tss_chroms == chrom)[0]
        if len(p_idx) == 0 or len(g_idx) == 0:
            continue
        ##
        tree = cKDTree(peak_mids[p_idx].reshape(-1, 1))
        neighbors = tree.query_ball_point(tss_pos[g_idx].reshape(-1, 1), r=max_distance)
        ##
        for gi, neigh in enumerate(neighbors):
            if len(neigh) == 0:
                continue
            gene_i = g_idx[gi]
            gene_id = gene_ids[gene_i]
            gene_tss = tss_pos[gene_i]
            ##
            local_peaks = p_idx[neigh]
            dists = np.abs(peak_mids[local_peaks] - gene_tss).astype(float)
            # Compute contact
            if method == 'exp':
                contact = np.exp(-dists / scale)
            elif method == 'power':
                contact = 1.0 / (np.power(dists, power) + 1e-9)
            else:
                raise ValueError("method must be 'exp' or 'power'")
            # Get activities per cell type
            activity = peak_matrix[local_peaks, :]  # shape (n_local_peaks, n_celltypes)
            # Multiply Activity × Contact for each cell type
            abc_raw = activity * contact[:, np.newaxis]  # broadcast contact
            ##
            if normalize_per_gene:
                denom = abc_raw.sum(axis=0, keepdims=True)  # per celltype
                denom[denom == 0] = 1
                abc_norm = abc_raw / denom
            else:
                abc_norm = abc_raw
            #
            n_local = len(local_peaks)
            # Build correctly aligned DataFrame
            df = pd.DataFrame({
                'gene_id': np.repeat(gene_id, n_local * n_celltypes),
                'chrom': np.repeat(chrom, n_local * n_celltypes),
                'peak_start': np.repeat(peaks_df.iloc[local_peaks]['start'].values, n_celltypes),
                'peak_end':   np.repeat(peaks_df.iloc[local_peaks]['end'].values, n_celltypes),
                'distance':   np.repeat(dists, n_celltypes),
                'celltype':   np.tile(celltypes, n_local),
                'activity':   activity.flatten(order='C'),    # row-major flatten (peak-major)
                'contact':    np.repeat(contact, n_celltypes),
                'ABC':        abc_norm.flatten(order='C')
            })
            ## Sanity check
            rows.append(df)
    ##
    abc_df = pd.concat(rows, ignore_index=True)
    ##Remove 0 activity rows
    abc_df = abc_df[abc_df['activity'] > 0].reset_index(drop=True)
    return abc_df


## --------------------------------------------------
## Compute ABC scores for spinal cord ATAC data
## --------------------------------------------------
import os
import anndata as ad
import pandas as pd
from tqdm import tqdm

analysis_dir =  "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis"

## --------------------------------------------------
## Load TSS data
## --------------------------------------------------
## Species genome annotation (GTF)
genome_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/anno_tables/ATAC/gencode"
species_gtf_files = {
    species: f"{genome_dir}/{species}_gencode.gtf.gz"
    for species in ["human", "macaque", "marmoset", "mouse"]
}

base_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/basal-ganglia/data/"
species_h5ads = {
    "human" : os.path.join(base_dir, "human", "crested_adata", "human_basalganglia_hmba_pre-print_crested_500bp.h5ad"),
    "macaque" : os.path.join(base_dir, "macaque", "crested_adata", "macaque_basalganglia_hmba_pre-print_crested_500bp.h5ad"),
    "marmoset": os.path.join(base_dir, "marmoset", "crested_adata", "marmoset_basalganglia_hmba_pre-print_crested_500bp.h5ad"),
}

## Process each species
for species in ["human", "macaque", "marmoset"]:
    ## --------------------------------------------------
    ## Load in accessibility data
    ## --------------------------------------------------
    adata = ad.read_h5ad(species_h5ads[species])
    ## Create peak dataframe
    peaks_df = pd.DataFrame(index=adata.var.index)
    peaks_df["chrom"] = adata.var["chr"].values
    peaks_df["start"] = adata.var["start"].values
    peaks_df["end"] = adata.var["end"].values
    ## --------------------------------------------------
    ## Load TSS data
    ## --------------------------------------------------
    ## Function to parse GTF attributes
    def parse_attributes(attributes_str):
        """Parse GTF attributes field into a dictionary."""
        attributes = {}
        for attr in attributes_str.strip().split(";"):
            if attr.strip():
                key_value = attr.strip().split(" ")
                if len(key_value) == 2:
                    key = key_value[0]
                    value = key_value[1].strip('"')
                    attributes[key] = value
        return attributes
    ## Load GTF files
    gtf = pd.read_csv(species_gtf_files[species], sep="\t", header=None, comment="#",
                    names=["Chromosome", "Source", "Feature", "Start", "End",
                            "Score", "Strand", "Frame", "Attributes"], compression='gzip')
    ## Keep only gene features
    transcripts = gtf[gtf["Feature"] == "gene"].copy()
    ## Parse attributes
    transcripts["Attributes"] = transcripts["Attributes"].apply(parse_attributes)
    transcripts["gene_id"] = transcripts["Attributes"].apply(lambda x: x.get("gene_id", None))
    transcripts["gene_name"] = transcripts["Attributes"].apply(lambda x: x.get("gene_name", None))
    transcripts["gene_type"] = transcripts["Attributes"].apply(lambda x: x.get("gene_type", None))
    ## Keep protein-coding genes
    if transcripts["gene_type"].str.contains("protein_coding").any():
        transcripts = transcripts[transcripts["gene_type"] == "protein_coding"].copy()
    ## Strand-aware TSS
    tss_df = transcripts.copy()
    tss_df["tss"] = tss_df.apply(lambda r: r["Start"] if r["Strand"] == "+" else r["End"], axis=1)
    tss_df["chrom"] = tss_df["Chromosome"]
    tss_df.set_index(tss_df.gene_name, inplace=True)
    ## --------------------------------------------------
    ## Compute ABC
    ## --------------------------------------------------
    ABC_scores = compute_ABC_scores(
        peaks_df,
        tss_df,
        sparse.csr_matrix(adata.X.T),              # scipy.sparse matrix (n_peaks x n_cells)
        celltypes=adata.obs_names.tolist(),  # list of celltype names
        method='power',             # 'exp' or 'power'
        scale=10000.0,            # for 'exp': scale in bp (decay distance); for 'power' not used
        power=0.75,                # for 'power': weight = 1/(d**power + eps)
        max_distance=500000,      # ignore peak-gene pairs beyond this (in bp)
    )
    ABC_scores.to_csv(os.path.join(analysis_dir, "atac", "abc-model", f"ABC_scores_{species}.csv"), index=False)

ABC_scores_f = ABC_scores.loc[ABC_scores.ABC > 0.02]
ABC_scores_f.to_csv(os.path.join(analysis_dir, "atac", "abc-model", f"ABC_scores_{species}_filtered_ABC_gt_0.02.csv"), index=False)