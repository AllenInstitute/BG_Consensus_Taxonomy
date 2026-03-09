import os, re, sys, glob
import anndata as ad
import pandas as pd
import keras
import crested

## Function to pad sequences equally on both sides
def pad_sequence(sequence, target_width):
    ## If sequence is to long cut out the middle
    if len(sequence) >= target_width:
        center = len(sequence) // 2
        start = int(center - (target_width / 2))
        end = int(center + (target_width / 2))
        return sequence[start:end]
    ## Otherwise pad the sequence with N's
    total_padding = target_width - len(sequence)
    left_padding = total_padding // 2
    right_padding = total_padding - left_padding
    return ('N' * left_padding) + sequence + ('N' * right_padding)

## Load enhancer library
# enhancer_library = pd.read_excel("/home/nelson.johansen/Analysis/HMBA_Genomics/Analysis/SeqModel/enhancer_tables/striatal_tools_supp_table1.xlsx")
enhancer_library = pd.read_excel("/home/nelson.johansen/Analysis/HMBA_Genomics/Analysis/SeqModel/enhancer_tables/VGT_master_tracker.xlsx", sheet_name="Enhancers")

enhancer_library = enhancer_library.loc[enhancer_library["Enhancer_ID"] == "AiE1519m", :].reset_index(drop=True)

## Filter for enhancers
enhancer_library["sequence"] = enhancer_library["Enhancer_sequence"]
enhancer_library['padded_sequence'] = enhancer_library["sequence"].apply(lambda seq: pad_sequence(seq, 2114)) ## Pad each sequence or cut out center to fit in 2114

# ## Filter to specific enhancer
# enhancer_library = enhancer_library.loc[enhancer_library["Target_brain_region"] == "CTX"]
# enhancer_library.reset_index(inplace=True, drop=True)

## --------------------------
## Load model
species_name = "human"
work_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling"

## -------------------------
## Setup genome
reference = os.path.join(work_dir, f"genomes/{species_name}")
if species_name == "macaque":
    reference = os.path.join(reference, "ncbi")

## Chromosome sizes
if os.path.exists(os.path.abspath(os.path.join(reference, f"{species_name}.chrom.sizes"))):
    chr_sizes = os.path.abspath(os.path.join(reference, f"{species_name}.chrom.sizes"))
else:
    chr_sizes = None

## Fasta and GTF paths
fasta_path = os.path.abspath(os.path.join(reference, f"fasta/genome.fa"))
gtf_path = os.path.abspath(os.path.join(reference, f"genes/genes.gtf.gz"))

##
genome = crested.Genome(
    fasta=fasta_path,
    annotation=gtf_path,
)
crested.register_genome(
    genome
)

## ---------------------------
## Load the pretrained CREsted model

model_name = "CrestedModelAPI--dilated-cnn-human-basalganglia--25-06-05-13-57-57"

## Check if model is finetuned or not.
model_dir = os.path.join(work_dir, "basal-ganglia", "models", f"crested/evogen-dnaseq-modeling/basal_ganglia/{model_name}")
if os.path.exists(model_dir + "_finetune"):
    model_dir = model_dir + "_finetune/checkpoints/"
else:
    model_dir = model_dir + "/checkpoints/"

## Find the model with highest training number this isn't IDEAL as the best model may not be the last one
best_model = max(glob.glob(os.path.join(model_dir, "*.keras")), key=lambda x: int(re.search(r'(\d+)\.keras$', x).group(1)))
print("Loading best model: ", best_model)

## Initalize model
model = keras.models.load_model(best_model, compile=False)

## ---------------------------
## Gather data model was trained on, really only needed to determine class (celltype) names
adata = ad.read_h5ad(os.path.join(work_dir, "basal-ganglia", "data", species_name, "crested_adata", f"{species_name}_basalganglia_hmba_pre-print_crested.h5ad"))

## Run predictions
seqs = enhancer_library.padded_sequence.to_list()
predictions = crested.tl.predict(seqs, model, batch_size=512, verbose=2)

prediction_result = pd.DataFrame(
    index=enhancer_library["Enhancer_ID"],
    columns=adata.obs_names,
    data=predictions  # or np.nan
)

## ---------------------------
## Setup output directories
figure_dir = os.path.join("/home/nelson.johansen/", "striatal_tools", "figures", f"{species_name}")
if not os.path.exists(figure_dir):
    os.makedirs(figure_dir)

## ---------------------------
## Run enhancer sequences through model
for enhancer_idx in range(len(enhancer_library)):
    ## Check if figure exists already
    figure_path = os.path.join(figure_dir, f"{species_name}_{enhancer_library['Enhancer_ID'][enhancer_idx]}_crested_all.png")
    if os.path.exists(figure_path):
        print(f"Figure for enhancer {enhancer_library['Enhancer_ID'][enhancer_idx]} already exists, skipping...")
        continue
    ##
    print(f"Scoring enhancer {enhancer_library['enhancer ID'][enhancer_idx]}")
    ## Calculate contribution scores
    scores, one_hot_encoded_sequences = crested.tl.contribution_scores(
        input=enhancer_library.padded_sequence[enhancer_idx],  # Pass a list of sequences
        target_idx = None,
        model = model,
        method = "integrated_grad",
    )
    ##
    crested.pl.patterns.contribution_scores(
        scores,
        one_hot_encoded_sequences,
        sequence_labels=[enhancer_library["enhancer ID"][enhancer_idx]],
        class_labels=adata.obs_names.to_list(),
        zoom_n_bases=len(enhancer_library["sequence"][enhancer_idx]), #(len(enhancer_library.iloc[enhancer_idx]["sequence"])//2), ## Sometimes we may want to zoom in on a part of the enhancer
        title=f"Species: {species_name} - Enhancer: {enhancer_library['enhancer ID'][enhancer_idx]}",
        save_path=f"{figure_dir}/{species_name}_{enhancer_library['enhancer ID'][enhancer_idx]}_crested_all.png",
    )

## ---------------------------
## Calculate contribution scores
scores, one_hot_encoded_sequences = crested.tl.contribution_scores(
    input=enhancer_library.padded_sequence.tolist(),  # Pass a list of sequences
    target_idx = None,
    model = model,
    method = "integrated_grad",
)

contribution_result = pd.DataFrame(
    index=enhancer_library["enhancer ID"],
    columns=adata.obs_names,
    data=0  # or np.nan
)

for enhancer_idx in range(enhancer_library.shape[0]):
    print(f"Calculating contribution scores for enhancer {enhancer_library['enhancer ID'][enhancer_idx]}")
    ## Figure zoom on enhancer
    zoom_n_bases=(len(enhancer_library.iloc[enhancer_idx]["sequence"])//2)
    center = int(scores.shape[2] / 2)
    start_idx = center - int(zoom_n_bases / 2)
    ## Filter one-hot
    seq_class_x = one_hot_encoded_sequences[0, start_idx : start_idx + zoom_n_bases,:]
    for class_idx in range(adata.obs_names.shape[0]):
        ## Gather scores from enhancer element in sequence
        seq_class_scores = scores[enhancer_idx, class_idx, start_idx : start_idx + zoom_n_bases, :]
        ## Create dataframe for integrated gradients
        intgrad_df = crested.pl.patterns._utils.grad_times_input_to_df(seq_class_x, seq_class_scores)
        intgrad_df['contribution_score'] = intgrad_df.max(axis=1)
        summary_stat = intgrad_df['contribution_score'].nlargest(10).median()
        ## Add to contribution result
        contribution_result.loc[enhancer_library["enhancer ID"][enhancer_idx], adata.obs_names[class_idx]] = summary_stat


## Plot heatmap of contribution scores
def plot_heatmap(data, title, normalization=None, save_path="."):
    """
    Plot a heatmap of the contribution scores.
    
    Parameters:
    - data: DataFrame with contribution scores.
    - title: Title for the heatmap.
    - save_path: Path to save the heatmap image.
    """
    import seaborn as sns
    import matplotlib.pyplot as plt
    ## Normalize data if specified
    if normalization == "min-max":
        ## Row-wise normalization, so each row is min-max normalized
        min_vals = data.min(axis=1)
        max_vals = data.max(axis=1)
        range_vals = max_vals - min_vals
        range_vals[range_vals == 0] = 1  # avoid division by zero
        data = (data.sub(min_vals, axis=0)).div(range_vals, axis=0)
    elif normalization == "probability":
        data = data.div(data.sum(axis=1) + 1e-12, axis=0)
    else:
        print("No normalization applied, using raw scores.")
    ## Figure setup
    plt.figure(figsize=(12, 10))
    # Plot clustered heatmap
    heat = sns.clustermap(
        data,
        cmap='viridis',
        method='average',      # linkage method for clustering
        metric='euclidean',    # distance metric
        col_cluster=False,
        row_cluster=True,      # cluster rows (enhancers)
        figsize=(12, 12),
        standard_scale=None,
        xticklabels=True,
        yticklabels=True,
    )
    plt.title(title)
    plt.setp(heat.ax_heatmap.get_xticklabels(), rotation=90)
    plt.setp(heat.ax_heatmap.get_yticklabels(), rotation=0)
    plt.xlabel('Cell Types')
    plt.ylabel('Enhancers')
    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    plt.close()

## Order data following consensus annotation
anno_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/AnnoTables"
consensus_anno = pd.read_csv(os.path.join(anno_dir, "consensus_anno_pre-print.csv"), index_col=0)
consensus_anno["Group"] = consensus_anno["Group"].str.replace(" ", "_").unique()
consensus_anno = consensus_anno.loc[consensus_anno["Group"].isin(prediction_result.columns.unique()), :]

prediction_result = prediction_result[consensus_anno["Group"].values]

plot_heatmap(
    prediction_result,
    normalization="min-max",
    title=f"{species_name} Enhancer Contribution Scores",
    save_path=os.path.join(figure_dir, f"{species_name}_enhancer_prediction_min-max_heatmap.png")
)
