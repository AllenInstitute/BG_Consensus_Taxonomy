
import anndata
import crested
import keras
import pysam
from tqdm import tqdm
import argparse
import os

def main(
    adata_path,
    genome_file,
    chrom_sizes,
    model_path,
    output_dir,
    top_k,
    manual_batching,
    manual_batch_size,
    gradient_method,
    filtering_method,
):
    # Load data
    adata = anndata.read_h5ad(adata_path)

    # Load and register genome
    genome = crested.Genome(genome_file, chrom_sizes)
    crested.register_genome(genome)

    # Load model
    model = keras.models.load_model(model_path, compile=False)

    # Fetch sequences
    fasta = pysam.FastaFile(genome_file)
    seqs = [
        fasta.fetch(chrom, start, end)
        for chrom, start, end in tqdm(zip(adata.var['chr'], adata.var['start'], adata.var['end']))
    ]

    # Predict and store
    predictions = crested.tl.predict(
        seqs,
        model
    )
    adata.layers["finetuned"] = predictions.T

    # Average with original data
    adata_combined = adata.copy()
    if filtering_method == 'combined':
        adata_combined.X = (adata_combined.X + adata_combined.layers["finetuned"]) / 2
    elif filtering_method == 'prediction':
        adata_combined.X = adata_combined.layers["finetuned"]
    elif filtering_method == 'atac':
        adata_combined.X = adata_combined.X
    else:
        raise ValueError("Invalid sorting method. Choose one of: 'combined', 'prediction', or 'atac'.")

    # Filter for top_k informative regions
    adata_filtered = adata_combined.copy()
    crested.pp.sort_and_filter_regions_on_specificity(
        adata_filtered, top_k=top_k, method="proportion"
    )

    # Compute contribution scores
    crested.tl.contribution_scores_specific(
        input=adata_filtered,
        target_idx=None,
        model=model,
        output_dir=output_dir,
        method=gradient_method,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run Crested contribution analysis pipeline.")
    parser.add_argument("--adata_path", type=str, required=True)
    parser.add_argument("--genome_file", type=str, required=True)
    parser.add_argument("--chrom_sizes", type=str, required=True)
    parser.add_argument("--model_path", type=str, required=True)
    parser.add_argument("--output_dir", type=str, required=True)
    parser.add_argument("--top_k", type=int, default=500)
    parser.add_argument("--manual_batching", action="store_true")
    parser.add_argument("--manual_batch_size", type=int, default=256000)
    parser.add_argument("--gradient_method", type=str, default='expected_integrated_grad')
    parser.add_argument("--filtering_method", type=str, default='combined')

    args = parser.parse_args()
    main(
        args.adata_path,
        args.genome_file,
        args.chrom_sizes,
        args.model_path,
        args.output_dir,
        args.top_k,
        args.manual_batching,
        args.manual_batch_size,
        args.gradient_method,
        args.filtering_method
    )
