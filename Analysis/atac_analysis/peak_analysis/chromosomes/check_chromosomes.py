from upsetplot import from_indicators, UpSet
import matplotlib.pyplot as plt
import pandas as pd
import os

os.chdir('/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/ATAC/human/Consensus_peaks/')

## Set custom resolution and size
plt.rcParams['figure.dpi'] = 300  # Resolution in dots per inch
plt.rcParams['figure.figsize'] = [10, 6]  # Width, height in inches

## Example binary matrix
peak_data = pd.read_csv("Group_by_peaks.csv", index_col=0) 

##
chromosomes = list(peak_data.index.str.extract(r'^(.*?):')[0].unique())

##
for chrom in chromosomes:
    print(f"Processing chromosome: {chrom}")
    ## Filter data for the current chromosome
    chrom_adata = peak_data[peak_data.index.str.startswith(chrom)]
    if chrom_adata.shape[0] <= 1:
        print(f"No data found for chromosome {chrom}. Skipping...")
        continue
    ## Bar plot of number of peaks per group
    group_counts = chrom_adata.sum(axis=0).sort_values(ascending=False)
    plt.figure(figsize=(10, 6))
    group_counts.plot(kind='bar', color='skyblue')
    plt.title(f'Number of Peaks per Group on {chrom}')
    plt.xlabel('Group')
    plt.ylabel('Number of Peaks')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    ## Save to PNG
    plt.savefig(f'/home/nelson.johansen/upset_plot_{chrom}_groups.png', dpi=300, bbox_inches='tight')
    plt.close()
    ## Bar plot of the number of groups each peak belongs to
    peak_group_counts = chrom_adata.sum(axis=1).value_counts().sort_index()
    plt.figure(figsize=(10, 6))
    peak_group_counts.plot(kind='bar', color='salmon')
    plt.title(f'Number of Groups per Peak on {chrom}')
    plt.xlabel('Number of Groups')
    plt.ylabel('Number of Peaks')
    plt.xticks(rotation=0)
    plt.tight_layout()
    ## Save to PNG
    plt.savefig(f'/home/nelson.johansen/upset_plot_{chrom}_peaks.png', dpi=300, bbox_inches='tight')
    plt.close()