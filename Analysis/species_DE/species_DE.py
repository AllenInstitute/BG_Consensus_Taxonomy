import os, sys
import xarray as xr

data_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/anno_tables/RNA/"

##
da = xr.open_dataset(os.path.join(data_dir, 'deseq2_bg.netcdf'))

##--------------------------------------------
species_de = {}
for contrast in ["Human-Macaque", "Human-Marmoset"]:
    lfc  = da['lfc'].sel({'contrast':contrast}).to_dataframe().reset_index()
    padj = da['padj'].sel({'contrast':contrast}).to_dataframe().reset_index()
    ##
    species_de[contrast] = lfc
    species_de[contrast]['padj'] = padj['padj']
    species_de[contrast] = species_de[contrast].loc[
        (species_de[contrast]['padj'] < 0.001) & (abs(species_de[contrast]['lfc']) > 2)
    ]

## Find human-specific DE genes (significant in both contrasts) for each Group term
human_specific_genes = {}
for group in da.coords["Group"].values:
    human_genes = set(species_de["Human-Macaque"].loc[species_de["Human-Macaque"]['Group'] == group]['gene']).intersection(
        set(species_de["Human-Marmoset"].loc[species_de["Human-Marmoset"]['Group'] == group]['gene'])
    )
    human_specific_genes[group] = human_genes
    print(f"{group}: {len(human_genes)} human-specific DE genes")

## Create DataFrame for human-specific DE genes
human_specific_df = pd.DataFrame(
    [(group, gene) for group, genes in human_specific_genes.items() for gene in genes],
    columns=['Group', 'gene']
)
human_specific_df.to_csv(os.path.join(data_dir, "human_specific_DE_genes_bg.csv"), index=False)

# conserved_da = da.min(dim=species_key, skipna=True) # can add .sel({'organism':['Human','Macaque','Marmoset']}) to do just primate (but looks the sameish)
# conserved_df = pd.DataFrame(
#     conserved_da.transpose("var", leaf_key).values,
#     index=genes,
#     columns=conserved_da.coords[leaf_key].data
# )

# top_genes = None
# if top_k is not None:
#     top_genes = {
#         g: conserved_df[g].dropna().sort_values(ascending=False).head(top_k)
#         for g in conserved_da.coords[leaf_key].data
#     }