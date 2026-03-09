import pandas as pd
import pyranges as pr
import os, sys, glob, re

## Reference species
species = "human"

## Macaque chrom alias
macaque_chrom_alias = pd.read_csv("/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/genomes/macaque/ncbi/rheMac10.chromAlias.txt", sep="\t", header=None, comment="#")

def extract_genome_and_species_string(path):
    match = re.search(r'To([A-Za-z0-9]+)(?:_([A-Za-z0-9]+))?', path)
    if match:
        genome = match.group(1).lower()
        species = match.group(2).lower() if match.group(2) else ""
        return f"{genome}_{species}" if species else genome
    return ""

species_naming = {
    "hg38": "human",
    "rhemac10": "macaque",
    "hmba_marmoset": "marmoset",
    "mm10": "mouse",
}

## Peak universe for species
peak_universe_species = pd.read_csv(f"/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/ATAC/{species}/Consensus_peaks/merged_peaks_with_names.bed", sep="\t", comment="#", header=None)

## Find liftover results
liftoverDir = f"/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/ATAC/{species}/Consensus_peaks"
liftoverResults = glob.glob(f"{liftoverDir}/merged_peaks_with_names_*.bed")

## Read in liftover results
liftoverDf = pd.DataFrame()
for file in liftoverResults:
    print(f"Reading {file}")
    liftover = pd.read_csv(file, sep="\t", comment="#", header=None)
    liftover.columns = ["Chromosome", "Start", "End", f"{species}_peak_name"]
    ## Identify the genome and species from the filename
    match = extract_genome_and_species_string(file)
    species_name = species_naming.get(match, "unknown")
    ## Add species information
    liftover["speciesFrom"] = species
    liftover["speciesTo"] = species_name
    ## Check if file is liftOver or unmapped
    if "unmapped" in file:
        liftover["liftOver"] = "unmapped"
    else:
        liftover["liftOver"] = "mapped"
    ## Combine into a results dataframe
    liftoverDf = pd.concat([liftoverDf, liftover], ignore_index=True)

## Update for macaque chrom alias
liftoverDf.Chromosome.replace(dict(zip(macaque_chrom_alias[2], macaque_chrom_alias[0])), inplace=True)
## Only keep chromosomes with chr1-chrY, ignore chrM and other non-standard contigs
liftoverDf = liftoverDf[liftoverDf["Chromosome"].str.match(r"^chr[0-9XY]+$")]
liftoverDf.reset_index(drop=True, inplace=True)

## ------------------- Identify orthologous / species specific peaks -------------------
##
human_df = liftoverDf[liftoverDf["speciesFrom"] == species].copy()
pivot = human_df.pivot_table(
    index=f"{species}_peak_name",
    columns="speciesTo",
    values="liftOver",
    aggfunc="first" ## Should be fine since there's only one entry per peak per species
)

## Used for accessologs, peak lookup table into common peak id space based on human peaks
# peak_lookup = human_df.copy()
# peak_lookup['species_peak_name'] = peak_lookup['Chromosome'] + ':' + peak_lookup['Start'].astype(str) + '-' + peak_lookup['End'].astype(str)
# peak_lookup = peak_lookup.pivot_table(index=f'{species}_peak_name', 
#     columns='speciesTo', 
#     values='species_peak_name',
#     aggfunc='first'
# )
# identical_rows = peak_lookup.apply(lambda row: len(set(row.dropna())) == 1, axis=1)
# peak_lookup = peak_lookup[~identical_rows].reset_index()
# peak_lookup.to_csv(liftoverDir + f"/peak_lookup_ref_{species}.csv", index=False)

## All species (excluding human) that were targeted
target_species = [s for s in pivot.columns if s != "human"]

## Boolean mask: All target species must be "unmapped"
human_specific_mask = (pivot[target_species] == "unmapped").all(axis=1)
human_specific_peaks = pivot[human_specific_mask].index.tolist()

## Identify the set of shared peaks mapped in all species
mapped_in_all_mask = (pivot[["macaque", "marmoset"]] == "mapped").all(axis=1)
primate_peaks = pivot[mapped_in_all_mask].index.tolist()

## Identify the set of shared peaks mapped in all species
mapped_in_all_mask = (pivot[target_species] == "mapped").all(axis=1)
mammalian_peaks = pivot[mapped_in_all_mask].index.tolist()

## ------------------- Activity information -------------------
## At this point we need to bring in cell type activity information
peak_table = pd.DataFrame()
for species_name in ["human", "macaque", "mouse"]:
    print(f"Reading peak table for {species_name}")
    ## Read in the peak table for the species
    species_peak_table = pd.read_csv(f"/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/ATAC/{species_name}/Consensus_peaks/Group_by_peaks.csv")
    species_peak_table["species"] = species_name
    ## Split out chrom:start-end into separate columns
    species_peak_table[["Chromosome", "start-end"]] = species_peak_table["Peaks"].str.split(":", expand=True)
    species_peak_table[["Start", "End"]] = species_peak_table["start-end"].str.split("-", expand=True)
    species_peak_table["Start"] = species_peak_table["Start"].astype(int)
    species_peak_table["End"] = species_peak_table["End"].astype(int)
    del species_peak_table["start-end"]
    ## macaque chrom alias
    if species_name == "macaque":
        species_peak_table["Chromosome"].replace(dict(zip(macaque_chrom_alias[2], macaque_chrom_alias[0])), inplace=True)
    ## Only keep chromosomes with chr1-chrY, ignore chrM and other non-standard contigs
    species_peak_table = species_peak_table[species_peak_table["Chromosome"].str.match(r"^chr[0-9XY]+$")]
    species_peak_table.reset_index(drop=True, inplace=True)
    ## No need to liftover reference species
    if species_name == species:
        ## Add human specific peaks
        species_peak_table[f"{species_name}_peak_name"] = species_peak_table["Peaks"]
        species_peak_table = species_peak_table[species_peak_table[f"{species_name}_peak_name"].isin(mammalian_peaks)]
        species_peak_table.reset_index(drop=True, inplace=True)
    else:
        ## Intersect with the liftover results
        gr = pr.PyRanges(species_peak_table)
        gr_liftover =  pr.PyRanges(liftoverDf[liftoverDf["speciesTo"] == species_name])
        ##
        intersected = gr.join(gr_liftover, report_overlap=True).df
        intersected = intersected.loc[intersected["Overlap"] > 50] ## 10% minimum overlap
        ## Only keep the best match for each species peak, multiple species peaks can map to the same species_name peak region
        intersected = (
            intersected
            .sort_values("Overlap", ascending=False)
            .drop_duplicates(subset=[f"{species}_peak_name"])
        )
        ## Filter to mammlian peaks
        intersected = intersected[intersected[f"{species}_peak_name"].isin(mammalian_peaks)]
        species_peak_table = intersected.reset_index(drop=True)
    ## Add to collection
    peak_table = pd.concat([peak_table, species_peak_table], ignore_index=True)
    print(f"Finished reading peak table for {species_name}, shape: {species_peak_table.shape}")

##
peak_lookup = peak_table[["Peaks", "human_peak_name", "species"]]
peak_lookup.to_csv(liftoverDir + f"/peak_lookup_ref_{species}.csv", index=False)

## Identify orthologous regions with peak across all species
peak_set = [group for _, group in peak_table.groupby('species')]
activity_conserved_peaks = set.intersection(*[set(df['human_peak_name']) for df in peak_set])

## Filter to only orthologous peak regions across all species
peak_table_common = peak_table[peak_table['human_peak_name'].isin(activity_conserved_peaks)].copy()
peak_table_common.reset_index(drop=True, inplace=True)

## Identify cell type vs. metadata columns
non_celltype_cols = ["Peaks", "species", "Chromosome", "Start", "End", "human_peak_name", 'Start_b', 'End_b', 'speciesFrom', 'speciesTo', 'liftOver', 'Overlap']
celltype_cols = [col for col in peak_table.columns if col not in non_celltype_cols]

## Step 2: Melt the DataFrame
peak_table_long = peak_table.melt(
    id_vars=[f"{species}_peak_name", "species"],
    value_vars=celltype_cols,
    var_name="celltype",
    value_name="is_present"
)

##
peak_table_long.fillna(False, inplace=True)

##
conserved_any = (
    peak_table_long[peak_table_long["is_present"]]
    .groupby(f"{species}_peak_name")["species"]
    .nunique()
    .reset_index(name="n_species_with_peak")
)
conserved_any["conserved_any_celltype"] = conserved_any["n_species_with_peak"] == 3

## Look for peaks conserved across species and cell types
present_df = peak_table_long[peak_table_long["is_present"]].copy()
celltype_sets = (
    present_df
    .groupby([f"{species}_peak_name", "species"])["celltype"]
    .apply(lambda x: frozenset(x))
    .reset_index(name="active_celltypes")
)

# Step 3: For each peak, find how many species have the *same* set of active cell types
# → Count how many times each unique set appears per peak
conserved_same_set = (
    celltype_sets
    .groupby([f"{species}_peak_name", "active_celltypes"])
    .size()
    .reset_index(name="n_species_with_this_set")
)

# Step 4: Keep only those with ≥2 species sharing the same active celltype set
conserved_peaks = (
    conserved_same_set[conserved_same_set["n_species_with_this_set"] >= 3]
    .human_peak_name
    .unique()
)

## ---------------- Plotting the conservation levels --------
## We will plot the conservation levels as nested circles, with the largest circle representing human peaks,
## and smaller circles representing primate and mammalian levels.
## The "Specific" circle will be placed between the human and primate levels.
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import numpy as np

## Data for plotting
labels = [
    ("Human", peak_universe_species.shape[0]),
    ("Primate Seq. Conserved", len(primate_peaks)),
    ("Mammalian Seq. Conserved", len(mammalian_peaks)),
    ("Mammalian Activity Conserved", len(activity_conserved_peaks)),
]
specific_label = ("Human Specific", len(human_specific_peaks))

## Convert count to radius  
def area_to_radius(area):
    return np.sqrt(area / np.pi)

## Get all radii proportional to the largest count (`species`)
max_area = labels[0][1]
radii = [area_to_radius(n / max_area) for _, n in labels]
specific_radius = area_to_radius(specific_label[1] / max_area)

## Plot setup
fig, ax = plt.subplots(figsize=(6, 6))
ax.set_aspect('equal')
ax.axis('off')

center = (0, 0)
colors = ["#377eb8", "#4bf8f8", "#82f463", "#b10aa6", "#a52222"]

## Draw nested conservation bubbles
for i, ((label, count), r, color) in enumerate(zip(labels, radii, colors)):
    circ = Circle(center, r, facecolor=color, alpha=0.8, edgecolor='black', linewidth=1.5)
    ax.add_patch(circ)

## Ensure "Specific" fits inside Human but outside Mammal 0
human_radius = radii[0]
mammal0_radius = radii[1]
margin = 0.01

min_d = mammal0_radius + specific_radius + margin
max_d = human_radius - specific_radius - margin

if min_d < max_d:
    raise ValueError("Specific circle cannot fit between Mammal 0 and Human with current proportions.")

## Place "Specific" at angle (e.g. 45 degrees) at midpoint distance
specific_distance = (min_d + max_d) / 2
angle_rad = np.pi / 4
x = specific_distance * np.cos(angle_rad)
y = specific_distance * np.sin(angle_rad)
specific_center = (x, y)

## Draw "Specific"
specific_color = "#FE0000"
specific_circle = Circle(specific_center, specific_radius/4, facecolor=specific_color, alpha=0.7, edgecolor='black', linewidth=1.5)
ax.add_patch(specific_circle)

## Finalize and save
plt.title("Human Peak Conservation and Specificity (Proportional, Unscaled)", fontsize=12)
plt.tight_layout()
plt.savefig("/home/nelson.johansen/nested_conservation_levels.pdf", dpi=600)