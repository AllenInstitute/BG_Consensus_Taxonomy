import pandas as pd

## Load in liftover results
halper_liftover = pd.read_csv("/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/ATAC/annotations/human_ref_liftover_HALPER_minMatch_0-5.tsv", sep="\t")
custom_liftover = pd.read_csv("/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/ATAC/annotations/human_ref_liftover_minMatch_0-5.tsv", sep="\t")
custom_liftover.rename(columns={"Macaca_mulatta_ID": "macaque_ID", 
                                "Macaca_mulatta_aligned_ID": "macaque_aligned_ID",
                                "Mus_musculus_ID": "mouse_ID",
                                "Mus_musculus_aligned_ID": "mouse_aligned_ID"}, inplace=True)

## Load in macaque chrom alias
macaque_chrom_alias = pd.read_csv("/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/genomes/macaque/ncbi/rheMac10.chromAlias.txt", sep="\t", header=None, comment="#", names=["alias", "type1", "type2"])
macaque_alias_dict = dict(zip(macaque_chrom_alias["type2"], macaque_chrom_alias["alias"]))

## For each liftover file update macaque to alias
halper_liftover.macaque_chrom = halper_liftover["macaque_ID"].str.split(":").str[0].map(macaque_alias_dict)
custom_liftover.macaque_chrom = custom_liftover["macaque_ID"].str.split(":").str[0].map(macaque_alias_dict)

## Now update the IDs with new chrom
halper_liftover.macaque_ID = halper_liftover.macaque_chrom + ":" + halper_liftover["macaque_ID"].str.split(":").str[1]
custom_liftover.macaque_ID = custom_liftover.macaque_chrom + ":" + custom_liftover["macaque_ID"].str.split(":").str[1]

##
custom_liftover.rename(columns={"human_ID": "human_ID_custom", 
                                "macaque_ID": "macaque_ID_custom", 
                                "macaque_aligned_ID": "macaque_aligned_ID_custom",
                                "mouse_ID": "mouse_ID_custom", 
                                "mouse_aligned_ID": "mouse_aligned_ID_custom",
                                "ortholog": "ortholog_custom"}, inplace=True)

## Common orthologs merge in just human, macaque, mouse IDs from custom
common_orthologs = pd.concat([
    halper_liftover,
    custom_liftover[["human_ID_custom", "macaque_ID_custom", "macaque_aligned_ID_custom", "mouse_ID_custom", "mouse_aligned_ID_custom", "ortholog_custom"]]
], axis=1, join="outer")

## Keep only orthologs with all three species
common_orthologs = common_orthologs[~common_orthologs["macaque_ID"].isna() & ~common_orthologs["macaque_ID_custom"].isna()].copy()
common_orthologs = common_orthologs[common_orthologs["ortholog"] & common_orthologs["ortholog_custom"]].copy()
common_orthologs.reset_index(drop=True, inplace=True)

## Create custom and halper pyranges then intersect the macaque_ID columns
halper_pr = pr.PyRanges(pd.DataFrame({"Chromosome":common_orthologs["macaque_aligned_ID"].str.split(":").str[0], 
                        "Start":common_orthologs["macaque_aligned_ID"].str.split(":").str[1].str.split("-").str[0].astype(int),
                        "End":common_orthologs["macaque_aligned_ID"].str.split(":").str[1].str.split("-").str[1].astype(int),
                        "ID":common_orthologs["human_ID"]})
)

custom_pr = pr.PyRanges(pd.DataFrame({"Chromosome":common_orthologs["macaque_aligned_ID_custom"].str.split(":").str[0],
                        "Start":common_orthologs["macaque_aligned_ID_custom"].str.split(":").str[1].str.split("-").str[0].astype(int),
                        "End":common_orthologs["macaque_aligned_ID_custom"].str.split(":").str[1].str.split("-").str[1].astype(int),
                        "ID":common_orthologs["human_ID"]})
)
intersect = halper_pr.join(custom_pr, suffix="_custom", report_overlap=True)
len(intersect.df.ID.unique())

## Write out some examples
common_orthologs.iloc[sample(list(common_orthologs.index), 10)].to_csv("/home/nelson.johansen/common_ortholog_examples.csv", index=False, header=True)

## Write out halper biased ortholog set
halper_liftover_orthologs = halper_liftover[halper_liftover["ortholog"]].copy()
halper_liftover_orthologs = halper_liftover_orthologs[~halper_liftover_orthologs["human_ID"].isin(common_orthologs["human_ID"])].reset_index(drop=True)
halper_liftover_orthologs.iloc[sample(list(halper_liftover_orthologs.index), 10)].to_csv("/home/nelson.johansen/halper_liftover_ortholog_examples.csv", index=False, header=True)

## Now custom biased ortholog set
custom_liftover_orthologs = custom_liftover[custom_liftover["ortholog_custom"]].copy()
custom_liftover_orthologs = custom_liftover_orthologs[~custom_liftover_orthologs["human_ID_custom"].isin(common_orthologs["human_ID"])].reset_index(drop=True)
custom_liftover_orthologs.iloc[sample(list(custom_liftover_orthologs.index), 10)].to_csv("/home/nelson.johansen/custom_liftover_ortholog_examples.csv", index=False, header=True)