## Loop through species human, macaque and mouse
# species=("human" "macaque" "marmoset")
# species_sci=("Homo_sapiens" "Macaca_mulatta" "Callithrix_jacchus")

species=("marmoset")
species_sci=("Callithrix_jacchus")

for i in "${!species[@]}"; do
    sp=${species[$i]}
    sp_sci=${species_sci[$i]}
    sbatch \
        -p celltypes \
        --array 1-892 \
        /home/nelson.johansen/Analysis/HMBA_Genomics/BasalGanglia/Analysis/atac_analysis/ortholog/liftover_cactus/ortholog_CACTUS/halper_map_peak_orthologs.sh \
	-b /allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/basal-ganglia/data/$sp/peaks/merged_peaks_updated_chrom.bed \
	-o /allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis/cactus/hal/$sp \
	-s $sp_sci \
	-t /allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis/cactus/hal/species.txt \
	-c /allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Matthew/genome/hal/zoonomia/447-mammalian-2022v1_replaceAIBSexceptHuman.hal \
	-f 'FALSE'
done


# sacct -u nelson.johansen --starttime=2025-09-01 \
#        --format=Elapsed -P | \
# awk -F'|' 'NR>1 {
#     split($1, t, ":");
#     hours = t[1]; minutes = t[2]; seconds = t[3];
#     total_hours += (hours + minutes/60 + seconds/3600)
# } END {print "Total wall-clock hours:", total_hours}'