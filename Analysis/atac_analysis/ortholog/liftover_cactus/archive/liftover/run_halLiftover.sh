## awk 'OFS="\t"{print $0, $1":"$2"-"$3}' merged_peaks.bed > merged_peaks_with_names.bed
## Mouse -> Human
halLiftover \
	/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Matthew/genome/hal/zoonomia/447-mammalian-2022v1_replaceAIBSexceptHuman.hal \
	Mus_musculus \
	/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/ATAC/mouse/Consensus_peaks/merged_peaks_with_names.bed \
	Homo_sapiens \
	/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/ATAC/hal/human/Mus_musculusToHomo_sapiens.bed

## Macaque -> Human
halLiftover \
	/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Matthew/genome/hal/zoonomia/447-mammalian-2022v1_replaceAIBSexceptHuman.hal \
	Macaca_mulatta \
	/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/ATAC/macaque/Consensus_peaks/merged_peaks_with_names.bed \
	Homo_sapiens \
	/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/ATAC/hal/human/Macaca_mulattaToHomo_sapiens.bed