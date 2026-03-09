#!/bin/bash

#SBATCH --job-name=halliftover
#SBATCH --time 24:00:00
#SBATCH --ntasks=4
#SBATCH --partition=celltypes
#SBATCH --error=logs/halliftover_%A_%a.err.txt
#SBATCH --output=logs/halliftover_%A_%a.out.txt

## Loop through 1:447 and call python script to make ortholog matrix for each species
# for i in {0..1}
# do
#         echo "Processing index $i"
#         ## Add SBATCH parameters as needed
#         sbatch --job-name=halliftover_$i --time=12:00:00 --ntasks=4 --mem 64GB --partition=celltypes --error=/home/nelson.johansen/logs/halliftover_${i}_%A.err.txt --output=/home/nelson.johansen/logs/halliftover_${i}_%A.out.txt \
#                 --wrap="python /home/nelson.johansen/repos/halLiftover-postprocessing/makePeakOrthologMatrix_sbatch.py \
#                 --bedFileName /allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/ATAC/human/Consensus_peaks/merged_peaks_with_names.bed \
#                 --orthologsFileNameListFileName /allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/ATAC/hal/human/archive/halper_files.txt \
#                 --orthologsFileNameListFileName_idx $i \
#                 --outputFileName /allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/ATAC/hal/human/ortholog_peak_matrix_"
# done

# find . -type f -name "*with_names*HALPER.narrowPeak.gz" ! -name "*Anc*" | sed 's|^\./||' > halper_files.txt

python /home/nelson.johansen/repos/halLiftover-postprocessing/makePeakOrthologMatrix.py \
        --bedFileName /allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/ATAC/human/Consensus_peaks/merged_peaks_with_names.bed \
        --orthologsFileNameListFileName /allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/ATAC/hal/human/archive/halper_files.txt \
        --outputFileName /allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/ATAC/hal/human/ortholog_peak_matrix.txt

## Count number of 1s in each row of the resulting ortholog matrix and write to a .txt file
while IFS= read -r line; do
  count=$(echo "$line" | grep -o "1" | wc -l)
  ## Write to a file
  ## First write the 1st column (peak name) and then the count
  species_name=$(echo "$line" | cut -f1)
  echo "$species_name,$count" >> ortholog_counts_per_species.txt
  echo "Number of 1s in this row: $count"
done < ortholog_peak_matrix.txt

## Remove the first line of the file (header)
sed -i '1d' ortholog_counts_per_species.txt

##
/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/ATAC/human/Consensus_peaks/merged_peaks_with_names.bed
zcat halLiftover.sFile.bed.gz | awk '{len += $3-$2} END{print len}'
zcat halLiftover.tFile.bed.gz | awk '{len += $3-$2} END{print len}'
fraction = tFile_length / sFile_length