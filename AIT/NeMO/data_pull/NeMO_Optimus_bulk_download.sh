#!/bin/bash

## Input CSV file
manifest_file="/home/nelson.johansen/Analysis/BG/NeMO/nemo_bican_manifest.csv"

## Output directory for downloaded files
output_dir="/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases"

## Read the CSV file line by line, assumption is no header
while IFS=',' read -r modality species bucket_name alignment_job_id library_id; do
    local_path="${output_dir}/${species}/NeMO/${modality}/${library_id}"
    sbatch /home/nelson.johansen/Analysis/BG/NeMO/NeMO_Optimus_bulk_download.slurm $modality $species $bucket_name $alignment_job_id $library_id $output_dir
    # if [ "$species" == "Human" ]; then
    #     ## Missed some files due to permisison issues, this will only pull data for libraries with no files in their directory.
    #     if [ "$(ls -A "$local_path" | wc -l)" -eq 0 ]; then
    #         echo "Downloading" $library_id $modality $bucket_name $alignment_job_id
    #        ## Download the selected file
    #        sbatch /home/nelson.johansen/Analysis/BG/NeMO/NeMO_Optimus_bulk_download.slurm $modality $species $bucket_name $alignment_job_id $library_id $output_dir
    #    fi
    #    #else
    #        #echo "Skipping" $library_id
    #    #fi
    #fi
done < "$manifest_file"
