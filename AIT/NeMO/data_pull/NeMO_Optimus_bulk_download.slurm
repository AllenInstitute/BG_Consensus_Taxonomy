#!/bin/bash
#SBATCH --job-name=NeMO_download    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=nelson.johansen@alleninstitute.org     # Where to send mail  
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=16gb                     # Job memory request (per node)
#SBATCH --time=04:00:00               # Time limit hrs:min:sec
#SBATCH --output=NeMO_download_%j.log   # Standard output and error log
#SBATCH --partition celltypes         # Partition used for processing
#SBATCH --tmp=16G                     # Request the amount of space your jobs needs on /scratch/fast
 
## Handle arguments
modality=$1 
species=$2 
bucket_name=$3 
alignment_job_id=$4 
library_id=$5
output_dir=$6

## Construct local path and create directory if it doesn't exist
local_path="${output_dir}/${species}/NeMO/${modality}/${library_id}"
mkdir -p "$local_path"

## Convert library_id to only contain lowercase alphanumeric characters but keep any special characters lilke "-" or "_"
# library_id=$(echo "$library_id" | tr '[:upper:]' '[:lower:]' | tr -cd '[:alnum:]-_')

## Construct the GCS file path for listing
gcs_list_path="gs://${bucket_name}/analysis/*${library_id}*.h5ad"
echo "Listing files: $gcs_list_path"

## List all matching files, if none found then skip.
matching_files=$(gsutil ls "$gcs_list_path")
if [ -z "$matching_files" ]; then
    echo "No files found for $gcs_list_path"
    exit 1
fi

## Filter files containing provided alignment_job_id
selected_file=$(echo "$matching_files" | grep "$alignment_job_id" | head -n 1)
if [ -z "$selected_file" ]; then
    # If no "JobId_" file exists, download first file
    selected_file=$(echo "$matching_files" | head -n 1)
else
    echo "'JobId_' file found: $selected_file"
fi

## Download the selected file
if gsutil cp "$selected_file" "$local_path/" 2>&1; then
    echo "Successfully downloaded $selected_file"
else
    echo "Failed to download $selected_file"
fi
