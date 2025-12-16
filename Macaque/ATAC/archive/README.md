## Snapatac2 processsing of Macaque ATAC (10X Multiome)

The pipeline for processing the Macaque ATAC data using snapatac2 follows this process:

----

`build_ATAC_study.py`: Identify alignment directories from cellranger containing fragment files per library. NOTE: This will change depending on how/who aligned the data.

* Input: Taxonomy .h5ad file from the Allen Institute
* Output: .csv table containing alignment id and directory information per library.

----

`snapatac2_preprocess.py/.sh`: Reading the .csv from the previous step to process the fragment file into snapatac2 .h5ad as well as computing basic QC metrics. The `.sh` script runs a SLURM JOB ARRAY to process the fragment files in parallel

* Input: `build_ATAC_study.py` .csv file.
* Output: snapatac2 `.h5ad` per library.

----

`snapatac2_anndataset_build.ipynb`: Concatenates the snapatac2 .h5ad file into an AnnDataSet object (.h5ads)

* Input: Location of snapatac2 `.h5ad` files.
* Output: Concatenated .h5ads file.

----

`create_bigwigs.py`: Combines the RNA and ATAC data files to derive per-annotation level bigwig files.

* Input: RNA .h5ad, ATAC .h5ads files.
* Output: .bw file per annotation in the level specified.

----

`call_peaks.py`: Calls consensus peaks across annotation level.

* Input: Subset .h5ads file with cell type labels.
* Output: Consensus peak .bed file.
  
 
