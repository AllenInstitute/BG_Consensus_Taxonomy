import os
import snapatac2 as snap

# Define the root directory to start walking
root_directory = '/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/data/human/ATAC/h5ad'

# Walk through the directory tree
for root, dirs, files in os.walk(root_directory):
    h5ad_files = {}
    for f in files:
        ## Check if the file is an h5ad file
        if f.endswith('.h5ad') & (not '.h5ads' in f) & (not '.h5ad.h5ad' in f) & (not 'concatenated' in f):
            print(root,f)
            file_path = os.path.join(root, f)
            print(f'Reading file: {file_path}')
            
            ## Read the h5ad file
            adata = snap.read(file_path,backed='r')
            h5ad_files[file_path] = adata
    ## Check if we have h5ad files to concatenate
    if h5ad_files:
        ## Concatenate all the h5ad files
        concatenated_adata = snap.AnnDataSet(
            adatas=[(k, h5ad_files[k]) for k in h5ad_files],
            filename=os.path.join(root, 'concatenated_102225.h5ads')
        )
        concatenated_adata.close()

for k in h5ad_files.keys():
    h5ad_files[k].close()
