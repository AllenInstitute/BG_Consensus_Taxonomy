#!/bin/bash

BASE_DIR=/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/genomes/

species=(hu mq ms mm) ## Human, Macaque, Marmoset, Mouse
files=(${BASE_DIR}/human/samap/species.cds.cleaned.fa \
       ${BASE_DIR}/macaque/ncbi/samap/species.cds.cleaned.fa \
       ${BASE_DIR}/marmoset/samap/species.cds.cleaned.fa \
       ${BASE_DIR}/mouse/samap/species.cds.cleaned.fa)

types=(nucl nucl nucl nucl) ## change to 'prot' if you have protein FASTAs

## loop over all unique pairs
for i in ${!species[@]}; do
    for j in $(seq $((i+1)) $((${#species[@]}-1))); do
        echo "Running SAMap: ${species[i]} vs ${species[j]}"
        bash /home/nelson.johansen/scripts/python/github/samap_directory/map_genes.sh \
            --tr1 ${files[i]} --t1 ${types[i]} --n1 ${species[i]} \
            --tr2 ${files[j]} --t2 ${types[j]} --n2 ${species[j]}
    done
done

## Some libraries aren't installed on worker nodes....
# ## loop over all unique pairs
# for i in ${!species[@]}; do
#     for j in $(seq $((i+1)) $((${#species[@]}-1))); do
#         job_name="SAMap_${species[i]}_${species[j]}"
#         echo "Submitting $job_name"

#         sbatch <<EOT
# #!/bin/bash
# #SBATCH --job-name=$job_name
# #SBATCH --output=logs/$job_name.%j.out
# #SBATCH --error=logs/$job_name.%j.err
# #SBATCH --time=12:00:00
# #SBATCH --partition=celltypes
# #SBATCH --mem=64G
# #SBATCH --cpus-per-task=8

# conda activate SAMap

# bash /home/nelson.johansen/scripts/python/github/samap_directory/map_genes.sh \
#     --tr1 ${files[i]} --t1 ${types[i]} --n1 ${species[i]} \
#     --tr2 ${files[j]} --t2 ${types[j]} --n2 ${species[j]}
# EOT

#     done
# done