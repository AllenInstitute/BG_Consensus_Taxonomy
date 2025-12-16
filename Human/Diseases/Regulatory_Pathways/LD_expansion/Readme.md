# LD Expansion Pipeline (PLINK)

This pipeline performs LD (Linkage Disequilibrium) expansion of significant SNPs using PLINK, with parallelization over multiple SLURM jobs.

---

## Overview

- Input: list of genome-wide significant SNPs (one per line)
- Output: single LD-expanded `.ld` file
- Parallelized across chunks using SLURM job arrays

---

## Step 1 — Split SNP list for parallelization

Run the following shell commands to split your list of genome-wide significant SNPs into `N` chunks (for parallel LD expansion):

```bash
# Number of jobs (chunks)
N=100
SNPLIST=/path/to/significant_snplist.txt
CHUNKDIR=/path/to/LD_chunks

mkdir -p $CHUNKDIR
split -d -n l/$N $SNPLIST $CHUNKDIR/snps_chunk_
```

---

## Step 2 — Run LD Expansion via SLURM Job Script

The expansion is parallelized using a SLURM job array. Each job processes a chunk of SNPs.

Run the following script to submit the job array:

```bash
sbatch ld_expansion.sh 
```

---

## Step 3 - Aggregate all LD results


```bash
OUTDIR=/path/to/LD_temp
FINAL=/path/to/LD_expanded_all.ld

# Use header from first file
head -n 1 $OUTDIR/LD_chunk_00.ld > $FINAL

# Concatenate the rest
cat $OUTDIR/LD_chunk_*.ld >> $FINAL
```



