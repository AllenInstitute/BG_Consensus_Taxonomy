#!/bin/bash
species="$1"
genome="$2"

bash /home/nelson.johansen/Analysis/CounterScreen/liftover.sh \
  "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/ATAC/${species}/Consensus_peaks/merged_peaks_with_names.bed" \
  "$genome" \
  "Hg38" \
  "0.5"

bash /home/nelson.johansen/Analysis/CounterScreen/liftover.sh \
  "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/ATAC/${species}/Consensus_peaks/merged_peaks_with_names.bed" \
  "$genome" \
  "RheMac10" \
  "0.5"

bash /home/nelson.johansen/Analysis/CounterScreen/liftover.sh \
  "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/ATAC/${species}/Consensus_peaks/merged_peaks_with_names.bed" \
  "$genome" \
  "HMBA_marmoset" \
  "0.5"

bash /home/nelson.johansen/Analysis/CounterScreen/liftover.sh \
  "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/ATAC/${species}/Consensus_peaks/merged_peaks_with_names.bed" \
  "$genome" \
  "Mm10" \
  "0.5"