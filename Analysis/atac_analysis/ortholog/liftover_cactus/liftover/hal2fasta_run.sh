#!/bin/bash
hal2maf \
  --refGenome Homo_sapiens \
  --refSequence chr8 \
  --start 121814459 \
  --length 500 \
  --noAncestors \
  /allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Matthew/genome/hal/zoonomia/447-mammalian-2022v1_replaceAIBSexceptHuman.hal \
  chr8_121814459_peak.maf

