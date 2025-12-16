#!/usr/bin/bash
#SBATCH --job-name=Full_MAGMA_pipeline
#SBATCH --partition celltypes 
#SBATCH --time=5:00:00
#SBATCH --mem=25G

# The 5 following scripts are adapted from Duncan et al. (2024)

# Step 1
echo "[1/5] Execution of step1.sh..."
bash Scripts/Gene_annotation_step1.sh gwas_config.env

# Step 2
echo "[2/5] Execution of step2.sh..."
bash Scripts/Gene_property_analysis_step2.sh gwas_config.env

# Step 3
echo "[3/5] Execution of analysis. 3..."
bash Scripts/run_script3.sh gwas_config.env

# Step 4
echo "[4/5] Execution of step4.sh..."
bash Scripts/conditional_analysis_step4.sh gwas_config.env

# Step 5
echo "[5/5] Execution of step5.R..."
bash Scripts/run_script5.sh gwas_config.env
