## Scripts to perform conserved and species-specific marker gene analysis

### Step 1 compute marker metrics (Tau / Auroc)
* 0_calc_auroc.py
* 0_calc_tau_beta.py

### Step 2 curate marker sets and marker .h5ad files for plotting
* 1_create_gene_sets.py

### Step 3: Plot using ComplexHeatmap and taxonomy info.
* 2_expression_plots.R

### Step 4-5: Species specific marker gene analysis (IN-PROGRESS, use with caution).
* 3_species_specific_analysis.py
* 4_expression_plots_per_species.R
