import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import mannwhitneyu
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42


# ------------------------------
# Load data
# ------------------------------
df = pd.read_csv("./all_peaks/df_for_features_updated.tsv", sep="\t")
df['target_group'] = df['status'].apply(lambda x: 'on-target' if x in ['high', 'medium', 'low'] else 'off-target')
df = df.rename(columns={'magnitude_specificity': 'ATAC_specificity'})

# ------------------------------
# Features to plot
# ------------------------------
features = ['ATAC_specificity', 'CREsted_specificity']
palette_colors = {'on-target': '#89CFF0', 'off-target': '#FDBB84'}  # light blue and orange

# Function to convert p-value to stars
def pval_to_star(p):
    if p < 0.001:
        return '***'
    elif p < 0.01:
        return '**'
    elif p < 0.05:
        return '*'
    else:
        return 'ns'  # not significant

plt.figure(figsize=(8,6))

for i, feature in enumerate(features):
    df_feat = df[['group','target_group', feature]].copy()
    df_feat = df_feat.rename(columns={feature:'Value'})
    
    positions = [i*2, i*2 + 0.8]  # on/off box positions
    
    # Draw boxes
    for j, target in enumerate(['on-target','off-target']):
        subset = df_feat[df_feat['target_group']==target]['Value']
        plt.boxplot(subset, positions=[positions[j]], widths=0.6, patch_artist=True,
                    boxprops=dict(facecolor=palette_colors[target], color='k'),
                    medianprops=dict(color='k', linewidth=2),
                    whiskerprops=dict(color='k'),
                    capprops=dict(color='k'),
                    flierprops=dict(marker='o', markersize=0, alpha=0))  # hide default fliers
    
    # Overlay points in front
    for j, target in enumerate(['on-target','off-target']):
        subset = df_feat[df_feat['target_group']==target]['Value']
        x = np.random.normal(positions[j], 0.08, size=len(subset))
        plt.scatter(x, subset.values, color='k', alpha=0.6, s=15, zorder=5)
    
    # Compute p-value and add star
    on_values = df_feat[df_feat['target_group']=='on-target']['Value']
    off_values = df_feat[df_feat['target_group']=='off-target']['Value']
    stat, pval = mannwhitneyu(on_values, off_values, alternative='two-sided')
    y_max = max(df_feat['Value']) * 1.05
    plt.text(np.mean(positions), y_max, f"p = {pval:.3e}", ha='center', fontsize=12, fontweight='bold')

    #star = pval_to_star(pval)
    #plt.text(np.mean(positions), y_max, star, ha='center', fontsize=12, fontweight='bold')

# X-axis labels
plt.xticks([0.4,2.4], features, fontsize=12)
plt.ylabel('Specificity', fontsize=12)
plt.xlim(-0.5,3.5)
plt.tight_layout()
plt.savefig("ATAC_CREsted_specificity_boxplot_on_off_points_front.pdf", dpi=300)
plt.show()

