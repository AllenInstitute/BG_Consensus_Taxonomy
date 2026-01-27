import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import TwoSlopeNorm
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42

# ------------------------------
# Load data
# ------------------------------
df = pd.read_csv("df_for_features.tsv", sep="\t")
df = df.rename(columns={"magnitude_specificity": "ATAC_specificity"})

# Map status to on/off
df['status'] = df['status'].str.lower()
df['on_off'] = df['status'].map({
    'high': 'on',
    'medium': 'on',
    'low': 'on',
    'none': 'off'
})

df['group'] = df['group'].str.replace('.', '-', regex=False)

# ------------------------------
# Group mapping
# ------------------------------
group_map = {
    "SN-VTR_CALB1_Dopa,SN-VTR_GAD2_Dopa,SN_SOX6_Dopa": "SN_TH_Dopamine",
    "STRd_D1_Matrix_MSN,STRd_D1_Striosome_MSN,STRv_D1_MSN": "STR_D1_MSN",
    "STRd_D2_Matrix_MSN,STRd_D2_Striosome_MSN,STRv_D2_MSN,STRd_D2_StrioMat_Hybrid_MSN": "STR_D2_MSN"
}
df['broad_group'] = df['group'].map(group_map).fillna(df['group'])

# ------------------------------
# Compute mean on/off ratios per broad group
# ------------------------------
groups = df['broad_group'].unique()
ratios = []

for g in groups:
    gdf = df[df['broad_group'] == g]

    on_crested = gdf.loc[gdf['on_off'] == 'on', 'CREsted_specificity'].mean()
    off_crested = gdf.loc[gdf['on_off'] == 'off', 'CREsted_specificity'].mean()
    crested_ratio = on_crested / off_crested if off_crested != 0 else np.nan

    on_atac = gdf.loc[gdf['on_off'] == 'on', 'ATAC_specificity'].mean()
    off_atac = gdf.loc[gdf['on_off'] == 'off', 'ATAC_specificity'].mean()
    atac_ratio = on_atac / off_atac if off_atac != 0 else np.nan

    ratios.append({
        'group': g,
        'CREsted': crested_ratio,
        'ATAC': atac_ratio,
        'diff': atac_ratio - crested_ratio
    })

ratios_df = pd.DataFrame(ratios)
ratios_df = ratios_df.sort_values('diff', ascending=False).reset_index(drop=True)

# ------------------------------
# Plot publication-ready
# ------------------------------
fig, ax = plt.subplots(figsize=(8, 4))

# Colormap for connecting lines and symbols
cmap = sns.color_palette("coolwarm", as_cmap=True)
norm = TwoSlopeNorm(vmin=ratios_df['diff'].min(), vcenter=0, vmax=ratios_df['diff'].max())

for _, row in ratios_df.iterrows():
    y_val = row['group']
    # Line colored by difference
    line_color = cmap(norm(row['diff']))
    ax.plot([row['CREsted'], row['ATAC']], [y_val, y_val], color=line_color, linewidth=2, zorder=0)
    # Symbols in same color as line
    ax.scatter(row['CREsted'], y_val, facecolors='none', color=line_color, marker='o', s=100, zorder=5)
    ax.scatter(row['ATAC'], y_val, color=line_color, marker='o', s=80, zorder=5)

# Vertical line at 1
ax.axvline(1, color='k', linestyle='--', alpha=0.5)

# X-axis, Y-axis
ax.set_xlabel("Mean On/Off Specificity Ratio", fontsize=12)
ax.set_ylabel("Group", fontsize=12)
ax.grid(axis='x', linestyle='--', alpha=0.4)
plt.tight_layout()

# Create neutral-colored legend manually
ax.legend([plt.Line2D([0], [0], marker='o', markerfacecolor='none', color='k', linestyle='None', markersize=8),
           plt.Line2D([0], [0], marker='o', color='k', linestyle='None', markersize=8)],
          ['CREsted', 'ATAC'], title='Feature', fontsize=10)

# Colorbar for line/symbol colors
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar = fig.colorbar(sm, ax=ax)
cbar.set_label("ATAC higher → warm|CREsted higher → cool", rotation=270, labelpad=15)

# ------------------------------
# Save publication-ready PDF
# ------------------------------
plt.savefig("CREsted_ATAC_ratio_per_broad_group_pub.pdf", dpi=300, bbox_inches='tight')
plt.show()

