import pandas as pd
import numpy as np
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import gridspec
import os
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42

# -------------------------------
# 1. Load and preprocess table
# -------------------------------
# df = pd.read_csv("crested_predictions_with_peakid.tsv", sep="\t")
# for ATAC reading this file
df = pd.read_csv("peaks_summary_all_bw.tsv", sep="\t") 
df['status'] = df['status'].fillna('none').str.lower()


# -------------------------------
# 5. Count status per group (High/Medium/Low/None)
# -------------------------------
statuses_order = ["high", "medium", "low", "none"]

status_table = (
    df.groupby(["group", "status"])
      .size()
      .unstack(fill_value=0)
)

# Ensure all 4 columns exist (even if some statuses are missing)
status_table = status_table.reindex(columns=statuses_order, fill_value=0)

# Add totals (optional but useful)
status_table["total"] = status_table.sum(axis=1)

# Sort groups by total peaks (optional)
status_table = status_table.sort_values("total", ascending=False)

print("\nCounts of peaks by Group x Status:")
print(status_table.to_string())

# Save to file (optional)
status_table.to_csv("ATAC_counts_by_group_status.tsv", sep="\t")
print("\nSaved table: ATAC_counts_by_group_status.tsv")

# Merge group names with friendly labels
group_mapping = {
    "STRd_D1_Matrix_MSN,STRd_D1_Striosome_MSN,STRv_D1_MSN": "STR_D1_MSN",
    "STRd_D2_Matrix_MSN,STRd_D2_Striosome_MSN,STRv_D2_MSN,STRd_D2_StrioMat_Hybrid_MSN": "STR_D2_MSN",
    "SN-VTR_CALB1_Dopa,SN-VTR_GAD2_Dopa,SN_SOX6_Dopa": "SN_TH_Dopamine"
}

for k, v in group_mapping.items():
    df.loc[df['group'] == k, 'group'] = v

# -------------------------------
# 2. Load color map from consensus annotation
# -------------------------------
anno_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/sarojaS/250506_bg_figure/merged_Fig_7and8/PyPeakRankr/"
cluster_meta = pd.read_excel(
    os.path.join(anno_dir, "HMBA_BG_consensus_annotation.xlsx"),
    sheet_name="consensus_anno_pre-print"
).drop_duplicates(subset=["Group"]).reset_index(drop=True)
#cluster_meta["Group"].replace(" ", "_", regex=True, inplace=True)

# Build dictionary: group name -> hex color
#color_map = dict(zip(cluster_meta["Group"], cluster_meta["color_hex_group"]))

# Some names in df['group'] might not match 'Group' exactly — normalize a bit
#df["group_norm"] = df["group"].str.replace(" ", "_")

# Assign colors using map; fallback to grey for missing groups
#df["color"] = df["group_norm"].map(color_map).fillna("lightgrey")


# 2. Assign auto-generated unique colors per group (no atlas colors)
#from matplotlib.cm import get_cmap
#import numpy as np

#groups = sorted(df['group'].unique())
#n_groups = len(groups)

## Use tab20 / tab20b / tab20c automatically depending on size
#if n_groups <= 20:
#    cmap = get_cmap("tab20")
#elif n_groups <= 40:
#    cmap = get_cmap("tab20b")
#else:
#    cmap = get_cmap("gist_ncar")  # fallback for very large palettes

#colors = [cmap(i / n_groups) for i in range(n_groups)]
#color_map = dict(zip(groups, colors))
#df["color"] = df["group"].map(color_map)

# -------------------------------
# Use representative group for COLOR ONLY
# -------------------------------
representative_map = {
    #"STRd_D1_Matrix_MSN,STRd_D1_Striosome_MSN,STRv_D1_MSN": "STRd_D1_Matrix_MSN",
    #"STRd_D2_Matrix_MSN,STRd_D2_Striosome_MSN,STRv_D2_MSN,STRd_D2_StrioMat_Hybrid_MSN": "STRd_D2_Matrix_MSN",
    #"SN-VTR_CALB1_Dopa,SN-VTR_GAD2_Dopa,SN_SOX6_Dopa": "SN_SOX6_Dopa"

    "STR_D1_MSN": "STRd_D1_Matrix_MSN",
    "STR_D2_MSN": "STRd_D2_Matrix_MSN",
    "SN_TH_Dopamine": "SN_SOX6_Dopa"


}

# Normalize atlas keys a bit
cluster_meta['Group'] = cluster_meta['Group'].str.replace(" ", "_")

# Build color map from atlas (Group → hex)
color_map = dict(zip(cluster_meta["Group"], cluster_meta["color_hex_group"]))

# Assign color via representative group
df["rep_group"] = df["group"]  # default: itself
for big, rep in representative_map.items():
    df.loc[df["group"] == big, "rep_group"] = rep

df["rep_group_norm"] = df["rep_group"].str.replace(" ", "_")

df["color"] = df["rep_group_norm"].map(color_map).fillna("lightgrey")

# -------------------------------
# Manual color override for specific group
# -------------------------------
df.loc[
    df["group"] == "STH_PVALB-PITX2_Glut",
    "color"
] = "#39e1e2"

df.loc[
    df["group"] == "SN_GATA3-PAX8_GABA",
    "color"
] = "#311a96"

# -------------------------------
# 3. t-SNE on numeric features
# -------------------------------
#meta_cols = ['region_id', 'chrom', 'start', 'end', 'group', 'status', 'peak_id']
#feature_cols = [c for c in df.columns if c not in meta_cols and c not in ['group_norm', 'color']]
#X = df[feature_cols].values

#tsne = TSNE(n_components=2, perplexity=30, random_state=42)
#X_emb = tsne.fit_transform(X)


# -------------------------------
# 3. t-SNE on numeric features
# -------------------------------
non_feature_cols = [
    #for ATAC file adding this line
     'chr', 'start', 'end', 'group', 'status', 'peak_id',
    #'region_id', 'chrom', 'start', 'end', 'group', 'status', 'peak_id',
    'group_norm', 'color', 'rep_group', 'rep_group_norm'
]

feature_cols = [c for c in df.columns if c not in non_feature_cols]
df[feature_cols] = df[feature_cols].apply(pd.to_numeric, errors='coerce')


X = df[feature_cols].values

tsne = TSNE(n_components=2, perplexity=30, random_state=42)
X_emb = tsne.fit_transform(X)



# -------------------------------
# 4. Combined figure: t-SNE + legend
# -------------------------------
fig = plt.figure(figsize=(14, 12))
gs = gridspec.GridSpec(1, 2, width_ratios=[4, 1], wspace=0.05)

# Main t-SNE plot
ax_tsne = fig.add_subplot(gs[0])
groups = df['group'].unique()
statuses = ['high', 'medium', 'low', 'none']

markers = {'high': 'o', 'medium': 'o', 'low': 'o', 'none': 'x'}
sizes = {'high': 120, 'medium': 90, 'low': 60, 'none': 80}

for g in groups:
    group_color = df.loc[df['group'] == g, 'color'].iloc[0]
    for s in statuses:
        idx = (df['group'] == g) & (df['status'] == s)
        if idx.sum() == 0:
            continue
        alpha_val = 0.8 if s != 'none' else 0.7
        ax_tsne.scatter(
            X_emb[idx, 0], X_emb[idx, 1],
            c=[group_color],
            marker=markers[s],
            s=sizes[s],
            alpha=alpha_val,
            edgecolors=group_color if markers[s] != 'x' else None,
            linewidths=0.8
        )


# -------------------------------
# 4a. Highlight OT_D1_ICj 'none' points with peak_id
# -------------------------------
highlight_idx = (df['group'] == 'OT_D1_ICj') & (df['status'] == 'none')
for i in np.where(highlight_idx)[0]:
    ax_tsne.text(
        X_emb[i, 0], X_emb[i, 1],
        df.loc[i, 'peak_id'],
        fontsize=7,
        color='red',
        alpha=0.8
    )


ax_tsne.set_xlabel("t-SNE 1", fontsize=12)
ax_tsne.set_ylabel("t-SNE 2", fontsize=12)
ax_tsne.set_title("t-SNE of Peaks by Group and Status", fontsize=14)
ax_tsne.set_aspect('equal', 'box')

# -------------------------------
# Legend panel
# -------------------------------
ax_legend = fig.add_subplot(gs[1])
ax_legend.axis('off')

legend_elements = []

# Groups first
for g in groups:
    c = df.loc[df['group'] == g, 'color'].iloc[0]
    legend_elements.append(Line2D([0], [0], marker='o', color='w',
                                  label=g, markerfacecolor=c,
                                  markersize=10, markeredgecolor=c))

# Status markers
status_labels = {"high": "High", "medium": "Medium", "low": "Low", "none": "None"}
for s in statuses:
    legend_elements.append(Line2D([0], [0], marker=markers[s], color='k',
                                  label=status_labels[s],
                                  markersize=np.sqrt(sizes[s]),
                                  linestyle='None', linewidth=1.5))

ax_legend.legend(handles=legend_elements, loc='center', frameon=False, fontsize=10)

# -------------------------------
# Save figure
# -------------------------------
output_path = "ATAC_tsne_all_peaks_group_status_colored.pdf"
plt.savefig(output_path, dpi=300, bbox_inches='tight')
plt.close()
print(f"Saved combined t-SNE + legend figure with color map: {output_path}")

