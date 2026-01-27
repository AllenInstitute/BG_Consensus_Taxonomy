import pandas as pd
from sklearn.ensemble import RandomForestClassifier
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42

# ------------------------------
# 1. Read file
# ------------------------------
df = pd.read_csv("./all_peaks/df_for_features_1210.tsv", sep="\t")

# Recode target: high/medium/low -> 1 (on-target), none -> 0 (off-target)
df['status_binary'] = df['status'].apply(lambda x: 0 if x.lower() == 'none' else 1)

# ------------------------------
# 2. Define features
# ------------------------------
exclude_cols = ['region_id', 'status', 'group', 'peak_id', 'status_binary']
feature_cols = [c for c in df.columns if c not in exclude_cols and df[c].dtype != 'object']

X = df[feature_cols]
y = df['status_binary']

# ------------------------------
# 3. Train Random Forest
# ------------------------------
rf = RandomForestClassifier(n_estimators=500, random_state=42)
rf.fit(X, y)

# ------------------------------
# 4. Feature importance
# ------------------------------
importances = pd.Series(rf.feature_importances_, index=feature_cols).sort_values(ascending=False)

# Determine color based on mean difference
colors = []
for f in importances.index:
    mean_on = df.loc[df['status_binary'] == 1, f].mean()
    mean_off = df.loc[df['status_binary'] == 0, f].mean()
    # Use seaborn warm/cool palettes
    if mean_on >= mean_off:
        colors.append(sns.color_palette("coolwarm", as_cmap=False)[0])  # cool for on-target
    else:
        colors.append(sns.color_palette("coolwarm", as_cmap=False)[-1])  # warm for off-target

# ------------------------------
# 5. Horizontal bar plot
# ------------------------------
plt.figure(figsize=(8, 6))
sns.barplot(
    x=importances.values,
    y=importances.index,
    palette=colors,
    orient='h'
)
plt.xlabel('Feature Importance', fontsize=10)
plt.ylabel('')
plt.xticks(fontsize=9)
plt.yticks(fontsize=9)
plt.title('Random Forest Feature Importance\n(Cool = higher on-target, Warm = higher off-target)', fontsize=11)
plt.tight_layout()
plt.savefig("RF_feature_importance_horizontal_seaborn_all_peaks.pdf", dpi=300)
plt.show()

