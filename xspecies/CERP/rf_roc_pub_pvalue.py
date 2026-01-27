import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_curve, auc
from scipy.ndimage import gaussian_filter1d
from scipy.stats import wilcoxon
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42

# -----------------------
# Load data
# -----------------------
df = pd.read_csv("all_peaks/df_for_features_1210.tsv", sep="\t")

# Map status to binary: high/medium/low -> 1, others -> 0
df['status_binary'] = df['status'].apply(lambda x: 1 if x.lower() in ['high', 'medium', 'low'] else 0)

# -----------------------
# Define features
# -----------------------
all_features = ['GC_content', 'PhyloP', 'TSS_distance', 'gene_tau_score',
                'CREsted_specificity', 'ATAC_specificity',
                'skewness', 'kurtosis', 'bimodality']

X_all = df[all_features].values
X_single = df[['ATAC_specificity']].values
y = df['status_binary'].values

# -----------------------
# Cross-validation setup
# -----------------------
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=35)

def compute_cv_roc(X, y, cv):
    """Compute mean ROC curve and per-fold AUCs using cross-validation."""
    tprs = []
    aucs = []
    mean_fpr = np.linspace(0, 1, 100)

    for train_idx, test_idx in cv.split(X, y):
        X_train, X_test = X[train_idx], X[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]

        clf = RandomForestClassifier(n_estimators=100, random_state=42)
        clf.fit(X_train, y_train)
        y_proba = clf.predict_proba(X_test)[:, 1]

        fpr, tpr, _ = roc_curve(y_test, y_proba)
        auc_score = auc(fpr, tpr)
        aucs.append(auc_score)

        tpr_interp = np.interp(mean_fpr, fpr, tpr)
        tpr_interp[0] = 0.0
        tprs.append(tpr_interp)

    mean_tpr = np.mean(tprs, axis=0)
    mean_tpr[-1] = 1.0
    mean_auc = np.mean(aucs)
    std_auc = np.std(aucs)

    return mean_fpr, mean_tpr, mean_auc, std_auc, aucs

# -----------------------
# Compute ROC curves
# -----------------------
fpr_single, tpr_single, mean_auc_single, std_auc_single, aucs_single = compute_cv_roc(X_single, y, cv)
fpr_all, tpr_all, mean_auc_all, std_auc_all, aucs_all = compute_cv_roc(X_all, y, cv)

# -----------------------
# Statistical test: Wilcoxon signed-rank
# -----------------------
stat, p_value = wilcoxon(aucs_all, aucs_single, alternative='greater')

# -----------------------
# Print results
# -----------------------
print("========================================")
print("Cross-validation AUC comparison")
print("========================================")
print(f"Single-feature (ATAC specificity): {mean_auc_single:.3f} ± {std_auc_single:.3f}")
print(f"All-features Random Forest:        {mean_auc_all:.3f} ± {std_auc_all:.3f}")
print("----------------------------------------")
print("Per-fold AUCs:")
for i, (a1, a2) in enumerate(zip(aucs_single, aucs_all), 1):
    print(f"Fold {i}: single = {a1:.3f}, all = {a2:.3f}, diff = {a2 - a1:.3f}")
print("----------------------------------------")
print(f"Wilcoxon signed-rank test (H1: all > single): p = {p_value:.4f}")
print("========================================")

# -----------------------
# Smooth TPR curves
# -----------------------
tpr_single_smooth = gaussian_filter1d(tpr_single, sigma=0.2)
tpr_all_smooth = gaussian_filter1d(tpr_all, sigma=0.2)

# -----------------------
# Plot ROC curves
# -----------------------
plt.figure(figsize=(6,6))

# single feature: grey
plt.plot(fpr_single, tpr_single_smooth, color='grey', lw=2, label='ATAC specificity')

# all features: black
plt.plot(fpr_all, tpr_all_smooth, color='black', lw=2, label='All features')

# Diagonal reference line
plt.plot([0, 1], [0, 1], 'k--', alpha=0.3)

# AUC text labels
plt.text(1.02, tpr_single_smooth[-1], f"AUC = {mean_auc_single:.2f} ± {std_auc_single:.2f}", 
         color='grey', fontsize=10, va='bottom', ha='left')
plt.text(1.02, tpr_all_smooth[-1]-0.05, f"AUC = {mean_auc_all:.2f} ± {std_auc_all:.2f}", 
         color='black', fontsize=10, va='bottom', ha='left')

# Add p-value text on the plot
plt.text(0.6, 0.2, f"Wilcoxon p = {p_value:.3e}", fontsize=10, ha='left')

plt.xlim([0.0, 1.1])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate', fontsize=12)
plt.ylabel('True Positive Rate', fontsize=12)
plt.title('ROC Curve with Cross-Validation', fontsize=14)
plt.legend(loc='lower right', fontsize=10)
plt.grid(alpha=0.3)

plt.tight_layout()
plt.savefig("ROC_curve_pub_pvalue.pdf", format='pdf', dpi=300)
plt.show()

