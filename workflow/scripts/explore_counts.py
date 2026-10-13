import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import linkage, dendrogram
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

expr = pd.read_csv(snakemake.input.tpm, sep="\t", index_col=0)
expr = expr.apply(pd.to_numeric, errors="coerce").fillna(0)

# Keep genes expressed in at least 2 samples
keep = (expr > 1).sum(axis=1) >= 2
expr = expr.loc[keep]
if expr.shape[0] < 10:
    raise ValueError("Too few expressed genes after filtering")

log_expr = np.log2(expr + 1)
log_expr.to_csv(snakemake.output.log_matrix, sep="\t")

# Top variable genes for PCA and correlation
variances = log_expr.var(axis=1).sort_values(ascending=False)
top_n = min(2000, len(variances))
top_genes = variances.head(top_n).index
mat = log_expr.loc[top_genes].T

# Remove zero-variance genes
mat = mat.loc[:, mat.std(axis=0) > 0]
if mat.shape[1] < 2:
    raise ValueError("Too few variable genes after filtering")

scaled = StandardScaler().fit_transform(mat)
n_components = min(5, scaled.shape[0], scaled.shape[1])
pca = PCA(n_components=n_components)
scores = pca.fit_transform(scaled)
variance_ratio = pca.explained_variance_ratio_

pca_df = pd.DataFrame(
    scores,
    index=mat.index,
    columns=[f"PC{i+1}" for i in range(n_components)],
)
pca_df.to_csv(snakemake.output.pca_scores, sep="\t")

fig, ax = plt.subplots(figsize=(7.5, 6.5))
ax.scatter(scores[:, 0], scores[:, 1], s=55, color="#2c7fb8")
for i, sample in enumerate(mat.index):
    ax.annotate(
        sample,
        (scores[i, 0], scores[i, 1]),
        xytext=(4, 3),
        textcoords="offset points",
        fontsize=7,
    )
ax.set_xlabel(f"PC1 ({variance_ratio[0]*100:.1f}%)")
ax.set_ylabel(f"PC2 ({variance_ratio[1]*100:.1f}%)")
ax.set_title("PCA of RNA-seq samples")
ax.axhline(0, color="grey", linewidth=0.5, alpha=0.5)
ax.axvline(0, color="grey", linewidth=0.5, alpha=0.5)
fig.tight_layout()
fig.savefig(snakemake.output.pca, dpi=300)
plt.close(fig)

corr = mat.T.corr(method="pearson")
g = sns.clustermap(
    corr,
    cmap="vlag",
    center=0,
    figsize=(10, 10),
    xticklabels=True,
    yticklabels=True,
)
g.fig.suptitle("Sample correlation (top variable genes)", y=1.02)
g.savefig(snakemake.output.correlation, dpi=300)
plt.close("all")

Z = linkage(scaled, method="ward")
pd.DataFrame(Z).to_csv(
    snakemake.output.linkage,
    sep="\t",
    header=["cluster_1", "cluster_2", "distance", "new_cluster"],
    index=False,
)

fig, ax = plt.subplots(figsize=(10, 5))
dendrogram(Z, labels=mat.index.tolist(), ax=ax)
ax.set_ylabel("Ward distance")
ax.set_title("Hierarchical clustering of RNA-seq samples")
plt.xticks(rotation=90)
fig.tight_layout()
fig.savefig(snakemake.output.clustering, dpi=300)
plt.close(fig)

with open(snakemake.output.summary, "w") as fh:
    fh.write(f"Genes after filtering: {expr.shape[0]}\n")
    fh.write(f"Top variable genes used: {top_n}\n")
    fh.write(f"Samples: {expr.shape[1]}\n")
    fh.write("\nPCA explained variance:\n")
    for i, value in enumerate(variance_ratio, start=1):
        fh.write(f"PC{i}: {value:.4f}\n")
    fh.write("\nNote: Group labels are not available in this run.\n")
