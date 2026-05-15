import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# =========================
# LOAD DATA
# =========================

df = pd.read_csv("../data/filtered.tsv.gz",
                 sep="\t",
                 compression="gzip")

# select GATA3 and XBP1 columns
# according to README:
# 4359 -> GATA3
# 4404 -> XBP1

df = df.iloc[:, [4359, 4404]]

df.columns = ["GATA3", "XBP1"]

# labels
labels = pd.read_csv("../data/class.tsv",
                     header=None,
                     names=["result"])

df["result"] = labels["result"]

# =========================
# PLOT 1
# GATA3 vs XBP1
# =========================

plt.figure(figsize=(6,6))

er_pos = df[df["result"] == 1]
er_neg = df[df["result"] == 0]

plt.scatter(er_pos["GATA3"],
            er_pos["XBP1"],
            c="green",
            label="ER+")

plt.scatter(er_neg["GATA3"],
            er_neg["XBP1"],
            c="red",
            label="ER-")

plt.xlabel("GATA3 Expression")
plt.ylabel("XBP1 Expression")

plt.title("GATA3 vs XBP1")

plt.legend()
plt.tight_layout()
plt.show()

# =========================
# PCA
# =========================

mat = df[["GATA3", "XBP1"]].values

# center data
mat_centered = mat - mat.mean(axis=0)

# covariance matrix
cov = np.cov(mat_centered, rowvar=False)

# eigendecomposition
eigenvalues, eigenvectors = np.linalg.eig(cov)

# sort by descending eigenvalue
sort = np.argsort(eigenvalues)[::-1]

eigenvalues = eigenvalues[sort]
eigenvectors = eigenvectors[:, sort]

# PC1
pc1 = eigenvectors[:, 0]

# projection onto PC1
pc1_scores = mat_centered @ pc1

# =========================
# PLOT 2
# Projection onto PC1
# =========================

pc1_er_pos = pc1_scores[df["result"] == 1]
pc1_er_neg = pc1_scores[df["result"] == 0]

plt.figure(figsize=(8,3))

plt.scatter(pc1_er_pos,
            np.zeros_like(pc1_er_pos),
            c="green",
            label="ER+")

plt.scatter(pc1_er_neg,
            np.zeros_like(pc1_er_neg),
            c="red",
            label="ER-")

plt.xlabel("Projection onto PC1")
plt.yticks([])

plt.title("Projection of Samples onto PC1")

plt.legend()
plt.tight_layout()
plt.show()
