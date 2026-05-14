from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR / "data"
RESULTS_DIR = BASE_DIR / "results"
RESULTS_DIR.mkdir(exist_ok=True)


def load_data():
    classes = pd.read_csv(DATA_DIR / "class.tsv", header=None, names=["class"])
    expression = pd.read_csv(DATA_DIR / "filtered.tsv.gz", sep="\t", compression="gzip")
    expression.columns = [str(column).strip() for column in expression.columns]
    probes = pd.read_csv(
        DATA_DIR / "columns.tsv.gz",
        sep="\t",
        comment="#",
        compression="gzip",
    )
    probes["ID"] = probes["ID"].astype(str)
    return classes, expression, probes


def plot_xbp1_gata3(classes, expression, probes):
    selected = probes[probes["GeneSymbol"].isin(["XBP1", "GATA3"])]
    gene_to_id = dict(zip(selected["GeneSymbol"], selected["ID"]))
    missing = {"XBP1", "GATA3"} - set(gene_to_id)
    if missing:
        raise ValueError(f"Missing gene IDs for: {', '.join(sorted(missing))}")

    points = pd.DataFrame(
        {
            "XBP1": expression[gene_to_id["XBP1"]],
            "GATA3": expression[gene_to_id["GATA3"]],
            "class": classes["class"],
        }
    )
    points.to_csv(RESULTS_DIR / "xbp1_gata3_expression.tsv", sep="\t", index=False)

    colors = points["class"].map({0: "#2f6fbb", 1: "#c9463d"})
    labels = {0: "ER- breast cancer", 1: "ER+ breast cancer"}

    fig, ax = plt.subplots(figsize=(6.2, 5.2), dpi=180)
    for klass, group in points.groupby("class"):
        ax.scatter(
            group["XBP1"],
            group["GATA3"],
            s=38,
            alpha=0.82,
            color=colors[group.index].iloc[0],
            edgecolor="white",
            linewidth=0.45,
            label=labels[klass],
        )
    ax.set_xlabel("XBP1 expression")
    ax.set_ylabel("GATA3 expression")
    ax.set_title("XBP1 vs GATA3 expression in 105 breast cancer samples")
    ax.legend(frameon=False, loc="best")
    ax.grid(alpha=0.18)
    fig.tight_layout()
    fig.savefig(RESULTS_DIR / "xbp1_gata3_scatter.png")
    plt.close(fig)
    return points


def run_pca(points):
    features = points[["XBP1", "GATA3"]].to_numpy(dtype=float)
    centered = features - features.mean(axis=0)
    scaled = centered / features.std(axis=0, ddof=1)

    _, singular_values, vt = np.linalg.svd(scaled, full_matrices=False)
    scores = scaled @ vt.T
    explained = singular_values**2 / np.sum(singular_values**2)

    result = pd.DataFrame(
        {
            "PC1": scores[:, 0],
            "PC2": scores[:, 1],
            "class": points["class"],
        }
    )
    result.to_csv(RESULTS_DIR / "pca_scores.tsv", sep="\t", index=False)

    rng = np.random.default_rng(102)
    y_jitter = result["class"].map({0: -0.08, 1: 0.08}).to_numpy()
    y_jitter = y_jitter + rng.normal(0, 0.018, len(result))
    colors = result["class"].map({0: "#2f6fbb", 1: "#c9463d"})

    fig, ax = plt.subplots(figsize=(6.8, 2.9), dpi=180)
    for klass, label in [(0, "ER- breast cancer"), (1, "ER+ breast cancer")]:
        mask = result["class"] == klass
        ax.scatter(
            result.loc[mask, "PC1"],
            y_jitter[mask],
            s=36,
            color=colors[mask].iloc[0],
            alpha=0.82,
            edgecolor="white",
            linewidth=0.4,
            label=label,
        )
    ax.axhline(0, color="#303030", linewidth=0.9, alpha=0.45)
    ax.set_yticks([])
    ax.set_xlabel(f"PC1 ({explained[0] * 100:.1f}% variance)")
    ax.set_title("Projection of XBP1/GATA3 expression onto PC1")
    ax.legend(frameon=False, loc="upper right")
    ax.grid(axis="x", alpha=0.18)
    fig.tight_layout()
    fig.savefig(RESULTS_DIR / "pc1_projection.png")
    plt.close(fig)


def main():
    classes, expression, probes = load_data()
    points = plot_xbp1_gata3(classes, expression, probes)
    run_pca(points)
    print(f"Wrote PCA results to {RESULTS_DIR}")


if __name__ == "__main__":
    main()
