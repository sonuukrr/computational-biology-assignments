# Gene Expression Analysis using PCA

This project explores the relationship between two biologically relevant genes, **GATA3** and **XBP1**, using gene expression data and applies **Principal Component Analysis (PCA)** to understand how well these genes separate different sample classes.

The goal is to visualize expression patterns and reduce the data into a simpler form that still captures the main variation between samples.

---

## Project Overview

This script performs two main tasks:

- Visualizes the expression relationship between **GATA3** and **XBP1**
- Applies **PCA** to reduce the 2D gene expression data into a single principal component for easier interpretation

The dataset contains labeled samples categorized into:

- **ER+ (Estrogen Receptor Positive)**
- **ER− (Estrogen Receptor Negative)**

These labels are used to color-code the visualizations for comparison.

---

## Files Required

Make sure the following files are available in the correct directory structure:

```bash
project/
│
├── data/
│   ├── filtered.tsv.gz
│   └── class.tsv
│
└── src/
    └── pca.py
````

### Dataset Description

**filtered.tsv.gz**
Compressed tab-separated dataset containing gene expression values.

**class.tsv**
Contains binary classification labels for each sample:

* `1` → ER+
* `0` → ER−

---

## Genes Used

This analysis focuses on:

* **GATA3** → Column index `4359`
* **XBP1** → Column index `4404`

These genes were selected from the full expression dataset for dimensionality reduction and visualization.

---

## Features

### 1. Scatter Plot of Gene Expression

The script first creates a scatter plot comparing:

* **GATA3 expression**
* **XBP1 expression**

Each point represents one sample.

Color coding:

* **Green** → ER+
* **Red** → ER−

This helps visually inspect whether the two classes show distinct clustering.

---

### 2. Principal Component Analysis (PCA)

The script then performs PCA manually using NumPy.

Steps include:

* Extracting the selected gene expression matrix
* Mean-centering the data
* Computing the covariance matrix
* Performing eigenvalue decomposition
* Sorting eigenvectors by explained variance
* Selecting the first principal component (PC1)
* Projecting all samples onto PC1

This reduces the two-dimensional gene expression space into a single axis that captures the maximum variance.

---

### 3. PC1 Projection Visualization

A second plot shows how all samples distribute along the first principal component.

This makes it easier to observe separation between:

* ER+ samples
* ER− samples

If the classes separate well along PC1, it suggests that these two genes contain meaningful discriminatory information.

---

## Installation

Install the required Python libraries:

```bash
pip install pandas matplotlib numpy
```

---

## How to Run

From the project directory:

```bash
python src/pca.py
```

---

## Output

Running the script generates two visualizations:

1. **GATA3 vs XBP1 scatter plot**
2. **Projection of samples onto PC1**

These plots help interpret gene expression patterns and dimensionality reduction results.

---

## Technologies Used

* Python
* Pandas
* NumPy
* Matplotlib

---

## Learning Outcome

This project is a simple hands-on implementation of PCA in bioinformatics.

It demonstrates:

* Working with gene expression datasets
* Feature selection from large biological datasets
* Manual PCA implementation without machine learning libraries
* Visualization of class-separated biological data

---
