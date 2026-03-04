# Essential-genes

A small machine-learning showcase using *Bacillus subtilis* essential gene data from [SubtiWiki](https://subtiwiki.uni-goettingen.de/).

## Overview

This project demonstrates two fundamental machine-learning techniques applied to a biology dataset:

| Technique | Target variable | Task |
|---|---|---|
| **Linear Regression** | `expression_level` | Predict continuous gene expression from gene length, GC content and interaction count |
| **Logistic Regression** | `is_essential` | Classify genes as essential (1) or non-essential (0) |

## Dataset

`data/bacillus_subtilis_genes.csv` – 78 genes (40 essential, 38 non-essential) with the following columns:

| Column | Description |
|---|---|
| `gene` | Gene name |
| `gene_length_bp` | Coding sequence length in base-pairs |
| `gc_content` | Fraction of G+C nucleotides |
| `expression_level` | Log-scale mRNA expression level |
| `num_interactions` | Number of known protein–protein interactions |
| `is_essential` | 1 = essential, 0 = non-essential |

Gene essentiality annotations are based on transposon-insertion and deletion studies catalogued in SubtiWiki.

## Requirements

```
pandas
scikit-learn
matplotlib
```

Install with:

```bash
pip install pandas scikit-learn matplotlib
```

## Usage

```bash
python analysis.py
```

The script prints model metrics to the terminal and saves two figures:

- `linear_regression_results.png` – actual vs. predicted expression levels
- `logistic_regression_confusion_matrix.png` – confusion matrix for essentiality classification

## Results

**Linear Regression**

- R² ≈ 0.90
- Number of protein–protein interactions is the strongest predictor of expression level.

**Logistic Regression**

- Accuracy ≈ 94 %
- Gene length and interaction count are the most informative features for distinguishing essential from non-essential genes.
