# Essential-genes

A small machine-learning showcase using *Bacillus subtilis* essential gene data from [SubtiWiki](https://subtiwiki.uni-goettingen.de/).

## Overview

This project demonstrates two fundamental machine-learning techniques applied to a biology dataset:

| Technique | Target variable | Task |
|---|---|---|
| **Linear Regression** | `expression_level` | Predict continuous gene expression from gene length, GC content and interaction count |
| **Logistic Regression** | `is_essential` | Classify genes as essential (1) or non-essential (0) |

## Dataset

`data/bacillus_subtilis_genes.csv` contains 78 *B. subtilis* 168 genes with the following columns:

| Column | Description | Source |
|---|---|---|
| `gene` | Gene name | SubtiWiki / genome annotation |
| `locus_tag` | BSU locus tag | SubtiWiki |
| `gene_length_bp` | CDS length in base-pairs | Derived from genome positions |
| `gc_content` | Fraction of G+C nucleotides in the CDS | Computed from genome sequence |
| `expression_level` | Mean log₂ mRNA expression level | SubtiWiki expression compendium |
| `num_interactions` | Number of curated protein–protein interactions | SubtiWiki interaction data |
| `is_essential` | 1 = essential, 0 = non-essential | SubtiWiki "Essential genes" category |

> **Note on the bundled CSV:** The gene names, locus tags, and essentiality labels are real and sourced from published B. subtilis research catalogued in SubtiWiki (Kobayashi et al. 2003; Koo et al. 2017). The numeric features (`gc_content`, `expression_level`, `num_interactions`) in the bundled file are **representative/illustrative values** based on the B. subtilis 168 reference genome and published literature, not values downloaded live from SubtiWiki. To replace the CSV with genuine SubtiWiki data, run `fetch_data.py` (see below).

## Fetching Real SubtiWiki Data

`fetch_data.py` queries SubtiWiki's REST API v4 and NCBI to download genuine data and regenerate the CSV:

```bash
pip install requests biopython pandas
python fetch_data.py
```

This will:
1. Fetch all gene positions from SubtiWiki (`/v4/api/gene/`)
2. Identify essential genes via SubtiWiki's "Essential genes" category
3. Compute per-gene GC content from the NCBI genome (accession AL009126.3)
4. Retrieve expression levels and interaction counts from SubtiWiki
5. Write the result to `data/bacillus_subtilis_genes.csv`

## Running the Analysis

```bash
pip install pandas scikit-learn matplotlib
python analysis.py
```

The script prints model metrics to the terminal and saves two figures:

- `linear_regression_results.png` – actual vs. predicted expression levels
- `logistic_regression_confusion_matrix.png` – confusion matrix for essentiality classification

## Results (on bundled representative dataset)

**Linear Regression**

- R² ≈ 0.90
- Number of protein–protein interactions is the strongest predictor of expression level.

**Logistic Regression**

- Accuracy ≈ 94 %
- Gene length and interaction count are the most informative features for distinguishing essential from non-essential genes.
