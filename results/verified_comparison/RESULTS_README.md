# Verified Geneformer vs Scanpy Comparison

Dataset: PBMC 3k, 2,671 cells.

- Geneformer V1 10M: 256-dimensional mean-pooled final-layer embeddings.
- Scanpy baseline: 2,000 highly variable genes and 50 PCA components.
- Metric: cosine silhouette score using existing cell-type annotations.
- Geneformer overall silhouette: 0.1075.
- Scanpy PCA overall silhouette: 0.2098.
- Difference (Geneformer minus PCA): -0.1023.

Interpretation:
In this dataset and configuration, the Scanpy PCA representation had
higher overall separation of the supplied cell-type annotations.
This is a descriptive result, not evidence of universal superiority.

Limitations:
- Only one dataset and one Geneformer embedding configuration were tested.
- Rare cell types have small sample sizes.
- Silhouette measures agreement with supplied labels, not general
  biological representation quality.
- UMAP is a visualization and should not be interpreted as a
  quantitative comparison of absolute coordinates.
