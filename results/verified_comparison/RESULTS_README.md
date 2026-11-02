# Verified Geneformer vs Scanpy Comparison

## Dataset and methods

- Dataset: PBMC 3k.
- Cells: 2,671.
- Cell types: 6.
- Geneformer: V1 10M, 256-dimensional mean-pooled second-to-last-layer embeddings.
- Scanpy baseline: 2,000 highly variable genes and 50 PCA components.
- Metric: cosine silhouette.
- Annotation reference: canonical `adata_annotated.h5ad`.

## Overall results

| Representation | Dimensions | Cosine silhouette |
|---|---:|---:|
| Geneformer V1 10M | 256 | 0.107173 |
| Scanpy PCA | 50 | 0.209903 |

Difference (Geneformer minus PCA): **-0.102731**.

In this dataset and under the tested configurations, Scanpy PCA
achieved higher overall separation of the supplied cell-type
annotations than Geneformer.

## Cell-type-specific results

| Cell type | Cells | Geneformer mean | PCA mean |
|---|---:|---:|---:|
| B cells | 347 | -0.035509 | 0.353437 |
| Candidate dendritic cells | 33 | 0.513240 | 0.485060 |
| Cytotoxic lymphocytes / NK-associated | 455 | 0.081706 | 0.204539 |
| Monocytes / myeloid | 638 | -0.017298 | 0.223300 |
| Platelet-associated | 15 | 0.429695 | 0.738441 |
| T cells | 1,183 | 0.210530 | 0.148263 |

These are descriptive mean silhouette values. Rare populations,
particularly platelet-associated cells and candidate dendritic cells,
have small sample sizes and require cautious interpretation.

## Annotation provenance audit

The canonical annotations were read from
`~/scrna-pbmc3k/data/adata_annotated.h5ad` and aligned to both
feature matrices by cell barcode.

All 2,671 barcodes matched in both set and order. The saved comparison
matrices had two cell-type assignments that differed from the
canonical annotations:

- `CCCTTACTGCAGTT-1`: saved cytotoxic/NK-associated; canonical T cells.
- `GAACCAACCACAAC-1`: saved cytotoxic/NK-associated; canonical T cells.

The remaining 2,669 cell-type assignments agreed. The original saved
cell-type summary is retained separately for audit purposes.

## Interpretation and limitations

The overall result is descriptive and applies only to this dataset
and these tested configurations. Silhouette measures geometric
separation relative to supplied labels; it does not independently
validate biological cell identities or establish universal
representation quality.

The six annotations are not independent biological ground truth.
Rare populations have small sample sizes. UMAP is used for
visualization and is not a quantitative comparison of absolute
coordinates.

The canonical metrics were recalculated from the saved feature
matrices using canonical cell-type labels and cosine silhouette.
