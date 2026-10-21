# PBMC 3K Geneformer Foundation Model Analysis

## Objective
Explore pretrained Geneformer V1 cell embeddings using a PBMC 3K single-cell RNA-seq dataset.

## Dataset
- Cells: 2,671
- Original Leiden clusters: 6
- Geneformer embedding dimensions: 256

## Workflow completed
1. Prepared the AnnData input for Geneformer.
2. Tokenized cells using Geneformer V1 resources.
3. Extracted pretrained cell embeddings.
4. Validated the embedding matrix and Leiden labels.
5. Generated a UMAP visualization of Geneformer embeddings.
6. Examined marker genes and refined preliminary cluster annotations.
7. Calculated exploratory silhouette scores.

## Preliminary results

| Representation | Dimensions | Distance | Silhouette |
|---|---:|---|---:|
| Geneformer V1 | 256 | Cosine | 0.1015 |
| Reconstructed PCA baseline | 30 | Cosine | 0.3044 |

## Preliminary biological annotations

| Leiden cluster | Cells | Annotation |
|---|---:|---|
| 0 | 1,183 | T cells |
| 1 | 638 | Myeloid cells / monocytes |
| 2 | 455 | NK / cytotoxic lymphocytes |
| 3 | 347 | B cells |
| 4 | 33 | Antigen-presenting cells; unresolved |
| 5 | 15 | Platelets |

## Limitations
- The PCA baseline was reconstructed and may not reproduce the original clustering representation exactly.
- The biological annotations were derived from the existing Leiden clusters and marker genes.
- Because annotations map one-to-one to Leiden labels, their silhouette score is not independent validation.
- The silhouette results are exploratory and do not establish that one representation is biologically superior.
- Cluster 2 requires further work to distinguish NK cells from cytotoxic T cells; cluster 4 remains unresolved.

## Next steps
- Evaluate cell neighborhoods in Geneformer embedding space.
- Investigate marker-informed biological coherence independently of the Leiden labels.
- Document the software versions and preprocessing decisions.
- Add the notebook, visualizations, and results to the GitHub repository.
