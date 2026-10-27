# Geneformer Reproducibility Notes

## Scope

This document records the Geneformer configuration used to generate the
PBMC single-cell embeddings in this repository. It distinguishes the
configuration observed during the successful Colab run from steps that
still need to be packaged and tested as a clean, end-to-end workflow.

The WSL `requirements.txt` file covers validation and visualization of
saved comparison results. It is not the dependency specification for
running Geneformer tokenization and embedding extraction.

## Dataset and input

- Dataset: PBMC 3k single-cell RNA-seq.
- Cells analyzed: 2,671.
- Input matrix dimensions: 2,671 cells by 35,635 genes.
- Input AnnData file used in Colab:
  `pbmc3k_verified_ensembl.h5ad`.
- Gene identifiers were mapped to Ensembl identifiers for Geneformer
  tokenization.
- Cell barcodes, Leiden labels, and preliminary cell-type annotations
  were retained for output validation and comparison.

The input file was stored in the Colab runtime at
`/content/geneformer_input/pbmc3k_verified_ensembl.h5ad`.
The source dataset download and all preprocessing steps must be
documented or packaged separately to reproduce this file from scratch.

## Model and dictionaries

- Model family: Geneformer V1.
- Model checkpoint directory: `Geneformer-V1-10M`.
- Model directory used in Colab: `/content/Geneformer/Geneformer-V1-10M`.
- Dictionary directory:
  `/content/Geneformer/geneformer/gene_dictionaries_30m/`.
- Dictionary files:
  - `token_dictionary_gc30M.pkl`
  - `gene_median_dictionary_gc30M.pkl`
  - `ensembl_mapping_dict_gc30M.pkl`

These filenames and runtime paths record the successful run. A reproducible
release should additionally record the Geneformer Git commit or release,
the exact model and dictionary source revisions, and checksums where
available. The model weights and dictionary files are not included in
this repository.

## Software environment

The successful Colab workflow used Python 3.12.13 and
`transformers==4.44.0`. The workflow also used Geneformer, Scanpy,
Hugging Face Datasets, and other supporting packages.

The project history records the following package versions during the
Colab work:

- `scanpy==1.12.4`
- `transformers==4.44.0`
- `huggingface_hub==0.36.2`
- `datasets==4.0.0`
- `loompy==3.0.8`

`tdigest` and `optuna` were installed during troubleshooting, but their
final versions were not recorded in the available run notes.

This list is a record of reported versions, not a complete tested
environment lockfile. Exact versions of all direct and transitive
dependencies, the PyTorch build, CUDA runtime, and the Geneformer source
revision still need to be captured from a successful runtime.

## Tokenization compatibility workaround

During tokenization, the installed pandas version caused an indexing
problem in Geneformer's `TranscriptomeTokenizer` implementation when
median expression values were indexed. The successful run used a
runtime-only compatibility patch to use `.iloc` for positional indexing
of the median values.

The patch was applied in the Colab runtime and was not saved as a
version-controlled source file. It will therefore not be applied
automatically by this repository.

Before publishing a fully reproducible tokenization script:

1. Inspect the exact Geneformer tokenizer implementation and pandas
   version in the intended environment.
2. Prefer a compatible, documented dependency combination if available.
3. If a patch remains necessary, implement it explicitly in a maintained
   script or documented compatibility layer.
4. Test tokenization from a fresh runtime and verify the output count
   and metadata alignment.

Do not assume that the runtime patch is present after restarting Colab.

## Embedding extraction configuration

The successful embedding extraction used:

- Model version: `V1`.
- Embedding mode: `cell`.
- Cell embedding style: `mean_pool`.
- Embedding layer: `-1` (final layer).
- Batch size: `16`.

Output file:
`pbmc3k_verified_geneformer_v1.csv`.

The exported table contained 2,671 cells and 256 embedding features,
along with the cell barcode, Leiden label, and cell-type annotation.

The successful output was checked for unique barcodes, barcode-set
agreement with the input, metadata alignment, finite values, missing
values, zero vectors, duplicate vectors, and zero-variance dimensions.

## Saved comparison results

Validated comparison artifacts are stored in
`results/verified_comparison/`.

The reported overall cosine silhouette scores were:

| Representation | Dimensions | Mean cosine silhouette |
|---|---:|---:|
| Geneformer V1 mean-pooled embeddings | 256 | 0.1075 |
| Scanpy PCA baseline | 50 | 0.2098 |

Difference (Geneformer minus PCA): `-0.1023`.

These are descriptive results for the tested configurations and supplied
cell-type annotations. They do not establish that one representation is
universally better. Rare cell populations have small sample sizes, and
the annotations are not independent biological ground truth.

## Current reproducibility status

The saved-results validation notebook can be executed against the
committed comparison artifacts. The complete pipeline from raw input
through preprocessing, Geneformer tokenization, embedding extraction,
and comparison has not yet been packaged and tested in a clean
environment.

Before claiming full end-to-end reproducibility, the project still needs:

- A documented source and retrieval procedure for the dataset.
- Reproducible preprocessing and Ensembl mapping steps.
- The Geneformer source revision and exact model/dictionary revisions.
- A tested dependency specification for the Geneformer runtime.
- A maintained solution to the tokenizer compatibility issue.
- An executable tokenization and embedding extraction script or notebook.
- A clean-runtime test and output validation procedure.
