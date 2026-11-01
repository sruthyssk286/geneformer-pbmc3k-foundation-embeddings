#!/usr/bin/env python
"""
Reproducible Geneformer V1 embedding workflow for the verified PBMC 3k dataset.

This script preserves the successful pandas-3-compatible tokenizer patch from
the original Colab notebook. It expects:
  1) an Ensembl-prepared AnnData file for tokenization;
  2) the verified AnnData file used to validate barcodes and labels;
  3) a local Geneformer checkout, V1-10M model directory, and V1 dictionaries.

Important:
- Run in an environment where the Geneformer dependencies are installed.
- The original notebook used pandas 3.x and required this runtime patch.
- Review the Geneformer model/dictionary license and version requirements.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from datasets import load_from_disk
from geneformer import EmbExtractor, TranscriptomeTokenizer
from geneformer.tokenizer import rank_genes, sum_ensembl_ids


def install_pandas_compatible_tokenizer_patch() -> None:
    """Patch positional Series indexing used by Geneformer's tokenizer."""
    def tokenize_anndata_fixed(
        self, adata_file_path, target_sum=10_000, file_format="h5ad"
    ):
        adata = sum_ensembl_ids(
            adata_file_path,
            self.collapse_gene_ids,
            self.gene_mapping_dict,
            self.gene_token_dict,
            self.custom_attr_name_dict,
            self.use_h5ad_index,
            file_format=file_format,
            chunk_size=self.chunk_size,
        )

        file_cell_metadata = (
            {key: [] for key in self.custom_attr_name_dict}
            if self.custom_attr_name_dict is not None
            else None
        )

        coding_miRNA_loc = np.where([
            self.genelist_dict.get(gene, False)
            for gene in adata.var["ensembl_id_collapsed"]
        ])[0]

        collapsed_ids = adata.var["ensembl_id_collapsed"]
        coding_ids = collapsed_ids.iloc[coding_miRNA_loc].tolist()

        norm_factor_vector = np.array([
            self.gene_median_dict[gene] for gene in coding_ids
        ])
        coding_miRNA_tokens = np.array([
            self.gene_token_dict[gene] for gene in coding_ids
        ])

        if "filter_pass" in adata.obs.columns:
            filter_pass_loc = np.where(
                adata.obs["filter_pass"].to_numpy() == 1
            )[0]
        else:
            print(
                f"{adata_file_path} has no 'filter_pass' column; "
                "tokenizing all cells."
            )
            filter_pass_loc = np.arange(adata.n_obs)

        tokenized_cells = []
        tokenized_counts = []

        for start in range(0, len(filter_pass_loc), self.chunk_size):
            idx = filter_pass_loc[start:start + self.chunk_size]
            batch = adata[idx]
            n_counts = batch.obs["n_counts"].to_numpy()[:, None]
            X_view = batch.X[:, coding_miRNA_loc]

            # Guard against zero-count cells rather than dividing by zero.
            if np.any(n_counts <= 0):
                raise ValueError(
                    "Found cells with non-positive n_counts during tokenization."
                )

            X_norm_unscaled = X_view / n_counts * target_sum
            X_norm = sp.csr_matrix(X_norm_unscaled / norm_factor_vector)
            X_norm_unscaled = sp.csr_matrix(X_norm_unscaled)

            for i in range(X_norm.shape[0]):
                tokenized_cells.append(
                    rank_genes(
                        X_norm[i].data,
                        coding_miRNA_tokens[X_norm[i].indices],
                    )
                )
                if self.keep_counts:
                    tokenized_counts.append(
                        rank_genes(
                            X_norm[i].data,
                            X_norm_unscaled[i].data,
                        )
                    )

            if file_cell_metadata is not None:
                for key in file_cell_metadata:
                    file_cell_metadata[key].extend(batch.obs[key].tolist())

        empty_indices = [
            i for i, sequence in enumerate(tokenized_cells)
            if len(sequence) == 0
        ]
        for i in reversed(empty_indices):
            del tokenized_cells[i]
            if self.keep_counts:
                del tokenized_counts[i]
            if file_cell_metadata is not None:
                for key in file_cell_metadata:
                    del file_cell_metadata[key][i]

        return tokenized_cells, file_cell_metadata, tokenized_counts

    TranscriptomeTokenizer.tokenize_anndata = tokenize_anndata_fixed
    print("Installed pandas-compatible tokenizer patch.")


def validate_tokenized_dataset(dataset_path: Path, expected_cells: int) -> None:
    ds = load_from_disk(str(dataset_path))
    required = {"input_ids", "cell_barcode", "leiden", "cell_type"}
    missing = required - set(ds.column_names)
    if missing:
        raise ValueError(f"Tokenized dataset is missing columns: {sorted(missing)}")
    if len(ds) != expected_cells:
        raise ValueError(
            f"Tokenized cell count ({len(ds)}) does not match input "
            f"({expected_cells})."
        )
    barcodes = ds["cell_barcode"]
    if any(value is None for value in barcodes):
        raise ValueError("Tokenized dataset contains missing cell barcodes.")
    if len(set(barcodes)) != len(barcodes):
        raise ValueError("Tokenized cell barcodes are not unique.")
    print(f"Tokenized dataset validated: {len(ds)} cells; columns={ds.column_names}")


def validate_embeddings(
    embedding_csv: Path,
    validation_h5ad: Path,
    output_dir: Path,
) -> dict:
    emb = pd.read_csv(embedding_csv)
    adata = sc.read_h5ad(validation_h5ad)
    embedding_cols = [str(i) for i in range(256)]
    required = {"cell_barcode", "leiden", "cell_type", *embedding_cols}
    missing = required - set(emb.columns)
    if missing:
        raise ValueError(f"Embedding CSV is missing columns: {sorted(missing)}")

    if not emb["cell_barcode"].is_unique or not adata.obs_names.is_unique:
        raise ValueError("Cell barcodes must be unique in both files.")
    if set(emb["cell_barcode"].astype(str)) != set(adata.obs_names.astype(str)):
        raise ValueError("Embedding barcodes do not match validation AnnData.")

    expected = adata.obs[["leiden", "cell_type"]].copy()
    expected.index = adata.obs_names.astype(str)
    observed = emb.copy()
    observed["cell_barcode"] = observed["cell_barcode"].astype(str)
    observed = observed.set_index("cell_barcode")
    expected = expected.loc[observed.index]

    if not np.array_equal(
        observed["leiden"].astype(str).to_numpy(),
        expected["leiden"].astype(str).to_numpy(),
    ):
        raise ValueError("Leiden metadata mismatch by barcode.")
    if not np.array_equal(
        observed["cell_type"].astype(str).to_numpy(),
        expected["cell_type"].astype(str).to_numpy(),
    ):
        raise ValueError("Cell-type metadata mismatch by barcode.")

    X = observed[embedding_cols].to_numpy(dtype=np.float64)
    if not np.isfinite(X).all():
        raise ValueError("Embedding matrix contains NaN or infinite values.")
    norms = np.linalg.norm(X, axis=1)
    if np.any(norms == 0):
        raise ValueError("Embedding matrix contains zero vectors.")
    if np.any(X.var(axis=0) == 0):
        raise ValueError("At least one embedding dimension has zero variance.")

    summary = {
        "n_cells": int(X.shape[0]),
        "embedding_dimensions": int(X.shape[1]),
        "all_values_finite": bool(np.isfinite(X).all()),
        "barcode_match": True,
        "leiden_match": True,
        "cell_type_match": True,
        "zero_vectors": int(np.sum(norms == 0)),
        "dimensions_with_zero_variance": int(np.sum(X.var(axis=0) == 0)),
        "duplicate_vectors": int(X.shape[0] - np.unique(X, axis=0).shape[0]),
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    summary_path = output_dir / "embedding_validation.json"
    import json
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print("Embedding validation passed:")
    print(json.dumps(summary, indent=2))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Tokenize PBMC AnnData and extract Geneformer V1 embeddings."
    )
    parser.add_argument("--input-dir", type=Path, required=True,
                        help="Directory containing the Ensembl-prepared .h5ad.")
    parser.add_argument("--input-identifier", default="pbmc3k_verified_ensembl_ready")
    parser.add_argument("--validation-h5ad", type=Path, required=True,
                        help="Verified symbol/metadata AnnData used for validation.")
    parser.add_argument("--geneformer-root", type=Path, required=True,
                        help="Path to the Geneformer checkout.")
    parser.add_argument("--model-dir", type=Path, required=True,
                        help="Path to Geneformer-V1-10M.")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--forward-batch-size", type=int, default=16)
    args = parser.parse_args()

    args = parser.parse_args()

    args.output_dir = args.output_dir.resolve()

    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        raise FileExistsError(
            f"Output directory is not empty: {args.output_dir}\n"
            "Choose a new output directory to avoid overwriting existing results."
        )

    gene_dict_dir = (
        args.geneformer_root / "geneformer" / "gene_dictionaries_30m"
    )
    dictionary_files = {
        "token dictionary": gene_dict_dir / "token_dictionary_gc30M.pkl",
        "gene median dictionary": gene_dict_dir / "gene_median_dictionary_gc30M.pkl",
        "gene mapping dictionary": gene_dict_dir / "ensembl_mapping_dict_gc30M.pkl",
    }
    dictionary_files = {
        "token dictionary": gene_dict_dir / "token_dictionary_gc30M.pkl",
        "gene median dictionary": gene_dict_dir / "gene_median_dictionary_gc30M.pkl",
        "gene mapping dictionary": gene_dict_dir / "ensembl_mapping_dict_gc30M.pkl",
    }
    for description, path in dictionary_files.items():
        if not path.is_file():
            raise FileNotFoundError(f"Missing {description}: {path}")
    if not args.model_dir.is_dir():
        raise FileNotFoundError(f"Model directory not found: {args.model_dir}")
    if not args.validation_h5ad.is_file():
        raise FileNotFoundError(f"Validation AnnData not found: {args.validation_h5ad}")

    source_files = list(args.input_dir.glob("*.h5ad"))
    target_input = args.input_dir / f"{args.input_identifier}.h5ad"
    if not target_input.is_file():
        raise FileNotFoundError(
            f"Expected tokenization input {target_input}; found "
            f"{[p.name for p in source_files]}"
        )

    verified = sc.read_h5ad(args.validation_h5ad)
    expected_cells = verified.n_obs
    if not {"leiden", "cell_type"}.issubset(verified.obs.columns):
        raise ValueError("Validation AnnData must contain leiden and cell_type.")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    tokenized_dir = args.output_dir / "tokenized"
    embedding_dir = args.output_dir / "embeddings"
    tokenized_dir.mkdir(parents=True, exist_ok=True)
    embedding_dir.mkdir(parents=True, exist_ok=True)

    install_pandas_compatible_tokenizer_patch()

    tokenizer = TranscriptomeTokenizer(
        custom_attr_name_dict={
            "cell_barcode": "cell_barcode",
            "leiden": "leiden",
            "cell_type": "cell_type",
        },
        nproc=1,
        model_version="V1",
        use_h5ad_index=True,
        token_dictionary_file=str(dictionary_files["token dictionary"]),
        gene_median_file=str(dictionary_files["gene median dictionary"]),
        gene_mapping_file=str(dictionary_files["gene mapping dictionary"]),
    )

    prefix = "pbmc3k_verified_ensembl_v1"
    tokenizer.tokenize_data(
        data_directory=str(args.input_dir),
        output_directory=str(tokenized_dir),
        output_prefix=prefix,
        file_format="h5ad",
        input_identifier=args.input_identifier,
    )
    dataset_path = tokenized_dir / f"{prefix}.dataset"
    if not dataset_path.exists():
        raise FileNotFoundError(f"Tokenizer did not create {dataset_path}")
    validate_tokenized_dataset(dataset_path, expected_cells)

    extractor = EmbExtractor(
        model_type="Pretrained",
        emb_mode="cell",
        cell_emb_style="mean_pool",
        emb_layer=-1,
        emb_label=["cell_barcode", "leiden", "cell_type"],
        max_ncells=None,
        forward_batch_size=args.forward_batch_size,
        nproc=1,
        model_version="V1",
    )
    extractor.extract_embs(
        model_directory=str(args.model_dir),
        input_data_file=str(dataset_path),
        output_directory=str(embedding_dir),
        output_prefix="pbmc3k_verified_geneformer_v1",
        output_torch_embs=True,
    )

    embedding_csv = embedding_dir / "pbmc3k_verified_geneformer_v1.csv"
    if not embedding_csv.is_file():
        raise FileNotFoundError(f"Embedding CSV not found: {embedding_csv}")
    validate_embeddings(embedding_csv, args.validation_h5ad, args.output_dir)


if __name__ == "__main__":
    main()
