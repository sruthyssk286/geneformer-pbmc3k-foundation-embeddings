
#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$HOME/geneformer-pbmc3k-foundation-embeddings"
NOTEBOOK="$PROJECT_ROOT/notebooks/geneformer_pbmc3k_analysis.ipynb"
OUTPUT_DIR="/tmp/geneformer_notebook_audit"

mkdir -p "$OUTPUT_DIR"

echo "Checking project..."
cd "$PROJECT_ROOT"

if [[ ! -f "$NOTEBOOK" ]]; then
    echo "ERROR: Notebook not found: $NOTEBOOK"
    exit 1
fi

echo "Checking Jupyter and notebook execution support..."
jupyter nbconvert --version
jupyter kernelspec list

echo "Executing notebook..."
jupyter nbconvert \
    --to notebook \
    --execute "$NOTEBOOK" \
    --output geneformer_pbmc3k_analysis.executed.ipynb \
    --output-dir "$OUTPUT_DIR" \
    --ExecutePreprocessor.timeout=600 \
    --ExecutePreprocessor.kernel_name=python3

echo
echo "NOTEBOOK EXECUTION COMPLETED SUCCESSFULLY"
echo "Executed notebook: $OUTPUT_DIR/geneformer_pbmc3k_analysis.executed.ipynb"
