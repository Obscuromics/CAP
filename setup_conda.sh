#!/bin/bash
# setup_conda.sh
# Description: Creates the conda environment and installs the BCT package.

set -e

ENV_NAME="cap-pipeline"
ENV_FILE="environment.yml"

# Check for submodules (TRASH2)
if [ -z "$(ls -A modules/TRASH_2 2>/dev/null)" ]; then
    echo "Initializing submodules..."
    git submodule update --init --recursive
fi

echo "Creating conda environment '$ENV_NAME' from $ENV_FILE..."
conda env create -f $ENV_FILE

echo "Activating environment..."
# Need to source conda.sh to use 'conda activate' in script, or use 'conda run'
# Assuming 'conda run' is available (newer conda versions)

echo "Setting up BCT (Bayesian Context Trees)..."
echo "Attempting to compile local C++ binary (faster)..."

# Try compilation first
if conda run -n $ENV_NAME make -C bin/src/BCT; then
    echo "✓ C++ binary compiled successfully."
else
    echo "⚠️  C++ compilation failed. Falling back to CRAN installation (slower)..."
    conda run -n $ENV_NAME Rscript install_bioc_packages.R
fi

echo "Setting permissions..."
chmod +x modules/TRASH_2/src/TRASH.R

echo "Setup complete! Activate the environment with:"
echo "conda activate $ENV_NAME"
