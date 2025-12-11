#!/bin/bash
# setup_conda.sh
# Description: Creates the conda environment and installs the BCT package.

set -e

ENV_NAME="cap-pipeline"
ENV_FILE="environment.yml"

echo "Creating conda environment '$ENV_NAME' from $ENV_FILE..."
conda env create -f $ENV_FILE

echo "Activating environment..."
# Need to source conda.sh to use 'conda activate' in script, or use 'conda run'
# Assuming 'conda run' is available (newer conda versions)

echo "Installing BCT package..."
conda run -n $ENV_NAME Rscript install_bioc_packages.R

echo "Setting permissions..."
chmod +x modules/TRASH_2/src/TRASH.R

echo "Setup complete! Activate the environment with:"
echo "conda activate $ENV_NAME"
