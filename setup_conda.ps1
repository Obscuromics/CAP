# setup_conda.ps1
# Description: Creates the conda environment and installs the BCT package.

$ErrorActionPreference = "Stop"

$ENV_NAME = "cap-pipeline"
$ENV_FILE = "environment.yml"

Write-Host "Creating conda environment '$ENV_NAME' from $ENV_FILE..."
conda env create -f $ENV_FILE

Write-Host "Installing BCT package..."
conda run -n $ENV_NAME Rscript install_bioc_packages.R

Write-Host "Setup complete! Activate the environment with:"
Write-Host "conda activate $ENV_NAME"
