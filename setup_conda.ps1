# setup_conda.ps1
# Description: Creates the conda environment and installs the BCT package.

$ErrorActionPreference = "Stop"

$ENV_NAME = "cap-pipeline"
$ENV_FILE = "environment.yml"

# Check for submodules (TRASH2)
if (-not (Test-Path "modules/TRASH_2/README.md")) {
    Write-Host "Initializing submodules..."
    git submodule update --init --recursive
}

Write-Host "Creating conda environment '$ENV_NAME' from $ENV_FILE..."
conda env create -f $ENV_FILE

Write-Host "Setting up BCT (Bayesian Context Trees)..."
Write-Host "Attempting to compile local C++ binary (faster)..."

$ErrorActionPreference = "Continue"
conda run -n $ENV_NAME make -C bin/src/BCT
$MakeExitCode = $LASTEXITCODE
$ErrorActionPreference = "Stop"

if ($MakeExitCode -eq 0) {
    Write-Host "✓ C++ binary compiled successfully."
} else {
    Write-Host "⚠️  C++ compilation failed. Falling back to CRAN installation (slower)..."
    conda run -n $ENV_NAME Rscript install_bioc_packages.R
}

Write-Host "Setup complete! Activate the environment with:"
Write-Host "conda activate $ENV_NAME"
