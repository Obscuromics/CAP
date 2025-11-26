#!/usr/bin/env Rscript
# Install Bioconductor packages that aren't available via conda

if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos='https://cloud.r-project.org/')
}

# Install MSA from Bioconductor
BiocManager::install("msa")

# Install BCT from CRAN
install.packages("BCT", repos='https://cloud.r-project.org/')

cat("✓ All additional packages installed successfully!\n")
