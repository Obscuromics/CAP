#!/usr/bin/env Rscript
# Install Bioconductor packages that aren't available via conda

# Install BCT from CRAN
if (!require("BCT", quietly = TRUE)) {
  install.packages("BCT", repos='https://cloud.r-project.org/')
}

cat("✓ All additional packages installed successfully!\n")
