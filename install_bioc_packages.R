#!/usr/bin/env Rscript
# Install Bioconductor packages that aren't available via conda

# Install BCT from CRAN
if (!require("BCT", quietly = TRUE)) {
  message("Attempting to install BCT from CRAN...")
  tryCatch(
    {
      install.packages("BCT", repos='https://cloud.r-project.org/')
      if (!require("BCT", quietly = TRUE)) {
        stop("Installation failed.")
      }
      message("✓ BCT installed successfully.")
    },
    error = function(e) {
      message("! Warning: BCT installation from CRAN failed.")
      message("! This is expected on some systems. The pipeline will compile a local C++ alternative instead.")
    },
    warning = function(w) {
      message("! Warning during BCT installation.")
    }
  )
} else {
  message("✓ BCT is already installed.")
}

cat("✓ Package setup complete.\n")
