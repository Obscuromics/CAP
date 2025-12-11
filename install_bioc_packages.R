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
      stop("! Critical Error: BCT installation from CRAN failed. Local compilation also failed, the pipeline cannot function without BCT.")
    },
    warning = function(w) {
      message("! Warning during BCT installation.")
    }
  )
} else {
  message("✓ BCT is already installed.")
}

cat("✓ Package setup complete.\n")
