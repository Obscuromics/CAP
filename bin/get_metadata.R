#!/usr/bin/env Rscript

# get_metadata.R
# Description: Extract metadata from assembly or read from provided metadata file

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: get_metadata.R <assembly.fasta> <output_metadata.csv> [-m <metadata_file>]")
}

assembly_fasta <- args[1]
output_metadata <- args[2]

# Check for optional metadata file argument
metadata_file <- NULL
if (length(args) >= 4 && args[3] == "-m") {
  metadata_file <- args[4]
}

# Load libraries
suppressMessages({library(Biostrings)})

# Load additional functions
source(file.path(Sys.getenv("WORKFLOW_DIR"), "bin", "auxfuns.R"))

# Load data or use provided metadata
if (!is.null(metadata_file)) {
  # Read metadata from provided file
  metadata_data <- read.csv(metadata_file)

  assembly.name <- basename(assembly_fasta)

  metadata_data <- metadata_data[metadata_data$assembly.name == assembly.name, ]

} else {
  # Generate metadata from assembly
  fasta <- readDNAStringSet(assembly_fasta)
  
  assembly.name <- basename(assembly_fasta)
  
  metadata_data <- data.frame(assembly.name = rep(assembly.name, length(fasta)),
                              chromosome.name = names(fasta),
                              size = width(fasta),
                              is.chr = rep(1, length(fasta)),
                              alternative_name = names(fasta))
  
  for(i in seq_len(nrow(metadata_data))) {
    metadata_data$chromosome.name[i] <- strsplit(metadata_data$chromosome.name[i], " ")[[1]][1]
  }
}

# Save output
write.csv(metadata_data, file = output_metadata, row.names = FALSE) 