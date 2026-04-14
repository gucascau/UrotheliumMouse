################################################################################
# Script 00: Install missing R packages
# Run this ONCE interactively before submitting the pipeline:
#   srun --partition=himem --mem=32G --cpus-per-task=4 --pty bash
#   module purge && module load GCC/9.3.0 OpenMPI/4.0.3 R/4.4.0
#   Rscript 00_install_packages.R
################################################################################

# BioConductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("DropletUtils", quietly = TRUE)) {
  message("Installing DropletUtils...")
  BiocManager::install("DropletUtils", ask = FALSE)
}

# GitHub-only packages (not on CRAN or Bioconductor)
if (!requireNamespace("remotes", quietly = TRUE))
  install.packages("remotes")

if (!requireNamespace("DoubletFinder", quietly = TRUE)) {
  message("Installing DoubletFinder from GitHub...")
  remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")
}

if (!requireNamespace("BPCells", quietly = TRUE)) {
  message("Installing BPCells from GitHub...")
  remotes::install_github("bnprks/BPCells/r")
}

# Confirm all packages load
pkgs <- c("Seurat", "harmony", "BPCells", "Matrix", "data.table", "dplyr",
          "ggplot2", "patchwork", "DropletUtils", "SingleCellExperiment",
          "zellkonverter", "DoubletFinder", "future")

results <- sapply(pkgs, requireNamespace, quietly = TRUE)
cat("\n===== Package check =====\n")
for (p in names(results)) {
  cat(sprintf("  %-25s %s\n", p, ifelse(results[p], "OK", "MISSING")))
}

if (all(results)) {
  cat("\nAll packages ready. You can now submit the pipeline.\n")
} else {
  cat("\nWARNING: Some packages still missing — do not submit yet.\n")
}
