################################################################################
# fix_lapack_packages.R
#
# Removes harmony2 (which requires LAPACK via arma::inv()) and installs
# harmony v1.2.4 (last stable v1 release) from the CRAN archive.
################################################################################

# ── 1. Remove current harmony (harmony2) ────────────────────────────────────
message("=== Removing current harmony package ===")
if (requireNamespace("harmony", quietly = TRUE)) {
  remove.packages("harmony")
  message("Removed.")
} else {
  message("harmony not found — skipping removal.")
}

# ── 2. Install harmony v1.2.4 from CRAN archive ──────────────────────────────
message("\n=== Installing harmony v1.2.4 from CRAN archive ===")
if (!requireNamespace("remotes", quietly = TRUE))
  install.packages("remotes", repos = "https://cloud.r-project.org")

remotes::install_version(
  "harmony",
  version = "1.2.4",
  repos   = "https://cloud.r-project.org"
)

# ── 3. Verify ────────────────────────────────────────────────────────────────
ver <- as.character(packageVersion("harmony"))
message("\n=== Done. Installed harmony version: ", ver, " ===")
if (startsWith(ver, "2")) {
  message("WARNING: still seeing harmony2 — check your R library path.")
} else {
  message("You can now resubmit submit_01_integrate.sh")
}
