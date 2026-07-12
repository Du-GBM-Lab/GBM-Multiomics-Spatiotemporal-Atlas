cran_packages <- c(
  "remotes",
  "BiocManager",
  "processx",
  "callr",
  "pkgbuild",
  "qs2",
  "rjags",
  "tidyr",
  "dplyr",
  "ggplot2",
  "patchwork"
)

missing_cran <- setdiff(cran_packages, rownames(installed.packages()))
if (length(missing_cran) > 0) {
  install.packages(missing_cran, repos = "https://cloud.r-project.org")
}

if (!requireNamespace("DoubletFinder", quietly = TRUE)) {
  remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")
}

if (!requireNamespace("infercnv", quietly = TRUE)) {
  if (!requireNamespace("rjags", quietly = TRUE)) {
    stop("rjags is required before installing infercnv. Install system-level JAGS 4 first.")
  }
  if (!grepl("JAGS", Sys.getenv("PATH"), ignore.case = TRUE)) {
    message("JAGS is not visible in PATH. Install JAGS 4.x and add it to PATH (e.g. C:/Program Files/JAGS/ on Windows) before installing infercnv.")
  }
  py <- Sys.getenv("PYTHON")
  if (!nzchar(py)) {
    message("Set the PYTHON environment variable to a Python 3 executable (used by infercnv/reticulate) if auto-detection fails.")
  }
  BiocManager::install("infercnv", ask = FALSE, update = FALSE)
}

cat("Package installation check completed.\n")
