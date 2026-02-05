# Shared utility functions for simulations
# Source this from simulation run scripts: source("R/utils.R")

#' Ensure a directory exists, creating it if necessary
safe_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
  invisible(path)
}

#' Save simulation results to both RDS and CSV
#'
#' @param results A data.frame of simulation results
#' @param dir Output directory
#' @param name Base filename (without extension)
save_results <- function(results, dir, name = "results") {
  safe_dir(dir)
  saveRDS(results, file.path(dir, paste0(name, ".rds")))
  write.csv(results, file.path(dir, paste0(name, ".csv")), row.names = FALSE)
  message("Results saved to ", dir)
  invisible(results)
}

#' Load simulation results from RDS
#'
#' @param dir Output directory
#' @param name Base filename (without extension)
load_results <- function(dir, name = "results") {
  path <- file.path(dir, paste0(name, ".rds"))
  if (!file.exists(path)) {
    stop("Results file not found: ", path)
  }
  readRDS(path)
}
