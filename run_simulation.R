#!/usr/bin/env Rscript
# Top-level simulation runner
# Usage: Rscript run_simulation.R <simulation_directory>
# Example: Rscript run_simulation.R simulations/example_sim

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript run_simulation.R <simulation_directory>")
}

sim_dir <- args[1]
run_script <- file.path(sim_dir, "run.R")

if (!file.exists(run_script)) {
  stop("run.R not found in: ", sim_dir)
}

message("=== Running simulation: ", sim_dir, " ===\n")
source(run_script)
message("\n=== Done ===")
