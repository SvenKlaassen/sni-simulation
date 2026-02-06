# Run script for density CI simulation
# Can be called directly or via: Rscript run_simulation.R simulations/density_ci

# Resolve paths relative to project root
sim_dir <- "simulations/density_ci"
source("R/utils.R")
source(file.path(sim_dir, "config.R"))
source(file.path(sim_dir, "dgp.R"))
source(file.path(sim_dir, "estimator.R"))

library(data.table)

# Setup ------------------------------------------------------------------------
set.seed(config$seed)
dgp_registry <- get_dgp_registry()
grid   <- config$grid
n_grid <- length(grid)

# Pre-compute true densities at grid points for each DGP
true_densities <- lapply(dgp_registry, function(dgp) dgp$density(grid))

# Build parameter grid (all combos of dgp x n x kernel x bw)
param_grid <- CJ(
  dgp_name  = config$dgp_names,
  n_obs     = config$n_values,
  kernel    = config$kernels,
  bw_method = config$bw_methods,
  sorted    = FALSE
)
n_combos <- nrow(param_grid)

message("Density CI simulation: ", n_combos, " parameter combos x ",
        config$n_reps, " reps = ", format(n_combos * config$n_reps, big.mark = ","),
        " total runs")

# Main loop --------------------------------------------------------------------
# Accumulate-then-summarize: keep running sums per grid point, not per-rep data
results_list <- vector("list", n_combos)
t_start <- proc.time()

for (cc in seq_len(n_combos)) {
  dgp_name <- param_grid$dgp_name[cc]
  n_obs    <- param_grid$n_obs[cc]
  kern     <- param_grid$kernel[cc]
  bw_m     <- param_grid$bw_method[cc]

  dgp    <- dgp_registry[[dgp_name]]
  f_true <- true_densities[[dgp_name]]

  message(sprintf("[%d/%d] dgp=%-10s n=%-5d kernel=%-13s bw=%s",
                  cc, n_combos, dgp_name, n_obs, kern, bw_m))

  # Accumulators
  cover_count <- numeric(n_grid)
  width_sum   <- numeric(n_grid)
  fhat_sum    <- numeric(n_grid)
  bw_sum      <- 0

  for (r in seq_len(config$n_reps)) {
    y   <- dgp$sample(n_obs)
    est <- estimate_density_ci(y, grid, kern, bw_m, config$alpha)

    covers      <- (est$ci_lower <= f_true) & (f_true <= est$ci_upper)
    cover_count <- cover_count + covers
    width_sum   <- width_sum + (est$ci_upper - est$ci_lower)
    fhat_sum    <- fhat_sum + est$f_hat
    bw_sum      <- bw_sum + est$bandwidth[1]
  }

  # Summarize over replications
  results_list[[cc]] <- data.table(
    dgp_name   = dgp_name,
    dgp_class  = dgp$class,
    n_obs      = n_obs,
    kernel     = kern,
    bw_method  = bw_m,
    x          = grid,
    f_true     = f_true,
    coverage   = cover_count / config$n_reps,
    mean_width = width_sum / config$n_reps,
    mean_fhat  = fhat_sum / config$n_reps,
    bias       = (fhat_sum / config$n_reps) - f_true,
    mean_bw    = bw_sum / config$n_reps
  )
}

elapsed <- (proc.time() - t_start)["elapsed"]
message(sprintf("\nSimulation completed in %.1f seconds", elapsed))

results <- rbindlist(results_list)

# Save results -----------------------------------------------------------------
results_dir <- file.path(sim_dir, "results")
save_results(results, results_dir)

# Print summary ----------------------------------------------------------------
summary_dt <- results[, .(
  mean_coverage = mean(coverage),
  mean_width    = mean(mean_width),
  mean_abs_bias = mean(abs(bias))
), by = .(dgp_class, n_obs)]

message("\n--- Summary by DGP class and sample size ---")
print(summary_dt)
