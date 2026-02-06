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

ci_methods <- if (!is.null(config$ci_methods)) config$ci_methods else "asymptotic"

# Build parameter grid (all combos of dgp x n x kernel x bw)
param_grid <- CJ(
  dgp_name  = config$dgp_names,
  n_obs     = config$n_values,
  kernel    = config$kernels,
  bw_method = config$bw_methods,
  ci_method = ci_methods,
  sorted    = FALSE
)

if (any(param_grid$ci_method %in% c("wc", "bc", "sni"))) {
  bad <- param_grid[kernel == "gaussian" & ci_method %in% c("wc", "bc", "sni")]
  if (nrow(bad) > 0) {
    message("Dropping gaussian kernel for wc/bc/sni CI methods (compact support required).")
    param_grid <- param_grid[!(kernel == "gaussian" & ci_method %in% c("wc", "bc", "sni"))]
  }
}
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
  ci_m     <- param_grid$ci_method[cc]

  dgp    <- dgp_registry[[dgp_name]]
  f_true <- true_densities[[dgp_name]]

  message(sprintf("[%d/%d] dgp=%-10s n=%-5d kernel=%-13s bw=%s ci=%s",
                  cc, n_combos, dgp_name, n_obs, kern, bw_m, ci_m))

  # Accumulators
  cover_count <- numeric(n_grid)
  width_sum   <- numeric(n_grid)
  fhat_sum    <- numeric(n_grid)
  bw_sum      <- 0
  bw2_sum     <- 0

  for (r in seq_len(config$n_reps)) {
    y   <- dgp$sample(n_obs)
    h1 <- compute_bandwidth_ci(
      y, bw_m, ci_m,
      config$sni_h_scale, config$sni_h_exponent,
      config$clt_h_scale, config$clt_h_exponent,
      config$fixed_h_scale, config$fixed_h_exponent
    )
    u_mat <- outer(y, grid, function(yy, xx) (yy - xx) / h1)
    k_mat <- kernel_eval(u_mat, kern)
    est <- estimate_density_ci(
      y, grid, kern, bw_m, config$alpha,
      ci_method   = ci_m,
      lipschitz_L = config$lipschitz_L,
      bc_h2_ratio = config$bc_h2_ratio,
      bc_u_grid_n = config$bc_u_grid_n,
      sni_t_multiplier = config$sni_t_multiplier,
      sni_u_grid_n     = config$sni_u_grid_n,
      sni_log_n_coeff  = config$sni_log_n_coeff,
      sni_h_scale      = config$sni_h_scale,
      sni_h_exponent   = config$sni_h_exponent,
      clt_h_scale      = config$clt_h_scale,
      clt_h_exponent   = config$clt_h_exponent,
      fixed_h_scale    = config$fixed_h_scale,
      fixed_h_exponent = config$fixed_h_exponent,
      h1               = h1,
      k_mat            = k_mat
    )

    covers      <- (est$ci_lower <= f_true) & (f_true <= est$ci_upper)
    cover_count <- cover_count + covers
    width_sum   <- width_sum + (est$ci_upper - est$ci_lower)
    fhat_sum    <- fhat_sum + est$f_hat
    bw_sum      <- bw_sum + est$bandwidth[1]
    if (ci_m == "bc") {
      bw2_sum <- bw2_sum + est$bandwidth2[1]
    }
  }

  # Summarize over replications
  results_list[[cc]] <- data.table(
    dgp_name   = dgp_name,
    dgp_class  = dgp$class,
    n_obs      = n_obs,
    kernel     = kern,
    bw_method  = bw_m,
    ci_method  = ci_m,
    x          = grid,
    f_true     = f_true,
    coverage   = cover_count / config$n_reps,
    mean_width = width_sum / config$n_reps,
    mean_fhat  = fhat_sum / config$n_reps,
    bias       = (fhat_sum / config$n_reps) - f_true,
    mean_bw    = bw_sum / config$n_reps,
    mean_bw2   = if (ci_m == "bc") bw2_sum / config$n_reps else NA_real_
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
), by = .(dgp_class, n_obs, ci_method)]

message("\n--- Summary by DGP class and sample size ---")
print(summary_dt)
