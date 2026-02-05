# Run script for example simulation
# Can be called directly or via: Rscript run_simulation.R simulations/example_sim

# Resolve paths relative to project root
sim_dir <- "simulations/example_sim"
source("R/utils.R")
source(file.path(sim_dir, "config.R"))
source(file.path(sim_dir, "dgp.R"))
source(file.path(sim_dir, "estimator.R"))

# Run simulation ---------------------------------------------------------------
set.seed(config$seed)
message("Running ", config$n_reps, " replications (n=", config$n_obs, ")...")

results <- data.frame(
  rep       = seq_len(config$n_reps),
  estimate  = numeric(config$n_reps),
  ci_lower  = numeric(config$n_reps),
  ci_upper  = numeric(config$n_reps),
  covers    = logical(config$n_reps)
)

for (i in seq_len(config$n_reps)) {
  y   <- generate_data(config$n_obs, config$mu_true, config$sigma)
  est <- estimate(y, alpha = config$alpha)

  results$estimate[i] <- est$estimate
  results$ci_lower[i] <- est$ci_lower
  results$ci_upper[i] <- est$ci_upper
  results$covers[i]   <- (est$ci_lower <= config$mu_true) & (config$mu_true <= est$ci_upper)
}

# Save results -----------------------------------------------------------------
results_dir <- file.path(sim_dir, "results")
save_results(results, results_dir)

# Print summary ----------------------------------------------------------------
coverage <- mean(results$covers)
bias     <- mean(results$estimate) - config$mu_true

message("\n--- Summary ---")
message("Coverage: ", sprintf("%.1f%%", coverage * 100),
        " (nominal: ", sprintf("%.1f%%", (1 - config$alpha) * 100), ")")
message("Bias:     ", sprintf("%.4f", bias))
message("Mean SE:  ", sprintf("%.4f", mean(results$ci_upper - results$ci_lower) / (2 * qt(1 - config$alpha / 2, config$n_obs - 1))))
