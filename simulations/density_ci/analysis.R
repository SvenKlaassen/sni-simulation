# Analysis script for density CI simulation
# Run after: Rscript run_simulation.R simulations/density_ci
# Usage:     Rscript simulations/density_ci/analysis.R

source("R/utils.R")
library(data.table)
library(ggplot2)

sim_dir     <- "simulations/density_ci"
results_dir <- file.path(sim_dir, "results")
results     <- load_results(results_dir)
fig_dir     <- file.path(results_dir, "figures")
safe_dir(fig_dir)

# Helpers ----------------------------------------------------------------------

save_plot <- function(p, name, width = 12, height = 8) {
  ggsave(file.path(fig_dir, paste0(name, ".png")), p,
         width = width, height = height, dpi = 150)
  message("Saved: ", name, ".png")
}

nominal <- 0.95  # 1 - alpha

# Plot 1: Pointwise coverage vs x ---------------------------------------------
# Faceted by dgp_class x n_obs, lines colored by kernel
# Fixed bw_method = "SJ"

plot_pointwise_coverage <- function(dt, bw_sel = "SJ") {
  sub <- dt[bw_method == bw_sel]
  p <- ggplot(sub, aes(x = x, y = coverage, color = kernel)) +
    geom_line(linewidth = 0.5) +
    geom_hline(yintercept = nominal, linetype = "dashed", color = "grey40") +
    facet_grid(dgp_class ~ n_obs, labeller = label_both) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(title = paste("Pointwise coverage (bw:", bw_sel, ")"),
         x = "Evaluation point x",
         y = "Coverage probability",
         color = "Kernel") +
    theme_minimal(base_size = 10) +
    theme(legend.position = "bottom")
  p
}

p1 <- plot_pointwise_coverage(results, "SJ")
save_plot(p1, "pointwise_coverage_SJ", width = 16, height = 10)

# Plot 2: Average coverage vs sample size --------------------------------------
# Aggregate mean coverage over grid, facet by dgp_class

plot_avg_coverage_vs_n <- function(dt) {
  agg <- dt[, .(avg_coverage = mean(coverage)),
            by = .(dgp_class, n_obs, kernel, bw_method)]
  p <- ggplot(agg, aes(x = n_obs, y = avg_coverage,
                        color = kernel, linetype = bw_method)) +
    geom_line() +
    geom_point(size = 1.5) +
    geom_hline(yintercept = nominal, linetype = "dashed", color = "grey40") +
    facet_wrap(~ dgp_class, ncol = 2) +
    scale_x_log10() +
    labs(title = "Average pointwise coverage vs sample size",
         x = "Sample size n (log scale)",
         y = "Average pointwise coverage",
         color = "Kernel",
         linetype = "Bandwidth") +
    theme_minimal(base_size = 10) +
    theme(legend.position = "bottom")
  p
}

p2 <- plot_avg_coverage_vs_n(results)
save_plot(p2, "avg_coverage_vs_n")

# Plot 3: CI width vs sample size ----------------------------------------------

plot_width_vs_n <- function(dt) {
  agg <- dt[, .(avg_width = mean(mean_width)),
            by = .(dgp_class, n_obs, kernel, bw_method)]
  p <- ggplot(agg, aes(x = n_obs, y = avg_width,
                        color = kernel, linetype = bw_method)) +
    geom_line() +
    geom_point(size = 1.5) +
    facet_wrap(~ dgp_class, ncol = 2) +
    scale_x_log10() +
    scale_y_log10() +
    labs(title = "Average CI width vs sample size",
         x = "Sample size n (log scale)",
         y = "Average CI width (log scale)",
         color = "Kernel",
         linetype = "Bandwidth") +
    theme_minimal(base_size = 10) +
    theme(legend.position = "bottom")
  p
}

p3 <- plot_width_vs_n(results)
save_plot(p3, "ci_width_vs_n")

# Plot 4: Coverage heatmap -----------------------------------------------------
# x vs n_obs, fill = coverage, per DGP

plot_coverage_heatmap <- function(dt, kern = "gaussian", bw_sel = "SJ") {
  sub <- dt[kernel == kern & bw_method == bw_sel]
  p <- ggplot(sub, aes(x = x, y = factor(n_obs), fill = coverage)) +
    geom_tile() +
    scale_fill_gradient2(low = "red", mid = "white", high = "blue",
                         midpoint = nominal, limits = c(0.5, 1)) +
    facet_wrap(~ dgp_class, ncol = 2) +
    labs(title = paste("Coverage heatmap:", kern, "kernel, bw:", bw_sel),
         x = "Evaluation point x",
         y = "Sample size n",
         fill = "Coverage") +
    theme_minimal(base_size = 10)
  p
}

p4 <- plot_coverage_heatmap(results, "gaussian", "SJ")
save_plot(p4, "coverage_heatmap_gaussian_SJ")

# Plot 5: Bias profile ---------------------------------------------------------
# Bias vs x, per DGP, colored by kernel

plot_bias_profile <- function(dt, bw_sel = "SJ") {
  sub <- dt[bw_method == bw_sel]
  p <- ggplot(sub, aes(x = x, y = bias, color = kernel)) +
    geom_line(linewidth = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    facet_grid(dgp_class ~ n_obs, labeller = label_both, scales = "free_y") +
    labs(title = paste("Bias profile (bw:", bw_sel, ")"),
         x = "Evaluation point x",
         y = "Bias (mean f_hat - f_true)",
         color = "Kernel") +
    theme_minimal(base_size = 10) +
    theme(legend.position = "bottom")
  p
}

p5 <- plot_bias_profile(results, "SJ")
save_plot(p5, "bias_profile_SJ", width = 16, height = 10)

message("\nAll figures saved to: ", fig_dir)
