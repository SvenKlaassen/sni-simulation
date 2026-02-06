# Configuration for density CI simulation
# Pointwise confidence intervals for kernel density estimation on [0,1]

config <- list(
  # Reproducibility
  seed    = 2024L,

  # Monte Carlo
  n_reps  = 200L,

  # Sample sizes
  n_values = c(200L, 1000L),

  # Significance level
  alpha   = 0.1,

  # Evaluation grid on [0,1], pulled inward to reduce boundary bias

  grid = seq(0.05, 0.95, length.out = 50),

  # Kernels (names accepted by stats::density)
  kernels = c("epanechnikov", "triangular"),

  # Bandwidth selectors: "nrd0", "SJ", "ucv", "fixed" ("fixed" uses fixed_h_scale * n^{-fixed_h_exponent})
  bw_methods = c("fixed"),

  # CI method options: "asymptotic", "wc", "bc", "sni", "clt"
  ci_methods = c("sni", "clt"),

  # Parameters for worst-case/bias-corrected CIs
  lipschitz_L  = 125,
  bc_h2_ratio  = 0.5,
  bc_u_grid_n  = 201,

  # Parameters for SNI/CLT CIs
  sni_t_multiplier = 3,
  sni_u_grid_n     = 41,
  sni_log_n_coeff  = 2,
  sni_h_scale      = NA_real_,
  sni_h_exponent   = NA_real_,
  clt_h_scale      = NA_real_,
  clt_h_exponent   = NA_real_,

  # Parameters for fixed bandwidth selection
  fixed_h_scale    = 1,
  fixed_h_exponent = 1 / 3,

  # DGP names (keys into dgp registry)
  dgp_names = c("lipschitz", "C1", "C2", "Cinfty")
)
