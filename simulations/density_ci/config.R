# Configuration for density CI simulation
# Pointwise confidence intervals for kernel density estimation on [0,1]

config <- list(
  # Reproducibility
  seed    = 2024L,

  # Monte Carlo
  n_reps  = 1000L,

  # Sample sizes
  n_values = c(50L, 100L, 200L, 500L, 1000L, 2000L),

  # Significance level
  alpha   = 0.05,

  # Evaluation grid on [0,1], pulled inward to reduce boundary bias

  grid    = seq(0.05, 0.95, length.out = 100),

  # Kernels (names accepted by stats::density)
  kernels = c("gaussian", "epanechnikov", "rectangular",
              "triangular", "biweight"),

  # Bandwidth selectors (passed to bw argument of density())
  bw_methods = c("nrd0", "SJ", "ucv"),

  # DGP names (keys into dgp registry)
  dgp_names = c("lipschitz", "C1", "C2", "Cinfty")
)
