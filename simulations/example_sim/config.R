# Configuration for example simulation: coverage of normal mean CI
# Modify these parameters to change the simulation setup

config <- list(
  # Simulation parameters
  n_reps    = 1000L,   # Number of Monte Carlo replications
  seed      = 42L,     # Base random seed for reproducibility

  # DGP parameters
  n_obs     = 100L,    # Sample size per replication
  mu_true   = 5,       # True mean
  sigma     = 2,       # True standard deviation

  # Estimator parameters
  alpha     = 0.05     # Significance level for confidence interval
)
