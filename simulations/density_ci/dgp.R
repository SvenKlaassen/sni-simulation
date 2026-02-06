# Data-generating processes for density CI simulation
# Four density classes on [0,1] with increasing smoothness

#' Registry of data-generating processes
#'
#' Each DGP is a list with:
#'   $name    - descriptive label
#'   $class   - smoothness class label
#'   $sample  - function(n) returning n draws on [0,1]
#'   $density - function(x) returning true density at x
#'
#' @return Named list of DGP objects
get_dgp_registry <- function() {

  # --- Lipschitz continuous (not C1): mixture of narrow triangles ---
  # Multiple sharp peaks increase local Lipschitz constants.
  tri_centers <- c(0.15, 0.35, 0.55, 0.75, 0.90)
  tri_widths  <- c(0.05, 0.04, 0.06, 0.04, 0.03)  # half-widths
  tri_weights <- c(0.25, 0.20, 0.25, 0.20, 0.10)

  tri_density <- function(x, center, half_width) {
    u <- (x - center) / half_width
    ifelse(abs(u) <= 1, (1 - abs(u)) / half_width, 0)
  }

  lipschitz_density <- function(x) {
    dens <- 0
    for (ii in seq_along(tri_centers)) {
      dens <- dens + tri_weights[ii] *
        tri_density(x, tri_centers[ii], tri_widths[ii])
    }
    dens
  }

  lipschitz_sample <- function(n) {
    comp <- sample(seq_along(tri_centers), size = n, replace = TRUE,
                   prob = tri_weights)
    u1 <- runif(n)
    u2 <- runif(n)
    tri_centers[comp] + tri_widths[comp] * (u1 - u2)
  }

  # --- C1: Beta(2,2) ---
  # f(x) = 6x(1-x), symmetric, smooth, well-understood baseline
  c1_density <- function(x) dbeta(x, 2, 2)
  c1_sample  <- function(n) rbeta(n, 2, 2)

  # --- C2: Mixture of Betas ---
  # 0.5 * Beta(3,3) + 0.5 * Beta(6,3)
  # Asymmetric density with bump, tests higher-curvature regions
  c2_density <- function(x) {
    0.5 * dbeta(x, 3, 3) + 0.5 * dbeta(x, 6, 3)
  }
  c2_sample <- function(n) {
    # Draw component indicator, then sample from each component
    comp <- rbinom(n, size = 1, prob = 0.5)
    ifelse(comp == 1, rbeta(n, 3, 3), rbeta(n, 6, 3))
  }

  # --- C-infinity: truncated Normal(0.5, 0.15) on [0,1] ---
  # Very smooth, near-Gaussian, best case for KDE
  mu_tn <- 0.5
  sd_tn <- 0.15
  F_lo  <- pnorm(0, mu_tn, sd_tn)
  F_hi  <- pnorm(1, mu_tn, sd_tn)
  C_tn  <- F_hi - F_lo

  cinfty_density <- function(x) dnorm(x, mu_tn, sd_tn) / C_tn
  cinfty_sample  <- function(n) {
    # Inverse CDF of truncated normal
    u <- runif(n)
    qnorm(F_lo + u * C_tn, mu_tn, sd_tn)
  }

  list(
    lipschitz = list(
      name = "Lipschitz (multi-triangle)", class = "Lipschitz",
      sample = lipschitz_sample, density = lipschitz_density
    ),
    C1 = list(
      name = "Beta(2,2)", class = "C1",
      sample = c1_sample, density = c1_density
    ),
    C2 = list(
      name = "Beta mixture", class = "C2",
      sample = c2_sample, density = c2_density
    ),
    Cinfty = list(
      name = "Trunc. Normal", class = "Cinfty",
      sample = cinfty_sample, density = cinfty_density
    )
  )
}
