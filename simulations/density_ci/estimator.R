# Kernel density estimation with pointwise confidence intervals

library(data.table)

#' Integrated squared kernel constant R(K) = int K^2(u) du
#'
#' @param kernel Character name of kernel (as accepted by stats::density)
#' @return Numeric value of R(K)
kernel_Rk <- function(kernel) {
  switch(kernel,
    gaussian     = 1 / (2 * sqrt(pi)),
    epanechnikov = 3 / 5,
    rectangular  = 1 / 2,
    triangular   = 2 / 3,
    biweight     = 5 / 7,
    stop("Unknown kernel: ", kernel)
  )
}

#' Estimate density and construct asymptotic pointwise CIs
#'
#' Uses stats::density() for KDE, then constructs normal CIs based on
#' the asymptotic variance: Var(f_hat(x)) ~ f(x) * R(K) / (n * h).
#' The true f(x) is replaced by the plug-in estimate f_hat(x).
#'
#' @param y         Numeric vector of observations on [0,1]
#' @param grid      Numeric vector of evaluation points
#' @param kernel    Character kernel name for stats::density()
#' @param bw_method Character bandwidth selector name (e.g. "nrd0", "SJ", "ucv")
#' @param alpha     Significance level (default 0.05)
#' @return data.table with columns: x, f_hat, se_hat, ci_lower, ci_upper, bandwidth
estimate_density_ci <- function(y, grid, kernel, bw_method, alpha = 0.05) {
  n <- length(y)

  # Compute bandwidth; bw.SJ and bw.ucv can fail on small/degenerate samples
  h <- tryCatch(
    do.call(paste0("bw.", bw_method), list(y)),
    error = function(e) bw.nrd0(y)
  )

  # Extend evaluation range beyond grid to avoid edge truncation
  dens <- density(y, bw = h, kernel = kernel,
                  from = min(grid) - 3 * h,
                  to   = max(grid) + 3 * h,
                  n    = 1024)

  # Interpolate f_hat onto the evaluation grid
  f_hat <- approx(dens$x, dens$y, xout = grid)$y

  # Asymptotic pointwise standard error
  Rk     <- kernel_Rk(kernel)
  se_hat <- sqrt(pmax(f_hat, 0) * Rk / (n * h))

  # Normal CI, clamped at 0 (density is non-negative)
  z        <- qnorm(1 - alpha / 2)
  ci_lower <- pmax(f_hat - z * se_hat, 0)
  ci_upper <- f_hat + z * se_hat

  data.table(
    x         = grid,
    f_hat     = f_hat,
    se_hat    = se_hat,
    ci_lower  = ci_lower,
    ci_upper  = ci_upper,
    bandwidth = h
  )
}
