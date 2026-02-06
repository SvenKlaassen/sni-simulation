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

#' Evaluate compact-support kernels K(u)
#'
#' @param u Numeric vector/matrix of scaled inputs
#' @param kernel Character kernel name
#' @return Numeric vector/matrix of kernel values
kernel_eval <- function(u, kernel) {
  switch(kernel,
    epanechnikov = ifelse(abs(u) <= 1, 0.75 * (1 - u^2), 0),
    rectangular  = ifelse(abs(u) <= 1, 0.5, 0),
    triangular   = ifelse(abs(u) <= 1, 1 - abs(u), 0),
    biweight     = ifelse(abs(u) <= 1, (15 / 16) * (1 - u^2)^2, 0),
    gaussian     = dnorm(u),
    stop("Unknown kernel: ", kernel)
  )
}

compute_bandwidth <- function(y, bw_method) {
  tryCatch(
    do.call(paste0("bw.", bw_method), list(y)),
    error = function(e) bw.nrd0(y)
  )
}

compute_bandwidth_ci <- function(y, bw_method, ci_method,
                                 sni_h_scale, sni_h_exponent,
                                 clt_h_scale, clt_h_exponent,
                                 fixed_h_scale, fixed_h_exponent) {
  n <- length(y)
  if (bw_method == "fixed") {
    if (is.na(fixed_h_scale) || is.na(fixed_h_exponent)) {
      stop("fixed bandwidth requires fixed_h_scale and fixed_h_exponent.")
    }
    h1 <- fixed_h_scale * n^(-fixed_h_exponent)
  } else {
    h1 <- compute_bandwidth(y, bw_method)
  }
  if (bw_method != "fixed" && ci_method == "sni" && !is.na(sni_h_scale) && !is.na(sni_h_exponent)) {
    h1 <- sni_h_scale * n^(-sni_h_exponent)
  }
  if (bw_method != "fixed" && ci_method == "clt" && !is.na(clt_h_scale) && !is.na(clt_h_exponent)) {
    h1 <- clt_h_scale * n^(-clt_h_exponent)
  }
  h1
}

kde_on_grid <- function(y, grid, h, kernel, pad = 3) {
  dens <- density(
    y, bw = h, kernel = kernel,
    from = min(grid) - pad * h,
    to   = max(grid) + pad * h,
    n    = 1024
  )
  approx(dens$x, dens$y, xout = grid)$y
}

kde_density <- function(y, grid, h, kernel, pad = 3, n = 1024) {
  density(
    y, bw = h, kernel = kernel,
    from = min(grid) - pad * h,
    to   = max(grid) + pad * h,
    n    = n
  )
}

compute_Rn <- function(v_hat, n, h, alpha) {
  log_u <- log(1.5 / alpha)
  denom <- n * h
  nu_n  <- if (sqrt(log_u / n) < 1) 1 / (1 - sqrt(log_u / n)) else Inf
  nu_n * (sqrt(2 * v_hat * log_u / denom) + 3.15 * log_u / denom)
}

compute_Rsni <- function(v_hat, n, h, alpha, t_multiplier = 3, log_n_coeff = 2) {
  t_val <- t_multiplier / alpha
  denom <- n * h
  sqrt(2 * v_hat * log(t_val) / denom + log_n_coeff * log(n) / (denom^2))
}

#' Estimate density and construct pointwise CIs
#'
#' Methods:
#'   - asymptotic: normal CI with plug-in variance
#'   - wc: worst-case CI for Lipschitz class
#'   - bc: bias-corrected CI using small-bandwidth bias estimate
#'   - sni: plug-in bias correction with SNI radius
#'   - clt: normal CI using sample variance of K((X-x)/h)/h
#'
#' @param y         Numeric vector of observations on [0,1]
#' @param grid      Numeric vector of evaluation points
#' @param kernel    Character kernel name for stats::density()
#' @param bw_method Character bandwidth selector name (e.g. "nrd0", "SJ", "ucv", "fixed")
#' @param alpha     Significance level (default 0.05)
#' @param ci_method Character CI method: "asymptotic", "wc", or "bc"
#' @param lipschitz_L Lipschitz constant for wc CI
#' @param bc_h2_ratio Ratio for secondary bandwidth in bc CI
#' @param bc_u_grid_n Number of integration points for bc CI
#' @param sni_t_multiplier Multiplier in t = multiplier / alpha for SNI radius
#' @param sni_u_grid_n Number of integration points for SNI bias correction
#' @param sni_log_n_coeff Coefficient for log(n) term in SNI radius
#' @param sni_h_scale Optional bandwidth scale for SNI (used with exponent)
#' @param sni_h_exponent Optional bandwidth exponent for SNI
#' @param clt_h_scale Optional bandwidth scale for CLT (used with exponent)
#' @param clt_h_exponent Optional bandwidth exponent for CLT
#' @param fixed_h_scale Fixed bandwidth scale (used when bw_method = "fixed")
#' @param fixed_h_exponent Fixed bandwidth exponent (used when bw_method = "fixed")
#' @return data.table with columns: x, f_hat, se_hat, ci_lower, ci_upper, bandwidth, bandwidth2
estimate_density_ci <- function(y, grid, kernel, bw_method, alpha = 0.05,
                                ci_method = "asymptotic",
                                lipschitz_L = NULL,
                                bc_h2_ratio = 0.5,
                                bc_u_grid_n = 201,
                                sni_t_multiplier = 3,
                                sni_u_grid_n = 41,
                                sni_log_n_coeff = 2,
                                sni_h_scale = NA_real_,
                                sni_h_exponent = NA_real_,
                                clt_h_scale = NA_real_,
                                clt_h_exponent = NA_real_,
                                fixed_h_scale = NA_real_,
                                fixed_h_exponent = NA_real_,
                                h1 = NULL,
                                k_mat = NULL) {
  n <- length(y)

  ci_method <- match.arg(ci_method, c("asymptotic", "wc", "bc", "sni", "clt"))

  # Compute bandwidth; bw.SJ and bw.ucv can fail on small/degenerate samples
  if (is.null(h1)) {
    h1 <- compute_bandwidth_ci(
      y, bw_method, ci_method,
      sni_h_scale, sni_h_exponent,
      clt_h_scale, clt_h_exponent,
      fixed_h_scale, fixed_h_exponent
    )
  }
  if (!is.null(k_mat)) {
    if (nrow(k_mat) != n || ncol(k_mat) != length(grid)) {
      stop("k_mat must be n x length(grid).")
    }
  }

  if (ci_method == "asymptotic") {
    if (is.null(k_mat)) {
      u_mat <- outer(y, grid, function(yy, xx) (yy - xx) / h1)
      k_mat <- kernel_eval(u_mat, kernel)
    }
    f_hat <- colMeans(k_mat) / h1
    Rk     <- kernel_Rk(kernel)
    se_hat <- sqrt(pmax(f_hat, 0) * Rk / (n * h1))

    z        <- qnorm(1 - alpha / 2)
    ci_lower <- pmax(f_hat - z * se_hat, 0)
    ci_upper <- f_hat + z * se_hat

    return(data.table(
      x          = grid,
      f_hat      = f_hat,
      se_hat     = se_hat,
      ci_lower   = ci_lower,
      ci_upper   = ci_upper,
      bandwidth  = h1,
      bandwidth2 = NA_real_
    ))
  }

  if (ci_method %in% c("wc", "bc", "sni") && kernel == "gaussian") {
    stop("wc/bc/sni CIs require compact-support kernels; drop gaussian.")
  }
  if (ci_method == "wc" && is.null(lipschitz_L)) {
    stop("wc CI requires lipschitz_L.")
  }
  if (ci_method == "sni" && (sni_u_grid_n %% 2) == 0) {
    stop("sni_u_grid_n must be odd so that u_grid includes 0.")
  }

  if (ci_method == "clt") {
    if (is.null(k_mat)) {
      u_mat <- outer(y, grid, function(yy, xx) (yy - xx) / h1)
      k_mat <- kernel_eval(u_mat, kernel)
    }

    f_hat <- colMeans(k_mat) / h1
    k_mean <- colMeans(k_mat)
    k_sq_mean <- colMeans(k_mat^2)
    s2 <- (n / (n - 1)) * pmax(k_sq_mean - k_mean^2, 0)
    se_hat <- sqrt(s2 / (n * h1^2))

    z        <- qnorm(1 - alpha / 2)
    ci_lower <- f_hat - z * se_hat
    ci_upper <- f_hat + z * se_hat

    return(data.table(
      x          = grid,
      f_hat      = f_hat,
      se_hat     = se_hat,
      ci_lower   = ci_lower,
      ci_upper   = ci_upper,
      bandwidth  = h1,
      bandwidth2 = NA_real_
    ))
  }

  # Kernel-matrix estimates at h1
  if (is.null(k_mat)) {
    u_mat <- outer(y, grid, function(yy, xx) (yy - xx) / h1)
    k_mat <- kernel_eval(u_mat, kernel)
  }
  f_hat <- colMeans(k_mat) / h1
  k_bar <- colMeans(k_mat)
  v_hat <- colSums((k_mat - matrix(k_bar, nrow = n, ncol = length(grid), byrow = TRUE))^2) / (n * h1)
  r_n   <- compute_Rn(v_hat, n, h1, alpha)

  if (ci_method == "wc") {
    ci_lower <- f_hat - lipschitz_L * h1 - r_n
    ci_upper <- f_hat + lipschitz_L * h1 + r_n

    return(data.table(
      x          = grid,
      f_hat      = f_hat,
      se_hat     = NA_real_,
      ci_lower   = ci_lower,
      ci_upper   = ci_upper,
      bandwidth  = h1,
      bandwidth2 = NA_real_
    ))
  }

  if (ci_method == "sni") {
    u_grid    <- seq(-1, 1, length.out = sni_u_grid_n)
    du        <- u_grid[2] - u_grid[1]
    k_u       <- kernel_eval(u_grid, kernel)

    # Use kernel-matrix estimates for SNI bias (no density interpolation)
    f_hat_km <- f_hat
    b_hat    <- numeric(length(grid))

    for (ii in seq_along(grid)) {
      x0 <- grid[ii]
      xq <- x0 + u_grid * h1
      u_mat_xq <- outer(y, xq, function(yy, xx) (yy - xx) / h1)
      k_mat_xq <- kernel_eval(u_mat_xq, kernel)
      f_hat_xq <- colMeans(k_mat_xq) / h1
      b_hat[ii] <- sum(k_u * (f_hat_xq - f_hat_km[ii])) * du
    }

    r_sni <- compute_Rsni(v_hat, n, h1, alpha, sni_t_multiplier, sni_log_n_coeff)
    ci_lower <- f_hat_km - b_hat - r_sni
    ci_upper <- f_hat_km - b_hat + r_sni

    return(data.table(
      x          = grid,
      f_hat      = f_hat_km,
      se_hat     = NA_real_,
      ci_lower   = ci_lower,
      ci_upper   = ci_upper,
      bandwidth  = h1,
      bandwidth2 = NA_real_
    ))
  }

  # Bias-corrected CI
  h2 <- h1 * bc_h2_ratio
  if (!is.finite(h2) || h2 <= 0) {
    stop("bc CI requires a positive h2; check bc_h2_ratio.")
  }

  u_mat2 <- outer(y, grid, function(yy, xx) (yy - xx) / h2)
  k_mat2 <- kernel_eval(u_mat2, kernel)
  ftilde_x0 <- colMeans(k_mat2) / h2

  u_grid     <- seq(-1, 1, length.out = bc_u_grid_n)
  du         <- u_grid[2] - u_grid[1]
  k_u        <- kernel_eval(u_grid, kernel)
  b_hat      <- numeric(length(grid))

  for (ii in seq_along(grid)) {
    x0 <- grid[ii]
    xq <- x0 + u_grid * h1
    u_mat_xq2 <- outer(y, xq, function(yy, xx) (yy - xx) / h2)
    k_mat_xq2 <- kernel_eval(u_mat_xq2, kernel)
    ftilde_vals <- colMeans(k_mat_xq2) / h2
    b_hat[ii] <- sum(k_u * (ftilde_vals - ftilde_x0[ii])) * du
  }

  ci_lower <- f_hat - b_hat - r_n
  ci_upper <- f_hat - b_hat + r_n

  data.table(
    x          = grid,
    f_hat      = f_hat,
    se_hat     = NA_real_,
    ci_lower   = ci_lower,
    ci_upper   = ci_upper,
    bandwidth  = h1,
    bandwidth2 = h2
  )
}
