# Estimator: sample mean with t-based confidence interval

#' Estimate the mean and construct a confidence interval
#'
#' @param y Numeric vector of observations
#' @param alpha Significance level (default 0.05 for 95% CI)
#' @return Named list with estimate, ci_lower, ci_upper
estimate <- function(y, alpha = 0.05) {
  n <- length(y)
  mu_hat <- mean(y)
  se <- sd(y) / sqrt(n)
  q <- qt(1 - alpha / 2, df = n - 1)

  list(
    estimate = mu_hat,
    ci_lower = mu_hat - q * se,
    ci_upper = mu_hat + q * se
  )
}
