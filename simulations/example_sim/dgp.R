# Data-generating process: draw from N(mu, sigma^2)

#' Generate a sample from the normal DGP
#'
#' @param n_obs Sample size
#' @param mu_true True population mean
#' @param sigma True population standard deviation
#' @return Numeric vector of length n_obs
generate_data <- function(n_obs, mu_true, sigma) {
  rnorm(n_obs, mean = mu_true, sd = sigma)
}
