#' @title KL Divergence for Two Continuous Prob Distributions
#' @param
#' @description
#' @details
#' @returns
#' @export

KLdivpqcont <- function(p_empdat, q_empdat, N = N) {
  # checks
  goodegg::assert_numeric(p_empdat)
  goodegg::assert_numeric(q_empdat)
  # Estimate the pdfs using kernel density estimation but bound by 0 and N
  p_pdf_estimate <- density(p_empdat, from = 0, to = N, kernel = "epanechnikov", bw = 1)
  q_pdf_estimate <- density(q_empdat, from = 0, to = N, kernel = "epanechnikov", bw = 1) # bw = 1 for ind

  # empiric pdf fxn for P and Q via interpolation for estimated pdfs
  p_pdf <- function(x) { approx(p_pdf_estimate$x, p_pdf_estimate$y, xout = x)$y }
  q_pdf <- function(x) { approx(q_pdf_estimate$x, q_pdf_estimate$y, xout = x)$y }

  # Define the integrand for the KL divergence
  kl_integrand <- function(x) {
    p <- p_pdf(x) + .Machine$double.eps
    q <- q_pdf(x) + .Machine$double.eps
    p * log(p / q)
    }

  # Calculate the KL divergence by numerical integration and bound by 0:N again
  kl_divergence <- cubature::adaptIntegrate(kl_integrand, lower = 0, upper = N)

  # out
  return(kl_divergence)
}



