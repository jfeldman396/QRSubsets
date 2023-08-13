#' Projected posterior distribution for regression coefficients, i.e. the posterior action
#'
#' Given draws from the posterior distribution of the conditional quantile function, project these
#' draws onto (a subset of) the covariates. This produces
#' many posterior draws for the regression coefficients,
#' which provides uncertainty quantification.
#' @param post_Q_tau \code{S x n} matrix of posterior samples of the conditional quantile function at each design point among \code{XX} covariate values
#' @param XX \code{n x p} matrix of covariates
#' @param sub_x vector of inclusion indicators for the \code{p} covariates;
#' the remaining coefficients will be fixed at zero

#' @return \code{post_beta}: the \code{S x p} matrix of
#' draws from the projected posterior distribution of the regression coefficients for any subset.
#' @export

proj_posterior = function(post_Q_tau, XX, sub_x = 1:ncol(XX)){

  # Get dimensions of the covariate matrix:
  n = nrow(XX); p = ncol(XX)
  S = nrow(post_Q_tau) # number of posterior simulations

  if(ncol(post_Q_tau) != n)
    stop('number of columns of post_Q_tau must equal the number of rows of XX')

  if(length(sub_x) > p)
    stop('length of sub_x must be less than or equal to p')

  # Subset the columns of XX:
  XX_nz = XX[, sub_x]

  # Storage: include zeros as needed!
  post_beta = array(0, c(S, p))
    post_beta[,sub_x] = tcrossprod(post_Q_tau, t(XX_nz)) %*%
      chol2inv(chol(crossprod(XX_nz)))

  return(post_beta)
}
