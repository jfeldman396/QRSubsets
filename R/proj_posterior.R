#' Projected predictive distribution for regression coefficients
#'
#' Given draws from the predictive distribution, project these
#' draws onto (a  subset of) the covariates. This produces
#' many predictive draws for the regression coefficients,
#' which provides uncertainty quantification.
#' @param post_y_pred \code{S x n} matrix of posterior predictive
#' at the given \code{XX} covariate values
#' @param XX \code{n x p} matrix of covariates
#' @param sub_x vector of inclusion indicators for the \code{p} covariates;
#' the remaining coefficients will be fixed at zero
#' @param use_ols logical; if TRUE, use ordinary least squares regression (default);
#' otherwise use logistic regression
#' @return \code{post_beta}: the \code{S x p} matrix of
#' draws from the projected predictive distribution of  the regression coefficients.
#' @export

proj_posterior = function(post_Q_tau, XX, sub_x = 1:ncol(XX), use_ols = TRUE){

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
  if(use_ols){
    post_beta[,sub_x] = tcrossprod(post_Q_tau, t(XX_nz)) %*%
      chol2inv(chol(crossprod(XX_nz)))
  } else{
    post_beta[,sub_x] = t(sapply(1:S, function(s){
      suppressWarnings(
        coef(
          glm(post_Q_tau[s,] ~ XX_nz - 1, family = binomial())
        )
      )
    }))
  }

  return(post_beta)
}
