#' Compute the posterior aggregated squared-error loss
#'
#' Use posterior draws of the conditional quantile function to compute the residual sum
#' of squares (RSS) between any subset of predictors and the model-fitted quantiles. These quantities
#' are computed given filtering of subsets from the branch and bound search algorithm, and the
#' and posterior samples of conditional quantile function from the overarching Bayesian model fit to
#' the data.
#'
#' @param post_Q_tau \code{S x n} matrix of S posterior samples of the conditional quantile function
#' \eqn{Q_{\tau}(x_i)} for each covariate value in \code{XX}
#' @param XX \code{n x p} matrix of covariates at which to evaluate subsets
#' @param indicators \code{L x p} matrix of inclusion indicators extracted from the branch and bound
#' search algorithm. Each row corresponds to a candidate subset, with the indicator determining whether \code{X_j} is included.
#' @return a list with four elements \code{post_loss}, \code{post_loss_raw},
#' \code{ref}, and \code{beta_hat_L}
#'
#' @details \code{post_loss} is an \code{S\times L} matrix including posterior samples of the difference in RSS between the anchor,
#' which is the posterior mean of the conditional quantile function under the unifying Bayesian model, and each subset from
#' the branch and bound filtering. This is used to determine the acceptable family.
#' \code{post_loss_raw} is an \code{S\times L} matrix with the posterior samples of the RSS for each subset in the branch and bound filtering.
#' \code{ref} is an \code{S x 1} array containing posterior samples of the RSS for the anchor action. \code{post_loss_raw} and \code{ref} are
#' the constituent terms to compute \code{post_loss}. Finally,
#' \code{beta_hat_L} is the optimal action for each subset in the branch and bound filtering
#' @export
#'

post_loss_l2 = function(post_Q_tau,
                        XX,
                        indicators

){

  # Get dimensions of the covariate matrix:
  n = nrow(XX); p = ncol(XX)

  # And some other dimensions
  S = nrow(post_Q_tau) # number of posterior simulations
  L = nrow(indicators)  # number of subsets to consider

  if(ncol(post_Q_tau) != n)
    stop('number of columns of post_Q_tau must equal the number of rows of XX')

  if(ncol(indicators) != p)
    stop('indicators must have the same number of columns as XX has rows')

  # Fitted values:
  Q_hat = colMeans(post_Q_tau)

  # Useful terms:
  XtX = crossprod(XX) # X'X for these testing points
  post_Q_tau2 = rowMeans(post_Q_tau^2) # squared predictive variables summed across n_out

  # Storage:
  post_loss = pred_loss1 =array(0, c(S, L))
  beta_hat = array(0, c(L,p))
  for(ell in 1:L){
    # Coefficients for this model:
    beta_ell = rep(0, p)
    beta_ell[indicators[ell,]] = coef(lm(Q_hat ~ XX[,indicators[ell,]] - 1))
    beta_hat[ell,] = beta_ell
    # And predictive loss:
    Q_pred<- XX%*%beta_ell

    post_loss[,ell] = post_Q_tau2+
      1/n*as.numeric(crossprod(Q_pred)) -
      2/n*tcrossprod(post_Q_tau, t(Q_pred))



  }

  # Reference: Q_hat and raw posterior loss
  post_loss_raw = post_loss

  ref<- post_Q_tau2 + 1/n*as.numeric(crossprod(Q_hat)) - 2/n*tcrossprod(post_Q_tau,t(Q_hat))


  # Percent difference in predictive loss relative to anchor Q_hat:
  post_loss = apply(post_loss, 2, function(ploss)
    100*(ploss - ref)/ref)


  return(list(
    post_loss = post_loss,
    post_loss_raw = post_loss_raw,
    ref = ref,
    beta_hat_L = beta_hat)
  )
}

