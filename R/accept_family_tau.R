
#' Curate Quantile Specific Acceptable Families
#'
#' Using posterior samples of the \eqn{\tau}th conditional quantile function, evaluate the posterior probability that each candidate
#' subset obtained from the branch and bound search algorithm has comparable predictive power to the anchoring action. The anchor is
#' chosen as the ``best" fit to the model-based conditional quantile function, which we set as the model-based posterior mean.
#'
#' @param post_Q_tau \code{S x n} matrix of S posterior samples of the conditional quantile function
#' \eqn{Q_{\tau}(x_i)} for each covariate value in \code{XX}
#' @param XX \code{n x p} matrix of covariates at which to evaluate subsets
#' @param indicators \code{L x p} matrix of inclusion indicators extracted from the branch and bound
#' @param eps_level slack variable \eqn{\epsilon}, used to determine member subsets in the acceptable family. Decrease to discover more acceptable subsets.
#' @param eta_level slack variable \eqn{\eta}, used to determine member subsets in the acceptable family. Increase to discover more acceptable subsets.
#' @return A list containing \code{all_accept} which are the indices of the subsets from \code{indicators} which are acceptable,
#' \code{ell_small}, which the index of the smallest acceptable subset, and \code{beta_hat_small} which is the optimal action (\eqn{\tau} quantile regression coefficients) for the smallest acceptable subset.
#'
#'
#' @details For any subset \eqn{S} obtained by the quantile-specific branch and bound search,
#' that subset is deemed acceptable if \eqn{P(D^{\tau}_{S, \hat{Q}}  \leq \eta) \geq \epsilon}, where \eqn{D^{\tau}_{S, \hat{Q}} = 100 \ \code{x} \ (L^{\tau}_{S}(\theta)  -L^{\tau}_{\hat{Q}}(\boldsymbol\theta) )/L^{\tau}_{\hat{Q}}(\theta) },
#'\eqn{\hat Q} is the posterior mean of \eqn{\tau}th conditional quantile function, \eqn{L} is the aggregated \eqn{L^2} loss, and \eqn{\theta} are the Bayesian regression model parameters, drawn from the posterior. The smallest acceptable subset
#' is then the smallest in cardanlity subset collected that maintains satisfactory predictive power under the criteria outlined
#' @export
#'
accept_family_tau = function(post_Q_tau,
                             XX,
                             indicators,
                             eps_level = 0.05,
                             eta_level = 0.00

){


  p_loss = post_loss_l2(post_Q_tau = post_Q_tau,
                        XX = XX,
                        indicators = indicators)
  post_loss = p_loss$post_loss
  insamp_loss = p_loss$insample_loss
  beta_hat = p_loss$beta_hat_L
  emp_loss = NULL

  # Indices of acceptable subsets:
  all_accept = which(colMeans(post_loss <= eta_level)
                     >= eps_level)

  prob_accept = colMeans(post_loss <= eta_level)[all_accept]


  # Subset sizes:
  subset_size = rowSums(indicators)

  # Minimize size of acceptable subset:
  min_size_accept = min(subset_size[all_accept])


  # Index of acceptable subsets with this size:
  ind_min_size_accept = which(subset_size[all_accept] == min_size_accept)

  # If more than one, select the minimum empirical loss (or expected predictive loss)
  if(length(ind_min_size_accept) > 1){
    ell_small = all_accept[ind_min_size_accept][which.max(prob_accept[ind_min_size_accept])]
  }else{
    ell_small = all_accept[ind_min_size_accept]
  }

  # Compute the coefficients for each index:

  beta_hat_small = beta_hat[ell_small,]


  return(
    list(
      all_accept = all_accept,
      beta_hat_small = beta_hat_small,
      ell_small = ell_small
    )
  )
}

