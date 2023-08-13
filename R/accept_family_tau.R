
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
#' @param eps_level slack variable, used to determine member subsets in the acceptable family
#' @param eta_level slack variable, used to determine member subsets in the acceptable family
#' @return \eqn{\tau}-specific acceptable families
#'
#' @details For any subset \code{S} obtained by the quantile-specific branch and bound search,
#' that subset is deemed acceptable if \eqn{P(D^{\tau}_{S, \hat{Q}}  \leq \eta) \geq \epsilon}, where \eqn{D^{\tau}_{S, \hat{Q}} = 100 \times (L^{\tau}_{S}(\boldsymbol\theta)  -L^{\tau}_{\hat{Q}}(\boldsymbol\theta) )/L^{\tau}_{\hat{Q}}(\boldsymbol\theta) }
#' with \eqn{\hat Q} the posterior mean of \eqn{\tau}th conditional quantile function, \eqn{L} is the aggregated \eqn{L^2} loss, and \code{\theta} are Bayesian model \code{\mathcal{M}} parameters, drawn from the posterior.
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

  #index of max D_post
  max_prob_accept<- which.max(colMeans(post_loss <= eta_level)[all_accept])


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

  ell_loss = all_accept[max_prob_accept]
  # Compute the coefficients for each index:

  beta_hat_small = beta_hat[ell_small,]
  beta_hat_loss = beta_hat[ell_loss,]


  return(
    list(
      all_accept = all_accept,
      beta_hat_small = beta_hat_small,
      beta_hat_loss = beta_hat_loss,
      ell_small = ell_small,
      ell_loss= ell_loss
    )
  )
}

