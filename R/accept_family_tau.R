
#' Curate Quantile Specific Acceptable Families
#'
#' @param post_Q_tau posterior samples of the conditional quantiles at \code{XX}
#' @param XX covariate values
#' @param indicators inclusion indicators from BBA algorithn
#' @param epsmm_level slack variable
#' @param eta_level slack variable
#' @return \eqn{\tau}-specific acceptable families
#' @export
#'
accept_family_tau = function(post_Q_tau,
                             XX,
                             indicators,
                             eps_level = 0.05,
                             eta_level = 0.00

){


  p_loss = post_loss_l2(post_Q = post_Q,
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

