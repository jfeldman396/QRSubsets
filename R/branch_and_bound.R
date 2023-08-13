#' Branch-and-bound algorithm for linear subset search
#'
#' Search for the "best" (according to residual sum of squares)
#' linear subsets of each size. The algorithm may collect
#' the \code{n_best} "best" subsets of each size, include or
#' exclude certain variables automatically, and apply
#' forward, backward, or exhaustive search.
#'
#' @param Q_hat_tau posterior mean of the \eqn{\tau}-conditional quantile function at \code{XX}
#' @param XX matrix of covariates
#' @param wts vector of observation weights (for weighted least squares)
#' @param n_best number of "best" subsets for each model size
#' @param to_include indices of covariates to include in *all* subsets
#' @param to_exclude indices of covariates to exclude from *all* subsets
#' @param searchtype use exhaustive search, forward selection, backward selection or sequential replacement to search
#' @return \code{inclusion_index}: the matrix of inclusion indicators (columns) for
#' each subset returned (rows)
#' @importFrom leaps regsubsets
#' @export

branch_and_bound = function(Q_hat_tau,
                            XX,
                            wts = NULL,
                            n_best = 15,
                            to_include = 1,
                            to_exclude = NULL,
                            searchtype = 'exhaustive'
){

  # Get dimensions of the covariate matrix:
  n = nrow(XX); p = ncol(XX)

  # Useful step:
  colnames(XX) = 1:p

  # Some basic checks:
  if(length(Q_hat_tau) !=n)
    stop('length of yy must equal the number of rows of XX')

  if(!is.null(to_include) && !is.null(to_exclude)){
    # Quick check:
    if(any(!is.na(match(to_include, to_exclude))))
      stop('Cannot include and exclude the same variables!')

    # Reindex to account for the excluded terms:
    to_include = match(to_include, (1:p)[-to_exclude])
  }

  # And a warning if the number of predictors is too large
  if(p - length(to_exclude) > 40 && searchtype == 'exhaustive'){
    warning("Inadvisable number of predictors for exhaustive search: slow computing! Try pre-screening to exclude some variables.")
  }

  if(!is.null(wts)){
    if(length(wts) != n)
      stop('wts must have length n')
  } else wts = rep(1,n)

  # Delete the excluded columns, if any:
  if(!is.null(to_exclude))
    XX = XX[,-to_exclude]

  # Branch-and-bound search:
  fit_all = leaps::regsubsets(x = XX, y = Q_hat_tau,
                       weights = wts,
                       nbest = n_best, nvmax = p,
                       method = searchtype,
                       intercept = FALSE,
                       really.big = TRUE,
                       force.in = to_include)

  # Indicator matrix of variable inclusions for each subset:
  temp_inclusion_index = summary(fit_all)$which

  # Inclusion index: need to adjust for excluded variables
  inclusion_index = matrix(FALSE,   # FALSE = excluded
                           nrow = nrow(temp_inclusion_index),
                           ncol = p)

  # Update the non-excluded entries:
  inclusion_index[,match(colnames(temp_inclusion_index),
                         1:p)] = temp_inclusion_index

  return(inclusion_index)
}
