
#' Variable importance for each quantile-specific acceptable family
#'
#' Given the candidate subsets and the indicators of any quantile-specific acceptable family
#' of subsets, compute for each variable the proportion of acceptable subsets in
#' which that variable appears. If specified, variable co-appearances
#' can be computed as reported as well.
#' @param indicators \code{L x p} matrix of inclusion indicators (booleans)
#' where each row denotes a candidate subset (from \code{branch_and_bound})
#' @param all_accept indices (i.e., rows of \code{indicators})
#' that correspond to the acceptable subsets (from \code{accept_family_tau})
#' @param co logical; if TRUE, compute and return the co-variable importances
#' @param xnames the names of the x-variables
#' @return a list with the variable importances \code{vi_inc} and the co-variable
#' importances \code{vi_co}
#' @export
var_imp = function(indicators, all_accept, co = TRUE, xnames = NULL){

  # Number of covariates:
  p = ncol(indicators)

  if(length(all_accept) > nrow(indicators))
    stop('Cannot have more acceptable subsets than candidate subsets!')

  if(!is.null(xnames)){
    if(length(xnames) != p)
      stop('xnames must have length p')
  } else xnames = paste(1:p)


  # For each variable, how often is it included in an acceptable subset?
  vi_inc = colMeans(indicators[all_accept,]); names(vi_inc) = xnames

  if(co){
    vi_co = diag(vi_inc, p)
    colnames(vi_co) = rownames(vi_co) = xnames
    for(j1 in 1:(p-1)){
      for(j2 in (j1 + 1):p){
        vi_co[j1,j2] = vi_co[j2,j1] =
          mean(indicators[all_accept, j1] & indicators[all_accept, j2])
      }
    }

  } else vi_co = NULL

  return(
    list(
      vi_inc = vi_inc,
      vi_co = vi_co
    )
  )
}
