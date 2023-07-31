

#' Compute the check loss between response \code{y} and quantile-predictions \code{qhat} for any quantile \code{tau}
#'
#' @param y response values
#' @param qhat quantile predictions
#' @param tau quantile of interest
#'
#' @return the average check loss across all quantile predictions
#' @export

check_loss<-function(y, qhat, tau){
  mu = y - qhat

  mean(mu*(tau - I(mu<0)))
}
