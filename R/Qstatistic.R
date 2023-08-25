#' Q statistic
#'
#' Cochran’s Q statistic
#'
#' @param b_group a vector containing effect size estimates
#' @param se_group a vector containing standard errors of effect size estimates
#' @param b_overall overall effect size estimate
#' 
#' @return list containing \itemize{
#' \item \code{Q} Cochran’s Q statistic
#' }
#' 
#' @export
#' 

Qstatistic <- function(b_group, se_group, b_overall) {
  Q = sum(1 / se_group**2 * (b_group - b_overall)**2)
  return(Q)
}
