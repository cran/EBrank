#' EBrank: Empirical Bayes Ranking
#'
#' Empirical Bayes ranking applicable to parallel-estimation settings where the estimated parameters are asymptotically unbiased and normal, with known standard errors.  A mixture normal prior for the parameter is estimated, subsequentially ranks for each parameter are simulated from the resulting posterior. Finally, experiments are ordered by expected posterior  rank, although computations minimizing other plausible rank-loss functions are also given.  
#' 
#' @section EBrank functions:
#' rankEM
#'
#' @docType package
#' @name EBrank
NULL
#> NULL