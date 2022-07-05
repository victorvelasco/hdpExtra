# posterior MCMC sampling for hdp
# returns an R list with two elements: 
# the updated hdp state as a list, and
# the vector of likelihoods for these iterations

#' @useDynLib hdpExtra hdpMultinomial_iterate
iterate <- function(hdplist, numiter, cpiter, allocations, numcp, verbosity){
  out <- .Call(hdpMultinomial_iterate, hdplist, numiter, cpiter,
               dolik=1, allocations, numcp, verbosity, PACKAGE="hdpExtra")
  return(out)
}
