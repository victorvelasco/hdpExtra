#' allocations is an N x S matrix, (N = number of mutations, S = number of posterior iterations)
#'
#' @param allocations An N x S matrix
#' @param Phi A list of S matrices. Each matrix of dimension V x K^(s)
#' @param allocations_best A vector of length N containing a representative partition
#' @param alpha The probability coverage of the intervals (by default 0.95)
#' @return Point estimates of the parameters characterising a representative partition
#' @export
hdp_postprocessing <- function(allocations, Phi, allocations_best, alpha = 0.95) {

  V <- nrow(Phi[[1]]) # Number of categories
  clabels <- unique(sort(allocations_best)) # Clusters in the partition allocations_best (indices)
  K <- length(clabels) # Number of clusters in the partition allocations_best
  S <- length(Phi) # Number of posterior iterations
  N <- length(allocations_best) # Total number of observations
  N_k <- unname(table(allocations_best)) # Cluster sizes in the best partition

  Phi_best <- array(0, dim = c(V, K, S))
  for (s in 1:S) {
    Phi_s <- Phi[[s]][, allocations[, s]] # Matrux of dimension V x N
    for (k in 1:K) {
      Phi_best[, k, s] <- rowMeans(
        as.matrix(Phi_s[, allocations_best == clabels[k]])
      )
    }
  }

  #### Phi_best_summary <- array(-1,
  ####                           c(V, 3, K),
  ####                           list(NULL, c("median", "lower", "upper"), NULL))
  #### for (k in 1:K) {
  ####   for (v in 1:V) {
  ####     Phi_best_summary[v, , k] <- quantile(Phi_best[v, k, ], c(0.5, (1-alpha)/2, (alpha+1)/2))
  ####   }
  #### }

  return(Phi_best)
}
