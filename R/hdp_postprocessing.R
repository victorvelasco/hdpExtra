#' Posterior draws the parameters of the signatures in a representative partition
#'
#' The posterior draws are stored in a 3D array of dimension V x K x S, where
#' V is the number of categories, K is the number of clusters (signatures) and
#' S is the number of posterior samples
#' @param hdp_extra_chains An HdpExtraChainMulti object
#' @param allocations_best A vector of length N containing a representative partition
#' @return 3D array with posterior samples of the signatures the representative partition
#' @export
hdp_postprocessing <- function(hdp_extra_chains, allocations_best) {

  allocations <- hdp_allocations(hdp_extra_chains)
  Phi <- hdp_signature_parameters(hdp_extra_chains)

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

  return(Phi_best)
}
