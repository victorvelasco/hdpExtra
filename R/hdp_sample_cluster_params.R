
hdp_sample_cluster_params <- function(classqq) {
  K <- ncol(classqq) # Number of clusters
  V <- nrow(classqq) # Number of categories
  Phi <- matrix(0, nrow = V, ncol = K)
  for (k in 1:K) {
    Phi[, k] <- MCMCpack::rdirichlet(1, classqq[, k] + 1)
  }

  return(Phi)
}
