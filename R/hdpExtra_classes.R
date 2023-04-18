#' An S4 class to store an MCMC chain
#'
#' @slot hdpChain An object of class hdp::hdpSampleChain for integration with the hdp package
#' @slot allocations A matrix whose colums are partitions (as vectors of labels) of the
#' observations sampled from the posterior distribution and rows index iterations
#' of the MCMC sampler
#' @slot cp_values A matrix containing the values of the concentration parameters drawn
#' across iterations of the sampler
#' @slot Phi A list whose elements are draws from the posterior over signature parameters
#' @slot niter Number of MCMC draws that are saved to memory
#' @slot burnin Number of MCMC draws that are discarded in the burnin period
#' @slot thin Thinning parameter
#' @slot nclust Number of clusters across iterations of the sampler
#' @importClassesFrom hdp hdpSampleChain
#' @export
setClass("HdpExtraChain",
  slots = c(
    hdpChain = "hdpSampleChain",
    cp_values = "matrix",
    allocations = "matrix",
    Phi = "list",
    niter = "numeric",
    burnin = "numeric",
    thin = "numeric",
    nclust = "numeric"
  )
)

#' An S4 class to store multiple MCMC chains
#'
#' @slot hdpChains An hdp::hdpSampleMulti object for integration with the hdp package
#' @slot cp_values A list of matrices containing the values of the concentration parameters drawn
#' across iterations of the sampler
#' @slot allocations A matrix where each column is a partition of the
#' observations sampled from the posterior distribution
#' @slot Phi A list whose elements are draws from the posterior over signature parameters
#' @slot niter Number of MCMC draws that are saved to memory
#' @slot burnin Number of MCMC draws that are discarded in the burnin period
#' @slot thin Thinning parameter
#' @slot nclust Number of clusters across iterations of the sampler
#' @importClassesFrom hdp hdpSampleMulti
#' @export
setClass("HdpExtraChainMulti",
  slots = c(
    hdpChains = "hdpSampleMulti",
    cp_values = "list",
    allocations = "matrix",
    Phi = "list",
    niter = "numeric",
    burnin = "numeric",
    thin = "numeric",
    nclust = "numeric"
  )
)


#' Creates an object of class HdpExtraChainMulti
#'
#' @param hdpExtraChains A list of objects of class HdpExtraChain
#' @export
HdpExtraChainMulti <- function(hdpExtraChains) {
  allocations = hdpExtraChains[[1]]@allocations
  Phi = hdpExtraChains[[1]]@Phi

  for (i in 2:length(hdpExtraChains)) {
    allocations <- cbind(allocations, hdpExtraChains[[i]]@allocations)
    Phi <- append(Phi, hdpExtraChains[[i]]@Phi)
  }

  new("HdpExtraChainMulti",
    hdpChains = hdp_multi_chain(lapply(1:length(hdpExtraChains), function(i) hdpExtraChains[[i]]@hdpChain)),
    cp_values = lapply(hdpExtraChains,
                  function(chain) coda::as.mcmc(chain@cp_values)),
    allocations = allocations,
    Phi = Phi,
    niter = sapply(hdpExtraChains, function(x) x@niter),
    burnin = sapply(hdpExtraChains, function(x) x@burnin),
    thin = sapply(hdpExtraChains, function(x) x@thin)
  )
}

#' @export
setGeneric("hdp_chain", function(x) standardGeneric("hdp_chain"))

#' @export
setMethod("hdp_chain", "HdpExtraChain", function(x) {
  x@hdpChain
})

#' @export
setMethod("hdp_chain", "HdpExtraChainMulti", function(x) {
  hdp::chains(x@hdpChains)
})

#' @export
setGeneric("hdp_allocations", function(x) standardGeneric("hdp_allocations"))

#' @export
setMethod("hdp_allocations", "HdpExtraChain", function(x) {
  x@allocations
})

#' @export
setMethod("hdp_allocations", "HdpExtraChainMulti", function(x) {
  x@allocations
})

#' @export
setGeneric("hdp_signature_parameters", function(x) standardGeneric("hdp_signature_parameters"))

#' @export
setMethod("hdp_signature_parameters", "HdpExtraChain", function(x) {
  x@Phi
})

#' @export
setMethod("hdp_signature_parameters", "HdpExtraChainMulti", function(x) {
  x@Phi
})
