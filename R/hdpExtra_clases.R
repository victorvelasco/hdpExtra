#' An S4 class to store an MCMC chain
#'
#' @slot hdpChain An hdpSampleChain object for integration with the hdp package
#' @slot allocations A matrix where each column is a partition of the
#' observations sampled from the posterior distribution
#' @slot Phi A list whose elements are draws from the posterior over signature parameters
#' @importClassesFrom hdp hdpSampleChain
#' @export
setClass("HdpExtraChain",
  slots = c(
    hdpChain = "hdpSampleChain",
    allocations = "matrix",
    Phi = "list"
  )
)

#' An S4 class to store multiple MCMC chains
#'
#' @slot hdpChains An hdpSampleMulti object for integration with the hdp package
#' @slot allocations A matrix where each column is a partition of the
#' observations sampled from the posterior distribution
#' @slot Phi A list whose elements are draws from the posterior over signature parameters
#' @importClassesFrom hdp hdpSampleMulti
#' @export
setClass("HdpExtraChainMulti",
  slots = c(
    hdpChains = "hdpSampleMulti",
    allocations = "matrix",
    Phi = "list"
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
    allocations = allocations,
    Phi = Phi
  )
}

#' @export
setGeneric("hdp_chain", function(x) standardGeneric("hdp_chain"))

#' @export
setMethod("hdp_chain", "HdpExtraChain", function(x) {
  x@hdpChain
})

#' @export
setMethod("hdp_chain", "HdpExtraChain", function(x) {
  hdp_multi_chain(x@hdpChain)
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
