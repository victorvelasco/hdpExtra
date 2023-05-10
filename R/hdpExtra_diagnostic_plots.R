#' Diagnostic plots for the concentration parameters
#'
#' @param hdpExtraChainMulti An object of class HdpExtraChainMulti
#' @param param_names Name of the concentration parameters, to display in
#' the diagnostic plots
#'
#' @export
hdpExtra_diagnostic_plots <- function(hdpExtraChainMulti, param_names) {
  list_cp_values <- coda::as.mcmc.list(hdpExtraChainMulti@cp_values)
  coda::varnames(list_cp_values) <- param_names
  list_cp_values <- ggmcmc::ggs(list_cp_values, keep_original_order = TRUE)

  list_nclust <- coda::as.mcmc.list(hdpExtraChainMulti@nclust)
  coda::varnames(list_nclust) <- "nclust"
  list_nclust <- ggmcmc::ggs(list_nclust, keep_original_order = TRUE)

  list(
    ggmcmc::ggs_traceplot(list_cp_values),
    ggmcmc::ggs_traceplot(list_nclust)
    # ggmcmc::ggs_compare_partial(ggs_object)
  )
}

#' Autocorrelation plots for the concentration parameters
#'
#' @param hdpExtraChainMulti An object of class HdpExtraChainMulti
#' @param param_names Name of the concentration parameters, to display in
#' the diagnostic plots
#'
#' @export
hdpExtra_acf <- function(hdpExtraChainMulti, param_names) {
  list_cp_values <- lapply(
    1:length(hdpExtraChainMulti@cp_values),
    function(i) coda::as.mcmc(hdpExtraChainMulti@cp_values[[i]][-(1:chlist@burnin[i]), ])
  )
  list_cp_values <- coda::as.mcmc.list(list_cp_values)
  coda::varnames(list_cp_values) <- param_names
  list_cp_values <- ggmcmc::ggs(list_cp_values, keep_original_order = TRUE)
  list(ggmcmc::ggs_autocorrelation(list_cp_values, nLags = 200))
}
