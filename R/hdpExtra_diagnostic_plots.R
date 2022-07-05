#' Diagnostic plots for the concentration parameters
#'
#' @param hdpExtraChainMulti An object of class HdpExtraChainMulti
#' @param param_names Name of the concentration parameters, to display in
#' the diagnostic plots
#'
#' @export
hdpExtra_diagnostic_plots <- function(hdpExtraChainMulti, param_names) {
  mcmc_list_object <- coda::as.mcmc.list(hdpExtraChainMulti@cp_values)
  coda::varnames(mcmc_list_object) <- param_names
  ggs_object <- ggmcmc::ggs(mcmc_list_object, keep_original_order = TRUE)
  list(
    ggmcmc::ggs_traceplot(ggs_object),
    ggmcmc::ggs_compare_partial(ggs_object)
  )
}
