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
