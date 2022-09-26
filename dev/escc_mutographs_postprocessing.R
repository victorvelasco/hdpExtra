library(hdp)
library(salso)

chains <- list()
for (i in 1:2) {
  chains[[i]] <- readRDS(paste0("output/escc/hdpExtra_chain", i, ".rds"))
}
chains <- HdpExtraChainMulti(chains)

pdf(file = "dev/plots/diagnostics.pdf")
hdpExtra_diagnostic_plots(chains, c("gamma", "alpha0", "alpha1"))
dev.off()

salso_obj <- readRDS("output/escc/salso_VI23.rds")
Phi_best <- hdp_postprocessing(hdp_extra_chains, salso_obj)
pdf(file = "dev/plots/signatures_VI23.pdf")
hdp_plot_sig_uncertainty(Phi_best, "plots")
dev.off()
