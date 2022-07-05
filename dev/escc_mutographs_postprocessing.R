library(hdp)
library(salso)

chains <- list()
for (i in 0:2) {
  chains[[i+1]] <- readRDS(paste0("dev/output/hdpExtra_chain", i, ".rds"))
}
chains <- HdpExtraChainMulti(chains)

hdpExtra_diagnostic_plots(chains)
