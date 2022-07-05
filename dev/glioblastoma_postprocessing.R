library(hdp)
library(salso)

chains <- list()
for (i in 1:4) {
  chains[[i]] <- readRDS(paste0("dev/output/glioblastoma/hdpExtra_chain", i, ".rds"))
}
chains <- HdpExtraChainMulti(chains)

hdpExtra_diagnostic_plots(chains, c("gamma","alpha"))
plot(1:2000,  chains[[1]]@cp_values[,1], 'l')
lines(1:2000, chains[[2]]@cp_values[,1], 'l', col = 'red')
lines(1:2000, chains[[3]]@cp_values[,1], 'l', col = 'red')

hdp_mut <- readRDS("dev/output/glioblastoma/hdp_mut.rds")

hdp_chain(chains)
