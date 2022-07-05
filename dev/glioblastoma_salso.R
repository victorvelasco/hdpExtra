library(hdpExtra)
library(hdp)
library(salso)

chains <- list()
for (i in 1:4) {
  chains[[i]] <- readRDS(paste0("output/glioblastoma/hdpExtra_chain", i, ".rds"))
}
chains <- HdpExtraChainMulti(chains)
rm(chains)
allocations <- hdp_allocations(chains)
allocations <- t(allocations)
salso_object <- salso(allocations, loss = VI(a= 2/3))

saveRDS(salso_object, "output/glioblastoma/salso_VI23.rds")
