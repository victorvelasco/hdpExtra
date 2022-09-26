library(hdp)
library(salso)
library(dplyr)

M <- read.csv("data/WGS_PCAWG.96.csv", header = TRUE, row.names = 1)
M <- M[, startsWith(colnames(M), "CNS.GBM")] # Keep glioblastoma samples only
M <- M %>% as.matrix %>% t %>% as.data.frame

chains <- list()
for (i in 1:4) {
  chains[[i]] <- readRDS(paste0("output/glioblastoma/hdpExtra_chain", i, ".rds"))
}
chains <- HdpExtraChainMulti(chains)

hdpExtra_diagnostic_plots(chains, c("gamma","alpha"))
plot(1:2000,  chains[[1]]@cp_values[,1], 'l')
lines(1:2000, chains[[2]]@cp_values[,1], 'l', col = 'red')
lines(1:2000, chains[[3]]@cp_values[,1], 'l', col = 'red')

for (i in c("1", "23", "25")) {
  salso_obj <- readRDS(paste0("output/glioblastoma/salso_VI", i, ".rds"))
  sort(table(salso_obj))
  hdp_extra_chains <- chains
  ## Phi_best <- hdp_postprocessing(hdp_extra_chains, salso_obj)
  ## hdp_plot_sig_uncertainty(Phi_best, "plots")


  J <- nrow(M)
  K <- max(salso_obj)

  nmuts <- rowSums(M)
  index_end <- cumsum(nmuts)
  index_beg <- c(0, index_end[1:(length(nmuts)-1)])+1
  names(index_beg) <- names(index_end)

  sort(table(salso_obj), decreasing = TRUE)

  exposures <- matrix(0, nrow = J, ncol = K)
  rownames(exposures) <- rownames(M)
  colnames(exposures) <- paste("Sig", 01:K)

  for (j in 1:J) {
    patient_allocs <- table(salso_obj[index_beg[j]:index_end[j]])
    exposures[j, as.integer(names(patient_allocs))] <- unname(patient_allocs)
  }

  write.csv(exposures, paste0("output/glioblastoma/exposures_VI", i, ".csv"))
}



Phi_best <- hdp_postprocessing(hdp_extra_chains, salso_obj)
hdp_plot_sig_uncertainty(Phi_best, "plots")
