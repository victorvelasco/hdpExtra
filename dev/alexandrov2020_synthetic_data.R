library(hdp)
library(coda)
library(salso)

M <- read.csv("/home/victor/phd/data/i.syn.pancreas/sp.sp/ground.truth.syn.catalog.csv", header = TRUE)
M <- M[, -c(1,2)]
rownames(M) <- colnames(mut_count)
sum(M)
M <- t(M)

# Initialise
hdp_mut <- hdp_init(ppindex = c(0, rep(1, nrow(M))), # index of parental nodes
                    cpindex = c(1, rep(2, nrow(M))), # index of concentration param
                    hh=rep(1, ncol(M)), # prior is uniform over the 96 mutation categories
                    alphaa=rep(1, 2), # shape hyperparams for five different CPs
                    alphab=rep(1, 2)) # rate hyperparams for five different CPs

# add data to leaf nodes (one per cancer sample, in row order of mut_count)
hdp_mut <- hdp_setdata(hdp_mut,
                       dpindex = 2:numdp(hdp_mut), # index of nodes to add data to
                       as.data.frame(M)) # input data (mutation counts, sample rows match up with specified dpindex)

hdp_mut

hdp_extra_chains <- vector("list", 2)
### for (i in 1:4){
###
###   # activate DPs, 10 initial components
###   hdp_activated <- dp_activate(hdp_mut, 1:numdp(hdp_mut), initcc=30, seed=i*200)
###
###   hdp_extra_chains[[i]] <- hdpExtra_posterior(hdp_activated,
###                                               burnin=100,
###                                               n=10,
###                                               space=10,
###                                               cpiter=1,
###                                               seed=i*1e3)
### }
for (i in 1:2) {
  hdp_extra_chains[[i]] <- readRDS(paste0("dev/output/alexandrov2020_syn_i/hdpExtra_chain", i, ".rds"))
}
hdp_extra_chains <- HdpExtraChainMulti(hdp_extra_chains)

Sys.time()
allocations_best <- salso(t(hdp_allocations(hdp_extra_chains)), loss = VI(a = 2/3))
Sys.time()
Phi_best <- hdp_postprocessing(hdp_extra_chains, allocations_best)
hdp_plot_sig_uncertainty(Phi_best, "plots")

hdpExtra_diagnostic_plots(hdp_extra_chains, c("gamma", "alpha"))


