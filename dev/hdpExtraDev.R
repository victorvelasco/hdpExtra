##### Script written with the purpose of developing hdpExtra
##### Not to be built into the package.

library(hdp)
library(coda)
library(salso)

dim(mut_count)
M <- mut_count[1:10, ]
M <- M + 1
sum(M)

# Initialise
hdp_mut <- hdp_init(ppindex = c(0, rep(1, nrow(M))), # index of parental nodes
                    cpindex = c(1, rep(2, nrow(M))), # index of concentration param
                    hh=rep(1, 96), # prior is uniform over the 96 mutation categories
                    alphaa=rep(1, 2), # shape hyperparams for five different CPs
                    alphab=rep(1, 2)) # rate hyperparams for five different CPs

# add data to leaf nodes (one per cancer sample, in row order of mut_count)
hdp_mut <- hdp_setdata(hdp_mut,
                       dpindex = 2:numdp(hdp_mut), # index of nodes to add data to
                       as.data.frame(M)) # input data (mutation counts, sample rows match up with specified dpindex)

hdp_mut

hdp_extra_chains <- vector("list", 4)
for (i in 1:4){

  # activate DPs, 10 initial components
  hdp_activated <- dp_activate(hdp_mut, 1:numdp(hdp_mut), initcc=10, seed=i*200)

  hdp_extra_chains[[i]] <- hdpExtra_posterior(hdp_activated,
                               burnin=50,
                               n=5,
                               space=200,
                               cpiter=3,
                               seed=i*1e3)
}

hdp_extra_chains <- HdpExtraChainMulti(hdp_extra_chains)

allocations_best <- salso(t(hdp_allocations(hdp_extra_chains)), loss = VI(a = 0.5))
Phi_best <- hdp_postprocessing(hdp_extra_chains, allocations_best)
hdp_plot_sig_uncertainty(Phi_best, "plots")


