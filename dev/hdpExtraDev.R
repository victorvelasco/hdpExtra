##### Script written with the purpose of developing hdpExtra
##### Not to be built into the package.

library(hdp)
library(hdpExtra)
library(coda)

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
                       M) # input data (mutation counts, sample rows match up with specified dpindex)

hdp_mut

chlist <- vector("list", 4)
allocations <- matrix(-1, nrow = sum(M), ncol = 0)
Phi <- list()
for (i in 1:4){

  # activate DPs, 10 initial components
  hdp_activated <- dp_activate(hdp_mut, 1:numdp(hdp_mut), initcc=10, seed=i*200)

  posterior_draws <- hdpExtra_posterior(hdp_activated,
                               burnin=50,
                               n=5,
                               space=200,
                               cpiter=3,
                               seed=i*1e3)
  chlist[[i]] <- posterior_draws$chains
  allocations <- cbind(allocations, posterior_draws$allocations)
  Phi <- append(Phi, posterior_draws$Phi)
}

allocations_best <- salso::salso(t(allocations), loss = "binder", maxNClusters = 4)
Phi_best <- hdp_postprocessing(allocations, Phi, allocations_best)
hdp_plot_sig_uncertainty(Phi_best, "plots")

class(Phi_best)
dim(Phi_best)


Phi_mcmc <- t(Phi_best[, 3, ])
dimnames(Phi_mcmc) <- list(
  NULL, channels
)
Phi_mcmc <- as.mcmc(Phi_mcmc)
Phi_mcmc_list <- as.mcmc.list(Phi_mcmc)

library(ggmcmc)


tail(ggs(Phi_mcmc_list, family = "*"))
