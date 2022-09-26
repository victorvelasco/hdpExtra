##### Script written with the purpose of developing hdpExtra
##### Not to be built into the package.

library(hdpExtra)
library(hdp)
library(coda)
library(salso)

dim(mut_count)
## M <- mut_count
## M <- M + 1
## sum(M)
M <- matrix(1, nrow = 7, ncol = 6)
M[1,1] <- M[1,2] <- 0
M[7, ] <- c(0,1,2,4,6,10)

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

hdp_extra_chains <- vector("list", 4)
for (i in 1:4){

  # activate DPs, 10 initial components
  hdp_activated <- dp_activate(hdp_mut, 1:numdp(hdp_mut), initcc=4, seed=i*200)

  hdp_extra_chains[[i]] <- hdpExtra_posterior(hdp_activated,
                               burnin=9,
                               n=12,
                               space=1,
                               cpiter=1,
                               seed=i*1e3)
}

hdp_extra_chains <- HdpExtraChainMulti(hdp_extra_chains)

allocations_best <- salso(t(hdp_allocations(hdp_extra_chains)), loss = VI(a = 2/3))
Phi_best <- hdp_postprocessing(hdp_extra_chains, allocations_best)
hdp_plot_sig_uncertainty(Phi_best, "plots")

hdpExtra_diagnostic_plots(hdp_extra_chains, c("gamma", "alpha"))
hdp_chains <- hdp_chain(hdp_extra_chains)
hdp_chains <- hdp_multi_chain(hdp_chains)
object.size(hdp_chains)/(1024**3)

par(mfrow=c(1,2), mar=c(4, 4, 2, 1))
p1 <- lapply(chains(hdp_chains), plot_lik, bty="L", start=1000)
p2 <- lapply(chains(hdp_chains), plot_numcluster, bty="L")
p3 <- lapply(chains(hdp_chains), plot_data_assigned, bty="L")
hdp_chains <- hdp_extract_components(hdp_chains)
hdp_chains

par(mfrow=c(1,1), mar=c(5, 4, 4, 2))
plot_comp_size(hdp_chains, bty="L")

trinuc_context <- sapply(strsplit(colnames(mut_count), '\\.'), `[`, 4)
group_factor <- as.factor(rep(c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"),
                              each=16))

# pick your colours
mut_colours <- c(RColorBrewer::brewer.pal(10, 'Paired')[seq(1,10,2)], 'grey70')

plot_comp_distn(hdp_chains, cat_names=trinuc_context,
                grouping=group_factor, col=mut_colours,
                col_nonsig="grey80", show_group_labels=TRUE)


plot_dp_comp_exposure(hdp_chains, 2:11, incl_nonsig = F,
                      col_comp = 1:11)
