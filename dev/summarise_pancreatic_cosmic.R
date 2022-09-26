library(hdp)
library(coda)
library(salso)
library(dplyr)

salso_path <- "output/glioblastoma/salso_VI25.rds"

chains_path <- "output/glioblastoma/hdpExtra_chain"

hdp_extra_chains <- vector("list", 4)
for (i in 1:4){
  hdp_extra_chains[[i]] <- readRDS(paste0(chains_path, i, ".rds"))
}

hdp_extra_chains <- HdpExtraChainMulti(hdp_extra_chains)

allocations_best <- readRDS(salso_path)
rm(list = (attr(allocations_best, "draws")))
Phi_best <- hdp_postprocessing(hdp_extra_chains, allocations_best)
hdp_plot_sig_uncertainty(Phi_best, "plots")

hdpExtra_diagnostic_plots(hdp_extra_chains, c("gamma", "alpha"))

acf(hdp_extra_chains@cp_values[[1]][10001:20000,2], lag.max = 100)

cpvalues <- hdp_extra_chains@cp_values
for (i in 1:4) {
  cpvalues[[i]] <- as.mcmc(cpvalues[[i]][-c(1:10000),])
}
cpvalues <- as.mcmc.list(cpvalues)
coda::gelman.diag(cpvalues)

my_acf <- acf(hdp_extra_chains@cp_values[[1]][10001:20000,1])
acl <- 2*sum(my_acf[[1]])-1
10000/acl

mut_example_multi <- hdp_extract_components(hdp_multi_chain(hdp_chain(hdp_extra_chains)))
par(mfrow=c(2,2), mar=c(4, 4, 2, 1))
p1 <- lapply(chains(mut_example_multi), plot_lik, bty="L", start=1000)
p2 <- lapply(chains(mut_example_multi), plot_numcluster, bty="L")
p3 <- lapply(chains(mut_example_multi), plot_data_assigned, bty="L")

mut_example_multi <- hdp_extract_components(mut_example_multi)
mut_example_multi

par(mfrow=c(1,1), mar=c(5, 4, 4, 2))
plot_comp_size(mut_example_multi, bty="L")

trinuc_context <- sapply(strsplit(colnames(mut_count), '\\.'), `[`, 4)
group_factor <- as.factor(rep(c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"),
                              each=16))

# pick your colours
mut_colours <- c(RColorBrewer::brewer.pal(10, 'Paired')[seq(1,10,2)], 'grey70')

plot_comp_distn(mut_example_multi, cat_names=trinuc_context,
                grouping=group_factor, col=mut_colours,
                col_nonsig="grey80", show_group_labels=TRUE)
plot_dp_comp_exposure(mut_example_multi, main_text="a",
                      dpindices=1+(1:41),
                      col=RColorBrewer::brewer.pal(12, "Set3"),
                      incl_nonsig=FALSE,
                      ylab_numdata = 'SNV count', ylab_exp = 'Signature exposure',
                      leg.title = 'Signature')
