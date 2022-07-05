library(hdp)
library(hdpExtra)
library(dplyr)

M <- read.table("data/ESCC_SBS96.txt", header = TRUE, row.names = 1) %>%
  as.matrix %>%
  t %>%
  as.data.frame

patient_country <- read.table("data/ESCC_patient_country.txt", header = TRUE)
patient_country <- patient_country[order(patient_country$country), ]
npatients_by_country <- patient_country %>%
  group_by(country) %>%
  summarise(n = n())

categories <- colnames(M)

M <- M[, order(substr(categories, 3, 5))]

countries <- patient_country$country %>% unique
ncountries <- length(countries)
npatients <- nrow(M)


##### hdp_mut <- hdp_init(ppindex = c(0, rep(1, ncountries), rep(1+(1:ncountries), npatients_by_country$n)), # index of parental nodes
#####                     cpindex = c(1, rep(2, ncountries), rep(2+(1:ncountries), npatients_by_country$n)), # index of concentration param
#####                     hh=rep(1, 96), # prior is uniform over the 96 mutation categories
#####                     alphaa=rep(1, 10), # shape hyperparams for five different CPs
#####                     alphab=rep(1, 10)) # rate hyperparams for five different CPs

hdp_mut <- hdp_init(ppindex = c(0, rep(1, ncountries), rep(2, npatients)), # index of parental nodes
                    cpindex = c(1, rep(2, ncountries), rep(3, npatients)), # index of concentration param
                    hh=rep(1, 96), # prior is uniform over the 96 mutation categories
                    alphaa=rep(1, 3), # shape hyperparams for five different CPs
                    alphab=rep(1, 3)) # rate hyperparams for five different CPs

hdp_mut <- hdp_setdata(hdp_mut,
                       dpindex = 10:numdp(hdp_mut), # index of nodes to add data to
                       M) # input data (mutation counts, sample rows match up with specified dpindex)

hdp_extra_chains <- vector("list", 4)
for (i in 1:4){

  # activate DPs, 10 initial components
  hdp_activated <- dp_activate(hdp_mut, 1:numdp(hdp_mut), initcc=50, seed=i*200)

  hdp_extra_chains[[i]] <- hdpExtra_posterior(hdp_activated,
                                              burnin=100,
                                              n=100,
                                              space=1,
                                              cpiter=1,
                                              seed=i*1e3)
  saveRDS(hdp_extra_chains[[i]], paste0("output/chain", i, "_escc_mutographs.rds"))
}

hdp_extra_chains <- HdpExtraChainMulti(hdp_extra_chains)
saveRDS(hdp_extra_chains, "output/chains_escc_mutographs.rds")
