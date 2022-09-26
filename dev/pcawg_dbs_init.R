library(dplyr)
library(hdpExtra)
library(hdp)

M <- read.csv("data/WGS_PCAWG.dinucs.csv", header = TRUE, row.names = 1) %>%
  as.matrix %>%
  t %>%
  as.data.frame

M <- M[rowSums(M) > 0, ]

tissues <- sapply(1:nrow(M), function(i) unlist(strsplit(rownames(M)[i], "[.]"))[[1]])

npatients_per_tissue <- unname(table(tissues))
ntissues <- length(npatients_per_tissue)

hdp_mut <- hdp_init(ppindex = c(0, rep(1, ntissues), rep(1+(1:ntissues), npatients_per_tissue)), # index of parental nodes
                    cpindex = c(1, rep(2, ntissues), rep(3, nrow(M))), # index of concentration param
                    hh=rep(1, ncol(M)), # prior is uniform over the 96 mutation categories
                    alphaa=rep(1, 3), # shape hyperparams for five different CPs
                    alphab=rep(1, 3)) # rate hyperparams for five different CPs
hdp_mut <- hdp_setdata(hdp_mut,
                       dpindex = (2+ntissues):numdp(hdp_mut), # index of nodes to add data to
                       M) # input data (mutation counts, sample rows match up with specified dpindex)

saveRDS(hdp_mut, "output/pcawg_dbs/hdp_mut.rds")
