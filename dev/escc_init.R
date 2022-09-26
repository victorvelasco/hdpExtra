library(dplyr)

M <- read.table("data/ESCC_SBS96.txt", header = TRUE, row.names = 1) %>%
  as.matrix %>%
  t %>%
  as.data.frame

patient_country <- read.table("data/ESCC_patient_country.txt", row.names = 1, header = TRUE)
# patient_country <- patient_country[order(patient_country$country), ]
npatients_by_country <- patient_country %>%
  group_by(country) %>%
  summarise(n = n())
categories <- colnames(M)

M <- M[, order(substr(categories, 3, 5))]
M <- cbind(patient_country, M)
M <- M[order(M$country), ]
M <- M[-c(1)]
countries <- patient_country$country %>% unique
ncountries <- length(countries)
npatients <- nrow(M)


hdp_mut <- hdp_init(ppindex = c(0, rep(1, ncountries), rep(1+(1:ncountries), npatients_by_country$n)), # index of parental nodes
                    cpindex = c(1, rep(2, ncountries), rep(3, nrow(M))), # index of concentration param
                    hh=rep(1, 96), # prior is uniform over the 96 mutation categories
                    alphaa=rep(10, 3), # shape hyperparams for five different CPs
                    alphab=rep(1, 3)) # rate hyperparams for five different CPs
hdp_mut <- hdp_setdata(hdp_mut,
                       dpindex = 10:numdp(hdp_mut), # index of nodes to add data to
                       M) # input data (mutation counts, sample rows match up with specified dpindex)

saveRDS(hdp_mut, "output/escc/hdp_mut.rds")
