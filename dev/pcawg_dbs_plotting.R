library(hdp)
library(salso)
library(dplyr)

M <- read.csv("data/WGS_PCAWG.dinucs.csv", header = TRUE, row.names = 1) %>%
  as.matrix %>%
  t %>%
  as.data.frame
M <- M[rowSums(M) > 0, ]

chains <- list()
for (i in 1:4) {
  chains[[i]] <- readRDS(paste0("output/pcawg_dbs/hdpExtra_chain", i, ".rds"))
}
chains <- HdpExtraChainMulti(chains)

### pdf(file = "dev/plots/diagnostics.pdf")
### hdpExtra_diagnostic_plots(chains, c("gamma", "alpha0", "alpha1"))
### dev.off()

salso_obj <- readRDS("output/pcawg_dbs/salso_VI1.rds")
sort(table(salso_obj))
Phi_best <- hdp_postprocessing(chains, salso_obj)
hdp_plot_sig_dbs_uncertainty(Phi_best, "plots")

J <- nrow(M)
K <- max(salso_obj)

nmuts <- rowSums(M)
index_end <- cumsum(nmuts)
index_beg <- c(0, index_end[1:(length(nmuts)-1)])+1
names(index_beg) <- names(index_end)

exposures <- matrix(0, nrow = J, ncol = K)
rownames(exposures) <- rownames(M)
colnames(exposures) <- paste("Sig", 01:K)

for (j in 1:J) {
  patient_allocs <- table(salso_obj[index_beg[j]:index_end[j]])
  exposures[j, as.integer(names(patient_allocs))] <- unname(patient_allocs)
}

weights <- diag(1/rowSums(exposures)) %*% exposures
colnames(weights) <- colnames(exposures)
rownames(weights) <- rownames(exposures)
tissues <- sapply(1:nrow(M), function(i) unlist(strsplit(rownames(M)[i], "[.]"))[[1]])

exposures_alexandrov <- read.csv("data/WGS_PCAWG_published_exposures.csv", header = TRUE, row.names = 1)
plot(
  exposures_alexandrov[, "DBS7"],
  exposures[, "Sig 9"]
)
abline(coef = c(0, 1))

exposures_nikzainal <- read.csv("data/icgc_nikzainal_published_exposures.csv", header = TRUE, row.names = 1)
plot(
  exposures_nikzainal[, "DBS26"],
  exposures[rownames(exposures_nikzainal), "Sig 9"],
  ylim = c(0, 350),
  ylab = "Exposure to Victor's Sig 9",
  xlab = "Exposure to Nik-Zainal's Sig DBS26 "
)
abline(coef = c(0, 1))



exposures_nikzainal <- read.csv("data/icgc_nikzainal_published_exposures.csv", header = TRUE, row.names = 1)
plot(
  exposures_nikzainal[, "DBS25"],
  exposures[rownames(exposures_nikzainal), "Sig 12"],
  cex = sample(6)/2,
  pch = 19
)
abline(coef = c(0, 1))

#### npatients_per_tissue <- unname(table(tissues))
#### ntissues <- length(npatients_per_tissue)
####
#### exposures[, "tissue"] <- tissues
#### exposures <- as.data.frame(exposures)
#### exposures %>%
####   group_by(tissue) %>%
####   summarise(across(everything(), sum))
####
#### knitr::kable(exposures %>%
####                group_by(tissue) %>%
####                summarise(across(everything(), sum)), "pipe")
#### exposures_by_tissue <- exposures %>%
####   group_by(tissue) %>%
####   summarise(across(everything(), sum)) %>%
####   select(where(~ !is.numeric(.x) || sum(.x) > 100))
####
#### library(ggplot2)
#### library(tidyr)
#### p1 <- exposures_by_tissue %>%
####   pivot_longer(., -tissue, names_to = "Signature", values_to = "Value") %>%
####   ggplot(aes(x = tissue, y = Value, fill = Signature)) +
####   geom_bar(stat = "identity", position = "fill") +
####   xlab("Tissue") + ylab("Proportion of mutations") +
####   theme(axis.text.x = element_text(angle = 90)) + ggtitle("Attributions of mutations to signatures") +
####   theme(text = element_text(size=12)) +
####   scale_fill_manual(values = pals::glasbey())
#### p1
