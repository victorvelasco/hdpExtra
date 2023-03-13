setwd("/home/victor/phd/software/hdpExtra/data/")
library(dplyr)

donors_specimens <- read.table("pcawg_donors_specimens.tsv", header = TRUE, row.names = 1)
head(donors_specimens)
dim(donors_specimens)

pcawg_exposures <- read.csv("WGS_PCAWG_published_exposures.csv", header = TRUE, row.names = 1)
head(unlist(lapply(strsplit(rownames(pcawg_exposures), "::"), function(x) x[2])))

pcawg_names <- cbind(
  rownames(pcawg_exposures),
  unlist(lapply(strsplit(rownames(pcawg_exposures), "::"), function(x) x[2]))
)
pcawg_names <- as.data.frame(pcawg_names)
rownames(pcawg_names) <- pcawg_names[, 2]

donors_in_pcawg <- donors_specimens[donors_specimens$icgc_specimen_id %in% unlist(lapply(strsplit(rownames(pcawg_exposures), "::"), function(x) x[2])), ]
rownames(donors_in_pcawg) <- donors_in_pcawg[, 1]
head(pcawg_names)
rownames(pcawg_names) <- pcawg_names[, 2]
pcawg_names$donor <- donors_in_pcawg[rownames(pcawg_names)]

nikzainal_exposures <- read.table("icgc_nikzainal_published_exposures.csv", header = TRUE)
rownames(nikzainal_exposures) <- nikzainal_exposures[, 1]
nikzainal_exposures <- nikzainal_exposures[-c(1,2,3)]
nikzainal_exposures[1:5,1:5]

nikzainal_exposures[, "donor_id"] <- donors_specimens[rownames(nikzainal_exposures), 1]
nikzainal_exposures <- nikzainal_exposures[!is.na(nikzainal_exposures$donor_id), ]
rownames(nikzainal_exposures) <- nikzainal_exposures[, "donor_id"]
nikzainal_exposures[, "donor_id_full"] <- pcawg_names[nikzainal_exposures$donor_id, 1]
rownames(nikzainal_exposures) <- nikzainal_exposures[, "donor_id_full"]
nikzainal_exposures <- nikzainal_exposures[, -c(40,41)]
head(nikzainal_exposures)
write.csv(nikzainal_exposures, "icgc_nikzainal_published_exposures.csv")
