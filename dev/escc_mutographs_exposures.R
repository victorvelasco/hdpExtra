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

salso_obj <- readRDS("output/escc/salso_VI23.rds")
tail(salso_obj %>% table %>% sort, 20)

salso_obj <- unclass(salso_obj)

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

have_opium_sig <- c("PD39988a", 
"PD40018a", 
"PD40335a", 
"PD42843a", 
"PD40854a", 
"PD40319a", 
"PD43793a", 
"PD39994a", 
"PD40021a", 
"PD39957a", 
"PD39972a", 
"PD40320a", 
"PD39987a", 
"PD40356a", 
"PD39974a", 
"PD40355a", 
"PD39978a", 
"PD39999a", 
"PD42845a", 
"PD39983a")

not_have_opium_sig <- c(
"PD45483a",
"PD45482a",
"PD44690a",
"PD44663a",
"PD44662a",
"PD44661a",
"PD44660a",
"PD44659a",
"PD43812a",
"PD43811a",
"PD43810a",
"PD43809a",
"PD43807a",
"PD43806a",
"PD43803a",
"PD43802a",
"PD43800a",
"PD43798a",
"PD43796a",
"PD43795a",
"PD43792a",
"PD43791a",
"PD43790a",
"PD43789a",
"PD43787a",
"PD43786a",
"PD43785a",
"PD43783a",
"PD43782a",
"PD43780a",
"PD43779a",
"PD43777a",
"PD43306a",
"PD43305a",
"PD43304a",
"PD42854a",
"PD42853a",
"PD42852a",
"PD42850a",
"PD42849a",
"PD42848a",
"PD42847a",
"PD42846a",
"PD42844a",
"PD42842a",
"PD42841a",
"PD42840a",
"PD42839a",
"PD42838a",
"PD41481a",
"PD41479a",
"PD41478a",
"PD41477a",
"PD41476a",
"PD41475a",
"PD41474a",
"PD41473a",
"PD41472a",
"PD41471a",
"PD41470a",
"PD41468a",
"PD41467a",
"PD41466a",
"PD40861a",
"PD40859a",
"PD40858a",
"PD40857a",
"PD40856a",
"PD40855a",
"PD40853a",
"PD40852a",
"PD40850a",
"PD40849a",
"PD40848a",
"PD40847a",
"PD40360a",
"PD40358a",
"PD40357a",
"PD40354a",
"PD40353a",
"PD40352a",
"PD40351a",
"PD40350a",
"PD40349a",
"PD40348a",
"PD40346a",
"PD40345a",
"PD40343a",
"PD40340a",
"PD40339a",
"PD40338a",
"PD40337a",
"PD40336a",
"PD40334a",
"PD40333a",
"PD40332a",
"PD40331a",
"PD40330a",
"PD40329a",
"PD40328a",
"PD40327a",
"PD40326a",
"PD40325a",
"PD40324a",
"PD40323a",
"PD40322a",
"PD40321a",
"PD40318a",
"PD40317a",
"PD40316a",
"PD40025a",
"PD40024a",
"PD40023a",
"PD40020a",
"PD40019a",
"PD40017a",
"PD40015a",
"PD40013a",
"PD40012a")

knitr::kable(exposures[have_opium_sig, as.integer(rev(names(tail(salso_obj %>% table %>% sort, 20))))], "pipe")

knitr::kable(exposures[not_have_opium_sig, as.integer(rev(names(tail(salso_obj %>% table %>% sort, 20))))], "pipe")

#### knitr::kable(cbind(patient_country$country, exposures[, as.integer(rev(names(tail(salso_obj %>% table %>% sort, 20))))]), "pipe")
