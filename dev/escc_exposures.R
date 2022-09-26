library(hdp)
library(salso)
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

nmuts_by_patient <- rowSums(M)

M <- M[, order(substr(categories, 3, 5))]
M <- cbind(patient_country, M)
M <- M[order(M$country), ]
M <- M[-c(1)]
countries <- patient_country$country %>% unique
ncountries <- length(countries)
npatients <- nrow(M)

salso_obj <- readRDS("output/escc/salso_binder1.rds")
salso_obj <- unclass(salso_obj)
object.size(salso_obj)/(1024**3)
length(salso_obj)
attr(salso_obj, "draws") <- NULL
table(salso_obj)

i <- 1 # Iterates over mutations in the dataset
J <- nrow(M) # Number of patients
exposures_matrix <- matrix(0, nrow = J, ncol = max(salso_obj))
for(j in 1:J) {
  exposures_patient <- table(salso_obj[i:(i+nmuts_by_patient[j]-1)])
  exposures_matrix[j, as.integer(names(exposures_patient))] <- exposures_patient
  i <- i + unname(nmuts_by_patient[j])
}
rownames(exposures_matrix) <- rownames(M)

opium <- c("PD39988a", "PD40018a", "PD40335a", "PD42843a", "PD40854a", "PD40319a", "PD43793a", "PD39994a", "PD40021a", "PD39957a", "PD39972a", "PD40320a", "PD39987a", "PD40356a", "PD39974a", "PD40355a", "PD39978a", "PD39999a", "PD42845a", "PD39983a")
exposures_matrix[opium, ]

no_opium <- c("PD39493a", "PD39523a", "PD39955a", "PD39958a", "PD39959a", "PD39960a", "PD39961a", "PD39962a", "PD39963a", "PD39964a", "PD39965a", "PD39966a", "PD39967a", "PD39969a", "PD39970a", "PD39971a", "PD39973a", "PD39975a", "PD39976a", "PD39977a", "PD39979a", "PD39981a", "PD39984a", "PD39985a", "PD39989a", "PD39990a", "PD39991a", "PD39992a", "PD39993a", "PD39996a", "PD39998a", "PD40000a", "PD40001a", "PD40004a", "PD40006a", "PD40007a", "PD40009a", "PD40010a", "PD40011a", "PD40012a", "PD40013a", "PD40015a", "PD40017a", "PD40019a", "PD40020a", "PD40023a", "PD40024a", "PD40025a", "PD40316a", "PD40317a", "PD40318a", "PD40321a", "PD40322a", "PD40323a", "PD40324a", "PD40325a", "PD40326a", "PD40327a", "PD40328a", "PD40329a", "PD40330a", "PD40331a", "PD40332a", "PD40333a", "PD40334a", "PD40336a", "PD40337a", "PD40338a", "PD40339a", "PD40340a", "PD40343a", "PD40345a", "PD40346a", "PD40348a", "PD40349a", "PD40350a", "PD40351a", "PD40352a", "PD40353a", "PD40354a", "PD40357a", "PD40358a", "PD40360a", "PD40847a", "PD40848a", "PD40849a", "PD40850a", "PD40852a", "PD40853a", "PD40855a", "PD40856a", "PD40857a", "PD40858a", "PD40859a", "PD40861a", "PD41466a", "PD41467a", "PD41468a", "PD41470a", "PD41471a", "PD41472a", "PD41473a", "PD41474a", "PD41475a", "PD41476a", "PD41477a", "PD41478a", "PD41479a", "PD41481a", "PD42838a", "PD42839a", "PD42840a", "PD42841a", "PD42842a", "PD42844a", "PD42846a", "PD42847a", "PD42848a", "PD42849a", "PD42850a", "PD42852a", "PD42853a", "PD42854a", "PD43304a", "PD43305a", "PD43306a", "PD43777a", "PD43779a", "PD43780a", "PD43782a", "PD43783a", "PD43785a", "PD43786a", "PD43787a", "PD43789a", "PD43790a", "PD43791a", "PD43792a", "PD43795a", "PD43796a", "PD43798a", "PD43800a", "PD43802a", "PD43803a", "PD43806a", "PD43807a", "PD43809a", "PD43810a", "PD43811a", "PD43812a", "PD44659a", "PD44660a", "PD44661a", "PD44662a", "PD44663a", "PD44690a", "PD45482a", "PD45483a")
exposures_matrix[no_opium, ]

