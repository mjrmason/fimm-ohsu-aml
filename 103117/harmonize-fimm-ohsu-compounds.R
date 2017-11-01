suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(synapseClient))

synapseLogin()

cat("Harmonize FIMM and OHSU compounds:\n")
cat("FIMM and OHSU will be harmonized by manually created map.\n")

## Load the map of drugs shared between OHSU and FIMM (FIMM_OHSU_Drug_Concentrations.xlsx)
synId <- "syn10163669"
obj <- synGet(synId, downloadFile = TRUE)
ohsu.fimm.drugs <- openxlsx:::read.xlsx(getFileLocation(obj), sheet=1)
ohsu.fimm.drugs <- ohsu.fimm.drugs[, !(colnames(ohsu.fimm.drugs) == "X1")]

## Save the mapping of FIMM and OHSU compounds in the folder
## FIMM_BEAT AML Collaboration/Files/Analyses/fimm-ohsu (syn11361089)
parentId <- "syn11361089"
path <- paste0("fimm-ohsu-drug-map.tsv")
write.table(file=path, ohsu.fimm.drugs, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
f <- File(path, parentId = parentId, synapseStore = TRUE)
synStore(f, executed = NULL, forceVersion = FALSE)

cat("Successfully created fimm-ohsu drug map.\n")