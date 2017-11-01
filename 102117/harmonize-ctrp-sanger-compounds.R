suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(synapseClient))
suppressPackageStartupMessages(library("foreach"))
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library("plyr"))

num.cores <- detectCores()
if(!is.na(num.cores) && (num.cores > 1)) {
  suppressPackageStartupMessages(library("doMC"))
  cat(paste("Registering ", num.cores-1, " cores.\n", sep=""))
  registerDoMC(cores=(num.cores-1))
}

synapseLogin()

cat("Harmonize CTRP and Sanger compounds:\n")
cat("CTRP and Sanger will be harmonized using PharmacoGx.\n")

## Use PharmacoGx to create a map between CTRP/CCLE and GDSC/Sanger
## No--instead we will use the published map between CTRP/CCLE and GSDC/Sanger

## Load the CTRP drugs
## CTRP folder on Synapse: syn5622707
## CTRP drug info: "v20.meta.per_compound.txt" (syn5632193)
synId <- "syn5632193"
obj <- synGet(synId, downloadFile = TRUE)
ctrp.compounds <- read.table(getFileLocation(obj), sep="\t", header=TRUE, as.is=TRUE, comment.char="", quote="")
ctrp.compounds <- ctrp.compounds[, c("cpd_name", "cpd_smiles", "master_cpd_id")]

## CTRPv2 pset in PharmacoGx has all of the ctrp compounds
##library(PharmacoGx)
##ctrp.pset <- downloadPSet("CTRPv2")
##table(ctrp.compounds$master_cpd_id %in% drugInfo(ctrp.pset)$master_cpd_id)

## Load the GDSC/Sanger drugs
## GDSC drug inf: "Screened_Compounds.xlsx"  (syn11275819)
synId <- "syn11275819"
obj <- synGet(synId, downloadFile = TRUE)
sanger.compounds <- openxlsx:::read.xlsx(getFileLocation(obj), sheet=1)

## But it's not clear what GDSC data set in PharmacoGx is the right one
## gdsc1000.pset <- downloadPSet("GDSC1000")
## table(sanger.compounds$Drug.ID %in% drugInfo(gdsc1000.pset)$DRUG.ID)

## Load the published map between CTRPv2 and GDSC
synId <- "syn11289136"
obj <- synGet(synId, downloadFile = TRUE)
ctrp.to.sanger.map <- openxlsx:::read.xlsx(getFileLocation(obj), sheet=1, startRow = 11)

## Make sure everything in the map is the CTRP and Sanger data sets
flag <- ctrp.to.sanger.map$GDSC1000.ids %in% sanger.compounds$Drug.ID
if(all(flag)) {
  cat("As expected all drugs in CTRP to Sanger map are in Sanger data\n")
} else {
  stop(paste0("Unexpectedly, the following drugs in CTRP to Sanger map are not in Sanger data: ",
       paste(ctrp.to.sanger.map$GDSC1000.ids[!flag], sep=","), "\n"))
}

flag <- ctrp.to.sanger.map$CTRP.master.cpd.id %in% ctrp.compounds$master_cpd_id
if(all(flag)) {
  cat("As expected all drugs in CTRP to Sanger map are in CTRP data\n")
} else {
  stop(paste0("Unexpectedly, the following drugs in CTRP to Sanger map are not in CTRP data: ",
       paste(ctrp.to.sanger.map$CTRP.master.cpd.id[!flag], sep=","), "\n"))
}

## The column loaded as "X2" has additional GDSC/Sanger ids for a drug.  i.e., is this entry
## is not NA, a drug has multiple mappings to CTRP in that row.  Drop such ambiguous drugs.
ctrp.to.sanger.map <- ctrp.to.sanger.map[is.na(ctrp.to.sanger.map$X2), c("CTRP.master.cpd.id", "GDSC1000.ids")]
colnames(ctrp.to.sanger.map) <- c("master_cpd_id", "Drug.ID")

ctrp.to.sanger.map <- merge(ctrp.to.sanger.map, ctrp.compounds)
ctrp.to.sanger.map <- merge(ctrp.to.sanger.map, sanger.compounds)
ctrp.to.sanger.map <- ctrp.to.sanger.map[!is.na(ctrp.to.sanger.map$cpd_smiles),]

## Save the mapping of CTRP to Sanger compounds in the folder
## FIMM_BEAT AML Collaboration/Files/Analyses/ctrp-sanger (syn11297344)
parentId <- "syn11297344"
path <- paste0("ctrp-sanger-drug-map.tsv")
write.table(file=path, ctrp.to.sanger.map, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
f <- File(path, parentId = parentId, synapseStore = TRUE)
synStore(f, executed = NULL, forceVersion = FALSE)

cat("Successfully created ctrp-sanger drug map.\n")