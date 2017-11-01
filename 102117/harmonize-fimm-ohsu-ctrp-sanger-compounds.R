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

cat("Harmonize FIMM, OHSU, CTRP, and Sanger compounds:\n")
cat("CTRP and Sanger will be harmonized using PharmacoGx.\n")
cat("FIMM and CTRP will be harmonized by converted FIMM InChi to SMILES (via webchem) and then comparing SMILES.\n")
cat("FIMM and OHSU will be harmonized by manually created map.\n")

## Robert Allaway's code to compare SMILES

## To run rcdk, needed to update Java JDK as described here:
## http://www.webupd8.org/2012/09/install-oracle-java-8-in-ubuntu-via-ppa.html
## Installed rJava from source as described here:
## install.packages("rJava", type = "source")
## Installed rcdk as described here:
## https://github.com/rajarshi/cdkr
## git clone https://github.com/rajarshi/cdkr.git
##   R CMD build rcdklibs
##   R CMD INSTALL rcdklibs_*gz
##   cd rcdkjar
##   ant clean jar
##   cd ../
##   R CMD build rcdk
##   R CMD INSTALL rcdk_*gz
library(dplyr)
library(rJava)
library(rcdk)
library(fingerprint)
library(parallel)

## input SMILES to fingerprint function
parseInputFingerprint <- function(input) {
  input.mol <- parse.smiles(input)
  fp.inp <- lapply(input.mol, get.fingerprint, type = "extended")
}

parseSingleInputFingerprint <- function(input) {
  input.mol <- NULL
  tryCatch({
    input.mol <- parse.smiles(input)
  }, error = function(e) { return(NULL) })
  fp.inp <- NULL
  tryCatch({
    fp.inp <- get.fingerprint(input.mol[[1]], type = "extended")
  }, error = function(e) { return(NULL) })
  fp.inp
}

## Load the map of drugs shared between OHSU and FIMM (FIMM_OHSU_Drug_Concentrations.xlsx)
synId <- "syn10163669"
obj <- synGet(synId, downloadFile = TRUE)
ohsu.fimm.drugs <- openxlsx:::read.xlsx(getFileLocation(obj), sheet=1)
ohsu.fimm.drugs <- ohsu.fimm.drugs[, !(colnames(ohsu.fimm.drugs) == "X1")]

## Load the FIMM annotations that have the InChi
synId <- "syn8434109"
obj <- synGet(synId, downloadFile = TRUE)
fimm.anno <- openxlsx:::read.xlsx(getFileLocation(obj), sheet=1)

## Merge the two tables
orig.sz <- nrow(ohsu.fimm.drugs)
ohsu.fimm.drugs <- merge(ohsu.fimm.drugs, fimm.anno, by.x = "FIMM_Batch_ID_Drug", by.y = "ID_Drug", all.x = TRUE)

## Convert the InChI to smiles using webchem
library(webchem)

## Use the chemical identifier resolver (cir) to resolve InChIKeys to SMILES
ohsu.fimm.drugs <- ohsu.fimm.drugs[!is.na(ohsu.fimm.drugs$InChI), ]
val <- cir_query(ohsu.fimm.drugs$InChI, representation = "smiles")
inchi.to.smiles <- ldply(names(val),
                         .fun = function(key) { 
                           df <- data.frame(smiles = val[[key]])
                           df$InChI <- key
                           df
                         })
inchi.to.smiles <- inchi.to.smiles[!is.na(inchi.to.smiles$InChI),]
inchi.to.smiles <- inchi.to.smiles[!is.na(inchi.to.smiles$smiles),]

## Merge the smiles to the OHSU/FIMM map
## Drop the dummy smiles column
ohsu.fimm.drugs <- ohsu.fimm.drugs[, !(grepl(colnames(ohsu.fimm.drugs), pattern="smiles", ignore.case=TRUE))]
## Drop some other extraneous columns
ohsu.fimm.drugs <- ohsu.fimm.drugs[, !(grepl(colnames(ohsu.fimm.drugs), pattern="conc", ignore.case=TRUE))]
ohsu.fimm.drugs <- ohsu.fimm.drugs[, !(grepl(colnames(ohsu.fimm.drugs), pattern="\\.y", ignore.case=TRUE))]
ohsu.fimm.drugs <- ohsu.fimm.drugs[, !(grepl(colnames(ohsu.fimm.drugs), pattern="\\.1", ignore.case=TRUE))]
colnames(ohsu.fimm.drugs) <- gsub(colnames(ohsu.fimm.drugs), pattern="\\.x", replacement="")
ohsu.fimm.drugs <- merge(ohsu.fimm.drugs, inchi.to.smiles, by.x = "InChI", by.y = "InChI")
ohsu.fimm.drugs <- ohsu.fimm.drugs[!is.na(ohsu.fimm.drugs$smiles),]

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

cat("Calculating drug fingerprints for CTRP/Sanger\n")
fp.ctrp <- parseInputFingerprint(as.character(ctrp.to.sanger.map$cpd_smiles))

cat("Calculating drug fingerprints for OHSU/FIMM\n")
fp.ohsu <- parseInputFingerprint(as.character(ohsu.fimm.drugs$smiles))

## compares every smiles on list one across list 2
## This calculates Tanimoto distance (the default)
cat("Calculating all pairwise Tanimoto distances between CTRP/Sanger drug fingerprints and OHSU/FIMM drug fingerprints\n")
sims <- mclapply(fp.ctrp, function(i) {
  sim <- sapply(fp.ohsu, function(j) {
    distance(i, j)
  })
  bar <- data.frame(sim=sim, match=names(sim))
  bar
}, mc.cores = detectCores())

## This is the Tanimoto cutoff used by the PharmacoGx paper
cutoff <- 0.95
sims.cutoff <- llply(sims, .parallel = FALSE, 
                     .fun = function(entry) {
                       flag <- entry$sim > cutoff
                       if(any(flag)) { return(entry$match[flag]) }
                       return(NULL)
                     })
flag <- unlist(lapply(sims.cutoff, is.null))
sims.cutoff <- sims.cutoff[!flag]
flag <- unlist(lapply(sims.cutoff, function(lst) length(lst) != 1))
sims.cutoff <- sims.cutoff[!flag]

struct.mapping <- ldply(1:length(sims.cutoff),
                           .fun = function(i) {
                             df <- data.frame(smiles2 = unlist(sims.cutoff[[i]]))
                             df$smiles1 <- names(sims.cutoff)[i]
                             df
                           })
struct.mapping <- struct.mapping[, c("smiles1", "smiles2")]
colnames(struct.mapping) <- c("ctrp.smiles", "ohsu.smiles")

flag <- struct.mapping$ohsu.smiles %in% ohsu.fimm.drugs$smiles
if(all(flag)) {
  cat("As expected, all mapped smiles are in the OHSU/FIMM data\n")
} else {
  stop(paste0("Unexpectedly, the following mapped smiles were not in the OHSU/FIMM data: ",
       paste(struct.mapping$ohsu.smiles[!flag], collapse=","), "\n"))
}

flag <- struct.mapping$ctrp.smiles %in% ctrp.to.sanger.map$cpd_smiles
if(all(flag)) {
  cat("As expected, all mapped smiles are in the CTRP/Sanger data\n")
} else {
  stop(paste0("Unexpectedly, the following mapped smiles were not in the CTRP/Sanger data: ",
       paste(struct.mapping$ctrp.smiles[!flag], collapse=","), "\n"))
}

grand.map <- merge(ctrp.to.sanger.map, struct.mapping, by.x = "cpd_smiles", by.y = "ctrp.smiles", all = FALSE)
grand.map <- merge(grand.map, ohsu.fimm.drugs, by.x = "ohsu.smiles", by.y = "smiles", all = FALSE)

## Save the mapping of CTRP to NTAP compounds in the folder
## FIMM_BEAT AML Collaboration/Files/Analyses/ctrp-sanger-fimm-ohsu (syn11288886)
parentId <- "syn11288886"
path <- paste0("ctrp-sanger-fimm-ohsu-drug-map.tsv")
write.table(file=path, grand.map, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
f <- File(path, parentId = parentId, synapseStore = TRUE)
synStore(f, executed = NULL, forceVersion = FALSE)

cat("Successfully created ctrp-sanger-fimm-ohsu drug map.\n")