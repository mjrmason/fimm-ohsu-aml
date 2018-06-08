suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(synapseClient))
suppressPackageStartupMessages(library("foreach"))
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("ggplot2"))
##suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("MonoInc"))
suppressPackageStartupMessages(library("data.table"))

source("../common/dss.R")
source("../common/drc-fit.R")
source("../common/plotting.R")

num.cores <- detectCores()
if(!is.na(num.cores) && (num.cores > 1)) {
  suppressPackageStartupMessages(library("doMC"))
  cat(paste("Registering ", num.cores-1, " cores.\n", sep=""))
  registerDoMC(cores=(num.cores-1))
}

synapseLogin()

## BEGIN setup

## Read in FIMM MCM raw drug response data
synId <- "syn8488824"
obj <- synGet(id=synId, downloadFile = TRUE)
file <- getFileLocation(obj)
drug.data.long.mcm <- read.xlsx(file)

## Read in FIMM CM raw drug response data (FIMM_DSRT_Data_CM.txt)
synId <- "syn12180768"
obj <- synGet(id=synId, downloadFile = TRUE)
file <- getFileLocation(obj)
## These data have the same columns and concentration ranges (hence in nM) as the
## FIMM MCM data ("syn8488824").  Just drop the extraneous "X" column.
## And "DSRT_SET" in these CM data correspond to "DRUG_SET" in MCM data.
drug.data.long.cm <- read.table(file, sep="\t", header=TRUE, as.is=TRUE)
colnames(drug.data.long.cm)[colnames(drug.data.long.cm) == "DSRT_SET"] <- "DRUG_SET"
drug.data.long.cm <- drug.data.long.cm[, !(colnames(drug.data.long.cm) == "X")]

## library(data.table)
## FIMM CM/validation data (raw log2 batch-corrected CPM)
synId <- "syn12176271"
obj <- synGet(id=synId, downloadFile = TRUE)
file <- getFileLocation(obj)
fimm.expr.cm <- read.table(file, header=TRUE, sep=",")
## fimm.expr.cm <- fread(file, data.table = FALSE)
## Make the columns samples
fimm.expr.cm <- t(fimm.expr.cm)

## FIMM MCM data (raw log2 batch-corrected CPM; RNASeq.CPM.log2.bc.csv)
synId <- "syn8270602"
obj <- synGet(id=synId, downloadFile = TRUE)
file <- getFileLocation(obj)
fimm.expr.mcm <- read.table(file, header=TRUE, sep=",")
## Make the columns samples
fimm.expr.mcm <- t(fimm.expr.mcm)

## CM
drug.cm.samples <- unique(drug.data.long.cm[, c("SCREEN_ID", "PATIENT_ID")])
colnames(drug.cm.samples) <- c("ORIG_SCREEN_ID", "PATIENT_ID")
drug.cm.samples$SCREEN_ID <- gsub(drug.cm.samples$ORIG_SCREEN_ID, pattern="_CM", replacement="")
drug.cm.samples$media <- "CM"
drug.cm.samples$Drug.assay <- TRUE

expr.cm.samples <- colnames(fimm.expr.cm)
expr.cm.samples <- data.frame(SCREEN_ID = expr.cm.samples)
expr.cm.samples$RNAseq.assay <- TRUE
expr.cm.samples$media <- "CM"

cm.samples <- merge(expr.cm.samples, drug.cm.samples, all=TRUE)

## MCM
drug.mcm.samples <- unique(drug.data.long.mcm[, c("SCREEN_ID", "PATIENT_ID")])
colnames(drug.mcm.samples) <- c("ORIG_SCREEN_ID", "PATIENT_ID")
drug.mcm.samples$SCREEN_ID <- drug.mcm.samples$ORIG_SCREEN_ID
drug.mcm.samples$media <- "MCM"
drug.mcm.samples$Drug.assay <- TRUE

expr.mcm.samples <- colnames(fimm.expr.mcm)
expr.mcm.samples <- data.frame(SCREEN_ID = expr.mcm.samples)
expr.mcm.samples$RNAseq.assay <- TRUE
expr.mcm.samples$media <- "MCM"

mcm.samples <- merge(expr.mcm.samples, drug.mcm.samples, all=TRUE)

fimm.samples <- rbind(cm.samples, mcm.samples)

fimm.samples$INFERRED_PATIENT_ID <- gsub(fimm.samples$SCREEN_ID, pattern="^([^_]+_[^_]+)_.+$", replacement="\\1")

write.table(file = "fimm.samples.tsv", fimm.samples, sep="\t", row.names=FALSE, col.names=TRUE, quote = FALSE)

