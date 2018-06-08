suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(synapseClient))

synapseLogin()

## BEGIN setup

## Read in FIMM CM raw drug response data (FIMM_DSRT_Data_CM.txt)
drug.synId <- "syn12180768"
## FIMM CM/validation data (raw log2 batch-corrected CPM)
expr.synId <- "syn12176271"
## Where to store subsetted drug and expr data
## FIMM Data/Processed Data
expr.parent.synId <- "syn8270577"
## FIMM Data/Processed Data
drug.parent.synId <- "syn8270577"

expr.prefix <- "fimm-cm-expr-data-with-drug-data"
drug.prefix <- "fimm-cm-drug-data-with-expr-data"

## END setup

obj <- synGet(id=drug.synId, downloadFile = TRUE)
file <- getFileLocation(obj)
## These data have the same columns and concentration ranges (hence in nM) as the
## FIMM MCM data ("syn8488824").  Just drop the extraneous "X" column.
## And "DSRT_SET" in these CM data correspond to "DRUG_SET" in MCM data.
drug.data.long <- read.table(file, sep="\t", header=TRUE, as.is=TRUE)
colnames(drug.data.long)[colnames(drug.data.long) == "DSRT_SET"] <- "DRUG_SET"
drug.data.long <- drug.data.long[, !(colnames(drug.data.long) == "X")]

obj <- synGet(id=expr.synId, downloadFile = TRUE)
file <- getFileLocation(obj)
fimm.expr <- read.table(file, header=TRUE, sep=",")
## Make the columns samples
fimm.expr <- t(fimm.expr)

drug.samples <- unique(drug.data.long[, c("SCREEN_ID", "PATIENT_ID")])
colnames(drug.samples) <- c("ORIG_SCREEN_ID", "PATIENT_ID")
drug.samples$SCREEN_ID <- gsub(drug.samples$ORIG_SCREEN_ID, pattern="_CM", replacement="")
drug.samples$Drug.assay <- TRUE

expr.samples <- colnames(fimm.expr)
expr.samples <- data.frame(SCREEN_ID = expr.samples)
expr.samples$RNAseq.assay <- TRUE

samples <- merge(expr.samples, drug.samples, all=TRUE)



