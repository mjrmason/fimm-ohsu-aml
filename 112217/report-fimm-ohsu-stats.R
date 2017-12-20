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

cat("Use DRC to fit L.4 (4-parameter logistic) and LL.4 (4-paramter log logistic) functions to FIMM data.\n")

## BEGIN setup

data.set <- "FIMM"
prefix <- "fimm"

## Read in FIMM raw drug response data
synId <- "syn8488824"
obj <- synGet(id=synId, downloadFile = TRUE)
file <- getFileLocation(obj)
drug.data.long <- read.xlsx(file)

## Load FIMM metadata (including relapse/refractory/diagnosis)
synId <- "syn8270594"
obj <- synGet(id=synId, downloadFile = TRUE, downloadLocation = ".")
file <- getFileLocation(obj)
fimm.metadata <- read.table(file, header=TRUE, sep=",")

## drug.data.long has one concentration/response per row.
## NB: concentrations are in nM.  Let's make that explicit.
drug.data.long$conc.nM <- drug.data.long$CONCENTRATION

## NB: PERCENT_INHBITION is ... inhibition!
## NB: PERCENT_INHIBITION really is a percent -- it runs from 0 to 100 (well, actually, slightly below and above).

conc.col <- "conc.nM"
response.col <- "PERCENT_INHIBITION"
response.is.viability <- FALSE

drug.screen.cols <- c("DRUG_ID", "DRUG_NAME", "DRUG_SET", "SCREEN_ID", "SAMPLE_DATE", "PATIENT_ID")
experimental.factors <- c("DRUG_NAME", "DRUG_SET", "SCREEN_ID", "PATIENT_ID")

all.cols <- c(drug.screen.cols, conc.col, response.col)

## Limit fo AML
fimm.aml.metadata <- subset(fimm.metadata, grepl(diagnosis, pattern="AML"))
drug.data.long <- subset(drug.data.long, (SCREEN_ID %in% rownames(fimm.aml.metadata)))
exclude <- list()

## Store in
## FIMM_BEAT AML Collaboration/Files/FIMM Data/Processed Data
parentId <- "syn8270577"

output.file <- paste0("fimm.l4.and.ll4.fits.112017.tsv")

## END setup

