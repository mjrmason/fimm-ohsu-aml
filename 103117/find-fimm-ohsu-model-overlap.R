## rm(list = ls())

suppressPackageStartupMessages(library(glmnet))
suppressPackageStartupMessages(library("randomForest"))
suppressPackageStartupMessages(library("e1071"))
suppressPackageStartupMessages(library("foreach"))
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(tidyr))

suppressPackageStartupMessages(library("minpack.lm"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("gridExtra"))
suppressPackageStartupMessages(library("caTools"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("pcaMethods"))
suppressPackageStartupMessages(library("VennDiagram"))
suppressPackageStartupMessages(library("caret"))

suppressPackageStartupMessages(library(ggbeeswarm))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(synapseClient))

suppressPackageStartupMessages(library("GSVA"))
# use GSA to read in a gmt file
suppressPackageStartupMessages(library("GSA"))

## Bring in common code
source("../common/data-preprocessing.R")
source("../common/models.R")
source("../common/plotting.R")

## Register parallelization
num.cores <- detectCores()
if(!is.na(num.cores) && (num.cores > 1)) {
  suppressPackageStartupMessages(library("doMC"))
  cat(paste("Registering ", num.cores-1, " cores.\n", sep=""))
  registerDoMC(cores=(num.cores-1))
}
num.processes <- num.cores - 1

## Log in to Synapse
synapseLogin()

## Begin setup (this will need to be changed as the data sets in the intersection change)

## Purpose:
## This file models 2 data sets (OHSU and FIMM)

## Data sets to be analyzed:
data.sets <- c("ohsu", "fimm")
data.sets.to.analyze <- data.sets

## The drug mapping file lists the drugs in common across all data sets with the drug.map.drug.id.cols
## providing the identifier for each drug for each respective data set.
## The drug mapping file may be created by a script such as 
## harmonize-fimm-ohsu-compounds.R
drug.mapping.synId <- "syn11361110" ## fimm-ohsu-drug-map.tsv
drug.map.drug.id.cols <- c("OHSU_DRUG_NAME", "FIMM_Batch_ID_Drug")
names(drug.map.drug.id.cols) <- data.sets

## The column name within the drug map file that will be used to identify drugs (across all data sets)
drug.name.col <- "FIMM_DRUG_NAME"

## These fit files are assumed to be for AUCs/DSS values calculated across shared concentration ranges across data sets.
## These files may be created by a script such as
## calculate-dss-fimm-ohsu.R
data.set.dss.fit.synIds <- c("syn11362885", "syn11362891")
names(data.set.dss.fit.synIds) <- data.sets
data.set.dss.fit.files <- c("ohsu-fimm-ohsu-outlier1-dss-t0.tsv", "fimm-fimm-ohsu-outlier1-dss-t0.tsv")

## Batch-corrected log2 cpm (or rma) expression files (NB: these exclude 2 outliers in OHSU)
## These files may be created by a script such as
## normalize-fimm-ohsu-expression-distributions.R
## Batch-corrected log2 cpm (or rma) expression files (NB: these exclude 2 outliers in OHSU)
data.set.expr.synIds <- c("syn11362256", "syn11362257")
names(data.set.expr.synIds) <- data.sets
data.set.expr.files <- c("ohsu-expr-fimm-ohsu-outlier1-combat.tsv", "fimm-expr-fimm-ohsu-outlier1-combat.tsv")

## The columns within each respective data set that identifies the drug (in the name space of that data set)
data.set.drug.id.cols <- list("ohsu" = "inhibitor", "fimm" = "DRUG_ID")

## The columns within each respective data set that identifies a screen (i.e., a particular replicate of a patient/cell line x drug)
screen.id.cols <- list("ohsu" = c("patient_id", "inhibitor", "lab_id", "SeqID", "replicant", "time_of_read"),
                       "fimm" = c("SCREEN_ID", "PATIENT_ID", "DRUG_ID"))

## The columns within each respective data set that identifies the patient/cell line/specimen in the drug fit files.
## These patients should match the patient identifiers/columns of the expression data.
## NB: in some cases, these columns may not exist in the drug fit files, but need to be merged in.
## In particular, "SeqID" (the identifier used in sequencing experiments) may not be included in the OHSU drug fit files.
## Instead, patients are identified by "lab_id"'s in the drug fit file.  Below, we merge the OHSU metadata/sample
## summary information that maps lab_id's to SeqID's.
patient.id.cols <- list("ohsu" = "SeqID", "fimm" = "SCREEN_ID")

## The columns within the drug fit files that are effectively annotations and are not dependent on the fit [e.g., using
## a 4-parameter logistic L.4 or 4-parameter log logistic LL.4 fit].
merge.cols <- list(
  "ohsu" = c(screen.id.cols[["ohsu"]], "min.conc", "max.conc", "all.concs.nM", "uniq.concs.nM"),
  "fimm" = c(screen.id.cols[["fimm"]], "min.conc", "max.conc", "all.concs.nM", "uniq.concs.nM"))

## Pull in the CTRP cell line metadata so we can see the cancer/histology type of each sample
synId <- "syn5632192"
obj <- synGet(synId, downloadFile = TRUE)
ctrp.cell.line.metadata <- read.table(getFileLocation(obj), sep="\t", header=TRUE, as.is=TRUE, comment.char="", quote="")

## A string that will be included in any output files.
## Here "foc-aml-s-aml" = "fimm-ohsu-ctrp-aml-s-aml"
file.prefix <- "fimm-ohsu-outlier-ds"

## End setup

cat(paste0("Model AUC ~ expr for data sets: ", paste(data.sets, collapse=","), "\n"))

## Load in the fit files
fits <- llply(data.set.dss.fit.synIds, .parallel = FALSE,
              .fun = function(synId) {
                       obj <- synGet(synId, downloadFile = TRUE)
                       file <- getFileLocation(obj)
                       read.table(file, sep="\t", header=TRUE, as.is=TRUE, comment.char="", quote="\"")
              })

## Merge metadata to fit files
cat("Merging metadata to fits\n")
for(ds in names(fits)) {
  switch(ds,
    "ohsu" = {
      ## Add the SeqID to the OHSU drug data.
      ## Read in the rna-seq summary table, which provides a map from the lab_id's used
      ## in the drug response data to the SeqID's used in the expression data.
      synId <- "syn10083817"
      obj <- synGet(synId, downloadFile = TRUE)
      file <- getFileLocation(obj)
      ohsu.rnaseq.sample.summary <- read.table(file, header=TRUE, sep="\t")
## Looks like the SeqID has already been included in the fit file.
##      fits[[ds]] <- merge(fits[[ds]], ohsu.rnaseq.sample.summary, by.x = "lab_id", by.y = "Original_LabID")
    },
    "fimm" = {
    },
    "ctrp" = {
      fits[[ds]] <- merge(fits[[ds]], ctrp.cell.line.metadata[, c("master_ccl_id", "ccl_name")])
    },
    "ntap" = {
    },
    "sanger" = {
    },
    {
      stop("Unrecognized data set\n")
    })
}

## Load the map of drugs between data sets
obj <- synGet(drug.mapping.synId, downloadFile = TRUE)
drug.map <- read.table(getFileLocation(obj), sep="\t", header=TRUE, quote="\"", comment.char="")

## Drop any drugs that show up in multiple rows of the drug map--which would indicate an ambiguous mapping.
my.dup <- function(x) duplicated(x, fromLast = TRUE) | duplicated(x, fromLast = FALSE)
for(col in drug.map.drug.id.cols) {
  flag <- !my.dup(drug.map[, col])
  drug.map <- drug.map[flag, ]
}

## Restrict the fits to those drugs in common across the data sets -- do so by merging.
cat("Restricting fits to those involving common drugs\n")
for(i in 1:length(fits)) {
  ds <- names(fits)[i]
  fits[[ds]] <- merge(fits[[ds]], drug.map, by.x = data.set.drug.id.cols[[ds]], by.y = drug.map.drug.id.cols[i], all = FALSE)
}

common.drugs <- Reduce("intersect", lapply(fits, function(fit) unique(fit[, drug.name.col])))
for(i in 1:length(fits)) {
  ds <- names(fits)[i]
  flag <- fits[[ds]][, drug.name.col] %in% common.drugs
  fits[[ds]] <- fits[[ds]][flag, ]
}

## Limit the drug map to common drugs (which, by design of the drug map file, should already be the case)
drug.name.tbl <- drug.map[drug.map[, drug.name.col] %in% common.drugs,]

common.drugs.by.data.set <- list()
for(i in 1:length(fits)) {
  ds <- names(fits)[i]
  common.drugs.by.data.set[[ds]] <- drug.name.tbl[, drug.name.col]
}

## Load in the expr files
exprs <- llply(data.set.expr.synIds, .parallel = FALSE,
              .fun = function(synId) {
                       obj <- synGet(synId, downloadFile = TRUE)
                       file <- getFileLocation(obj)
                       read.table(file, sep="\t", header=TRUE, as.is=TRUE, comment.char="", quote="\"")
              })

## Sanitize some of the column headers 
for(nm in names(exprs)) {
  switch(nm,
    "ohsu" = {
      ## Convert the column names from X20.00347 to 20-00347
      colnames(exprs[[nm]]) <- gsub("X(.*)\\.(.*)", "\\1-\\2", colnames(exprs[[nm]]))
    },
    "fimm" = {
    },
    "ctrp" = {
    },
    "ntap" = {
    },
    "sanger" = {
      colnames(exprs[[nm]]) <- gsub("X(.*)", "\\1", colnames(exprs[[nm]]))
    },
    {
      stop("Unrecognized data set\n")
    })
}

common.genes <- Reduce("intersect", lapply(exprs, rownames))

## Restrict the matrices to those genes in common across data sets
for(i in 1:length(exprs)) {
  ds <- names(exprs)[i]
  exprs[[ds]] <- exprs[[ds]][common.genes, ]
}

## Sanity check: compare Sanger and CTRP drug AUC and IC50 values.
if(("sanger" %in% names(fits)) && ("ctrp" %in% names(fits))) {
  cat("Comparing Sanger and CTRP\n")
  synId <- "syn11289147"
  obj <- synGet(synId, downloadFile = TRUE)
  gdsc.to.ctrp.cell.line.map <- openxlsx:::read.xlsx(getFileLocation(obj), sheet = 1, startRow = 14)
  gdsc.to.ctrp.cell.line.map <- gdsc.to.ctrp.cell.line.map[, c("GDSC1000.cosmic.id", "CTRP.master.ccl.id")]
  ft <- fits[["sanger"]]
  ft <- merge(ft, gdsc.to.ctrp.cell.line.map, by.x = "COSMIC_ID", by.y = "GDSC1000.cosmic.id", all = FALSE)
  ft <- merge(ft, fits[["ctrp"]], by.x = c("CTRP.master.ccl.id", "cpd_name"), by.y = c("master_ccl_id", "cpd_name"), suffixes = c(".sanger", ".ctrp"))
  g1 <- plot.r2(ft$dss.auc.ll4.ctrp, ft$dss.auc.ll4.sanger)
  g1 <- plot.r2(ft$dss.auc.l4.ctrp, ft$dss.auc.l4.sanger)
  g1 <- g1 + xlab("CTRP AUC")
  g1 <- g1 + ylab("Sanger AUC")
  g2 <- plot.r2(ft$e.ll4.ctrp, ft$e.ll4.sanger)
  g2 <- g2 + xlab("CTRP IC50")
  g2 <- g2 + ylab("Sanger IC50")
  grid.arrange(g1, g2)

  if(FALSE) {
  gdsc.to.ctrp.cell.line.map <- read.table(getFileLocation(synGet("syn9988099", downloadFile = TRUE)), sep="\t", header=TRUE, quote="\"", comment.char="")
  ft <- fits[["sanger"]]
  ft <- merge(ft, gdsc.to.ctrp.cell.line.map, by.x = "COSMIC_ID", by.y = "COSMIC_ID", all = FALSE)
  ft$CCLE.name <- gsub(ft$CCLE.name, pattern="^([^_]+).*", replacement="\\1")
  ft <- ft[!is.na(ft$CCLE.name),]
  ft <- merge(ft, fits[["ctrp"]], by.x = c("CCLE.name", "cpd_name"), by.y = c("ccl_name.l4", "cpd_name"), suffixes = c(".sanger", ".ctrp"))

  flag <- (ft$dss.auc.ll4.ctrp < 2*10^5) & (ft$dss.auc.ll4.sanger < 2*10^5)
  g1 <- plot.r2(ft$dss.auc.ll4.ctrp[flag], ft$dss.auc.ll4.sanger[flag])
  g1 <- g1 + xlab("CTRP AUC")
  g1 <- g1 + ylab("Sanger AUC")
  flag <- (ft$e.ll4.ctrp < 50000) & (ft$e.ll4.sanger < 50000)
  g2 <- plot.r2(ft$e.ll4.ctrp[flag], ft$e.ll4.sanger[flag])
  g2 <- g2 + xlab("CTRP IC50")
  g2 <- g2 + ylab("Sanger IC50")
  grid.arrange(g1, g2)
  }
}

## Download biocarta, kegg, reactome, and hallmark data sets that have been translated to ensembl ids

## Hallmark gene sets: The list in this R object is called hallmark.gene.sets
synId <- "syn10644516"
load(getFileLocation(synGet(synId)))

## KEGG gene sets: The list in this R object is called kegg.gene.sets
synId <- "syn10644519"
load(getFileLocation(synGet(synId)))

## Biocarta gene sets: The list in this R object is called biocarta.gene.sets
synId <- "syn10644517"
load(getFileLocation(synGet(synId)))

## Reactome gene sets: The list in this R object is called reactome.gene.sets
synId <- "syn10644520"
load(getFileLocation(synGet(synId)))

## BEGIN modeling setup

## Do analysis at pathway- (hallmark, biocarta, kegg, and reactome) as well as at gene-level
gene.sets <- list("gene" = "gene", "hallmark" = hallmark.gene.sets, "biocarta" = biocarta.gene.sets, "kegg" = kegg.gene.sets, "reactome" = reactome.gene.sets)
## Do analysis only at gene-level
gene.sets <- list("gene" = "gene")

## Filter the gene expression data to restrict to highly-expressed genes
expr.filt.matrices <- list()
for(ds in names(exprs)) {
  iqrs <- unlist(apply(exprs[[ds]], 1, IQR))
  expr.filt.matrices[[ds]] <- exprs[[ds]][iqrs > median(iqrs[iqrs > 0]),]
}

study.genes <- (lapply(exprs, rownames))
study.filt.genes <- (lapply(expr.filt.matrices, rownames))
common.genes <- Reduce(intersect, study.filt.genes)

## END modeling setup

## BEGIN Define training and test data sets

## Based on whether the data sets are listed in data.sets.to.analyze (defined above), the following code will define

## Training sets:
## 60% of OHSU (stratified by drug)
## CTRP

## Test sets:
## 40% of OHSU (stratified by drug)
## FIMM
## Sanger
## NTAP

train.set.names <- list()
train.dss.args <- list()
train.expr.args <- list()
train.genomic.args <- list()
train.clinical.args <- list()
train.common.drugs <- list()
train.drugs <- list()
train.drug.cols <- list()
train.patient.cols <- list()
train.response.cols <- list()

test.set.names <- list()
test.dss.args <- list()
test.expr.args <- list()
test.genomic.args <- list()
test.clinical.args <- list()
test.common.drugs <- list()
test.drug.cols <- list()
test.patient.cols <- list()
test.response.cols <- list()

ds <- "ohsu"
if(ds %in% data.sets.to.analyze) {

  nxt.indx <- length(train.set.names) + 1
  train.set.names[[nxt.indx]] <- "ohsu"
  train.dss.args[[nxt.indx]] <- fits[[ds]]
  train.expr.args[[nxt.indx]] <- expr.filt.matrices[[ds]][common.genes, ] 
  train.genomic.args[nxt.indx] <- list(NULL)
  train.clinical.args[nxt.indx] <- list(NULL)
  train.common.drugs[[nxt.indx]] <- common.drugs.by.data.set[[ds]]
  train.drugs[[nxt.indx]] <- common.drugs.by.data.set[[ds]]
  train.drug.cols[[nxt.indx]] <- drug.name.col
  train.patient.cols[[nxt.indx]] <- patient.id.cols[[ds]]
  train.response.cols[[nxt.indx]] <- "dss.auc.ll4"

  half.index <- createDataPartition(fits[[ds]][, data.set.drug.id.cols[[ds]]], p = .5, list = FALSE)
  fits.set1 <- fits[[ds]][half.index, ]
  fits.set2 <- fits[[ds]][-half.index, ]

  nxt.indx <- length(train.set.names) + 1
  train.set.names[[nxt.indx]] <- "ohsu.set1"
  train.dss.args[[nxt.indx]] <- fits.set1
  train.expr.args[[nxt.indx]] <- expr.filt.matrices[[ds]][common.genes, ] 
  train.genomic.args[nxt.indx] <- list(NULL)
  train.clinical.args[nxt.indx] <- list(NULL)
  train.common.drugs[[nxt.indx]] <- common.drugs.by.data.set[[ds]]
  train.drugs[[nxt.indx]] <- common.drugs.by.data.set[[ds]]
  train.drug.cols[[nxt.indx]] <- drug.name.col
  train.patient.cols[[nxt.indx]] <- patient.id.cols[[ds]]
  train.response.cols[[nxt.indx]] <- "dss.auc.ll4"

  nxt.indx <- length(train.set.names) + 1
  train.set.names[[nxt.indx]] <- "ohsu.set2"
  train.dss.args[[nxt.indx]] <- fits.set2
  train.expr.args[[nxt.indx]] <- expr.filt.matrices[[ds]][common.genes, ] 
  train.genomic.args[nxt.indx] <- list(NULL)
  train.clinical.args[nxt.indx] <- list(NULL)
  train.common.drugs[[nxt.indx]] <- common.drugs.by.data.set[[ds]]
  train.drugs[[nxt.indx]] <- common.drugs.by.data.set[[ds]]
  train.drug.cols[[nxt.indx]] <- drug.name.col
  train.patient.cols[[nxt.indx]] <- patient.id.cols[[ds]]
  train.response.cols[[nxt.indx]] <- "dss.auc.ll4"


} ## if("ohsu" %in% data.sets.to.analyze) 

ds <- "fimm"
if(ds %in% data.sets.to.analyze) {

  nxt.indx <- length(train.set.names) + 1
  train.set.names[[nxt.indx]] <- ds
  train.dss.args[[nxt.indx]] <- fits[[ds]]
  train.expr.args[[nxt.indx]] <- expr.filt.matrices[[ds]][common.genes, ]
  train.genomic.args[nxt.indx] <- list(NULL)
  train.clinical.args[nxt.indx] <- list(NULL)
  train.common.drugs[[nxt.indx]] <- common.drugs.by.data.set[[ds]]
  train.drugs[[nxt.indx]] <- common.drugs.by.data.set[[ds]]
  train.drug.cols[[nxt.indx]] <- drug.name.col
  train.patient.cols[[nxt.indx]] <- patient.id.cols[[ds]]
  train.response.cols[[nxt.indx]] <- "dss.auc.ll4"

  half.index <- createDataPartition(fits[[ds]][, data.set.drug.id.cols[[ds]]], p = .5, list = FALSE)
  fits.set1 <- fits[[ds]][half.index, ]
  fits.set2 <- fits[[ds]][-half.index, ]

  nxt.indx <- length(train.set.names) + 1
  train.set.names[[nxt.indx]] <- "fimm.set1"
  train.dss.args[[nxt.indx]] <- fits.set1
  train.expr.args[[nxt.indx]] <- expr.filt.matrices[[ds]][common.genes, ] 
  train.genomic.args[nxt.indx] <- list(NULL)
  train.clinical.args[nxt.indx] <- list(NULL)
  train.common.drugs[[nxt.indx]] <- common.drugs.by.data.set[[ds]]
  train.drugs[[nxt.indx]] <- common.drugs.by.data.set[[ds]]
  train.drug.cols[[nxt.indx]] <- drug.name.col
  train.patient.cols[[nxt.indx]] <- patient.id.cols[[ds]]
  train.response.cols[[nxt.indx]] <- "dss.auc.ll4"

  nxt.indx <- length(train.set.names) + 1
  train.set.names[[nxt.indx]] <- "fimm.set2"
  train.dss.args[[nxt.indx]] <- fits.set2
  train.expr.args[[nxt.indx]] <- expr.filt.matrices[[ds]][common.genes, ] 
  train.genomic.args[nxt.indx] <- list(NULL)
  train.clinical.args[nxt.indx] <- list(NULL)
  train.common.drugs[[nxt.indx]] <- common.drugs.by.data.set[[ds]]
  train.drugs[[nxt.indx]] <- common.drugs.by.data.set[[ds]]
  train.drug.cols[[nxt.indx]] <- drug.name.col
  train.patient.cols[[nxt.indx]] <- patient.id.cols[[ds]]
  train.response.cols[[nxt.indx]] <- "dss.auc.ll4"

} ## if("fimm" %in% data.sets.to.analyze) 

ds <- "ctrp"
if(ds %in% data.sets.to.analyze) {

  nxt.indx <- length(train.set.names) + 1
  train.set.names[[nxt.indx]] <- ds
  train.dss.args[[nxt.indx]] <- fits[[ds]]
  train.expr.args[[nxt.indx]] <- expr.filt.matrices[[ds]][common.genes, ]
  train.genomic.args[nxt.indx] <- list(NULL)
  train.clinical.args[nxt.indx] <- list(NULL)
  train.common.drugs[[nxt.indx]] <- common.drugs.by.data.set[[ds]]
  train.drugs[[nxt.indx]] <- common.drugs.by.data.set[[ds]]
  train.drug.cols[[nxt.indx]] <- drug.name.col
  train.patient.cols[[nxt.indx]] <- patient.id.cols[[ds]]
  train.response.cols[[nxt.indx]] <- "dss.auc.ll4"

} ## if("ctrp" %in% data.sets.to.analyze) 

ds <- "ntap"
if(ds %in% data.sets.to.analyze) {

  nxt.indx <- length(train.set.names) + 1
  train.set.names[[nxt.indx]] <- ds
  train.dss.args[[nxt.indx]] <- fits[[ds]]
  train.expr.args[[nxt.indx]] <- expr.filt.matrices[[ds]][common.genes, ]
  train.genomic.args[nxt.indx] <- list(NULL)
  train.clinical.args[nxt.indx] <- list(NULL)
  train.common.drugs[[nxt.indx]] <- common.drugs.by.data.set[[ds]]
  train.drugs[[nxt.indx]] <- common.drugs.by.data.set[[ds]]
  train.drug.cols[[nxt.indx]] <- drug.name.col
  train.patient.cols[[nxt.indx]] <- patient.id.cols[[ds]]
  train.response.cols[[nxt.indx]] <- "dss.auc.ll4"
} ## if("ntap" %in% data.sets.to.analyze) 

ds <- "sanger"
if(ds %in% data.sets.to.analyze) {

  nxt.indx <- length(train.set.names) + 1
  train.set.names[[nxt.indx]] <- ds
  train.dss.args[[nxt.indx]] <- fits[[ds]]
  train.expr.args[[nxt.indx]] <- expr.filt.matrices[[ds]][common.genes, ]
  train.genomic.args[nxt.indx] <- list(NULL)
  train.clinical.args[nxt.indx] <- list(NULL)
  train.common.drugs[[nxt.indx]] <- common.drugs.by.data.set[[ds]]
  train.drugs[[nxt.indx]] <- common.drugs.by.data.set[[ds]]
  train.drug.cols[[nxt.indx]] <- drug.name.col
  train.patient.cols[[nxt.indx]] <- patient.id.cols[[ds]]
  train.response.cols[[nxt.indx]] <- "dss.auc.ll4"
} ## if("sanger" %in% data.sets.to.analyze) 

## END Define training and test data sets

## Create the table of number of samples (to use for downsampling) vs drug.
## To do so, look for the non-NA samples in common between expression and drug
## This is only used if we do downsampling.
num.samples.per.drug <- drug.name.tbl
for(train.indx in 1:length(train.set.names)) {
   train.set <- train.set.names[[train.indx]]
   l <- prepare.drug.response.and.expr.matrices(train.dss.args[[train.indx]], train.expr.args[[train.indx]], drugs = train.common.drugs[[train.indx]], 
                                                drug.col = train.drug.cols[[train.indx]], patient.col = train.patient.cols[[train.indx]], 
                                                response.col = train.response.cols[[train.indx]])

   drc <- l[["drc.df"]]
   expr <- l[["expr.df"]]
   samples <- colnames(drc)

   if(!is.null(expr)) {
     samples <- intersect(samples, colnames(expr))
   }

   col <- paste0(train.set, ".train.num.samples")
   num.samples.per.drug[, col] <- NA
   for(i in 1:nrow(num.samples.per.drug)) {
     if(is.na(num.samples.per.drug[i, train.drug.cols[[train.indx]]])) { next }
     vec <- drc[as.character(num.samples.per.drug[i, train.drug.cols[[train.indx]]]),]
     flag <- !is.na(vec)
     drug.samples <- intersect(samples, names(vec)[flag])
     num.samples.per.drug[i, col] <- length(unique(drug.samples))
   }
}

if(length(test.set.names) > 0) {
  for(test.indx in 1:length(test.set.names)) {
     test.set <- test.set.names[[test.indx]]
     l <- prepare.drug.response.and.expr.matrices(test.dss.args[[test.indx]], test.expr.args[[test.indx]], drugs = test.common.drugs[[test.indx]], 
                                                  drug.col = test.drug.cols[[test.indx]], patient.col = test.patient.cols[[test.indx]], 
                                                  response.col = test.response.cols[[test.indx]])
  
     drc <- l[["drc.df"]]
     expr <- l[["expr.df"]]
     samples <- colnames(drc)
  
     if(!is.null(expr)) {
       samples <- intersect(samples, colnames(expr))
     }
     col <- paste0(test.set, ".test.num.samples")
     num.samples.per.drug[, col] <- NA
     for(i in 1:nrow(num.samples.per.drug)) {
       if(is.na(num.samples.per.drug[i, test.drug.cols[[test.indx]]])) { next }
       vec <- drc[as.character(num.samples.per.drug[i, test.drug.cols[[test.indx]]]),]
       flag <- !is.na(vec)
       drug.samples <- intersect(samples, names(vec)[flag])
       num.samples.per.drug[i, col] <- length(unique(drug.samples))
     }
  }
}

## The column named 'num.samples' will be used to define the number of samples (per drug) used for downsampling.
## Define that here, separately for both training and testing.
min.cnt.col <- "num.samples"
num.train.samples.per.drug <- num.samples.per.drug
num.test.samples.per.drug <- num.samples.per.drug

## The number of samples to use for downsampling for training should be the minimum of all samples across all training sets
## (defined on a per-drug basis)
cnt.cols <- grepl(colnames(num.train.samples.per.drug), pattern="train.num.samples")
num.train.samples.per.drug[, min.cnt.col] <- unlist(apply(num.train.samples.per.drug[, cnt.cols, drop=FALSE], 1, min))

## The number of samples to use for downsampling for testing should be the minimum of all samples across all testing sets
## (defined on a per-drug basis)
cnt.cols <- grepl(colnames(num.test.samples.per.drug), pattern="test.num.samples")
num.test.samples.per.drug[, min.cnt.col] <- unlist(apply(num.test.samples.per.drug[, cnt.cols, drop=FALSE], 1, min))

## BEGIN modeling

cat("Training and testing with AUC\n")

## Set the response to use AUC based on 4-parameter log logistic fit (LL.4)
train.response.cols <- rep("dss.auc.ll4", length(train.response.cols))
test.response.cols <- rep("dss.auc.ll4", length(test.response.cols))

## Train on each of the training sets and apply those models to each of the test data sets.
## At both the training and test stages downsample to the number of samples (which vary by drug)
## in the num.{train,test}.samples.per.drug table.
## num.iterations determines the number of down-sampling.
## with.replacement = FALSE indicates downsampling, whereas = TRUE indicates bootstrapping.
## NB: can return fits by setting return.fits = TRUE, but this will take a lot of memory (and disk space when doing a save.image)
res.aucs <- list()
do.downsampling <- FALSE
use.rf <- TRUE
if(do.downsampling) {
  res.aucs[[1]] <- train.and.test.crossproduct.with.downsampling(gene.sets = gene.sets, drug.name.tbl = drug.name.tbl, train.set.names = train.set.names, train.dss.args = train.dss.args, 
                                         train.expr.args = train.expr.args, train.genomic.args = train.genomic.args, train.clinical.args = train.clinical.args, 
                                         train.common.drugs = train.common.drugs, train.drugs = train.drugs, train.drug.cols = train.drug.cols, train.patient.cols = train.patient.cols, 
                                         train.response.cols = train.response.cols, test.set.names = test.set.names, test.dss.args = test.dss.args, test.expr.args = test.expr.args, 
                                         test.genomic.args = test.genomic.args, test.clinical.args = test.clinical.args, test.common.drugs = test.common.drugs,
                                         test.drug.cols = test.drug.cols, test.patient.cols = test.patient.cols, test.response.cols = test.response.cols, 
                                         seed = 1234, use.rf = use.rf, use.svm = FALSE, use.mean = FALSE, num.processes = num.processes, 
                                         num.train.samples.per.drug = num.train.samples.per.drug, num.test.samples.per.drug = num.test.samples.per.drug,
                                         num.iterations = 100, with.replacement = FALSE, return.fits = FALSE) 
} else { 
  ## Alternately, can do modeling without downsampling
  seeds <- c(1234, 4321)
  for(i in 1:length(seeds)) {
    seed <- seeds[i]
    cat(paste0("Training for seed ", seed, "\n"))
    res.aucs[[i]] <- train.and.test.crossproduct(gene.sets = gene.sets, drug.name.tbl = drug.name.tbl, train.set.names = train.set.names, train.dss.args = train.dss.args, 
                                           train.expr.args = train.expr.args, train.genomic.args = train.genomic.args, train.clinical.args = train.clinical.args, 
                                           train.common.drugs = train.common.drugs, train.drugs = train.drugs, train.drug.cols = train.drug.cols, train.patient.cols = train.patient.cols, 
                                           train.response.cols = train.response.cols, test.set.names = test.set.names, test.dss.args = test.dss.args, test.expr.args = test.expr.args, 
                                           test.genomic.args = test.genomic.args, test.clinical.args = test.clinical.args, test.common.drugs = test.common.drugs,
                                           test.drug.cols = test.drug.cols, test.patient.cols = test.patient.cols, test.response.cols = test.response.cols, 
                                           seed = seed, use.rf = use.rf, use.svm = TRUE, use.mean = FALSE, num.processes = num.processes)
  }
}
cat("Done training and testing with AUC\n")

save.image(".Rdata.overlap.fimm.ohsu")

do.shifted.auc <- FALSE
if(do.shifted.auc) {
  ## Set the response to use "shifted" AUC based on 4-parameter log logistic fit (LL.4)
  train.response.cols <- rep("shifted.auc.ll4", length(train.response.cols))
  test.response.cols <- rep("shifted.auc.ll4", length(test.response.cols))
  res.shifted.auc <- train.and.test.crossproduct.with.downsampling(gene.sets = gene.sets, drug.name.tbl = drug.name.tbl, train.set.names = train.set.names, train.dss.args = train.dss.args, 
                                         train.expr.args = train.expr.args, train.genomic.args = train.genomic.args, train.clinical.args = train.clinical.args, 
                                         train.common.drugs = train.common.drugs, train.drugs = train.drugs, train.drug.cols = train.drug.cols, train.patient.cols = train.patient.cols, 
                                         train.response.cols = train.response.cols, test.set.names = test.set.names, test.dss.args = test.dss.args, test.expr.args = test.expr.args, 
                                         test.genomic.args = test.genomic.args, test.clinical.args = test.clinical.args, test.common.drugs = test.common.drugs,
                                         test.drug.cols = test.drug.cols, test.patient.cols = test.patient.cols, test.response.cols = test.response.cols, 
                                         seed = 1234, use.rf = use.rf, use.svm = FALSE, use.mean = FALSE, num.processes = num.processes,
                                         num.train.samples.per.drug = num.train.samples.per.drug, num.test.samples.per.drug = num.test.samples.per.drug,
                                         num.iterations = 100, with.replacement = FALSE, return.fits = FALSE) 
  cat("Done training and testing with shifted AUC\n")

##  save.image(".Rdata.overlap.fimm.ohsu")
}

trn.tbl <- ldply(train.set.names,
             .fun = function(train.set) {
                      ldply(names(gene.sets),
                        .fun = function(gene.set) {
                                 ldply(res.aucs[[1]][["all.fits"]][[train.set]][[gene.set]],
                                   .fun = function(lst) {
                                     ldply(lst, 
                                       .fun = function(df) {
                                         ## Only pull out the regression and rf results
                                         if((df$model != "glmnet") && (df$model != "rf")) { return(NULL) }
                                         coeffs <- NULL
                                         switch(df$model,
                                           "glmnet" = {
                                             coeffs <- coefficients(df$fit)
                                           },
                                           "rf" = {
                                             col <- colnames(importance(df$fit))[1]
                                             if(!grepl(col, pattern="IncMSE")) { stop(paste0("Was expected ", col, " to be %IncMSE\n")) }
                                             coeffs <- importance(df$fit)[,1,drop=FALSE]
                                           },
                                           { stop(paste0("Unknown model ", df$model, "\n")) })
                                         ns <- rownames(coeffs)
                                         coeffs <- as.vector(coeffs)
                                         names(coeffs) <- ns
                                         coeffs <- coeffs[!grepl(pattern="Intercept", names(coeffs))]
                                         n.features <- length(coeffs) 
                                         coeffs <- coeffs[order(names(coeffs))]
                                         coeff.vals <- paste(coeffs, collapse=",")
                                         coeffs <- paste(names(coeffs), collapse=",")
                                         ret <- data.frame(train.set = train.set, gene.set = gene.set, train.drug = df$train.drug,
                                                           alpha = df$alpha, model = df$model, n.features = n.features, coeffs = coeffs, coeff.vals = coeff.vals)
                                         ret
                                       })
                                   })
                        })
             })
trn.tbl$model <- as.character(trn.tbl$model)
trn.tbl[(trn.tbl$model == "glmnet") & (trn.tbl$alpha == 0),"model"] <- "ridge"
trn.tbl[(trn.tbl$model == "glmnet") & (trn.tbl$alpha == 1),"model"] <- "lasso"
trn.tbl <- subset(trn.tbl, model != "glmnet")
trn.tbl$train.drug <- as.character(trn.tbl$train.drug)
trn.tbl$model <- as.character(trn.tbl$model)
trn.tbl$coeffs <- as.character(trn.tbl$coeffs)

trn.tbl.lasso <- subset(trn.tbl, model == "lasso")
trn.tbl.ridge <- subset(trn.tbl, model == "ridge")
trn.tbl.rf <- subset(trn.tbl, model == "rf")

## save.image(".Rdata.overlap.fimm.ohsu")

cat("Done subsetting table by model\n")

## Looking at correlation of coefficients
## method = c("pearson", "spearman", "kendall")
find.overlap.of.dense.predictors <- function(tbl.input, method) {
  ret.tbl <- c()
  for(gs in unique(tbl.input$gene.set)) {
    tmp <- subset(tbl.input, gene.set == gs)
    tr.sets <- as.character(unique(tmp$train.set))
    for(i in 1:(length(tr.sets)-1)) {
      tr1 <- tr.sets[i]
      tmp1 <- subset(tmp, train.set == tr1)
      merge.col <- drug.name.col
      tmp1 <- merge(drug.name.tbl, tmp1, by.x = merge.col, by.y = "train.drug")
      for(j in (i+1):length(tr.sets)) {
        tr2 <- tr.sets[j]
        if(((tr1 == "ohsu") || (tr2 == "ohsu")) && (grepl(pattern="ohsu.set1", paste0(tr1,tr2)) || grepl(pattern="ohsu.set2", paste0(tr1,tr2)))) { next }
        if(((tr1 == "fimm") || (tr2 == "fimm")) && (grepl(pattern="fimm.set1", paste0(tr1,tr2)) || grepl(pattern="fimm.set2", paste0(tr1,tr2)))) { next }
        tmp2 <- subset(tmp, train.set == tr2)
        merge.col <- drug.name.col
        tmp2 <- merge(drug.name.tbl, tmp2, by.x = merge.col, by.y = "train.drug")
  
        m <- merge(tmp1, tmp2, by = merge.col, suffixes = c(".1", ".2"))
        m <- m[, c(merge.col, "n.features.1", "coeffs.1", "coeffs.2", "coeff.vals.1", "coeff.vals.2", "train.set.1", "train.set.2", "gene.set.1")]
        m$id <- 1:nrow(m)
        hyp <- ddply(m, .variables = "id",
               .fun = function(r) {
                        coeff.vals.1 <- as.numeric(unlist(strsplit(as.character(r$coeff.vals.1), split=",")))
                        coeffs.1 <- (unlist(strsplit(as.character(r$coeffs.1), split=",")))
                        names(coeff.vals.1) <- coeffs.1
                        coeffs.1 <- coeff.vals.1
                        coeff.vals.2 <- as.numeric(unlist(strsplit(as.character(r$coeff.vals.2), split=",")))
                        coeffs.2 <- (unlist(strsplit(as.character(r$coeffs.2), split=",")))
                        names(coeff.vals.2) <- coeffs.2
                        coeffs.2 <- coeff.vals.2
                        inter <- intersect(names(coeffs.1), names(coeffs.2))
                        ct <- cor.test(coeffs.1[inter], coeffs.2[inter], method = method)
                        pval <- ct$p.value
                        cor <- ct$estimate
                        vec <- c(r$id[1], cor, pval)
                        names(vec) <- c("id", "cor", "pval")
                        vec
                      })
        hyp <- merge(hyp, m, by = "id")
        ret.tbl <- rbind(ret.tbl, hyp)
      }
    }
  }
  ret.tbl
}

## transform = c("none", "abs", "negate")
## top = # of features to consider in overall
## top = 0 -> use all
find.overlap.of.sparse.predictors <- function(tbl.input, transform, top) {
  ret.tbl <- c()
  for(gs in unique(tbl.input$gene.set)) {
    tmp <- subset(tbl.input, gene.set == gs)
    tr.sets <- as.character(unique(tmp$train.set))
    for(i in 1:(length(tr.sets)-1)) {
      tr1 <- tr.sets[i]
      tmp1 <- subset(tmp, train.set == tr1)
      merge.col <- drug.name.col
      tmp1 <- merge(drug.name.tbl, tmp1, by.x = merge.col, by.y = "train.drug")
      for(j in (i+1):length(tr.sets)) {
        tr2 <- tr.sets[j]
        if(((tr1 == "ohsu") || (tr2 == "ohsu")) && (grepl(pattern="ohsu.set1", paste0(tr1,tr2)) || grepl(pattern="ohsu.set2", paste0(tr1,tr2)))) { next }
        if(((tr1 == "fimm") || (tr2 == "fimm")) && (grepl(pattern="fimm.set1", paste0(tr1,tr2)) || grepl(pattern="fimm.set2", paste0(tr1,tr2)))) { next }
        tmp2 <- subset(tmp, train.set == tr2)
        merge.col <- drug.name.col
        tmp2 <- merge(drug.name.tbl, tmp2, by.x = merge.col, by.y = "train.drug")

        m <- merge(tmp1, tmp2, by = merge.col, suffixes = c(".1", ".2"))
        m <- m[, c(merge.col, "n.features.1", "coeffs.1", "coeffs.2", "coeff.vals.1", "coeff.vals.2", "train.set.1", "train.set.2", "gene.set.1")]
        m$id <- 1:nrow(m)
        hyp <- ddply(m, .variables = "id", 
               .fun = function(r) {
                        n.f <- as.numeric(r$n.features.1)
                        coeff.vals.1 <- as.numeric(unlist(strsplit(as.character(r$coeff.vals.1), split=",")))
                        coeffs.1 <- (unlist(strsplit(as.character(r$coeffs.1), split=",")))
                        names(coeff.vals.1) <- coeffs.1
                        coeffs.1 <- coeff.vals.1
                        coeff.vals.2 <- as.numeric(unlist(strsplit(as.character(r$coeff.vals.2), split=",")))
                        coeffs.2 <- (unlist(strsplit(as.character(r$coeffs.2), split=",")))
                        names(coeff.vals.2) <- coeffs.2
                        coeffs.2 <- coeff.vals.2
                        switch(transform,
                          "drop.zero" = { 
                             coeffs.1 <- coeffs.1[coeffs.1 > 0]
                             coeffs.2 <- coeffs.2[coeffs.2 > 0]
                          },
                          "none" = { },
                          "abs" = {
                             coeffs.1 <- abs(coeffs.1)
                             coeffs.2 <- abs(coeffs.2)
                          },
                          "negate" = {
                             coeffs.1 <- - coeffs.1
                             coeffs.2 <- - coeffs.2
                          })
                        if(top != 0) {
                          coeffs.1 <- sort(coeffs.1)
                          coeffs.2 <- sort(coeffs.2)
                          coeffs.1 <- coeffs.1[1:min(top, length(coeffs.1))]
                          coeffs.2 <- coeffs.2[1:min(top, length(coeffs.2))]
                        }
                        n.1 <- length(coeffs.1)
                        n.2 <- length(coeffs.2)
                        both <- intersect(names(coeffs.1), names(coeffs.2))
                        n.both <- length(both)
                        both <- paste(both, collapse=",")
                        ## pval <- phyper(q = n.both - 1, m = n.1, n = n.f - n.1, k = n.2, lower.tail = FALSE) 
                        fet <- fisher.test(cbind(c(n.both, n.1), c(n.2, n.f - n.1 - n.2 + n.both)), alternative = "greater")
                        pval <- fet$p.value
                        vec <- c(r$id[1], n.1, n.2, n.both, both, pval)
                        names(vec) <- c("id", "n.1", "n.2", "n.both", "both", "pval")
                        vec
                      })
        hyp <- merge(hyp, m, by = "id")
        ret.tbl <- rbind(ret.tbl, hyp)
      }
    }
  }
  ret.tbl
}

cat("Summarizing ridge results\n")
all.ridge <- find.overlap.of.dense.predictors(trn.tbl.ridge, method = "spearman")
all.ridge$pval <- as.numeric(all.ridge$pval)
cat("Done summarizing ridge results\n")

if(use.rf) {
  cat("Summarizing rf results\n")
  ## all.rf <- find.overlap.of.dense.predictors(trn.tbl.rf, method = "spearman")
  all.rf <- find.overlap.of.sparse.predictors(trn.tbl.rf, transform = "negate", top = 100)
  all.rf$pval <- as.numeric(all.rf$pval)
  cat("Done summarizing rf results\n")
}

lst <- list("ridge" = all.ridge)
all.cross.sigs <- list()
for(nm in names(lst)) {

  all.sig <- subset(lst[[nm]], pval < 0.05 & cor > 0)
  cross.sig <- subset(all.sig, train.set.1 == "ohsu" & train.set.2 == "fimm")
  ohsu.sig <- subset(all.sig, (train.set.1 == "ohsu.set1" & train.set.2 == "ohsu.set2") | (train.set.2 == "ohsu.set1" & train.set.1 == "ohsu.set2"))
  fimm.sig <- subset(all.sig, (train.set.1 == "fimm.set1" & train.set.2 == "fimm.set2") | (train.set.2 == "fimm.set1" & train.set.1 == "fimm.set2"))
  all.cross.sigs[[nm]] <- cross.sig

  tbl <- as.data.frame(table(cross.sig$gene.set.1))
  colnames(tbl) <- c("Feature Set", "Drugs w/ Correlated Features")
  pdf(paste0(file.prefix, "-fimm-vs-ohsu-", nm, "-feature-overlap.pdf"))
  plot.table(tbl, main = "OHSU (train) vs FIMM (train)")
  d <- dev.off()

  tbl <- as.data.frame(table(ohsu.sig$gene.set.1))
  colnames(tbl) <- c("Feature Set", "Drugs w/ Correlated Features")
  pdf(paste0(file.prefix, "-ohsu1-vs-ohsu2-", nm, "-feature-overlap.pdf"))
  plot.table(tbl, main = "OHSU 1 (train) vs OHSU 2 (train)")
  d <- dev.off()

  tbl <- as.data.frame(table(fimm.sig$gene.set.1))
  colnames(tbl) <- c("Feature Set", "Drugs w/ Correlated Features")
  pdf(paste0(file.prefix, "-fimm1-vs-fimm2-", nm, "-feature-overlap.pdf"))
  plot.table(tbl, main = "FIMM 1 (train) vs FIMM 2 (train)")
  d <- dev.off()
}

## save.image(".Rdata.overlap.fimm.ohsu")

all <- find.overlap.of.sparse.predictors(trn.tbl.lasso, transform = "drop.zero", top = 0)
all$pval <- as.numeric(all$pval)
cat("Done summarizing lasso results\n")

## save.image(".Rdata.overlap.fimm.ohsu")

lst <- list("lasso" = all)
if(use.rf) {
  lst <- list("lasso" = all, "rf" = all.rf)
}
all.cross.sigs <- list()
for(nm in names(lst)) {
  all.sig <- subset(lst[[nm]], pval < 0.05)
  cross.sig <- subset(all.sig, train.set.1 == "ohsu" & train.set.2 == "fimm")
  all.cross <- subset(lst[[nm]], train.set.1 == "ohsu" & train.set.2 == "fimm")
  ohsu.sig <- subset(all.sig, (train.set.1 == "ohsu.set1" & train.set.2 == "ohsu.set2") | (train.set.2 == "ohsu.set1" & train.set.1 == "ohsu.set2"))
  fimm.sig <- subset(all.sig, (train.set.1 == "fimm.set1" & train.set.2 == "fimm.set2") | (train.set.2 == "fimm.set1" & train.set.1 == "fimm.set2"))
  all.cross.sigs[[nm]] <- cross.sig

  table(cross.sig$gene.set.1)
  table(ohsu.sig$gene.set.1)
  table(fimm.sig$gene.set.1)

  tbl <- as.data.frame(table(cross.sig$gene.set.1))
  colnames(tbl) <- c("Feature Set", "Drugs w/ Overlapping Features")
  pdf(paste0(file.prefix, "-ohsu-vs-fimm-lasso-feature-overlap.pdf"))
  plot.table(tbl, main = "OHSU (train) vs FIMM (train)")
  d <- dev.off()

  if(nm == "lasso") { 
    s <- subset(all.cross, n.1 > 0 & n.2 > 0)
    tmp <- s[,c("FIMM_DRUG_NAME", "train.set.1", "n.1", "train.set.2", "n.2", "gene.set.1", "n.both")]
    write.table(file = paste0(file.prefix, "-ohsu-vs-fimm-lasso-feature-overlap.tsv"), tmp, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  }

  tbl <- as.data.frame(table(ohsu.sig$gene.set.1))
  colnames(tbl) <- c("Feature Set", "Drugs w/ Overlapping Features")
  pdf(paste0(file.prefix, "-ohsu1-vs-ohsu2-lasso-feature-overlap.pdf"))
  plot.table(tbl, main = "OHSU 1 (train) vs OHSU 2 (train)")
  d <- dev.off()

  tbl <- as.data.frame(table(fimm.sig$gene.set.1))
  colnames(tbl) <- c("Feature Set", "Drugs w/ Overlapping Features")
  pdf(paste0(file.prefix, "-fimm1-vs-fimm2-lasso-feature-overlap.pdf"))
  plot.table(tbl, main = "FIMM 1 (train) vs FIMM 2 (train)")
  d <- dev.off()
}


save.image(".Rdata.overlap.fimm.ohsu")
stop("successfully completed")
