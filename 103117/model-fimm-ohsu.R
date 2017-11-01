
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
## gene.sets <- list("gene" = "gene", "hallmark" = hallmark.gene.sets, "biocarta" = biocarta.gene.sets, "kegg" = kegg.gene.sets, "reactome" = reactome.gene.sets)
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
  set.seed(1234)

  ## Use 60% of samples (stratified by drug) for training and 40% for test.
  train.index <- createDataPartition(fits[[ds]][, data.set.drug.id.cols[[ds]]], p = .6, list = FALSE)
  ohsu.dss.orig.train <- fits[[ds]][train.index, ]
  ohsu.dss.orig.test <- fits[[ds]][-train.index, ]

  nxt.indx <- length(train.set.names) + 1
  train.set.names[[nxt.indx]] <- "ohsu.train"
  train.dss.args[[nxt.indx]] <- ohsu.dss.orig.train
  train.expr.args[[nxt.indx]] <- expr.filt.matrices[[ds]][common.genes, ] 
  train.genomic.args[nxt.indx] <- list(NULL)
  train.clinical.args[nxt.indx] <- list(NULL)
  train.common.drugs[[nxt.indx]] <- common.drugs.by.data.set[[ds]]
  train.drugs[[nxt.indx]] <- common.drugs.by.data.set[[ds]]
  train.drug.cols[[nxt.indx]] <- drug.name.col
  train.patient.cols[[nxt.indx]] <- patient.id.cols[[ds]]
  train.response.cols[[nxt.indx]] <- "dss.auc.ll4"

  nxt.indx <- length(test.set.names) + 1
  test.set.names[[nxt.indx]] <- "ohsu.test"
  test.dss.args[[nxt.indx]] <- ohsu.dss.orig.test
  test.expr.args[[nxt.indx]] <- expr.filt.matrices[[ds]][common.genes, ]
  test.genomic.args[nxt.indx] <- list(NULL)
  test.clinical.args[nxt.indx] <- list(NULL)
  test.common.drugs[[nxt.indx]] <- common.drugs.by.data.set[[ds]]
  test.drug.cols[[nxt.indx]] <- drug.name.col
  test.patient.cols[[nxt.indx]] <- patient.id.cols[[ds]]
  test.response.cols[[nxt.indx]] <- "dss.auc.ll4"

  rm(ohsu.dss.orig.train)
  rm(ohsu.dss.orig.test)
} ## if("ohsu" %in% data.sets.to.analyze) 

ds <- "fimm"
if(ds %in% data.sets.to.analyze) {

  nxt.indx <- length(test.set.names) + 1
  test.set.names[[nxt.indx]] <- ds
  test.dss.args[[nxt.indx]] <- fits[[ds]]
  test.expr.args[[nxt.indx]] <- expr.filt.matrices[[ds]][common.genes, ]
  test.genomic.args[nxt.indx] <- list(NULL)
  test.clinical.args[nxt.indx] <- list(NULL)
  test.common.drugs[[nxt.indx]] <- common.drugs.by.data.set[[ds]]
  test.drug.cols[[nxt.indx]] <- drug.name.col
  test.patient.cols[[nxt.indx]] <- patient.id.cols[[ds]]
  test.response.cols[[nxt.indx]] <- "dss.auc.ll4"
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

  nxt.indx <- length(test.set.names) + 1
  test.set.names[[nxt.indx]] <- ds
  test.dss.args[[nxt.indx]] <- fits[[ds]]
  test.expr.args[[nxt.indx]] <- expr.filt.matrices[[ds]][common.genes, ]
  test.genomic.args[nxt.indx] <- list(NULL)
  test.clinical.args[nxt.indx] <- list(NULL)
  test.common.drugs[[nxt.indx]] <- common.drugs.by.data.set[[ds]]
  test.drug.cols[[nxt.indx]] <- drug.name.col
  test.patient.cols[[nxt.indx]] <- patient.id.cols[[ds]]
  test.response.cols[[nxt.indx]] <- "dss.auc.ll4"
} ## if("ntap" %in% data.sets.to.analyze) 

ds <- "sanger"
if(ds %in% data.sets.to.analyze) {

  nxt.indx <- length(test.set.names) + 1
  test.set.names[[nxt.indx]] <- ds
  test.dss.args[[nxt.indx]] <- fits[[ds]]
  test.expr.args[[nxt.indx]] <- expr.filt.matrices[[ds]][common.genes, ]
  test.genomic.args[nxt.indx] <- list(NULL)
  test.clinical.args[nxt.indx] <- list(NULL)
  test.common.drugs[[nxt.indx]] <- common.drugs.by.data.set[[ds]]
  test.drug.cols[[nxt.indx]] <- drug.name.col
  test.patient.cols[[nxt.indx]] <- patient.id.cols[[ds]]
  test.response.cols[[nxt.indx]] <- "dss.auc.ll4"
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
     vec <- drc[num.samples.per.drug[i, train.drug.cols[[train.indx]]],]
     flag <- !is.na(vec)
     drug.samples <- intersect(samples, names(vec)[flag])
     num.samples.per.drug[i, col] <- length(unique(drug.samples))
   }
}

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
     vec <- drc[num.samples.per.drug[i, test.drug.cols[[test.indx]]],]
     flag <- !is.na(vec)
     drug.samples <- intersect(samples, names(vec)[flag])
     num.samples.per.drug[i, col] <- length(unique(drug.samples))
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

save.image(".Rdata.focs")

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
res.auc <- NULL
do.downsampling <- FALSE
if(do.downsampling) {
  res.auc <- train.and.test.crossproduct.with.downsampling(gene.sets = gene.sets, drug.name.tbl = drug.name.tbl, train.set.names = train.set.names, train.dss.args = train.dss.args, 
                                         train.expr.args = train.expr.args, train.genomic.args = train.genomic.args, train.clinical.args = train.clinical.args, 
                                         train.common.drugs = train.common.drugs, train.drugs = train.drugs, train.drug.cols = train.drug.cols, train.patient.cols = train.patient.cols, 
                                         train.response.cols = train.response.cols, test.set.names = test.set.names, test.dss.args = test.dss.args, test.expr.args = test.expr.args, 
                                         test.genomic.args = test.genomic.args, test.clinical.args = test.clinical.args, test.common.drugs = test.common.drugs,
                                         test.drug.cols = test.drug.cols, test.patient.cols = test.patient.cols, test.response.cols = test.response.cols, 
                                         seed = 1234, use.rf = FALSE, use.svm = FALSE, use.mean = FALSE, num.processes = num.processes, 
                                         num.train.samples.per.drug = num.train.samples.per.drug, num.test.samples.per.drug = num.test.samples.per.drug,
                                         num.iterations = 100, with.replacement = FALSE, return.fits = FALSE) 
} else { 
  ## Alternately, can do modeling without downsampling
  res.auc <- train.and.test.crossproduct(gene.sets = gene.sets, drug.name.tbl = drug.name.tbl, train.set.names = train.set.names, train.dss.args = train.dss.args, 
                                         train.expr.args = train.expr.args, train.genomic.args = train.genomic.args, train.clinical.args = train.clinical.args, 
                                         train.common.drugs = train.common.drugs, train.drugs = train.drugs, train.drug.cols = train.drug.cols, train.patient.cols = train.patient.cols, 
                                         train.response.cols = train.response.cols, test.set.names = test.set.names, test.dss.args = test.dss.args, test.expr.args = test.expr.args, 
                                         test.genomic.args = test.genomic.args, test.clinical.args = test.clinical.args, test.common.drugs = test.common.drugs,
                                         test.drug.cols = test.drug.cols, test.patient.cols = test.patient.cols, test.response.cols = test.response.cols, 
                                         seed = 1234, use.rf = FALSE, use.svm = FALSE, use.mean = FALSE, num.processes = num.processes)
}
## save.image(".Rdata")
cat("Done training and testing with AUC\n")

## save.image(".Rdata.focs")

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
                                         seed = 1234, use.rf = FALSE, use.svm = FALSE, use.mean = FALSE, num.processes = num.processes,
                                         num.train.samples.per.drug = num.train.samples.per.drug, num.test.samples.per.drug = num.test.samples.per.drug,
                                         num.iterations = 100, with.replacement = FALSE, return.fits = FALSE) 
  ## save.image(".Rdata")
  cat("Done training and testing with shifted AUC\n")

  save.image(".Rdata.focs")
}

## Sanitize the above results list into a table that has columns for 
## train.set, test.set, gene.set (i.e., "gene", "kegg", "biocarta", etc.), response (e.g., "AUC", "IC50"),
## metric ("pearson" or "spearman" correlation), val (the correlation), and pval (the associated p-value)
## Note that this table will have one row for each random (down) sample.
tbl.full <- ldply(train.set.names,
             .fun = function(train.set) {
                      ldply(test.set.names, 
                            .fun = function(test.set) {
                                     ldply(names(gene.sets),
                                           .fun = function(gene.set) {
                                                    tmp.auc <- extract.results(res.auc[["all.comparisons"]][[train.set]][[test.set]][[gene.set]])
                                                    tmp.auc$train.set <- train.set
                                                    tmp.auc$test.set <- test.set
                                                    tmp.auc$gene.set <- gene.set
                                                    tmp.auc$response <- "AUC"
                                                    if(do.shifted.auc) {
                                                      tmp.shifted.auc <- extract.results(res.shifted.auc[["all.comparisons"]][[train.set]][[test.set]][[gene.set]])
                                                      tmp.shifted.auc$train.set <- train.set
                                                      tmp.shifted.auc$test.set <- test.set
                                                      tmp.shifted.auc$gene.set <- gene.set
                                                      tmp.shifted.auc$response <- "shifted-AUC"
                                                      return(rbind(tmp.auc, tmp.shifted.auc))
                                                    }
                                                    return(tmp.auc)                         
                                                  }) 
                                   }) 
                    })

cat("Extracted results\n")
save.image(".Rdata.focs")

## Create a table that summarizes downsamplings by taking the entry corresponding to the median correlation.
if(do.downsampling) { 
  vars <- colnames(tbl.full)
  vars <- vars[!(vars %in% c("val", "pval"))]
  tbl <- ddply(tbl.full, .variables = vars, 
               .fun = function(df) {
                 if(all(is.na(df$val) | is.na(df$pval))) {
                   vec <- c(NA, NA)
                   names(vec) <- c("val", "pval")
                   return(vec)
                 }
                 df <- df[!is.na(df$val) & !is.na(df$pval),,drop=FALSE]
                 df <- df[order(as.numeric(df$val), decreasing=FALSE),,drop=FALSE]
                 indx <- max(1,floor(nrow(df)/2))
                 vec <- c(df$val[indx], df$pval[indx])
                 names(vec) <- c("val", "pval")
                 vec
               })
} else {
  tbl <- tbl.full
}
cat("Summarized results\n")

## For each train-test combination create a vector of correlations (one for each drug) of predicted vs actual in corr.lists
## There is a corresponding list of pvalues in pval.lists.
corr.lists <- list()
pval.lists <- list()
for(test.set.name in unique(tbl$test.set)) {
  corr.lists[[test.set.name]] <- list()
  pval.lists[[test.set.name]] <- list()
  for(mdl in unique(tbl$model)) {
    corr.lists[[test.set.name]][[mdl]] <- list()
    pval.lists[[test.set.name]][[mdl]] <- list()
    for(resp in unique(tbl$response)) {
      corr.lists[[test.set.name]][[mdl]][[resp]] <- list()
      pval.lists[[test.set.name]][[mdl]][[resp]] <- list()
      mets <- c("pearson", "spearman")
      for(met in mets[mets %in% tbl$metric]) {
        corr.lists[[test.set.name]][[mdl]][[resp]][[met]] <- list()
        pval.lists[[test.set.name]][[mdl]][[resp]][[met]] <- list()
        for(gs in unique(tbl$gene.set)) {
          corr.lists[[test.set.name]][[mdl]][[resp]][[met]][[gs]] <- list()
          pval.lists[[test.set.name]][[mdl]][[resp]][[met]][[gs]] <- list()
          t.sub <- subset(tbl, (test.set == test.set.name) & (model == mdl) & (response == resp) & (gene.set == gs) & (metric == met))
          na.drugs <- t.sub$test.drug[is.na(t.sub$val)]
          t.sub <- t.sub[!(t.sub$test.drug %in% na.drugs),]
          test.drugs <- unique(t.sub$test.drug)
          for(train.set.name in unique(tbl$train.set)) {
            tr.sub <- subset(t.sub, train.set == train.set.name) 
            vec <- tr.sub$val
            names(vec) <- tr.sub$test.drug
            pvec <- tr.sub$pval
            names(pvec) <- tr.sub$test.drug
            corr.lists[[test.set.name]][[mdl]][[resp]][[met]][[gs]][[train.set.name]] <- vec[test.drugs]
            pval.lists[[test.set.name]][[mdl]][[resp]][[met]][[gs]][[train.set.name]] <- pvec[test.drugs]
          }
        }
      }
    }
  }
}

## Make violin plots of correlations (model prediction of drug response vs actual drug response) faceted by
## training/test data set and type of feature (gene, biocart, hallmark, kegg, and reactome).
tbl$gene.set <- factor(tbl$gene.set, levels = c("gene", "biocarta", "hallmark", "kegg", "reactome"))
tbl$train.set <- factor(tbl$train.set, levels = c("ctrp", "ctrp.all", "ctrp.heme", "ctrp.aml", "ohsu.train", "ntap"))
for(resp in unique(tbl$response)) {
  for(met in c("pearson", "spearman")) {
    for(mdl in c("ridge", "lasso")) {
      sub <- subset(tbl, metric == met & model == mdl & response == resp)
      g <- ggplot(data = sub, aes(x = gene.set, y = val))
      g <- g + facet_grid(train.set ~ test.set)
      g <- g + geom_violin()
##    g <- g + geom_beeswarm()
      g <- g + geom_boxplot(width = 0.5)
      g <- g + geom_hline(yintercept = 0, linetype = "dashed")
      g <- g + ylab(paste0(capwords(met), " Correlation"))
      g <- g + ggtitle(paste0(resp, " ", capwords(met), " ", mdl))
      g <- g + ylim(c(-1,1))
      file <- paste0(file.prefix, "-", resp, "-", mdl, "-", met, "-vs-gene-sets.pdf")
      pdf(file)
      print(g)
      d <- dev.off()
##    glist[[length(glist)+1]] <- g
    }
  }
##  do.call("grid.arrange", glist)
}

## Make volcano plots of correlations (model prediction of drug response vs actual drug response) faceted by
## training/test data set and type of feature (gene, biocart, hallmark, kegg, and reactome).
for(resp in unique(tbl$response)) {
  for(met in c("pearson", "spearman")) {
    for(mdl in c("ridge", "lasso")) {
      file <- paste0(file.prefix, "-", resp, "-", mdl, "-", met, "-vs-gene-sets-volcano.pdf")
      pdf(file)
      for(gs in unique(tbl$gene.set)) {
        sub <- subset(tbl, metric == met & model == mdl & response == resp & gene.set == gs)
        tmp <- sub
        pvalue <- rep(">= 0.05", length(tmp$pval))
        pvalue[!is.na(tmp$pval) & (tmp$pval < 0.05)] <- "< 0.05"
        tmp$pvalue <- pvalue
        tmp$pval <- -log10(tmp$pval)
        ## Training on OHSU and testing on OHSU leads to very low pvalues--artificially
        ## adjust them so that they do not set the scale for the other comparisons.
##        flag <- sub$train.set == "ohsu.train" & sub$test.set == "ohsu.test"
##        if(any(flag)) { tmp$pval[flag] <- tmp$pval[flag] / max(tmp$pval[flag], na.rm=TRUE) }
##        if(any(flag)) { tmp$pval[flag] <- tmp$pval[flag] / 10 }
        g <- ggplot(data = tmp, aes(x = val, y = pval, colour = pvalue))
        g <- g + facet_grid(train.set ~ test.set)
        g <- g + geom_point()
        g <- g + xlab(paste0(capwords(met), " Correlation"))
        g <- g + geom_hline(yintercept = -log10(0.05), linetype = "dashed")
        g <- g + geom_vline(xintercept = 0, linetype = "dashed")
        g <- g + ylab("-Log10 p-value")
        g <- g + ggtitle(paste0(resp, " ", capwords(met), " ", mdl, " ", gs))
        print(g)
      }
      d <- dev.off()
##    glist[[length(glist)+1]] <- g
    }
  }
##  do.call("grid.arrange", glist)
}

cat("Created volcano plots\n")

save.image(".Rdata.focs")
cat("Successfully exiting\n")
q()

panel.diff <- function (x, y, col = par("col"), bg = NA, pch = par("pch"), 
                              cex = 1, col.regres = "red", ...) 
{ 
  ok <- is.finite(x) & is.finite(y) 
print(x-y)
vec <- x - y
  points(1:length(vec), vec, pch = pch, col = col, bg = bg, cex = cex) 
##  abline(h = 0, col = col.regres, ...) 
} 

save.image(".Rdata.focs")

if(length(unique(tbl$train.set)) > 1) { 
## Create pairwise correlation plots
for(test.set.name in unique(tbl$test.set)) {
  for(mdl in c("lasso", "ridge")) {
    for(resp in unique(tbl$response)) {
      pdf(paste(file.prefix, test.set.name, mdl, resp, "corr.pdf", sep="-"))
      mets <- c("pearson", "spearman")
      for(met in mets[mets %in% tbl$metric]) {
        for(gs in unique(tbl$gene.set)) {
          mat <- do.call("cbind", corr.lists[[test.set.name]][[mdl]][[resp]][[met]][[gs]])
          cols <- NULL
          if("Class.explained" %in% colnames(drug.name.tbl)) {
            cls <- drug.name.tbl$Class.explained
            names(cls) <- drug.name.tbl[, drug.name.col]
            cols <- factor(cls[names(corr.lists[[test.set.name]][[mdl]][[resp]][[met]][[gs]][[1]])])
          }
##        pdf(paste(test.set.name, mdl, resp, gs, "corr.pdf", sep="-"))
          pairs(mat, lower.panel=panel.regression, upper.panel=panel.cor, main = paste(test.set.name, mdl, resp, met, gs, sep=" "), col = cols, oma=c(4,4,6,12))
          par(xpd=TRUE)
          legend(0.85, 0.7, legend = levels(cols), fill = c(1:length(levels(cols))), cex = 0.7)
##        d <- dev.off()
        }
      }
      d <- dev.off()
    }
  }
}

for(test.set.name in unique(tbl$test.set)) {
  for(mdl in c("lasso", "ridge")) {
    for(resp in unique(tbl$response)) {
      pdf(paste(file.prefix, test.set.name, mdl, resp, "diff.pdf", sep="-"))
      mets <- c("pearson", "spearman")
      for(met in mets[mets %in% tbl$metric]) {
        for(gs in unique(tbl$gene.set)) {
          mat <- do.call("cbind", corr.lists[[test.set.name]][[mdl]][[resp]][[met]][[gs]])
          cols <- NULL
          if("Class.explained" %in% colnames(drug.name.tbl)) {
            cls <- drug.name.tbl$Class.explained
            names(cls) <- drug.name.tbl[, drug.name.col]
            cols <- factor(cls[names(corr.lists[[test.set.name]][[mdl]][[resp]][[met]][[gs]][[1]])])
          }
          len <- length(corr.lists[[test.set.name]][[mdl]][[resp]][[met]][[gs]])
          glist <- list()
          for(i in 1:len) {
            for(j in 1:len) {
              g <- NULL
              if(i != j) {
                x <- y <- NULL
                if(i > j) {
                  x <- corr.lists[[test.set.name]][[mdl]][[resp]][[met]][[gs]][[i]]
                  y <- corr.lists[[test.set.name]][[mdl]][[resp]][[met]][[gs]][[j]]
                  df <- data.frame(x = x, y = y, col = cols)
                } else {
                  x <- corr.lists[[test.set.name]][[mdl]][[resp]][[met]][[gs]][[j]]
                  y <- corr.lists[[test.set.name]][[mdl]][[resp]][[met]][[gs]][[i]]
                  df <- data.frame(x = 1:length(x), y = x - y, col = cols)
                }
                if(i < j) { g <- ggplot(df, aes(x = x, y = y, colour = col)) } else { g <- ggplot(df, aes(x = x, y = y)) }
                g <- g + geom_point()
                g <- g + theme(legend.position="none")
                if(i < j) {
                  g <- g + geom_hline(yintercept = 0)
                } else {
                  g <- g + geom_smooth(method='lm')
                  lims <- c(min(min(df$x, na.rm=TRUE), min(df$y, na.rm=TRUE)), max(max(df$x, na.rm=TRUE), max(df$y, na.rm=TRUE)))
                  g <- g + ggtitle(lm_eqn_bq(df))
                  g <- g + geom_abline(intercept = 0, slope = 1, linetype="dashed")
                  ## g <- g + geom_text(x = lims[1] + 0.2, y = lims[1] + 0.1, label = lm_eqn(df), parse=TRUE)
                }
              } else {
                g <- ggplot() + annotation_custom(textGrob(names(corr.lists[[test.set.name]][[mdl]][[resp]][[met]][[gs]])[i]),
                                 xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) 
              }
              g <- g + xlab("") + ylab("")
              if(i < j) {
                g <- g + ylab(paste0(names(corr.lists[[test.set.name]][[mdl]][[resp]][[met]][[gs]])[j], " - ", names(corr.lists[[test.set.name]][[mdl]][[resp]][[met]][[gs]])[i]))
              }
              if(i > j) {
                g <- g + xlab(names(corr.lists[[test.set.name]][[mdl]][[resp]][[met]][[gs]])[i])
                g <- g + ylab(names(corr.lists[[test.set.name]][[mdl]][[resp]][[met]][[gs]])[j])
              }
  
              glist[[length(glist)+1]] <- g
            }
          }
          do.call("grid.arrange", c(glist, "top" = paste(test.set.name, mdl, resp, met, gs, sep=" ")))
        }
      }
      d <- dev.off()
    }
  }
}

} ## if(length(unique(tbl$train.set)) > 1) 

stop("successfully completed")
