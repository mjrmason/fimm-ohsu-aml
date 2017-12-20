rm(list = ls())

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
source("../common/utils.R")

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
source("fimm-ohsu-setup-120817.R")
file.prefix <- "ohsu-bootstrap-120817"
rdata.file <- ".Rdata.ohsu.bootstrap.120817"

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

## Drop any fits that are flagged as outliers
fct.postfixes <- list("LL.4" = ".ll4", "L.4" = ".l4")
fct <- "LL.4"
for(ds in names(fits)) {
    exclude.cols <- colnames(fits[[ds]])[grepl(pattern=paste0("exclude", fct.postfixes[[fct]]), x=colnames(fits[[ds]]))]
    any.excluded <- unlist(apply(fits[[ds]][, exclude.cols], 1, function(row) any(row)))
    orig.nrow <- nrow(fits[[ds]])
    fits[[ds]] <- fits[[ds]][!any.excluded, ]
    new.nrow <- nrow(fits[[ds]])
    cat(paste0(ds, ": filtered ", orig.nrow - new.nrow, " of the original ", orig.nrow, " fits.\n"))
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
for(ds in names(fits)) {
  fits[[ds]] <- merge(fits[[ds]], drug.map, by.x = data.set.drug.id.cols[[ds]], by.y = drug.map.drug.id.cols[[ds]], all = FALSE)
}

common.drugs <- Reduce("intersect", lapply(fits, function(fit) unique(fit[, drug.name.col])))
for(ds in names(fits)) {
  flag <- fits[[ds]][, drug.name.col] %in% common.drugs
  fits[[ds]] <- fits[[ds]][flag, ]
}

## Limit the drug map to common drugs (which, by design of the drug map file, should already be the case)
drug.name.tbl <- drug.map[drug.map[, drug.name.col] %in% common.drugs,]

common.drugs.by.data.set <- list()
for(ds in names(fits)) {
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
## common.genes <- intersect(common.genes, gene.subset)

## Restrict the matrices to those genes in common across data sets
for(ds in names(fits)) {
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
expr.filt.matrices <- exprs
if(FALSE) {
expr.filt.matrices <- list()
for(ds in names(exprs)) {
  iqrs <- unlist(apply(exprs[[ds]], 1, IQR))
  expr.filt.matrices[[ds]] <- exprs[[ds]][iqrs > median(iqrs[iqrs > 0]),]
}
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
  train.patient.cols[[nxt.indx]] <- expr.patient.id.cols[[ds]]
  train.response.cols[[nxt.indx]] <- "dss.auc.ll4"

  set.seed(1234)
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
  train.patient.cols[[nxt.indx]] <- expr.patient.id.cols[[ds]]
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
  train.patient.cols[[nxt.indx]] <- expr.patient.id.cols[[ds]]
  train.response.cols[[nxt.indx]] <- "dss.auc.ll4"


} ## if("ohsu" %in% data.sets.to.analyze) 

ds <- "fimm"
if(ds %in% data.sets.to.analyze) {

  if(FALSE) {
  nxt.indx <- length(test.set.names) + 1
  test.set.names[[nxt.indx]] <- ds
  test.dss.args[[nxt.indx]] <- fits[[ds]]
  test.expr.args[[nxt.indx]] <- expr.filt.matrices[[ds]][common.genes, ]
  test.genomic.args[nxt.indx] <- list(NULL)
  test.clinical.args[nxt.indx] <- list(NULL)
  test.common.drugs[[nxt.indx]] <- common.drugs.by.data.set[[ds]]
  test.drug.cols[[nxt.indx]] <- drug.name.col
  test.patient.cols[[nxt.indx]] <- expr.patient.id.cols[[ds]]
  test.response.cols[[nxt.indx]] <- "dss.auc.ll4"

  nxt.indx <- length(train.set.names) + 1
  train.set.names[[nxt.indx]] <- ds
  train.dss.args[[nxt.indx]] <- fits[[ds]]
  train.expr.args[[nxt.indx]] <- expr.filt.matrices[[ds]][common.genes, ]
  train.genomic.args[nxt.indx] <- list(NULL)
  train.clinical.args[nxt.indx] <- list(NULL)
  train.common.drugs[[nxt.indx]] <- common.drugs.by.data.set[[ds]]
  train.drugs[[nxt.indx]] <- common.drugs.by.data.set[[ds]]
  train.drug.cols[[nxt.indx]] <- drug.name.col
  train.patient.cols[[nxt.indx]] <- expr.patient.id.cols[[ds]]
  train.response.cols[[nxt.indx]] <- "dss.auc.ll4"
  }

  if(FALSE) {
  set.seed(1234)
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
  train.patient.cols[[nxt.indx]] <- expr.patient.id.cols[[ds]]
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
  train.patient.cols[[nxt.indx]] <- expr.patient.id.cols[[ds]]
  train.response.cols[[nxt.indx]] <- "dss.auc.ll4"
  }

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
  train.patient.cols[[nxt.indx]] <- expr.patient.id.cols[[ds]]
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
  train.patient.cols[[nxt.indx]] <- expr.patient.id.cols[[ds]]
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
  train.patient.cols[[nxt.indx]] <- expr.patient.id.cols[[ds]]
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
glmnet.alphas <- seq(from = 0, to = 1, by = 0.1)
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
  seeds <- c(1234)
  for(i in 1:length(seeds)) {
    seed <- seeds[i]
    cat(paste0("Training for seed ", seed, "\n"))
    res.aucs[[i]] <- train.and.test.crossproduct(gene.sets = gene.sets, drug.name.tbl = drug.name.tbl, train.set.names = train.set.names, train.dss.args = train.dss.args, 
                                           train.expr.args = train.expr.args, train.genomic.args = train.genomic.args, train.clinical.args = train.clinical.args, 
                                           train.common.drugs = train.common.drugs, train.drugs = train.drugs, train.drug.cols = train.drug.cols, train.patient.cols = train.patient.cols, 
                                           train.response.cols = train.response.cols, test.set.names = test.set.names, test.dss.args = test.dss.args, test.expr.args = test.expr.args, 
                                           test.genomic.args = test.genomic.args, test.clinical.args = test.clinical.args, test.common.drugs = test.common.drugs,
                                           test.drug.cols = test.drug.cols, test.patient.cols = test.patient.cols, test.response.cols = test.response.cols, 
                                           seed = seed, use.rf = use.rf, use.svm = FALSE, use.mean = FALSE, num.processes = num.processes, alphas = glmnet.alphas,
                                           keep.forest = FALSE)
  }
}
cat("Done training and testing with AUC\n")

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
                                                    tmp.auc <- extract.results(res.aucs[[1]][["all.comparisons"]][[train.set]][[test.set]][[gene.set]])
                                                    tmp.auc$train.set <- train.set
                                                    tmp.auc$test.set <- test.set
                                                    tmp.auc$gene.set <- gene.set
                                                    tmp.auc$response <- "AUC"

                                                    tmp.features <- extract.features(res.aucs[[1]][["all.fits"]][[train.set]][[gene.set]], s = "lambda.min")
                                                    tmp.auc <- merge(tmp.auc, tmp.features, all.x = TRUE, all.y = FALSE)

                                                    return(tmp.auc)                         
                                                  }) 
                                   }) 
                    })


trn.tbl <- extract.table(res.aucs, train.set.names, gene.sets, s = "lambda.min") 
trn.tbl.1se <- extract.table(res.aucs, train.set.names, gene.sets, s = "lambda.1se") 

trn.tbl.glmnet <- subset(trn.tbl, model %in% c("glmnet", "lasso", "ridge"))
trn.tbl <- subset(trn.tbl, model != "glmnet")

trn.tbl.lasso <- subset(trn.tbl, model == "lasso")
trn.tbl.ridge <- subset(trn.tbl, model == "ridge")
trn.tbl.rf <- subset(trn.tbl, model == "rf")

list.sizes <- unique(c(seq(from=10,to=2000,by=10), seq(from=2000, to=10000, by=100)))
list.sizes <- c(seq(from=10,to=1000, by=10), 10,100,1000,2000,5000,7000,7500,7600,length(common.genes)-100,
                length(common.genes)-10,length(common.genes)-1, length(common.genes),length(common.genes)+1,8000) 
list.sizes <- list.sizes[list.sizes <= length(common.genes)]
list.sizes <- sort(list.sizes)
names(list.sizes) <- list.sizes

tbl.rf.cross <- subset(tbl.full, model == "rf" & train.set == "ohsu" & test.set == "fimm" & metric == "pearson")
colnames(tbl.rf.cross)[colnames(tbl.rf.cross) == "pval"] <- "corr.pval"
colnames(tbl.rf.cross)[colnames(tbl.rf.cross) == "val"] <- "corr"

num.rf.sigs <- NULL
if(use.rf) {

  num.rf.sigs <- ldply(list.sizes, .parallel = TRUE,
                       .fun = function(sz) {
                                tmp <- find.overlap.of.sparse.predictors(trn.tbl.rf, universe = common.genes, transform = "none", top = sz, verbose = FALSE)
                                tmp <- tmp[, c("FIMM_DRUG_NAME", "train.set.1", "train.set.2", "n.1", "n.2", "n.both", "n.neither", "pval", "top", "gene.set.1")]
                                tmp <- subset(tmp, (train.set.1 == "ohsu") & (train.set.2 == "fimm"))
                                tmp
                              })

  tmp <- merge(num.rf.sigs, tbl.rf.cross, by.x = c("FIMM_DRUG_NAME", "gene.set.1"), by.y = c("test.drug", "gene.set"))
  tmp$top <- as.numeric(tmp$top)
  tmp$n.both <- as.numeric(tmp$n.both)
  tmp$pval <- as.numeric(tmp$pval)
  pdf(paste0(file.prefix, "-rf-top-feature-overlap.pdf"))
  d_ply(tmp, .variables = c("FIMM_DRUG_NAME"),
        .fun = function(df) {
                 flag <- !is.na(df$n.both)
                 corr <- df$corr[1]; corr.pval <- df$corr.pval[1]
##                 if((corr < 0) || (corr.pval > 0.05)) { return() }
##                 if(!any(df$pval < 0.05)) { return() }
##                 plot(df$top[flag], df$n.both[flag], main = paste0(df$FIMM_DRUG_NAME[1], ": corr = ", signif(df$corr[1], digits=2), " pval = ", 
##                      signif(df$corr.pval[1], digits=2)))
                 g1 <- ggplot(data = df[flag,]) + geom_point(aes(x = top, y = n.both)) + xlab("Num Top Features") + ylab("Num Overlapping")
                 g1 <- g1 + ggtitle(paste0(df$FIMM_DRUG_NAME[1], ": corr = ", signif(df$corr[1], digits=2), " pval = ", signif(df$corr.pval[1], digits=2)))
                 g1 <- g1 + geom_abline(intercept = 0, slope = 1)
            
                 g2 <- ggplot(data = df[flag,]) + geom_point(aes(x = top, y = -log10(pval))) + xlab("Num Top Features") + ylab("Overlap -log10 p-value")
                 g2 <- g2 + geom_hline(yintercept = -log10(0.05))
                 grid.arrange(g1, g2)
        })
  d <- dev.off()
}


cat("Done.\n")

save.image(rdata.file)

q(status=0)


cat("Making scatter plots of FIMM vs OHSU betas RF\n")
## Make a scatter plot of FIMM vs OHSU _ridge_ coefficients
pdf(paste0(file.prefix, "-fimm-vs-ohsu-", "-rf-feature-overlap-correlation.pdf"))
ds1 <- "ohsu"
ds2 <- "fimm"
drugs <- train.drugs[[1]]

##tbl <- read.table(file="fimm-ohsu-validate-trained-models-112217-one-ohsu-AUC-ridge-pearson-vs-gene-sets-all-pos-and-sig.xls", sep="\t", header=TRUE, as.is=TRUE)

## drugs <- c("AZD1152-HQPA", "Cabozantinib", "Cytarabine", "Erlotinib", "Foretinib", "Gefitinib", "Pictilisib", "Ruboxistaurin", "Sorafenib", "Sunitinib")
## drugs <- c("Cabozantinib", "Ruboxistaurin", "Sorafenib")


for(train.drug in drugs) {
  drug1 <- train.drug
  drug2 <- train.drug
  scatterplot(ds1, drug1, ds2, drug2, save.to.pdf = FALSE, model = "rf", alpha = NA)
}
d <- dev.off()
cat("Done making scatter plots of FIMM vs OHSU betas RF\n")

validate.prefix <- "fimm-ohsu-validate-trained-models-112217-one-ohsu"
tbl <- read.table(file=paste0(validate.prefix, "-AUC-ridge-pearson-vs-gene-sets-all.xls"), sep="\t", header=TRUE, as.is=TRUE)
pos.flag <- !is.na(tbl$fimm.val) & (tbl$fimm.val > 0)
sig.flag <- !is.na(tbl$fimm.pval) & (tbl$fimm.pval < 0.05)
drugs <- tbl$train.drug[pos.flag]
sig.drugs <- tbl$train.drug[sig.flag & pos.flag]
pos.drugs <- tbl$train.drug[pos.flag]

names(drugs) <- drugs
fimm.ohsu.common.markers <- llply(drugs, .fun = function(drug) common.markers("fimm", drug, "ohsu", drug, alpha = 0, model = "glmnet", top = 25))
## fimm.ohsu.common.markers <- llply(drugs, .fun = function(drug) common.markers("fimm", drug, "ohsu", drug, alpha = 1, model = "glmnet", top = 20))
tmp <- ldply(fimm.ohsu.common.markers, .fun = function(lst) data.frame(gene = lst))
colnames(tmp) <- c("drug", "gene")
library(mygene)
baz <- queryMany(tmp$gene, scopes="ensembl.gene", fields=c("symbol"), species="human")
fimm.ohsu.common.markers <- merge(tmp, unique(as.data.frame(baz[,c("query","symbol")])), by.x = "gene", by.y = "query")
mat <- spread(fimm.ohsu.common.markers[, c("drug", "gene", "symbol")], key = drug, value = gene)
flag <- grepl(pattern="^AC", x = rownames(mat)) | grepl(pattern="^LINC", x = rownames(mat)) | grepl(pattern="^LOC", x = rownames(mat))
rownames(mat) <- mat[, 1]
mat <- as.matrix(mat[, -1])
mat[!is.na(mat)] <- 1
mat[is.na(mat)] <- 0
class(mat) <- "numeric"
library(gplots)

pdf(paste0(file.prefix, "-fimm-vs-ohsu-", "-ridge-pearson-sig-drugs-feature-heatmap.pdf"))
col.flag <- colnames(mat) %in% sig.drugs
tmp <- mat[, col.flag]
flag <- rowSums(tmp) > 0
heatmap(tmp[flag,], Rowv = NA, Colv = NA, scale = "none", main = "Significant (p < 0.05), positively correlated drugs")
d <- dev.off()

pdf(paste0(file.prefix, "-fimm-vs-ohsu-", "-ridge-pearson-sig-drugs-feature-heatmap.pdf"))
row.flag <- rowSums(mat) > 1
tmp <- mat[row.flag, , drop=FALSE]
flag <- colSums(tmp) > 0
heatmap(tmp[, flag, drop=FALSE], Rowv = NA, Colv = NA, scale = "none", main = "Positively correlated drugs; genes occurring at least twice")
d <- dev.off()

cat("Making scatter plots of FIMM vs OHSU betas glmnet\n")
## Make a scatter plot of FIMM vs OHSU _ridge_ coefficients
pdf(paste0(file.prefix, "-fimm-vs-ohsu-", "-ridge-feature-overlap-correlation.pdf"))
ds1 <- "ohsu"
ds2 <- "fimm"
for(train.drug in drugs) {
  drug1 <- train.drug
  drug2 <- train.drug
  scatterplot(ds1, drug1, ds2, drug2, save.to.pdf = FALSE, model = "glmnet", alpha = 0, main = train.drug)
}
d <- dev.off()
cat("Done making scatter plots of FIMM vs OHSU betas ridge\n")


save.image(rdata.file)

source("continue.R")


cat("Done saving\n")
q(status = 0)

foo <- extract.fit(res.aucs[[1]][["all.fits"]][["ohsu"]][["gene"]], train.drug = "PD184352", alpha = 0, model = "glmnet")
## bar <- extract.fit(res.aucs[[1]][["all.fits"]][["ohsu"]][["gene"]], train.drug = "Alvocidib", alpha = 0, model = "glmnet")
bar <- extract.fit(res.aucs[[1]][["all.fits"]][["ohsu"]][["gene"]], train.drug = "Sunitinib", alpha = 0, model = "glmnet")
## bar <- extract.fit(res.aucs[[1]][["all.fits"]][["fimm"]][["gene"]], train.drug = "PD184352", alpha = 0, model = "glmnet")

ds1 <- "ohsu"
gene1 <- "PD184352"
ds2 <- "fimm"
gene2 <- "PD184352"

scatterplot(ds1, gene1, ds2, gene2, save.to.pdf = TRUE)

ds1 <- "ohsu"
gene1 <- "PD184352"
ds2 <- "ohsu"
gene2 <- "Sunitinib"

scatterplot(ds1, gene1, ds2, gene2, save.to.pdf = TRUE)


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


trn.tbl.1se.lasso <- subset(trn.tbl.1se, model == "lasso")
trn.tbl.1se.ridge <- subset(trn.tbl.1se, model == "ridge")
trn.tbl.1se.rf <- subset(trn.tbl.1se, model == "rf")

## save.image(".Rdata.overlap.fimm.ohsu")

cat("Done subsetting table by model\n")

cat("Summarizing ridge results\n")
all.ridge <- find.overlap.of.dense.predictors(trn.tbl.ridge, method = "spearman")
all.ridge$pval <- as.numeric(all.ridge$pval)
cat("Done summarizing ridge results\n")

if(use.rf) {
  cat("Summarizing rf results\n")
  ## all.rf <- find.overlap.of.dense.predictors(trn.tbl.rf, method = "spearman")
  all.rf <- find.overlap.of.sparse.predictors(trn.tbl.rf, transform = "none", top = 100)
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

all.1se <- find.overlap.of.sparse.predictors(trn.tbl.1se.lasso, transform = "drop.zero", top = 0)
all.1se$pval <- as.numeric(all.1se$pval)
cat("Done summarizing lasso results\n")

## save.image(".Rdata.overlap.fimm.ohsu")

lst <- list("lasso" = all, "lasso.1se" = all.1se)
if(use.rf) {
  lst <- list("lasso" = all, "lasso.1se" = all.1se, "rf" = all.rf)
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
  pdf(paste0(file.prefix, "-ohsu-vs-fimm-", nm , "-feature-overlap.pdf"))
  plot.table(tbl, main = "OHSU (train) vs FIMM (train)")
  d <- dev.off()

  if(grepl(x=nm, pattern="lasso")) {  
    s <- subset(all.cross, n.1 > 0 & n.2 > 0)
    tmp <- s[,c("FIMM_DRUG_NAME", "train.set.1", "n.1", "train.set.2", "n.2", "gene.set.1", "n.both")]
    write.table(file = paste0(file.prefix, "-ohsu-vs-fimm-", nm , "-feature-overlap.tsv"), tmp, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  }

  tbl <- as.data.frame(table(ohsu.sig$gene.set.1))
  colnames(tbl) <- c("Feature Set", "Drugs w/ Overlapping Features")
  pdf(paste0(file.prefix, "-ohsu1-vs-ohsu2-", nm , "-feature-overlap.pdf"))
  plot.table(tbl, main = "OHSU 1 (train) vs OHSU 2 (train)")
  d <- dev.off()

  tbl <- as.data.frame(table(fimm.sig$gene.set.1))
  colnames(tbl) <- c("Feature Set", "Drugs w/ Overlapping Features")
  pdf(paste0(file.prefix, "-fimm1-vs-fimm2-", nm , "-feature-overlap.pdf"))
  plot.table(tbl, main = "FIMM 1 (train) vs FIMM 2 (train)")
  d <- dev.off()
}

save.image(rdata.file)
stop("successfully completed")
