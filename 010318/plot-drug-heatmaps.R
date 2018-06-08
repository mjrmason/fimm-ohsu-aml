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
suppressPackageStartupMessages(library(mygene))

suppressPackageStartupMessages(library("GSVA"))
# use GSA to read in a gmt file
suppressPackageStartupMessages(library("GSA"))

library(pcaMethods)
library(rJava)
library(rcdk)
library(fingerprint)
library(plyr)
library(tidyverse)

library(dynamicTreeCut)
library(MCL)
library(clusteval)

library(WeightedCluster)
library(igraph)


## Log in to Synapse
synapseLogin()

## Bring in common code
source("../common/data-preprocessing.R")
source("../common/models.R")
source("../common/plotting.R")
source("../common/utils.R")
source("heatmap-utils.R")

## Register parallelization
num.cores <- detectCores()
if(!is.na(num.cores) && (num.cores > 1)) {
  suppressPackageStartupMessages(library("doMC"))
  cat(paste("Registering ", num.cores-1, " cores.\n", sep=""))
  registerDoMC(cores=(num.cores-1))
}
num.processes <- num.cores - 1

## Begin setup (this will need to be changed as the data sets in the intersection change)

source("fimm-ohsu-setup-120817.R")

if(TRUE) {
## Overwrite the expression synIds to exclude outliers (as before), but _not_ to restrict to highly variable/cancer genes
## Batch-corrected log2 cpm (or rma) expression files (NB: these exclude 2 outliers in OHSU)
data.set.expr.synIds <- list("ohsu" = "syn11362256", "fimm" = "syn11362257")
## names(data.set.expr.synIds) <- data.sets
expr.folder.synIds <- list("ohsu" = "syn11361089", "fimm" = "syn11361089")
data.set.expr.files <- list("ohsu" = "ohsu-expr-fimm-ohsu-outlier1-combat.tsv", "fimm" = "fimm-expr-fimm-ohsu-outlier1-combat.tsv")
}

rdata.file <- ".Rdata.model.mean.response.120817"

## Over-write some of the variables
## These fit files are assumed to be for AUCs/DSS values calculated across shared concentration ranges across data sets.
## These files may be created by a script such as
## calculate-dss-fimm-ohsu.R
## data.set.dss.fit.synIds <- c("syn11362885", "syn11362891")
## names(data.set.dss.fit.synIds) <- data.sets
## data.set.dss.fit.files <- c("ohsu-fimm-ohsu-outlier1-dss-t0.tsv", "fimm-fimm-ohsu-outlier1-dss-t0.tsv")

## A string that will be included in any output files.
file.prefix <- "fimm-ohsu-validate-trained-models-112217-one-ohsu-no-filtering-all-mean-response"

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
for(ds in names(fits)) {
  fits[[ds]] <- merge(fits[[ds]], drug.map, by.x = data.set.drug.id.cols[[ds]], by.y = drug.map.drug.id.cols[[ds]], all = FALSE)
}

common.drugs <- Reduce("intersect", lapply(fits, function(fit) unique(fit[, drug.name.col])))
for(i in 1:length(fits)) {
  ds <- names(fits)[i]
  flag <- fits[[ds]][, drug.name.col] %in% common.drugs
  fits[[ds]] <- fits[[ds]][flag, ]
}

cat(paste0("Number of drugs in common: ", length(unique(common.drugs)), "\n"))

## Drop any fits that are flagged as outliers
if(TRUE) {
fct.postfixes <- list("LL.4" = ".ll4", "L.4" = ".l4")
fct <- "LL.4"
for(ds in names(fits)) {
    exclude.cols <- colnames(fits[[ds]])[grepl(pattern=paste0("exclude", fct.postfixes[[fct]]), x=colnames(fits[[ds]]))]
    any.excluded <- unlist(apply(fits[[ds]][, exclude.cols], 1, function(row) any(row)))
    orig.nrow <- nrow(fits[[ds]])
    fits[[ds]] <- fits[[ds]][!any.excluded, ]
    new.nrow <- nrow(fits[[ds]])
    cat(paste0(ds, ": filtered ", orig.nrow - new.nrow, " of the original ", orig.nrow, " fits, leaving ", new.nrow, " fits.\n"))
}
}

## Limit the drug map to common drugs (which, by design of the drug map file, should already be the case)
drug.name.tbl <- drug.map[drug.map[, drug.name.col] %in% common.drugs,]

file <- paste0(file.prefix, "-drug-name-tbl.tsv")
write.table(file = file, drug.name.tbl, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

cat("Breakdown of drug classes for common drugs:\n")
print(table(drug.name.tbl[, "Class.explained"]))

common.drugs.by.data.set <- list()
for(i in 1:length(fits)) {
  ds <- names(fits)[i]
  common.drugs.by.data.set[[ds]] <- drug.name.tbl[, drug.name.col]
}

## Load in the expr files
cat("Reading expression data sets\n")
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
if(FALSE) {
use.subset <- TRUE
  genes <- read.table("../common-resources/Geeleher-genes.tsv", sep="\t", header=TRUE, stringsAsFactors = FALSE)$symbol
  gene.subset <- genes[grepl(genes, pattern="ENSG")]
  genes.to.translate <- genes[!grepl(genes, pattern="ENSG")]
  trns <- symbols.to.ensg.mapping(genes.to.translate)
  gene.subset <- unique(c(gene.subset, trns$ensg))
}
if(FALSE) {
  genes <- read.table("../common-resources/gene_names.tsv", sep="\t", header=FALSE, stringsAsFactors = FALSE)$V1
  gene.subset <- c(gene.subset, genes[grepl(genes, pattern="ENSG")])
  genes.to.translate <- genes[!grepl(genes, pattern="ENSG")]
  trns <- symbols.to.ensg.mapping(genes.to.translate)
  gene.subset <- unique(c(gene.subset, trns$ensg))
}

if(use.subset) {
  common.genes <- intersect(common.genes, gene.subset)
}

cat(paste0("Num expression features: ", length(unique(common.genes)), "\n"))

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

all.gene.sets <- c(hallmark.gene.sets, kegg.gene.sets, biocarta.gene.sets, reactome.gene.sets)
ensg.to.gene.sets <- invert.list(all.gene.sets)

## BEGIN modeling setup

## Do analysis at pathway- (hallmark, biocarta, kegg, and reactome) as well as at gene-level
gene.sets <- list("gene" = "gene", "hallmark" = hallmark.gene.sets, "biocarta" = biocarta.gene.sets, "kegg" = kegg.gene.sets, "reactome" = reactome.gene.sets)
## Do analysis only at gene-level
## gene.sets <- list("gene" = "gene")

## Filter the gene expression data to restrict to highly-expressed genes
## We have already filtered before getting here.
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
  set.seed(1234)

  ## Use 60% of samples (stratified by drug) for training and 40% for test.
  train.index <- createDataPartition(fits[[ds]][, data.set.drug.id.cols[[ds]]], p = .6, list = FALSE)
  ohsu.dss.orig.train <- fits[[ds]][train.index, ]
  ohsu.dss.orig.test <- fits[[ds]][-train.index, ]

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

  nxt.indx <- length(test.set.names) + 1
  test.set.names[[nxt.indx]] <- "ohsu"
  test.dss.args[[nxt.indx]] <- fits[[ds]]
  test.expr.args[[nxt.indx]] <- expr.filt.matrices[[ds]][common.genes, ] 
  test.genomic.args[nxt.indx] <- list(NULL)
  test.clinical.args[nxt.indx] <- list(NULL)
  test.common.drugs[[nxt.indx]] <- common.drugs.by.data.set[[ds]]
  test.drug.cols[[nxt.indx]] <- drug.name.col
  test.patient.cols[[nxt.indx]] <- expr.patient.id.cols[[ds]]
  test.response.cols[[nxt.indx]] <- "dss.auc.ll4"

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

  if(FALSE) {
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
  }

  nxt.indx <- length(test.set.names) + 1
  test.set.names[[nxt.indx]] <- "ohsu.set2"
  test.dss.args[[nxt.indx]] <- fits.set2
  test.expr.args[[nxt.indx]] <- expr.filt.matrices[[ds]][common.genes, ] 
  test.genomic.args[nxt.indx] <- list(NULL)
  test.clinical.args[nxt.indx] <- list(NULL)
  test.common.drugs[[nxt.indx]] <- common.drugs.by.data.set[[ds]]
  test.drug.cols[[nxt.indx]] <- drug.name.col
  test.patient.cols[[nxt.indx]] <- expr.patient.id.cols[[ds]]
  test.response.cols[[nxt.indx]] <- "dss.auc.ll4"

  nxt.indx <- length(test.set.names) + 1
  test.set.names[[nxt.indx]] <- "ohsu.set1"
  test.dss.args[[nxt.indx]] <- fits.set1
  test.expr.args[[nxt.indx]] <- expr.filt.matrices[[ds]][common.genes, ] 
  test.genomic.args[nxt.indx] <- list(NULL)
  test.clinical.args[nxt.indx] <- list(NULL)
  test.common.drugs[[nxt.indx]] <- common.drugs.by.data.set[[ds]]
  test.drug.cols[[nxt.indx]] <- drug.name.col
  test.patient.cols[[nxt.indx]] <- expr.patient.id.cols[[ds]]
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

  if(TRUE) {
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

  if(FALSE) {
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

  nxt.indx <- length(test.set.names) + 1
  test.set.names[[nxt.indx]] <- "fimm.set2"
  test.dss.args[[nxt.indx]] <- fits.set2
  test.expr.args[[nxt.indx]] <- expr.filt.matrices[[ds]][common.genes, ] 
  test.genomic.args[nxt.indx] <- list(NULL)
  test.clinical.args[nxt.indx] <- list(NULL)
  test.common.drugs[[nxt.indx]] <- common.drugs.by.data.set[[ds]]
  test.drug.cols[[nxt.indx]] <- drug.name.col
  test.patient.cols[[nxt.indx]] <- expr.patient.id.cols[[ds]]
  test.response.cols[[nxt.indx]] <- "dss.auc.ll4"
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
  test.patient.cols[[nxt.indx]] <- expr.patient.id.cols[[ds]]
  test.response.cols[[nxt.indx]] <- "dss.auc.ll4"
} ## if("sanger" %in% data.sets.to.analyze) 

## END Define training and test data sets

cutoff <- 0.1

genes <- read.table("../common-resources/Geeleher-genes.tsv", sep="\t", header=TRUE, stringsAsFactors = FALSE)$symbol
gene.subset <- genes[grepl(genes, pattern="ENSG")]
genes.to.translate <- genes[!grepl(genes, pattern="ENSG")]
trns <- symbols.to.ensg.mapping(genes.to.translate)
gene.subset <- unique(c(gene.subset, trns$ensg))

for(train.indx in 1:length(train.set.names)) {
   train.expr.args[[train.indx]] <- train.expr.args[[train.indx]][common.genes, ]
}

for(test.indx in 1:length(test.set.names)) {
   test.expr.args[[test.indx]] <- test.expr.args[[test.indx]][common.genes, ]
}

## BEGIN modeling

train.indx <- 1
train.set <- train.set.names[[train.indx]]
l <- prepare.drug.response.and.expr.matrices(train.dss.args[[train.indx]], train.expr.args[[train.indx]], drugs = train.common.drugs[[train.indx]], 
                                             drug.col = train.drug.cols[[train.indx]], patient.col = train.patient.cols[[train.indx]], 
                                             response.col = train.response.cols[[train.indx]])
ohsu.drc <- l[["drc.df"]]

test.indx <- 1
test.set <- test.set.names[[test.indx]]
l <- prepare.drug.response.and.expr.matrices(test.dss.args[[test.indx]], test.expr.args[[test.indx]], drugs = test.common.drugs[[test.indx]], 
                                             drug.col = test.drug.cols[[test.indx]], patient.col = test.patient.cols[[test.indx]], 
                                             response.col = test.response.cols[[test.indx]])
fimm.drc <- l[["drc.df"]]


drug.metadata <- drug.map[, "Mechanism.Targets", drop=F]
drug.metadata$Mechanism.Targets <- as.character(drug.metadata$Mechanism.Targets)
rownames(drug.metadata) <- drug.map$FIMM_DRUG_NAME
tbl <- as.data.frame(table(drug.metadata$Mechanism.Targets))
tbl <- tbl[tbl$Freq == 1,]
drug.metadata$Mechanism.Targets[!is.na(drug.metadata$Mechanism.Targets) & (drug.metadata$Mechanism.Targets %in% as.character(tbl$Var1))] <- NA
drug.metadata$Mechanism.Targets	<- gsub(drug.metadata$Mechanism.Targets, pattern=" inhibitor", replacement="")

orig.ohsu.metadata <- ohsu.metadata
orig.fimm.metadata <- fimm.metadata

frac.na <- 0.25
tmp <- ohsu.drc
flag <- unlist(apply(tmp, 1, function(row) all(is.na(row))))
tmp <- tmp[!flag,]
##flag <- unlist(apply(tmp, 1, function(row) length(which(is.na(row))) > frac.na * length(row)))
##tmp <- tmp[!flag,]
flag <- unlist(apply(tmp, 2, function(col) all(is.na(col))))
tmp <- tmp[, !flag]
##flag <- unlist(apply(tmp, 2, function(col) length(which(is.na(col))) > frac.na * length(col)))
##tmp <- tmp[, !flag]
tmp.scaled <- t(scale(t(tmp)))

subset.and.zscore.matrix <- function(df, row.cutoff = 20, col.cutoff = 20) {
##  frac.na <- 0.25
##  cutoff <- frac.na * nrow(df)
##  cutoff <- 20
  tmp <- df
  flag <- unlist(apply(tmp, 1, function(row) all(is.na(row))))
  tmp <- tmp[!flag,]
  flag <- unlist(apply(tmp, 1, function(row) length(which(!is.na(row))) > row.cutoff))
  tmp <- tmp[flag,]
  flag <- unlist(apply(tmp, 2, function(col) all(is.na(col))))
  tmp <- tmp[, !flag]
  flag <- unlist(apply(tmp, 2, function(col) length(which(!is.na(col))) > col.cutoff))
  tmp <- tmp[, flag]
  tmp.scaled <- t(scale(t(tmp)))
  tmp.scaled
}

## Calculate the correlation of each drug with the GLDS
ohsu.glds <- colMeans(tmp.scaled, na.rm=TRUE)

indices <- 1:nrow(tmp.scaled)
names(indices) <- row.names(tmp.scaled)
drug.glds.corr <- ldply(indices,
                        .fun = function(i) {
                                 drug.resp <- tmp.scaled[i, ]
                                 ct <- cor.test(drug.resp, ohsu.glds)
                                 data.frame(corr = ct$estimate, pval = ct$p.value)
                        })
colnames(drug.glds.corr) <- c("drug", "corr", "pval")
drug.glds.corr <- drug.glds.corr[order(drug.glds.corr$corr, decreasing=FALSE), ]

g <- ggplot()
g <- g + geom_point(data = drug.glds.corr, aes(x = corr, y = -log10(pval)))
df <- subset(drug.glds.corr, pval > 0.1)
## g <- g + geom_text(data = df, aes(x = corr, y = -log10(pval), label = drug))



plot(drug.glds.corr$corr, -log10(drug.glds.corr$pval))



safe.t.test <- function(...) {
  res <- tryCatch({t.test(...)}, error = function(e) { data.frame(p.value = 1) })
  res
}

##cutoff <- 10^-5
##sensitive <- unlist(apply(tmp.scaled, 2, function(col) safe.t.test(col, alternative="greater")$p.value < cutoff))
##insensitive <- unlist(apply(tmp.scaled, 2, function(col) safe.t.test(col, alternative="less")$p.value < cutoff))
sensitive <- unlist(apply(tmp.scaled, 2, function(col) log(safe.t.test(col, alternative="greater")$p.value)))
insensitive <- unlist(apply(tmp.scaled, 2, function(col) log(safe.t.test(col, alternative="less")$p.value)))
df <- data.frame(sensitive = sensitive, insensitive = insensitive)
rownames(df) <- colnames(tmp.scaled)
rownames(df) <- unlist(lapply(rownames(df), function(nm) gsub("X(.*)\\.(.*)", "\\1-\\2", nm)))
ohsu.metadata <- merge(orig.ohsu.metadata, df, by = "row.names", all = TRUE)
rownames(ohsu.metadata) <- ohsu.metadata$Row.names
ohsu.metadata <- ohsu.metadata[, !(colnames(ohsu.metadata) %in% "Row.names")]

hc <- hclust(na.dist(t(tmp.scaled)))
k <- 4
cat(paste0("\nRidge quant cut k = ", k, "\n"))
## print(cutree(hc, k = k))
clusters <- t(t(cutree(hc, k=k)))
colnames(clusters) <- "cluster"
rownames(clusters) <- unlist(lapply(rownames(clusters), function(nm) gsub("X(.*)\\.(.*)", "\\1-\\2", nm)))
ohsu.metadata <- merge(clusters, ohsu.metadata, by = "row.names", all = TRUE)
rownames(ohsu.metadata) <- ohsu.metadata$Row.names
ohsu.metadata <- ohsu.metadata[, !(colnames(ohsu.metadata) %in% "Row.names")]

tmp <- ddply(ohsu.metadata, "cluster", .fun = function(df) mean(df$sensitive))
sensitive.cluster <- tmp$cluster[!is.na(tmp$V1) & (tmp$V1 == min(tmp$V1, na.rm=TRUE))]
tmp <- ddply(ohsu.metadata, "cluster", .fun = function(df) mean(df$insensitive))
insensitive.cluster <- tmp$cluster[!is.na(tmp$V1) & (tmp$V1 == min(tmp$V1, na.rm=TRUE))]
ohsu.metadata$status <- NA
ohsu.metadata$status[!is.na(ohsu.metadata$cluster) & (ohsu.metadata$cluster == sensitive.cluster)] <- "sensitive"
ohsu.metadata$status[!is.na(ohsu.metadata$cluster) & (ohsu.metadata$cluster == insensitive.cluster)] <- "insensitive"

ohsu.sensitive.samples <- data.frame(sample = rownames(ohsu.metadata)[!is.na(ohsu.metadata$status) & (ohsu.metadata$status == "sensitive")])
ohsu.not.sensitive.samples <- data.frame(sample = rownames(ohsu.metadata)[is.na(ohsu.metadata$status) | (ohsu.metadata$status != "sensitive")])
ohsu.intermediate.sensitive.samples <- data.frame(sample = rownames(ohsu.metadata)[is.na(ohsu.metadata$status)])
ohsu.insensitive.samples <- data.frame(sample = rownames(ohsu.metadata)[!is.na(ohsu.metadata$status) & (ohsu.metadata$status == "insensitive")])
write.table(file="ohsu.sensitive.samples.tsv", ohsu.sensitive.samples, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
write.table(file="ohsu.not.sensitive.samples.tsv", ohsu.not.sensitive.samples, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
write.table(file="ohsu.intermediate.sensitive.samples.tsv", ohsu.intermediate.sensitive.samples, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
write.table(file="ohsu.insensitive.samples.tsv", ohsu.insensitive.samples, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)

## Load the original OHSU drug data (i.e., both restricting to same scale as FIMM and to overlapping compounds)
data.set.drug.fit.synIds <- list("ohsu" = "syn11470803", "fimm" = "syn11471586")
data.set.drug.fit.files <- c("ohsu.l4.and.ll4.fits.112017.tsv", "fimm.l4.and.ll4.fits.112017.tsv")

orig.fits <- llply(data.set.drug.fit.synIds, .parallel = FALSE,
              .fun = function(synId) {
                       obj <- synGet(synId, downloadFile = TRUE)
                       file <- getFileLocation(obj)
                       read.table(file, sep="\t", header=TRUE, as.is=TRUE, comment.char="", quote="\"")
              })

## Drop any fits that are flagged as outliers
if(TRUE) {
fct.postfixes <- list("LL.4" = ".ll4", "L.4" = ".l4")
fct <- "LL.4"
for(ds in names(orig.fits)) {
    exclude.cols <- colnames(orig.fits[[ds]])[grepl(pattern=paste0("exclude", fct.postfixes[[fct]]), x=colnames(orig.fits[[ds]]))]
    any.excluded <- unlist(apply(orig.fits[[ds]][, exclude.cols], 1, function(row) any(row)))
    orig.nrow <- nrow(orig.fits[[ds]])
    orig.fits[[ds]] <- orig.fits[[ds]][!any.excluded, ]
    new.nrow <- nrow(orig.fits[[ds]])
    cat(paste0(ds, ": filtered ", orig.nrow - new.nrow, " of the original ", orig.nrow, " fits, leaving ", new.nrow, " fits.\n"))
}
}

##source("../common/drug-mapping-functions/mapping_helpers.R")
##source("../common/drug-mapping-functions/drug_mapping.R")



##converts SMILES string to fingerprint
parseInputFingerprint <- function(input) {
  input.mol <- parse.smiles(input)
  lapply(input.mol, do.typing)
  lapply(input.mol, do.aromaticity)
  lapply(input.mol, do.isotopes)
  fp.inp <- lapply(input.mol, get.fingerprint, type = "extended")
}


drug.structures.tbl <- read.table("../common/drug-mapping-functions/drug_structures.txt", comment.char="", sep="\t", header=TRUE, as.is=TRUE)
drug.structures.tbl <- merge(drug.structures.tbl, drug.map[, c("OHSU_DRUG_NAME", "FIMM_DRUG_NAME")], by.x = "drug", by.y = "FIMM_DRUG_NAME")

cat("Calculating drug fingerprints for common drugs\n")
fp.common <- parseInputFingerprint(as.character(drug.structures.tbl$smiles))
rownames(drug.structures.tbl) <-  drug.structures.tbl$smiles
names(fp.common) <- drug.structures.tbl[names(fp.common), "OHSU_DRUG_NAME"]

## compares every smiles on list one across list 2
cat("Calculating all pairwise distances between drugs\n")
drug.struct.dist <- ldply(fp.common, .fun = function(i) {
  sim <- sapply(fp.common, function(j) {
    1 - distance(i, j, method = "tanimoto")
  })
  sim
})
rownames(drug.struct.dist) <- drug.struct.dist$.id
drug.struct.dist <- drug.struct.dist[, !(colnames(drug.struct.dist) == ".id")]

orig.ohsu.drc <- prepare.drug.response.matrix(orig.fits[["ohsu"]], drug.col = "inhibitor", patient.col = "patient_id", response.col = "auc.ll4")
orig.ohsu.drc.scaled <- subset.and.zscore.matrix(orig.ohsu.drc)
ohsu.drug.response.dist <- as.matrix(na.dist(orig.ohsu.drc.scaled))

tmp <- merge(orig.fits[["fimm"]], drug.map[, c("FIMM_DRUG_NAME", "OHSU_DRUG_NAME")], by.x = "DRUG_NAME", by.y = "FIMM_DRUG_NAME")
orig.fimm.drc <- prepare.drug.response.matrix(tmp, drug.col = "OHSU_DRUG_NAME", patient.col = "PATIENT_ID", response.col = "auc.ll4")
orig.fimm.drc.scaled <- subset.and.zscore.matrix(orig.fimm.drc)
fimm.drug.response.dist <- as.matrix(na.dist(orig.fimm.drc.scaled))


common <- intersect(rownames(drug.struct.dist), rownames(ohsu.drug.response.dist))
tmp1 <- ohsu.drug.response.dist[common, common]
tmp2 <- drug.struct.dist[common, common]

x <- tmp1[upper.tri(tmp1, diag=FALSE)]
y <- tmp2[upper.tri(tmp2, diag=FALSE)]
ct <- cor.test(x, y)

common <- intersect(rownames(fimm.drug.response.dist), rownames(ohsu.drug.response.dist))
tmp1 <- ohsu.drug.response.dist[common, common]
tmp2 <- fimm.drug.response.dist[common, common]

x <- tmp1[upper.tri(tmp1, diag=FALSE)]
y <- tmp2[upper.tri(tmp2, diag=FALSE)]
ct <- cor.test(x, y)



mat <- orig.ohsu.drc.scaled[common,]
hc <- hclust(na.dist(mat))
ct <- cutreeDynamic(hc, distM = as.matrix(na.dist(mat)), minClusterSize = 1)
df <- data.frame(drug = labels(hc), cluster = ct)
ohsu.drug.clusters <- dlply(df, .variables = "cluster", .fun = function(tbl) as.character(tbl$drug))

mat <- orig.fimm.drc.scaled[common,]
hc <- hclust(na.dist(mat))
ct <- cutreeDynamic(hc, distM = as.matrix(na.dist(mat)), minClusterSize = 1)
df <- data.frame(drug = labels(hc), cluster = ct)
fimm.drug.clusters <- dlply(df, .variables = "cluster", .fun = function(tbl) as.character(tbl$drug))

clusters <- list("ohsu" = ohsu.drug.clusters, "fimm" = fimm.drug.clusters)

jaccard.index <- function(set1, set2) {
  if(length(intersect(set1, set2)) == 0) { return(0) }
  length(intersect(set1, set2))/length(union(set1, set2))
}

## Perform meta-clustering using the Jaccard similarity as a distance metric/edge weight
meta.mcl_ <- function(cls, ...) {
  edges <- ldply(cls,
                .fun = function(cluster1) {
                         df.c1 <- ldply(cls,
                                        .fun = function(cluster2) { 
                                                 sim <- jaccard.index(cluster1, cluster2)
                                                 data.frame(jaccard = sim)
                                        })
                         colnames(df.c1) <- c("node2", "jaccard")   
                         df.c1 
                })
  colnames(edges) <- c("node1", "node2", "jaccard")
  ## Remove self-loops
##  flag <- (edges$method1 == edges$method2) & (edges$cluster1 == edges$cluster2)
##  edges <- edges[!flag,]

  adj.matrix <- as.matrix(get.adjacency(graph.data.frame(edges[, c("node1", "node2", "jaccard")]), attr="jaccard"))
  mcl.out <- mcl(adj.matrix, ...)
  df <- data.frame(meta.cluster = 0, name = colnames(adj.matrix))
  if("Cluster" %in% names(mcl.out)) {
    df$meta.cluster <- mcl.out$Cluster
  }
  df
}

meta.mcl <- function(clusters, ...) {
  cls <- unlist(clusters, recursive = FALSE)
  meta.mcl_(cls)
}

df <- meta.mcl(clusters, inflation = 1, addLoops = FALSE)

calc.mcl.weighted.silhouette <- function(clusters, n.iters = 1000, frac.to.downsample = 0.8, ...) {
  cls <- unlist(clusters, recursive = FALSE)
  mcl.all.out <- meta.mcl_(cls, ...)

  all.entries <- llply(clusters, .fun = function(lst) unname(unlist(lst)))
  common <- Reduce(intersect, all.entries)
  for(i in 1:length(all.entries)) {
    if(any(!(union(all.entries[[i]], common) %in% intersect(all.entries[[i]], common)))) {
      cat(paste0("clusters[[", i, "]] != common: ", union(setdiff(all.entries[[i]], common), setdiff(common, all.entries[[i]])), "\n"))
    }
  }

  iters <- 1:n.iters
  mcl.downsampled.clusters <- 
    llply(iters, .parallel = TRUE,
          .fun = function(i) {
                   set.seed(i)
                   downsampled.items <- common[sample.int(n = length(common), size = floor(frac.to.downsample * length(common)))]
                   downsampled.clusters <- llply(clusters,
                                                 .fun = function(method) {
                                                          llply(method,
                                                                .fun = function(cluster) {
                                                                         ds <- cluster[cluster %in% downsampled.items]
                                                                         ds
                                                                 })
                                                 })
                   ds.cls <- unlist(downsampled.clusters, recursive = FALSE)
                   meta.mcl_(ds.cls, ...)
          })

  ## Create a list, where each entry is a list of items in one of the (meta-)clusters
  cluster.members <- unlist(llply(mcl.downsampled.clusters,
                                  .fun = function(mcl.res) {
                                           dlply(mcl.res[mcl.res$meta.cluster != 0,,drop=FALSE], .variables = "meta.cluster", 
                                                 .fun = function(df) paste(df$name, collapse = ","))
                                         }))

  ## Create a consensus matrix, where entry i,j is the fraction of the iterations in which
  ## clusters i and j are clustered together
  all.clusters <- names(cls)
  cons.matrix <- as.data.frame(matrix(data = 0, nrow = length(all.clusters), ncol = length(all.clusters)))
  rownames(cons.matrix) <- all.clusters
  colnames(cons.matrix) <- all.clusters
  for(i in 1:length(cluster.members)) {
    items <- unlist(strsplit(cluster.members[[i]], split=",[ ]*"))
    if(length(items) < 2) { next }
    for(j in 1:(length(items)-1)) {
      itemj <- items[j]
      for(k in (j+1):length(items)) {
        itemk <- items[k]
        cons.matrix[itemj, itemk] <- cons.matrix[itemj, itemk] + 1
        cons.matrix[itemk, itemj] <- cons.matrix[itemk, itemj] + 1
      }
    }
  }
  cons.matrix <- cons.matrix / n.iters

}

stop("stop")






## jaccard dissimilarity
## mcl clustering

hm <- Heatmap(mat, clustering_distance_rows = function(x) na.dist(x), clustering_distance_columns = function(x) na.dist(x), na_col = "black")

## TODO
## - how to cut a cluster (in wgcna)
## - cluster by drug response
## - consensus clustering across ohsu and fimm
## - look at tanimoto distance within cluster
## - look at other annotations
## - predict all within cluster, with covariate for each drug
## - look at overlap with drug covariate model
## - look at prediction accuracy with drug covariate model (limit post hoc to drugs in common)

## hm <- Heatmap(mat, clustering_distance_rows = function(x) na.dist(x), clustering_distance_columns = function(x) na.dist(x), na_col = "black")

pdf(paste0(file.prefix, "-ohsu-all-heatmap.pdf"))
plot.heatmap(tmp.scaled, ohsu.metadata, drug.metadata)
d <- dev.off()

tmp <- fimm.drc
flag <- unlist(apply(tmp, 1, function(row) all(is.na(row))))
tmp <- tmp[!flag,]
flag <- unlist(apply(tmp, 1, function(row) length(which(is.na(row))) > frac.na * length(row)))
tmp <- tmp[!flag,]
flag <- unlist(apply(tmp, 2, function(col) all(is.na(col))))
tmp <- tmp[, !flag]
flag <- unlist(apply(tmp, 2, function(col) length(which(is.na(col))) > frac.na * length(col)))
tmp <- tmp[, !flag]
tmp.scaled <- t(scale(t(tmp)))

##cutoff <- 10^-5
##sensitive <- unlist(apply(tmp.scaled, 2, function(col) safe.t.test(col, alternative="greater")$p.value < cutoff))
##insensitive <- unlist(apply(tmp.scaled, 2, function(col) safe.t.test(col, alternative="less")$p.value < cutoff))
sensitive <- unlist(apply(tmp.scaled, 2, function(col) log(safe.t.test(col, alternative="greater")$p.value)))
insensitive <- unlist(apply(tmp.scaled, 2, function(col) log(safe.t.test(col, alternative="less")$p.value)))
df <- data.frame(sensitive = sensitive, insensitive = insensitive)
rownames(df) <- colnames(tmp.scaled)
fimm.metadata <- merge(orig.fimm.metadata, df, by = "row.names", all = TRUE)
rownames(fimm.metadata) <- fimm.metadata$Row.names
fimm.metadata <- fimm.metadata[, !(colnames(fimm.metadata) %in% "Row.names")]

hc <- hclust(na.dist(t(tmp.scaled)))
k <- 3
cat(paste0("\nRidge quant cut k = ", k, "\n"))
## print(cutree(hc, k = k))
clusters <- t(t(cutree(hc, k=k)))
colnames(clusters) <- "cluster"
rownames(clusters) <- unlist(lapply(rownames(clusters), function(nm) gsub("X(.*)\\.(.*)", "\\1-\\2", nm)))
fimm.metadata <- merge(clusters, fimm.metadata, by = "row.names", all = TRUE)
rownames(fimm.metadata) <- fimm.metadata$Row.names
fimm.metadata <- fimm.metadata[, !(colnames(fimm.metadata) %in% "Row.names")]

tmp <- ddply(fimm.metadata, "cluster", .fun = function(df) mean(df$sensitive))
sensitive.cluster <- tmp$cluster[!is.na(tmp$V1) & (tmp$V1 == min(tmp$V1, na.rm=TRUE))]
tmp <- ddply(fimm.metadata, "cluster", .fun = function(df) mean(df$insensitive))
insensitive.cluster <- tmp$cluster[!is.na(tmp$V1) & (tmp$V1 == min(tmp$V1, na.rm=TRUE))]
fimm.metadata$status <- NA
fimm.metadata$status[!is.na(fimm.metadata$cluster) & (fimm.metadata$cluster == sensitive.cluster)] <- "sensitive"
fimm.metadata$status[!is.na(fimm.metadata$cluster) & (fimm.metadata$cluster == insensitive.cluster)] <- "insensitive"

fimm.sensitive.samples <- data.frame(sample = rownames(fimm.metadata)[!is.na(fimm.metadata$status) & (fimm.metadata$status == "sensitive")])
fimm.not.sensitive.samples <- data.frame(sample = rownames(fimm.metadata)[is.na(fimm.metadata$status) | (fimm.metadata$status != "sensitive")])
fimm.intermediate.sensitive.samples <- data.frame(sample = rownames(fimm.metadata)[is.na(fimm.metadata$status)])
fimm.insensitive.samples <- data.frame(sample = rownames(fimm.metadata)[!is.na(fimm.metadata$status) & (fimm.metadata$status == "insensitive")])
write.table(file="fimm.sensitive.samples.tsv", fimm.sensitive.samples, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
write.table(file="fimm.not.sensitive.samples.tsv", fimm.not.sensitive.samples, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
write.table(file="fimm.intermediate.sensitive.samples.tsv", fimm.intermediate.sensitive.samples, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
write.table(file="fimm.insensitive.samples.tsv", fimm.insensitive.samples, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)

pdf(paste0(file.prefix, "-fimm-all-heatmap.pdf"))
plot.heatmap(tmp.scaled, fimm.metadata, drug.metadata)
d <- dev.off()

save.image(rdata.file)

cat("Exiting\n")
q(status = 0)

