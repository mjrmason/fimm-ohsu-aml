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
suppressPackageStartupMessages(library(GGally))

suppressPackageStartupMessages(library("GSVA"))
# use GSA to read in a gmt file
suppressPackageStartupMessages(library("GSA"))

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
##data.set.dss.fit.synIds <- c("syn11362885", "syn11362891")
##names(data.set.dss.fit.synIds) <- data.sets
##data.set.dss.fit.files <- c("ohsu-fimm-ohsu-outlier1-dss-t0.tsv", "fimm-fimm-ohsu-outlier1-dss-t0.tsv")

## A string that will be included in any output files.
file.prefix <- "venetoclax"

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

use.subset <- FALSE
if(FALSE) {
aml.logsdon.gene.tbl <- openxlsx:::read.xlsx("../common-resources/41467_2017_2465_MOESM7_ESM.xlsx", sheet = 1, startRow = 3)
aml.logsdon.gene.tbl <- aml.logsdon.gene.tbl[order(as.numeric(aml.logsdon.gene.tbl$SCORE), decreasing=TRUE),]
## Take the top 500, even though the plateau in scores is around 3000
aml.logsdon.genes <- unique(aml.logsdon.gene.tbl$GENE[1:500])
aml.logsdon.genes <- symbols.to.ensg.mapping(aml.logsdon.genes)
aml.logsdon.genes <- na.omit(aml.logsdon.genes$ensg)

## Take all Cancer genes, not just those with Cancer Type == LAML or confidence = A (High), B (Medium), or C (low)
cancer.genes <- openxlsx::read.xlsx("../common-resources/TableS2A.xlsx", sheet = 1, startRow = 3)
cancer.genes <- unique(cancer.genes$Gene)
cancer.genes <- symbols.to.ensg.mapping(cancer.genes)
cancer.genes <- na.omit(cancer.genes$ensg)

## From Suleiman and Ammad
genes <- read.table("../common-resources/gene_names.tsv", sep="\t", header=FALSE, stringsAsFactors = FALSE)$V1
gene.subset <- genes[grepl(genes, pattern="ENSG")]
genes.to.translate <- genes[!grepl(genes, pattern="ENSG")]
trns <- symbols.to.ensg.mapping(genes.to.translate)
gene.subset <- unique(c(gene.subset, trns$ensg))

print(length(aml.logsdon.genes))
print(length(cancer.genes))
print(length(gene.subset))
gene.subset <- unique(union(aml.logsdon.genes, union(cancer.genes, gene.subset)))
print(length(gene.subset))
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

if(FALSE) {
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
}
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

if(FALSE) {
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
}

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

train.indx <- 1
train.set <- train.set.names[[train.indx]]
l <- prepare.drug.response.and.expr.matrices(train.dss.args[[train.indx]], train.expr.args[[train.indx]], drugs = train.common.drugs[[train.indx]], 
                                             drug.col = train.drug.cols[[train.indx]], patient.col = train.patient.cols[[train.indx]], 
                                             response.col = train.response.cols[[train.indx]])
ohsu.drc <- l[["drc.df"]]
ohsu.expr <- l[["expr.df"]]


train.indx <- 2
train.set <- train.set.names[[train.indx]]
l <- prepare.drug.response.and.expr.matrices(train.dss.args[[train.indx]], train.expr.args[[train.indx]], drugs = train.common.drugs[[train.indx]], 
                                             drug.col = train.drug.cols[[train.indx]], patient.col = train.patient.cols[[train.indx]], 
                                             response.col = train.response.cols[[train.indx]])
fimm.drc <- l[["drc.df"]]
fimm.expr <- l[["expr.df"]]

drug.metadata <- drug.map[, "Mechanism.Targets", drop=F]
drug.metadata$Mechanism.Targets <- as.character(drug.metadata$Mechanism.Targets)
rownames(drug.metadata) <- drug.map$FIMM_DRUG_NAME
tbl <- as.data.frame(table(drug.metadata$Mechanism.Targets))
tbl <- tbl[tbl$Freq == 1,]
drug.metadata$Mechanism.Targets[!is.na(drug.metadata$Mechanism.Targets) & (drug.metadata$Mechanism.Targets %in% as.character(tbl$Var1))] <- NA
drug.metadata$Mechanism.Targets	<- gsub(drug.metadata$Mechanism.Targets, pattern=" inhibitor", replacement="")

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

## Do analysis at pathway- (hallmark, biocarta, kegg, and reactome) as well as at gene-level
gene.sets <- list("kegg" = kegg.gene.sets, "hallmark" = hallmark.gene.sets, "biocarta" = biocarta.gene.sets, "reactome" = reactome.gene.sets)


library(GO.db)
library(org.Hs.eg.db)

ontologies <- c("MF", "CC", "BP")
names(ontologies) <- ontologies
xx <- as.list(org.Hs.egGO2EG)
go.lists <- llply(ontologies, 
                  .fun = function(ont) {
                           terms <- xx[Ontology(names(xx)) == ont]
                           names(terms) <- unlist(lapply(names(terms), function(x) Term(x)))
                           llply(terms, 
                                 .fun = function(term) {
                                          genes <- as.character(unlist(mget(as.character(unlist(term)),org.Hs.egENSEMBL)))
                                          na.omit(genes)
                                 })
                  })

for(nm in names(go.lists)) {
  gene.sets[[nm]] <- go.lists[[nm]]
}

library(GSVA)
source("interpret-mean-response-functions.R")

fimm.corrs <- interpret.mean.response(fimm.drc, fimm.expr, prefix = "fimm-mean-response", gene.sets)

ohsu.corrs <- interpret.mean.response(ohsu.drc, ohsu.expr, prefix = "ohsu-mean-response", gene.sets)

save.image(".interpret.Rdata")

nms <- names(fimm.corrs)
names(nms) <- nms
both.corrs <- llply(nms, 
                    .fun = function(nm) {
                             m <- merge(fimm.corrs[[nm]], ohsu.corrs[[nm]], by = "term", suffixes = c(".fimm", ".ohsu"))
                             m
                    })

bm <- ensg.to.sym.mapping(both.corrs[["gene"]]$term) 

gene.tbl <- merge(both.corrs[["gene"]], bm, by.x = "term", by.y = "ensg", all = FALSE)
gene.tbl$cor.mean <- unlist(apply(gene.tbl[, c("cor.fimm", "cor.ohsu")], 1, function(row) mean(as.numeric(row))))
gene.tbl$rank.fimm <- rank(gene.tbl$pval.fimm)
gene.tbl$rank.ohsu <- rank(gene.tbl$pval.ohsu)
gene.tbl$rank.mean <- unlist(apply(gene.tbl[,c("rank.fimm", "rank.ohsu")], 1, function(row) mean(as.numeric(row))))
## Create files for GSEA
write.table(file="fimm-mean-response-gene-correlations.rnk", gene.tbl[,c("symbol", "cor.fimm")], quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(file="fimm-ohsu-mean-response-gene-correlations.rnk", gene.tbl[,c("symbol", "cor.mean")], quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(file="ohsu-mean-response-gene-correlations.rnk", gene.tbl[,c("symbol", "cor.ohsu")], quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

mdr.genes <- c("ABCB1", "ABCG2", "ABCC1", "ABCC2", "MVP")
flag <- gene.tbl$symbol %in% mdr.genes
plot(gene.tbl$cor.fimm[flag], gene.tbl$cor.ohsu[flag])
text(gene.tbl$cor.fimm[flag], gene.tbl$cor.ohsu[flag], gene.tbl$symbol[flag])
gene.tbl[flag,]
gene.tbl[flag,2:10]

baz <- queryMany(both.corrs[["gene"]]$term, scopes="ensembl.gene", fields=c("go"), species="human")
indices <- 1:nrow(baz)
names(indices) <- baz$query
term.types <- c("go.BP", "go.CC", "go.MF")
names(term.types) <- term.types
cat("End query\n")
ensg.to.pathway.list <- llply(term.types, .fun = function(col) {
                              llply(indices, .fun = function(indx) { 
                                unlist(lapply(baz[indx,col], function(df) df$term)) }) })


both.sig.corrs <- llply(both.corrs,
                        .fun = function(df) {
                                 flag <- !is.na(df$pval.fimm) & !is.na(df$pval.ohsu) & (df$pval.fimm < 0.05) & 
                                              (df$pval.ohsu < 0.05) & sign(df$cor.fimm) == sign(df$cor.ohsu)
                                 df[flag, ,drop=F]
                        })

## orig.ohsu.metadata <- ohsu.metadata
## orig.fimm.metadata <- fimm.metadata

frac.na <- 0.25
tmp <- ohsu.drc
flag <- unlist(apply(tmp, 1, function(row) all(is.na(row))))
tmp <- tmp[!flag,]
flag <- unlist(apply(tmp, 2, function(col) all(is.na(col))))
tmp <- tmp[, !flag]
tmp.scaled <- t(scale(t(tmp)))

venetoclax.feature.syms <- c("BCL2", "CD14", "FLT3", "MAFB", "LRP1")
venetoclax.feature.gene.tbl <- symbols.to.ensg.mapping(venetoclax.feature.syms)
rownames(venetoclax.feature.gene.tbl) <-  venetoclax.feature.gene.tbl$ensg


if(FALSE) {
colnames(ohsu.expr) <- unlist(lapply(colnames(ohsu.expr), function(nm) gsub("X(.*)\\.(.*)", "\\1-\\2", nm)))
if(any(venetoclax.feature.gene.tbl$ensg %in% rownames(ohsu.expr))) {
  genes <- venetoclax.feature.gene.tbl$ensg[venetoclax.feature.gene.tbl$ensg %in% rownames(ohsu.expr)]
  ohsu.metadata <- ohsu.metadata[rownames(ohsu.metadata) %in% colnames(ohsu.expr),]
  for(gene in genes) {
    ohsu.metadata[,gene] <- as.numeric(ohsu.expr[gene,rownames(ohsu.metadata)])
  }
}
}

colnames(tmp.scaled) <- gsub("X(.*)\\.(.*)", "\\1-\\2", colnames(tmp.scaled))
colnames(ohsu.expr) <- gsub("X(.*)\\.(.*)", "\\1-\\2", colnames(ohsu.expr))

samples <- colnames(tmp.scaled)
## samples <- intersect(colnames(tmp.scaled), rownames(ohsu.metadata))
bar <- tmp.scaled[, samples]
bar <- rbind(bar, colMeans(bar, na.rm=TRUE))
rownames(bar)[nrow(bar)] <- "mean.response"
if(any(venetoclax.feature.gene.tbl$ensg %in% rownames(ohsu.expr))) {
  genes <- venetoclax.feature.gene.tbl$ensg[venetoclax.feature.gene.tbl$ensg %in% rownames(ohsu.expr)]
  for(gene in genes) {
    bar <- rbind(bar, as.numeric(ohsu.expr[gene,samples]))
    rownames(bar)[nrow(bar)] <- gene
  }
}

bar <- bar[, !is.na(bar["Venetoclax",])]
flag <- grepl(rownames(bar), pattern="ENSG") | rownames(bar) %in% c("Venetoclax", "mean.response")
mat <- bar[flag,order(bar["Venetoclax",])]
flag <- rownames(mat) %in% rownames(venetoclax.feature.gene.tbl)
if(any(flag)) {
  rownames(mat)[flag] <- venetoclax.feature.gene.tbl[rownames(mat)[flag], "symbol"]
}

pdf(paste0("ohsu-venetoclax-heatmap.pdf"))
heatmap(mat, Rowv = NA, Colv = NA, main = "OHSU")
## heatmap(bar[,order(bar["Venetoclax",])], Rowv = NA, Colv = NA)
d <- dev.off()

pdf(paste0("ohsu-venetoclax-ggpairs.pdf"))
g <- ggpairs(as.data.frame(t(mat)), title = "FIMM")
print(g)
d <- dev.off()

pdf(paste0("ohsu-venetoclax-dendro.pdf"))
hc <- hclust(as.dist(1-abs(cor(t(mat), method="pearson"))), method = "ward.D")
plot(hc, main = "OHSU")
d <- dev.off()



## samples <- intersect(colnames(tmp.scaled), rownames(ohsu.metadata))
samples <- colnames(tmp.scaled)
bar <- tmp.scaled[, samples]
bar <- rbind(bar, colMeans(bar, na.rm=TRUE))
rownames(bar)[nrow(bar)] <- "mean.response"
if(any(venetoclax.feature.gene.tbl$ensg %in% rownames(fimm.expr))) {
  genes <- venetoclax.feature.gene.tbl$ensg[venetoclax.feature.gene.tbl$ensg %in% rownames(fimm.expr)]
  for(gene in genes) {
    bar <- rbind(bar, as.numeric(fimm.expr[gene,samples]))
    rownames(bar)[nrow(bar)] <- gene
  }
}

bar <- bar[, !is.na(bar["Venetoclax",])]
flag <- grepl(rownames(bar), pattern="ENSG") | rownames(bar) %in% c("Venetoclax", "mean.response")
mat <- bar[flag,order(bar["Venetoclax",])]
flag <- rownames(mat) %in% rownames(venetoclax.feature.gene.tbl)
if(any(flag)) {
  rownames(mat)[flag] <- venetoclax.feature.gene.tbl[rownames(mat)[flag], "symbol"]
}

pdf(paste0("fimm-venetoclax-heatmap.pdf"))
heatmap(mat, Rowv = NA, Colv = NA, main = "FIMM")
## heatmap(bar[,order(bar["Venetoclax",])], Rowv = NA, Colv = NA)
d <- dev.off()

pdf(paste0("fimm-venetoclax-ggpairs.pdf"))
g <- ggpairs(as.data.frame(t(mat)), title = "FIMM")
print(g)
d <- dev.off()

pdf(paste0("fimm-venetoclax-dendro.pdf"))
hc <- hclust(as.dist(1-abs(cor(t(mat), method="pearson"))), method = "ward.D")
plot(hc, main = "FIMM")
d <- dev.off()

load(".Rdata.overlap.120817")

train.set <- "ohsu"
test.set <- "fimm"
gene.set <- "gene"

## ridge
alpha <- 0
model <- "glmnet"
nm <- "ridge"
foo <- extract.predicted.actual(res.aucs[["mean.feature"]][["all.comparisons"]][[train.set]][[test.set]][[gene.set]],
                                train.drug = "Venetoclax", alpha = alpha, model = model, s = NA)
g1 <- plot.r2(foo$predicted, foo$actual)
g1 <- g1 + ggtitle(paste0("Venetoclax ", train.set, " vs ", test.set, " (", nm, ": ", gene.set, ")"))
g1 <- g1 + ylab("Actual Venetoclax Response")
g1 <- g1 + xlab("Predicted Venetoclax Response")
pdf(paste0("Venetoclax-predicted-vs-actual-", train.set, "-vs-", test.set, "-", nm, ".pdf"))
print(g1)
d <- dev.off()

## lasso
alpha <- 1
model <- "glmnet"
nm <- "lasso"
foo <- extract.predicted.actual(res.aucs[["mean.feature"]][["all.comparisons"]][[train.set]][[test.set]][[gene.set]],
                                train.drug = "Venetoclax", alpha = alpha, model = model, s = NA)
g1 <- plot.r2(foo$predicted, foo$actual)
g1 <- g1 + ggtitle(paste0("Venetoclax ", train.set, " vs ", test.set, " (", nm, ": ", gene.set, ")"))
g1 <- g1 + ylab("Actual Venetoclax Response")
g1 <- g1 + xlab("Predicted Venetoclax Response")
pdf(paste0("Venetoclax-predicted-vs-actual-", train.set, "-vs-", test.set, "-", nm, ".pdf"))
print(g1)
d <- dev.off()

## rf
alpha <- NA
model <- "rf"
nm <- "rf"
foo <- extract.predicted.actual(res.aucs[["mean.feature"]][["all.comparisons"]][[train.set]][[test.set]][[gene.set]],
                                train.drug = "Venetoclax", alpha = alpha, model = model, s = NA)

g1 <- plot.r2(foo$predicted, foo$actual)
g1 <- g1 + ggtitle(paste0("Venetoclax ", train.set, " vs ", test.set, " (", model, ": ", gene.set, ")"))
g1 <- g1 + ylab("Actual Venetoclax Response")
g1 <- g1 + xlab("Predicted Venetoclax Response")
pdf(paste0("Venetoclax-predicted-vs-actual-", train.set, "-vs-", test.set, "-", nm, ".pdf"))
print(g1)
d <- dev.off()