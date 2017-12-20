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
suppressPackageStartupMessages(library("pcaMethods"))

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
rdata.file <- ".Rdata.validate.120817.one.ohsu.no.filtering.all.mean.response"

## Over-write some of the variables
## These fit files are assumed to be for AUCs/DSS values calculated across shared concentration ranges across data sets.
## These files may be created by a script such as
## calculate-dss-fimm-ohsu.R
data.set.dss.fit.synIds <- c("syn11362885", "syn11362891")
names(data.set.dss.fit.synIds) <- data.sets
data.set.dss.fit.files <- c("ohsu-fimm-ohsu-outlier1-dss-t0.tsv", "fimm-fimm-ohsu-outlier1-dss-t0.tsv")

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
if(FALSE) {
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

baz <- queryMany(common.genes, scopes="ensembl.gene", fields=c("symbol","go"), species="human")
indices <- 1:nrow(baz)
names(indices) <- baz$query
ensg.to.go.gene.sets <- llply(indices, 
                      .fun = function(indx) {
                               unlist(llply(c("go.BP", "go.CC", "go.MF"), .fun = function(col) unlist(lapply(baz[indx,col], function(df) df$term))))
                             })
go.gene.sets <- invert.list(ensg.to.go.gene.sets)
baz <- baz[!is.na(baz$symbol),]
ensg.to.symbols <- as.list(baz$symbol)
names(ensg.to.symbols) <- baz$query


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

train.indx <- 1
train.set <- train.set.names[[train.indx]]
l <- prepare.drug.response.and.expr.matrices(train.dss.args[[train.indx]], train.expr.args[[train.indx]], drugs = train.common.drugs[[train.indx]], 
                                             drug.col = train.drug.cols[[train.indx]], patient.col = train.patient.cols[[train.indx]], 
                                             response.col = train.response.cols[[train.indx]])
ohsu.drc <- l[["drc.df"]]
ohsu.expr <- l[["expr.df"]]

test.indx <- 1
test.set <- test.set.names[[test.indx]]
l <- prepare.drug.response.and.expr.matrices(test.dss.args[[test.indx]], test.expr.args[[test.indx]], drugs = test.common.drugs[[test.indx]], 
                                             drug.col = test.drug.cols[[test.indx]], patient.col = test.patient.cols[[test.indx]], 
                                             response.col = test.response.cols[[test.indx]])
fimm.drc <- l[["drc.df"]]
fimm.expr <- l[["expr.df"]]


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
flag <- unlist(apply(tmp, 1, function(row) length(which(is.na(row))) > frac.na * length(row)))
tmp <- tmp[!flag,]
flag <- unlist(apply(tmp, 2, function(col) all(is.na(col))))
tmp <- tmp[, !flag]
flag <- unlist(apply(tmp, 2, function(col) length(which(is.na(col))) > frac.na * length(col)))
tmp <- tmp[, !flag]
ohsu.drc.filtered <- tmp
tmp.scaled <- t(scale(t(ohsu.drc.filtered)))
ohsu.drc.filtered.scaled <- tmp.scaled
ohsu.mean.response <- colMeans(ohsu.drc.filtered.scaled, na.rm = TRUE)

library(mice)

impute.pca <- function(mat, mean.response, main.prefix, impute.row = TRUE) {
  mice.imputation.methods <- c("pmm", "midastouch", "sample", "cart", "rf", "mean", "norm", "norm.nob", "norm.boot", "norm.predict", "quadratic", "ri")
  names(mice.imputation.methods) <- mice.imputation.methods
  imps <- llply(mice.imputation.methods, .parallel = TRUE, .fun = function(imp) {
    mat.impute <- NULL
    if(impute.row) {
      mi <- tryCatch({mice(t(mat), method = imp, m = 1)}, error = function(e) { return(NULL) })
      if(!is.null(mi)) { mat.impute <- t(complete(mi)) }
    } else {
      mi <- tryCatch({mice(mat, method = imp, m = 1)}, error = function(e) { return(NULL) })
      if(!is.null(mi)) { mat.impute <- complete(mi) }
    }
    mat.impute
  })
  print(names(imps))
  for(imp in mice.imputation.methods) {
    mat.impute <- imps[[imp]]
    if(!is.null(mat.impute)) {
      pr <- tryCatch({prcomp(t(mat.impute), scale = FALSE, center = FALSE)}, error = function(e) { return(NULL) })
      if(!is.null(pr)) {
        pcs <- pr$x
        var <- summary(pr)$importance[2, ]
        cm <- intersect(names(mean.response), rownames(pcs))
        pcs <- pcs[cm, ]
        resp <- mean.response[cm]
        resp.rank <- rank(abs(resp))
        col <- unlist(lapply(resp, function(r) ifelse(r < 0, "black", "steelblue2")))
        percent.nas <- 100 * length(which(is.na(mat))) / (nrow(mat) * ncol(mat))
        main <- paste0(main.prefix, " ", imp, " imputation: ", round(percent.nas, digits=1), "% NAs")
        str <- "row"
        if(!impute.row) { str <- "col" }
        pdf(paste0(file.prefix, "-", main.prefix, "-imp-", imp, "-pca-", str, ".pdf"))
        par(mfrow=c(3,2), oma=c(0,0,2,0))
        plot(var*100, xlab = "Index", ylab = "Percent Variance Explained")
        for(i in seq(from=1, to=9, by=2)) {
          j <- i + 1
          symbols(pcs[,i], pcs[,j], circles = resp.rank, inches=1/20, bg = col, xlab = paste0("PC", i, " (", 100 * round(as.numeric(var[i]), digits=2), "%)"),
                  ylab = paste0("PC", j, " (", 100 * round(as.numeric(var[j]), digits=2), "%)"))
        }
        title(main, outer=TRUE)
        d <- dev.off()
      }
    }
  }
}

## for pca: samples/observations in rows and variables as columns
do.pca <- function(mat, mean.response, main.prefix) {
  nPcs <- min(20, min(nrow(mat), ncol(mat)))
  pca.methods <- c("nipals", "bpca", "ppca", "svdImpute")
  names(pca.methods) <- pca.methods
  pc.res <- llply(pca.methods, .parallel = TRUE, .fun = function(method) {
    res <- tryCatch({pca(t(mat), method = method, nPcs = nPcs, scale = "none", center = FALSE)}, error = function(e) { return(NULL) })
    res    
  })
  for(method in pca.methods) {
    res <- pc.res[[method]]
    if(!is.null(res)) {
      pcs <- scores(res)
##      var <- sDev(res)^2
      var <- 100 * res@R2
      cm <- intersect(names(mean.response), rownames(pcs))
      pcs <- pcs[cm, ]
      resp <- mean.response[cm]
      resp.rank <- rank(abs(resp))
      col <- unlist(lapply(resp, function(r) ifelse(r < 0, "black", "steelblue2")))
      percent.nas <- 100 * length(which(is.na(mat))) / (nrow(mat) * ncol(mat))
      main <- paste0(main.prefix, " pca method = ", method, ": ", round(percent.nas, digits=1), "% NAs")
      str <- "row"
      pdf(paste0(file.prefix, "-", main.prefix, "-pca-method-", method, ".pdf"))
      par(mfrow=c(3,2), oma=c(0,0,2,0))
      plot(var, xlab = "Index", ylab = "Percent Variance Explained")
      for(i in seq(from=1, to=9, by=2)) {
        j <- i + 1
        symbols(pcs[,i], pcs[,j], circles = resp.rank, inches=1/20, bg = col, xlab = paste0("PC", i, " (", round(as.numeric(var[i]), digits=2), "%)"),
                ylab = paste0("PC", j, " (", round(as.numeric(var[j]), digits=2), "%)"))
      }
      title(main, outer=TRUE)
      d <- dev.off()
    }
  }
}

mat <- ohsu.drc.filtered.scaled
mean.response <- ohsu.mean.response
main.prefix <- "ohsu"
do.pca(mat, mean.response, main.prefix)
impute.pca(mat, mean.response, main.prefix, impute.row = TRUE)
impute.pca(mat, mean.response, main.prefix, impute.row = FALSE)


##cutoff <- 10^-5
##sensitive <- unlist(apply(tmp.scaled, 2, function(col) t.test(col, alternative="greater")$p.value < cutoff))
##insensitive <- unlist(apply(tmp.scaled, 2, function(col) t.test(col, alternative="less")$p.value < cutoff))
sensitive <- unlist(apply(tmp.scaled, 2, function(col) log(t.test(col, alternative="greater")$p.value)))
insensitive <- unlist(apply(tmp.scaled, 2, function(col) log(t.test(col, alternative="less")$p.value)))
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

tmp <- fimm.drc
flag <- unlist(apply(tmp, 1, function(row) all(is.na(row))))
tmp <- tmp[!flag,]
flag <- unlist(apply(tmp, 1, function(row) length(which(is.na(row))) > frac.na * length(row)))
tmp <- tmp[!flag,]
flag <- unlist(apply(tmp, 2, function(col) all(is.na(col))))
tmp <- tmp[, !flag]
flag <- unlist(apply(tmp, 2, function(col) length(which(is.na(col))) > frac.na * length(col)))
tmp <- tmp[, !flag]
fimm.drc.filtered <- tmp
tmp.scaled <- t(scale(t(fimm.drc.filtered)))
fimm.drc.filtered.scaled <- tmp.scaled
fimm.mean.response <- colMeans(fimm.drc.filtered.scaled, na.rm = TRUE)

mat <- fimm.drc.filtered.scaled
mean.response <- fimm.mean.response
main.prefix <- "fimm"
do.pca(mat, mean.response, main.prefix)
impute.pca(mat, mean.response, main.prefix, impute.row = TRUE)
impute.pca(mat, mean.response, main.prefix, impute.row = FALSE)

stop("stop")

##cutoff <- 10^-5
##sensitive <- unlist(apply(tmp.scaled, 2, function(col) t.test(col, alternative="greater")$p.value < cutoff))
##insensitive <- unlist(apply(tmp.scaled, 2, function(col) t.test(col, alternative="less")$p.value < cutoff))
sensitive <- unlist(apply(tmp.scaled, 2, function(col) log(t.test(col, alternative="greater")$p.value)))
insensitive <- unlist(apply(tmp.scaled, 2, function(col) log(t.test(col, alternative="less")$p.value)))
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

cat("Exiting\n")
q(status = 0)

