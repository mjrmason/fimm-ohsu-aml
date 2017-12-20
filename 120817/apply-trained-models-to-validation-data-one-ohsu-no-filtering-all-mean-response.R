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

## Train on each of the training sets and apply those models to each of the test data sets.
## At both the training and test stages downsample to the number of samples (which vary by drug)
## in the num.{train,test}.samples.per.drug table.
## num.iterations determines the number of down-sampling.
## with.replacement = FALSE indicates downsampling, whereas = TRUE indicates bootstrapping.
## NB: can return fits by setting return.fits = TRUE, but this will take a lot of memory (and disk space when doing a save.image)
res.aucs <- list()
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
  cat("Randoming OHSU expr\n")

  set.seed(1234)
  random.train.expr.args <- train.expr.args
  for(nm in 1:length(random.train.expr.args)) {
    indices <- sample.int(ncol(random.train.expr.args[[nm]]))
    samples <- colnames(random.train.expr.args[[nm]])
    random.train.expr.args[[nm]] <- random.train.expr.args[[nm]][, indices]
    colnames(random.train.expr.args[[nm]]) <- samples
  }
  res.aucs[["random"]] <- train.and.test.crossproduct(gene.sets = gene.sets, drug.name.tbl = drug.name.tbl, train.set.names = train.set.names, train.dss.args = train.dss.args, 
                                         train.expr.args = random.train.expr.args, train.genomic.args = train.genomic.args, train.clinical.args = train.clinical.args, 
                                         train.common.drugs = train.common.drugs, train.drugs = train.drugs, train.drug.cols = train.drug.cols, train.patient.cols = train.patient.cols, 
                                         train.response.cols = train.response.cols, test.set.names = test.set.names, test.dss.args = test.dss.args, test.expr.args = test.expr.args, 
                                         test.genomic.args = test.genomic.args, test.clinical.args = test.clinical.args, test.common.drugs = test.common.drugs,
                                         test.drug.cols = test.drug.cols, test.patient.cols = test.patient.cols, test.response.cols = test.response.cols, 
                                         seed = 1234, use.rf = TRUE, use.svm = FALSE, use.mean = FALSE, num.processes = num.processes, alphas = c(0, 1), include.mean.response.as.feature = FALSE,
                                         subtract.mean.response = TRUE, keep.forest = TRUE)

  cat("Subtracting mean response\n")
  res.aucs[["subtract.mean"]] <- train.and.test.crossproduct(gene.sets = gene.sets, drug.name.tbl = drug.name.tbl, train.set.names = train.set.names, train.dss.args = train.dss.args, 
                                         train.expr.args = train.expr.args, train.genomic.args = train.genomic.args, train.clinical.args = train.clinical.args, 
                                         train.common.drugs = train.common.drugs, train.drugs = train.drugs, train.drug.cols = train.drug.cols, train.patient.cols = train.patient.cols, 
                                         train.response.cols = train.response.cols, test.set.names = test.set.names, test.dss.args = test.dss.args, test.expr.args = test.expr.args, 
                                         test.genomic.args = test.genomic.args, test.clinical.args = test.clinical.args, test.common.drugs = test.common.drugs,
                                         test.drug.cols = test.drug.cols, test.patient.cols = test.patient.cols, test.response.cols = test.response.cols, 
                                         seed = 1234, use.rf = TRUE, use.svm = FALSE, use.mean = FALSE, num.processes = num.processes, alphas = c(0, 1), include.mean.response.as.feature = FALSE,
                                         subtract.mean.response = TRUE, keep.forest = TRUE)
  cat("Including mean response\n")
  res.aucs[["mean.feature"]] <- train.and.test.crossproduct(gene.sets = gene.sets, drug.name.tbl = drug.name.tbl, train.set.names = train.set.names, train.dss.args = train.dss.args, 
                                         train.expr.args = train.expr.args, train.genomic.args = train.genomic.args, train.clinical.args = train.clinical.args, 
                                         train.common.drugs = train.common.drugs, train.drugs = train.drugs, train.drug.cols = train.drug.cols, train.patient.cols = train.patient.cols, 
                                         train.response.cols = train.response.cols, test.set.names = test.set.names, test.dss.args = test.dss.args, test.expr.args = test.expr.args, 
                                         test.genomic.args = test.genomic.args, test.clinical.args = test.clinical.args, test.common.drugs = test.common.drugs,
                                         test.drug.cols = test.drug.cols, test.patient.cols = test.patient.cols, test.response.cols = test.response.cols, 
                                         seed = 1234, use.rf = TRUE, use.svm = FALSE, use.mean = FALSE, num.processes = num.processes, alphas = c(0, 1), include.mean.response.as.feature = TRUE,
                                         subtract.mean.response = FALSE, keep.forest = TRUE)
  cat("No transform")
  res.aucs[["no.transform"]] <- train.and.test.crossproduct(gene.sets = gene.sets, drug.name.tbl = drug.name.tbl, train.set.names = train.set.names, train.dss.args = train.dss.args, 
                                         train.expr.args = train.expr.args, train.genomic.args = train.genomic.args, train.clinical.args = train.clinical.args, 
                                         train.common.drugs = train.common.drugs, train.drugs = train.drugs, train.drug.cols = train.drug.cols, train.patient.cols = train.patient.cols, 
                                         train.response.cols = train.response.cols, test.set.names = test.set.names, test.dss.args = test.dss.args, test.expr.args = test.expr.args, 
                                         test.genomic.args = test.genomic.args, test.clinical.args = test.clinical.args, test.common.drugs = test.common.drugs,
                                         test.drug.cols = test.drug.cols, test.patient.cols = test.patient.cols, test.response.cols = test.response.cols, 
                                         seed = 1234, use.rf = TRUE, use.svm = FALSE, use.mean = FALSE, num.processes = num.processes, alphas = c(0, 1), include.mean.response.as.feature = FALSE,
                                         subtract.mean.response = FALSE, keep.forest = TRUE)
}
cat("Done training and testing with AUC\n")

## save.image(rdata.file)

## Sanitize the above results list into a table that has columns for 
## train.set, test.set, gene.set (i.e., "gene", "kegg", "biocarta", etc.), response (e.g., "AUC", "IC50"),
## metric ("pearson" or "spearman" correlation), val (the correlation), and pval (the associated p-value)
## Note that this table will have one row for each random (down) sample.
tbls <- llply(res.aucs, .fun = function(res.auc) res.list.to.table(res.auc))
tbls <- llply(names(tbls), .fun = function(nm) { tbls[[nm]]$pt.response <- nm; return(tbls[[nm]]) } )
tbl.full <- as.data.frame(do.call("rbind", tbls))

cat("Extracted results\n")
## save.image(rdata.file)

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
tbl$train.set <- factor(tbl$train.set, levels = c("ctrp", "ctrp.all", "ctrp.heme", "ctrp.aml", "ohsu", "ntap"))
tbl$pt.response <- factor(tbl$pt.response, levels = c("random", "no.transform", "subtract.mean", "mean.feature"))
tbl$model <- factor(tbl$model, levels = c("lasso", "ridge", "rf"))
for(resp in unique(tbl$response)) {
  for(met in c("pearson", "spearman")) {
##    for(gset in unique(tbl$gene.set)) { }
##    for(mdl in c("ridge", "lasso", "rf")) { }
      sub <- subset(tbl, metric == met & response == resp & train.set == "ohsu" & test.set == "fimm")
      g <- ggplot(data = sub, aes(x = model, y = val))
      g <- g + facet_grid(gene.set ~ pt.response)
      g <- g + geom_violin()
##    g <- g + geom_beeswarm()
      g <- g + geom_boxplot(width = 0.5)
      g <- g + geom_hline(yintercept = 0, linetype = "dashed")
      g <- g + ylab(paste0(capwords(met), " Correlation"))
      g <- g + ggtitle(paste0(resp, " ", capwords(met), " "))
      g <- g + ylim(c(-1,1))
      file <- paste0(file.prefix, "-", resp, "-", met, "-vs-gene-sets.pdf")
      pdf(file)
      print(g)
      d <- dev.off()
##    glist[[length(glist)+1]] <- g
##    { }
  }
##  do.call("grid.arrange", glist)
}

## Look at performance of those lasso results w/ vs w/o mean.resp (but always require some feature)
for(mdl in unique(tbl$model)) {
  for(resp in unique(tbl$response)) {
    for(met in c("pearson", "spearman")) {
      for(gset in unique(tbl$gene.set)) {
        sub <- subset(tbl, metric == met & response == resp & gene.set == gset & train.set == "ohsu" & test.set == "fimm" & model == mdl & pt.response %in% c("mean.feature", "random"))
        sub$includes.mean.feature <- grepl(pattern="mean.resp", sub$features)
        g <- ggplot(data = sub, aes(x = includes.mean.feature, y = val))
        g <- g + facet_grid(gene.set ~ pt.response)
        g <- g + geom_violin()
        g <- g + geom_boxplot(width = 0.5)
        g <- g + geom_beeswarm()
        g <- g + geom_hline(yintercept = 0, linetype = "dashed")
        g <- g + ylab(paste0(capwords(met), " Correlation"))
        g <- g + ggtitle(paste0(mdl, " ", resp, " ", capwords(met), " "))
        g <- g + ylim(c(-1,1))
        file <- paste0(file.prefix, "-", resp, "-", met, "-", mdl, "-", gset, "-incl-mean-vs-gene-sets.pdf")
        pdf(file)
        print(g)
        d <- dev.off()
      }
    }
  }
}

sub <- subset(tbl, metric == "pearson" & response == "AUC" & train.set == "ohsu" & test.set == "fimm" & model == "lasso" & pt.response == "mean.feature" & 
              gene.set == "gene")
sub$includes.mean.feature <- grepl(pattern="mean.resp", sub$features)

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
flag <- unlist(apply(tmp, 1, function(row) length(which(is.na(row))) > frac.na * length(row)))
tmp <- tmp[!flag,]
flag <- unlist(apply(tmp, 2, function(col) all(is.na(col))))
tmp <- tmp[, !flag]
flag <- unlist(apply(tmp, 2, function(col) length(which(is.na(col))) > frac.na * length(col)))
tmp <- tmp[, !flag]
tmp.scaled <- t(scale(t(tmp)))

library(pcaMethods)

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

pdf(paste0(file.prefix, "-ohsu-with-mean-feature-heatmap.pdf"))
plot.heatmap(tmp.scaled[rownames(tmp.scaled) %in% sub$train.drug[sub$includes.mean.feature == TRUE],], ohsu.metadata, drug.metadata)
d <- dev.off()

pdf(paste0(file.prefix, "-ohsu-without-mean-feature-heatmap.pdf"))
plot.heatmap(tmp.scaled[rownames(tmp.scaled) %in% sub$train.drug[sub$includes.mean.feature == FALSE],], ohsu.metadata, drug.metadata)
d <- dev.off()

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

pdf(paste0(file.prefix, "-fimm-all-heatmap.pdf"))
plot.heatmap(tmp.scaled, fimm.metadata, drug.metadata)
d <- dev.off()

rm(res.aucs)
save.image(rdata.file)

cat("Exiting\n")
q(status = 0)

save.image(rdata.file)

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
        flag <- !is.na(tmp$pval) & (tmp$pval > -log10(0.05))
        add.labels <- FALSE
        if(add.labels && any(flag)) {
          labels <- tmp[flag,]
          g <- g + geom_text(data = labels, aes(x = val, y = pval, label = train.drug))
        }
        g <- g + ggtitle(paste0(resp, " ", capwords(met), " ", mdl, " ", gs))
        print(g)
      }
      d <- dev.off()
##    glist[[length(glist)+1]] <- g
    }
  }
##  do.call("grid.arrange", glist)
}



annotate.table <- function(tbl, map, from.col, to.col, do.sort = FALSE) {
  tbl[, to.col] <- unlist(llply(tbl[, from.col], .parallel = FALSE,
                                .fun = function(key.str) {
                                         keys <- unlist(strsplit(key.str, split=",[ ]*"))
                                         values <- unlist(lapply(keys, function(key) map[[key]]))
                                         values <- unique(values)
                                         values <- values[!is.null(values)]
                                         if(do.sort) { values <- sort(values) }
                                         paste(values, collapse=", ")
                                }))
  tbl
}

annotate.table.sig.overlap <- function(tbl, gene.sets, from.col, to.col, universe, do.sort = FALSE) {
  tbl[, to.col] <- unlist(llply(tbl[, from.col], .parallel = FALSE,
                                .fun = function(key.str) {
                                         keys <- unlist(strsplit(key.str, split=",[ ]*"))
                                         genes1 <- intersect(keys, universe)
                                         pvals <- llply(gene.sets, .parallel = FALSE,
                                                      .fun = function(gene.set) {
                                                               genes2 <- intersect(gene.set, universe)
                                                               both <- intersect(genes1, genes2)
                                                               n.both <- length(both)
                                                               n.1 <- length(setdiff(genes1, genes2))
                                                               n.2 <- length(setdiff(genes2, genes1))
                                                               n.neither <- length(setdiff(universe, union(genes1, genes2)))
                                                               mat <- cbind(c(n.both, n.1), c(n.2, n.neither))
                                                               fet <- fisher.test(mat, alternative = "greater")
                                                               pval <- fet$p.value
                                                               pval
                                                      })   
                                         pvals <- unlist(pvals)
                                         qvals <- p.adjust(pvals, method = "BH")
                                         values <- names(gene.sets)[!is.na(qvals) & (qvals < .20)]
                                         if(do.sort) { values <- sort(values) }
                                         paste(values, collapse=", ")
                                }))
  tbl
}


## For a given training set, output all drugs that are positive correlated across all test sets
for(resp in unique(tbl$response)) {
  for(met in c("pearson", "spearman")) {
    for(mdl in c("ridge", "lasso")) {
      for(gs in unique(tbl$gene.set)) {
        for(tr.set in unique(tbl$train.set)) {
          sub <- subset(tbl, metric == met & model == mdl & response == resp & gene.set == gs & train.set == tr.set)
          base.tr.set <- gsub(x = tr.set, pattern = ".train", replacement = "")
          base.tr.set <- gsub(x = base.tr.set, pattern = ".train", replacement = "")
          ## sub <- merge(sub, drug.map, by.x = "train.drug", by.y = drug.map.drug.id.cols[[base.tr.set]], all.x = TRUE)
          test.sets <- unique(sub$test.set)
          names(test.sets) <- test.sets
          test.tbls <- llply(test.sets,
                             .fun = function(ts) {
                                      ts.sub <- subset(sub, test.set == ts)
                                      ts.sub
                             })
          for(ts in test.sets) {
            colnames(test.tbls[[ts]])[grepl(colnames(test.tbls[[ts]]), pattern="n.test")] <- paste0(ts, ".n.test")
            colnames(test.tbls[[ts]])[grepl(colnames(test.tbls[[ts]]), pattern="pval")] <- paste0(ts, ".pval")
            colnames(test.tbls[[ts]])[grepl(colnames(test.tbls[[ts]]), pattern="^val$")] <- paste0(ts, ".val")
          }
          test.tbl <- Reduce(function(...) merge(..., all = FALSE, by = c("train.drug", "train.set", "features", "alpha", "model", "gene.set", "s", "metric", "response", "test.drug", "n.train")), test.tbls)
##          test.tbl <- Reduce(function(...) merge(..., all = FALSE), test.tbls)
          test.tbl <- merge(test.tbl, drug.map[ ,!grepl(pattern="conc", colnames(drug.map), ignore.case=TRUE)], by.x = "train.drug", by.y = drug.name.col, all.x = TRUE)
          pval.cols <- colnames(test.tbl)[grepl(colnames(test.tbl), pattern="pval")]
          val.cols <- colnames(test.tbl)[grepl(colnames(test.tbl), pattern="\\.val")]

          pval.cutoff <- 0.05
          fdr.cutoff <- 0.20
          for(pval.col in pval.cols) {
            fdr.col <- gsub(pval.col, pattern="pval", replacement="fdr")
            val.col <- gsub(pval.col, pattern="pval", replacement="val")
            test.tbl[, fdr.col] <- p.adjust(test.tbl[, pval.col], method = "BH")
            test.set.name <- gsub(pval.col, pattern=".pval", replacement="")
            pos.flag <- !is.na(test.tbl[, val.col]) & (test.tbl[, val.col] > 0)
            pval.flag <- !is.na(test.tbl[, pval.col]) & (test.tbl[, pval.col] < pval.cutoff) 
            fdr.flag <- !is.na(test.tbl[, fdr.col]) & (test.tbl[, fdr.col] < fdr.cutoff) 
            cat(paste0("Num sig for ", paste(resp, met, mdl, gs, tr.set, test.set.name, collapse=" "), ": p < ", pval.cutoff, " (n = ", length(which(pval.flag)), "); p < ", pval.cutoff, " and pos corr (n = ", length(which(pval.flag & pos.flag)), 
                       "); FDR < ", fdr.cutoff, " (n = ", length(which(fdr.flag)), "); FDR < ", fdr.cutoff, " and pos corr (n = ", length(which(fdr.flag & pos.flag)), ")\n")) 
          }
##          flag <- unlist(apply(test.tbl[, pval.cols, drop=FALSE], 1, function(row) any(is.na(row))))
##          test.tbl <- test.tbl[!flag, ,drop=FALSE]

##          sig.tbl <- test.tbl[all.pos.and.sig.flag,,drop=FALSE ]
##          sig.tbl <- test.tbl
          ## Annotate lasso results
          if(mdl == "lasso-foo") {
             ## features column has comma separated list of ensembl ids
             test.tbl <- annotate.table(test.tbl, ensg.to.gene.sets, "features", "pathways", do.sort = TRUE)
             universe <- intersect(common.genes, names(ensg.to.gene.sets))
             test.tbl <- annotate.table.sig.overlap(test.tbl, all.gene.sets, "features", "sig.pathways", universe = universe, do.sort = TRUE)
             test.tbl <- annotate.table(test.tbl, ensg.to.go.gene.sets, "features", "go", do.sort = TRUE)
             universe <- intersect(common.genes, names(ensg.to.go.gene.sets))
             test.tbl <- annotate.table.sig.overlap(test.tbl, go.gene.sets, "features", "sig.go", universe = universe, do.sort = TRUE)
             test.tbl <- annotate.table(test.tbl, ensg.to.symbols, "features", "symbol", do.sort = TRUE)
##             test.tbl <- test.tbl[, c("test.drug", "fimm.val", "fimm.pval", "Mechanism.Targets", "Gene.Targets", "Class.explained", "symbol", "pathways", "sig.pathways", "go", "sig.go", "features")]
          }
cat(paste0("writing ", file, "\n"))
          all.pos.and.sig.flag <- unlist(apply(test.tbl[, pval.cols, drop=FALSE], 1, function(row) all(!is.na(row)) && all(row < 0.05))) &
                                  unlist(apply(test.tbl[, val.cols, drop=FALSE], 1, function(row) all(!is.na(row)) && all(row > 0)))
          none.sig.flag <- unlist(apply(test.tbl[, pval.cols, drop=FALSE], 1, function(row) all(is.na(row)) || all(row[!is.na(row)] > 0.05)))
          file <- paste0(file.prefix, "-", resp, "-", mdl, "-", met, "-vs-gene-sets-all.xls")
          ## Report mean, sd, and 95% CI for n.train, n.test, and correlation of all, of sig (p < 0.05), and of sig (fdr < 20%)
          col <- "n.train"
          flag <- rep(TRUE, nrow(test.tbl))
          flag <- !is.na(test.tbl[, col])
          ci <- quantile(test.tbl[flag, col], probs = c(0.025, 1 - 0.025), na.rm = TRUE)
          print(ci)
          cat(paste0(col, " for ", paste(resp, met, mdl, gs, tr.set, test.set.name, collapse=" "), ": mean = ", mean(test.tbl[flag,col]), " sd = ", sd(test.tbl[flag, col]), " 95% CI = ", as.numeric(ci[1]), " - ", as.numeric(ci[2]), "\n"))
          for(pval.col in pval.cols) {
            fdr.col <- gsub(pval.col, pattern="pval", replacement="fdr")

            n.test.col <- gsub(pval.col, pattern="pval", replacement="n.test")
            val.col <- gsub(pval.col, pattern="pval", replacement="val")
            for(col in c("n.train", n.test.col, val.col)) {

              flag <- !is.na(test.tbl[, pval.col])
              if(any(flag)) {
                ci <- quantile(test.tbl[flag, col], probs = c(0.025, 1 - 0.025), na.rm = TRUE)
                cat(paste0(col, " for ", paste(resp, met, mdl, gs, tr.set, test.set.name, collapse=" "), ": mean = ", mean(test.tbl[flag,col]), " sd = ", sd(test.tbl[flag, col]), " 95% CI = ", as.numeric(ci[1]), " - ", as.numeric(ci[2]), "\n"))
              }

              flag <- !is.na(test.tbl[, pval.col]) & (test.tbl[, pval.col] < pval.cutoff)
              if(any(flag)) {
                ci <- quantile(test.tbl[flag, col], probs = c(0.025, 1 - 0.025), na.rm = TRUE)
                cat(paste0(col, " (p < ", pval.cutoff, ") for ", paste(resp, met, mdl, gs, tr.set, test.set.name, collapse=" "), ": mean = ", mean(test.tbl[flag,col]), " sd = ", sd(test.tbl[flag, col]), " 95% CI = ", as.numeric(ci[1]), " - ", as.numeric(ci[2]), "\n"))
              }

              flag <- !is.na(test.tbl[, fdr.col]) & (test.tbl[, fdr.col] < fdr.cutoff)
              if(any(flag)) {
                ci <- quantile(test.tbl[flag, col], probss = c(0.025, 1 - 0.025), na.rm = TRUE)
                cat(paste0(col, " (fdr < ", fdr.cutoff, ") for ", paste(resp, met, mdl, gs, tr.set, test.set.name, collapse=" "), ": mean = ", mean(test.tbl[flag,col]), " sd = ", sd(test.tbl[flag, col]), " 95% CI = ", as.numeric(ci[1]), " - ", as.numeric(ci[2]), "\n"))
              }
            }

          }
          
          write.table(file = file, test.tbl, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
          file <- paste0(file.prefix, "-", resp, "-", mdl, "-", met, "-vs-gene-sets-all-pos-and-sig.xls")
          write.table(file = file, test.tbl[all.pos.and.sig.flag, ,drop=FALSE], row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
          file <- paste0(file.prefix, "-", resp, "-", mdl, "-", met, "-vs-gene-sets-none-sig.xls")
          write.table(file = file, test.tbl[none.sig.flag, ,drop=FALSE], row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
        }
      }
    }
  }
}

cat("Created volcano plots\n")

save.image(rdata.file)
cat("Successfully exiting\n")
q(status = 0)


get.coef <- function(fit, s) {
  cf <- coefficients(fit, s = s)
  cf <- as.matrix(cf)
  cf[cf[,1] != 0, ]
}


