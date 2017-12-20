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

## Look at overlap across independent fits that have noise due to 
## randomness in glmnet.

## Do this using original response data for FIMM and OHSU and
## drc-reprocessed data.

## Original FIMM data
## DSS.csv (syn8270589)
synId <- "syn8270589"
obj <- synGet(synId, downloadFile = TRUE)
fimm.drc.orig <- read.table(getFileLocation(obj), sep=",", header=TRUE, as.is=TRUE, comment.char="", quote="\"")

colnames(fimm.drc.orig) <- fimm.drc.orig[1,]
## rownames(fimm.drc.orig) <- fimm.drc.orig[,1]

drug.cols <- grepl(pattern="FIMM", colnames(fimm.drc.orig))
sample.rows <- grepl(pattern="FHRB", rownames(fimm.drc.orig))

## samples are rows; drugs are columns
fimm.drc.orig <- fimm.drc.orig[sample.rows,drug.cols]
## Invert so samples are columns
fimm.drc.orig <- t(fimm.drc.orig)

## fimm-log-cpm-gene-expr.tsv (syn11288849)
synId <- "syn11288849"
obj <- synGet(synId, downloadFile = TRUE)
file <- getFileLocation(obj)
fimm.expr <- read.table(file, sep="\t", header=TRUE, as.is=TRUE, comment.char="", quote="\"")

## Original OHSU drug response data 
## interpreted_inhibitor_results_2016_10_10.txt (syn7440451)
synId <- "syn7440451"
obj <- synGet(synId, downloadFile = TRUE)
file <- getFileLocation(obj)
ohsu.drc.orig <- read.table(file, sep="\t", header=TRUE, as.is=TRUE, comment.char="", quote="\"")

synId <- "syn10083817"
obj <- synGet(synId, downloadFile = TRUE)
file <- getFileLocation(obj)
ohsu.rnaseq.sample.summary <- read.table(file, header=TRUE, sep="\t")

ohsu.drc.orig <- merge(ohsu.drc.orig, ohsu.rnaseq.sample.summary, by.x = "lab_id", by.y = "Original_LabID")

collapse.replicates.to.median <- function(drc.df, drug.col, patient.col, response.col) {
  drc.df <- drc.df[, c(drug.col, patient.col, response.col)]
  drc.df <- ddply(drc.df, c(drug.col, patient.col), .fun = function(df) median(df[, response.col]))
  colnames(drc.df) <- c(drug.col, patient.col, response.col)
  drc.df
}

ohsu.drc.orig <- collapse.replicates.to.median(ohsu.drc.orig, drug.col = "drug", patient.col = "SeqID", response.col = "Area_under_the_curve")

ohsu.drc.orig.mat <- prepare.drug.response.matrix(ohsu.drc.orig, drug.col = "drug", patient.col = "SeqID", response.col = "Area_under_the_curve")

## Read in the OHSU expression data
## path <- "ohsu.expr.tsv"
synId <- "syn10083723"
obj <- synGet(synId, downloadFile = TRUE)
file <- getFileLocation(obj)
ohsu.expr <- read.table(file, header=TRUE, sep="\t")
## Convert the column names from X20.00347 to 20-00347
colnames(ohsu.expr) <- gsub("X(.*)\\.(.*)", "\\1-\\2", colnames(ohsu.expr))

library(caret)

set.seed(1234)
indices <- createDataPartition(ohsu.drc.orig[, "drug"], p = 0.5, list = FALSE)
drc1 <- ohsu.drc.orig[indices, ]
drc2 <- ohsu.drc.orig[-indices, ]

## Stratified bootstrapping
library(dplyr)
tbl <- ohsu.drc.orig %>% 
         group_by(drug) %>% 
         sample_frac(size = 1, replace=T)
bar <- ohsu.drc.orig
bar$indx <- 1:nrow(bar)
foo <- merge(tbl, bar)
indices <- foo$indx
bs.drc <- ohsu.drc.orig[indices, ]

x.arg <- ohsu.expr
train.drc <- prepare.drug.response.matrix(drc1, drug.col = "drug", patient.col = "SeqID", response.col = "Area_under_the_curve")
train.drugs <- rownames(train.drc)
fits <- list()
fits[[1]] <- train.model_(x.arg, train.drc, train.drugs, seed = i, use.rf = FALSE, use.svm = FALSE, use.ridge = TRUE, use.mean = FALSE, calc.coefficient.pvals = FALSE)

train.drc <- prepare.drug.response.matrix(drc2, drug.col = "drug", patient.col = "SeqID", response.col = "Area_under_the_curve")
train.drugs <- rownames(train.drc)
fits[[2]] <- train.model_(x.arg, train.drc, train.drugs, seed = i, use.rf = FALSE, use.svm = FALSE, use.ridge = TRUE, use.mean = FALSE, calc.coefficient.pvals = FALSE)

output.model.diffs(ohsu.fits)

num.iters <- 3

x.arg <- fimm.expr
train.drc <- fimm.drc.orig
train.drugs <- rownames(fimm.drc.orig)

fimm.fits <- list()
for(i in 1:num.iters) {
  fimm.fits[[i]] <- train.model_(x.arg, train.drc, train.drugs, seed = i, use.rf = FALSE, use.svm = FALSE, use.ridge = TRUE, use.mean = FALSE, calc.coefficient.pvals = FALSE)
}

x.arg <- ohsu.expr
train.drc <- ohsu.drc.orig.mat
train.drugs <- rownames(train.drc)

ohsu.fits <- list()
for(i in 1:num.iters) {
  print(i)
  ohsu.fits[[i]] <- train.model_(x.arg, train.drc, train.drugs, seed = i, use.rf = FALSE, use.svm = FALSE, use.ridge = TRUE, use.mean = FALSE, calc.coefficient.pvals = FALSE)
}

output.model.diffs <- function(fits) {
  num <- length(fits)
  for(i in 1:(num-1)) {
    for(j in (i+1):num) {
      for(drug.indx in 1:length(fits[[i]])) {
        model.indx <- 2
        if(class(fits[[i]][[drug.indx]][[model.indx]]$fit) != "cv.glmnet") { next }
        if(class(fits[[j]][[drug.indx]][[model.indx]]$fit) != "cv.glmnet") { next }
        if(class(fits[[i]][[drug.indx]][[model.indx]]$fit) != "cv.glmnet") { 
          genesi <- c("Intercept")
        } else {
          coeffsi <- coefficients(fits[[i]][[drug.indx]][[model.indx]]$fit) 
          genesi <- rownames(coeffsi)[coeffsi[,1] > 0]
        }
        if(class(fits[[j]][[drug.indx]][[model.indx]]$fit) != "cv.glmnet") { 
          genesj <- c("Intercept")
        } else {
          coeffsj <- coefficients(fits[[j]][[drug.indx]][[model.indx]]$fit)
          genesj <- rownames(coeffsj)[coeffsi[,1] > 0]
        }
        genesi <- genesi[!grepl(pattern="Intercept", genesi)]
        genesj <- genesj[!grepl(pattern="Intercept", genesj)]
        if(!setequal(genesi, genesj)) {
          cat(paste0("Diff (drug = ", drug.indx, ", ", i , " = ", length(genesi), ", ", j, " = ", length(genesj), "): ", 
              paste(union(setdiff(genesi, genesj), setdiff(genesj, genesi)), collapse=", "), "\n"))
        } else {
          if(length(genesi) > 0) {
            cat(paste0("Same (drug = ", drug.indx, ", ", i, ", ", j , " = ", length(genesi), "): ", paste(genesi, collapse=", "), "\n"))
          }
        }
      }
    }
  }
}

output.model.diffs(ohsu.fits)

for(i in 1:num.iters) {
  model.indx <- 2
  for(drug.indx in 1:length(fits[[i]])) {
print(class(fits[[i]][[drug.indx]][[model.indx]]$fit))
    if(class(fits[[i]][[drug.indx]][[model.indx]]$fit) != "cv.glmnet") { next }
    coeffs <- coefficients(fits[[i]][[drug.indx]][[model.indx]]$fit)
    tmp <- as.vector(coeffs[,1])
    names(tmp) <- names(coeffs)
    coeffs <- tmp
    print(names(coeffs))
  }
}

