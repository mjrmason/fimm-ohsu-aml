library(synapseClient)
synapseLogin()

suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(Biobase))

library(psych)
library(plyr)
library(dplyr)

suppressPackageStartupMessages(library("foreach"))
suppressPackageStartupMessages(library("parallel"))
library(sva)
library(preprocessCore)
suppressPackageStartupMessages(library("xlsx"))
library(edgeR)
library(ggplot2)
library(gridExtra)

num.cores <- detectCores()
if(!is.na(num.cores) && (num.cores > 1)) {
  suppressPackageStartupMessages(library("doMC"))
  cat(paste("Registering ", num.cores-1, " cores.\n", sep=""))
  registerDoMC(cores=(num.cores-1))
}

symbols.to.ensg.mapping <- function(symbols) {
  # Get a mapping from gene symbol to ensembl id
  ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
  bm <- getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id'), 
              filters = 'hgnc_symbol', 
              values = symbols, 
              mart = ensembl)
  names(bm) <- c("symbol", "ensg")
  bm <- bm[!(bm$ensg %in% c("")),]
  bm
}

cancer.genes <- read.xlsx("../common-resources/TableS2A.xlsx", sheetIndex = 1, startRow = 3)
cancer.genes <- unique(cancer.genes$Gene)
cancer.genes <- symbols.to.ensg.mapping(cancer.genes)
cancer.genes <- na.omit(cancer.genes$ensg)

## Begin setup (this will need to be changed as the data sets in the intersection change)

data.sets <- c("ohsu", "fimm")
data.set.expr.synIds <- c("syn10083723", "syn11288849")
names(data.set.expr.synIds) <- data.sets
data.set.expr.files <- c("ohsu.expr.tsv", "fimm-log-cpm-gene-expr.tsv")

patient.subsets <- list("ohsu" = NULL, "fimm" = NULL)

## The synapseIds of folders in which to store the batch-corrected data for each respective data set.
## FIMM_BEAT AML Collaboration/Files/Analyses/fimm-ohsu (syn11361089)
parentId <- "syn11361089"
output.folder.synIds <- list("ohsu" = "syn11361089", "fimm" = "syn11361089")

## A string that will be included in the batch-corrected expression output file name.
## Here "foc-aml-s-aml" = "fimm-ohsu-ctrp-aml-s-aml"
bc.file.prefix <- "fimm-ohsu-outlier1"

## End setup

cat(paste0("Normalize/batch correct expression of genes shared across data sets: ", paste(data.sets, collapse=","), "\n"))

## Load in the expr files
exprs <- llply(data.set.expr.synIds, .parallel = FALSE,
              .fun = function(synId) {
                       obj <- synGet(synId, downloadFile = TRUE)
                       file <- getFileLocation(obj)
                       read.table(file, sep="\t", header=TRUE, as.is=TRUE, comment.char="", quote="\"")
              })

common.genes <- Reduce("intersect", lapply(exprs, rownames))

## Restrict the matrices to those genes in common across data sets
for(i in 1:length(exprs)) {
  ds <- names(exprs)[i]
  exprs[[ds]] <- exprs[[ds]][common.genes, ]
}

## Fix column names
for(ds in names(exprs)) {
  switch(ds, 
    "ohsu" = {
      ## Convert the column names from X20.00347 to 20-00347
      colnames(exprs[[ds]]) <- gsub("X(.*)\\.(.*)", "\\1-\\2", colnames(exprs[[ds]]))
    },
    "sanger" = {
      colnames(exprs[[ds]]) <- gsub("X(.*)", "\\1", colnames(exprs[[ds]]))
    },
    {
    })
}

## Restrict to samples of interest
for(ds in names(exprs)) {
  if(!is.null(patient.subsets[[ds]])) {
    flag <- !(patient.subsets[[ds]] %in% colnames(exprs[[ds]]))
    if(any(flag)) {
      cat(paste0("The following samples were not in the expression matrix for data set ", ds, ": ", paste(patient.subsets[[ds]][flag], collapse=","), "\n"))
    }
    flag <- colnames(exprs[[ds]]) %in% patient.subsets[[ds]]
    orig.sz <- ncol(exprs[[ds]])
    exprs[[ds]] <- exprs[[ds]][, flag]
    cat(paste0("Limiting ", ds, " data set to ", ncol(exprs[[ds]]), " samples from ", orig.sz, ".\n"))
  }
}

## Merge the matrices.
all.expr <- do.call("cbind", unname(exprs))

source("../common/corr-of-corrs.R")
source("../common/plotting.R")

## Batch correct based on data set
batch.data.sets <- unlist(llply(names(exprs), .parallel = FALSE, .fun = function(nm) rep(nm, ncol(exprs[[nm]]))))
names(batch.data.sets) <- colnames(all.expr)
batch <- as.factor(batch.data.sets)
names(batch) <- colnames(all.expr)

cat("Calculating PCA of original data\n")
all.pca <- prcomp(t(all.expr))

file <- paste0(bc.file.prefix, "-pca-with-outliers.png")
png(file)
## par(mfrow=c(2,1))
## g1 <- varplot(all.pca, cols=as.numeric(as.factor(as.vector(batch.data.sets))), Main = "Uncorrected")
cols <- as.vector(batch.data.sets)
cols <- as.factor(cols)
g1 <- varplot(all.pca, cols=cols, Main = "Uncorrected (w/ Outliers)")
g1 <- g1 + coord_fixed()
print(g1)
dev.off()

x.sort <- sort(abs(all.pca$x[,2]), decreasing=TRUE)
outliers <- names(x.sort)[1:2]
if(!(all(outliers %in% c("14-00800", "20-00062")))) {
  stop(paste0("Was expecting outliers: ", paste(c("14-00800", "20-00062"), collapse=","), " but got: ", paste(outliers, collapse=","), "\n"))
}

## There are two outliers in the PCAs (in both uncorrected and combat-corrected): 14-00800 and 20-00062.
## See if they have exceptional RIN scores or other metrics
## Load the dashboard to get these.

## Read in the latest rna-seq dashboard
## BeatAML_rnaseq_2017_06_09_public_dashboard.xlsx (syn10008793)
## As stated in the README (datawave-2-q3_rnaseq_README.pdf),
## this is a cumulative dashboard that covers all current releases (i.e., DR1 and DR2)
rna.obj <- synGet("syn10008793", downloadFile=TRUE)
file <- getFileLocation(rna.obj)

## RIN scores are in sheet 2
tbl <- openxlsx:::read.xlsx(file, sheet=2)
metrics <- c("RIN_Score")
for(metric in metrics) {
  all.vals <- as.numeric(tbl[,metric])
  all.vals <- all.vals[!is.na(all.vals)]
  num.vals <- length(all.vals)
  cat(paste0(metric, "\n"))
  for(outlier in outliers) {
    val <- as.numeric(subset(tbl, SeqID == outlier)[, metric])
    if(length(val) > 1) { stop("Was only expecting a single value\n") }
    num.less <- length(which(all.vals < val))
    cat(paste0(outlier, " has ", metric, ": ", val, " (", format(100 * num.less / num.vals, digits = 2), " percentile)\n"))
  }
  cat("\n")
}

## Look at some sequencing metrics -- which are on sheet 5
tbl <- openxlsx:::read.xlsx(file, sheet=5)
metrics <- colnames(tbl)
metrics <- metrics[!(metrics %in% c("SampleGroup", "SeqID", "Original_LabID", "Outliers", "FlowCell"))]
for(metric in metrics) {
  all.vals <- as.numeric(tbl[,metric])
  all.vals <- all.vals[!is.na(all.vals)]
  num.vals <- length(all.vals)
  cat(paste0(metric, "\n"))
  for(outlier in outliers) {
    val <- as.numeric(subset(tbl, SeqID == outlier)[, metric])
    if(length(val) > 1) { stop("Was only expecting a single value\n") }
    num.less <- length(which(all.vals < val))
    cat(paste0(outlier, " has ", metric, ": ", val, " (", format(100 * num.less / num.vals, digits = 2), " percentile)\n"))
  }
  cat("\n")
}

## The above shows that these samples have low RIN (< 8th percentile), low aligned (< 7th percentile), and high unmapped (> 94th percentile).

## Drop the outliers and redo pca

cat("Calculating PCA of original data after dropping outliers\n")
all.expr <- all.expr[, !(colnames(all.expr) %in% outliers)]
batch <- batch[colnames(all.expr)]
batch.data.sets <- batch.data.sets[colnames(all.expr)]

all.pca <- prcomp(t(all.expr))

cat("Fitting combat\n")
combat.fit <- fit.ComBat(all.expr, batch = batch, mod = model.matrix(~1, data = batch))

cat("Applying combat\n")
all.expr.combat <- apply.ComBat(all.expr, batch = batch, mod = model.matrix(~1, data = batch), 
                                B.hat = combat.fit[["B.hat"]], gamma.star = combat.fit[["gamma.star"]], delta.star = combat.fit[["delta.star"]])

cat("Calculating PCA of batch-corrected data\n")
all.combat.pca <- prcomp(t(all.expr.combat))

## write.table(all.combat.pca$x[, 1:2], row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)

## Store the combat-corrected expression files
for(ds in data.sets) {
  flag <- as.vector(batch.data.sets) == ds
  expr <- all.expr.combat[, flag]
  file <- paste0(ds, "-expr-", bc.file.prefix, "-combat.tsv")
  write.table(file = file, expr, sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
  f <- File(file, parentId = output.folder.synIds[[ds]], synapseStore = TRUE)
  synStore(f, executed = NULL, forceVersion = FALSE)
}

tmp <- all.expr.combat
offset <- min(tmp)
if(offset < 0) {
  tmp <- tmp - offset
}
mu <- unlist(apply(tmp, 1, mean, na.rm=TRUE))
std.dev <- unlist(apply(tmp, 1, sd, na.rm=TRUE))
file <- paste0(bc.file.prefix, "-mu-vs-sd.png")
png(file)
smoothScatter(mu, std.dev, xlab = "Mean Gene Expression", ylab = "SD Gene Expression", main = paste0("Combat-normalized ", paste(data.sets, collapse=",")))
dev.off()

cv <- std.dev / mu
names(cv) <- names(std.dev)
scv <- sort(cv, decreasing = TRUE)
file <- paste0(bc.file.prefix, "-cv.png")
png(file)
indx <- floor(length(scv)/2)
cutoff <- scv[indx]
plot(density(cv), main = paste0("Combat-normalized ", paste(data.sets, collapse=",")), xlab="CV", ylab="Density")
abline(v = cutoff)
high.cv.genes <- names(cv)[cv >= cutoff]
dev.off()

file <- paste0(bc.file.prefix, "-pca-corrected.png")
png(file)
## par(mfrow=c(2,1))
## g1 <- varplot(all.pca, cols=as.numeric(as.factor(as.vector(batch.data.sets))), Main = "Uncorrected")
cols <- as.vector(batch.data.sets)
cols <- as.factor(cols)
g1 <- varplot(all.pca, cols=cols, Main = "Uncorrected")
g1 <- g1 + coord_fixed()
##legend("top", legend = unique(as.factor(as.vector(batch.data.sets))), fill = unique(as.factor(as.vector(batch.data.sets))), ncol=2)
##legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA, legend = unique(as.factor(as.vector(batch.data.sets))), fill = unique(as.factor(batch.data.sets)), ncol=2)
## g2 <- varplot(all.combat.pca, cols=as.numeric(as.factor(as.vector(batch.data.sets))), Main = paste0("Study Corrected (", paste(data.sets, collapse=","), ")"))
g2 <- varplot(all.combat.pca, cols=cols, Main = paste0("Study Corrected (", paste(data.sets, collapse=","), ")"))
## legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA, legend = unique(as.factor(as.vector(batch.data.sets))), fill = unique(as.factor(as.vector(batch.data.sets)))
##legend("top", legend = unique(as.factor(as.vector(batch.data.sets))), fill = unique(as.factor(as.vector(batch.data.sets))), ncol=2)
g2 <- g2 + coord_fixed()
grid.arrange(g1,g2)
dev.off()

save.image(".Rdata")

do.correlation.analysis <- function(all.expr, data.sets, gene.subset, data.set.names, postfix) {

  flags <- lapply(data.set.names, function(name) data.sets == name)
  names(flags) <- data.set.names
  genes <- gene.subset

  cat(paste0("Doing correlation of correlation for ", postfix, "\n"))
  all.data.corr <- do.corrs.of.corrs(all.expr, data.sets, genes, data.set.names, flags, postfix = postfix)
  rm(all.data.corr)

  n.samples.by.data.set <- (lapply(flags, function(flag) length(which(flag))))
  names(n.samples.by.data.set) <- data.set.names
  n.samples <- min(unlist(n.samples.by.data.set))
  cat(paste0("Downsampling ", 
      paste(unlist(data.set.names, function(ds) paste0(ds, " (", n.samples.by.data.set[[ds]], ")")), collapse=" "), 
      " to ", n.samples, "\n"))

  for(i.sample in 1:10) {
    set.seed(i.sample)
    cat(paste0("Doing downsample for ", postfix, "-ds-", i.sample, "\n"))
    indices.lst <- lapply(data.set.names, function(ds) sample.int(n = n.samples.by.data.set[[ds]], size = n.samples, replace = FALSE))
    names(indices.lst) <- data.set.names
    ds.expr <- lapply(data.set.names, 
                       function(ds) {
                         flag <- flags[[ds]]
                         tmp <- all.expr[, flag]
                         tmp <- tmp[gene.subset, indices.lst[[ds]]]
                         tmp              
                       })
    ds.expr <- do.call("cbind", ds.expr)
    ds.data.sets <- do.call("c", lapply(data.set.names, function(ds) rep(ds, n.samples)))
    ds.flags <- lapply(data.set.names, function(ds) ds.data.sets == ds)
  
    ds.data.corr <- do.corrs.of.corrs(ds.expr, ds.data.sets, genes, data.set.names, ds.flags, postfix = paste0(postfix, "-ds-", i.sample))
    rm(ds.data.corr)
  }
}

## Do correlation analysis subsetting only to common genes.
data.set.names <- as.list(data.sets)
do.correlation.analysis(all.expr = all.expr, data.sets = as.vector(batch.data.sets), gene.subset = common.genes, data.set.names = data.set.names, postfix = paste0(bc.file.prefix, "-common"))

## Do correlation analysis subsetting to genes with high covariance.
do.correlation.analysis(all.expr = all.expr, data.sets = as.vector(batch.data.sets), gene.subset = high.cv.genes, data.set.names = data.set.names, postfix = paste0(bc.file.prefix, "-high-cv"))

cat("Successfully completed normalization")
