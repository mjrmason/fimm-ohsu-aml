## rm(list = ls())

suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(synapseClient))

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

## Batch-corrected log2 cpm (or rma) expression files (NB: these exclude 2 outliers in OHSU)
orig.data.set.expr.synIds <- list("ohsu" = "syn11362256", "fimm" = "syn11362257")
orig.data.set.expr.files <- c("ohsu-expr-fimm-ohsu-outlier1-combat.tsv", "fimm-expr-fimm-ohsu-outlier1-combat.tsv")

## Load in the expr files
exprs <- llply(orig.data.set.expr.synIds, .parallel = FALSE,
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

expr.filt.matrices <- list()
for(ds in names(exprs)) {
  iqrs <- unlist(apply(exprs[[ds]], 1, IQR))
  expr.filt.matrices[[ds]] <- exprs[[ds]][iqrs > median(iqrs[iqrs > 0]),]
}

study.genes <- (lapply(exprs, rownames))
study.filt.genes <- (lapply(expr.filt.matrices, rownames))
common.genes <- Reduce(intersect, study.filt.genes)

### common.genes <- Reduce("intersect", lapply(exprs, rownames))

## Restrict the matrices to those genes in common across data sets
for(i in 1:length(expr.filt.matrices)) {
  ds <- names(expr.filt.matrices)[i]
  expr.filt.matrices[[ds]] <- expr.filt.matrices[[ds]][common.genes, ]
}

## Merge the matrices
all.expr <- do.call("cbind", unname(expr.filt.matrices))

get.highly.variable.genes <- function(mat, top = 200, ...) {

  tmp <- mat
  offset <- min(tmp)
  if(offset < 0) {
    tmp <- tmp - offset
  }
  mu <- unlist(apply(tmp, 1, mean, na.rm=TRUE))
  std.dev <- unlist(apply(tmp, 1, sd, na.rm=TRUE))
  cv <- std.dev / mu
  
  ## Limit to expressed genes
  plot(density(mu), ...)
  flag <- mu > 10
  cv <- cv[flag]
  
  cv <- sort(cv, decreasing = TRUE)
  cutoff <- cv[top]
  plot(density(cv))
  abline(v=cutoff)

  highly.variable.genes <- names(cv)[cv > cutoff]
  highly.variable.genes
}

all.variable <- get.highly.variable.genes(all.expr, main = "all")
fimm.variable <- get.highly.variable.genes(expr.filt.matrices[["fimm"]], main = "fimm")
ohsu.variable <- get.highly.variable.genes(expr.filt.matrices[["ohsu"]], main = "ohsu")

fimm.mu <- unlist(apply(expr.filt.matrices[["fimm"]], 1, mean, na.rm = TRUE))
plot(density(fimm.mu), main = "fimm")

ohsu.mu <- unlist(apply(expr.filt.matrices[["ohsu"]], 1, mean, na.rm = TRUE))
plot(density(ohsu.mu), main = "ohsu")

length(intersect(fimm.variable, ohsu.variable))
length(intersect(fimm.variable, all.variable))
length(intersect(ohsu.variable, all.variable))

## Also include the cancer genes
suppressPackageStartupMessages(library("xlsx"))
source("../common/utils.R")
cancer.genes <- read.xlsx("../common-resources/TableS2A.xlsx", sheetIndex = 1, startRow = 3)
cancer.genes <- unique(cancer.genes)
cancer.genes <- symbols.to.ensg.mapping(cancer.genes)
cancer.genes <- na.omit(cancer.genes$ensg)

table(cancer.genes %in% common.genes)

gene.subset <- unique(c(cancer.genes, all.variable, fimm.variable, ohsu.variable))
gene.subset <- gene.subset[gene.subset %in% common.genes]

for(ds in names(exprs)) {
  exprs[[ds]] <- exprs[[ds]][gene.subset, ]
}

source("fimm-ohsu-setup-112217.R")
for(ds in names(exprs)) {
  expr <- exprs[[ds]]
  file <- data.set.expr.files[[ds]]
  write.table(file = file, expr, sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
  f <- File(file, parentId = expr.folder.synIds[[ds]], synapseStore = TRUE)
  synStore(f, executed = NULL, forceVersion = FALSE)

}