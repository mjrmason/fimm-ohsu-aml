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

suppressPackageStartupMessages(library(mygene))
suppressPackageStartupMessages(library(gplots))

## Log in to Synapse
synapseLogin()

## Bring in common code
source("../common/data-preprocessing.R")
source("../common/models.R")
source("../common/plotting.R")
source("../common/utils.R")
##source("../common/drug-mapping-functions/mapping_helpers.R")
##source("../common/drug-mapping-functions/drug_mapping.R")

## Register parallelization
num.cores <- detectCores()
if(!is.na(num.cores) && (num.cores > 1)) {
  suppressPackageStartupMessages(library("doMC"))
  cat(paste("Registering ", num.cores-1, " cores.\n", sep=""))
  registerDoMC(cores=(num.cores-1))
}
num.processes <- num.cores - 1


## Begin setup (this will need to be changed as the data sets in the intersection change)

## Purpose:
## This file models 2 data sets (OHSU and FIMM)
source("fimm-ohsu-setup-040218.R")
source("process-drc-and-expr.R")

file.prefix <- "fimm-ohsu-overlap-120817"
rdata.file <- ".Rdata.overlap.120817"

## End setup

common.genes <- Reduce("intersect", lapply(exprs, rownames))

if(use.subset) {
  common.genes <- intersect(common.genes, gene.subset)
}

## Restrict the matrices to those genes in common across data sets
for(ds in names(fits)) {
  cat(paste0("Subset ", ds, "\n"))
  exprs[[ds]] <- exprs[[ds]][common.genes, ]
  cat(paste0("End subset ", ds, "\n"))
}

## Filter the gene expression data to restrict to highly-expressed genes
expr.filt.matrices <- exprs
if(FALSE) {
  expr.filt.matrices <- list()
  for(ds in names(exprs)) {
    iqrs <- unlist(apply(exprs[[ds]], 1, IQR))
    expr.filt.matrices[[ds]] <- exprs[[ds]][iqrs > median(iqrs[iqrs > 0]),]
  }
}

study.genes <- lapply(exprs, rownames)
study.filt.genes <- lapply(expr.filt.matrices, rownames)
common.genes <- Reduce(intersect, study.filt.genes)

## Restrict the matrices to those genes in common across data sets
for(ds in names(fits)) {
  cat(paste0("Subset ", ds, "\n"))
  expr.filt.matrices[[ds]] <- expr.filt.matrices[[ds]][common.genes, ]
  cat(paste0("End subset ", ds, "\n"))
}

venetoclax.feature.syms <- c("BCL2", "CD14", "FLT3", "MAFB", "LRP1")
venetoclax.feature.gene.tbl <- symbols.to.ensg.mapping(venetoclax.feature.syms)
venetoclax.feature.gene.tbl$ensg <- as.character(venetoclax.feature.gene.tbl$ensg)
ensg <- venetoclax.feature.gene.tbl$ensg

cat(paste0("Length of genes: ", length(common.genes), "\n"))

## Do analysis at pathway- (hallmark, biocarta, kegg, and reactome) as well as at gene-level
gene.sets <- list("gene" = "gene", "kegg" = kegg.gene.sets, "hallmark" = hallmark.gene.sets, "biocarta" = biocarta.gene.sets, "reactome" = reactome.gene.sets)
## Do analysis only at gene-level
## gene.sets <- list("gene" = "gene")

cat("Begin query\n")
baz <- queryMany(common.genes, scopes="ensembl.gene", fields=c("symbol","go"), species="human")
indices <- 1:nrow(baz)
names(indices) <- baz$query
term.types <- c("go.BP", "go.CC", "go.MF")
names(term.types) <- term.types
cat("End query\n")
ensg.to.pathway.list <- llply(term.types, .fun = function(col) {
                              llply(indices, .fun = function(indx) { 
                                unlist(lapply(baz[indx,col], function(df) df$term)) }) })
for(nm in names(gene.sets)[names(gene.sets) != "gene"]) {
  cat(paste0("Invert ", nm, "\n"))
  ensg.to.pathway.list[[nm]] <- invert.list(gene.sets[[nm]])
  cat(paste0("End invert ", nm, "\n"))
}

## BEGIN modeling setup


## END modeling setup

## BEGIN Define training and test data sets

## Based on whether the data sets are listed in data.sets.to.analyze (defined above), the following code will define
source("define-training-and-test-sets.R")

## END Define training and test data sets

## The column named 'num.samples' will be used to define the number of samples (per drug) used for downsampling.
## Define that here, separately for both training and testing.
min.cnt.col <- "num.samples"

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

use.rf <- TRUE
glmnet.alphas <- seq(from = 0, to = 1, by = 0.1)


  keep.forest <- TRUE
  use.cforest <- FALSE

  do.random <- TRUE
  do.subtract.mean.response <- FALSE
  do.model.mean.feature <- TRUE

  cat("Modeling mean response\n")
  mean.response.train.clinical.args <- train.clinical.args
  ## Read in the modeled mean response
  for(idx in 1:length(mean.response.train.clinical.args)) {
    nm1 <- train.set.names[[idx]]
    nm2 <- nm1
    mean.resp.file <- paste0("fimm-ohsu-validate-trained-models-112217-one-ohsu-no-filtering-all-mean-response-", nm1, "-vs-", nm2, "-ridge-gene-prediction.tsv")
    mean.resp.file <- paste0("fimm-ohsu-validate-trained-models-112217-one-ohsu-no-filtering-all-mean-response-", nm1, "-vs-", nm2, "-rf-gene-prediction.tsv")
    mean.resp.file <- paste0("fimm-ohsu-validate-trained-models-112217-one-ohsu-no-filtering-all-mean-response-", nm1, "-vs-", nm2, "-lasso-gene-prediction.tsv")
    mean.resp.file <- paste0("model-glds-", nm1, "-vs-", nm2, "-lasso-gene-prediction.tsv")
    print(mean.resp.file)
    mean.resp <- read.table(file = mean.resp.file, sep="\t", header=TRUE)
    if(grepl(nm2, pattern="ohsu")) {
      ## Convert the column names from X20.00347 to 20-00347
##      colnames(mean.resp) <- gsub("X(.*)\\.(.*)", "\\1-\\2", colnames(mean.resp))
    }
    rownames(mean.resp) <- "mean.resp"
    mean.response.train.clinical.args[[idx]] <- mean.resp
  }

  mean.response.test.clinical.args <- test.clinical.args
  ## Read in the modeled mean response
  for(idx in 1:length(mean.response.test.clinical.args)) {
    nm2 <- test.set.names[[idx]]
    nm1 <- "ohsu"
##    nm1 <- nm2
    mean.resp.file <- paste0("fimm-ohsu-validate-trained-models-112217-one-ohsu-no-filtering-all-mean-response-", nm1, "-vs-", nm2, "-ridge-gene-prediction.tsv")
    mean.resp.file <- paste0("fimm-ohsu-validate-trained-models-112217-one-ohsu-no-filtering-all-mean-response-", nm1, "-vs-", nm2, "-rf-gene-prediction.tsv")
    mean.resp.file <- paste0("fimm-ohsu-validate-trained-models-112217-one-ohsu-no-filtering-all-mean-response-", nm1, "-vs-", nm2, "-lasso-gene-prediction.tsv")
    mean.resp.file <- paste0("model-glds-", nm1, "-vs-", nm2, "-lasso-gene-prediction.tsv")
    print(mean.resp.file)
    mean.resp <- read.table(file = mean.resp.file, sep="\t", header=TRUE)
    if(grepl(nm2, pattern="ohsu")) {
      ## Convert the column names from X20.00347 to 20-00347
##      colnames(mean.resp) <- gsub("X(.*)\\.(.*)", "\\1-\\2", colnames(mean.resp))
    }
    rownames(mean.resp) <- "mean.resp"
    mean.response.test.clinical.args[[idx]] <- mean.resp
  }



  if(do.model.mean.feature) {
  res.aucs[["model.mean.feature"]] <- train.and.test.crossproduct(gene.sets = gene.sets, drug.name.tbl = drug.name.tbl, train.set.names = train.set.names, train.dss.args = train.dss.args, 
                                         train.expr.args = train.expr.args, train.genomic.args = train.genomic.args, train.clinical.args = mean.response.train.clinical.args, 
                                         train.common.drugs = train.common.drugs, train.drugs = train.drugs, train.drug.cols = train.drug.cols, train.patient.cols = train.patient.cols, 
                                         train.response.cols = train.response.cols, test.set.names = test.set.names, test.dss.args = test.dss.args, test.expr.args = test.expr.args, 
                                         test.genomic.args = test.genomic.args, test.clinical.args = mean.response.test.clinical.args, test.common.drugs = test.common.drugs,
                                         test.drug.cols = test.drug.cols, test.patient.cols = test.patient.cols, test.response.cols = test.response.cols, 
                                         seed = 1234, use.rf = TRUE, use.cforest = use.cforest, use.svm = FALSE, use.mean = FALSE, num.processes = num.processes, alphas = c(0, 1), include.mean.response.as.feature = FALSE,
                                         subtract.mean.response = FALSE, keep.forest = keep.forest)
  }

  set.seed(1234)
  random.train.expr.args <- train.expr.args
  for(nm in 1:length(random.train.expr.args)) {
    shuffle.samples <- FALSE
    if(shuffle.samples) {
      cat(paste0("Shuffling samples for ", train.set.names[[nm]], "\n"))
      indices <- sample.int(ncol(random.train.expr.args[[nm]]))
      samples <- colnames(random.train.expr.args[[nm]])
      random.train.expr.args[[nm]] <- random.train.expr.args[[nm]][, indices]
      colnames(random.train.expr.args[[nm]]) <- samples
    } else {
      cat(paste0("Shuffling genes for ", train.set.names[[nm]], "\n"))
      indices <- sample.int(nrow(random.train.expr.args[[nm]]))
      genes <- rownames(random.train.expr.args[[nm]])
      random.train.expr.args[[nm]] <- random.train.expr.args[[nm]][indices, ]
      rownames(random.train.expr.args[[nm]]) <- genes
    }
  }
  random.train.dss.args <- train.dss.args
  if(do.random) {
  res.aucs[["random"]] <- train.and.test.crossproduct(gene.sets = gene.sets, drug.name.tbl = drug.name.tbl, train.set.names = train.set.names, train.dss.args = random.train.dss.args, 
                                         train.expr.args = random.train.expr.args, train.genomic.args = train.genomic.args, train.clinical.args = train.clinical.args, 
                                         train.common.drugs = train.common.drugs, train.drugs = train.drugs, train.drug.cols = train.drug.cols, train.patient.cols = train.patient.cols, 
                                         train.response.cols = train.response.cols, test.set.names = test.set.names, test.dss.args = test.dss.args, test.expr.args = test.expr.args, 
                                         test.genomic.args = test.genomic.args, test.clinical.args = test.clinical.args, test.common.drugs = test.common.drugs,
                                         test.drug.cols = test.drug.cols, test.patient.cols = test.patient.cols, test.response.cols = test.response.cols, 
                                         seed = 1234, use.rf = TRUE, use.cforest = use.cforest, use.svm = FALSE, use.mean = FALSE, num.processes = num.processes, alphas = c(0, 1), include.mean.response.as.feature = FALSE,
                                         subtract.mean.response = FALSE, keep.forest = keep.forest)
  }

  cat("Subtracting mean response\n")

  if(do.subtract.mean.response) {
  res.aucs[["subtract.mean"]] <- train.and.test.crossproduct(gene.sets = gene.sets, drug.name.tbl = drug.name.tbl, train.set.names = train.set.names, train.dss.args = train.dss.args, 
                                         train.expr.args = train.expr.args, train.genomic.args = train.genomic.args, train.clinical.args = train.clinical.args, 
                                         train.common.drugs = train.common.drugs, train.drugs = train.drugs, train.drug.cols = train.drug.cols, train.patient.cols = train.patient.cols, 
                                         train.response.cols = train.response.cols, test.set.names = test.set.names, test.dss.args = test.dss.args, test.expr.args = test.expr.args, 
                                         test.genomic.args = test.genomic.args, test.clinical.args = test.clinical.args, test.common.drugs = test.common.drugs,
                                         test.drug.cols = test.drug.cols, test.patient.cols = test.patient.cols, test.response.cols = test.response.cols, 
                                         seed = 1234, use.rf = TRUE, use.cforest = use.cforest, use.svm = FALSE, use.mean = FALSE, num.processes = num.processes, alphas = c(0, 1), include.mean.response.as.feature = FALSE,
                                         subtract.mean.response = TRUE, keep.forest = keep.forest)
  }

  empty.train.expr.args <- train.expr.args
  ## Read in the modeled mean response
  for(idx in 1:length(empty.train.expr.args)) {
    empty.train.expr.args[[idx]] <- NA
  }

  empty.test.expr.args <- test.expr.args
  ## Read in the modeled mean response
  for(idx in 1:length(empty.test.expr.args)) {
    empty.test.expr.args[[idx]] <- NA
  }

  cat("Including mean response only\n")
  res.aucs[["mean.feature.only"]] <- train.and.test.crossproduct(gene.sets = gene.sets, drug.name.tbl = drug.name.tbl, train.set.names = train.set.names, train.dss.args = train.dss.args, 
                                         train.expr.args = empty.train.expr.args, train.genomic.args = train.genomic.args, train.clinical.args = train.clinical.args, 
                                         train.common.drugs = train.common.drugs, train.drugs = train.drugs, train.drug.cols = train.drug.cols, train.patient.cols = train.patient.cols, 
                                         train.response.cols = train.response.cols, test.set.names = test.set.names, test.dss.args = test.dss.args, test.expr.args = empty.test.expr.args, 
                                         test.genomic.args = test.genomic.args, test.clinical.args = test.clinical.args, test.common.drugs = test.common.drugs,
                                         test.drug.cols = test.drug.cols, test.patient.cols = test.patient.cols, test.response.cols = test.response.cols, 
                                         seed = 1234, use.rf = TRUE, use.svm = FALSE, use.mean = FALSE, num.processes = num.processes, alphas = c(0, 1), include.mean.response.as.feature = TRUE,
                                         subtract.mean.response = FALSE, keep.forest = TRUE)


  cat("Including mean response\n")
  res.aucs[["mean.feature"]] <- train.and.test.crossproduct(gene.sets = gene.sets, drug.name.tbl = drug.name.tbl, train.set.names = train.set.names, train.dss.args = train.dss.args, 
                                         train.expr.args = train.expr.args, train.genomic.args = train.genomic.args, train.clinical.args = train.clinical.args, 
                                         train.common.drugs = train.common.drugs, train.drugs = train.drugs, train.drug.cols = train.drug.cols, train.patient.cols = train.patient.cols, 
                                         train.response.cols = train.response.cols, test.set.names = test.set.names, test.dss.args = test.dss.args, test.expr.args = test.expr.args, 
                                         test.genomic.args = test.genomic.args, test.clinical.args = test.clinical.args, test.common.drugs = test.common.drugs,
                                         test.drug.cols = test.drug.cols, test.patient.cols = test.patient.cols, test.response.cols = test.response.cols, 
                                         seed = 1234, use.rf = TRUE, use.cforest = use.cforest, use.svm = FALSE, use.mean = FALSE, num.processes = num.processes, alphas = c(0, 1), include.mean.response.as.feature = TRUE,
                                         subtract.mean.response = FALSE, keep.forest = keep.forest)

  cat("Unpenalized mean response\n")
  res.aucs[["unpenalized.mean.feature"]] <- train.and.test.crossproduct(gene.sets = gene.sets, drug.name.tbl = drug.name.tbl, train.set.names = train.set.names, train.dss.args = train.dss.args, 
                                         train.expr.args = train.expr.args, train.genomic.args = train.genomic.args, train.clinical.args = train.clinical.args, 
                                         train.common.drugs = train.common.drugs, train.drugs = train.drugs, train.drug.cols = train.drug.cols, train.patient.cols = train.patient.cols, 
                                         train.response.cols = train.response.cols, test.set.names = test.set.names, test.dss.args = test.dss.args, test.expr.args = test.expr.args, 
                                         test.genomic.args = test.genomic.args, test.clinical.args = test.clinical.args, test.common.drugs = test.common.drugs,
                                         test.drug.cols = test.drug.cols, test.patient.cols = test.patient.cols, test.response.cols = test.response.cols, 
                                         seed = 1234, use.rf = FALSE, use.svm = FALSE, use.mean = FALSE, num.processes = num.processes, alphas = c(0, 1), include.mean.response.as.feature = TRUE,
                                         subtract.mean.response = FALSE, keep.forest = TRUE, penalize.mean.response = FALSE)


  cat("Including unpenalized modeled mean feature\n")
  res.aucs[["unpenalized.model.mean.feature"]] <- train.and.test.crossproduct(gene.sets = gene.sets, drug.name.tbl = drug.name.tbl, train.set.names = train.set.names, train.dss.args = train.dss.args, 
                                         train.expr.args = train.expr.args, train.genomic.args = train.genomic.args, train.clinical.args = mean.response.train.clinical.args, 
                                         train.common.drugs = train.common.drugs, train.drugs = train.drugs, train.drug.cols = train.drug.cols, train.patient.cols = train.patient.cols, 
                                         train.response.cols = train.response.cols, test.set.names = test.set.names, test.dss.args = test.dss.args, test.expr.args = test.expr.args, 
                                         test.genomic.args = test.genomic.args, test.clinical.args = mean.response.test.clinical.args, test.common.drugs = test.common.drugs,
                                         test.drug.cols = test.drug.cols, test.patient.cols = test.patient.cols, test.response.cols = test.response.cols, 
                                         seed = 1234, use.rf = TRUE, use.svm = FALSE, use.mean = FALSE, num.processes = num.processes, alphas = c(0, 1), include.mean.response.as.feature = FALSE,
                                         subtract.mean.response = FALSE, keep.forest = TRUE, penalize.mean.response = FALSE)


  cat("No transform")
  res.aucs[["no.transform"]] <- train.and.test.crossproduct(gene.sets = gene.sets, drug.name.tbl = drug.name.tbl, train.set.names = train.set.names, train.dss.args = train.dss.args, 
                                         train.expr.args = train.expr.args, train.genomic.args = train.genomic.args, train.clinical.args = train.clinical.args, 
                                         train.common.drugs = train.common.drugs, train.drugs = train.drugs, train.drug.cols = train.drug.cols, train.patient.cols = train.patient.cols, 
                                         train.response.cols = train.response.cols, test.set.names = test.set.names, test.dss.args = test.dss.args, test.expr.args = test.expr.args, 
                                         test.genomic.args = test.genomic.args, test.clinical.args = test.clinical.args, test.common.drugs = test.common.drugs,
                                         test.drug.cols = test.drug.cols, test.patient.cols = test.patient.cols, test.response.cols = test.response.cols, 
                                         seed = 1234, use.rf = TRUE, use.cforest = use.cforest, use.svm = FALSE, use.mean = FALSE, num.processes = num.processes, alphas = c(0, 1), include.mean.response.as.feature = FALSE,
                                         subtract.mean.response = FALSE, keep.forest = keep.forest)

cat("Done training and testing with AUC\n")

cat(paste0("Saving rdata file: ", rdata.file, "\n"))
save.image(rdata.file)
cat(paste0("Done saving rdata file: ", rdata.file, "\n"))

## Sanitize the above results list into a table that has columns for 
## train.set, test.set, gene.set (i.e., "gene", "kegg", "biocarta", etc.), response (e.g., "AUC", "IC50"),
## metric ("pearson" or "spearman" correlation), val (the correlation), and pval (the associated p-value)
## Note that this table will have one row for each random (down) sample.
tbls <- llply(res.aucs, .fun = function(res.auc) res.list.to.table(res.auc))
tbls <- llply(names(tbls), .fun = function(nm) { tbls[[nm]]$pt.response <- nm; return(tbls[[nm]]) } )
tbl.full <- as.data.frame(do.call("rbind", tbls))

tbls <- llply(res.aucs, .fun = function(res.auc) extract.table.2(res.auc))
tbls <- llply(names(tbls), .fun = function(nm) { tbls[[nm]]$pt.response <- nm; return(tbls[[nm]]) } )
trn.tbl <- as.data.frame(do.call("rbind", tbls))

max.lengths <- list()
max.lengths[["gene"]] <- length(common.genes)
for(gset in names(gene.sets)[names(gene.sets) != "gene"]) {
  max.lengths[[gset]] <- length(gene.sets[[gset]])
}

universes <- list()
universes[["gene"]] <- common.genes
for(gset in names(gene.sets)[names(gene.sets) != "gene"]) {
  universes[[gset]] <- names(gene.sets[[gset]])
}


cat(paste0("Saving rdata file: ", rdata.file, "again\n"))
save.image(rdata.file)
cat(paste0("Done saving rdata file: ", rdata.file, "\n"))

cat("Sourcing clean-foo\n")
source("clean-foo.R")
cat("Done sourcing clean-foo\n")

cat(paste0("Saving rdata file: ", rdata.file, "\n"))
save.image(rdata.file)
cat(paste0("Done saving rdata file: ", rdata.file, "\n"))

cat("Sourcing compare-feature-overlaps\n")
source("compare-feature-overlaps.R")
cat("Done sourcing compare-feature-overlaps\n")

cat(paste0("Saving rdata file: ", rdata.file, "\n"))
save.image(rdata.file)
cat(paste0("Done saving rdata file: ", rdata.file, "\n"))

cat("Sourcing explore-venetoclax.R")
source("explore-venetoclax.R")
cat("Done sourcing explore-venetoclax.R")

cat("Sourcing plot-validation-performance.R")
source("plot-validation-performance.R")
cat("Done sourcing plot-validation-performance.R")

cat("Sourcing plot-venetoclax-model-and-gene-correlations.R")
source("plot-venetoclax-model-and-gene-correlations.R")
cat("Done sourcing plot-venetoclax-model-and-gene-correlations.R")




cat("Done in find-model-overlap.R\n")
