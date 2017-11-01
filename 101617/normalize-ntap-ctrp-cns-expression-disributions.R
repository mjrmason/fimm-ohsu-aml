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

num.cores <- detectCores()
if(!is.na(num.cores) && (num.cores > 1)) {
  suppressPackageStartupMessages(library("doMC"))
  cat(paste("Registering ", num.cores-1, " cores.\n", sep=""))
  registerDoMC(cores=(num.cores-1))
}


fastCor = function(x, Split=1000, Method = "spearman")
{
  if(Method == "spearman"){x        = apply(x,2,rank)}
  groups   = split(1:nrow(x), ceiling(seq_along(1:nrow(x))/min(Split, nrow(x))))
  intraCor = ldply(groups, function(gs){cor(t(x[gs,]), t(x), use = "pairwise")},.parallel=T)[,-1]; 
  rownames(intraCor) = colnames(intraCor)
  return(intraCor)
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

cancer.genes <- read.xlsx("TableS2A.xlsx", sheetIndex = 1, startRow = 3)
cancer.genes <- unique(cancer.genes$Gene)
cancer.genes <- symbols.to.ensg.mapping(cancer.genes)
cancer.genes <- na.omit(cancer.genes$ensg)

## Begin setup (this will need to be changed as the data sets in the intersection change)
data.sets <- c("ctrp", "ntap")
drug.mapping.synId <- "syn11244433"
drug.map.drug.id.cols <- c("ctrp.id", "ntap.smiles")
data.set.expr.synIds <- c("syn11261484", "syn11261464")
names(data.set.expr.synIds) <- data.sets
data.set.expr.files <- c("ctrpv2-log-cpm-gene-expr.tsv", "ntap-log-cpm-gene-expr.tsv")
data.set.id.cols <- c("master_cpd_id", "SMI")
screen.id.cols <- list("ctrp" = c("experiment_id", "master_cpd_id", "master_ccl_id"), "ntap" = c("SMI", "sampleIdentifier"))
patient.id.cols <- list("ctrp" = "ccl_name", "ntap" = "sampleIdentifier")
heatmap.label.cols  <- list("ctrp" = "ccle_primary_hist", "ntap" = NULL)
heatmap.label.cols  <- list("ctrp" = "ccle_primary_site", "ntap" = NULL)

## Pull in the cell line metadata so we can see the cancer/histology type of each sample
synId <- "syn5632192"
obj <- synGet(synId, downloadFile = TRUE)
ctrp.cell.line.metadata <- read.table(getFileLocation(obj), sep="\t", header=TRUE, as.is=TRUE, comment.char="", quote="")

## Subset the cell lines to those of interest
cl.subset <- ctrp.cell.line.metadata[!is.na(ctrp.cell.line.metadata$ccle_primary_site),]
cl.subset <- subset(ctrp.cell.line.metadata, ccle_primary_site == "central_nervous_system")

patient.id.cols <- list("ctrp" = "ccl_name", "ntap" = "sampleIdentifier")
patient.subsets <- list("ctrp" = cl.subset$ccl_name, "ntap" = NULL)

## The synapseIds of folders in which to store the batch-corrected data for each respective data set.
## Here, we will store both in Files/Analysis/Drug Sensitivity/NCATS + CTRP Drug Response Modeling
output.folder.synIds <- list("ctrp" = "syn11244429", "ntap" = "syn11244429")

## A string that will be included in the batch-corrected expression output file name.
## Here "nc" refers to "ntap-ctrp-cns"; cns = central_nervous_system
bc.file.prefix <- "nc-cns"


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
batch <- as.factor(batch.data.sets)
cat("Fitting combat\n")
combat.fit <- fit.ComBat(all.expr, batch = batch, mod = model.matrix(~1, data = batch))

cat("Applying combat\n")
all.expr.combat <- apply.ComBat(all.expr, batch = batch, mod = model.matrix(~1, data = batch), 
                                B.hat = combat.fit[["B.hat"]], gamma.star = combat.fit[["gamma.star"]], delta.star = combat.fit[["delta.star"]])

cat("Calculating PCA of original data\n")
all.pca <- prcomp(t(all.expr))

cat("Calculating PCA of batch-corrected data\n")
all.combat.pca <- prcomp(t(all.expr.combat))

## Store the combat-corrected expression files
for(ds in data.sets) {
  flag <- batch.data.sets == ds
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
par(mfrow=c(2,1))
varplot(all.pca, cols=as.numeric(as.factor(batch.data.sets)), Main = "Uncorrected")
legend("top", legend = unique(as.factor(batch.data.sets)), fill = unique(as.factor(batch.data.sets)), ncol=2)
varplot(all.combat.pca, cols=as.numeric(as.factor(batch.data.sets)), Main = paste0("Study Corrected (", paste(data.sets, collapse=","), ")"))
legend("top", legend = unique(as.factor(batch.data.sets)), fill = unique(as.factor(batch.data.sets)), ncol=2)
dev.off()

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

  for(i.sample in 1:2) {
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
do.correlation.analysis(all.expr = all.expr, data.sets = batch.data.sets, gene.subset = common.genes, data.set.names = data.set.names, postfix = paste0(bc.file.prefix, "-common"))

## Do correlation analysis subsetting to genes with high covariance.
do.correlation.analysis(all.expr = all.expr, data.sets = batch.data.sets, gene.subset = high.cv.genes, data.set.names = data.set.names, postfix = paste0(bc.file.prefix, "-high-cv"))

stop("Successfully completed normalization")
