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

## Read in the OHSU expression data
## path <- "ohsu.expr.tsv"
synId <- "syn10083723"
obj <- synGet(synId, downloadFile = TRUE)
file <- getFileLocation(obj)
ohsu.expr <- read.table(file, header=TRUE, sep="\t")
## Convert the column names from X20.00347 to 20-00347
colnames(ohsu.expr) <- gsub("X(.*)\\.(.*)", "\\1-\\2", colnames(ohsu.expr))

## Read in the FIMM expression data (RNASeq.CPM.log2.bc.csv)
## Load (batch corrected) FIMM expression data
obj <- synGet(id="syn8270602", downloadFile = TRUE)
file <- getFileLocation(obj)
fimm.expr <- read.table(file, header=TRUE, sep=",")
## Make the columns samples
fimm.expr <- t(fimm.expr)

## Read in the CTRP expression data
## Read in the CTRP expression data (this is counts)
## CCLE_RNAseq_081117.reads.tsv
synId <- "syn10808299"
obj <- synGet(synId, downloadFile = TRUE)
file <- getFileLocation(obj)
ctrp.cnt.expr <- read.table(file, header=TRUE, sep="\t")

## Convert expr to cpm using edger and log2-transform
ctrp.log.cpm.expr <- cpm(ctrp.cnt.expr, log=TRUE)

synId <- "syn5632192"
obj <- synGet(synId, downloadFile = TRUE)
ctrp.cell.line.metadata <- read.table(getFileLocation(obj), sep="\t", header=TRUE, as.is=TRUE, comment.char="", quote="")

## Samples in CTRP do not have "."s and are all in caps
all.ccls <- colnames(ctrp.log.cpm.expr)
heme.ccls <- all.ccls[all.ccls %in% ctrp.cell.line.metadata[ctrp.cell.line.metadata$ccle_primary_hist == "haematopoietic_neoplasm", "ccl_name"]]
aml.ccls <- all.ccls[all.ccls %in% ctrp.cell.line.metadata[ctrp.cell.line.metadata$ccle_hist_subtype_1 == "acute_myeloid_leukaemia", "ccl_name"]]

all.expr.cols <- intersect(colnames(ctrp.log.cpm.expr), all.ccls)
heme.expr.cols <- intersect(colnames(ctrp.log.cpm.expr), heme.ccls)
aml.expr.cols <- intersect(colnames(ctrp.log.cpm.expr), aml.ccls)

common.genes <- intersect(rownames(fimm.expr), rownames(ohsu.expr))
common.genes <- intersect(common.genes, rownames(ctrp.log.cpm.expr)) 
cat(paste0("Number of genes common to all 3 data sets ", length(common.genes), "\n"))

source("corr-of-corrs.R")

## Batch correct based on data set: FIMM and OHSU
all.expr <- cbind(fimm.expr[common.genes, ],
                  ohsu.expr[common.genes, ],
                  ctrp.log.cpm.expr[common.genes, aml.expr.cols])

data.sets <- c(rep("fimm", ncol(fimm.expr)), rep("ohsu", ncol(ohsu.expr)), rep("ctrp-aml", ncol(ctrp.log.cpm.expr[, aml.expr.cols])))

batch.data.sets <- data.sets
batch <- as.factor(batch.data.sets)
combat.fit <- fit.ComBat(all.expr, batch = batch, mod = model.matrix(~1, data = batch))

all.expr.combat <- apply.ComBat(all.expr, batch = batch, mod = model.matrix(~1, data = batch), 
                                B.hat = combat.fit[["B.hat"]], gamma.star = combat.fit[["gamma.star"]], delta.star = combat.fit[["delta.star"]])

all.pca <- prcomp(t(all.expr))
all.combat.pca <- prcomp(t(all.expr.combat))

## Store the combat-corrected expression files

## Store the OHSU combat-corrected expression data in Files/BEAT AML Data/Expression Data (synId = "syn7433787")
## Store the FIMM combat-corrected expression data in Files/FIMM Data/Processed Data (synId = "syn8270577")
## Store the CTRP combat-corrected expression data in Files/CTRPv2 Data (synId = "syn10288724")
parentIds <- c("syn7433787", "syn8270577", "syn10288724")
tags <- c("ohsu", "fimm", "ctrp-aml")
for(i in 1:length(parentIds)) {
  flag <- data.sets == tags[i]
  expr <- all.expr.combat[, flag]
  file <- paste0(tags[i], "-expr-foc-aml-combat.tsv")
  write.table(file = file, expr, sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
  f <- File(file, parentId = parentIds[i], synapseStore = TRUE)
  synStore(f, executed = NULL, forceVersion = FALSE)
}

tmp <- all.expr.combat
offset <- min(tmp)
if(offset < 0) {
  tmp <- tmp - offset
}
mu <- unlist(apply(tmp, 1, mean, na.rm=TRUE))
std.dev <- unlist(apply(tmp, 1, sd, na.rm=TRUE))
png("foc-mu-vs-sd.png")
smoothScatter(mu, std.dev, xlab = "Mean Gene Expression", ylab = "SD Gene Expression", main = "Combat-normalized FIMM, OHSU, and CTRP AML")
dev.off()

cv <- std.dev / mu
names(cv) <- names(std.dev)
scv <- sort(cv, decreasing = TRUE)
png("foc-cv.png")
cutoff <- scv[8000]
plot(density(cv), main = "Combat-normalized FIMM, OHSU, and CTRP AML", xlab="CV", ylab="Density", xlim=c(0,0.2))
abline(v = cutoff)
high.cv.genes <- names(cv)[cv >= cutoff]
dev.off()

png("pca-foc-corrected.png")
par(mfrow=c(2,1))
varplot(all.pca, cols=as.numeric(as.factor(batch.data.sets)), Main = "Uncorrected")
legend("top", legend = unique(as.factor(batch.data.sets)), fill = unique(as.factor(batch.data.sets)), ncol=2)
varplot(all.combat.pca, cols=as.numeric(as.factor(batch.data.sets)), Main = "Study Corrected (FIMM, OHSU, CTRP AML)")
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
data.set.names <- list("ohsu", "fimm", "ctrp-aml")
do.correlation.analysis(all.expr = all.expr, data.sets = data.sets, gene.subset = common.genes, data.set.names = data.set.names, postfix = "foc-aml-common") 

## Do correlation analysis subsetting to genes with high covariance.
do.correlation.analysis(all.expr = all.expr, data.sets = data.sets, gene.subset = high.cv.genes, data.set.names = data.set.names, postfix = "foc-aml-high-cv") 

save.image(".Rdata")

ctrp.aml.expr <- ctrp.log.cpm.expr[, aml.expr.cols]

set.seed(1234)
fimm1.indices <- sample.int(size = floor(ncol(fimm.expr)/2), n = ncol(fimm.expr), replace = FALSE) 
ohsu1.indices <- sample.int(size = floor(ncol(ohsu.expr)/2), n = ncol(ohsu.expr), replace = FALSE) 
ctrp1.indices <- sample.int(size = floor(ncol(ctrp.aml.expr)/2), n = ncol(ctrp.aml.expr), replace = FALSE) 

sub.expr <- cbind(fimm.expr[common.genes, fimm1.indices],
                  fimm.expr[common.genes, -fimm1.indices],
                  ohsu.expr[common.genes, ohsu1.indices],
                  ohsu.expr[common.genes, -ohsu1.indices],
                  ctrp.aml.expr[common.genes, ctrp1.indices],
                  ctrp.aml.expr[common.genes, -ctrp1.indices])

sub.data.sets <- c(rep("fimm1", ncol(fimm.expr[, fimm1.indices])), rep("fimm2", ncol(fimm.expr[, -fimm1.indices])), 
                   rep("ohsu1", ncol(ohsu.expr[, ohsu1.indices])), rep("ohsu2", ncol(ohsu.expr[, -ohsu1.indices])), 
                   rep("ctrp.aml1", ncol(ctrp.aml.expr[, ctrp1.indices])), rep("ctrp.aml2", ncol(ctrp.aml.expr[, -ctrp1.indices])))
sub.data.set.names <- list("fimm1", "fimm2", "ohsu1", "ohsu2", "ctrp.aml1", "ctrp.aml2")

do.correlation.analysis(all.expr = sub.expr, data.sets = sub.data.sets, gene.subset = common.genes, data.set.names = sub.data.set.names, postfix = "foc-aml-sub-common") 

## Do correlation analysis subsetting to genes with high covariance.
do.correlation.analysis(all.expr = sub.expr, data.sets = sub.data.sets, gene.subset = high.cv.genes, data.set.names = sub.data.set.names, postfix = "foc-aml-sub-high-cv") 

save.image(".Rdata")

stop("Successfully completed normalization")
