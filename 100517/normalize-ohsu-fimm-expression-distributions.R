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

common.genes <- intersect(rownames(fimm.expr), rownames(ohsu.expr))
cat(paste0("Number of genes common to FIMM and OHSU data sets ", length(common.genes), "\n"))

source("corr-of-corrs.R")

## Batch correct based on data set: FIMM and OHSU
all.expr <- cbind(fimm.expr[common.genes, ],
                  ohsu.expr[common.genes, ])

data.sets <- c(rep("fimm", ncol(fimm.expr)), rep("ohsu", ncol(ohsu.expr)))

all.res <- do.combat(all.expr, as.factor(data.sets), postfix = "fo-all")

batch.data.sets <- c(rep("fimm", ncol(fimm.expr)), rep("ohsu", ncol(ohsu.expr)))
batch <- as.factor(batch.data.sets)
combat.fit <- fit.ComBat(all.expr, batch = batch, mod = model.matrix(~1, data = batch))

all.expr.combat <- apply.ComBat(all.expr, batch = batch, mod = model.matrix(~1, data = batch), 
                                B.hat = combat.fit[["B.hat"]], gamma.star = combat.fit[["gamma.star"]], delta.star = combat.fit[["delta.star"]])

all.pca <- prcomp(t(all.expr))
all.combat.pca <- prcomp(t(all.expr.combat))

## Store the combat-corrected expression files

## Store the OHSU combat-corrected expression data in Files/BEAT AML Data/Expression Data (synId = "syn7433787")
## Store the FIMM combat-corrected expression data in Files/Processed Data (synId = "syn8270577")
parentIds <- c("syn7433787", "syn8270577")
tags <- c("ohsu", "fimm")
for(i in 1:length(parentIds)) {
  flag <- data.sets == tags[i]
  expr <- all.expr.combat[, flag]
  file <- paste0(tags[i], "-expr-fo-combat.tsv")
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
png("fo-mu-vs-sd.png")
smoothScatter(mu, std.dev, xlab = "Mean Gene Expression", ylab = "SD Gene Expression", main = "Combat-normalized FIMM and OHSU")
dev.off()

cv <- std.dev / mu
names(cv) <- names(std.dev)
scv <- sort(cv, decreasing = TRUE)
png("fo-cv.png")
cutoff <- scv[8000]
plot(density(cv), main = "Study Corrected CV", xlab="CV", ylab="Density", xlim=c(0,0.2))
abline(v = cutoff)
high.cv.genes <- names(cv)[cv >= cutoff]
dev.off()

png("pca-fo-corrected.png")
par(mfrow=c(2,1))
varplot(all.pca, cols=as.numeric(as.factor(batch.data.sets)), Main = "Uncorrected")
legend("top", legend = unique(as.factor(batch.data.sets)), fill = unique(as.factor(batch.data.sets)), ncol=2)
varplot(all.combat.pca, cols=as.numeric(as.factor(batch.data.sets)), Main = "Study Corrected (FIMM, OHSU)")
legend("bottom", legend = unique(as.factor(batch.data.sets)), fill = unique(as.factor(batch.data.sets)), ncol=2)
dev.off()

do.correlation.analysis <- function(fimm.expr, ohsu.expr, gene.subset, postfix) {
  all.expr <- cbind(fimm.expr[gene.subset, ],
                    ohsu.expr[gene.subset, ])

  data.sets <- c(rep("fimm", ncol(fimm.expr)), rep("ohsu", ncol(ohsu.expr)))
  names <- list("ohsu", "fimm")

  fimm.flag <- data.sets == "fimm"
  ohsu.flag <- data.sets == "ohsu"

  flags <- list(ohsu.flag, fimm.flag)

  genes <- gene.subset

  cat(paste0("Doing correlation of correlation for ", postfix, "\n"))
  all.data.corr <- do.corrs.of.corrs(all.expr, data.sets, genes, names, flags, postfix = postfix)
  rm(all.data.corr)

  fimm.n.samples <- ncol(fimm.expr)
  ohsu.n.samples <- ncol(ohsu.expr)
  n.samples <- pmin(fimm.n.samples, ohsu.n.samples)
  cat(paste0("Downsampling FIMM (", fimm.n.samples, "), OHSU (", ohsu.n.samples, ") to ", n.samples, " samples\n"))

  for(i.sample in 1:10) {
    set.seed(i.sample)
    cat(paste0("Doing downsample for ", postfix, "-ds-", i.sample, "\n"))
    fimm.indices <- sample.int(n = fimm.n.samples, size = n.samples, replace = FALSE)
    ohsu.indices <- sample.int(n = ohsu.n.samples, size = n.samples, replace = FALSE)

    all.expr <- cbind(fimm.expr[gene.subset, fimm.indices],
                      ohsu.expr[gene.subset, ohsu.indices])

    data.sets <- c(rep("fimm", ncol(fimm.expr[, fimm.indices])), rep("ohsu", ncol(ohsu.expr[, ohsu.indices])))
    names <- list("ohsu", "fimm")

    fimm.flag <- data.sets == "fimm"
    ohsu.flag <- data.sets == "ohsu"

    flags <- list(ohsu.flag, fimm.flag)
  
    ds.data.corr <- do.corrs.of.corrs(all.expr, data.sets, genes, names, flags, postfix = paste0(postfix, "-ds-", i.sample))
    rm(ds.data.corr)
  }
}

## Do correlation analysis subsetting only to common genes.
do.correlation.analysis(fimm.expr = fimm.expr, ohsu.expr = ohsu.expr, gene.subset = common.genes, postfix = "fo-common") 

## Do correlation analysis subsetting to genes with high covariance.
do.correlation.analysis(fimm.expr = fimm.expr, ohsu.expr = ohsu.expr, gene.subset = high.cv.genes, postfix = "fo-high-cv") 

stop("stop")
