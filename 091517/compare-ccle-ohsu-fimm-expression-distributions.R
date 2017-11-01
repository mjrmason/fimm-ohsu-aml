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

## Read in the CTRP expression data (this is counts)
## CCLE_RNAseq_081117.reads.tsv
synId <- "syn10808299"
obj <- synGet(synId, downloadFile = TRUE)
file <- getFileLocation(obj)
ctrp.cnt.expr <- read.table(file, header=TRUE, sep="\t")

synId <- "syn5616092"
obj <- synGet(synId, downloadFile = TRUE)
file <- getFileLocation(obj)
ctrp.expr <- read.table(file, header=TRUE, sep=",")
rownames(ctrp.expr) <- ctrp.expr[,1]
ctrp.expr <- ctrp.expr[,-1]

min.expr <- min(ctrp.expr[ctrp.expr != 0])
ctrp.log2.expr <- ctrp.expr
ctrp.log2.expr[ctrp.log2.expr == 0] <- min.expr
ctrp.log2.expr <- log2(ctrp.log2.expr)

## Convert expr to cpm using edger and log2-transform
ctrp.log.cpm.expr <- cpm(ctrp.cnt.expr, log=TRUE)

if(FALSE) {
iqrs <- unlist(apply(ctrp.log.cpm.expr, 1, IQR))
ctrp.expr.filt <- ctrp.log.cpm.expr[iqrs > median(iqrs[iqrs > 0]),]

iqrs <- unlist(apply(ctrp.log2.expr, 1, IQR))
ctrp.log2.expr.filt <- ctrp.log2.expr[iqrs > median(iqrs[iqrs > 0]),]

iqrs <- unlist(apply(ohsu.expr, 1, IQR))
ohsu.expr.filt <- ohsu.expr[iqrs > median(iqrs[iqrs > 0]),]

iqrs <- unlist(apply(fimm.expr, 1, IQR))
fimm.expr.filt <- fimm.expr[iqrs > median(iqrs[iqrs > 0]),]

common.genes <- intersect(rownames(ctrp.expr.filt), rownames(ohsu.expr.filt))
common.genes <- intersect(common.genes, rownames(fimm.expr.filt)) 
common.genes <- intersect(common.genes, rownames(ctrp.log2.expr.filt)) 
}

common.genes <- intersect(rownames(ctrp.expr), rownames(ohsu.expr))
common.genes <- intersect(common.genes, rownames(fimm.expr)) 
common.genes <- intersect(common.genes, rownames(ctrp.log.cpm.expr)) 
cat(paste0("Number of genes common to all 3 data sets ", length(common.genes), "\n"))



## common.genes <- intersect(common.genes, cancer.genes)
## cat(paste0("Number of genes common to all 3 data sets and implicated in cancer ", length(common.genes), "\n"))

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

##covC = (apply(amtCtrp[,aml.expr.cols],1,sd)/apply(amtCtrp[,aml.expr.cols],1,mean))
##z = covC[common.genes]
##plot(x,z,pch = 19, cex = .25)
##smoothScatter(x,z,pch = 19, cex = .25)

## NB: only using AML in CTRP below
mid.indx <- floor(length(aml.expr.cols)/2)
##ctrp.aml.log.cpm.expr <- ctrp.log.cpm.expr[, aml.expr.cols]
ctrp.aml.log.cpm.expr <- ctrp.log.cpm.expr[, aml.expr.cols[1:mid.indx]]
ctrp.aml2.log.cpm.expr <- ctrp.log.cpm.expr[, aml.expr.cols[(mid.indx+1):length(aml.expr.cols)]]
ctrp.non.aml.log.cpm.expr <- ctrp.log.cpm.expr[, !(colnames(ctrp.log.cpm.expr) %in% aml.expr.cols)]
ctrp.non.aml.heme.log.cpm.expr <- ctrp.log.cpm.expr[, heme.expr.cols]
ctrp.non.aml.heme.log.cpm.expr <- ctrp.non.aml.heme.log.cpm.expr[, !(colnames(ctrp.non.aml.heme.log.cpm.expr) %in% aml.expr.cols)]

## Shift by - offset before doing this
cov.ctrp.aml <- apply(ctrp.aml.log.cpm.expr, 1, sd) / apply(ctrp.aml.log.cpm.expr, 1, mean)
cov.ctrp.aml2 <- apply(ctrp.aml2.log.cpm.expr, 1, sd) / apply(ctrp.aml2.log.cpm.expr, 1, mean)
cov.fimm <- apply(fimm.expr, 1, sd) / apply(fimm.expr, 1, mean)
cov.ohsu <- apply(ohsu.expr, 1, sd) / apply(ohsu.expr, 1, mean)

inter <- intersect(intersect(rownames(fimm.expr), rownames(ohsu.expr)), rownames(ctrp.log.cpm.expr))
plot(apply(ohsu.expr[inter, ], 1, sd), apply(ohsu.expr[inter, ], 1, mean))
plot(apply(fimm.expr[inter, ], 1, sd), apply(fimm.expr[inter, ], 1, mean))

smoothScatter(cov.ohsu[common.genes], cov.fimm[common.genes])

source("corr-of-corrs.R")


## Batch correct based on data set: FIMM, OHSU, or CTRP (ALL)
all.expr <- cbind(fimm.expr[common.genes, ],
                  ohsu.expr[common.genes, ],
                  ctrp.log.cpm.expr[common.genes, ])

data.sets <- c(rep("fimm", ncol(fimm.expr)), rep("ohsu", ncol(ohsu.expr)), rep("ctrp.all", ncol(ctrp.log.cpm.expr)))

all.res <- do.combat(all.expr, as.factor(data.sets), postfix = "foc-all")

sub.data.sets <- data.sets
flag <- ( sub.data.sets == "ctrp.all" ) & ( colnames(all.expr) %in% aml.expr.cols )
sub.data.sets[flag] <- "ctrp.aml"
flag <- ( sub.data.sets == "ctrp.all" ) & ( colnames(all.expr) %in% heme.expr.cols )
sub.data.sets[flag] <- "ctrp.heme"

varplot(all.res[["pca"]], cols=as.numeric(as.factor(sub.data.sets)))
legend("top", legend = unique(as.factor(sub.data.sets)), fill = unique(as.factor(sub.data.sets)), ncol=2)
varplot(all.res[["combat.pca"]], cols=as.numeric(as.factor(sub.data.sets)))
legend("top", legend = unique(as.factor(sub.data.sets)), fill = unique(as.factor(sub.data.sets)), ncol=2)

varplot2(all.res[["pca"]], cols=as.numeric(as.factor(sub.data.sets)))
legend("top", legend = unique(as.factor(sub.data.sets)), fill = unique(as.factor(sub.data.sets)), ncol=2)
varplot2(all.res[["combat.pca"]], cols=as.numeric(as.factor(sub.data.sets)))
legend("top", legend = unique(as.factor(sub.data.sets)), fill = unique(as.factor(sub.data.sets)), ncol=2)

## Batch correct based on data set: FIMM, OHSU, or CTRP (AML) -- only include CTRP AML
aml.expr <- cbind(fimm.expr[common.genes, ],
                  ohsu.expr[common.genes, ],
                  ctrp.log.cpm.expr[common.genes, aml.expr.cols])

aml.data.sets <- c(rep("fimm", ncol(fimm.expr)), rep("ohsu", ncol(ohsu.expr)), rep("ctrp.aml", ncol(ctrp.log.cpm.expr[, aml.expr.cols])))

batch.aml.data.sets <- c(rep("fimm", ncol(fimm.expr)), rep("ohsu", ncol(ohsu.expr)), rep("ctrp", ncol(ctrp.log.cpm.expr[, aml.expr.cols])))
batch <- as.factor(batch.aml.data.sets)
combat.fit <- fit.ComBat(aml.expr, batch = batch, mod = model.matrix(~1, data = batch))

batch.data.sets <- c(rep("fimm", ncol(fimm.expr)), rep("ohsu", ncol(ohsu.expr)), rep("ctrp", ncol(ctrp.log.cpm.expr)))
batch <- as.factor(batch.data.sets)
all.expr.combat <- apply.ComBat(all.expr, batch = batch, mod = model.matrix(~1, data = batch), 
                                B.hat = combat.fit[["B.hat"]], gamma.star = combat.fit[["gamma.star"]], delta.star = combat.fit[["delta.star"]])
all.expr.combat2 <- apply.ComBat2(all.expr, batch = batch, mod = model.matrix(~1, data = batch), 
                                  gamma.star = combat.fit[["gamma.star"]], delta.star = combat.fit[["delta.star"]])
stand.mean <- combat.fit[["stand.mean"]][,1] %*% t(rep(1, ncol(all.expr)))
all.expr.combat3 <- apply.ComBat3(all.expr, batch = batch, mod = model.matrix(~1, data = batch), 
                                  B.hat = combat.fit[["B.hat"]], var.pooled = combat.fit[["var.pooled"]], stand.mean = stand.mean,
                                  gamma.star = combat.fit[["gamma.star"]], delta.star = combat.fit[["delta.star"]])

all.pca <- prcomp(t(all.expr))
all.combat.pca <- prcomp(t(all.expr.combat))
all.combat2.pca <- prcomp(t(all.expr.combat2))
all.combat3.pca <- prcomp(t(all.expr.combat3))

tmp <- all.expr.combat
offset <- min(tmp)
if(offset < 0) {
  tmp <- tmp - offset
}
mu <- unlist(apply(tmp, 1, mean, na.rm=TRUE))
std.dev <- unlist(apply(tmp, 1, sd, na.rm=TRUE))
png("foc-mu-vs-sd.png")
smoothScatter(mu, std.dev, xlab = "Mean Gene Expression", ylab = "SD Gene Expression")
dev.off()

cv <- std.dev / mu
names(cv) <- names(std.dev)
png("foc-cv.png")
plot(density(cv), main = "Study Corrected CV", xlab="CV", ylab="Density", xlim=c(0,0.2))
abline(v = 0.02)
high.cv.genes <- names(cv)[cv > 0.02]
dev.off()

png("pca-foc-ctrp-aml-corrected.png")
par(mfrow=c(2,1))
varplot(all.pca, cols=as.numeric(as.factor(sub.data.sets)), Main = "Uncorrected")
legend("top", legend = unique(as.factor(sub.data.sets)), fill = unique(as.factor(sub.data.sets)), ncol=2)
##varplot(all.combat.pca, cols=as.numeric(as.factor(sub.data.sets)))
##legend("top", legend = unique(as.factor(sub.data.sets)), fill = unique(as.factor(sub.data.sets)), ncol=2)
##varplot(all.combat2.pca, cols=as.numeric(as.factor(sub.data.sets)))
##legend("top", legend = unique(as.factor(sub.data.sets)), fill = unique(as.factor(sub.data.sets)), ncol=2)
varplot(all.combat3.pca, cols=as.numeric(as.factor(sub.data.sets)), Main = "Study Corrected (FIMM, OHSU, CTRP AML)")
legend("top", legend = unique(as.factor(sub.data.sets)), fill = unique(as.factor(sub.data.sets)), ncol=2)
dev.off()

png("pca-foc-ctrp-all-corrected.png")
par(mfrow=c(2,1))
varplot(all.pca, cols=as.numeric(as.factor(sub.data.sets)), Main = "Uncorrected")
legend("top", legend = unique(as.factor(sub.data.sets)), fill = unique(as.factor(sub.data.sets)), ncol=2)
varplot(all.combat.pca, cols=as.numeric(as.factor(sub.data.sets)), Main = "Study Corrected (FIMM, OHSU, CTRP ALL)")
legend("top", legend = unique(as.factor(sub.data.sets)), fill = unique(as.factor(sub.data.sets)), ncol=2)
dev.off()

do.correlation.analysis <- function(fimm.expr, ohsu.expr, ctrp.expr, ctrp.aml.expr, ctrp.aml2.expr, ctrp.non.aml.expr, ctrp.non.aml.heme.expr, gene.subset, postfix) {
  all.expr <- cbind(fimm.expr[gene.subset, ],
                    ohsu.expr[gene.subset, ],
                    ctrp.expr[gene.subset, ],
                    ctrp.aml.expr[gene.subset, ],
                    ctrp.aml2.expr[gene.subset, ],
                    ctrp.non.aml.expr[gene.subset, ],
                    ctrp.non.aml.heme.expr[gene.subset, ])

  data.sets <- c(rep("fimm", ncol(fimm.expr)), rep("ohsu", ncol(ohsu.expr)), rep("ctrp.all", ncol(ctrp.expr)), rep("ctrp.aml", ncol(ctrp.aml.expr)), rep("ctrp.aml2", ncol(ctrp.aml2.expr)), 
                 rep("ctrp.non.aml", ncol(ctrp.non.aml.expr)), rep("ctrp.non.aml.heme", ncol(ctrp.non.aml.heme.expr)))
  names <- list("ohsu", "fimm", "ctrp.all", "ctrp.aml", "ctrp.aml2", "ctrp.non.aml", "ctrp.non.aml.heme")

  fimm.flag <- data.sets == "fimm"
  ohsu.flag <- data.sets == "ohsu"
  ctrp.flag <- data.sets == "ctrp.all"
  ctrp.aml.flag <- data.sets == "ctrp.aml"
  ctrp.aml2.flag <- data.sets == "ctrp.aml2"
  ctrp.non.aml.flag <- data.sets == "ctrp.non.aml"
  ctrp.non.aml.heme.flag <- data.sets == "ctrp.non.aml.heme"

  flags <- list(ohsu.flag, fimm.flag, ctrp.flag, ctrp.aml.flag, ctrp.aml2.flag, ctrp.non.aml.flag, ctrp.non.aml.heme.flag)

  genes <- gene.subset

  cat(paste0("Doing correlation of correlation for ", postfix, "\n"))
  all.data.corr <- do.corrs.of.corrs(all.expr, data.sets, genes, names, flags, postfix = postfix)
  rm(all.data.corr)

  fimm.n.samples <- ncol(fimm.expr)
  ohsu.n.samples <- ncol(ohsu.expr)
  ctrp.n.samples <- ncol(ctrp.expr)
  ctrp.aml.n.samples <- ncol(ctrp.aml.expr)
  ctrp.aml2.n.samples <- ncol(ctrp.aml2.expr)
  ctrp.non.aml.n.samples <- ncol(ctrp.non.aml.expr)
  ctrp.non.aml.heme.n.samples <- ncol(ctrp.non.aml.heme.expr)
  n.samples <- pmin(fimm.n.samples, ohsu.n.samples, ctrp.n.samples, ctrp.aml.n.samples, ctrp.non.aml.n.samples, ctrp.non.aml.heme.n.samples)
  cat(paste0("Downsampling FIMM (", fimm.n.samples, "), OHSU (", ohsu.n.samples, "), CTRP ALL (", ctrp.n.samples, ") CTRP AML (", ctrp.aml.n.samples, "), and CTRP non-AML (", ctrp.non.aml.n.samples, ") to ", n.samples, " samples\n"))

  for(i.sample in 1:10) {
    set.seed(i.sample)
    cat(paste0("Doing downsample for ", postfix, "-ds-", i.sample, "\n"))
    fimm.indices <- sample.int(n = fimm.n.samples, size = n.samples, replace = FALSE)
    ohsu.indices <- sample.int(n = ohsu.n.samples, size = n.samples, replace = FALSE)
    ctrp.indices <- sample.int(n = ctrp.n.samples, size = n.samples, replace = FALSE)
    ctrp.aml.indices <- sample.int(n = ctrp.aml.n.samples, size = n.samples, replace = FALSE)
    ctrp.aml2.indices <- sample.int(n = ctrp.aml2.n.samples, size = n.samples, replace = FALSE)
    ctrp.non.aml.indices <- sample.int(n = ctrp.non.aml.n.samples, size = n.samples, replace = FALSE)
    ctrp.non.aml.heme.indices <- sample.int(n = ctrp.non.aml.heme.n.samples, size = n.samples, replace = FALSE)

    all.expr <- cbind(fimm.expr[gene.subset, fimm.indices],
                      ohsu.expr[gene.subset, ohsu.indices],
                      ctrp.expr[gene.subset, ctrp.indices],
                      ctrp.aml.expr[gene.subset, ctrp.aml.indices],
                      ctrp.aml2.expr[gene.subset, ctrp.aml2.indices],
                      ctrp.non.aml.expr[gene.subset, ctrp.non.aml.indices],
                      ctrp.non.aml.heme.expr[gene.subset, ctrp.non.aml.heme.indices])

    data.sets <- c(rep("fimm", ncol(fimm.expr[, fimm.indices])), rep("ohsu", ncol(ohsu.expr[, ohsu.indices])), rep("ctrp.all", ncol(ctrp.expr[, ctrp.indices])), 
                   rep("ctrp.aml", ncol(ctrp.aml.expr[, ctrp.aml.indices])), rep("ctrp.aml2", ncol(ctrp.aml2.expr[, ctrp.aml2.indices])), 
                   rep("ctrp.non.aml", ncol(ctrp.non.aml.expr[, ctrp.non.aml.indices])), rep("ctrp.non.aml.heme", ncol(ctrp.non.aml.heme.expr[, ctrp.non.aml.heme.indices])))
    names <- list("ohsu", "fimm", "ctrp.all", "ctrp.aml", "ctrp.aml2", "ctrp.non.aml", "ctrp.non.aml.heme")

    fimm.flag <- data.sets == "fimm"
    ohsu.flag <- data.sets == "ohsu"
    ctrp.flag <- data.sets == "ctrp.all"
    ctrp.aml.flag <- data.sets == "ctrp.aml"
    ctrp.aml2.flag <- data.sets == "ctrp.aml2"
    ctrp.non.aml.flag <- data.sets == "ctrp.non.aml"
    ctrp.non.aml.heme.flag <- data.sets == "ctrp.non.aml.heme"

    flags <- list(ohsu.flag, fimm.flag, ctrp.flag, ctrp.aml.flag, ctrp.aml2.flag, ctrp.non.aml.flag, ctrp.non.aml.heme.flag)
  
    ds.data.corr <- do.corrs.of.corrs(all.expr, data.sets, genes, names, flags, postfix = paste0(postfix, "-ds-", i.sample))
    rm(ds.data.corr)
  }
}

## Do correlation analysis subsetting only to common genes.
do.correlation.analysis(fimm.expr = fimm.expr, ohsu.expr = ohsu.expr, ctrp.expr = ctrp.log.cpm.expr, 
                        ctrp.aml.expr = ctrp.aml.log.cpm.expr, ctrp.aml2.expr = ctrp.aml2.log.cpm.expr, 
                        ctrp.non.aml.expr = ctrp.non.aml.log.cpm.expr, ctrp.non.aml.heme.expr = ctrp.non.aml.heme.log.cpm.expr, 
                        gene.subset = common.genes, postfix = "common") 

## Do correlation analysis subsetting to genes with cov > 0.02.
do.correlation.analysis(fimm.expr = fimm.expr, ohsu.expr = ohsu.expr, ctrp.expr = ctrp.log.cpm.expr, 
                        ctrp.aml.expr = ctrp.aml.log.cpm.expr, ctrp.aml2.expr = ctrp.aml2.log.cpm.expr, 
                        ctrp.non.aml.expr = ctrp.non.aml.log.cpm.expr, ctrp.non.aml.heme.expr = ctrp.non.aml.heme.log.cpm.expr, 
                        gene.subset = high.cv.genes, postfix = "high.cv") 

stop("stop")

## In the PCA, CTRP samples are divided into 3 clusters AML/heme, a second cluster near it, and a larger third one.
## Let's find out what the second cluster is comprised of--these are those samples that have a value < 0
## in PC1.
samples <- rownames(all.combat3.pca$x)[(all.combat3.pca$x[,1]) < 0]
samples <- samples[!(samples %in% aml.expr.cols)]
samples <- samples[!(samples %in% heme.expr.cols)]
sample.metadata <- subset(ctrp.cell.line.metadata, ccl_name %in% samples)
## All of these samples are "haematopoietic_and_lymphoid_tissue"
table(sample.metadata$ccle_primary_site)
## All are them are lymphoid_neoplasm
table(sample.metadata$ccle_primary_hist)
## These are mostly lymphomomas
table(sample.metadata$ccle_hist_subtype_1)
heme.sample.metadata <- subset(ctrp.cell.line.metadata, ccl_name %in% heme.expr.cols)
## Heme are chronic and acute leukemia
table(heme.sample.metadata$ccle_hist_subtype_1)
##               acute_myeloid_leukaemia blast_phase; chronic_myeloid_leukaemia              chronic_myeloid_leukaemia 
##                                    32                                     12                                      2 
##            essential_thrombocythaemia 
##                                     1                   

aml.res <- do.combat(aml.expr, as.factor(aml.data.sets), postfix = "foc-aml")

varplot(aml.res[["pca"]], cols=as.numeric(as.factor(aml.data.sets)))
legend("top", legend = unique(as.factor(aml.data.sets)), fill = unique(as.factor(aml.data.sets)), ncol=2)
varplot(aml.res[["combat.pca"]], cols=as.numeric(as.factor(aml.data.sets)))
legend("top", legend = unique(as.factor(aml.data.sets)), fill = unique(as.factor(aml.data.sets)), ncol=2)

varplot2(aml.res[["pca"]], cols=as.numeric(as.factor(aml.data.sets)))
legend("top", legend = unique(as.factor(aml.data.sets)), fill = unique(as.factor(aml.data.sets)), ncol=2)
varplot2(aml.res[["combat.pca"]], cols=as.numeric(as.factor(aml.data.sets)))
legend("top", legend = unique(as.factor(aml.data.sets)), fill = unique(as.factor(aml.data.sets)), ncol=2)


all.expr <- cbind(fimm.expr[common.genes, ],
                  ohsu.expr[common.genes, ],
                  ctrp.log.cpm.expr[common.genes, ],
                  ctrp.aml.log.cpm.expr[common.genes, ],
                  ctrp.aml2.log.cpm.expr[common.genes, ],
                  ctrp.non.aml.log.cpm.expr[common.genes, ],
                  ctrp.non.aml.heme.log.cpm.expr[common.genes, ])

data.sets <- c(rep("fimm", ncol(fimm.expr)), rep("ohsu", ncol(ohsu.expr)), rep("ctrp.all", ncol(ctrp.log.cpm.expr)), rep("ctrp.aml", ncol(ctrp.aml.log.cpm.expr)), rep("ctrp.aml2", ncol(ctrp.aml2.log.cpm.expr)), 
               rep("ctrp.non.aml", ncol(ctrp.non.aml.log.cpm.expr)), rep("ctrp.non.aml.heme", ncol(ctrp.non.aml.heme.log.cpm.expr)))
names <- list("ohsu", "fimm", "ctrp.all", "ctrp.aml", "ctrp.aml2", "ctrp.non.aml", "ctrp.non.aml.heme")

fimm.flag <- data.sets == "fimm"
ohsu.flag <- data.sets == "ohsu"
ctrp.aml.flag <- data.sets == "ctrp.aml"
ctrp.aml2.flag <- data.sets == "ctrp.aml2"
ctrp.non.aml.flag <- data.sets == "ctrp.non.aml"
ctrp.non.aml.heme.flag <- data.sets == "ctrp.non.aml.heme"

flags <- list(ohsu.flag, fimm.flag, ctrp.aml.flag, ctrp.aml2.flag, ctrp.non.aml.flag, ctrp.non.aml.heme.flag)

genes <- common.genes

all.data.corr <- do.corrs.of.corrs(all.expr, data.sets, genes, names, flags, postfix = "all-data-")

## Downsample to the size of the smallest data set
fimm.n.samples <- ncol(fimm.expr)
ohsu.n.samples <- ncol(ohsu.expr)
ctrp.aml.n.samples <- ncol(ctrp.aml.log.cpm.expr)
ctrp.aml2.n.samples <- ncol(ctrp.aml2.log.cpm.expr)
ctrp.non.aml.n.samples <- ncol(ctrp.non.aml.log.cpm.expr)
ctrp.non.aml.heme.n.samples <- ncol(ctrp.non.aml.heme.log.cpm.expr)
n.samples <- pmin(fimm.n.samples, ohsu.n.samples, ctrp.aml.n.samples, ctrp.non.aml.n.samples, ctrp.non.aml.heme.n.samples)
cat(paste0("Downsampling FIMM (", fimm.n.samples, "), OHSU (", ohsu.n.samples, "), CTRP AML (", ctrp.aml.n.samples, "), and CTRP non-AML (", ctrp.non.aml.n.samples, ") to ", n.samples, " samples\n"))
n.genes <- 5000

corrs <- list()
for(i.sample in 1:1) {
  set.seed(i.sample)

  fimm.indices <- sample.int(n = fimm.n.samples, size = n.samples, replace = FALSE)
  ohsu.indices <- sample.int(n = ohsu.n.samples, size = n.samples, replace = FALSE)
  ctrp.aml.indices <- sample.int(n = ctrp.aml.n.samples, size = n.samples, replace = FALSE)
  ctrp.aml2.indices <- sample.int(n = ctrp.aml2.n.samples, size = n.samples, replace = FALSE)
  ctrp.non.aml.indices <- sample.int(n = ctrp.non.aml.n.samples, size = n.samples, replace = FALSE)
  ctrp.non.aml.heme.indices <- sample.int(n = ctrp.non.aml.heme.n.samples, size = n.samples, replace = FALSE)

  do.downsample <- FALSE
  if(!do.downsample) {
    fimm.indices <- 1:fimm.n.samples
    ohsu.indices <- 1:ohsu.n.samples
    ctrp.aml.indices <- 1:ctrp.aml.n.samples
    ctrp.aml2.indices <- 1:ctrp.aml2.n.samples
    ctrp.non.aml.indices <- 1:ctrp.non.aml.n.samples
    ctrp.non.aml.heme.indices <- 1:ctrp.non.aml.heme.n.samples
  }

  all.expr <- cbind(fimm.expr[common.genes, fimm.indices],
                    ohsu.expr[common.genes, ohsu.indices],
                    ctrp.aml.log.cpm.expr[common.genes, ctrp.aml.indices],
                    ctrp.aml2.log.cpm.expr[common.genes, ctrp.aml2.indices],
                    ctrp.non.aml.log.cpm.expr[common.genes, ctrp.non.aml.indices],
                    ctrp.non.aml.heme.log.cpm.expr[common.genes, ctrp.non.aml.heme.indices])

  data.sets <- c(rep("fimm", n.samples), rep("ohsu", n.samples), rep("ctrp.aml", n.samples), rep("ctrp.aml2", n.samples), rep("ctrp.non.aml", n.samples), rep("ctrp.non.aml.heme", n.samples))

  all.expr <- cbind(fimm.expr[common.genes, fimm.indices],
                    ohsu.expr[common.genes, ohsu.indices])

  data.sets <- c(rep("fimm", n.samples), rep("ohsu", n.samples))
  
  study.expressed.genes <- llply(unique(data.sets), .fun = function(data.set) {
                                   flag <- data.sets == data.set
                                   avg.log.cpm <- unlist(apply(all.expr[, flag], 1, function(row) mean(row, na.rm=TRUE)))
                                   expressed.genes <- rownames(all.expr)[avg.log.cpm > 0]
                                   expressed.genes
                                 })
  expressed.genes <- Reduce(intersect, study.expressed.genes)

  study.variable.genes <- llply(unique(data.sets), .fun = function(data.set) { 
                                                     flag <- data.sets == data.set
                                                     sd.log.cpm <- unlist(apply(all.expr[, flag], 1, function(row) sd(row, na.rm=TRUE)))
                                                     var.genes <- rownames(all.expr); var.genes <- var.genes[order(sd.log.cpm, decreasing=TRUE)]
                                                     var.genes[1:min(length(var.genes),n.genes)]
                                                   })
  variable.genes <- Reduce(intersect, study.variable.genes)

  fimm.flag <- data.sets == "fimm"
  ohsu.flag <- data.sets == "ohsu"
  ctrp.aml.flag <- data.sets == "ctrp.aml"
  ctrp.aml2.flag <- data.sets == "ctrp.aml2"
  ctrp.non.aml.flag <- data.sets == "ctrp.non.aml"
  ctrp.non.aml.heme.flag <- data.sets == "ctrp.non.aml.heme"

  cat(paste0("Doing downsample iteration ", i.sample, "\n"))

  ## corrs[[length(corrs)+1]] <- do.corrs.of.corrs(all.expr, data.sets, common.genes, ohsu.flag, fimm.flag, ctrp.flag, postfix = paste0("ds-", i.sample))
  ## corrs[[length(corrs)+1]] <- do.corrs.of.corrs(all.expr, data.sets, common.genes, c("ohsu", "fimm", "ctrp"), c(ohsu.flag, fimm.flag, ctrp.flag), postfix = paste0("ds-", i.sample))
  genes <- intersect(common.genes, expressed.genes)
  genes <- intersect(genes, variable.genes)
  print(length(genes))
##  genes <- intersect(genes, cancer.genes)
  corrs[[length(corrs)+1]] <- do.corrs.of.corrs(all.expr, data.sets, genes, list("ohsu", "fimm", "ctrp.aml", "ctrp.aml2", "ctrp.non.aml", "ctrp.non.aml.heme"), 
                                                list(ohsu.flag, fimm.flag, ctrp.aml.flag, ctrp.aml2.flag, ctrp.non.aml.flag, ctrp.non.aml.heme.flag), postfix = paste0("ds-", i.sample))
}

stop("stop")


inter <- intersect(intersect(rownames(fimm.expr), rownames(ohsu.expr)), rownames(ctrp.log.cpm.expr))
all.expr <- cbind(fimm.expr[inter,],
                  ohsu.expr[inter,],
                  ctrp.log.cpm.expr[inter, ])
cat(paste0("all.expr = ", nrow(all.expr) , " x ", ncol(all.expr), "\n"))
data.sets <- c(rep("fimm", ncol(fimm.expr)), rep("ohsu", ncol(ohsu.expr)), rep("ctrp", ncol(ctrp.log.cpm.expr[inter, ])))

all.expr.combat <- do.combat(all.expr, data.sets)
plot(apply(all.expr.combat, 1, sd), apply(all.expr.combat, 1, mean))

stop("stop")

batch           <- as.factor(data.sets)
modcombat       <- model.matrix(~1,data=batch)

all.expr.combat  <- ComBat(all.expr, batch = batch, mod=modcombat)
offset <- min(all.expr.combat)
if(offset < 0) {
  all.expr.combat <- all.expr.combat - offset + 0.001
}

q()

nq.expr <- normalize.quantiles(as.matrix(all.expr))

nq.combat  <- ComBat(nq.expr, batch = batch, mod=modcombat)
offset <- min(nq.combat)
if(offset < 0) {
  nq.combat <- nq.combat - offset + 0.001
}

nq <- normalize.quantiles(as.matrix(all.expr.combat))

offset <- min(all.expr)
if(offset < 0) {
  all.expr <- all.expr - offset + 0.001
}


combat.covs <- list()
covs <- list()
nq.covs <- list()
nq.combat.covs <- list()
for(ds in c("fimm", "ohsu", "ctrp")) {
  flag <- data.sets == ds
  combat.covs[[ds]] <- apply(all.expr.combat[, flag], 1, sd) / apply(all.expr.combat[, flag], 1, mean) 
  covs[[ds]] <- apply(all.expr[, flag], 1, sd) / apply(all.expr[, flag], 1, mean) 
  nq.covs[[ds]] <- apply(nq[, flag], 1, sd) / apply(nq[, flag], 1, mean) 
  nq.combat.covs[[ds]] <- apply(nq.combat[, flag], 1, sd) / apply(nq.combat[, flag], 1, mean) 
}

if(FALSE) {
smoothScatter(covs[["fimm"]], covs[["ctrp"]])
smoothScatter(combat.covs[["fimm"]], combat.covs[["ctrp"]])
smoothScatter(nq.covs[["fimm"]], nq.covs[["ctrp"]])
smoothScatter(nq.combat.covs[["fimm"]], nq.combat.covs[["ctrp"]])


## To see mu vs sd plot
mean <- unlist(apply(ohsu.expr[common.genes,], 1, mean, na.rm=TRUE))
sd <- unlist(apply(ohsu.expr[common.genes,], 1, sd, na.rm=TRUE))
smoothScatter(mean, sd)
}

fimm.flag <- data.sets == "fimm"
ohsu.flag <- data.sets == "ohsu"
ctrp.flag <- data.sets == "ctrp"
flag <- data.sets %in% c("fimm", "ohsu")

if(FALSE) {
mean <- unlist(apply(all.expr[, flag], 1, mean, na.rm=TRUE))
sd <- unlist(apply(all.expr[, flag], 1, sd, na.rm=TRUE))
smoothScatter(mean, sd)
}

if(FALSE) {
plot(density(unlist(all.expr.combat[,fimm.flag])))
lines(density(unlist(all.expr.combat[,ohsu.flag])))
lines(density(unlist(all.expr.combat[,ctrp.flag])))

plot(density(unlist(all.expr[,fimm.flag])))
lines(density(unlist(all.expr[,ohsu.flag])))
lines(density(unlist(all.expr[,ctrp.flag])))

plot(density(unlist(nq.expr[,fimm.flag])))
lines(density(unlist(nq.expr[,ohsu.flag])))
lines(density(unlist(nq.expr[,ctrp.flag])))

all.expr <- cbind(fimm.expr[common.genes,] - min(fimm.expr[common.genes,]) + 0.001, 
                  ohsu.expr[common.genes,] - min(ohsu.expr[common.genes,]) + 0.001, 
                  ctrp.log.cpm.expr[common.genes, aml.expr.cols] - min(ctrp.log.cpm.expr[common.genes,]) + 0.001)
data.sets <- c(rep("fimm", ncol(fimm.expr)), rep("ohsu", ncol(ohsu.expr)), rep("ctrp", ncol(ctrp.log.cpm.expr[common.genes, aml.expr.cols])))

nq <- normalize.quantiles(as.matrix(all.expr))

covs <- list()
nq.covs <- list()
for(ds in c("fimm", "ohsu", "ctrp")) {
  flag <- data.sets == ds
  covs[[ds]] <- apply(all.expr[, flag], 1, sd) / apply(all.expr[, flag], 1, mean) 
  nq.covs[[ds]] <- apply(nq[, flag], 1, sd) / apply(nq[, flag], 1, mean) 
}

lmf <- lm(nq.covs[["fimm"]] ~ nq.covs[["ctrp"]])
summary(lmf)

plot(covs[["fimm"]], covs[["ctrp"]])
plot(nq.covs[["fimm"]], nq.covs[["ctrp"]])
}  ## if(FALSE)

## expr.mats <- list("ohsu" = all.expr[common.genes, ohsu.flag], "fimm" = all.expr[common.genes, fimm.flag], "ctrp" = all.expr[common.genes, ctrp.flag])
expr.mats <- list("ohsu" = all.expr.combat[common.genes, ohsu.flag], "fimm" = all.expr.combat[common.genes, fimm.flag], "ctrp" = all.expr.combat[common.genes, ctrp.flag])

intraCors <- list()
for(i in 1:length(expr.mats))
{
  print(i)
  tempEx      <- expr.mats[[i]]; 
  intraCors[[i]] <- fastCor(tempEx, Split = 1500,Method = "spearman")
}

names(intraCors) <- names(expr.mats)

png("intraStudyCors.png")
## plot(density(fisherz(unlist(as.dist(intraCors[[1]]))),na.rm=T), main = "Intra Study Correlations", xlab="Fisher z transform of\nPearson Correlations", ylim= c(0,2.5),xlim=c(-1,5))
plot(density(fisherz(unlist(as.dist(intraCors[[1]]))),na.rm=T), main = "Intra Study Correlations", xlab="Fisher z transform of\nPearson Correlations")
labels <- names(intraCors)
n.labels <- length(labels)
legend(2.5,2.5,legend=labels, fill=1:n.labels,border=1:n.labels, box.col="white")
for(i in 2:length(intraCors)){lines(density(fisherz(unlist(as.dist(intraCors[[i]]))),na.rm=T), col = i); print(i)}
dev.off()

corMat     <- matrix(NA, choose(length(common.genes),2), length(expr.mats))

for(i in 1:length(intraCors)) {
  print(i)
  corMat[,i]  <- as.vector(as.dist(intraCors[[i]]))
}
rm(intraCors)
gc()

save.image(".Rdata")

# Use Fisher's z transform to remove hard (-1,1) bounds on correlation distributions
corMat <- apply(corMat,2, fisherz)

corMat[corMat == Inf]  <- NA
cors = cor(corMat, use = "pairwise")
rm(corMat);gc()
colnames(cors) <- labels
rownames(cors) <- labels

cors

write.table(cors, file="corOfCors.txt",sep="\t", col.names=T, row.names=T,quote=F)


if(FALSE) {
ohsu.expr.filt <- ohsu.expr.filt[common.genes,]
ctrp.expr.filt <- ctrp.expr.filt[common.genes,]
ctrp.log2.expr.filt <- ctrp.log2.expr.filt[common.genes,]
fimm.expr.filt <- fimm.expr.filt[common.genes,]

ctrp.expr.un <- unlist(ctrp.expr.filt)
ctrp.log2.expr.un <- unlist(ctrp.log2.expr.filt)
ohsu.expr.un <- unlist(ohsu.expr.filt)
fimm.expr.un <- unlist(fimm.expr.filt)

plot(density(ctrp.expr.un)); lines(density(ohsu.expr.un)); lines(density(fimm.expr.un)); lines(density(ctrp.log2.expr.un))
}