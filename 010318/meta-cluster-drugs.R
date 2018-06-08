rm(list = ls())

suppressPackageStartupMessages(library("foreach"))
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library(tidyr))

suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("gridExtra"))
suppressPackageStartupMessages(library("caTools"))
suppressPackageStartupMessages(library("data.table"))

suppressPackageStartupMessages(library(ggbeeswarm))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(synapseClient))
suppressPackageStartupMessages(library(mygene))
suppressPackageStartupMessages(library(openxlsx))

library(rJava)
library(rcdk)
library(fingerprint)
library(plyr)
library(tidyverse)

library(dynamicTreeCut)
library(MCL)
library(clusteval)

library(WeightedCluster)
library(igraph)

library("dendextend")
## library("dendextendRcpp")


## Log in to Synapse
synapseLogin()

## Bring in common code
source("../common/data-preprocessing.R")
source("../common/models.R")
source("../common/plotting.R")
source("../common/utils.R")
source("heatmap-utils.R")

## These are defnied in heatmap-utils
orig.ohsu.metadata <- ohsu.metadata
orig.fimm.metadata <- fimm.metadata




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

if(TRUE) {
## Overwrite the expression synIds to exclude outliers (as before), but _not_ to restrict to highly variable/cancer genes
## Batch-corrected log2 cpm (or rma) expression files (NB: these exclude 2 outliers in OHSU)
data.set.expr.synIds <- list("ohsu" = "syn11362256", "fimm" = "syn11362257")
## names(data.set.expr.synIds) <- data.sets
expr.folder.synIds <- list("ohsu" = "syn11361089", "fimm" = "syn11361089")
data.set.expr.files <- list("ohsu" = "ohsu-expr-fimm-ohsu-outlier1-combat.tsv", "fimm" = "fimm-expr-fimm-ohsu-outlier1-combat.tsv")
}

rdata.file <- ".Rdata.model.mean.response.120817"

## Over-write some of the variables
## These fit files are assumed to be for AUCs/DSS values calculated across shared concentration ranges across data sets.
## These files may be created by a script such as
## calculate-dss-fimm-ohsu.R
## data.set.dss.fit.synIds <- c("syn11362885", "syn11362891")
## names(data.set.dss.fit.synIds) <- data.sets
## data.set.dss.fit.files <- c("ohsu-fimm-ohsu-outlier1-dss-t0.tsv", "fimm-fimm-ohsu-outlier1-dss-t0.tsv")

## A string that will be included in any output files.
file.prefix <- "fimm-ohsu-validate-trained-models-112217-one-ohsu-no-filtering-all-mean-response"

## End setup

## Load the map of drugs between data sets
obj <- synGet(drug.mapping.synId, downloadFile = TRUE)
drug.map <- read.table(getFileLocation(obj), sep="\t", header=TRUE, quote="\"", comment.char="")

## Drop any drugs that show up in multiple rows of the drug map--which would indicate an ambiguous mapping.
my.dup <- function(x) duplicated(x, fromLast = TRUE) | duplicated(x, fromLast = FALSE)
for(col in drug.map.drug.id.cols) {
  flag <- !my.dup(drug.map[, col])
  drug.map <- drug.map[flag, ]
}

subset.and.zscore.matrix <- function(df, row.frac.cutoff = 0.25, col.frac.cutoff = 0.25) {
##  frac.na <- 0.25
##  cutoff <- frac.na * nrow(df)
##  cutoff <- 20
  tmp <- df
  flag <- unlist(apply(tmp, 1, function(row) all(is.na(row))))
  tmp <- tmp[!flag,]
  flag <- unlist(apply(tmp, 1, function(row) length(which(!is.na(row))) < ( row.frac.cutoff * length(row) )))
  tmp <- tmp[!flag,]
  flag <- unlist(apply(tmp, 2, function(col) all(is.na(col))))
  tmp <- tmp[, !flag]
  flag <- unlist(apply(tmp, 2, function(col) length(which(!is.na(col))) < ( col.frac.cutoff * length(col) )))
  tmp <- tmp[, !flag]
  tmp.scaled <- t(scale(t(tmp)))
  tmp.scaled
}

safe.t.test <- function(...) {
  res <- tryCatch({t.test(...)}, error = function(e) { data.frame(p.value = 1) })
  res
}

process.fits <- function(fits, drug.map) {
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


  ## Drop any fits that are flagged as outliers
  if(TRUE) {
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
  fits
}

orig.fits <- llply(data.set.drug.fit.synIds, .parallel = FALSE,
              .fun = function(synId) {
                       obj <- synGet(synId, downloadFile = TRUE)
                       file <- getFileLocation(obj)
                       read.table(file, sep="\t", header=TRUE, as.is=TRUE, comment.char="", quote="\"")
              })

fits <- llply(data.set.dss.fit.synIds, .parallel = FALSE,
              .fun = function(synId) {
                       obj <- synGet(synId, downloadFile = TRUE)
                       file <- getFileLocation(obj)
                       read.table(file, sep="\t", header=TRUE, as.is=TRUE, comment.char="", quote="\"")
              })


orig.fits <- process.fits(orig.fits, drug.map)
fits <- process.fits(fits, drug.map)

##source("../common/drug-mapping-functions/mapping_helpers.R")
##source("../common/drug-mapping-functions/drug_mapping.R")

##converts SMILES string to fingerprint
parseInputFingerprint <- function(input) {
  input.mol <- parse.smiles(input)
  lapply(input.mol, do.typing)
  lapply(input.mol, do.aromaticity)
  lapply(input.mol, do.isotopes)
  fp.inp <- lapply(input.mol, get.fingerprint, type = "extended")
}

drug.structures.tbl <- read.table("../common/drug-mapping-functions/drug_structures.txt", comment.char="", sep="\t", header=TRUE, as.is=TRUE)
drug.structures.tbl <- merge(drug.structures.tbl, drug.map[, c("OHSU_DRUG_NAME", "FIMM_DRUG_NAME")], by.x = "drug", by.y = "FIMM_DRUG_NAME")

cat("Calculating drug fingerprints for common drugs\n")
fp.common <- parseInputFingerprint(as.character(drug.structures.tbl$smiles))
rownames(drug.structures.tbl) <-  drug.structures.tbl$smiles
names(fp.common) <- drug.structures.tbl[names(fp.common), "OHSU_DRUG_NAME"]

safe.cor.test <- function(...) {
  ct <- tryCatch({cor.test(...)}, error = function(e) { data.frame(estimate = NA, p.value = NA) })
  ct
}

plot.drug.heatmap <- function(mat, drug.tbl, show.col.names = TRUE, show.row.names = TRUE) {
  if(!show.col.names) { colnames(mat) <- NULL }
  if(!show.row.names) { rownames(mat) <- NULL }
  hm <- Heatmap(mat, clustering_distance_rows = function(x) na.dist(x), clustering_distance_columns = function(x) na.dist(x), na_col = "black", cluster_rows = FALSE)
  row.annotation <- data.frame(cluster = as.character(drug.tbl$cluster))
  rownames(row.annotation) <- drug.tbl$OHSU_DRUG_NAME
  cols <- rainbow(length(unique(row.annotation$cluster)))
  names(cols) <- unique(row.annotation$cluster)
  row.col.list <- list("cluster" = cols)
  ra <- rowAnnotation(df = row.annotation, col = row.col.list, show_annotation_name = TRUE)
  hm <- hm + ra
  hm
}

## Plot correlation of each drug with GLDS
plot.drug.correlation.with.glds <- function(mat) {
  glds <- colMeans(mat, na.rm=TRUE)
  indices <- 1:nrow(mat)
  names(indices) <- row.names(mat)
  drug.glds.corr <- ldply(indices,
                          .fun = function(i) {
                                   drug.resp <- mat[i, ]
                                   ct <- safe.cor.test(drug.resp, glds)
                                   data.frame(corr = ct$estimate, pval = ct$p.value)
                          })
  colnames(drug.glds.corr) <- c("drug", "corr", "pval")
  drug.glds.corr <- na.omit(drug.glds.corr)
  drug.glds.corr <- drug.glds.corr[order(drug.glds.corr$corr, decreasing=FALSE), ]

  g <- ggplot()
  g <- g + geom_point(data = drug.glds.corr, aes(x = corr, y = -log10(pval)))
  df <- subset(drug.glds.corr, pval > 0.1)
  df <- subset(drug.glds.corr, corr < 0)
  df <- drug.glds.corr[1:min(5, nrow(drug.glds.corr)),]
  df$drug <- unlist(lapply(df$drug, function(str) gsub(str, pattern="([^(]+)[ ]*\\(.*", replacement="\\1")))
  g <- g + xlab("Pearson Correlation")
  g <- g + xlim(c(min(min(drug.glds.corr$corr, na.rm=TRUE) - 0.1, -0.2), max(drug.glds.corr$corr, na.rm=TRUE) + 0.1))
  g <- g + geom_text(data = df, aes(x = corr, y = -log10(pval), label = drug), position = position_jitter(width=0, height=1), hjust = 1)
  list("g" = g, "drug.glds.corr" = drug.glds.corr)
}

## compares every smiles on list one across list 2
cat("Calculating all pairwise distances between drugs\n")
drug.struct.dist <- ldply(fp.common, .fun = function(i) {
  sim <- sapply(fp.common, function(j) {
    1 - distance(i, j, method = "tanimoto")
  })
  sim
})

rownames(drug.struct.dist) <- drug.struct.dist$.id
drug.struct.dist <- drug.struct.dist[, !(colnames(drug.struct.dist) == ".id")]

orig.ohsu.drc <- prepare.drug.response.matrix(orig.fits[["ohsu"]], drug.col = "inhibitor", patient.col = "patient_id", response.col = "auc.ll4")
orig.ohsu.drc.scaled <- subset.and.zscore.matrix(orig.ohsu.drc, col.frac.cutoff = 0.25, row.frac.cutoff = 0.25)
orig.ohsu.drug.response.dist <- as.matrix(na.dist(orig.ohsu.drc.scaled))

## tmp <- merge(orig.fits[["fimm"]], drug.map[, c("FIMM_DRUG_NAME", "OHSU_DRUG_NAME")], by.x = "DRUG_NAME", by.y = "FIMM_DRUG_NAME")
orig.fimm.drc <- prepare.drug.response.matrix(orig.fits[["fimm"]], drug.col = "OHSU_DRUG_NAME", patient.col = "PATIENT_ID", response.col = "auc.ll4")
orig.fimm.drc.scaled <- subset.and.zscore.matrix(orig.fimm.drc, col.frac.cutoff = 0.25, row.frac.cutoff = 0.25)
orig.fimm.drug.response.dist <- as.matrix(na.dist(orig.fimm.drc.scaled))

ohsu.drc <- prepare.drug.response.matrix(fits[["ohsu"]], drug.col = "inhibitor", patient.col = "patient_id", response.col = "auc.ll4")
ohsu.drc.scaled <- subset.and.zscore.matrix(ohsu.drc, col.frac.cutoff = 0.25, row.frac.cutoff = 0.25)
ohsu.drug.response.dist <- as.matrix(na.dist(ohsu.drc.scaled))

## tmp <- merge(fits[["fimm"]], drug.map[, c("FIMM_DRUG_NAME", "OHSU_DRUG_NAME")], by.x = "DRUG_NAME", by.y = "FIMM_DRUG_NAME")
fimm.drc <- prepare.drug.response.matrix(fits[["fimm"]], drug.col = "OHSU_DRUG_NAME", patient.col = "PATIENT_ID", response.col = "auc.ll4")
fimm.drc.scaled <- subset.and.zscore.matrix(fimm.drc, col.frac.cutoff = 0.25, row.frac.cutoff = 0.25)
fimm.drug.response.dist <- as.matrix(na.dist(fimm.drc.scaled))

orig.common.drugs <- Reduce("intersect", lapply(list(orig.ohsu.drc.scaled, orig.fimm.drc.scaled), function(mat) rownames(mat)))
common.drugs <- Reduce("intersect", lapply(list(ohsu.drc.scaled, fimm.drc.scaled), function(mat) rownames(mat)))

pdf("ohsu-all-drug-glds-corr.pdf")
ret <- plot.drug.correlation.with.glds(orig.ohsu.drc.scaled)
g <- ret$g
g <- g + ggtitle("OHSU: (All) Drug correlation with GLDS")
print(g)
d <- dev.off()

pdf("ohsu-common-drug-glds-corr.pdf")
ret <- plot.drug.correlation.with.glds(orig.ohsu.drc.scaled[orig.common.drugs, ])
g <- ret$g
common.ohsu.glds.corr <- ret$drug.glds.corr
g <- g + ggtitle("OHSU: (Common) Drug correlation with GLDS")
print(g)
d <- dev.off()

pdf("fimm-all-drug-glds-corr.pdf")
ret <- plot.drug.correlation.with.glds(orig.fimm.drc.scaled)
g <- ret$g
g <- g + ggtitle("FIMM: (All) Drug correlation with GLDS")
print(g)
d <- dev.off()

pdf("fimm-common-drug-glds-corr.pdf")
ret <- plot.drug.correlation.with.glds(orig.fimm.drc.scaled[orig.common.drugs, ])
g <- ret$g
common.fimm.glds.corr <- ret$drug.glds.corr
g <- g + ggtitle("FIMM: (Common) Drug correlation with GLDS")
print(g)
d <- dev.off()

common.glds.corr <- merge(common.ohsu.glds.corr, common.fimm.glds.corr, by = "drug", suffixes = c(".ohsu", ".fimm"))
g <- ggplot(data = common.glds.corr)
g <- g + geom_point(aes(x = corr.ohsu, y = corr.fimm))
g <- g + xlab("OHSU Pearson") + ylab("FIMM Pearson")
df <- subset(common.glds.corr, (corr.ohsu < 0) | (corr.fimm < 0))
df$drug <- unlist(lapply(df$drug, function(str) gsub(str, pattern="([^(]+)[ ]*\\(.*", replacement="\\1")))
df2 <- data.frame(x = common.glds.corr$corr.ohsu, y = common.glds.corr$corr.fimm)
g <- g + geom_text(data = df, aes(x = corr.ohsu, y = corr.fimm, label = drug), position = position_jitter(width=0, height=0), hjust = 0)
g <- g + geom_smooth(data = df2, aes(x = x, y = y), method='lm')
g <- g + geom_text(x = 0, y = 0.75, label = lm_eqn(df2), parse=TRUE)

pdf("fimm-ohsu-correlations-with-glds.pdf")
print(g)
d <- dev.off()

common <- intersect(rownames(drug.struct.dist), rownames(orig.ohsu.drug.response.dist))
tmp1 <- orig.ohsu.drug.response.dist[common, common]
tmp2 <- drug.struct.dist[common, common]

x <- tmp1[upper.tri(tmp1, diag=FALSE)]
y <- tmp2[upper.tri(tmp2, diag=FALSE)]
ct <- cor.test(x, y)

common <- intersect(rownames(orig.fimm.drug.response.dist), rownames(orig.ohsu.drug.response.dist))
tmp1 <- orig.ohsu.drug.response.dist[common, common]
tmp2 <- orig.fimm.drug.response.dist[common, common]

x <- tmp1[upper.tri(tmp1, diag=FALSE)]
y <- tmp2[upper.tri(tmp2, diag=FALSE)]
ct <- cor.test(x, y)

## pdf(paste0(file.prefix, "-ohsu-all-heatmap.pdf"))
## plot.heatmap(tmp.scaled, ohsu.metadata, drug.metadata)
## d <- dev.off()

mat <- orig.ohsu.drc.scaled
hm <- Heatmap(mat, clustering_distance_rows = function(x) na.dist(x), clustering_distance_columns = function(x) na.dist(x), na_col = "black")
ohsu.row.order <- rownames(mat)[row_order(hm)[[1]]]
ohsu.col.order <- colnames(mat)[column_order(hm)]

mat <- orig.fimm.drc.scaled
hm <- Heatmap(mat, clustering_distance_rows = function(x) na.dist(x), clustering_distance_columns = function(x) na.dist(x), na_col = "black")
fimm.row.order <- rownames(mat)[row_order(hm)[[1]]]
fimm.col.order <- colnames(mat)[column_order(hm)]

mats <- list("orig-ohsu-common-drugs" = orig.ohsu.drc.scaled[orig.common.drugs,], 
             "ohsu-common-drugs" = ohsu.drc.scaled[common.drugs,],
             "orig-ohsu-all-drugs" = orig.ohsu.drc.scaled,
             "orig-fimm-common-drugs" = orig.fimm.drc.scaled[orig.common.drugs,], 
             "fimm-common-drugs" = fimm.drc.scaled[common.drugs,],
             "orig-fimm-all-drugs" = orig.fimm.drc.scaled,
             "orig-fimm-ohsu-order-common-drugs" = orig.fimm.drc.scaled[orig.common.drugs,], 
             "fimm-ohsu-order-common-drugs" = fimm.drc.scaled[common.drugs,],
             "orig-fimm-ohsu-order-all-drugs" = orig.fimm.drc.scaled)

mat.rows <- list("orig-ohsu-common-drugs" = ohsu.row.order,
             "ohsu-common-drugs" = ohsu.row.order,
             "orig-ohsu-all-drugs" = ohsu.row.order,
             "orig-fimm-common-drugs" = fimm.row.order,
             "fimm-common-drugs" = fimm.row.order,
             "orig-fimm-all-drugs" = fimm.row.order,
             "orig-fimm-ohsu-order-common-drugs" = ohsu.row.order, 
             "fimm-ohsu-order-common-drugs" = ohsu.row.order,
             "orig-fimm-ohsu-order-all-drugs" = ohsu.row.order)


mat.cols <- list("orig-ohsu-common-drugs" = ohsu.col.order,
             "ohsu-common-drugs" = ohsu.col.order,
             "orig-ohsu-all-drugs" = ohsu.col.order,
             "orig-fimm-common-drugs" = fimm.col.order,
             "fimm-common-drugs" = fimm.col.order,
             "orig-fimm-all-drugs" = fimm.col.order,
             "orig-fimm-ohsu-order-common-drugs" = fimm.col.order, 
             "fimm-ohsu-order-common-drugs" = fimm.col.order,
             "orig-fimm-ohsu-order-all-drugs" = fimm.col.order)


for(nm in names(mats)) {
  pdf(paste0(nm, "-heatmap.pdf"))
  mat <- mats[[nm]]
##  hm <- Heatmap(mat, clustering_distance_rows = function(x) na.dist(x), clustering_distance_columns = function(x) na.dist(x), na_col = "black")
  rows <- mat.rows[[nm]][mat.rows[[nm]] %in% rownames(mat)]
  cols <- mat.cols[[nm]][mat.cols[[nm]] %in% colnames(mat)]
  mat <- mat[rows, cols]
  colnames(mat) <- NULL
  hm <- Heatmap(mat, cluster_rows = FALSE, cluster_columns = FALSE, na_col = "black")
  print(hm)
  d <- dev.off()
}

for(nm in c("ohsu-common-drugs", "fimm-common-drugs")) {
  pdf(paste0(nm, "-hc-heatmap.pdf"))
  mat <- mats[[nm]]
  colnames(mat) <- NULL
  hm <- Heatmap(mat, clustering_distance_rows = function(x) na.dist(x), clustering_distance_columns = function(x) na.dist(x), na_col = "black")
  print(hm)
  d <- dev.off()
}

plot.drug.clusters <- function(mat, k.char, drug.map, hc.file = NULL) {
  ## ct <- cutreeDynamic(hc, distM = as.matrix(na.dist(mat)), minClusterSize = 1)
  ## df <- data.frame(drug = labels(hc), cluster = ct)
  ## df <- data.frame(drug = labels(hc), cluster = ct)
  ## Plot the tree to see manually where we should make the cut
  ## plot(hc)
  if(!is.null(hc.file)) { pdf(hc.file) }
  dend <- mat %>% na.dist %>% hclust %>% as.dendrogram
  # plot + color the dend's branches before, based on 3 clusters:
  ## dend %>% color_branches(k=3) %>% plot(horiz=FALSE)
  dend %>% plot(horiz=FALSE)
  ## add horiz rect
  #dend %>% rect.dendrogram(k=3,horiz=TRUE)
  ## add horiz (well, vertical) line:
  #abline(v = heights_per_k.dendrogram(dend)["3"] + .6, lwd = 2, lty = 2, col = "blue")
  h <- unname(heights_per_k.dendrogram(dend)[k.char])
  abline(h = h, lwd = 2, lty = 2, col = "blue")
  if(!is.null(hc.file)) { d <- dev.off() }
  hc <- hclust(na.dist(mat))
  ct <- cutree(hc, k=as.numeric(k.char))
  df <- data.frame(drug = names(ct), cluster = unname(ct))
  drug.clusters <- dlply(df, .variables = "cluster", .fun = function(tbl) as.character(tbl$drug))
  drug.tbl <- merge(drug.map, df, by.x = "OHSU_DRUG_NAME", by.y = "drug")
  rownames(mat) <- unlist(lapply(rownames(mat), function(str) gsub(str, pattern="([^(]+)[ ]+\\(.*", replacement="\\1")))
  drug.tbl$OHSU_DRUG_NAME <- unlist(lapply(drug.tbl$OHSU_DRUG_NAME, function(str) gsub(str, pattern="([^(]+)[ ]+\\(.*", replacement="\\1")))
  rownames(drug.tbl) <- as.character(drug.tbl$OHSU_DRUG_NAME)
  common <- intersect(rownames(mat), rownames(drug.tbl))
  drug.tbl <- drug.tbl[common,]
  drug.tbl <- drug.tbl[order(drug.tbl$cluster),]                                       
  common <- rownames(drug.tbl)
  mat <- mat[common, ]
  hm <- plot.drug.heatmap(mat, drug.tbl, show.col.names = FALSE)
  list("heatmap" = hm, "clusters" = drug.clusters)
}


mat <- orig.ohsu.drc.scaled[orig.common.drugs,]
hc <- hclust(na.dist(mat))
ct <- cutreeDynamic(hc, distM = as.matrix(na.dist(mat)), minClusterSize = 1)
df <- data.frame(drug = labels(hc), cluster = ct)
orig.ohsu.drug.clusters <- dlply(df, .variables = "cluster", .fun = function(tbl) as.character(tbl$drug))

mat <- orig.fimm.drc.scaled[orig.common.drugs,]
hc <- hclust(na.dist(mat))
## ct <- cutreeDynamic(hc, distM = as.matrix(na.dist(mat)), minClusterSize = 1)
## ct <- cutreeHybrid(hc, distM = as.matrix(na.dist(mat)), minClusterSize = 1)
## df <- data.frame(drug = labels(hc), cluster = ct)
## Plot the tree to see manually where we should make the cut
## plot(hc)

mat <- ohsu.drc.scaled[common.drugs,]
ret <- plot.drug.clusters(mat, "12", drug.map, hc.file = "ohsu-common-drugs-hc.pdf")
ohsu.drug.clusters <- ret$clusters
pdf("ohsu-common-drugs-clustered-heatmap.pdf")
print(ret$heatmap)
d <- dev.off()

mat <- fimm.drc.scaled[common.drugs,]
ret <- plot.drug.clusters(mat, "12", drug.map, hc.file = "fimm-common-drugs.hc.pdf")
fimm.drug.clusters <- ret$clusters
pdf("fimm-common-drugs-clustered-heatmap.pdf")
print(ret$heatmap)
d <- dev.off()

jaccard.index <- function(set1, set2) {
  if(length(intersect(set1, set2)) == 0) { return(0) }
  length(intersect(set1, set2))/length(union(set1, set2))
}

## Perform meta-clustering using the Jaccard similarity as a distance metric/edge weight
meta.mcl_ <- function(cls, ...) {
  edges <- ldply(cls,
                .fun = function(cluster1) {
                         df.c1 <- ldply(cls,
                                        .fun = function(cluster2) { 
                                                 sim <- jaccard.index(cluster1, cluster2)
                                                 data.frame(jaccard = sim)
                                        })
                         colnames(df.c1) <- c("node2", "jaccard")   
                         df.c1 
                })
  colnames(edges) <- c("node1", "node2", "jaccard")
  ## Remove self-loops
##  flag <- (edges$method1 == edges$method2) & (edges$cluster1 == edges$cluster2)
##  edges <- edges[!flag,]

  adj.matrix <- as.matrix(get.adjacency(graph.data.frame(edges[, c("node1", "node2", "jaccard")]), attr="jaccard"))
  mcl.out <- mcl(adj.matrix, ...)
  df <- data.frame(meta.cluster = 0, name = colnames(adj.matrix))
  if("Cluster" %in% names(mcl.out)) {
    df$meta.cluster <- mcl.out$Cluster
  }
  df
}

meta.mcl <- function(clusters, ...) {
  cls <- unlist(clusters, recursive = FALSE)
  meta.mcl_(cls, ...)
}

## df <- meta.mcl(clusters, inflation = 1, addLoops = FALSE)

source("calc-mcl.R")

## Change format from item to cluster assigments to a list of all items in a cluster
get.cluster.items <- function(res, cluster.col, item.col) {
  cluster.members <- dlply(res[res[, cluster.col] != 0,,drop=FALSE], .variables = cluster.col, 
                           .fun = function(df) paste(df[, item.col], collapse = ","))
  cluster.members
}

## Expand clusters to items within that cluster
expand.cluster.to.items <- function(meta.cluster, cluster.to.item.map) {
  clusters <- unlist(strsplit(meta.cluster, split=",[ ]*"))
  items <- unlist(llply(clusters, .fun = function(cls) cluster.to.item.map[[cls]]))
  paste(sort(items), collapse=",")
}

## Finally, find drugs that occur twice in a cluster--i.e., are included from both data sets
find.consistently.coclustered.drugs <- function(cluster.of.drugs) {
  clusters.of.consistent.drugs <- llply(clusters.of.drugs, .parallel = FALSE,
                                         .fun = function(cluster) {
                                                  items <- unlist(strsplit(cluster, split=",[ ]*"))
                                                  if(!(any(duplicated(items)))) { return(NA) }
                                                  items <- items[duplicated(items)]
                                                  items
                                         })
  clusters.of.consistent.drugs <- clusters.of.consistent.drugs[!is.na(clusters.of.consistent.drugs)]

  ## Limit to cases in which there are at least 2 consistent drugs in a cluster
  clusters.of.consistent.drugs <- llply(clusters.of.consistent.drugs,
                                         .fun = function(cluster) {
                                                  items <- cluster
                                                  if(length(items) == 1) { return(NA) }
                                                  items
                                         })                                   
  clusters.of.consistent.drugs <- clusters.of.consistent.drugs[!is.na(clusters.of.consistent.drugs)]
  clusters.of.consistent.drugs

}    


inflations <- seq(from=1, to=10, by=0.1)
names(inflations) <- inflations

## clusters <- list("ohsu" = orig.ohsu.drug.clusters, "fimm" = orig.fimm.drug.clusters)
clusters <- list("ohsu" = ohsu.drug.clusters, "fimm" = fimm.drug.clusters)

sils <- ldply(inflations, 
              .fun = function(inflation) {
                       data.frame(sil = calc.mcl.weighted.silhouette(clusters, n.iters = 1000, frac.to.downsample = 0.8, addLoops = FALSE, inflation = inflation, allow1 = TRUE))
                     })
colnames(sils) <- c("inflation", "sil")

save.image(".Rdata.inflate")
stop("stop")

inflations2 <- seq(from=1.4, to=2.0, by=0.1)
inflations2 <- seq(from=1.4, to=1.8, by=0.4)
names(inflations2) <- inflations2

sils2 <- ldply(inflations2, 
              .fun = function(inflation) {
                       data.frame(sil = calc.mcl.weighted.silhouette(clusters, n.iters = 1000, frac.to.downsample = 0.8, addLoops = FALSE, inflation = inflation, allow1 = TRUE))
                     })
colnames(sils2) <- c("inflation", "sil")


cls <- unlist(clusters, recursive = FALSE)
cluster.to.item.map <- cls
mcl.all.out <- meta.mcl_(cls, addLoops = FALSE, allow1 = TRUE, inflation = 1.5)

clusters.of.clusters <- get.cluster.items(mcl.all.out, "meta.cluster", "name")
clusters.of.drugs <- unlist(llply(unlist(clusters.of.clusters), .fun = function(meta.cluster) expand.cluster.to.items(meta.cluster, cluster.to.item.map)))
clusters.of.consistent.drugs <- find.consistently.coclustered.drugs(cluster.of.drugs) 



## Calculate structural distances for intra- and inter-cluster drugs
df <- ldply(clusters.of.consistent.drugs, .parallel = FALSE,
            .fun = function(cluster1) {
                     cls1.res <- ldply(clusters.of.consistent.drugs, .parallel = FALSE,
                                       .fun = function(cluster2) {
                                                cls2.res <- ldply(cluster1, .parallel = FALSE,
                                                                  .fun = function(drug1) {
                                                                           drg1.res <- ldply(cluster2, .parallel = FALSE,
                                                                                             .fun = function(drug2) {
                                                                                                      data.frame(drug2 = drug2, distance =  drug.struct.dist[drug1, drug2])
                                                                                                    })
                                                                           colnames(drg1.res) <- c("drug2", "distance")
                                                                           drg1.res <- cbind(drug1 = drug1, drg1.res)
                                                                           drg1.res
                                                                         })
                                                cls2.res
                                              })
                     colnames(cls1.res) <- c("cluster2", "drug1", "drug2", "distance")
                     cls1.res
                   })
colnames(df)[1] <- "cluster1"
df <- df[df$drug1 != df$drug2,]
df$relation <- unlist(apply(df[, c("cluster1", "cluster2")], 1, function(row) ifelse(row[1] == row[2], "Intra-Cluster", "Inter-Cluster")))

g <- ggplot(data = df, aes(x = relation, y = distance))
g <- g + geom_violin()
## g <- g + geom_beeswarm()
g <- g + geom_boxplot(width = 0.5)
g <- g + ylab("Tanimoto Distance")
pdf("cluster-relation-vs-tanimoto.pdf")
print(g)
d <- dev.off()

save(clusters.of.consistent.drugs, file="clusters.of.consistent.drugs.Rd")

## Expand back out to format in which a drug is assigned to a cluster
consistent.drug.tbl <- ldply(clusters.of.consistent.drugs,
                                       .fun = function(cluster) {
                                                items <- cluster
                                                data.frame(drug = items)
                                       })
colnames(consistent.drug.tbl) <- c("cluster", "drug")

consistent.drug.tbl <- merge(drug.map, consistent.drug.tbl, by.x = "OHSU_DRUG_NAME", by.y = "drug")
consistent.drug.tbl <- consistent.drug.tbl[order(consistent.drug.tbl$cluster),]                                       
             
consistent.drug.tbl[, c("OHSU_DRUG_NAME", "Mechanism.Targets.1", "cluster")]

mat <- fimm.drc.scaled[as.character(consistent.drug.tbl$OHSU_DRUG_NAME),]   
hm <- plot.drug.heatmap(mat, consistent.drug.tbl, show.col.names = FALSE)
pdf("fimm-consistently-co-clustered-drugs-heatmap.pdf")
print(hm)
d <- dev.off()

mat <- ohsu.drc.scaled[as.character(consistent.drug.tbl$OHSU_DRUG_NAME),]   
hm <- plot.drug.heatmap(mat, consistent.drug.tbl, show.col.names = FALSE)
pdf("ohsu-consistently-co-clustered-drugs-heatmap.pdf")
print(hm)
d <- dev.off()



## Look at intra- vs inter-correlations for these drugs (in both data sets)
## Modify code to do drug classes
## Can we model each class in OHSU and apply to FIMM -- plot resp vs predicted, coloring by drug
## Plot correlation vs mean response


## TODO
## - cluster by drug response
## - look at other annotations
## - predict all within cluster, with covariate for each drug
## - look at overlap with drug covariate model
## - look at prediction accuracy with drug covariate model (limit post hoc to drugs in common)

## hm <- Heatmap(mat, clustering_distance_rows = function(x) na.dist(x), clustering_distance_columns = function(x) na.dist(x), na_col = "black")

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
##sensitive <- unlist(apply(tmp.scaled, 2, function(col) safe.t.test(col, alternative="greater")$p.value < cutoff))
##insensitive <- unlist(apply(tmp.scaled, 2, function(col) safe.t.test(col, alternative="less")$p.value < cutoff))
sensitive <- unlist(apply(tmp.scaled, 2, function(col) log(safe.t.test(col, alternative="greater")$p.value)))
insensitive <- unlist(apply(tmp.scaled, 2, function(col) log(safe.t.test(col, alternative="less")$p.value)))
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

save.image(rdata.file)

cat("Exiting\n")
q(status = 0)

