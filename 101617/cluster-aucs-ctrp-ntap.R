suppressPackageStartupMessages(library("foreach"))
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library(tidyr))

suppressPackageStartupMessages(library("minpack.lm"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("gridExtra"))
suppressPackageStartupMessages(library("caTools"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("maftools"))
suppressPackageStartupMessages(library("pcaMethods"))

library(ggbeeswarm)
library(plyr)
library(stringr)
library(synapseClient)

num.cores <- detectCores()
if(!is.na(num.cores) && (num.cores > 1)) {
  suppressPackageStartupMessages(library("doMC"))
  cat(paste("Registering ", num.cores-1, " cores.\n", sep=""))
  registerDoMC(cores=(num.cores-1))
}

synapseLogin()

## Begin setup (this will need to be changed as the data sets in the intersection change)
data.sets <- c("ctrp", "ntap")
drug.mapping.synId <- "syn11244433"
drug.map.drug.id.cols <- c("ctrp.id", "ntap.smiles")
## These fit files are assumed to be for AUCs/DSS values calculated across shared concentration ranges across data sets.
data.set.dss.fit.synIds <- c("syn11257776", "syn11257777")
names(data.set.dss.fit.synIds) <- data.sets
data.set.dss.fit.files <- c("ctrp-cn-dss-t0.tsv", "ntap-cn-dss-t0.tsv")
data.set.id.cols <- c("master_cpd_id", "SMI")
screen.id.cols <- list("ctrp" = c("experiment_id", "master_cpd_id", "master_ccl_id"), "ntap" = c("SMI", "sampleIdentifier"))
patient.id.cols <- list("ctrp" = "ccl_name", "ntap" = "sampleIdentifier")
heatmap.label.cols  <- list("ctrp" = "ccle_primary_hist", "ntap" = NULL)
heatmap.label.cols  <- list("ctrp" = "ccle_primary_site", "ntap" = NULL)
## Columns whose values are not dependent on the (L.4 vs LL.4) fit.  
merge.cols <- list(
  "ctrp" = c("master_cpd_id", "experiment_id", "ccl_name", "ccl_availability", "ccle_primary_site", "ccle_primary_hist", "ccle_hist_subtype_1", "min.conc", "max.conc", "all.concs.nM", "uniq.concs.nM"),
  "ntap" = c("SMI", "SID", "NAME", "sampleIdentifier", "min.conc", "max.conc", "all.concs.nM", "uniq.concs.nM"))
## The column that can be used to uniquely identify the drug (across all data sets) after merging the drug map
drug.name.col <- "ctrp.name"
## End setup

cat(paste0("Cluster (ranked) AUCs for shared concentration ranges of drugs shared across data sets: ", paste(data.sets, collapse=","), "\n"))

## Load in the fit files
fits <- llply(data.set.dss.fit.synIds, .parallel = FALSE,
              .fun = function(synId) {
                       obj <- synGet(synId, downloadFile = TRUE)
                       file <- getFileLocation(obj)
                       read.table(file, sep="\t", header=TRUE, as.is=TRUE, comment.char="", quote="\"")
              })

## Load the map of drugs between data sets
obj <- synGet(drug.mapping.synId, downloadFile = TRUE)
drug.map <- read.table(getFileLocation(obj), sep="\t", header=TRUE, quote="\"", comment.char="")

## Add a dummy id to the drug.map
## drug.map$dummyId <- 1:nrow(drug.map)

## Restrict the fits to those drugs in common across the data sets -- do so by merging.
for(i in 1:length(fits)) {
  ds <- names(fits)[i]
  fits[[ds]] <- merge(fits[[ds]], drug.map, by.x = data.set.id.cols[i], by.y = drug.map.drug.id.cols[i], all = FALSE)
}

## Convert all fits to a drug (AUC) by patient matrix.
source("../common/data-preprocessing.R")

response.matrices <- list()
for(ds in names(fits)) {
  print(ds)
  response.matrices[[ds]] <- prepare.drug.response.matrix(fits[[ds]], drug.col = drug.name.col, patient.col = patient.id.cols[[ds]], response.col = "dss.auc.ll4", drugs = NULL)
}

## Make sure matrices all hold the same drugs in the same order 
common.drugs <- Reduce("intersect", lapply(response.matrices, rownames))
response.matrices <- llply(response.matrices, .fun = function(mat) mat[common.drugs, ])

fs <- names(fits)
names(fs) <- fs
patient.to.label.map <- llply(fs, .parallel = FALSE,
                              .fun = function(nm) {
                                if(!is.null(heatmap.label.cols[[nm]])) {
                                  tmp <- unique(fits[[nm]][, c(patient.id.cols[[nm]], heatmap.label.cols[[nm]])])
                                  colnames(tmp) <- c("pt", "label")
                                  mp <- tmp$label
                                  names(mp) <- tmp$pt
                                  return(mp)
                                }
                                tmp <- unique(fits[[nm]][, patient.id.cols[[nm]], drop = TRUE])
                                mp <- rep(nm, length(tmp))
                                names(mp) <- tmp
                                mp
                              }) 

## Merge the matrices.
response.matrix <- do.call("cbind", response.matrices)
labels <- unlist(llply(names(response.matrices), .fun = function(nm) patient.to.label.map[[nm]][colnames(response.matrices[[nm]])]))
short.labels <- unlist(lapply(labels, function(lb) substr(lb, 1, 10)))

## Rank each drug independently for each patient (patient = col).
rank.matrix <- apply(response.matrix, 2, rank)

## Cluster the ranked AUCs to see how and whether the data sets are interleaved.
##rownames(rank.matrix) <- NULL
##colnames(rank.matrix) <- NULL
do.plot <- FALSE
if(do.plot) {
  library(RColorBrewer)
  ## col.map <- data.frame(labels = unique(labels), color = I(brewer.pal(length(unique(labels)), name = 'Dark2')))
  col.map <- data.frame(labels = unique(labels), color = I(rainbow(length(unique(labels)))))
  rownames(col.map) <- col.map$labels
  cols <- col.map[labels,c("color")]

  heatmap(rank.matrix, scale = "none", ColSideColors = cols)
  col.map$short.labels <- unlist(laply(col.map$labels, function(lb) substr(lb, 1, 10)))
  legend(x = "topright", legend = col.map$short.labels, col = col.map$color, pch = 21, bty = 'n', xjust = 1)

  hr <- hclust(as.dist(1-cor(rank.matrix, method="pearson")), method="ward.D"); 
  hc <- hclust(as.dist(1-cor(t(rank.matrix), method="pearson")), method="ward.D"); 
  plot(hr, labels = short.labels, cex = 0.5)
  hor.dendro <- as.dendrogram(hr)
}

hr <- hclust(as.dist(1-cor(rank.matrix, method="pearson")), method="ward.D"); 
hor.dendro <- as.dendrogram(hr)


lbs <- labels(hor.dendro)
lbs <- lbs[grepl(lbs, pattern="ctrp")]
lbs <- gsub(lbs, pattern="ctrp.", replacement="")
all.tbl <- as.data.frame(table(patient.to.label.map[["ctrp"]][lbs]))
rownames(all.tbl) <- all.tbl$Var1
ag <- "autonomic_ganglia"
cns <- "central_nervous_system"
cns.all.num <- all.tbl[cns, "Freq"]
ag.all.num <- all.tbl[ag, "Freq"]
if(is.na(cns.all.num)) { cns.all.num <- 0 }
if(is.na(ag.all.num)) { ag.all.num <- 0 }
n.all.tot <- length(lbs)
cat(paste0("Proportion of ", ag, " and ", cns, " in all data = ", ag.all.num, " + ", cns.all.num, " = ", (ag.all.num + cns.all.num) / n.all.tot))

## Just traversed the dendrogram until I found the following cluster (hd) holding ntap
## First, confirm that ntap is in this cluster
## NTAP cludster is hor.dendro[[1]][[1]]
## Immediate subtree containing the ntap subtree in [[1]][[1]][[2]]
hd <- hor.dendro[[1]][[1]][[2]]
lbs <- labels(hd)
if(!any(grepl(pattern="ntap", x=lbs))) {
  cat("This subtree was supposed to contain ntap samples, but does not!\n")
  cat(paste0(paste(lbs, collapse=","), "\n"))
}  else {
  cat("As expected, this subtree contains ntap samples\n")
}
lbs <- lbs[grepl(lbs, pattern="ctrp")]
lbs <- gsub(lbs, pattern="ctrp.", replacement="")
tbl <- as.data.frame(table(patient.to.label.map[["ctrp"]][lbs]))
rownames(tbl) <- tbl$Var1
cns.num <- tbl[cns, "Freq"]
ag.num <- tbl[ag, "Freq"]
if(is.na(cns.num)) { cns.num <- 0 }
if(is.na(ag.num)) { ag.num <- 0 }
n.tot <- length(lbs)
cat(paste0("Proportion of ", ag, " and ", cns, " in closest branch = ", ag.num, " + ", cns.num, " = ", (ag.num + cns.num) / n.tot))
mat <- rbind(c(cns.num+ag.num,n.tot-(cns.num+ag.num)), c(cns.all.num+ag.all.num,n.all.tot-(cns.all.num+ag.all.num)))
ft <- fisher.test(mat)
print(ft)

lbs <- labels(hor.dendro[[1]][[1]])
lbs <- lbs[grepl(lbs, pattern="ctrp")]
lbs <- gsub(lbs, pattern="ctrp.", replacement="")
tbl <- as.data.frame(table(patient.to.label.map[["ctrp"]][lbs]))
rownames(tbl) <- tbl$Var
cns.num <- tbl[cns, "Freq"]
ag.num <- tbl[ag, "Freq"]
if(is.na(cns.num)) { cns.num <- 0 }
if(is.na(ag.num)) { ag.num <- 0 }
n.tot <- length(lbs)
cat(paste0("Proportion of ", ag, " and ", cns, " in closest two branches = ", ag.num, " + ", cns.num, " = ", (ag.num + cns.num) / n.tot))
mat <- rbind(c(cns.num+ag.num,n.tot-(cns.num+ag.num)), c(cns.all.num+ag.all.num,n.all.tot-(cns.all.num+ag.all.num)))
ft <- fisher.test(mat)
print(ft)

## From above, it looks like their may be a mild tendency for ntap samples to cluster with autonomic ganglia and central nervous system samples.
## But below shows that this is really driven by autonomic ganglia.

## Look at the over-representation of each cancer type within the parental subtree of the ntap subtree
hd <- hor.dendro[[1]][[1]][[2]]
lbs <- labels(hd)
if(!any(grepl(pattern="ntap", x=lbs))) {
  cat("This subtree was supposed to contain ntap samples, but does not!\n")
  cat(paste0(paste(lbs, collapse=","), "\n"))
}  else {
  cat("As expected, this subtree contains ntap samples\n")
}
lbs <- lbs[grepl(lbs, pattern="ctrp")]
lbs <- gsub(lbs, pattern="ctrp.", replacement="")
tbl <- as.data.frame(table(patient.to.label.map[["ctrp"]][lbs]))
n.tot <- length(lbs)
rownames(tbl) <- tbl$Var1
for(typ in rownames(all.tbl)) {
##  cat(paste0("Testing over-representation of ", typ, "\n"))
  typ.all.num <- all.tbl[typ, "Freq"]
  typ.num <- tbl[typ, "Freq"]
  if(is.na(typ.all.num)) { typ.all.num <- 0 }
  if(is.na(typ.num)) { typ.num <- 0 }
  mat <- rbind(c(typ.num,n.tot-(typ.num)), c(typ.all.num,n.all.tot-(typ.all.num)))
  ft <- fisher.test(mat)
  ## estimate is odds ratio
  if((ft$estimate > 1) && (ft$p.value < 0.25)) {
    cat(paste0("Over-represented proportion of ", typ, " in parental subtree = ", typ.num, " / ", n.tot, " = ", typ.num / n.tot, " pval : ", ft$p.value, "\n"))
  } else if((ft$estimate < 1) && (ft$p.value < 0.25)) {
    cat(paste0("Under-represented proportion of ", typ, " in parental subtree = ", typ.num, " / ", n.tot, " = ", typ.num / n.tot, " pval : ", ft$p.value, "\n"))
  }
}

hd <- hor.dendro[[1]][[1]]
lbs <- labels(hd)
if(!any(grepl(pattern="ntap", x=lbs))) {
  cat("This subtree was supposed to contain ntap samples, but does not!\n")
  cat(paste0(paste(lbs, collapse=","), "\n"))
}  else {
  cat("As expected, this subtree contains ntap samples\n")
}
lbs <- lbs[grepl(lbs, pattern="ctrp")]
lbs <- gsub(lbs, pattern="ctrp.", replacement="")
tbl <- as.data.frame(table(patient.to.label.map[["ctrp"]][lbs]))
n.tot <- length(lbs)
rownames(tbl) <- tbl$Var1
for(typ in rownames(all.tbl)) {
##  cat(paste0("Testing over-representation of ", typ, "\n"))
  typ.all.num <- all.tbl[typ, "Freq"]
  typ.num <- tbl[typ, "Freq"]
  if(is.na(typ.all.num)) { typ.all.num <- 0 }
  if(is.na(typ.num)) { typ.num <- 0 }
  mat <- rbind(c(typ.num,n.tot-(typ.num)), c(typ.all.num,n.all.tot-(typ.all.num)))
  ft <- fisher.test(mat)
  ## estimate is odds ratio
  if((ft$estimate > 1) && (ft$p.value < 0.25)) {
    cat(paste0("Over-represented roportion of ", typ, " in grandparental subtree = ", typ.num, " / ", n.tot, " = ", typ.num / n.tot, " pval : ", ft$p.value, "\n"))
  } else if((ft$estimate < 1) && (ft$p.value < 0.25)) {
    cat(paste0("Under-represented roportion of ", typ, " in grandparental subtree = ", typ.num, " / ", n.tot, " = ", typ.num / n.tot, " pval : ", ft$p.value, "\n"))
  }
}



