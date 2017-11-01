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

heatmap.hash <- list("ctrp" = ctrp.cell.line.metadata, "ntap" = NULL)

## Columns whose values are not dependent on the (L.4 vs LL.4) fit.  
merge.cols <- list(
  "ctrp" = c("master_cpd_id", "experiment_id", "ccl_name", "ccl_availability", "ccle_primary_site", "ccle_primary_hist", "ccle_hist_subtype_1", "min.conc", "max.conc", "all.concs.nM", "uniq.concs.nM"),
  "ntap" = c("SMI", "SID", "NAME", "sampleIdentifier", "min.conc", "max.conc", "all.concs.nM", "uniq.concs.nM"))
## The column that can be used to uniquely identify the drug (across all data sets) after merging the drug map
drug.name.col <- "ctrp.name"
## End setup


cat(paste0("Cluster (ranked) gene expr for genes shared across data sets: ", paste(data.sets, collapse=","), "\n"))

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

nms <- names(exprs)
names(nms) <- nms
patient.to.label.map <- llply(nms, .parallel = FALSE,
                              .fun = function(nm) {
                                if(!is.null(heatmap.label.cols[[nm]])) {
                                  tmp <- heatmap.hash[[nm]][, c(patient.id.cols[[nm]], heatmap.label.cols[[nm]])]
                                  colnames(tmp) <- c("pt", "label")
                                  mp <- tmp$label
                                  names(mp) <- tmp$pt
                                  return(mp)
                                }
                                tmp <- unique(colnames(exprs[[nm]]))
                                mp <- rep(nm, length(tmp))
                                names(mp) <- tmp
                                mp
                              }) 

## Merge the matrices.
expr.matrix <- do.call("cbind", exprs)
labels <- unlist(llply(names(exprs), .fun = function(nm) patient.to.label.map[[nm]][colnames(exprs[[nm]])]))
short.labels <- unlist(lapply(labels, function(lb) substr(lb, 1, 10)))

## Rank each drug independently for each patient (patient = col).
rank.matrix <- apply(expr.matrix, 2, rank)

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

## This is ntap cluster -- found manually
## labels(hor.dendro[[2]][[2]][[2]][[2]][[1]])

## Just traversed the dendrogram until I found the cluster holding ntap
lbs <- labels(hor.dendro)
lbs <- lbs[grepl(lbs, pattern="ctrp")]
lbs <- gsub(lbs, pattern="ctrp.", replacement="")
translated <- patient.to.label.map[["ctrp"]][lbs]
all.tbl <- as.data.frame(table(translated))
rownames(all.tbl) <- all.tbl$translated
ag <- "autonomic_ganglia"
cns <- "central_nervous_system"
cns.all.num <- all.tbl[cns, "Freq"]
ag.all.num <- all.tbl[ag, "Freq"]
if(is.na(cns.all.num)) { cns.all.num <- 0 }
if(is.na(ag.all.num)) { ag.all.num <- 0 }
n.all.tot <- length(lbs)
cat(paste0("Proportion of ", ag, " and ", cns, " in all data = ", ag.all.num, " + ", cns.all.num, " = ", (ag.all.num + cns.all.num) / n.all.tot))

hd <- hor.dendro[[2]][[2]][[2]][[2]][[2]][[1]]
plot(hd)
lbs <- labels(hd)
lbs <- lbs[grepl(lbs, pattern="ctrp")]
lbs <- gsub(lbs, pattern="ctrp.", replacement="")
translated <- patient.to.label.map[["ctrp"]][lbs]
tbl <- as.data.frame(table(translated))
rownames(tbl) <- tbl$translated
cns.num <- tbl[cns, "Freq"]
ag.num <- tbl[ag, "Freq"]
if(is.na(cns.num)) { cns.num <- 0 }
if(is.na(ag.num)) { ag.num <- 0 }
n.tot <- length(lbs)
cat(paste0("Proportion of ", ag, " and ", cns, " in closest branch = ", ag.num, " + ", cns.num, " = ", (ag.num + cns.num) / n.tot))
mat <- rbind(c(cns.num+ag.num,n.tot-(cns.num+ag.num)), c(cns.all.num+ag.all.num,n.all.tot-(cns.all.num+ag.all.num)))
ft <- fisher.test(mat)
print(ft)

hd <- hor.dendro[[2]][[2]][[2]][[2]]
plot(hd)
lbs <- labels(hd)
lbs <- lbs[grepl(lbs, pattern="ctrp")]
lbs <- gsub(lbs, pattern="ctrp.", replacement="")
translated <- patient.to.label.map[["ctrp"]][lbs]
tbl <- as.data.frame(table(translated))
rownames(tbl) <- tbl$translated
cns.num <- tbl[cns, "Freq"]
ag.num <- tbl[ag, "Freq"]
if(is.na(cns.num)) { cns.num <- 0 }
if(is.na(ag.num)) { ag.num <- 0 }
n.tot <- length(lbs)
cat(paste0("Proportion of ", ag, " and ", cns, " in large neighboring branch = ", ag.num, " + ", cns.num, " = ", (ag.num + cns.num) / n.tot))
mat <- rbind(c(cns.num+ag.num,n.tot-(cns.num+ag.num)), c(cns.all.num+ag.all.num,n.all.tot-(cns.all.num+ag.all.num)))
ft <- fisher.test(mat)
print(ft)

## From above, it looks like their may be a mild tendency for ntap samples to cluster with autonomic ganglia and central nervous system samples.

## Look at the over-representation of each cancer type within the parental subtree of the ntap subtree
hd <- hor.dendro[[2]][[2]][[2]][[2]]
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
  if((ft$estimate > 1) && (ft$p.value < 0.05)) {
    cat(paste0("Over-represented roportion of ", typ, " in parental subtree = ", typ.num, " / ", n.tot, " = ", typ.num / n.tot, " pval : ", ft$p.value, "\n"))
  } else if((ft$estimate < 1) && (ft$p.value < 0.25)) {
    cat(paste0("Under-represented roportion of ", typ, " in parental subtree = ", typ.num, " / ", n.tot, " = ", typ.num / n.tot, " pval : ", ft$p.value, "\n"))
  }
}



