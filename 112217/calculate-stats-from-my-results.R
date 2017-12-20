rdata.file <- ".Rdata.validate.112217.one.ohsu.no.filtering"
load(rdata.file)

## rm(list = ls())

source("../common/models.R")
source("../common/utils.R")
library(plyr)
library(dplyr)
library(gplots)
library(tidyr)
library(synapseClient)

prefix <- "observed"

synapseLogin()

train.indx <- 1
train.set <- train.set.names[[train.indx]]
l <- prepare.drug.response.and.expr.matrices(train.dss.args[[train.indx]], train.expr.args[[train.indx]], drugs = train.common.drugs[[train.indx]], 
                                             drug.col = train.drug.cols[[train.indx]], patient.col = train.patient.cols[[train.indx]], 
                                             response.col = train.response.cols[[train.indx]])
ohsu.drc.orig <- l[["drc.df"]]



## Get the OHSU covariates
synId <- "syn8149174"
obj <- synGet(id=synId, downloadFile = TRUE)
file <- getFileLocation(obj)
tmp <- read.xlsx(file, sheet = 1)
ohsu.metadata <- tmp[, c("patient_id", "gender", "lab_id", "specimen_type")]

tmp <- read.xlsx(file, sheet = 3)
tmp <- tmp[grepl(pattern="FAB", ignore.case=TRUE, tmp$lab_type),]
tmp <- unique(tmp[, c("patient_id", "lab_result")])
colnames(tmp) <- c("patient_id", "fab")
ohsu.metadata <- merge(ohsu.metadata, tmp, by = "patient_id", all = TRUE)

tmp <- read.xlsx(file, sheet = 5)
tmp <- tmp[as.numeric(tmp$age) > 0,]
tmp <- tmp[, c("patient_id", "lab_id", "age")]
tmp <- ddply(unique(tmp[, c("patient_id", "lab_id", "age")]), c("patient_id", "lab_id"), .fun = function(df) { df$age <- mean(df$age); return(unique(df)) })
ohsu.metadata <- merge(ohsu.metadata, tmp[, c("patient_id", "lab_id", "age")], by = c("patient_id", "lab_id"), all = TRUE)

tmp <- read.xlsx(file, sheet = 5)
tmp <- tmp[as.numeric(tmp$age) > 0,]
tmp <- tmp[, c("patient_id", "lab_id", "response")]
ohsu.simplified.response.tbl <- ddply(tmp[!is.na(tmp$response),], c("patient_id", "lab_id"), 
                                      .fun = function(df) {
                                        if(any(grepl(x=df$response, pattern="Complete Response"))) { return("Complete Response")}
                                        if(any(df$response == "Hematologic CR")) { return("Hematologic CR")}
                                        if(any(df$response == "Refractory")) { return("Refractory")}
                                        if(any(df$response == "Supportive/Palliative Care")) { return("Supportive/Palliative Care")}
                                        if(any(df$response == "Unknown")) { return("Unknown")}
                                        warning(paste0("What is ", df$response, "\n"))
                                      })
colnames(ohsu.simplified.response.tbl) <- c("patient_id", "lab_id", "response")
ohsu.metadata <- merge(ohsu.metadata, ohsu.simplified.response.tbl, by = c("patient_id", "lab_id"), all = TRUE)


## Remove NAs
ohsu.metadata <- ohsu.metadata[!is.na(ohsu.metadata$patient_id),]

## Map lab_id to seq_ids
synId <- "syn10083817"
obj <- synGet(synId, downloadFile = TRUE)
file <- getFileLocation(obj)
ohsu.rnaseq.sample.summary <- read.table(file, header=TRUE, sep="\t")
ohsu.metadata <- merge(ohsu.metadata, ohsu.rnaseq.sample.summary[, c("Original_LabID", "SeqID")], by.x = "lab_id", by.y = "Original_LabID", all = FALSE)

## One patient has 2 different FAB subtypes
ohsu.metadata <- ohsu.metadata[!(ohsu.metadata$SeqID == "20-00347"),]
rownames(ohsu.metadata) <- ohsu.metadata$SeqID
ohsu.metadata <- ohsu.metadata[, !(colnames(ohsu.metadata) %in% c("lab_id", "patient_id", "SeqID"))]

library(ComplexHeatmap)
library(circlize)

ohsu.metadata <- ohsu.metadata[, c("gender", "specimen_type", "fab", "response", "age")]

na.dist <- function(x,...) {
 t.dist <- dist(x,...)
 t.dist <- as.matrix(t.dist)
 t.limit <- 1.1*max(t.dist,na.rm=T)
 t.dist[is.na(t.dist)] <- t.limit
 t.dist <- as.dist(t.dist)
 return(t.dist)
}

plot.heatmap <- function(mat, annotations) {

   colnames(mat) <- gsub("X(.*)\\.(.*)", "\\1-\\2", colnames(mat))
   common <- intersect(rownames(annotations), colnames(mat))
   annotations <- annotations[common, ]
   mat <- mat[, common]

   col.list <- list()
   for(col in c("gender", "specimen_type", "fab", "age", "response")) {
     flag <- is.na(annotations[, col]) | (annotations[, col] == "NA")
     annotations[flag, col] <- NA
   ##  annotations[flag,col] <- "NA"
     cols <- NULL
     vals <- annotations[, col]
     vals.no.na <- vals[!is.na(vals)]
   ##  vals.no.na <- vals[vals != "NA"]
     if(col != "age") {
       cols <- rainbow(length(unique(vals.no.na)))
       names(cols) <- unique(vals.no.na)
       col.list[[col]] <- cols
     } else {
       annotations[, col] <- as.numeric(annotations[, col])
       vals <- as.numeric(vals)
       cols <- colorRamp2(c(min(0, floor(min(vals, na.rm=TRUE))), ceiling(max(vals, na.rm=TRUE))), c("white", "red"))
     }
   ##  col.list[[col]] <- cols
   }

   ## ha <- HeatmapAnnotation(df = annotations, col = col.list, na_col = "white", show_annotation_name = TRUE)

   colnames(mat) <- NULL
   ## ha <- HeatmapAnnotation(df = annotations[, !(colnames(annotations) %in% c("age"))], resp = anno_boxplot(mat), age = anno_points(annotations$age), col = col.list, na_col = "black", show_annotation_name = TRUE)
   ha <- HeatmapAnnotation(df = annotations[, !(colnames(annotations) %in% c("age"))], age = anno_points(annotations$age), col = col.list, na_col = "black", show_annotation_name = TRUE)
   ha_boxplot <- HeatmapAnnotation(response = anno_boxplot(mat, axis = TRUE, outline = FALSE), show_annotation_name = TRUE)
   Heatmap(mat, top_annotation = ha, clustering_distance_rows = function(x) na.dist(x), clustering_distance_columns = function(x) na.dist(x), na_col = "black", bottom_annotation = ha_boxplot, bottom_annotation_height = unit(3, "cm"))
}

plot.ordered.based.on.entropy <- function(mat) {
  tmp.scaled <- t(scale(t(mat)))
  ## heatmap.2(tmp.scaled, scale="none", trace="none", Rowv = TRUE, Colv = TRUE, dist = na.dist)
  ranked <- ldply(1:nrow(tmp.scaled), .fun = function(i) rank(tmp.scaled[i,])/max(rank(tmp.scaled[i,])))
  entropy <- unlist(apply(ranked, 2, function(col) -sum((col)*log(col))))
##  entropy <- unlist(apply(ranked, 2, function(col) max(col) - min(col)))
  nms <- names(sort(entropy))
  heatmap.2(tmp.scaled[, nms], scale="none", trace="none", Rowv = TRUE, Colv = FALSE, dist = na.dist)
}

stop("stop")

tmp <- ohsu.drc.orig
## tmp <- tmp[rownames(tmp) %in% rownames(mat),]
flag <- unlist(apply(tmp, 1, function(row) all(is.na(row))))
tmp <- tmp[!flag,]
flag <- unlist(apply(tmp, 2, function(col) all(is.na(col))))
tmp.scaled <- t(scale(t(tmp[,!flag])))
all.mean <- colMeans(tmp.scaled, na.rm=TRUE)
pdf(paste0(prefix, "-ohsu-all.pdf"))
plot.heatmap(tmp.scaled, ohsu.metadata)
d <- dev.off()

tbl.all <- read.table("fimm-ohsu-validate-trained-models-112217-one-ohsu-no-filtering-AUC-ridge-pearson-vs-gene-sets-all.xls", sep="\t", header=TRUE)
fdr.col <- "fimm.fdr"
cor.col <- "fimm.val"
drug.col <- "train.drug"
tbl.all <- tbl.all[order(as.numeric(tbl.all[, fdr.col]), decreasing=FALSE),]

## look at the top quantile
n.quant <- floor(0.1 * nrow(tbl.all))

tbl.quant <- tbl.all[1:n.quant,]
tmp <- ohsu.drc.orig
drugs <- tbl.quant$train.drug
drugs <- c("Sirolimus", "MGCD-265", "Dasatinib", "Trametinib", "Selumetinib", "Tanespimycin", "Sorafenib", "Cabozantinib", "Sunitinib")
drugs <- c("Venetoclax", "Erlotinib", "Selumetinib", "Trametinib", "Sorafenib", "Sirolimus", "PD184352", "Foretinib", "Cabozantinib")
tmp <- tmp[rownames(tmp) %in% drugs,]
flag <- unlist(apply(tmp, 1, function(row) all(is.na(row))))
tmp <- tmp[!flag,]
flag <- unlist(apply(tmp, 2, function(col) all(is.na(col))))
tmp.scaled <- t(scale(t(tmp[,!flag])))
validated.mean <- colMeans(tmp.scaled, na.rm=TRUE)
pdf(paste0(prefix, "-ohsu-validated.pdf"))
plot.heatmap(tmp.scaled, ohsu.metadata)
d <- dev.off()

pdf(paste0(prefix, "-mean-response.pdf"))
plot(validated.mean, all.mean[names(validated.mean)], xlab = "OHSU (FIMM-Validated)", ylab = "OHSU (All)")
d <- dev.off()

stop("stop")

heatmap.2(tmp.scaled, scale="none", trace="none", Rowv = TRUE, Colv = TRUE, dist = na.dist)
ranked <- ldply(1:nrow(tmp.scaled), .fun = function(i) rank(tmp.scaled[i,])/max(rank(tmp.scaled[i,])))
entropy <- unlist(apply(ranked, 2, function(col) -sum((col)*log(col))))
entropy <- unlist(apply(ranked, 2, function(col) max(col) - min(col)))
nms <- names(sort(entropy))[1:100]
heatmap.2(tmp.scaled[, nms], scale="none", trace="none", Rowv = TRUE, Colv = TRUE, dist = na.dist)

heatmap.2(tmp.scaled[, nms], scale="none", trace="none", Rowv = TRUE, Colv = TRUE, dist = na.dist)

entropy		       <- unlist(apply(ranked,		 2, function(col) max(col) - min(col)))
nms <- names(sort(entropy)[200:302])

df <- data.frame(sample = nms)
df$sample <- gsub("X(.*)\\.(.*)", "\\1-\\2", df$sample)
write.table(file="ohsu-variably-sensitive-samples.tsv", df, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)

heatmap.2(tmp.scaled[, nms], scale="none", trace="none", Rowv = TRUE, Colv = TRUE, dist = na.dist)

## Load FIMM metadata (including relapse/refractory/diagnosis)
synId <- "syn8270594"
obj <- synGet(id=synId, downloadFile = TRUE, downloadLocation = ".")
file <- getFileLocation(obj)
fimm.metadata <- read.table(file, header=TRUE, sep=",")
disease.stage <- fimm.metadata$diseases.stage
names(disease.stage) <- rownames(fimm.metadata)

library(scales)
dat <- data.frame(sample = names(disease.stage), disease.stage = as.character(disease.stage))
n <- nlevels(disease.stage)
dat.col <- data.frame(stage = unique(as.character(disease.stage)),
            col =brewer_pal()(n))  ## you can also use rainbow(n)
dat <- merge(dat, dat.col, by.x = "disease.stage", by.y = "stage")
rownames(dat) <- dat$sample

do.stats <- function(tbl, cor.col, str) {
  mean.cor <- mean(tbl[, cor.col], na.rm=TRUE)
  sd.cor <- sd(tbl[, cor.col], na.rm=TRUE)
  ci.cor <- quantile(tbl[, cor.col], na.rm=TRUE, probs = c(0.025, 1 - 0.025))
  min.cor <- min(tbl[, cor.col], na.rm=TRUE)
  max.cor <- max(tbl[, cor.col], na.rm=TRUE)
  cat(paste0(paste(str, "n:", nrow(tbl), "mean corr:", mean.cor, "sd:", sd.cor, "min:", min.cor, "max:", max.cor, "95% CI:", ci.cor[1], "-", ci.cor[2], sep=" "), "\n"))
}

get.tbl.of.responses <- function(drugs) {
  tbl <- ldply(drugs, 
               .fun = function(drug) {
                        foo <- extract.predicted.actual(res.auc[["all.comparisons"]][["ohsu"]][["fimm"]][["gene"]], train.drug = drug, alpha = 1, model = "glmnet", s = "lambda.min")
                        df <- data.frame(drug = drug, sample = foo$sample.names, actual = foo$actual)
                        df
                      })
  tbl
}

fdr.cutoff <- 0.2

tbl.all <- read.table("fimm-ohsu-validate-trained-models-112217-one-ohsu-no-filtering-AUC-ridge-pearson-vs-gene-sets-all.xls", sep="\t", header=TRUE)
fdr.col <- "fimm.fdr"
cor.col <- "fimm.val"
drug.col <- "train.drug"
tbl.all <- tbl.all[order(as.numeric(tbl.all[, fdr.col]), decreasing=FALSE),]

## look at the top quantile
n.quant <- floor(0.1 * nrow(tbl.all))

tbl.quant <- tbl.all[1:n.quant,]
do.stats(tbl.quant, cor.col, "Ridge top quant:")
sig.flag <- !is.na(tbl.all[, fdr.col]) & (tbl.all[, fdr.col] < fdr.cutoff)
pos.flag <- !is.na(tbl.all[, cor.col]) & (tbl.all[, cor.col] > 0)
do.stats(tbl.all[sig.flag & pos.flag, ], cor.col, "Ridge sig and pos:")

## Cluster the actual drug responses for all drugs and the top quantile
## Use WGCNA to cut tree
## library(dynamicTreeCut)

tbl <- get.tbl.of.responses(as.character(tbl.quant$train.drug))
mat <- spread(tbl, sample, actual)
rownames(mat) <- mat[,1]
mat <- mat[,-1]
## Cluster the drugs (manually)
hc <- hclust(dist(mat))
k <- 4
cat(paste0("\nRidge quant cut k = ", k, "\n"))
## print(cutree(hc, k = k))
clusters <- invert.list(as.list(cutree(hc, k = k)))
l_ply(clusters, .fun = function(cluster) { cat(paste0("Cluster: ", paste(sort(cluster), collapse=", "), " ; mechanisms = ", paste(sort(unique(tbl.all$Mechanism.Targets[tbl.all$train.drug %in% cluster])), collapse="; "), "\n")) })
cat("\n")

## Cluster the patients
## hc <- hclust(dist(t(mat)))
## cutree(hc, k = 4)

pdf("ridge-top-quant-heatmap.pdf")
## heatmap(as.matrix(mat))
heatmap.2(as.matrix(mat), scale = "none", trace = "none", margins = c(6, 8), ColSideColors = as.character(dat[colnames(mat), "col"]))
d <- dev.off()

tbl <- get.tbl.of.responses(as.character(tbl.all$train.drug))
mat <- spread(tbl, sample, actual)
rownames(mat) <- mat[,1]
mat <- mat[,-1]
pdf("ridge-all-heatmap.pdf")
## heatmap(as.matrix(mat))
heatmap.2(as.matrix(mat), scale = "none", trace = "none", margins = c(6, 8))
d <- dev.off()

tbl.all.ridge <- tbl.all
tbl.quant.ridge <-  tbl.quant
min.cor.ridge <- min(tbl.quant.ridge[, cor.col])

## Repeat for Suleiman's results

tbl.all <- read.table("TR-Ridge.tsv", sep="\t", header=TRUE)
fdr.col <- "TR.fdr"
cor.col <- "TR.Pearson.Correlation"
drug.col <- "train.drug"
tbl.all[,drug.col] <- rownames(tbl.all)
tbl.all <- merge(tbl.all, tbl.all.ridge[, c("train.drug", "Mechanism.Targets", "Class.explained")], by = "train.drug")
tbl.all <- tbl.all[order(as.numeric(tbl.all[, fdr.col]), decreasing=FALSE),]

## look at the top quantile
n.quant <- floor(0.1 * nrow(tbl.all))

tbl.quant <- tbl.all[1:n.quant,]
tbl.min.cor <- tbl.all[tbl.all[, cor.col] >= min.cor.ridge, ]

sig.flag <- !is.na(tbl.all[, fdr.col]) & (tbl.all[, fdr.col] < fdr.cutoff)
pos.flag <- !is.na(tbl.all[, cor.col]) & (tbl.all[, cor.col] > 0)
do.stats(tbl.all[sig.flag & pos.flag, ], cor.col, "TR sig and pos:")
do.stats(tbl.quant, cor.col, "TR top quant:")
do.stats(tbl.min.cor, cor.col, "TR min cor:")


## Drugs in common with ridge
both.drugs <- intersect(tbl.quant$train.drug, tbl.quant.ridge$train.drug)
tr.only.drugs <- setdiff(tbl.quant$train.drug, tbl.quant.ridge$train.drug)
ridge.only.drugs <- setdiff(tbl.quant.ridge$train.drug, tbl.quant$train.drug)
neither.drugs <- setdiff(tbl.all$train.drug, unique(c(both.drugs, tr.only.drugs, ridge.only.drugs)))
mat <- cbind(c(length(both.drugs), length(tr.only.drugs)), c(length(ridge.only.drugs), length(neither.drugs)))
print(mat)
ft <- fisher.test(mat, alternative = "greater")

cat(paste0("Drugs in common between ridge and TR: ", paste(intersect(tbl.quant$train.drug, tbl.quant.ridge$train.drug), collapse=", "), "\n"))
cat(paste0("Drugs unique to TR: ", paste(setdiff(tbl.quant$train.drug, tbl.quant.ridge$train.drug), collapse=", "), "\n"))
cat(paste0("Drugs unique to ridge: ", paste(setdiff(tbl.quant.ridge$train.drug, tbl.quant$train.drug), collapse=", "), "\n"))

cat(paste0("Drugs in common between ridge and TR (min cor): ", paste(intersect(tbl.min.cor$train.drug, tbl.quant.ridge$train.drug), collapse=", "), "\n"))
cat(paste0("Drugs unique to TR (min cor): ", paste(setdiff(tbl.min.cor$train.drug, tbl.quant.ridge$train.drug), collapse=", "), "\n"))
cat(paste0("Drugs unique to ridge (min cor): ", paste(setdiff(tbl.quant.ridge$train.drug, tbl.min.cor$train.drug), collapse=", "), "\n"))

## Plot volcano for TR


## Cluster the actual drug responses for all drugs and the top quantile
## Use WGCNA to cut tree
## library(dynamicTreeCut)

tbl <- get.tbl.of.responses(as.character(tbl.quant$train.drug))
mat <- spread(tbl, sample, actual)
rownames(mat) <- mat[,1]
mat <- mat[,-1]

pdf("tr-top-quant-heatmap.pdf")
## heatmap(as.matrix(mat))
heatmap.2(as.matrix(mat), scale = "none", trace = "none", margins = c(6, 8))
d <- dev.off()

## Cluster the drugs (manually)
hc <- hclust(dist(mat))
k <- 4
cat(paste0("\nTR quant cut k = ", k, "\n"))
## print(cutree(hc, k = k))
clusters <- invert.list(as.list(cutree(hc, k = k)))
l_ply(clusters, .fun = function(cluster) { cat(paste0("Cluster: ", paste(sort(cluster), collapse=", "), " ; mechanisms = ", paste(sort(unique(tbl.all$Mechanism.Targets[tbl.all$train.drug %in% cluster])), collapse="; "), "\n")) })
cat("\n")

tbl <- get.tbl.of.responses(as.character(tbl.min.cor$train.drug))
mat <- spread(tbl, sample, actual)
rownames(mat) <- mat[,1]
mat <- mat[,-1]

pdf("tr-min-cor-heatmap.pdf")
## heatmap(as.matrix(mat))
heatmap.2(as.matrix(mat), scale = "none", trace = "none", margins = c(6, 8))
d <- dev.off()

## Cluster the drugs (manually)
hc <- hclust(dist(mat))
k <- 7
cat(paste0("\nTR quant cut k = ", k, "\n"))
## print(cutree(hc, k = k))
clusters <- invert.list(as.list(cutree(hc, k = k)))
l_ply(clusters, .fun = function(cluster) { cat(paste0("Cluster: ", paste(sort(cluster), collapse=", "), " ; mechanisms = ", paste(sort(unique(tbl.all$Mechanism.Targets[tbl.all$train.drug %in% cluster])), collapse="; "), "\n")) })
cat("\n")

## Cluster the patients
## hc <- hclust(dist(t(mat)))
## cutree(hc, k = 4)

tbl <- get.tbl.of.responses(as.character(tbl.all$train.drug))
mat <- spread(tbl, sample, actual)
rownames(mat) <- mat[,1]
mat <- mat[,-1]
pdf("tr-all-heatmap.pdf")
## heatmap(as.matrix(mat))
heatmap.2(as.matrix(mat), scale = "none", trace = "none", margins = c(6, 8))
d <- dev.off()

