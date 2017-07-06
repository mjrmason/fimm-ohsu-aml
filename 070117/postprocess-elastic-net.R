library(synapseClient)
library(data.table)
library(ggplot2)
library(gridExtra)

suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(caret))
suppressPackageStartupMessages(library(glmnet))
library( ReporteRs )
library(scales)

suppressPackageStartupMessages(library("openxlsx"))

synapseLogin()

synId <- "syn8488824"
obj <- synGet(id=synId, downloadFile = TRUE)
file <- getFileLocation(obj)
fimm.raw.dr <- read.xlsx(file)

## path <- "ohsu.dss.t0.tsv"
synId <- "syn10083494"
obj <- synGet(synId, downloadFile = TRUE)
file <- getFileLocation(obj)
ohsu.dss.t0 <- as.data.frame(fread(file))
## ohsu.dss.t0$ic50.ll4 <- remove.outliers(ohsu.dss.t0$e.ll4)
## ohsu.dss.t0$ic50.l4 <- remove.outliers(ohsu.dss.t0$e.l4)

synId <- "syn10083488"
obj <- synGet(synId, downloadFile = TRUE)
file <- getFileLocation(obj)
fimm.dss.t0 <- as.data.frame(fread(file))

load(".Rdata.ohsu")

## Calculate variable importance a la caret
## i.e., scaled absolute magnitude excluding intercept
## caret::getModelInfo("glmnet")$glmnet$varImp
varImpGLMnet <- function(table, coef.name.col = "coef.name", coef.value.col = "coef.value") {
  tbl <- table[!grepl(pattern="Intercept", table[,coef.name.col]),,drop = FALSE]
  tbl[,coef.value.col] <- abs(tbl[,coef.value.col]) * 100 / max(abs(tbl[,coef.value.col]))
  o <- order(tbl[,coef.value.col], decreasing = TRUE)
  tbl <- tbl[o,]
  tbl
}

## Create a drug vs gene heatmap of average importance

translate.ensg.to.symbol <- function(ids) {
  library(biomaRt)
  # Get a mapping from Ensembl id to symbol
  ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
  bm <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol'),
              filters = 'ensembl_gene_id', 
              values = ids, 
              mart = ensembl)
  names(bm) <- c("ID", "TRANSLATION")
  bm <- bm[!(bm$TRANSLATION %in% c("")),]
  bm
}

make.drug.vs.gene.matrix <- function(lst, coef.value.col = "coef.value", coef.name.col = "coef.name") {
  ## Take the ridge results
  alpha.indx <- 2
  num.bootstraps <- length(lst[[1]])
  num.drugs <- length(lst)
  all.coeffs <- ldply(lst, .parallel = FALSE,
                     .fun = function(lst.drug) {
                               ldply(lst.drug, .parallel = FALSE,
                                     .fun = function(lst.boot) {
                                              if(!is.null(lst.boot$coeffs[[alpha.indx]])) {
                                                tbl <- varImpGLMnet(lst.boot$coeffs[[alpha.indx]])
                                                if(nrow(tbl) >= 1) {
                                                  tbl$drug <- lst.boot$drug
                                                }
                                                tbl
                                              }
                                     })
                     })
  ## Calculate the mean variable importance of each gene within each drug
  ## i.e., sum and divide by the number of bootstraps
  all.coeffs <- ddply(all.coeffs, .variables = c(coef.name.col, "drug"),
                      .parallel = FALSE,
                      .fun = function(df) { 
                                vec <- c(as.character(df[1,coef.name.col]), sum(df[,coef.value.col]) / num.bootstraps, df[1,"drug"])
                                names(vec) <- c(coef.name.col, coef.value.col, "drug")
                                vec
                      })
  all.coeffs[,coef.value.col] <- as.numeric(all.coeffs[,coef.value.col])
  all.coeffs
  
  ## Spread into a gene vs drug matrix 
  all.coeffs <- spread_(all.coeffs, key="drug", value=coef.value.col, fill=0)
  
  ## Translate ENSG gene ids to symbols
  bm <- translate.ensg.to.symbol(all.coeffs[, coef.name.col])
  
  ## If any are duplicated, just keep the ENSG name
  flag <- duplicated(bm$TRANSLATION, fromLast = TRUE) | duplicated(bm$TRANSLATION, fromLast = FALSE)
  bm$TRANSLATION[flag] <- bm$ID[flag]
  bm <- unique(bm)
    
  all.coeffs <- merge(bm, all.coeffs, by.x = "ID", by.y = coef.name.col)

  rownames(all.coeffs) <- all.coeffs$TRANSLATION
  all.coeffs <- all.coeffs[,!(colnames(all.coeffs) %in% c("ID", "TRANSLATION", coef.name.col))]
  all.coeffs
}

get.number.genes.in.model <- function(lst, coef.value.col = "coef.value", coef.name.col = "coef.name") {
  ## Take the ridge results
  alpha.indx <- 2
  num.bootstraps <- length(lst[[1]])
  num.drugs <- length(lst)
  all.coeffs <- ldply(lst, .parallel = FALSE,
                      .fun = function(lst.drug) {
                        ldply(lst.drug, .parallel = FALSE,
                              .fun = function(lst.boot) {
                                num.genes <- 0
                                if(!is.null(lst.boot$coeffs[[alpha.indx]])) {
                                  tbl <- varImpGLMnet(lst.boot$coeffs[[alpha.indx]])
                                  if(nrow(tbl) >= 1) {
                                    tbl$drug <- lst.boot$drug
                                    num.genes <- nrow(tbl)
                                  }
                                  vec <- c(num.genes = num.genes, drug = lst.boot$drug)
                                  vec
                                }
                              })
                      })
  all.coeffs
}  

num.genes.fimm.orig <- get.number.genes.in.model(fimm.lst)
num.genes.ohsu.orig <- get.number.genes.in.model(ohsu.lst)

source("analysis")
fimm.pp <- post.process.fits(fimm.lst, "fimm-")
ohsu.pp <- post.process.fits(ohsu.lst, "fimm-")

orig.fimm.coefs <- make.drug.vs.gene.matrix(fimm.lst)
orig.ohsu.coefs <- make.drug.vs.gene.matrix(ohsu.lst)

fimm.coefs <- orig.fimm.coefs
ohsu.coefs <- orig.ohsu.coefs

## Translate FIMM columns
df <- data.frame(drug = colnames(fimm.coefs))
m <- merge(df, ohsu.fimm.drugs, by.x = "drug", by.y = "ID_Drug")
rownames(m) <- m$drug
m <- m[colnames(fimm.coefs),]
colnames(fimm.coefs) <- m$ID_Drug.ohsu
common.cols <- intersect(colnames(fimm.coefs), colnames(ohsu.coefs))
fimm.coefs <- fimm.coefs[,common.cols]
ohsu.coefs <- ohsu.coefs[,common.cols]

## Keep only those genes that have an average variable importance of at 
## least 25 in _some_ drug in either data set.
fimm.genes <- rownames(fimm.coefs)[unlist(lapply(rownames(fimm.coefs), function(name) any(fimm.coefs[name,] > 40)))]
ohsu.genes <- rownames(ohsu.coefs)[unlist(lapply(rownames(ohsu.coefs), function(name) any(ohsu.coefs[name,] > 40)))]
important.genes <- unique(c(fimm.genes, ohsu.genes))
## Set any undefined to zero
if(!all(important.genes %in% rownames(fimm.coefs))) {
  row.names <- important.genes[!(important.genes %in% rownames(fimm.coefs))]
  mat <- matrix(data = 0, nrow=length(row.names), ncol=ncol(fimm.coefs))
  rownames(mat) <- row.names
  colnames(mat) <- colnames(fimm.coefs)
  fimm.coefs <- rbind(fimm.coefs, mat)
}
if(!all(important.genes %in% rownames(ohsu.coefs))) {
  row.names <- important.genes[!(important.genes %in% rownames(ohsu.coefs))]
  mat <- matrix(data = 0, nrow=length(row.names), ncol=ncol(ohsu.coefs))
  rownames(mat) <- row.names
  colnames(mat) <- colnames(ohsu.coefs)
  ohsu.coefs <- rbind(ohsu.coefs, mat)
}
fimm.filtered.coefs <- fimm.coefs[important.genes,]
ohsu.filtered.coefs <- ohsu.coefs[important.genes,]

## Don't reorder rows or columns
hc <- hclust(dist(as.matrix(ohsu.filtered.coefs)))
fimm.filtered.coefs <- fimm.filtered.coefs[rev(hc$labels[hc$order]),]
ohsu.filtered.coefs <- ohsu.filtered.coefs[rev(hc$labels[hc$order]),]

lmat = rbind(c(0,3),c(2,1),c(0,4))
lwid = c(0.1,4)
lhei = c(0.1,4,0.5)
pdf("fimm-heatmap.pdf")
heatmap.2(t(as.matrix(fimm.filtered.coefs)), scale="none", trace="none", density.info="histogram", dendrogram="none", Rowv = FALSE, Colv = FALSE, lmat = lmat, lwid = lwid, lhei = lhei, margins = c(1,20), main = "FIMM")
d <- dev.off()
pdf("ohsu-heatmap.pdf")
heatmap.2(t(as.matrix(ohsu.filtered.coefs)), scale="none", trace="none", density.info="histogram", dendrogram="none", Rowv = FALSE, Colv = FALSE, lmat = lmat, lwid = lwid, lhei = lhei, margins = c(1,20), main = "OHSU")
d <- dev.off()
## heatmap.2(t(as.matrix(ohsu.filtered.coefs)), scale="none", trace="none", density.info="histogram", dendrogram="none", Rowv = FALSE, Colv = FALSE)


##plot(density(fimm.dss.t0$gof.l4[fimm.dss.t0$b.l4 < 0]), col="blue")
##lines(density(fimm.dss.t0$gof.l4[fimm.dss.t0$b.l4 > 0]))

##plot(density(ohsu.dss.t0$gof.l4[ohsu.dss.t0$b.l4 < 0]), col="blue")
##lines(density(ohsu.dss.t0$gof.l4[ohsu.dss.t0$b.l4 > 0]))

