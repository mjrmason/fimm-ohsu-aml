
library(glmnet)
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
## suppressPackageStartupMessages(library("GGally"))
suppressPackageStartupMessages(library("VennDiagram"))
suppressPackageStartupMessages(library("caret"))

library(ggbeeswarm)
library(plyr)
library(stringr)
library(synapseClient)

source("../../common/data-preprocessing.R")
source("../../common/models.R")
source("../../common/plotting.R")

num.cores <- detectCores()
if(!is.na(num.cores) && (num.cores > 1)) {
  suppressPackageStartupMessages(library("doMC"))
  cat(paste("Registering ", num.cores-1, " cores.\n", sep=""))
  registerDoMC(cores=(num.cores-1))
}

synapseLogin()

lm_eqn <- function(df){
    m <- lm(y ~ x, df);
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
         list(a = format(coef(m)[1], digits = 2), 
              b = format(coef(m)[2], digits = 2), 
             r2 = format(summary(m)$r.squared, digits = 3)))
    eq <- substitute(italic(r)^2~"="~r2, 
             list(r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(eq));                 
}

lm_eqn_bq <- function(df){
    m <- lm(y ~ x, df);
    r2 = format(summary(m)$r.squared, digits = 3)
    eq <- bquote(italic(r)^2~"="~.(r2)) 
    eq
}

suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(Biobase))

ensg.to.sym.mapping <- function(gene.ids) {
  # Get a mapping from ensembl id to hugo symbol
  ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")
  bm <- getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id'), 
              filters = 'ensembl_gene_id', 
              values = gene.ids,
              mart = ensembl)
  names(bm) <- c("symbol", "ensg")
  bm <- bm[!(bm$symbol %in% c("")),]
  bm <- bm[!(bm$ensg %in% c("")),]
  bm
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

symbols.to.ensg.mapping <- function(symbols, ensembl, curl = NULL) {
  # Get a mapping from gene symbol to ensembl id
  bm <- getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id'), 
              filters = 'hgnc_symbol', 
              values = symbols, 
              mart = ensembl,
              curl = curl)
  names(bm) <- c("symbol", "ensg")
  bm <- bm[!(bm$ensg %in% c("")),]
  bm
}

## Begin setup (this will need to be changed as the data sets in the intersection change)
data.sets <- c("ohsu", "ctrp", "sanger", "fimm")
data.sets.to.analyze <- data.sets
drug.mapping.synId <- "syn11290968" ## ctrp-sanger-fimm-ohsu-drug-map.tsv
drug.map.drug.id.cols <- c("OHSU_DRUG_NAME", "master_cpd_id", "Drug.ID", "FIMM_Batch_ID_Drug")

## These fit files are assumed to be for AUCs/DSS values calculated across shared concentration ranges across data sets.
data.set.dss.fit.synIds <- c("syn11309113", "syn11309114", "syn11309115", "syn11309116")
names(data.set.dss.fit.synIds) <- data.sets
data.set.dss.fit.files <- c("ohsu-foc-aml-s-aml-outlier1-dss-t0.tsv", "ctrp-foc-aml-s-aml-outlier1-dss-t0.tsv", "sanger-foc-aml-s-aml-outlier1-dss-t0.tsv", "fimm-foc-aml-s-aml-outlier1-dss-t0.tsv")

## Batch-corrected log2 cpm (or rma) expression files (NB: these exclude 2 outliers in OHSU)
data.set.expr.synIds <- c("syn11307673", "syn11307674", "syn11307678", "syn11307680")
names(data.set.expr.synIds) <- data.sets
data.set.expr.files <- c("ohsu-expr-foc-aml-s-aml-outlier1-combat.tsv", "ctrp-expr-foc-aml-s-aml-outlier1-combat.tsv", "sanger-expr-foc-aml-s-aml-outlier1-combat.tsv", "fimm-expr-foc-aml-s-aml-outlier1-combat.tsv")

data.set.drug.id.cols <- list("ohsu" = "inhibitor", "ctrp" = "master_cpd_id", "sanger" = "DRUG_ID", "fimm" = "DRUG_ID")
screen.id.cols <- list("ohsu" = c("patient_id", "inhibitor", "lab_id", "SeqID", "replicant", "time_of_read"),
                       "ctrp" = c("experiment_id", "master_cpd_id", "master_ccl_id", "ccl_name"), 
                       "sanger" = c("IC_RESULTS_ID", "COSMIC_ID", "CELL_LINE_NAME", "DRUG_ID"),
                       "fimm" = c("SCREEN_ID", "PATIENT_ID", "DRUG_ID"))

## Columns whose values are not dependent on the (L.4 vs LL.4) fit.  
merge.cols <- list(
  "ohsu" = c(screen.id.cols[["ohsu"]], "min.conc", "max.conc", "all.concs.nM", "uniq.concs.nM"),
  "ctrp" = c(screen.id.cols[["ctrp"]], "min.conc", "max.conc", "all.concs.nM", "uniq.concs.nM"),
  "sanger" = c(screen.id.cols[["sanger"]], "min.conc", "max.conc", "all.concs.nM", "uniq.concs.nM"),
  "fimm" = c(screen.id.cols[["fimm"]], "min.conc", "max.conc", "all.concs.nM", "uniq.concs.nM"))

## The synapseIds of folders in which to store the fits for each respective data set.
## Here, we will store all in FIMM_BEAT AML Collaboration/Files/Analyses/ctrp-sanger-fimm-ohsu (syn11288886)
fit.folder.synIds <- list("ohsu" = "syn11288886", "ctrp" = "syn11288886", "sanger" = "syn11288886", "fimm" = "syn11288886")
## A string that will be included in the fit file name.
## Here "foc-aml-s-aml" = "fimm-ohsu-ctrp-aml-s-aml"
fit.file.prefix <- "foc-aml-s-aml-outlier"

## Pull in the CTRP cell line metadata so we can see the cancer/histology type of each sample
synId <- "syn5632192"
obj <- synGet(synId, downloadFile = TRUE)
ctrp.cell.line.metadata <- read.table(getFileLocation(obj), sep="\t", header=TRUE, as.is=TRUE, comment.char="", quote="")

## Subset the cell lines to those of interest
ccle.subset <- ctrp.cell.line.metadata[!is.na(ctrp.cell.line.metadata$ccle_hist_subtype_1),]
ccle.subset <- subset(ccle.subset, ccle_hist_subtype_1 == "acute_myeloid_leukaemia")

## Read in the Sanger cell line metadata so we can see the cancer type of each sample (Cell_Lines_Details.xlsx)
synId <- "syn11275809"
sanger.cell.line.metadata <- openxlsx:::read.xlsx(getFileLocation(synGet(synId)), sheet=1)
sanger.subset <- sanger.cell.line.metadata[!is.na(sanger.cell.line.metadata$`Cancer.Type.(matching.TCGA.label)`),]
sanger.subset <- sanger.subset[sanger.subset$`Cancer.Type.(matching.TCGA.label)` == "LAML", ]

patient.id.cols <- list("ohsu" = "SeqID", "ctrp" = "ccl_name", "sanger" = "COSMIC_ID", "fimm" = "SCREEN_ID")
patient.subsets <- list("ohsu" = NULL, "ctrp" = ccle.subset$ccl_name, "sanger" = sanger.subset$COSMIC.identifier, "fimm" = NULL)

drug.name.col <- "cpd_name"

## A string that will be included in the batch-corrected expression output file name.
## Here "foc-aml-s-aml" = "fimm-ohsu-ctrp-aml-s-aml"
file.prefix <- fit.file.prefix

## End setup


cat(paste0("Model AUC ~ expr for data sets: ", paste(data.sets, collapse=","), "\n"))

## Load in the fit files
fits <- llply(data.set.dss.fit.synIds, .parallel = FALSE,
              .fun = function(synId) {
                       obj <- synGet(synId, downloadFile = TRUE)
                       file <- getFileLocation(obj)
                       read.table(file, sep="\t", header=TRUE, as.is=TRUE, comment.char="", quote="\"")
              })

for(ds in names(fits)) {
  switch(ds,
    "ohsu" = {
      ## Add the SeqID to the OHSU drug data.
      ## Read in the rna-seq summary table, which provides a map from the lab_id's used
      ## in the drug response data to the SeqID's used in the expression data.
      synId <- "syn10083817"
      obj <- synGet(synId, downloadFile = TRUE)
      file <- getFileLocation(obj)
      ohsu.rnaseq.sample.summary <- read.table(file, header=TRUE, sep="\t")
      fits[[ds]] <- merge(fits[[ds]], ohsu.rnaseq.sample.summary, by.x = "lab_id", by.y = "Original_LabID")
    },
    "fimm" = {
    },
    "ctrp" = {
      fits[[ds]] <- merge(fits[[ds]], ctrp.cell.line.metadata[, c("master_ccl_id", "ccl_name")])
    },
    "ntap" = {
    },
    "sanger" = {
    },
    {
      stop("Unrecognized data set\n")
    })
}

## Load the map of drugs between data sets
obj <- synGet(drug.mapping.synId, downloadFile = TRUE)
drug.map <- read.table(getFileLocation(obj), sep="\t", header=TRUE, quote="\"", comment.char="")

my.dup <- function(x) duplicated(x, fromLast = TRUE) | duplicated(x, fromLast = FALSE)
for(col in drug.map.drug.id.cols) {
  flag <- !my.dup(drug.map[, col])
  drug.map <- drug.map[flag, ]
}

## Add a dummy id to the drug.map
## drug.map$dummyId <- 1:nrow(drug.map)

## Restrict the fits to those drugs in common across the data sets -- do so by merging.
for(i in 1:length(fits)) {
  ds <- names(fits)[i]
  fits[[ds]] <- merge(fits[[ds]], drug.map, by.x = data.set.drug.id.cols[[ds]], by.y = drug.map.drug.id.cols[i], all = FALSE)
}

common.drugs <- Reduce("intersect", lapply(fits, function(fit) unique(fit[, drug.name.col])))
for(i in 1:length(fits)) {
  ds <- names(fits)[i]
  flag <- fits[[ds]][, drug.name.col] %in% common.drugs
  fits[[ds]] <- fits[[ds]][flag, ]
}

drug.name.tbl <- drug.map[drug.map[, drug.name.col] %in% common.drugs,]

common.drugs.by.data.set <- list()
for(i in 1:length(fits)) {
  ds <- names(fits)[i]
##  common.drugs.by.data.set[[ds]] <- drug.name.tbl[, drug.map.drug.id.cols[i]]
  common.drugs.by.data.set[[ds]] <- drug.name.tbl[, drug.name.col]
}

if(("sanger" %in% names(fits)) && ("ctrp" %in% names(fits))) {
  synId <- "syn11289147"
  obj <- synGet(synId, downloadFile = TRUE)
  gdsc.to.ctrp.cell.line.map <- openxlsx:::read.xlsx(getFileLocation(obj), sheet = 1, startRow = 14)
  gdsc.to.ctrp.cell.line.map <- gdsc.to.ctrp.cell.line.map[, c("GDSC1000.cosmic.id", "CTRP.master.ccl.id")]
  ft <- fits[["sanger"]]
  ft <- merge(ft, gdsc.to.ctrp.cell.line.map, by.x = "COSMIC_ID", by.y = "GDSC1000.cosmic.id", all = FALSE)
  ft <- merge(ft, fits[["ctrp"]], by.x = c("CTRP.master.ccl.id", "cpd_name"), by.y = c("master_ccl_id", "cpd_name"), suffixes = c(".sanger", ".ctrp"))
  g1 <- plot.r2(ft$dss.auc.ll4.ctrp, ft$dss.auc.ll4.sanger)
  g1 <- plot.r2(ft$dss.auc.l4.ctrp, ft$dss.auc.l4.sanger)
  g1 <- g1 + xlab("CTRP AUC")
  g1 <- g1 + ylab("Sanger AUC")
  g2 <- plot.r2(ft$e.ll4.ctrp, ft$e.ll4.sanger)
  g2 <- g2 + xlab("CTRP IC50")
  g2 <- g2 + ylab("Sanger IC50")
  grid.arrange(g1, g2)

  gdsc.to.ctrp.cell.line.map <- read.table(getFileLocation(synGet("syn9988099", downloadFile = TRUE)), sep="\t", header=TRUE, quote="\"", comment.char="")
  ft <- fits[["sanger"]]
  ft <- merge(ft, gdsc.to.ctrp.cell.line.map, by.x = "COSMIC_ID", by.y = "COSMIC_ID", all = FALSE)
  ft$CCLE.name <- gsub(ft$CCLE.name, pattern="^([^_]+).*", replacement="\\1")
  ft <- ft[!is.na(ft$CCLE.name),]
  ft <- merge(ft, fits[["ctrp"]], by.x = c("CCLE.name", "cpd_name"), by.y = c("ccl_name.l4", "cpd_name"), suffixes = c(".sanger", ".ctrp"))

  flag <- (ft$dss.auc.ll4.ctrp < 2*10^5) & (ft$dss.auc.ll4.sanger < 2*10^5)
  g1 <- plot.r2(ft$dss.auc.ll4.ctrp[flag], ft$dss.auc.ll4.sanger[flag])
  g1 <- g1 + xlab("CTRP AUC")
  g1 <- g1 + ylab("Sanger AUC")
  flag <- (ft$e.ll4.ctrp < 50000) & (ft$e.ll4.sanger < 50000)
  g2 <- plot.r2(ft$e.ll4.ctrp[flag], ft$e.ll4.sanger[flag])
  g2 <- g2 + xlab("CTRP IC50")
  g2 <- g2 + ylab("Sanger IC50")
  grid.arrange(g1, g2)

}


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

for(nm in names(exprs)) {
  switch(nm,
    "ohsu" = {
      ## Convert the column names from X20.00347 to 20-00347
      colnames(exprs[[nm]]) <- gsub("X(.*)\\.(.*)", "\\1-\\2", colnames(exprs[[nm]]))
    },
    "fimm" = {
    },
    "ctrp" = {
    },
    "ntap" = {
    },
    "sanger" = {
      colnames(exprs[[nm]]) <- gsub("X(.*)", "\\1", colnames(exprs[[nm]]))
    },
    {
      stop("Unrecognized data set\n")
    })
}

save.image(".Rdata.focs")

## BEGIN modeling setup


## Filter the gene expression data to restrict to highly-expressed genes
expr.filt.matrices <- list()
for(ds in names(exprs)) {
  iqrs <- unlist(apply(exprs[[ds]], 1, IQR))
  expr.filt.matrices[[ds]] <- exprs[[ds]][iqrs > median(iqrs[iqrs > 0]),]
}

study.genes <- (lapply(exprs, rownames))
study.filt.genes <- (lapply(expr.filt.matrices, rownames))
common.genes <- Reduce(intersect, study.filt.genes)

## END modeling setup

## BEGIN modeling

## Do analysis at pathway-, rather than gene-, level
suppressPackageStartupMessages(library("GSVA"))
# use GSA to read in a gmt file
suppressPackageStartupMessages(library("GSA"))

## Download biocarta, kegg, reactome, and hallmark data sets that have been translated to ensembl ids

## Hallmark gene sets: The list in this R object is called hallmark.gene.sets
synId <- "syn10644516"
load(getFileLocation(synGet(synId)))

## KEGG gene sets: The list in this R object is called kegg.gene.sets
synId <- "syn10644519"
load(getFileLocation(synGet(synId)))

## Biocarta gene sets: The list in this R object is called biocarta.gene.sets
synId <- "syn10644517"
load(getFileLocation(synGet(synId)))

## Reactome gene sets: The list in this R object is called reactome.gene.sets
synId <- "syn10644520"
load(getFileLocation(synGet(synId)))

## Do some bootstrap analysis of the ctrp data
ctrp.drug.col <- "master_cpd_id"
ctrp.patient.col <- "ccl_name"

## Assemble training and test sets

## Use ctrp for training and ntap for test

gene.sets <- list("gene" = "gene", "hallmark" = hallmark.gene.sets, "biocarta" = biocarta.gene.sets, "kegg" = kegg.gene.sets, "reactome" = reactome.gene.sets)

train.set.names <- list()
train.dss.args <- list()
train.expr.args <- list()
train.genomic.args <- list()
train.clinical.args <- list()
train.common.drugs <- list()
train.drugs <- list()
train.drug.cols <- list()
train.patient.cols <- list()
train.response.cols <- list()

test.set.names <- list()
test.dss.args <- list()
test.expr.args <- list()
test.genomic.args <- list()
test.clinical.args <- list()
test.common.drugs <- list()
test.drug.cols <- list()
test.patient.cols <- list()
test.response.cols <- list()

ohsu.drug.col <- "inhibitor"
ohsu.patient.col <- "SeqID"

ds <- "ohsu"
if(ds %in% data.sets.to.analyze) {
  set.seed(1234)

  train.index <- createDataPartition(fits[[ds]][, data.set.drug.id.cols[[ds]]], p = .6, list = FALSE)
  ohsu.dss.orig.train <- fits[[ds]][train.index, ]
  ohsu.dss.orig.test <- fits[[ds]][-train.index, ]

  nxt.indx <- length(train.set.names) + 1
  train.set.names[[nxt.indx]] <- "ohsu.train"
  train.dss.args[[nxt.indx]] <- ohsu.dss.orig.train
  train.expr.args[[nxt.indx]] <- expr.filt.matrices[[ds]][common.genes, ] 
  train.genomic.args[nxt.indx] <- list(NULL)
  train.clinical.args[nxt.indx] <- list(NULL)
  train.common.drugs[[nxt.indx]] <- common.drugs.by.data.set[[ds]]
  train.drugs[[nxt.indx]] <- common.drugs.by.data.set[[ds]]
  train.drug.cols[[nxt.indx]] <- drug.name.col
  train.patient.cols[[nxt.indx]] <- patient.id.cols[[ds]]
  train.response.cols[[nxt.indx]] <- "dss.auc.ll4"

  nxt.indx <- length(test.set.names) + 1
  test.set.names[[nxt.indx]] <- "ohsu.test"
  test.dss.args[[nxt.indx]] <- ohsu.dss.orig.test
  test.expr.args[[nxt.indx]] <- expr.filt.matrices[[ds]][common.genes, ]
  test.genomic.args[nxt.indx] <- list(NULL)
  test.clinical.args[nxt.indx] <- list(NULL)
  test.common.drugs[[nxt.indx]] <- common.drugs.by.data.set[[ds]]
  test.drug.cols[[nxt.indx]] <- drug.name.col
  test.patient.cols[[nxt.indx]] <- patient.id.cols[[ds]]
  test.response.cols[[nxt.indx]] <- "dss.auc.ll4"

  rm(ohsu.dss.orig.train)
  rm(ohsu.dss.orig.test)
} ## if("ohsu" %in% data.sets.to.analyze) 

fimm.drug.col <- "DRUG_ID"
fimm.patient.col <- "SCREEN_ID"

ds <- "fimm"
if(ds %in% data.sets.to.analyze) {

  nxt.indx <- length(test.set.names) + 1
  test.set.names[[nxt.indx]] <- ds
  test.dss.args[[nxt.indx]] <- fits[[ds]]
  test.expr.args[[nxt.indx]] <- expr.filt.matrices[[ds]][common.genes, ]
  test.genomic.args[nxt.indx] <- list(NULL)
  test.clinical.args[nxt.indx] <- list(NULL)
  test.common.drugs[[nxt.indx]] <- common.drugs.by.data.set[[ds]]
  test.drug.cols[[nxt.indx]] <- drug.name.col
  test.patient.cols[[nxt.indx]] <- patient.id.cols[[ds]]
  test.response.cols[[nxt.indx]] <- "dss.auc.ll4"
} ## if("fimm" %in% data.sets.to.analyze) 

ctrp.drug.col <- "master_cpd_id"
ctrp.patient.col <- "ccl_name"

if("ctrp" %in% data.sets.to.analyze) {

  nxt.indx <- length(train.set.names) + 1
  train.set.names[[nxt.indx]] <- "ctrp"
  train.dss.args[[nxt.indx]] <- fits[["ctrp"]]
  train.expr.args[[nxt.indx]] <- expr.filt.matrices[["ctrp"]][common.genes, ]
  train.genomic.args[nxt.indx] <- list(NULL)
  train.clinical.args[nxt.indx] <- list(NULL)
  train.common.drugs[[nxt.indx]] <- common.drugs.by.data.set[["ctrp"]]
  train.drugs[[nxt.indx]] <- common.drugs.by.data.set[["ctrp"]]
  train.drug.cols[[nxt.indx]] <- drug.name.col
  train.patient.cols[[nxt.indx]] <- patient.id.cols[["ctrp"]]
  train.response.cols[[nxt.indx]] <- "dss.auc.ll4"

} ## if("ctrp" %in% data.sets.to.analyze) 

if("ntap" %in% data.sets.to.analyze) {

  nxt.indx <- length(test.set.names) + 1
  test.set.names[[nxt.indx]] <- "ntap"
  test.dss.args[[nxt.indx]] <- fits[["ntap"]]
  test.expr.args[[nxt.indx]] <- expr.filt.matrices[["ntap"]][common.genes, ]
  test.genomic.args[nxt.indx] <- list(NULL)
  test.clinical.args[nxt.indx] <- list(NULL)
  test.common.drugs[[nxt.indx]] <- common.drugs.by.data.set[["ntap"]]
  test.drug.cols[[nxt.indx]] <- drug.name.col
  test.patient.cols[[nxt.indx]] <- patient.id.cols[["ntap"]]
  test.response.cols[[nxt.indx]] <- "dss.auc.ll4"
} ## if("ntap" %in% data.sets.to.analyze) 

ds <- "sanger"
if(ds %in% data.sets.to.analyze) {

  nxt.indx <- length(test.set.names) + 1
  test.set.names[[nxt.indx]] <- ds
  test.dss.args[[nxt.indx]] <- fits[[ds]]
  test.expr.args[[nxt.indx]] <- expr.filt.matrices[[ds]][common.genes, ]
  test.genomic.args[nxt.indx] <- list(NULL)
  test.clinical.args[nxt.indx] <- list(NULL)
  test.common.drugs[[nxt.indx]] <- common.drugs.by.data.set[[ds]]
  test.drug.cols[[nxt.indx]] <- drug.name.col
  test.patient.cols[[nxt.indx]] <- patient.id.cols[[ds]]
  test.response.cols[[nxt.indx]] <- "dss.auc.ll4"
} ## if("sanger" %in% data.sets.to.analyze) 

num.processes <- num.cores - 1

cat("Training and testing with AUC\n")
train.response.cols <- rep("dss.auc.ll4", length(train.response.cols))
test.response.cols <- rep("dss.auc.ll4", length(test.response.cols))
res.auc <- train.and.test.crossproduct(gene.sets = gene.sets, drug.name.tbl = drug.name.tbl, train.set.names = train.set.names, train.dss.args = train.dss.args, 
                                       train.expr.args = train.expr.args, train.genomic.args = train.genomic.args, train.clinical.args = train.clinical.args, 
                                       train.common.drugs = train.common.drugs, train.drugs = train.drugs, train.drug.cols = train.drug.cols, train.patient.cols = train.patient.cols, 
                                       train.response.cols = train.response.cols, test.set.names = test.set.names, test.dss.args = test.dss.args, test.expr.args = test.expr.args, 
                                       test.genomic.args = test.genomic.args, test.clinical.args = test.clinical.args, test.common.drugs = test.common.drugs,
                                       test.drug.cols = test.drug.cols, test.patient.cols = test.patient.cols, test.response.cols = test.response.cols, 
                                       seed = 1234, use.rf = FALSE, use.svm = FALSE, num.processes = num.processes) 
## save.image(".Rdata")
cat("Done training and testing with AUC\n")

save.image(".Rdata.focs")

do.shifted.auc <- FALSE
if(do.shifted.auc) {
train.response.cols <- rep("shifted.auc.ll4", length(train.response.cols))
test.response.cols <- rep("shifted.auc.ll4", length(test.response.cols))
res.shifted.auc <- train.and.test.crossproduct(gene.sets = gene.sets, drug.name.tbl = drug.name.tbl, train.set.names = train.set.names, train.dss.args = train.dss.args, 
                                       train.expr.args = train.expr.args, train.genomic.args = train.genomic.args, train.clinical.args = train.clinical.args, 
                                       train.common.drugs = train.common.drugs, train.drugs = train.drugs, train.drug.cols = train.drug.cols, train.patient.cols = train.patient.cols, 
                                       train.response.cols = train.response.cols, test.set.names = test.set.names, test.dss.args = test.dss.args, test.expr.args = test.expr.args, 
                                       test.genomic.args = test.genomic.args, test.clinical.args = test.clinical.args, test.common.drugs = test.common.drugs,
                                       test.drug.cols = test.drug.cols, test.patient.cols = test.patient.cols, test.response.cols = test.response.cols, 
                                       seed = 1234, use.rf = FALSE, use.svm = FALSE, num.processes = num.processes) 
## save.image(".Rdata")
cat("Done training and testing with shifted AUC\n")

save.image(".Rdata.focs")
}

## Now aggregate all of those plots into a violin (showing points)
## See plot.correlation.vs.model for how to create a violin plot

tbl <- ldply(train.set.names,
             .fun = function(train.set) {
                      ldply(test.set.names, 
                            .fun = function(test.set) {
                                     ldply(names(gene.sets),
                                           .fun = function(gene.set) {
                                                    tmp.auc <- extract.results(res.auc[["all.comparisons"]][[train.set]][[test.set]][[gene.set]])
                                                    tmp.auc$train.set <- train.set
                                                    tmp.auc$test.set <- test.set
                                                    tmp.auc$gene.set <- gene.set
                                                    tmp.auc$response <- "AUC"
                                                    if(do.shifted.auc) {
                                                      tmp.shifted.auc <- extract.results(res.shifted.auc[["all.comparisons"]][[train.set]][[test.set]][[gene.set]])
                                                      tmp.shifted.auc$train.set <- train.set
                                                      tmp.shifted.auc$test.set <- test.set
                                                      tmp.shifted.auc$gene.set <- gene.set
                                                      tmp.shifted.auc$response <- "shifted-AUC"
                                                      return(rbind(tmp.auc, tmp.shifted.auc))
                                                    }
                                                    return(tmp.auc)                         
                                                  }) 
                                   }) 
                    })

## For each train-test combination create a vector of correlations (one for each drug) of predicted vs actual
corr.lists <- list()
pval.lists <- list()
for(test.set.name in unique(tbl$test.set)) {
  corr.lists[[test.set.name]] <- list()
  pval.lists[[test.set.name]] <- list()
  for(mdl in unique(tbl$model)) {
    corr.lists[[test.set.name]][[mdl]] <- list()
    pval.lists[[test.set.name]][[mdl]] <- list()
    for(resp in unique(tbl$response)) {
      corr.lists[[test.set.name]][[mdl]][[resp]] <- list()
      pval.lists[[test.set.name]][[mdl]][[resp]] <- list()
      mets <- c("pearson", "spearman")
      for(met in mets[mets %in% tbl$metric]) {
        corr.lists[[test.set.name]][[mdl]][[resp]][[met]] <- list()
        pval.lists[[test.set.name]][[mdl]][[resp]][[met]] <- list()
        for(gs in unique(tbl$gene.set)) {
          corr.lists[[test.set.name]][[mdl]][[resp]][[met]][[gs]] <- list()
          pval.lists[[test.set.name]][[mdl]][[resp]][[met]][[gs]] <- list()
          t.sub <- subset(tbl, (test.set == test.set.name) & (model == mdl) & (response == resp) & (gene.set == gs) & (metric == met))
          na.drugs <- t.sub$test.drug[is.na(t.sub$val)]
          t.sub <- t.sub[!(t.sub$test.drug %in% na.drugs),]
          test.drugs <- unique(t.sub$test.drug)
          for(train.set.name in unique(tbl$train.set)) {
            tr.sub <- subset(t.sub, train.set == train.set.name) 
            vec <- tr.sub$val
            names(vec) <- tr.sub$test.drug
            pvec <- tr.sub$pval
            names(pvec) <- tr.sub$test.drug
            corr.lists[[test.set.name]][[mdl]][[resp]][[met]][[gs]][[train.set.name]] <- vec[test.drugs]
            pval.lists[[test.set.name]][[mdl]][[resp]][[met]][[gs]][[train.set.name]] <- pvec[test.drugs]
          }
        }
      }
    }
  }
}

cat("Created correlation list\n")
save.image(".Rdata.focs")

capwords <- function(s, strict = FALSE) {
    cap <- function(s) paste(toupper(substring(s, 1, 1)),
                  {s <- substring(s, 2); if(strict) tolower(s) else s},
                             sep = "", collapse = " " )
    sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}

## Make a plot of correlation vs gene, biocart, hallmark, kegg, and reactome 
tbl$gene.set <- factor(tbl$gene.set, levels = c("gene", "biocarta", "hallmark", "kegg", "reactome"))
tbl$train.set <- factor(tbl$train.set, levels = c("ctrp", "ctrp.all", "ctrp.heme", "ctrp.aml", "ohsu.train", "ntap"))
for(resp in unique(tbl$response)) {
  for(met in c("pearson", "spearman")) {
    for(mdl in c("ridge", "lasso")) {
      sub <- subset(tbl, metric == met & model == mdl & response == resp)
      g <- ggplot(data = sub, aes(x = gene.set, y = val))
      g <- g + facet_grid(train.set ~ test.set)
      g <- g + geom_violin()
##    g <- g + geom_beeswarm()
      g <- g + geom_boxplot(width = 0.5)
      g <- g + geom_hline(yintercept = 0, linetype = "dashed")
      g <- g + ylab(paste0(capwords(met), " Correlation"))
      g <- g + ggtitle(paste0(resp, " ", capwords(met), " ", mdl))
      g <- g + ylim(c(-1,1))
      file <- paste0(file.prefix, "-", resp, "-", mdl, "-", met, "-vs-gene-sets.pdf")
      pdf(file)
      print(g)
      d <- dev.off()
##    glist[[length(glist)+1]] <- g
    }
  }
##  do.call("grid.arrange", glist)
}

for(resp in unique(tbl$response)) {
  for(met in c("pearson", "spearman")) {
    for(mdl in c("ridge", "lasso")) {
      file <- paste0(file.prefix, "-", resp, "-", mdl, "-", met, "-vs-gene-sets-volcano.pdf")
      pdf(file)
      for(gs in unique(tbl$gene.set)) {
        sub <- subset(tbl, metric == met & model == mdl & response == resp & gene.set == gs)
        tmp <- sub
        pvalue <- rep(">= 0.05", length(tmp$pval))
        pvalue[!is.na(tmp$pval) & (tmp$pval < 0.05)] <- "< 0.05"
        tmp$pvalue <- pvalue
        tmp$pval <- -log10(tmp$pval)
        flag <- sub$train.set == "ohsu.train" & sub$test.set == "ohsu.test"
##        if(any(flag)) { tmp$pval[flag] <- tmp$pval[flag] / max(tmp$pval[flag], na.rm=TRUE) }
        if(any(flag)) { tmp$pval[flag] <- tmp$pval[flag] / 10 }
        g <- ggplot(data = tmp, aes(x = val, y = pval, colour = pvalue))
        g <- g + facet_grid(train.set ~ test.set)
        g <- g + geom_point()
        g <- g + xlab(paste0(capwords(met), " Correlation"))
        g <- g + geom_hline(yintercept = -log10(0.05), linetype = "dashed")
        g <- g + geom_vline(xintercept = 0, linetype = "dashed")
        g <- g + ylab("-Log10 p-value")
        g <- g + ggtitle(paste0(resp, " ", capwords(met), " ", mdl, " ", gs))
        print(g)
      }
      d <- dev.off()
##    glist[[length(glist)+1]] <- g
    }
  }
##  do.call("grid.arrange", glist)
}

panel.diff <- function (x, y, col = par("col"), bg = NA, pch = par("pch"), 
                              cex = 1, col.regres = "red", ...) 
{ 
  ok <- is.finite(x) & is.finite(y) 
print(x-y)
vec <- x - y
  points(1:length(vec), vec, pch = pch, col = col, bg = bg, cex = cex) 
##  abline(h = 0, col = col.regres, ...) 
} 

save.image(".Rdata.focs")

if(length(unique(tbl$train.set)) > 1) { 
## Create pairwise correlation plots
for(test.set.name in unique(tbl$test.set)) {
  for(mdl in c("lasso", "ridge")) {
    for(resp in unique(tbl$response)) {
      pdf(paste(file.prefix, test.set.name, mdl, resp, "corr.pdf", sep="-"))
      mets <- c("pearson", "spearman")
      for(met in mets[mets %in% tbl$metric]) {
        for(gs in unique(tbl$gene.set)) {
          mat <- do.call("cbind", corr.lists[[test.set.name]][[mdl]][[resp]][[met]][[gs]])
          cols <- NULL
          if("Class.explained" %in% colnames(drug.name.tbl)) {
            cls <- drug.name.tbl$Class.explained
            names(cls) <- drug.name.tbl[, drug.name.col]
            cols <- factor(cls[names(corr.lists[[test.set.name]][[mdl]][[resp]][[met]][[gs]][[1]])])
          }
##        pdf(paste(test.set.name, mdl, resp, gs, "corr.pdf", sep="-"))
          pairs(mat, lower.panel=panel.regression, upper.panel=panel.cor, main = paste(test.set.name, mdl, resp, met, gs, sep=" "), col = cols, oma=c(4,4,6,12))
          par(xpd=TRUE)
          legend(0.85, 0.7, legend = levels(cols), fill = c(1:length(levels(cols))), cex = 0.7)
##        d <- dev.off()
        }
      }
      d <- dev.off()
    }
  }
}

for(test.set.name in unique(tbl$test.set)) {
  for(mdl in c("lasso", "ridge")) {
    for(resp in unique(tbl$response)) {
      pdf(paste(file.prefix, test.set.name, mdl, resp, "diff.pdf", sep="-"))
      mets <- c("pearson", "spearman")
      for(met in mets[mets %in% tbl$metric]) {
        for(gs in unique(tbl$gene.set)) {
          mat <- do.call("cbind", corr.lists[[test.set.name]][[mdl]][[resp]][[met]][[gs]])
          cols <- NULL
          if("Class.explained" %in% colnames(drug.name.tbl)) {
            cls <- drug.name.tbl$Class.explained
            names(cls) <- drug.name.tbl[, drug.name.col]
            cols <- factor(cls[names(corr.lists[[test.set.name]][[mdl]][[resp]][[met]][[gs]][[1]])])
          }
          len <- length(corr.lists[[test.set.name]][[mdl]][[resp]][[met]][[gs]])
          glist <- list()
          for(i in 1:len) {
            for(j in 1:len) {
              g <- NULL
              if(i != j) {
                x <- y <- NULL
                if(i > j) {
                  x <- corr.lists[[test.set.name]][[mdl]][[resp]][[met]][[gs]][[i]]
                  y <- corr.lists[[test.set.name]][[mdl]][[resp]][[met]][[gs]][[j]]
                  df <- data.frame(x = x, y = y, col = cols)
                } else {
                  x <- corr.lists[[test.set.name]][[mdl]][[resp]][[met]][[gs]][[j]]
                  y <- corr.lists[[test.set.name]][[mdl]][[resp]][[met]][[gs]][[i]]
                  df <- data.frame(x = 1:length(x), y = x - y, col = cols)
                }
                if(i < j) { g <- ggplot(df, aes(x = x, y = y, colour = col)) } else { g <- ggplot(df, aes(x = x, y = y)) }
                g <- g + geom_point()
                g <- g + theme(legend.position="none")
                if(i < j) {
                  g <- g + geom_hline(yintercept = 0)
                } else {
                  g <- g + geom_smooth(method='lm')
                  lims <- c(min(min(df$x, na.rm=TRUE), min(df$y, na.rm=TRUE)), max(max(df$x, na.rm=TRUE), max(df$y, na.rm=TRUE)))
                  g <- g + ggtitle(lm_eqn_bq(df))
                  g <- g + geom_abline(intercept = 0, slope = 1, linetype="dashed")
                  ## g <- g + geom_text(x = lims[1] + 0.2, y = lims[1] + 0.1, label = lm_eqn(df), parse=TRUE)
                }
              } else {
                g <- ggplot() + annotation_custom(textGrob(names(corr.lists[[test.set.name]][[mdl]][[resp]][[met]][[gs]])[i]),
                                 xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) 
              }
              g <- g + xlab("") + ylab("")
              if(i < j) {
                g <- g + ylab(paste0(names(corr.lists[[test.set.name]][[mdl]][[resp]][[met]][[gs]])[j], " - ", names(corr.lists[[test.set.name]][[mdl]][[resp]][[met]][[gs]])[i]))
              }
              if(i > j) {
                g <- g + xlab(names(corr.lists[[test.set.name]][[mdl]][[resp]][[met]][[gs]])[i])
                g <- g + ylab(names(corr.lists[[test.set.name]][[mdl]][[resp]][[met]][[gs]])[j])
              }
  
              glist[[length(glist)+1]] <- g
            }
          }
          do.call("grid.arrange", c(glist, "top" = paste(test.set.name, mdl, resp, met, gs, sep=" ")))
        }
      }
      d <- dev.off()
    }
  }
}

} ## if(length(unique(tbl$train.set)) > 1) 

stop("successfully completed")
