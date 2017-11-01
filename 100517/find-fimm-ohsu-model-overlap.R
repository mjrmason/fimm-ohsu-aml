
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

num.cores <- detectCores()
if(!is.na(num.cores) && (num.cores > 1)) {
  suppressPackageStartupMessages(library("doMC"))
  cat(paste("Registering ", num.cores-1, " cores.\n", sep=""))
  registerDoMC(cores=(num.cores-1))
}

synapseLogin()

cat("Compare features in common between FIMM and OHSU\n")

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

## BEGIN read in preprocessed OHSU, FIMM, and CTRVP2 data from Synapse


## Update as change data sets
data.sets.to.analyze <- c("ohsu", "fimm")

## Update as change data sets
processed.file.names <- c("ohsu.fo.common.drugs.dss.t0.tsv", "fimm.fo.common.drugs.dss.t0.tsv")
names(processed.file.names) <- c("ohsu", "fimm")
## parentId is the processed data folder of the respective data sets
parentIds <- c("syn10083332", "syn8270577")
names(parentIds) <- c("ohsu", "fimm")

## Read in the processed DSS files and store in dss.orig list
dss.orig <- list()
for(i in 1:length(processed.file.names)) {
  parentId <- parentIds[i]
  tbl <- synQuery(paste0("SELECT name,id FROM file WHERE parentId==\'", parentId, "\'"))
  file.name <- processed.file.names[i]
  synId <- tbl[tbl$file.name == file.name, "file.id"]
  obj <- synGet(synId, downloadFile = TRUE)
  file <- getFileLocation(obj)
  ## Read in the t=0 thresholded DSS values
  dss.orig[[names(processed.file.names)[i]]] <- as.data.frame(fread(file))
}

## Read in the expression matrices (possibly post-combat corrected)

## Update as change data sets
## NB: these have been combat-corrected for the comparison of interest.  
## e.g., if comparing OHSU and FIMM, only these 2 should have been input to combat.
## If comparing OHSU, FIMM, and CTRP AML, these 3 should have been input to combat.
## Below are the raw expression data--before being corrected by combat.
## raw.paths <- c("ohsu.expr.tsv", "RNASeq.CPM.log2.bc.csv", "CCLE_RNA-seq_tpm_matrix.csv")
## raw.expression.synIds <- c("syn10083723", "syn8270602", "syn5616092")
## paths <- c("ohsu-expr-fo-combat.tsv", "fimm-expr-fo-combat.tsv")
expression.synIds <- c("syn11055512", "syn11055519")
names(expression.synIds) <- c("ohsu", "fimm")

expr.matrices <- list()

for(i in 1:length(expression.synIds)) {
  synId <- expression.synIds[i]
  obj <- synGet(synId, downloadFile = TRUE)
  file <- getFileLocation(obj)
  expr <- read.table(file, header=TRUE, sep="\t")
  switch(names(expression.synIds)[i],
    "ohsu" = {
      ## Convert the column names from X20.00347 to 20-00347
      colnames(expr) <- gsub("X(.*)\\.(.*)", "\\1-\\2", colnames(expr))
    },
    "fimm" = {
    },
    "ctrp" = {
    },
    {
      die("Unrecognized data set\n")
    })
  expr.matrices[[names(expression.synIds)[i]]] <- expr
}

## Read in the OHSU diagnosis info in the raw drug file and drug outcome file
## and limit samples to aml
if("ohsu" %in% data.sets.to.analyze) {
  aml.diagnosis <- "ACUTE MYELOID LEUKAEMIA (AML) AND RELATED PRECURSOR NEOPLASMS"

  ## Read in the OHSU patient diagnosis table (corresponding to the drug response data).
  ## Diagnosis_Labs_Treatments_Outcomes_2017_01_12.xlsx (syn8149174)
  synId <- "syn8149174"
  obj <- synGet(synId, downloadFile=TRUE)
  ohsu.diagnosis.tbl <- openxlsx:::read.xlsx(getFileLocation(obj), sheet=1)
  ohsu.diagnosis.tbl <- subset(ohsu.diagnosis.tbl, most_recent_diagnosis == aml.diagnosis)
  ohsu.diagnosis.tbl <- subset(ohsu.diagnosis.tbl, diagnosis_at_time_of_specimen_acquisition == aml.diagnosis)

  ## Read in the raw OHSU inhibitor data 
  ## inhibitor_data_points_2017_01_12.txt
  ## This latest file seems to have data from the prior releases
  synId <- "syn8149180"
  inh.obj <- synGet(synId, downloadFile=TRUE)
  ohsu.inh.tbl <- fread(getFileLocation(inh.obj))
  ohsu.inh.tbl <- as.data.frame(ohsu.inh.tbl)
  ohsu.inh.tbl <- subset(ohsu.inh.tbl, diagnosis == aml.diagnosis)
  ohsu.inh.tbl <- unique(ohsu.inh.tbl[, c("patient_id", "lab_id", "inhibitor", "diagnosis")])

  ## Limit the drugs to those that come from aml samples
  dss.orig[["ohsu"]] <- subset(dss.orig[["ohsu"]], lab_id %in% ohsu.diagnosis.tbl$lab_id)
  dss.orig[["ohsu"]] <- subset(dss.orig[["ohsu"]], lab_id %in% ohsu.inh.tbl$lab_id)

  ## Add the SeqID to the OHSU drug data.
  ## Read in the rna-seq summary table, which provides a map from the lab_id's used
  ## in the drug response data to the SeqID's used in the expression data.
  synId <- "syn10083817"
  obj <- synGet(synId, downloadFile = TRUE)
  file <- getFileLocation(obj)
  ohsu.rnaseq.sample.summary <- read.table(file, header=TRUE, sep="\t")

  dss.orig[["ohsu"]] <- merge(dss.orig[["ohsu"]], ohsu.rnaseq.sample.summary, by.x = "lab_id", by.y = "Original_LabID")
}

if("fimm" %in% data.sets.to.analyze) {
  ## Read in the FIMM Metadata and ensure that we only analyze AML data
  synId <- "syn8270594"
  obj <- synGet(id=synId, downloadFile = TRUE)
  file <- getFileLocation(obj)
  fimm.metadata <- read.table(file, header=TRUE, sep=",")
  fimm.metadata <- subset(fimm.metadata, grepl(diagnosis, pattern="AML"))

  expr.matrices[["fimm"]] <- expr.matrices[["fimm"]][, colnames(expr.matrices[["fimm"]]) %in% rownames(fimm.metadata)]
}

## Read in the CTRP cell line metadata so we can limit (or not) cell lines by disease type
## CTRP cell line metadata: "v20.meta.per_cell_line.txt" (syn5632192)
all.ccls <- NULL
heme.ccls <- NULL
aml.ccls <- NULL
if("ctrp" %in% data.sets.to.analyze) {
  synId <- "syn5632192"
  obj <- synGet(synId, downloadFile = TRUE)
  ctrp.cell.line.metadata <- read.table(getFileLocation(obj), sep="\t", header=TRUE, as.is=TRUE, comment.char="", quote="")

  ## Samples in CTRP do not have "."s and are all in caps
  colnames(expr.matrices[["ctrp"]]) <- gsub(x=colnames(expr.matrices[["ctrp"]]), pattern="\\.", replacement="")
  colnames(expr.matrices[["ctrp"]]) <- toupper(colnames(expr.matrices[["ctrp"]]))
  dss.orig[["ctrp"]] <- merge(dss.orig[["ctrp"]], ctrp.cell.line.metadata[, c("master_ccl_id", "ccl_name")])
  ## all.ccls <- intersect(colnames(expr.matrices[["ctrp"]]), dss.orig[["ctrp"]][, "ccl_name"])
  all.ccls <- unique(dss.orig[["ctrp"]][, "ccl_name"])
  heme.ccls <- all.ccls[all.ccls %in% ctrp.cell.line.metadata[ctrp.cell.line.metadata$ccle_primary_hist == "haematopoietic_neoplasm", "ccl_name"]]
  aml.ccls <- all.ccls[all.ccls %in% ctrp.cell.line.metadata[ctrp.cell.line.metadata$ccle_hist_subtype_1 == "acute_myeloid_leukaemia", "ccl_name"]]
}

## END read in preprocessed OHSU, FIMM, and CTRVP2 data from Synapse

## BEGIN define drugs common to FIMM, OHSU, and CTRP

## Subset to analysis to drugs common to FIMM, OHSU, and CTRP
cols <- c()
new.col.names <- c()
if("ohsu" %in% data.sets.to.analyze) {
  cols <- c(cols, "inhibitor")
  new.col.names <- c(new.col.names, "inhibitor")
}
if("fimm" %in% data.sets.to.analyze) {
  cols <- c(cols, c("FIMM_Batch_ID_Drug.ll4", "FIMM_DRUG_NAME.ll4"))
  new.col.names <- c(new.col.names, "DRUG_ID", "DRUG_NAME")
}
if("ctrp" %in% data.sets.to.analyze) {
  cols <- c(cols, "master_cpd_id.ll4")
  new.col.names <- c(new.col.names, "master_cpd_id")
}

common.drugs <- dss.orig[[names(dss.orig)[1]]][, cols]
colnames(common.drugs) <- new.col.names

if("ohsu" %in% data.sets.to.analyze) {
  common.drugs <- subset(common.drugs, inhibitor %in% dss.orig[["ohsu"]][, "inhibitor"])
}
if("fimm" %in% data.sets.to.analyze) {
  common.drugs <- subset(common.drugs, DRUG_ID %in% dss.orig[["fimm"]][, "DRUG_ID"])
}
if("ctrp" %in% data.sets.to.analyze) {
  common.drugs <- subset(common.drugs, master_cpd_id %in% dss.orig[["ctrp"]][, "master_cpd_id"])
}

drug.name.tbl <- unique(common.drugs)
common.drugs.by.data.set <- list()
if("ohsu" %in% data.sets.to.analyze) {
  common.drugs.by.data.set[["ohsu"]] <- unique(drug.name.tbl[, "inhibitor"])
}
if("fimm" %in% data.sets.to.analyze) {
  common.drugs.by.data.set[["fimm"]] <- unique(drug.name.tbl[, "DRUG_ID"])
}
if("ctrp" %in% data.sets.to.analyze) {
  common.drugs.by.data.set[["ctrp"]] <- unique(drug.name.tbl[, "master_cpd_id"])
}

## Download FIMM_OHSU_Drug_Concentrations.xlsx (which has drug targets)
if("ohsu" %in% data.sets.to.analyze) {
  synId <- "syn10163669"
  drug.target.tbl <- openxlsx:::read.xlsx(getFileLocation(synGet(synId)), sheet=1)
  tmp <- drug.target.tbl[, c("OHSU_DRUG_NAME", "Gene.Targets", "Class.explained", "Mechanism.Targets")]
  colnames(tmp) <- c("OHSU_DRUG_NAME", "Gene.Targets", "Class", "Mechanism")
  drug.name.tbl <- merge(drug.name.tbl, tmp, by.x = "inhibitor", by.y = "OHSU_DRUG_NAME")
}

## END define drugs common to FIMM, OHSU, and CTRP

## BEGIN modeling setup

source("../common/data-preprocessing.R")
source("models.R")
source("../common/plotting.R")

## Define CTRP AML, Heme, and All cell lines.  AML is a subset of Heme is a subset of all.
all.expr.cols <- NULL
heme.expr.cols <- NULL
aml.expr.cols <- NULL
if("ctrp" %in% data.sets.to.analyze) {
  all.expr.cols <- intersect(colnames(expr.matrices[["ctrp"]]), all.ccls)
  heme.expr.cols <- intersect(colnames(expr.matrices[["ctrp"]]), heme.ccls)
  aml.expr.cols <- intersect(colnames(expr.matrices[["ctrp"]]), aml.ccls)

  ## Log-transform the CTRP data (FIMM and OHSU are already in log space)
  min.expr <- min(expr.matrices[["ctrp"]][expr.matrices[["ctrp"]] != 0])
  ctrp.log2.expr <- expr.matrices[["ctrp"]]
  ctrp.log2.expr[ctrp.log2.expr == 0] <- min.expr
  ctrp.log2.expr <- log2(ctrp.log2.expr)
  expr.matrices[["ctrp"]] <- ctrp.log2.expr
  rm(ctrp.log2.expr)
  gc()
}

## Filter the gene expression data to restrict to highly-expressed genes
expr.filt.matrices <- list()
for(ds in names(expr.matrices)) {
  iqrs <- unlist(apply(expr.matrices[[ds]], 1, IQR))
  expr.filt.matrices[[ds]] <- expr.matrices[[ds]][iqrs > median(iqrs[iqrs > 0]),]
}

study.genes <- (lapply(expr.matrices, rownames))
study.filt.genes <- (lapply(expr.filt.matrices, rownames))
common.genes <- Reduce(intersect, study.filt.genes)

cat(paste0(length(common.genes), " common across studies.\n"))

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
if("ctrp" %in% data.sets.to.analyze) {
  ctrp.drug.col <- "master_cpd_id"
  ctrp.patient.col <- "ccl_name"

  ## Define the drug response variable.
  ## Use the AUC value computing using a log-logistic (LL.4) fit.  
  ## "dss." implies that AUC was calculated over the concentration range common
  ## to all three data sets.
  train.response.col <- "dss.auc.ll4"

  boot.res <- list()

  cat("Bootstrapping CTRP AML\n")
  ret <- bootstrap.model(train.dss.arg = dss.orig[["ctrp"]], train.expr.arg = expr.filt.matrices[["ctrp"]][common.genes, aml.expr.cols], train.common.drugs = common.drugs.by.data.set[["ctrp"]], train.drugs = common.drugs.by.data.set[["ctrp"]], 
                         train.drug.col = ctrp.drug.col, train.patient.col = ctrp.patient.col, train.response.col = train.response.col, num.bootstraps = 100)
  boot.res[["ctrp.aml"]] <- ret
  save.image(".Rdata")

  cat("Bootstrapping CTRP ALL\n")
  ret <- bootstrap.model(train.dss.arg = dss.orig[["ctrp"]], train.expr.arg = expr.filt.matrices[["ctrp"]][common.genes, ], train.common.drugs = common.drugs.by.data.set[["ctrp"]], train.drugs = common.drugs.by.data.set[["ctrp"]], 
                         train.drug.col = ctrp.drug.col, train.patient.col = ctrp.patient.col, train.response.col = train.response.col, num.bootstraps = 100)
  boot.res[["ctrp.all"]] <- ret
  save.image(".Rdata")

  nms <- unlist(lapply(boot.res[["ctrp.aml"]], function(x) x[[1]]$train.drug))
  names(boot.res[["ctrp.aml"]]) <- as.character(nms)
  names(nms) <- nms
  ctrp.aml.ranks <- llply(nms,
                          .parallel = TRUE,
                          .fun = function(indx) {
                                   lst <- boot.res[["ctrp.aml"]][[as.character(indx)]]
                                   flag <- unlist(lapply(lst, function(x) !is.null(x$coeffs) && nrow(x$coeffs) > 0))
                                   lst <- lst[flag]
                                   foo <- lapply(lst, function(x) { df <- data.frame(coeff = as.numeric(x$coeffs), gene = rownames(x$coeffs)); df$rank <- rank(- 1 * abs(df$coeff)); df })
                                   bar <- Reduce(function(x, y) merge(x, y, all=TRUE, by = "gene"), foo)
                                   bar$gene <- as.character(bar$gene)
                                   baz <- ldply(1:nrow(bar), .fun = function(i) { gene <- bar$gene[i]; rank.cols <- which(grepl(pattern="rank", colnames(bar))); grand.rank <- mean(as.numeric(bar[i,rank.cols]), na.rm=TRUE); c(gene, grand.rank) })
                                   colnames(baz) <- c("gene", "grand.rank")
                                   baz$grand.rank <- as.numeric(baz$grand.rank)
                                   baz <- baz[order(baz$grand.rank),]
                                   baz
                                 })

  nms <- unlist(lapply(boot.res[["ctrp.all"]], function(x) x[[1]]$train.drug))
  names(boot.res[["ctrp.all"]]) <- as.character(nms)
  names(nms) <- nms
  ctrp.all.ranks <- llply(nms,
                          .fun = function(indx) {
                                   lst <- boot.res[["ctrp.all"]][[as.character(indx)]]
                                   flag <- unlist(lapply(lst, function(x) !is.null(x$coeffs) && nrow(x$coeffs) > 0))
                                   lst <- lst[flag]
                                   foo <- lapply(lst, function(x) { df <- data.frame(coeff = as.numeric(x$coeffs), gene = rownames(x$coeffs)); df$rank <- rank(- 1 * abs(df$coeff)); df })
                                   bar <- Reduce(function(x, y) merge(x, y, all=TRUE, by = "gene"), foo)
                                   bar$gene <- as.character(bar$gene)
                                   baz <- ldply(1:nrow(bar), .fun = function(i) { gene <- bar$gene[i]; rank.cols <- which(grepl(pattern="rank", colnames(bar))); grand.rank <- mean(as.numeric(bar[i,rank.cols]), na.rm=TRUE); c(gene, grand.rank) })
                                   colnames(baz) <- c("gene", "grand.rank")
                                   baz$grand.rank <- as.numeric(baz$grand.rank)
                                   baz <- baz[order(baz$grand.rank),]
                                   baz
                                 })
## save.image(".Rdata")

  ctrp.aml.ranks.anno <- annotate.ranked.list(ctrp.aml.ranks)
  ctrp.all.ranks.anno <- annotate.ranked.list(ctrp.all.ranks)

  ## HALLMARK
  ## http://software.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/6.0/h.all.v6.0.symbols.gmt
  synId <- "syn10507487"
  obj <- synGet(synId, downloadFile = TRUE)
  file <- getFileLocation(obj)
  gsea.gene.sets <- GSA.read.gmt(file)

  ## BIOCARTA
  ## http://software.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/6.0/c2.cp.biocarta.v6.0.symbols.gmt
  synId <- "syn10507483"
  obj <- synGet(synId, downloadFile = TRUE)
  file <- getFileLocation(obj)
  gsea.gene.sets <- GSA.read.gmt(file)

  ctrp.aml.anno <- merge(drug.name.tbl, ctrp.aml.ranks.anno, by.x = "master_cpd_id", by.y = "drug")
  ctrp.all.anno <- merge(drug.name.tbl, ctrp.all.ranks.anno, by.x = "master_cpd_id", by.y = "drug")

  foo <- annotate.with.overlapping.gene.sets(ctrp.aml.anno, gsea.gene.sets, "biocarta", "genes.symbol", "Gene.Targets")
  foo$overlapping.targets <- unlist(apply(foo[, c("genes.symbol", "Gene.Targets")], 1, 
                                          function(row) {
                                            gs1 <- unlist(strsplit(row[1], split=",[ ]*"))
                                            gs2 <- unlist(strsplit(row[2], split=",[ ]*"))
                                            overlap <- ""
                                            if(!is.na(gs1) && !is.na(gs2) && any(gs1 %in% gs2)) { overlap <- gs1[gs1 %in% gs2] }
                                            overlap
                                          }))
  foo <- foo[, c("DRUG_NAME", "Mechanism", "Class", "Gene.Targets", "genes.symbol", "ranks", "genes.ensg", "biocarta", "overlapping.targets")]
  ctrp.aml.anno <- foo

  foo <- annotate.with.overlapping.gene.sets(ctrp.all.anno, gsea.gene.sets, "biocarta", "genes.symbol", "Gene.Targets")
  foo$overlapping.targets <- unlist(apply(foo[, c("genes.symbol", "Gene.Targets")], 1, 
                                          function(row) {
                                            gs1 <- unlist(strsplit(row[1], split=",[ ]*"))
                                            gs2 <- unlist(strsplit(row[2], split=",[ ]*"))
                                            overlap <- ""
                                            if(!is.na(gs1) && !is.na(gs2) && any(gs1 %in% gs2)) { overlap <- gs1[gs1 %in% gs2] }
                                            overlap
                                          }))
  foo <- foo[, c("DRUG_NAME", "Mechanism", "Class", "Gene.Targets", "genes.symbol", "ranks", "genes.ensg", "biocarta", "overlapping.targets")]
  ctrp.all.anno <- foo

  write.table(file="ctrp.aml.anno.tsv", ctrp.aml.anno, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
  write.table(file="ctrp.all.anno.tsv", ctrp.all.anno, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
} ## if("ctrp" %in% data.sets.to.analyze) 

## Assemble training and test sets

## Split the OHSU data into training (60%) and test (40%) -- do this on a per-drug basis because of the
## missing data.
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

if("ohsu" %in% data.sets.to.analyze) {
  set.seed(1234)
  train.index <- createDataPartition(dss.orig[["ohsu"]][, "inhibitor"], p = .6, list = FALSE)
  ohsu.dss.orig.train <- dss.orig[["ohsu"]][train.index, ]
  ohsu.dss.orig.test <- dss.orig[["ohsu"]][-train.index, ]

  nxt.indx <- length(train.set.names) + 1
  train.set.names[[nxt.indx]] <- "ohsu.train"
  train.dss.args[[nxt.indx]] <- ohsu.dss.orig.train
  train.expr.args[[nxt.indx]] <- expr.filt.matrices[["ohsu"]][common.genes, ] 
  train.genomic.args[nxt.indx] <- list(NULL)
  train.clinical.args[nxt.indx] <- list(NULL)
  train.common.drugs[[nxt.indx]] <- common.drugs.by.data.set[["ohsu"]]
  train.drugs[[nxt.indx]] <- common.drugs.by.data.set[["ohsu"]]
  train.drug.cols[[nxt.indx]] <- ohsu.drug.col
  train.patient.cols[[nxt.indx]] <- ohsu.patient.col
  train.response.cols[[nxt.indx]] <- "dss.auc.ll4"

  nxt.indx <- length(test.set.names) + 1
  test.set.names[[nxt.indx]] <- "ohsu.test"
  test.dss.args[[nxt.indx]] <- ohsu.dss.orig.test
  test.expr.args[[nxt.indx]] <- expr.filt.matrices[["ohsu"]][common.genes, ]
  test.genomic.args[nxt.indx] <- list(NULL)
  test.clinical.args[nxt.indx] <- list(NULL)
  test.common.drugs[[nxt.indx]] <- common.drugs.by.data.set[["ohsu"]]
  test.drug.cols[[nxt.indx]] <- ohsu.drug.col
  test.patient.cols[[nxt.indx]] <- ohsu.patient.col
  test.response.cols[[nxt.indx]] <- "dss.auc.ll4"

  rm(ohsu.dss.orig.train)
  rm(ohsu.dss.orig.test)
} ## if("ohsu" %in% data.sets.to.analyze) 

fimm.drug.col <- "DRUG_ID"
fimm.patient.col <- "SCREEN_ID"

if("fimm" %in% data.sets.to.analyze) {

  nxt.indx <- length(test.set.names) + 1
  test.set.names[[nxt.indx]] <- "fimm"
  test.dss.args[[nxt.indx]] <- dss.orig[["fimm"]]
  test.expr.args[[nxt.indx]] <- expr.filt.matrices[["fimm"]][common.genes, ]
  test.genomic.args[nxt.indx] <- list(NULL)
  test.clinical.args[nxt.indx] <- list(NULL)
  test.common.drugs[[nxt.indx]] <- common.drugs.by.data.set[["fimm"]]
  test.drug.cols[[nxt.indx]] <- fimm.drug.col
  test.patient.cols[[nxt.indx]] <- fimm.patient.col
  test.response.cols[[nxt.indx]] <- "dss.auc.ll4"
} ## if("fimm" %in% data.sets.to.analyze) 

if("ctrp" %in% data.sets.to.analyze) {
  ctrp.drug.col <- "master_cpd_id"
  ctrp.patient.col <- "ccl_name"

  nxt.indx <- length(train.set.names) + 1
  train.set.names[[nxt.indx]] <- "ctrp.all"
  train.dss.args[[nxt.indx]] <- dss.orig[["ctrp"]]
  train.expr.args[[nxt.indx]] <- expr.filt.matrices[["ctrp"]][common.genes, ] 
  train.genomic.args[nxt.indx] <- list(NULL)
  train.clinical.args[nxt.indx] <- list(NULL)
  train.common.drugs[[nxt.indx]] <- common.drugs.by.data.set[["ctrp"]]
  train.drugs[[nxt.indx]] <- common.drugs.by.data.set[["ctrp"]]
  train.drug.cols[[nxt.indx]] <- ctrp.drug.col
  train.patient.cols[[nxt.indx]] <- ctrp.patient.col
  train.response.cols[[nxt.indx]] <- "dss.auc.ll4"

  nxt.indx <- length(train.set.names) + 1
  train.set.names[[nxt.indx]] <- "ctrp.heme"
  train.dss.args[[nxt.indx]] <- dss.orig[["ctrp"]]
  train.expr.args[[nxt.indx]] <- expr.filt.matrices[["ctrp"]][common.genes, heme.expr.cols] 
  train.genomic.args[nxt.indx] <- list(NULL)
  train.clinical.args[nxt.indx] <- list(NULL)
  train.common.drugs[[nxt.indx]] <- common.drugs.by.data.set[["ctrp"]]
  train.drugs[[nxt.indx]] <- common.drugs.by.data.set[["ctrp"]]
  train.drug.cols[[nxt.indx]] <- ctrp.drug.col
  train.patient.cols[[nxt.indx]] <- ctrp.patient.col
  train.response.cols[[nxt.indx]] <- "dss.auc.ll4"

  nxt.indx <- length(train.set.names) + 1
  train.set.names[[nxt.indx]] <- "ctrp.aml"
  train.dss.args[[nxt.indx]] <- dss.orig[["ctrp"]]
  train.expr.args[[nxt.indx]] <- expr.filt.matrices[["ctrp"]][common.genes, aml.expr.cols] 
  train.genomic.args[nxt.indx] <- list(NULL)
  train.clinical.args[nxt.indx] <- list(NULL)
  train.common.drugs[[nxt.indx]] <- common.drugs.by.data.set[["ctrp"]]
  train.drugs[[nxt.indx]] <- common.drugs.by.data.set[["ctrp"]]
  train.drug.cols[[nxt.indx]] <- ctrp.drug.col
  train.patient.cols[[nxt.indx]] <- ctrp.patient.col
  train.response.cols[[nxt.indx]] <- "dss.auc.ll4"

} ## if("ctrp" %in% data.sets.to.analyze) 

num.processes <- num.cores - 1

cat("Training with AUC\n")
empty.list <- list()

train.only.set.names <- list()
train.only.dss.args <- list()
train.only.expr.args <- list()
train.only.genomic.args <- list()
train.only.clinical.args <- list()
train.only.common.drugs <- list()
train.only.drugs <- list()
train.only.drug.cols <- list()
train.only.patient.cols <- list()
train.only.response.cols <- list()

## Train all ohsu, and two halfs of the ohsu data set
if("ohsu" %in% data.sets.to.analyze) {

  nxt.indx <- length(train.only.set.names) + 1
  train.only.set.names[[nxt.indx]] <- "ohsu"
  train.only.dss.args[[nxt.indx]] <- dss.orig[["ohsu"]]
  train.only.expr.args[[nxt.indx]] <- expr.filt.matrices[["ohsu"]][common.genes, ] 
  train.only.genomic.args[nxt.indx] <- list(NULL)
  train.only.clinical.args[nxt.indx] <- list(NULL)
  train.only.common.drugs[[nxt.indx]] <- common.drugs.by.data.set[["ohsu"]]
  train.only.drugs[[nxt.indx]] <- common.drugs.by.data.set[["ohsu"]]
  train.only.drug.cols[[nxt.indx]] <- ohsu.drug.col
  train.only.patient.cols[[nxt.indx]] <- ohsu.patient.col
  train.only.response.cols[[nxt.indx]] <- "dss.auc.ll4"

  set.seed(1234)
  half.index <- createDataPartition(dss.orig[["ohsu"]][, "inhibitor"], p = .6, list = FALSE)

  nxt.indx <- length(train.only.set.names) + 1
  train.only.set.names[[nxt.indx]] <- "ohsu.set1"
  train.only.dss.args[[nxt.indx]] <- dss.orig[["ohsu"]][half.index, ]
  train.only.expr.args[[nxt.indx]] <- expr.filt.matrices[["ohsu"]][common.genes, ] 
  train.only.genomic.args[nxt.indx] <- list(NULL)
  train.only.clinical.args[nxt.indx] <- list(NULL)
  train.only.common.drugs[[nxt.indx]] <- common.drugs.by.data.set[["ohsu"]]
  train.only.drugs[[nxt.indx]] <- common.drugs.by.data.set[["ohsu"]]
  train.only.drug.cols[[nxt.indx]] <- ohsu.drug.col
  train.only.patient.cols[[nxt.indx]] <- ohsu.patient.col
  train.only.response.cols[[nxt.indx]] <- "dss.auc.ll4"

  nxt.indx <- length(train.only.set.names) + 1
  train.only.set.names[[nxt.indx]] <- "ohsu.set2"
  train.only.dss.args[[nxt.indx]] <- dss.orig[["ohsu"]][-half.index, ]
  train.only.expr.args[[nxt.indx]] <- expr.filt.matrices[["ohsu"]][common.genes, ] 
  train.only.genomic.args[nxt.indx] <- list(NULL)
  train.only.clinical.args[nxt.indx] <- list(NULL)
  train.only.common.drugs[[nxt.indx]] <- common.drugs.by.data.set[["ohsu"]]
  train.only.drugs[[nxt.indx]] <- common.drugs.by.data.set[["ohsu"]]
  train.only.drug.cols[[nxt.indx]] <- ohsu.drug.col
  train.only.patient.cols[[nxt.indx]] <- ohsu.patient.col
  train.only.response.cols[[nxt.indx]] <- "dss.auc.ll4"
}

if("fimm" %in% data.sets.to.analyze) {
  nxt.indx <- length(train.only.set.names) + 1
  train.only.set.names[[nxt.indx]] <- "fimm"
  train.only.dss.args[[nxt.indx]] <- dss.orig[["fimm"]]
  train.only.expr.args[[nxt.indx]] <- expr.filt.matrices[["fimm"]][common.genes, ] 
  train.only.genomic.args[nxt.indx] <- list(NULL)
  train.only.clinical.args[nxt.indx] <- list(NULL)
  train.only.common.drugs[[nxt.indx]] <- common.drugs.by.data.set[["fimm"]]
  train.only.drugs[[nxt.indx]] <- common.drugs.by.data.set[["fimm"]]
  train.only.drug.cols[[nxt.indx]] <- fimm.drug.col
  train.only.patient.cols[[nxt.indx]] <- fimm.patient.col
  train.only.response.cols[[nxt.indx]] <- "dss.auc.ll4"

  set.seed(1234)
  half.index <- createDataPartition(dss.orig[["fimm"]][, "DRUG_ID"], p = .6, list = FALSE)

  nxt.indx <- length(train.only.set.names) + 1
  train.only.set.names[[nxt.indx]] <- "fimm.set1"
  train.only.dss.args[[nxt.indx]] <- dss.orig[["fimm"]][half.index, ]
  train.only.expr.args[[nxt.indx]] <- expr.filt.matrices[["fimm"]][common.genes, ] 
  train.only.genomic.args[nxt.indx] <- list(NULL)
  train.only.clinical.args[nxt.indx] <- list(NULL)
  train.only.common.drugs[[nxt.indx]] <- common.drugs.by.data.set[["fimm"]]
  train.only.drugs[[nxt.indx]] <- common.drugs.by.data.set[["fimm"]]
  train.only.drug.cols[[nxt.indx]] <- fimm.drug.col
  train.only.patient.cols[[nxt.indx]] <- fimm.patient.col
  train.only.response.cols[[nxt.indx]] <- "dss.auc.ll4"

  nxt.indx <- length(train.only.set.names) + 1
  train.only.set.names[[nxt.indx]] <- "fimm.set2"
  train.only.dss.args[[nxt.indx]] <- dss.orig[["fimm"]][-half.index, ]
  train.only.expr.args[[nxt.indx]] <- expr.filt.matrices[["fimm"]][common.genes, ] 
  train.only.genomic.args[nxt.indx] <- list(NULL)
  train.only.clinical.args[nxt.indx] <- list(NULL)
  train.only.common.drugs[[nxt.indx]] <- common.drugs.by.data.set[["fimm"]]
  train.only.drugs[[nxt.indx]] <- common.drugs.by.data.set[["fimm"]]
  train.only.drug.cols[[nxt.indx]] <- fimm.drug.col
  train.only.patient.cols[[nxt.indx]] <- fimm.patient.col
  train.only.response.cols[[nxt.indx]] <- "dss.auc.ll4"
}

train.response.cols <- rep("dss.auc.ll4", length(train.response.cols))
res.trn <- train.and.test.crossproduct(gene.sets = gene.sets, drug.name.tbl = drug.name.tbl, train.set.names = train.only.set.names, 
                                       train.dss.args = train.only.dss.args, train.expr.args = train.only.expr.args, 
                                       train.genomic.args = train.only.genomic.args, train.clinical.args = train.only.clinical.args, 
                                       train.common.drugs = train.only.common.drugs, train.drugs = train.only.drugs, train.drug.cols = train.only.drug.cols, 
                                       train.patient.cols = train.only.patient.cols, train.response.cols = train.only.response.cols, 
                                       test.set.names = empty.list, test.dss.args = empty.list, test.expr.args = empty.list, 
                                       test.genomic.args = empty.list, test.clinical.args = empty.list, test.common.drugs = empty.list,
                                       test.drug.cols = empty.list, test.patient.cols = empty.list, test.response.cols = empty.list, 
                                       seed = 1234, use.rf = FALSE, use.svm = FALSE, use.ridge = TRUE, use.mean = FALSE, num.processes = num.processes) 
save.image(".Rdata")

trn.tbl <- ldply(train.only.set.names,
             .fun = function(train.set) {
                      ldply(names(gene.sets),
                        .fun = function(gene.set) {
                                 ldply(res.trn[["all.fits"]][[train.set]][[gene.set]],
                                   .fun = function(lst) {
                                     ldply(lst, 
                                       .fun = function(df) {
                                         coeffs <- NULL
                                         coeffs <- coefficients(df$fit)
                                         ns <- rownames(coeffs)
                                         coeffs <- as.vector(coeffs)
                                         n.features <- length(coeffs) - 1 ## don't count intercept
                                         names(coeffs) <- ns
                                         coeffs <- coeffs[!grepl(pattern="Intercept", names(coeffs))]
                                         coeffs <- coeffs[order(names(coeffs))]
                                         coeff.vals <- paste(coeffs, collapse=",")
                                         coeffs <- paste(names(coeffs), collapse=",")
                                         ret <- data.frame(train.set = train.set, gene.set = gene.set, train.drug = df$train.drug,
                                                           alpha = df$alpha, model = df$model, n.features = n.features, coeffs = coeffs, coeff.vals = coeff.vals)
                                         ret
                                       })
                                   })
                        })
             })
trn.tbl$model <- as.character(trn.tbl$model)
trn.tbl[(trn.tbl$model == "glmnet") & (trn.tbl$alpha == 0),"model"] <- "ridge"
trn.tbl[(trn.tbl$model == "glmnet") & (trn.tbl$alpha == 1),"model"] <- "lasso"
trn.tbl <- subset(trn.tbl, model != "glmnet")
trn.tbl$train.drug <- as.character(trn.tbl$train.drug)
trn.tbl$model <- as.character(trn.tbl$model)
trn.tbl$coeffs <- as.character(trn.tbl$coeffs)

trn.tbl.lasso <- subset(trn.tbl, model == "lasso")
trn.tbl.ridge <- subset(trn.tbl, model == "ridge")

save.image(".Rdata")

## Looking at correlation of coefficients
cat("Summarizing ridge results\n")
all.ridge <- c()
for(gs in unique(trn.tbl.ridge$gene.set)) {
  tmp <- subset(trn.tbl.ridge, gene.set == gs)
  tr.sets <- as.character(unique(tmp$train.set))
  for(i in 1:(length(tr.sets)-1)) {
    tr1 <- tr.sets[i]
    tmp1 <- subset(tmp, train.set == tr1)
    merge.col <- "DRUG_ID"
    if(grepl(pattern = "ohsu", tr1)) { merge.col <- "inhibitor" }
    tmp1 <- merge(drug.name.tbl, tmp1, by.x = merge.col, by.y = "train.drug")
    for(j in (i+1):length(tr.sets)) {
      tr2 <- tr.sets[j]
      if(((tr1 == "ohsu") || (tr2 == "ohsu")) && (grepl(pattern="ohsu.set1", paste0(tr1,tr2)) || grepl(pattern="ohsu.set2", paste0(tr1,tr2)))) { next }
      if(((tr1 == "fimm") || (tr2 == "fimm")) && (grepl(pattern="fimm.set1", paste0(tr1,tr2)) || grepl(pattern="fimm.set2", paste0(tr1,tr2)))) { next }
      tmp2 <- subset(tmp, train.set == tr2)
      merge.col <- "DRUG_ID"
      if(grepl(pattern = "ohsu", tr2)) { merge.col <- "inhibitor" }
      tmp2 <- merge(drug.name.tbl, tmp2, by.x = merge.col, by.y = "train.drug")

      m <- merge(tmp1, tmp2, by = "inhibitor", suffixes = c(".1", ".2"))
      m <- m[, c("inhibitor", "n.features.1", "coeffs.1", "coeffs.2", "coeff.vals.1", "coeff.vals.2", "train.set.1", "train.set.2", "gene.set.1")]
      m$id <- 1:nrow(m)
      hyp <- ddply(m, .variables = "id",
             .fun = function(r) {
                      coeff.vals.1 <- as.numeric(unlist(strsplit(as.character(r$coeff.vals.1), split=",")))
                      coeffs.1 <- (unlist(strsplit(as.character(r$coeffs.1), split=",")))
                      names(coeff.vals.1) <- coeffs.1
                      coeffs.1 <- coeff.vals.1
                      coeff.vals.2 <- as.numeric(unlist(strsplit(as.character(r$coeff.vals.2), split=",")))
                      coeffs.2 <- (unlist(strsplit(as.character(r$coeffs.2), split=",")))
                      names(coeff.vals.2) <- coeffs.2
                      coeffs.2 <- coeff.vals.2
                      inter <- intersect(names(coeffs.1), names(coeffs.2))
                      ct <- cor.test(coeffs.1[inter], coeffs.2[inter])
                      pval <- ct$p.value
                      cor <- ct$estimate
                      vec <- c(r$id[1], cor, pval)
                      names(vec) <- c("id", "cor", "pval")
		      vec
                    })
      hyp <- merge(hyp, m, by = "id")
      all.ridge <- rbind(all.ridge, hyp)
    }
  }
}
all.ridge$pval <- as.numeric(all.ridge$pval)
cat("Done summarizing ridge results\n")

all.ridge.sig <- subset(all.ridge, pval < 0.05 & cor > 0)
cross.ridge.sig <- subset(all.ridge.sig, train.set.1 == "ohsu" & train.set.2 == "fimm")
ohsu.ridge.sig <- subset(all.ridge.sig, (train.set.1 == "ohsu.set1" & train.set.2 == "ohsu.set2") | (train.set.2 == "ohsu.set1" & train.set.1 == "ohsu.set2"))
fimm.ridge.sig <- subset(all.ridge.sig, (train.set.1 == "fimm.set1" & train.set.2 == "fimm.set2") | (train.set.2 == "fimm.set1" & train.set.1 == "fimm.set2"))

tbl <- as.data.frame(table(cross.ridge.sig$gene.set.1))
colnames(tbl) <- c("Feature Set", "Drugs w/ Correlated Features")
pdf("ohsu-vs-fimm-ridge-feature-overlap.pdf")
plot.table(tbl, main = "OHSU (train) vs FIMM (test)")
d <- dev.off()

tbl <- as.data.frame(table(ohsu.ridge.sig$gene.set.1))
colnames(tbl) <- c("Feature Set", "Drugs w/ Correlated Features")
pdf("ohsu1-vs-ohsu2-ridge-feature-overlap.pdf")
plot.table(tbl, main = "OHSU 1 (train) vs OHSU 2 (test)")
d <- dev.off()

tbl <- as.data.frame(table(fimm.ridge.sig$gene.set.1))
colnames(tbl) <- c("Feature Set", "Drugs w/ Correlated Features")
pdf("fimm1-vs-fimm2-ridge-feature-overlap.pdf")
plot.table(tbl, main = "FIMM 1 (train) vs FIMM 2 (test)")
d <- dev.off()


save.image(".Rdata")

all <- c()
cat("Summarizing lasso results\n")
for(gs in unique(trn.tbl.lasso$gene.set)) {
  tmp <- subset(trn.tbl.lasso, gene.set == gs)
  tr.sets <- as.character(unique(tmp$train.set))
  for(i in 1:(length(tr.sets)-1)) {
    tr1 <- tr.sets[i]
    tmp1 <- subset(tmp, train.set == tr1)
    merge.col <- "DRUG_ID"
    if(grepl(pattern = "ohsu", tr1)) { merge.col <- "inhibitor" }
    tmp1 <- merge(drug.name.tbl, tmp1, by.x = merge.col, by.y = "train.drug")
    for(j in (i+1):length(tr.sets)) {
      tr2 <- tr.sets[j]
      if(((tr1 == "ohsu") || (tr2 == "ohsu")) && (grepl(pattern="ohsu.set1", paste0(tr1,tr2)) || grepl(pattern="ohsu.set2", paste0(tr1,tr2)))) { next }
      if(((tr1 == "fimm") || (tr2 == "fimm")) && (grepl(pattern="fimm.set1", paste0(tr1,tr2)) || grepl(pattern="fimm.set2", paste0(tr1,tr2)))) { next }
      tmp2 <- subset(tmp, train.set == tr2)
      merge.col <- "DRUG_ID"
      if(grepl(pattern = "ohsu", tr2)) { merge.col <- "inhibitor" }
      tmp2 <- merge(drug.name.tbl, tmp2, by.x = merge.col, by.y = "train.drug")

      m <- merge(tmp1, tmp2, by = "inhibitor", suffixes = c(".1", ".2"))
      m <- m[, c("inhibitor", "n.features.1", "coeffs.1", "coeffs.2", "train.set.1", "train.set.2", "gene.set.1")]
      m$id <- 1:nrow(m)
      hyp <- ddply(m, .variables = "id",
             .fun = function(r) {
                      n.f <- as.numeric(r$n.features.1)
                      coeff.vals.1 <- as.numeric(unlist(strsplit(as.character(r$coeff.vals.1), split=",")))
                      coeffs.1 <- (unlist(strsplit(as.character(r$coeffs.1), split=",")))
                      names(coeff.vals.1) <- coeffs.1
                      coeffs.1 <- coeff.vals.1
                      coeff.vals.2 <- as.numeric(unlist(strsplit(as.character(r$coeff.vals.2), split=",")))
                      coeffs.2 <- (unlist(strsplit(as.character(r$coeffs.2), split=",")))
                      names(coeff.vals.2) <- coeffs.2
                      coeffs.2 <- coeff.vals.2
                      coeffs.1 <- coeffs.1[coeffs.1 > 0]
                      coeffs.2 <- coeffs.2[coeffs.2 > 0]
                      n.1 <- length(unlist(strsplit(coeffs.1, split=",")))
                      n.2 <- length(unlist(strsplit(coeffs.2, split=",")))
                      both <- intersect(unlist(strsplit(coeffs.1, split=",")), unlist(strsplit(coeffs.2, split=",")))
                      n.both <- length(both)
                      both <- paste(both, collapse=",")
                      pval <- phyper(q = n.both - 1, m = n.1, n = n.f - n.1, k = n.2, lower.tail = FALSE) 
                      vec <- c(r$id[1], n.1, n.2, n.both, both, pval)
                      names(vec) <- c("id", "n.1", "n.2", "n.both", "both", "pval")
		      vec
                    })
      hyp <- merge(hyp, m, by = "id")
      all <- rbind(all, hyp)
    }
  }
}
all$pval <- as.numeric(all$pval)
cat("Done summarizing lasso results\n")

save.image(".Rdata")

all.sig <- subset(all, pval < 0.05)
cross.sig <- subset(all.sig, train.set.1 == "ohsu" & train.set.2 == "fimm")
ohsu.sig <- subset(all.sig, (train.set.1 == "ohsu.set1" & train.set.2 == "ohsu.set2") | (train.set.2 == "ohsu.set1" & train.set.1 == "ohsu.set2"))
fimm.sig <- subset(all.sig, (train.set.1 == "fimm.set1" & train.set.2 == "fimm.set2") | (train.set.2 == "fimm.set1" & train.set.1 == "fimm.set2"))
sig.overlap <- ddply(all.sig, .variables = c("gene.set.1", "train.set.1", "train.set.2"),
      .fun = function(df) { 
               num <- length(which(df$pval < 0.05))
               df <- unique(df[, c("gene.set.1", "train.set.1", "train.set.2")])
               df$num <- num
               df
             })
table(cross.sig$gene.set.1)
table(ohsu.sig$gene.set.1)
table(fimm.sig$gene.set.1)

tbl <- as.data.frame(table(cross.sig$gene.set.1))
colnames(tbl) <- c("Feature Set", "Tests w/ Overlapping Features")
pdf("ohsu-vs-fimm-feature-overlap.pdf")
plot.table(tbl, main = "OHSU (train) vs FIMM (test)")
d <- dev.off()

tbl <- as.data.frame(table(ohsu.sig$gene.set.1))
colnames(tbl) <- c("Feature Set", "Tests w/ Overlapping Features")
pdf("ohsu1-vs-ohsu2-feature-overlap.pdf")
plot.table(tbl, main = "OHSU 1 (train) vs OHSU 2 (test)")
d <- dev.off()

tbl <- as.data.frame(table(fimm.sig$gene.set.1))
colnames(tbl) <- c("Feature Set", "Tests w/ Overlapping Features")
pdf("fimm1-vs-fimm2-feature-overlap.pdf")
plot.table(tbl, main = "FIMM 1 (train) vs FIMM 2 (test)")
d <- dev.off()



response.col <- "dss.auc.ll4"
training.res <- list()
if("ohsu" %in% data.sets.to.analyze) {
  cat("Train on OHSU\n")
  ret <- train.model(train.dss.arg = dss.orig[["ohsu"]], train.expr.arg = expr.filt.matrices[["ohsu"]][common.genes, ], 
                     train.genomic.arg = NULL, train.clinical.arg = NULL, train.common.drugs = common.drugs.by.data.set[["ohsu"]], 
                     train.drugs = common.drugs.by.data.set[["ohsu"]], train.drug.col = ohsu.drug.col,
                     train.patient.col = ohsu.patient.col, train.response.col = response.col, seed = 1234, use.rf = FALSE, use.svm = FALSE)
  training.res[["ohsu"]] <- ret
}

if("ohsu" %in% data.sets.to.analyze){ 
  cat("Train on fimm\n")
  ret <- train.model(train.dss.arg = dss.orig[["fimm"]], train.expr.arg = expr.filt.matrices[["fimm"]][common.genes, ], 
                     train.genomic.arg = NULL, train.clinical.arg = NULL, train.common.drugs = common.drugs.by.data.set[["fimm"]], 
                     train.drugs = common.drugs.by.data.set[["fimm"]], train.drug.col = fimm.drug.col,
                     train.patient.col = fimm.patient.col, train.response.col = response.col, seed = 1234, use.rf = FALSE, use.svm = FALSE)
  training.res[["fimm"]] <- ret
}

## Create a table of drug by feature for lasso
for(indx1 in 1:(length(names(training.res))-1)) {
  train.set1 <- names(training.res)[indx1]
  merge.col <- NULL
  switch(train.set1,
         "fimm" = { merge.col <- "DRUG_ID" },
         "ohsu" = { merge.col <- "inhibitor" },
         { die(paste0("Unknown train set ", train.set1, "\n")) })
  tbl1 <- ldply(training.res[[train.set1]], 
                .fun = function(lst) {
                  for(i in 1:length(lst)) {
                    if((lst[[i]]$model == "glmnet") && (lst[[i]]$alpha == 1)) {
                      coeffs <- coefficients(lst[[i]]$fit)
                      ns <- rownames(coeffs)
                      coeffs <- as.vector(coeffs)
                      n.features <- length(coeffs) - 1 ## don't count intercept
                      names(coeffs) <- ns
                      coeffs <- coeffs[coeffs > 0]
                      coeffs <- coeffs[!grepl(pattern="Intercept", names(coeffs))]
                      vec <- c(lst[[i]]$train.drug, paste(names(coeffs), collapse=","), n.features)
                      names(vec) <- c(merge.col, "features", "n.features")
                      return(vec)
                    }
                  }
                })
  tbl1$train.set <- train.set1
  tbl1 <- merge(tbl1, drug.name.tbl, by = merge.col)
  for(indx2 in (indx1+1):length(names(training.res))) {
    train.set2 <- names(training.res)[indx2]
    merge.col <- NULL
    switch(train.set2,
           "fimm" = { merge.col <- "DRUG_ID" },
           "ohsu" = { merge.col <- "inhibitor" },
           { die(paste0("Unknown train set ", train.set2, "\n")) })
    tbl2 <- ldply(training.res[[train.set2]], 
                  .fun = function(lst) {
                    for(i in 1:length(lst)) {
                      if((lst[[i]]$model == "glmnet") && (lst[[i]]$alpha == 1)) {
                        coeffs <- coefficients(lst[[i]]$fit)
                        ns <- rownames(coeffs)
                        coeffs <- as.vector(coeffs)
                        n.features <- length(coeffs) - 1 ## don't count intercept
                        names(coeffs) <- ns
                        coeffs <- coeffs[coeffs > 0]
                        coeffs <- coeffs[!grepl(pattern="Intercept", names(coeffs))]
                        vec <- c(lst[[i]]$train.drug, paste(names(coeffs), collapse=","), n.features)
                        names(vec) <- c(merge.col, "features")
                        return(vec)
                      }
                    }
                  })
    tbl2$train.set <- train.set2
    tbl2 <- merge(tbl2, drug.name.tbl, by = merge.col)
    mtbl <- merge(tbl1, tbl2, by = "inhibitor", suffixes = c(paste0(".", train.set1), paste0(".", train.set2)))
  }
}

mtbl$id <- 1:nrow(mtbl)
hyp <- ddply(mtbl, .variables = "id",
             .fun = function(r) {
                      n.f <- as.numeric(r$n.features)
                      n.fimm <- length(unlist(strsplit(r$features.fimm, split=",")))
                      n.ohsu <- length(unlist(strsplit(r$features.ohsu, split=",")))
                      n.both <- length(intersect(unlist(strsplit(r$features.fimm, split=",")), unlist(strsplit(r$features.ohsu, split=","))))
                      pval <- phyper(q = n.both - 1, m = n.fimm, n = n.f - n.fimm, k = n.ohsu, lower.tail = FALSE) 
                      vec <- c(r$id[1], n.fimm, n.ohsu, n.both, pval)
                      names(vec) <- c("id", "n.fimm", "n.ohsu", "n.both", "pval")
		      vec
                    })
hyp <- merge(hyp, mtbl, by = "id")

mtbl$pval <- unlist(apply(mtbl[, c("n.features", "features.fimm", "features.ohsu")], 1, 
                    function(row) {
                      n.f <- as.numeric(row[1])
                      n.fimm <- length(unlist(strsplit(row[2], split=",")))
                      n.ohsu <- length(unlist(strsplit(row[3], split=",")))
                      n.both <- length(intersect(unlist(strsplit(row[2], split=",")), unlist(strsplit(row[3], split=","))))
print(n.both)
                      phyper(q = n.both - 1, m = n.fimm, n = n.f - n.fimm, k = n.ohsu, lower.tail = FALSE) 
                    }))

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

cat("Training and testing with DSS2\n")
train.response.cols <- rep("dss2.ll4", length(train.response.cols))
test.response.cols <- rep("dss2.ll4", length(test.response.cols))
res.dss2 <- train.and.test.crossproduct(gene.sets = gene.sets, drug.name.tbl = drug.name.tbl, train.set.names = train.set.names, train.dss.args = train.dss.args, 
                                        train.expr.args = train.expr.args, train.genomic.args = train.genomic.args, train.clinical.args = train.clinical.args, 
                                        train.common.drugs = train.common.drugs, train.drugs = train.drugs, train.drug.cols = train.drug.cols, train.patient.cols = train.patient.cols, 
                                        train.response.cols = train.response.cols, test.set.names = test.set.names, test.dss.args = test.dss.args, test.expr.args = test.expr.args, 
                                        test.genomic.args = test.genomic.args, test.clinical.args = test.clinical.args, test.common.drugs = test.common.drugs,
                                        test.drug.cols = test.drug.cols, test.patient.cols = test.patient.cols, test.response.cols = test.response.cols, 
                                        seed = 1234, use.rf = FALSE, use.svm = FALSE, num.processes = num.processes) 
cat("Done training and testing with DSS2\n")

save.image(".Rdata")

tbl <- ldply(train.set.names,
             .fun = function(train.set) {
                      ldply(test.set.names, 
                            .fun = function(test.set) {
                                     ldply(names(gene.sets),
                                           .fun = function(gene.set) {
                                                    tmp.dss2 <- extract.results(res.dss2[["all.comparisons"]][[train.set]][[test.set]][[gene.set]])
                                                    tmp.dss2$train.set <- train.set
                                                    tmp.dss2$test.set <- test.set
                                                    tmp.dss2$gene.set <- gene.set
                                                    tmp.dss2$response <- "DSS2"
                                                    tmp.auc <- extract.results(res.auc[["all.comparisons"]][[train.set]][[test.set]][[gene.set]])
                                                    tmp.auc$train.set <- train.set
                                                    tmp.auc$test.set <- test.set
                                                    tmp.auc$gene.set <- gene.set
                                                    tmp.auc$response <- "AUC"
                                                    rbind(tmp.auc, tmp.dss2)                         
                                                  }) 
                                   }) 
                    })

save.image(".Rdata")

## For each train-test combination create a vector of correlations (one for each drug) of predicted vs actual
corr.lists <- list()
for(test.set.name in unique(tbl$test.set)) {
  corr.lists[[test.set.name]] <- list()
  for(mdl in unique(tbl$model)) {
    corr.lists[[test.set.name]][[mdl]] <- list()
    for(resp in unique(tbl$response)) {
      corr.lists[[test.set.name]][[mdl]][[resp]] <- list()
      for(gs in unique(tbl$gene.set)) {
        corr.lists[[test.set.name]][[mdl]][[resp]][[gs]] <- list()
        t.sub <- subset(tbl, (test.set == test.set.name) & (model == mdl) & (response == resp) & (gene.set == gs))
        na.drugs <- t.sub$test.drug[is.na(t.sub$val)]
        t.sub <- t.sub[!(t.sub$test.drug %in% na.drugs),]
        test.drugs <- unique(t.sub$test.drug)
        for(train.set.name in unique(tbl$train.set)) {
          tr.sub <- subset(t.sub, train.set == train.set.name) 
          vec <- tr.sub$val
          names(vec) <- tr.sub$test.drug
          corr.lists[[test.set.name]][[mdl]][[resp]][[gs]][[train.set.name]] <- vec[test.drugs]
        }
      }
    }
  }
}

capwords <- function(s, strict = FALSE) {
    cap <- function(s) paste(toupper(substring(s, 1, 1)),
                  {s <- substring(s, 2); if(strict) tolower(s) else s},
                             sep = "", collapse = " " )
    sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}

## Make a plot of correlation vs gene, biocart, hallmark, kegg, and reactome 
tbl$gene.set <- factor(tbl$gene.set, levels = c("gene", "biocarta", "hallmark", "kegg", "reactome"))
tbl$train.set <- factor(tbl$train.set, levels = c("ctrp.all", "ctrp.heme", "ctrp.aml", "ohsu.train"))
for(resp in c("DSS2", "AUC")) {
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
      file <- paste0(resp, "-", mdl, "-", met, "-vs-gene-sets-fo.pdf")
      pdf(file)
      print(g)
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

if(length(unique(tbl$train.set)) > 1) { 
## Create pairwise correlation plots
for(test.set.name in unique(tbl$test.set)) {
  for(mdl in c("lasso", "ridge")) {
    for(resp in unique(tbl$response)) {
      pdf(paste(test.set.name, mdl, resp, "corr-fo.pdf", sep="-"))
      for(gs in unique(tbl$gene.set)) {
        mat <- do.call("cbind", corr.lists[[test.set.name]][[mdl]][[resp]][[gs]])
        cols <- NULL
        if("Class" %in% colnames(drug.name.tbl)) {
          cls <- drug.name.tbl$Class
          switch(test.set.name,
            "fimm" = { names(cls) <- drug.name.tbl$DRUG_ID },
            "ohsu.test" = { names(cls) <- drug.name.tbl$inhibitor },
             { die(paste0("Unknown test set ", test.set.name, "\n")) }) 
          cols <- factor(cls[names(corr.lists[[test.set.name]][[mdl]][[resp]][[gs]][[1]])])
        }
##        pdf(paste(test.set.name, mdl, resp, gs, "corr.pdf", sep="-"))
        pairs(mat, lower.panel=panel.regression, upper.panel=panel.cor, main = paste(test.set.name, mdl, resp, gs, sep=" "), col = cols, oma=c(4,4,6,12))
        par(xpd=TRUE)
        legend(0.85, 0.7, legend = levels(cols), fill = c(1:length(levels(cols))), cex = 0.7)
##        d <- dev.off()
      }
      d <- dev.off()
    }
  }
}

for(test.set.name in unique(tbl$test.set)) {
  for(mdl in c("lasso", "ridge")) {
    for(resp in unique(tbl$response)) {
      pdf(paste(test.set.name, mdl, resp, "diff-fo.pdf", sep="-"))
      for(gs in unique(tbl$gene.set)) {
        mat <- do.call("cbind", corr.lists[[test.set.name]][[mdl]][[resp]][[gs]])
        cols <- NULL
        if("Class" %in% colnames(drug.name.tbl)) {
          cls <- drug.name.tbl$Class
          switch(test.set.name,
            "fimm" = { names(cls) <- drug.name.tbl$DRUG_ID },
            "ohsu.test" = { names(cls) <- drug.name.tbl$inhibitor },
             { die(paste0("Unknown test set ", test.set.name, "\n")) }) 
          cols <- factor(cls[names(corr.lists[[test.set.name]][[mdl]][[resp]][[gs]][[1]])])
        }
        len <- length(corr.lists[[test.set.name]][[mdl]][[resp]][[gs]])
        glist <- list()
        for(i in 1:len) {
          for(j in 1:len) {
            x <- corr.lists[[test.set.name]][[mdl]][[resp]][[gs]][[i]]
            y <- corr.lists[[test.set.name]][[mdl]][[resp]][[gs]][[j]]
            g <- NULL
            if(i != j) {
              if(i > j) {
                df <- data.frame(x = x, y = y, col = cols)
              } else {
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
                ## g <- g + geom_text(x = lims[1] + 0.2, y = lims[1] + 0.1, label = lm_eqn(df), parse=TRUE)
              }
            } else {
              g <- ggplot() + annotation_custom(textGrob(names(corr.lists[[test.set.name]][[mdl]][[resp]][[gs]])[i]),
                               xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) 
            }
            g <- g + xlab("") + ylab("")
            glist[[length(glist)+1]] <- g
          }
        }
        do.call("grid.arrange", c(glist, "top" = paste(test.set.name, mdl, resp, gs, sep=" ")))
      }
      d <- dev.off()
    }
  }
}

} ## if(length(unique(tbl$train.set)) > 1)

stop("successfully completed")

test.set.name <- "fimm"
mdl <- "ridge"
resp <- "AUC"
gs <- "gene"
x <- corr.lists[[test.set.name]][[mdl]][[resp]][[gs]][["ctrp.all"]]
y <- corr.lists[[test.set.name]][[mdl]][[resp]][[gs]][["ctrp.aml"]]
df.fimm <- data.frame(name = names(x), fimm = y - x, fimm.ctrp.aml = y, fimm.ctrp.all = x)
df.fimm <- merge(drug.name.tbl, df.fimm, by.y = "name", by.x = "DRUG_ID")

test.set.name <- "ohsu.test"
x <- corr.lists[[test.set.name]][[mdl]][[resp]][[gs]][["ctrp.all"]]
y <- corr.lists[[test.set.name]][[mdl]][[resp]][[gs]][["ctrp.aml"]]
df.ohsu <- data.frame(name = names(x), ohsu = y - x, ohsu.ctrp.aml = y, ohsu.ctrp.all = x)
df <- merge(df.ohsu, df.fimm, by.x = "name", by.y = "inhibitor")

g <- ggplot(data = df)
g <- g + geom_point(aes(x = fimm, y = ohsu))
g <- g + geom_hline(yintercept = 0)
g <- g + geom_vline(xintercept = 0)
g <- g + xlab("CTRP AML - CTRP ALL (FIMM Test)")
g <- g + ylab("CTRP AML - CTRP ALL (OHSU Test)")
print(g)

## Drugs with predictions positively correlated with actual response in AML and with greater correlations in AML than in ALL
subset(df[, c("name", "Mechanism", "Class", "fimm", "ohsu", "fimm.ctrp.aml", "fimm.ctrp.all", "ohsu.ctrp.aml", "ohsu.ctrp.all")], 
       (fimm > 0) & (ohsu > 0) & (fimm.ctrp.aml > 0) & (ohsu.ctrp.aml > 0))

## Drugs with predictions positively correlated with actual response in ALL and with greater correlations in ALL than in AML
subset(df[, c("name", "Mechanism", "Class", "fimm", "ohsu", "fimm.ctrp.aml", "fimm.ctrp.all", "ohsu.ctrp.aml", "ohsu.ctrp.all")], 
       (fimm < 0) & (ohsu < 0) & (fimm.ctrp.all > 0) & (ohsu.ctrp.all > 0))


num.fimm.samples <- length(unique(intersect(fimm.dss.orig$SCREEN_ID, colnames(fimm.expr))))

fisher.z.transform <- TRUE
for(met in c("spearman", "mse", "pearson")) {
  all.results <- rbind(td, te, tf, tg, th)
  ## all.results <- subset(all.results, metric == "spearman")
  all.results <- subset(all.results, metric == met)
  all.results$train.set <- factor(all.results$train.set, levels = c(ctrp.all.label, ctrp.heme.label, ctrp.aml.label, ohsu.label, ohsu.label.down))
  tmp <- all.results
  cat(paste0(met, "\n"))
  if(met != "mse") {
    ## Mean model can not be evaluated by correlation
    tmp <- tmp[tmp$model != "mean",]
  }
  if(fisher.z.transform && (met != "mse")) {
    cat("Fisher z transforming\n")
    tmp$val <- 0.5 * log(1 + tmp$val) - 0.5 * log(1 - tmp$val)
    tmp$val <- tmp$val * sqrt(tmp$n.test - 3)
  }
  g <- ggplot(data = tmp, aes(x = train.set, y = val))
  g <- g + geom_violin()
  g <- g + geom_boxplot(width = 0.1)
  g <- g + facet_wrap(~ model)
  g <- g + xlab(paste0("Training Data Set [Test Data Set = FIMM (n=", num.fimm.samples, ")]"))
  g <- g + ylab(met)
  g <- g + theme(axis.text.x=element_text(angle = 45, hjust = 1))
  pdf(paste0("ctrp-ohsu-fimm-boxplot-", met, "-fo.pdf"))
  print(g)
  d <- dev.off()

  g <- ggplot(data = tmp, aes(x = train.set, y = val))
  g <- g + geom_violin()
  g <- g + geom_beeswarm()
  g <- g + facet_wrap(~ model)
  g <- g + xlab(paste0("Training Data Set [Test Data Set = FIMM (n=", num.fimm.samples, ")]"))
  g <- g + ylab(met)
  g <- g + theme(axis.text.x=element_text(angle = 45, hjust = 1))
  pdf(paste0("ctrp-ohsu-fimm-beeswarm-", met, "-fo.pdf"))
  print(g)
  d <- dev.off()
}

for(mdl in c("lasso", "ridge", "mean")) {
  for(met in c("spearman", "pearson")) {
    all.results <- rbind(td, te, tf, tg, th, ti)
    sub <- subset(all.results, (model == mdl) & (metric == met))
    sub <- sub[order(sub$val, decreasing=TRUE),]
    sub$drug.name <- sub$test.drug
    sub$drug.name <- factor(sub$drug.name, levels = sub[sub$cmp == fimm.label, "drug.name"])
    glist <- dlply(sub, .variables = "cmp",
                        .fun = function(df) {
                          g <- ggplot(df, aes(x = drug.name, y = val))
                          g <- g + geom_point()
                          g <- g + ylab(met)
                          g <- g + ggtitle(paste0(mdl, ": ", df$cmp[1]))
                          g
                        })
    do.call("grid.arrange", glist)
  }
}

for(mdl in c("lasso", "ridge", "mean")) {
  for(met in c("spearman", "pearson")) {
    datasets <- c("ctrp.all", "ctrp.heme", "ctrp.aml", "ohsu", "fimm")
    glist <- llply(datasets[datasets != "fimm"],
                   .fun = function(dataset) {
                     train.vs.train <- paste0(dataset, ".vs.", dataset)
                     train.err <- extract.results(rets[[train.vs.train]])
                     train.err <- subset(train.err, (model == mdl) & (metric == met))
                     train.vs.test <- paste0(dataset, ".vs.", "fimm")
                     val.err <- extract.results(rets[[train.vs.test]])
                     val.err <- subset(val.err, (model == mdl) & (metric == met))
                     tbl.tmp <- merge(train.err, val.err, by = c("train.drug"), suffixes = c(".train", ".test"))
                     df <- tbl.tmp[, c("val.train", "val.test")]
                     colnames(df) <- c("x", "y")
                     g <- ggplot(df, aes(x = x, y = y))
                     g <- g + geom_point()
                     g <- g + ylab(train.vs.test)
                     g <- g + xlab(train.vs.train)
                     g <- g + ggtitle(paste0(mdl, " ", met, ": ", train.vs.test))
                     g <- g + geom_smooth(method='lm')
                     lims <- c(min(min(df$x, na.rm=TRUE), min(df$y, na.rm=TRUE)), max(max(df$x, na.rm=TRUE), max(df$y, na.rm=TRUE)))
                     ## g <- g + xlim(c(-1,1))
                     ## g <- g + ylim(c(-1,1))
                     g <- g + xlim(lims)
                     g <- g + ylim(lims)
                     ## g <- g + geom_text(x = min(df$x, na.rm=TRUE) + 0.2, y = min(df$y, na.rm=TRUE) + 0.1, label = lm_eqn(df), parse=TRUE)
                     ## g <- g + geom_text(x = -0.8, y = -0.9, label = lm_eqn(df), parse=TRUE)
                     g <- g + geom_text(x = lims[1] + 0.2, y = lims[1] + 0.1, label = lm_eqn(df), parse=TRUE)
                     g
                   })
    pdf(paste0(mdl, "-", met, "-train-vs-test-correlations-fo.pdf"))
    do.call("grid.arrange", glist)
    d <- dev.off()
  }
}

all.results <- rbind(td, te, tf, tg)
## all.results <- subset(all.results, metric == "spearman")
all.results <- subset(all.results, metric == "mse")
all.results$train.set <- factor(all.results$train.set, levels = c(ctrp.all.label, ctrp.heme.label, ctrp.aml.label, ohsu.label))
g <- ggplot(data = all.results, aes(x = train.set, y = val))
g <- g + geom_violin()
g <- g + geom_beeswarm()
g <- g + facet_wrap(~ model)
g <- g + xlab(paste0("Training Data Set [Test Data Set = FIMM (n=", num.fimm.samples, ")]"))
g <- g + ylab("MSE")
g <- g + theme(axis.text.x=element_text(angle = 45, hjust = 1))
pdf("ctrp-ohsu-fimm-beeswarm-mse-fo.pdf")
print(g)
d <- dev.off()


## Look at consistency of features across training sets in lasso fits

## lasso is index 2 in our results above
model.indx <- 2

all.training.sets <- c("ohsu", "fimm", "ctrp.all", "ctrp.heme", "ctrp.aml")

## Column in drug.name.tbl that holds the drug identifier for each of the above training sets
drug.id.col <- c("inhibitor", "DRUG_ID", "master_cpd_id", "master_cpd_id", "master_cpd_id")

all.training.sets <- c("ohsu", "fimm", "ctrp.aml")
drug.id.col <- c("inhibitor", "DRUG_ID", "master_cpd_id")

all.training.sets <- c("ohsu", "fimm")
drug.id.col <- c("inhibitor", "DRUG_ID")

drug.gene.matrices <- 
  llply(1:length(all.training.sets),
        .fun = function(i) {
                 train.set <- all.training.sets[i]
                 df <- ldply(training.res[[train.set]], 
                             .fun = function(l) {
                                      lasso.fit <- l[[model.indx]]$fit
                                      drug <- l[[model.indx]]$train.drug
                                      if(is.na(lasso.fit)) { return(NULL) }
                                      coef <- as.matrix(coefficients(lasso.fit))
                                      coef <- coef[coef[,1] > 0,,drop=F]
                                      if(nrow(coef) == 0) { return(NULL) }
                                      features <- rownames(coef)
                                      return(data.frame(drug = rep(drug, length(features)), feature = features))
                                    })
                 df <- df[!grepl(df$feature, pattern="Intercept"),]
                 df <- merge(df, drug.name.tbl, by.x = "drug", by.y = drug.id.col[i])
                 df <- df[, c("DRUG_NAME", "feature")]
                 df <- unique(df)
                 ## Convert ENSG id to Hugo symbol
                 mp <- ensg.to.sym.mapping(unique(df$feature))
                 df <- merge(df, mp, by.x = "feature", by.y = "ensg")
                 df <- df[, c("DRUG_NAME", "symbol")]
                 colnames(df) <- c("drug", "feature")
                 df
                 table(df$feature, df$drug)
               })

rows <- unique(unlist(lapply(drug.gene.matrices, rownames)))
cols <- unique(unlist(lapply(drug.gene.matrices, colnames)))


hm <- matrix(data = 0, nrow = length(rows), ncol = length(cols), dimnames = list(rows, cols))
for(i in 1:length(drug.gene.matrices)) {
  mat <- drug.gene.matrices[[i]]
  rs <- rownames(mat)
  cs <- colnames(mat)
  hm[rs, cs] <- hm[rs, cs] + mat[rs, cs]
}


suppressPackageStartupMessages(library(heatmap.2))

pdf("fimm-ohsu-gene-drug-lasso-matrix.pdf")
my_palette <- colorRampPalette(c("white", "red", "black"))(n = 5)
heatmap.2(hm, trace = "none", col = my_palette)
d <- dev.off()

drug.cnts <- apply(hm, 2, function(v) sum(v[v > 0]))
gene.cnts <- apply(hm, 1, function(v) sum(v[v > 0]))

data <- data.frame(x = names(drug.cnts), y = drug.cnts)
data$x <- factor(data$x)
g <- ggplot(data = data, aes(x = x, y = y))
g <- g + geom_point()
g <- g + xlab("Drug")
g <- g + ylab("Num features")
g <- g + theme(axis.text.x = element_text(angle = 90))
print(g)

data <- data.frame(y = names(gene.cnts), x = gene.cnts)
data <- subset(data, x > 2)
data <- data[order(data$x, decreasing = TRUE),]
data$y <- factor(data$y)
g <- ggplot(data = data, aes(x = x, y = y))
g <- g + geom_point()
g <- g + ylab("Gene")
g <- g + xlab("Occurrences in Drug Model")
g <- g + theme(axis.text.y = element_text(size = 10))
## g <- g + theme(axis.text.x = element_text(angle = 90))
pdf("fimm-ohsu-gene-occurrence.pdf")
print(g)
d <- dev.off()