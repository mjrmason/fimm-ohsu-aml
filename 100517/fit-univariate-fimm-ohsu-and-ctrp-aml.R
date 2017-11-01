
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

cat("Do gene vs drug univariate fits for OHSU, FIMM, and CTRP AML\n")

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
data.sets.to.analyze <- c("ohsu", "fimm", "ctrp")

## Update as change data sets
processed.file.names <- c("ohsu.foc.common.drugs.dss.t0.tsv", "fimm.foc.common.drugs.dss.t0.tsv", "ctrp.foc.common.drugs.dss.t0.tsv")
names(processed.file.names) <- c("ohsu", "fimm", "ctrp")
## parentId is the processed data folder of the respective data sets
parentIds <- c("syn10083332", "syn8270577", "syn10288724")
names(parentIds) <- c("ohsu", "fimm", "ctrp")

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
## paths <- c("ohsu-expr-foc-aml-combat.tsv", "fimm-expr-foc-aml-combat.tsv", "ctrp-aml-expr-foc-aml-combat.tsv")
expression.synIds <- c("syn11154802", "syn11154803", "syn11154806")
names(expression.synIds) <- c("ohsu", "fimm", "ctrp")

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
## NB: the ctrp expression data is limited to AML
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

uni.fits <- list()
names(expression.synIds) <- c("ohsu", "fimm", "ctrp")
drug.cols <- c("inhibitor", "foo", "bar")
drug.pt.cols <- c("SeqID", "baz", "baf")
uni.fits <- llply(1:length(names(expression.synIds)), 
                  .parallel = TRUE,
                  .fun = function(i) {
  data.set <- names(expression.synIds)[i]
  fits <- c()
  drugs <- unique(dss.orig[[data.set]][, drug.cols[i]])
  genes <- rownames(expr.matrices[[data.set]])
  seq.pts <- colnames(expr.matrices[[data.set]])
  for(drug in drugs) {
    flag <- dss.orig[[data.set]][, drug.cols[i]] == drug
    drug.vals <- dss.orig[[data.set]][flag, "dss.auc.ll4"]
    names(drug.vals) <- dss.orig[[data.set]][flag, drug.pt.cols[i]]
    pts <- intersect(names(drug.vals), seq.pts)
    for(gene in genes) {
      print(c(data.set, drug, gene))
      gene.expr <- expr.matrices[[data.set]][gene, pts]
      drug.resp <- drug.vals[pts]
      flag <- !is.na(gene.expr) & !is.na(drug.resp)
      gene.expr <- gene.expr[flag]
      drug.resp <- drug.resp[flag]
      df <- data.frame(x = gene.expr, y = drug.resp)
      lmf <- lm(y ~ x, data = df)
      f <- summary(lmf)$fstatistic
      p <- pf(f[1],f[2],f[3],lower.tail=F)
      vec <- c(drug, gene, as.numeric(p))
      names(vec) <- c("drug", "gene", "pval")
      fits <- rbind(fits, vec)
    }
  }
  fits })

save.image(".Rdata")

stop("Successfully completed\n")


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

do.bootstrap <- FALSE
if((do.bootstrap) && ("ctrp" %in% data.sets.to.analyze)) {

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

ctrp.drug.col <- "master_cpd_id"
ctrp.patient.col <- "ccl_name"

if("ctrp" %in% data.sets.to.analyze) {

  ## The expression data for CTRP is just the AML samples--so don't do heme or all
  if(FALSE) { 
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
  } ## FALSE

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

do.dss2 <- FALSE
if(do.dss2) { 
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
}

save.image(".Rdata")

## Plot a few examples of predicted vs actual 
plot.r2 <- function(x, y) {
  df <- data.frame(x = x, y = y)
  g <- ggplot(df, aes(x = x, y = y))
  g <- g + geom_point()
  g <- g + theme(legend.position="none")
  g <- g + geom_smooth(method='lm')
  x.min <- min(df$x, na.rm=TRUE)
  x.max <- max(df$x, na.rm=TRUE)
  y.min <- min(df$y, na.rm=TRUE)
  y.max <- max(df$y, na.rm=TRUE)
  g <- g + geom_text(x = x.min + 0.5 * (x.max - x.min), y = y.min + 1 * (y.max - y.min), label = lm_eqn(df), parse=TRUE)
  ## g <- g + ggtitle(lm_eqn_bq(df))
  g
}

## Let's plot a few of the ridge results--extract all of them
ridge.res <- llply(res.auc[["all.comparisons"]][["ohsu.train"]][["fimm"]][["gene"]],
                   .fun = function(lst) {
                     for(entry in lst) {
                       if((entry$model == "glmnet") && (entry$alpha == 0)) { return(entry) }
                     }
                   })

## Calculate the correlations between predicted and actual to pick out a good one
cors <- unlist(llply(ridge.res, .fun = function(lst) as.numeric(cor.test(as.numeric(lst$predicted[,1]), lst$actual)$estimate)))
best.indx <- which(cors == max(cors, na.rm=TRUE))
indx <- best.indx
g <- plot.r2(as.numeric(ridge.res[[indx]]$predicted[,1]), ridge.res[[indx]]$actual)
g <- g + xlab("Predicted") + ylab("Actual")
g <- g + ggtitle(ridge.res[[indx]]$train.drug)
file <- paste0(ridge.res[[indx]]$train.drug, "-ridge-ohsu-vs-fimm-pred-vs-actual.pdf")
pdf(file)
print(g)
d <- dev.off()

indx <- 1
g <- plot.r2(as.numeric(ridge.res[[indx]]$predicted[,1]), ridge.res[[indx]]$actual)
g <- g + xlab("Predicted") + ylab("Actual")
g <- g + ggtitle(ridge.res[[indx]]$train.drug)
file <- paste0(ridge.res[[indx]]$train.drug, "-ridge-ohsu-vs-fimm-pred-vs-actual.pdf")
pdf(file)
print(g)
d <- dev.off()

indx <- 2
g <- plot.r2(as.numeric(ridge.res[[indx]]$predicted[,1]), ridge.res[[indx]]$actual)
g <- g + xlab("Predicted") + ylab("Actual")
g <- g + ggtitle(ridge.res[[indx]]$train.drug)
file <- paste0(ridge.res[[indx]]$train.drug, "-ridge-ohsu-vs-fimm-pred-vs-actual.pdf")
pdf(file)
print(g)
d <- dev.off()

df <- data.frame(cors = cors, x = "gene features")
g <- ggplot(data = df, aes(x = x, y = cors))
g <- g + geom_violin()
g <- g + xlab("") + ylab("Pearson's Correlation")
g <- g + geom_beeswarm()
g <- g + ggtitle("Ridge: OHSU.train vs FIMM")
file <- paste0("gene-ridge-ohsu-vs-fimm-pred-vs-actual.pdf")
pdf(file)
print(g)
d <- dev.off()

## Now aggregate all of those plots into a violin (showing points)
## See plot.correlation.vs.model for how to create a violin plot


tbl <- ldply(train.set.names,
             .fun = function(train.set) {
                      ldply(test.set.names, 
                            .fun = function(test.set) {
                                     ldply(names(gene.sets),
                                           .fun = function(gene.set) {
                                                    tmp.dss2 <- NULL
                                                    if(do.dss2) {
                                                      tmp.dss2 <- extract.results(res.dss2[["all.comparisons"]][[train.set]][[test.set]][[gene.set]])
                                                      tmp.dss2$train.set <- train.set
                                                      tmp.dss2$test.set <- test.set
                                                      tmp.dss2$gene.set <- gene.set
                                                      tmp.dss2$response <- "DSS2"
                                                    }
                                                    tmp.auc <- extract.results(res.auc[["all.comparisons"]][[train.set]][[test.set]][[gene.set]])
                                                    tmp.auc$train.set <- train.set
                                                    tmp.auc$test.set <- test.set
                                                    tmp.auc$gene.set <- gene.set
                                                    tmp.auc$response <- "AUC"
                                                    if(do.dss2) {
                                                      return(rbind(tmp.auc, tmp.dss2))
                                                    }
                                                    return(tmp.auc)                         
                                                  }) 
                                   }) 
                    })

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
      file <- paste0(resp, "-", mdl, "-", met, "-vs-gene-sets.pdf")
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
      pdf(paste(test.set.name, mdl, resp, "corr.pdf", sep="-"))
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
      pdf(paste(test.set.name, mdl, resp, "diff.pdf", sep="-"))
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
            g <- NULL
            if(i != j) {
              x <- y <- NULL
              if(i > j) {
                x <- corr.lists[[test.set.name]][[mdl]][[resp]][[gs]][[i]]
                y <- corr.lists[[test.set.name]][[mdl]][[resp]][[gs]][[j]]
                df <- data.frame(x = x, y = y, col = cols)
              } else {
                x <- corr.lists[[test.set.name]][[mdl]][[resp]][[gs]][[j]]
                y <- corr.lists[[test.set.name]][[mdl]][[resp]][[gs]][[i]]
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
              g <- ggplot() + annotation_custom(textGrob(names(corr.lists[[test.set.name]][[mdl]][[resp]][[gs]])[i]),
                               xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) 
            }
            g <- g + xlab("") + ylab("")
            if(i < j) {
              g <- g + ylab(paste0(names(corr.lists[[test.set.name]][[mdl]][[resp]][[gs]])[j], " - ", names(corr.lists[[test.set.name]][[mdl]][[resp]][[gs]])[i]))
            }
            if(i > j) {
              g <- g + xlab(names(corr.lists[[test.set.name]][[mdl]][[resp]][[gs]])[i])
              g <- g + ylab(names(corr.lists[[test.set.name]][[mdl]][[resp]][[gs]])[j])
            }

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
