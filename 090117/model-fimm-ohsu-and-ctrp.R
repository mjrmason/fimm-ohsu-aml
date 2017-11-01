
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

cat("Compare de novo vs in intro analyses\n")
cat("i.e., compare the robustness of two fits:\n")
cat("(1) fit to OHSU and validate against FIMM\n")
cat("(1) fit to CTRP and validate against FIMM\n")

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

## parentId is the OHSU processed data folder
parentId <- "syn10083332"
tbl <- synQuery(paste0("SELECT name,id FROM file WHERE parentId==\'", parentId, "\'"))

## Read in the t=0 thresholded DSS values
file.name <- "ohsu.foc.common.drugs.dss.t0.tsv"
synId <- tbl[tbl$file.name == file.name, "file.id"]
obj <- synGet(synId, downloadFile = TRUE)
file <- getFileLocation(obj)
ohsu.dss.orig <- as.data.frame(fread(file))

## parentId is the FIMM processed data folder
parentId <- "syn8270577"
tbl <- synQuery(paste0("SELECT name,id FROM file WHERE parentId==\'", parentId, "\'"))

## Read in the t=0 thresholded DSS values
file.name <- "fimm.foc.common.drugs.dss.t0.tsv"
synId <- tbl[tbl$file.name == file.name, "file.id"]
obj <- synGet(synId, downloadFile = TRUE)
file <- getFileLocation(obj)
fimm.dss.orig <- as.data.frame(fread(file))

## parentId is the CTRP processed data folder
parentId <- "syn10288724"
tbl <- synQuery(paste0("SELECT name,id FROM file WHERE parentId==\'", parentId, "\'"))

## Read in the t=0 thresholded DSS values
file.name <- "ctrp.foc.common.drugs.dss.t0.tsv"
synId <- tbl[tbl$file.name == file.name, "file.id"]
obj <- synGet(synId, downloadFile = TRUE)
file <- getFileLocation(obj)
ctrp.dss.orig <- as.data.frame(fread(file))

## Read in the expression matrices
## Read in the OHSU expression data
## path <- "ohsu.expr.tsv"
synId <- "syn10083723"
obj <- synGet(synId, downloadFile = TRUE)
file <- getFileLocation(obj)
ohsu.expr <- read.table(file, header=TRUE, sep="\t")
## Convert the column names from X20.00347 to 20-00347
colnames(ohsu.expr) <- gsub("X(.*)\\.(.*)", "\\1-\\2", colnames(ohsu.expr))

## Read in the OHSU genomic data
## NB: a 1 indicates that the gene had a non-synonymous or splice site mutation.
## path <- "ohsu.genomic.ns.ss.tsv"
## Columns correspond to samples.
synId <- "syn10220117"
obj <- synGet(id=synId, downloadFile = TRUE)
file <- getFileLocation(obj)
ohsu.genomic <- read.table(file, header=TRUE, sep="\t")
colnames(ohsu.genomic) <- gsub("X(.*)\\.(.*)", "\\1-\\2", colnames(ohsu.genomic))

## Read in the OHSU diagnosis info in the raw drug file and drug outcome file
## and limit samples to aml
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
ohsu.dss.orig <- subset(ohsu.dss.orig, lab_id %in% ohsu.diagnosis.tbl$lab_id)
ohsu.dss.orig <- subset(ohsu.dss.orig, lab_id %in% ohsu.inh.tbl$lab_id)

## Add the SeqID to the OHSU drug data.
## Read in the rna-seq summary table, which provides a map from the lab_id's used
## in the drug response data to the SeqID's used in the expression data.
synId <- "syn10083817"
obj <- synGet(synId, downloadFile = TRUE)
file <- getFileLocation(obj)
ohsu.rnaseq.sample.summary <- read.table(file, header=TRUE, sep="\t")

ohsu.dss.orig <- merge(ohsu.dss.orig, ohsu.rnaseq.sample.summary, by.x = "lab_id", by.y = "Original_LabID")

## Read in the FIMM expression data (RNASeq.CPM.log2.bc.csv)
## Load (batch corrected) FIMM expression data
obj <- synGet(id="syn8270602", downloadFile = TRUE)
file <- getFileLocation(obj)
fimm.expr <- read.table(file, header=TRUE, sep=",")
## Make the columns samples
fimm.expr <- t(fimm.expr)

## Read in the FIMM genomic data
## The "bin" file simply is 0/1 matrix indicating whether the sample (row) has a mutation
## in the gene (column).  These use gene symbols.
## NB: a 1 indicates that the gene had a non-synonymous or splice site mutation.
## path <- "fimm.genomic.ns.ss.tsv"
## Columns correspond to genes.
synId <- "syn10220160"
obj <- synGet(id=synId, downloadFile = TRUE)
file <- getFileLocation(obj)
fimm.genomic <- read.table(file, header=TRUE, sep="\t")
## Make the columns samples
fimm.genomic <- t(fimm.genomic)

## Read in the FIMM Metadata and ensure that we only analyze AML data
synId <- "syn8270594"
obj <- synGet(id=synId, downloadFile = TRUE)
file <- getFileLocation(obj)
fimm.metadata <- read.table(file, header=TRUE, sep=",")
fimm.metadata <- subset(fimm.metadata, grepl(diagnosis, pattern="AML"))

fimm.genomic <- fimm.genomic[, colnames(fimm.genomic) %in% rownames(fimm.metadata)]
fimm.expr <- fimm.expr[, colnames(fimm.expr) %in% rownames(fimm.metadata)]

## Read in the CTRP expression data
synId <- "syn5616092"
obj <- synGet(synId, downloadFile = TRUE)
file <- getFileLocation(obj)
ctrp.expr <- read.table(file, header=TRUE, sep=",")
rownames(ctrp.expr) <- ctrp.expr[,1]
ctrp.expr <- ctrp.expr[,-1]

## Read in the CTRP cell line metadata so we can limit (or not) cell lines by disease type
## CTRP cell line metadata: "v20.meta.per_cell_line.txt" (syn5632192)
synId <- "syn5632192"
obj <- synGet(synId, downloadFile = TRUE)
ctrp.cell.line.metadata <- read.table(getFileLocation(obj), sep="\t", header=TRUE, as.is=TRUE, comment.char="", quote="")

## Samples in CTRP do not have "."s and are all in caps
colnames(ctrp.expr) <- gsub(x=colnames(ctrp.expr), pattern="\\.", replacement="")
colnames(ctrp.expr) <- toupper(colnames(ctrp.expr))
ctrp.dss.orig <- merge(ctrp.dss.orig, ctrp.cell.line.metadata[, c("master_ccl_id", "ccl_name")])
## all.ccls <- intersect(colnames(ctrp.expr), ctrp.dss.orig$ccl_name)
all.ccls <- unique(ctrp.dss.orig$ccl_name)
heme.ccls <- all.ccls[all.ccls %in% ctrp.cell.line.metadata[ctrp.cell.line.metadata$ccle_primary_hist == "haematopoietic_neoplasm", "ccl_name"]]
aml.ccls <- all.ccls[all.ccls %in% ctrp.cell.line.metadata[ctrp.cell.line.metadata$ccle_hist_subtype_1 == "acute_myeloid_leukaemia", "ccl_name"]]

## END read in preprocessed OHSU, FIMM, and CTRVP2 data from Synapse

## BEGIN define drugs common to FIMM, OHSU, and CTRP

## Subset to analysis to drugs common to FIMM, OHSU, and CTRP
common.drugs <- ohsu.dss.orig[, c("inhibitor", "FIMM_Batch_ID_Drug.ll4", "FIMM_DRUG_NAME.ll4", "Mechanism.Targets.ll4", "Class.explained.ll4", "master_cpd_id.ll4")]
colnames(common.drugs) <- c("inhibitor", "DRUG_ID", "DRUG_NAME", "Mechanism", "Class", "master_cpd_id")
common.drugs <- subset(common.drugs, DRUG_ID %in% fimm.dss.orig$DRUG_ID)
common.drugs <- subset(common.drugs, master_cpd_id %in% ctrp.dss.orig$master_cpd_id)

un <- unique(common.drugs)
fimm.common.drugs <- un$DRUG_ID
ohsu.common.drugs <- un$inhibitor
ctrp.common.drugs <- un$master_cpd_id
drug.name.tbl <- un

fimm.drugs <- un$DRUG_ID
ohsu.drugs <- un$inhibitor
ctrp.drugs <- un$master_cpd_id

ohsu.query.drugs <- ohsu.common.drugs
fimm.query.drugs <- fimm.common.drugs
ctrp.query.drugs <- ctrp.common.drugs

## END define drugs common to FIMM, OHSU, and CTRP

## BEGIN modeling setup

source("../common/data-preprocessing.R")
source("models.R")
source("../common/plotting.R")

## Define CTRP AML, Heme, and All cell lines.  AML is a subset of Heme is a subset of all.
all.expr.cols <- intersect(colnames(ctrp.expr), all.ccls)
heme.expr.cols <- intersect(colnames(ctrp.expr), heme.ccls)
aml.expr.cols <- intersect(colnames(ctrp.expr), aml.ccls)

## Log-transform the CTRP data (FIMM and OHSU are already in log space)
min.expr <- min(ctrp.expr[ctrp.expr != 0])
ctrp.log2.expr <- ctrp.expr
ctrp.log2.expr[ctrp.log2.expr == 0] <- min.expr
ctrp.log2.expr <- log2(ctrp.log2.expr)

## Define the drug response variable.
## Use the AUC value computing using a log-logistic (LL.4) fit.  
## "dss." implies that AUC was calculated over the concentration range common
## to all three data sets.
test.response.col <- "dss.auc.ll4"
train.response.col <- test.response.col
ctrp.response.col <- test.response.col

## Filter the gene expression data to restrict to highly-expressed genes
iqrs <- unlist(apply(ctrp.log2.expr, 1, IQR))
ctrp.expr.filt <- ctrp.log2.expr[iqrs > median(iqrs[iqrs > 0]),]

iqrs <- unlist(apply(ohsu.expr, 1, IQR))
ohsu.expr.filt <- ohsu.expr[iqrs > median(iqrs[iqrs > 0]),]

iqrs <- unlist(apply(fimm.expr, 1, IQR))
fimm.expr.filt <- fimm.expr[iqrs > median(iqrs[iqrs > 0]),]

## Restrict analysis to genes common to all 3 data sets
common.genes <- intersect(rownames(ctrp.expr.filt), rownames(ohsu.expr.filt))
common.genes <- intersect(common.genes, rownames(fimm.expr.filt))

## plot(density(unlist(na.omit(ohsu.expr.filt[common.genes, ]))))
## lines(density(unlist(na.omit(ctrp.expr.filt[common.genes, ]))))
## lines(density(unlist(na.omit(fimm.expr.filt[common.genes, ]))))

## plot(density(na.omit(ohsu.dss.orig[, "auc.ll4"])))
## lines(density(na.omit(fimm.dss.orig[, "auc.ll4"])))
## lines(density(na.omit(ctrp.dss.orig[, ctrp.response.col])))

## Split the OHSU data into training (60%) and test (40%) -- do this on a per-drug basis because of the
## missing data.
set.seed(1234)
train.index <- createDataPartition(ohsu.dss.orig$inhibitor, p = .6, list = FALSE)
ohsu.dss.orig.train <- ohsu.dss.orig[train.index, ]
ohsu.dss.orig.test <- ohsu.dss.orig[-train.index, ]

num.ohsu.samples <- length(unique(intersect(ohsu.dss.orig$SeqID, colnames(ohsu.expr))))
num.ohsu.train.samples <- length(unique(intersect(ohsu.dss.orig.train$SeqID, colnames(ohsu.expr))))
num.ohsu.test.samples <- length(unique(intersect(ohsu.dss.orig.test$SeqID, colnames(ohsu.expr))))
rm(ohsu.dss.orig)

num.fimm.samples <- length(unique(intersect(fimm.dss.orig$SCREEN_ID, colnames(fimm.expr))))
ctrp.samples <- unique(intersect(ctrp.dss.orig$ccl_name, colnames(ctrp.log2.expr)))
num.ctrp.all.samples <- length(ctrp.samples)
num.ctrp.heme.samples <- length(intersect(ctrp.samples, heme.expr.cols))
num.ctrp.aml.samples <- length(intersect(ctrp.samples, aml.expr.cols))

## END modeling setup

## BEGIN modeling

ohsu.drug.col <- "inhibitor"
ohsu.patient.col <- "SeqID"
fimm.drug.col <- "DRUG_ID"
fimm.patient.col <- "SCREEN_ID"
ctrp.drug.col <- "master_cpd_id"
ctrp.patient.col <- "ccl_name"

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

## ras.sig.genes <- data.frame(symbols = c("NPTX2", "LMAN2L", "IL33", "DUSP6", "MYL4", "PHLDA1", "EYA2", "PHIP", "KCNN4", "HSPB3", "SKP1", "LAMC2", "HOXB6", "SKAP1", "RNF39", "SLC17A1", "BANP", "TBXAS1", "KCNG1", "TRPV1", "ZBTB40", "EDEM2", "SNX3", "NINL", "HDLBP", "MYBPC1", "S100A6", "SYNJ1", "GABBR1", "APOBEC1", "FAM60A", "ZBTB48", "HOXB7", "ASGR1", "SYCP1", "EREG", "CDKAL1", "BNC1", "MFSD11", "C10orf68", "KRT81", "LRRC6", "IFNGR2", "KRT24", "STAMBP", "PRUNE2", "GABRQ", "TNNI1", "CAPN9"))
## ras.sig.genes <- na.omit(symbols.to.ensg.mapping(ras.sig.genes$symbols))


## save.image(".Rdata")

gene.sets <- list("gene" = "gene", "hallmark" = hallmark.gene.sets, "biocarta" = biocarta.gene.sets, "kegg" = kegg.gene.sets, "reactome" = reactome.gene.sets)
train.set.names <- list("ctrp.all", "ctrp.heme", "ctrp.aml", "ohsu.train")
train.dss.args <- list(ctrp.dss.orig, ctrp.dss.orig, ctrp.dss.orig, ohsu.dss.orig.train)
train.expr.args <- list(ctrp.expr.filt[common.genes, ], ctrp.expr.filt[common.genes, heme.expr.cols], ctrp.expr.filt[common.genes, aml.expr.cols], ohsu.expr.filt[common.genes, ])
train.genomic.args <- list(NULL, NULL, NULL, NULL)
train.clinical.args <- list(NULL, NULL, NULL, NULL)
train.common.drugs <- list(ctrp.common.drugs, ctrp.common.drugs, ctrp.common.drugs, ohsu.common.drugs)
train.drugs <- list(ctrp.query.drugs, ctrp.query.drugs, ctrp.query.drugs, ohsu.query.drugs)
train.drug.cols <- list(ctrp.drug.col, ctrp.drug.col, ctrp.drug.col, ohsu.drug.col)
train.patient.cols <- list(ctrp.patient.col, ctrp.patient.col, ctrp.patient.col, ohsu.patient.col)
train.response.cols <- list("dss.auc.ll4", "dss.auc.ll4", "dss.auc.ll4", "dss.auc.ll4")

if(FALSE) {
n.drugs <- 3
gene.sets <- list("gene" = "gene", "hallmark" = hallmark.gene.sets)
train.set.names <- list("ctrp.aml", "ohsu.train")
train.dss.args <- list(ctrp.dss.orig, ohsu.dss.orig.train)
train.expr.args <- list(ctrp.expr.filt[common.genes, aml.expr.cols], ohsu.expr.filt[common.genes, ])
train.genomic.args <- list(NULL, NULL)
train.clinical.args <- list(NULL, NULL)
train.common.drugs <- list(ctrp.common.drugs, ohsu.common.drugs)
train.drugs <- list(ctrp.query.drugs[1:n.drugs], ohsu.query.drugs[1:n.drugs])
train.drug.cols <- list(ctrp.drug.col, ohsu.drug.col)
train.patient.cols <- list(ctrp.patient.col, ohsu.patient.col)
train.response.cols <- list("dss.auc.ll4", "dss.auc.ll4")
}

test.set.names <- list("ohsu.test", "fimm")
test.dss.args <- list(ohsu.dss.orig.test, fimm.dss.orig)
test.expr.args <- list(ohsu.expr.filt[common.genes, ], fimm.expr.filt[common.genes, ])
test.genomic.args <- list(NULL, NULL)
test.clinical.args <- list(NULL, NULL)
test.common.drugs <- list(ohsu.common.drugs, fimm.common.drugs)
test.drug.cols <- list(ohsu.drug.col, fimm.drug.col)
test.patient.cols <- list(ohsu.patient.col, fimm.patient.col)
test.response.cols <- list("dss.auc.ll4", "dss.auc.ll4")
num.processes <- num.cores - 1

cat("Training and testing with AUC\n")
train.response.cols <- list("dss.auc.ll4", "dss.auc.ll4", "dss.auc.ll4", "dss.auc.ll4")
test.response.cols <- list("dss.auc.ll4", "dss.auc.ll4")
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
train.response.cols <- list("dss2.ll4", "dss2.ll4", "dss2.ll4", "dss2.ll4")
test.response.cols <- list("dss2.ll4", "dss2.ll4")
res.dss2 <- train.and.test.crossproduct(gene.sets = gene.sets, drug.name.tbl = drug.name.tbl, train.set.names = train.set.names, train.dss.args = train.dss.args, 
                                        train.expr.args = train.expr.args, train.genomic.args = train.genomic.args, train.clinical.args = train.clinical.args, 
                                        train.common.drugs = train.common.drugs, train.drugs = train.drugs, train.drug.cols = train.drug.cols, train.patient.cols = train.patient.cols, 
                                        train.response.cols = train.response.cols, test.set.names = test.set.names, test.dss.args = test.dss.args, test.expr.args = test.expr.args, 
                                        test.genomic.args = test.genomic.args, test.clinical.args = test.clinical.args, test.common.drugs = test.common.drugs,
                                        test.drug.cols = test.drug.cols, test.patient.cols = test.patient.cols, test.response.cols = test.response.cols, 
                                        seed = 1234, use.rf = FALSE, use.svm = FALSE, num.processes = num.processes) 
cat("Done training and testing with DSS2\n")

save.image(".Rdata")

stop("successfully completed")

tbl <- ldply(train.set.names,
             .fun = function(train.set) {
                      ldply(test.set.names, 
                            .fun = function(test.set) {
                                     ldply(names(gene.sets),
                                           .fun = function(gene.set) {
                                                    tmp.dss2 <- extract.results(res.dss2[[train.set]][[test.set]][[gene.set]])
                                                    tmp.dss2$train.set <- train.set
                                                    tmp.dss2$test.set <- test.set
                                                    tmp.dss2$gene.set <- gene.set
                                                    tmp.dss2$response <- "DSS2"
                                                    tmp.auc <- extract.results(res.auc[[train.set]][[test.set]][[gene.set]])
                                                    tmp.auc$train.set <- train.set
                                                    tmp.auc$test.set <- test.set
                                                    tmp.auc$gene.set <- gene.set
                                                    tmp.auc$response <- "AUC"
                                                    rbind(tmp.auc, tmp.dss2)                         
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

## Create pairwise correlation plots
for(test.set.name in unique(tbl$test.set)) {
  for(mdl in c("lasso", "ridge")) {
    for(resp in unique(tbl$response)) {
      pdf(paste(test.set.name, mdl, resp, "corr.pdf", sep="-"))
      for(gs in unique(tbl$gene.set)) {
        mat <- do.call("cbind", corr.lists[[test.set.name]][[mdl]][[resp]][[gs]])
        cols <- NULL
        cls <- drug.name.tbl$Class
        if(test.set.name == "fimm") {
          names(cls) <- drug.name.tbl$DRUG_ID
        } else {
          names(cls) <- drug.name.tbl$inhibitor
        }
        cols <- factor(cls[names(corr.lists[[test.set.name]][[mdl]][[resp]][[gs]][[1]])])
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

panel.diff <- function (x, y, col = par("col"), bg = NA, pch = par("pch"), 
                              cex = 1, col.regres = "red", ...) 
{ 
  ok <- is.finite(x) & is.finite(y) 
print(x-y)
vec <- x - y
  points(1:length(vec), vec, pch = pch, col = col, bg = bg, cex = cex) 
##  abline(h = 0, col = col.regres, ...) 
} 

for(test.set.name in unique(tbl$test.set)) {
  for(mdl in c("lasso", "ridge")) {
    for(resp in unique(tbl$response)) {
      pdf(paste(test.set.name, mdl, resp, "diff.pdf", sep="-"))
      for(gs in unique(tbl$gene.set)) {
        mat <- do.call("cbind", corr.lists[[test.set.name]][[mdl]][[resp]][[gs]])
        cols <- NULL
        cls <- drug.name.tbl$Class
        if(test.set.name == "fimm") {
          names(cls) <- drug.name.tbl$DRUG_ID
        } else {
          names(cls) <- drug.name.tbl$inhibitor
        }
        len <- length(corr.lists[[test.set.name]][[mdl]][[resp]][[gs]])
        cols <- factor(cls[names(corr.lists[[test.set.name]][[mdl]][[resp]][[gs]][[1]])])
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
              g <- ggplot(df, aes(x = x, y = y, colour = col))
              g <- g + geom_point()
              g <- g + theme(legend.position="none")
              if(i < j) {
                g <- g + geom_hline(yintercept = 0)
              }
            } else {
              g <- qplot(1,1) + annotation_custom(textGrob(names(corr.lists[[test.set.name]][[mdl]][[resp]][[gs]])[i]),
                               xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) 
            }
            glist[[length(glist)+1]] <- g
          }
        }
        do.call("grid.arrange", glist)
##        pdf(paste(test.set.name, mdl, resp, gs, "diff.pdf", sep="-"))
        pairs(mat, lower.panel=panel.diff, upper.panel = panel.diff, main = paste(test.set.name, mdl, resp, gs, sep=" "), col = cols, oma=c(4,4,6,12))
        par(xpd=TRUE)
        legend(0.85, 0.7, legend = levels(cols), fill = c(1:length(levels(cols))), cex = 0.7)
##        d <- dev.off()
      }
      d <- dev.off()
    }
  }
}

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

## Drugs with predictions positively correlated with actual response in AML and with greater correlations in AML than in ALL
subset(df[, c("name", "Mechanism", "Class", "fimm", "ohsu", "fimm.ctrp.aml", "fimm.ctrp.all", "ohsu.ctrp.aml", "ohsu.ctrp.all")], (fimm > 0) & (ohsu > 0) & (fimm.ctrp.aml > 0) & (ohsu.ctrp.aml > 0))

## Drugs with predictions positively correlated with actual response in ALL and with greater correlations in ALL than in AML
subset(df[, c("name", "Mechanism", "Class", "fimm", "ohsu", "fimm.ctrp.aml", "fimm.ctrp.all", "ohsu.ctrp.aml", "ohsu.ctrp.all")], (fimm < 0) & (ohsu < 0) & (fimm.ctrp.all > 0) & (ohsu.ctrp.all > 0))

print(g)


## Create pairwise correlation plots and then look at diffs


## Add abline
## Consider coloring by targeted


capwords <- function(s, strict = FALSE) {
    cap <- function(s) paste(toupper(substring(s, 1, 1)),
                  {s <- substring(s, 2); if(strict) tolower(s) else s},
                             sep = "", collapse = " " )
    sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}

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
      file <- paste0(resp, "-", mdl, "-", met, "-vs-gene-sets.pdf")
      pdf(file)
      print(g)
      d <- dev.off()
##    glist[[length(glist)+1]] <- g
    }
  }
##  do.call("grid.arrange", glist)
}



## WORKING
## - convert biocarta, etc.
## - run auc (all)
## - run dss2 (all)
## - combine into plot

stop("stop")

rets <- list()

cat("Train on OHSU; Test on FIMM\n")



ret <- model.train.and.test(train.dss.arg = ohsu.dss.orig, train.expr.arg = ohsu.expr.filt[common.genes, ], train.genomic.arg = NULL, train.clinical.arg = NULL, 
                            train.common.drugs = ohsu.common.drugs, train.drugs = ohsu.query.drugs, train.drug.col = "inhibitor", 
                            train.patient.col = "SeqID", train.response.col = train.response.col, 
                            test.dss.arg = fimm.dss.orig, test.expr.arg = fimm.expr.filt[common.genes, ], test.genomic.arg = NULL, test.clinical.arg = NULL, 
                            test.common.drugs = fimm.common.drugs, test.drugs = fimm.query.drugs, test.drug.col = "DRUG_ID", 
                            test.patient.col = "SCREEN_ID", test.response.col = test.response.col, seed = 1234, use.rf = FALSE, use.svm = FALSE)  
rets[["ohsu.vs.fimm"]] <- ret

cat("Train on OHSU; Test on OHSU\n")
ret <- model.train.and.test(train.dss.arg = ohsu.dss.orig, train.expr.arg = ohsu.expr.filt[common.genes, ], train.genomic.arg = NULL, train.clinical.arg = NULL, 
                            train.common.drugs = ohsu.common.drugs, train.drugs = ohsu.query.drugs, train.drug.col = ohsu.drug.col,
                            train.patient.col = ohsu.patient.col, train.response.col = train.response.col, 
                            test.dss.arg = ohsu.dss.orig, test.expr.arg = ohsu.expr.filt[common.genes, ], test.genomic.arg = NULL, test.clinical.arg = NULL, 
                            test.common.drugs = ohsu.common.drugs, test.drugs = ohsu.query.drugs, test.drug.col = ohsu.drug.col,
                            test.patient.col = ohsu.patient.col, test.response.col = test.response.col, seed = 1234, use.rf = FALSE, use.svm = FALSE)  
rets[["ohsu.vs.ohsu"]] <- ret

cat("Train on fimm; Test on fimm\n")
ret <- model.train.and.test(train.dss.arg = fimm.dss.orig, train.expr.arg = fimm.expr.filt[common.genes, ], train.genomic.arg = NULL, train.clinical.arg = NULL, 
                            train.common.drugs = fimm.common.drugs, train.drugs = fimm.query.drugs, train.drug.col = fimm.drug.col,
                            train.patient.col = fimm.patient.col, train.response.col = train.response.col, 
                            test.dss.arg = fimm.dss.orig, test.expr.arg = fimm.expr.filt[common.genes, ], test.genomic.arg = NULL, test.clinical.arg = NULL, 
                            test.common.drugs = fimm.common.drugs, test.drugs = fimm.query.drugs, test.drug.col = fimm.drug.col,
                            test.patient.col = fimm.patient.col, test.response.col = test.response.col, seed = 1234, use.rf = FALSE, use.svm = FALSE)  
rets[["fimm.vs.fimm"]] <- ret

cat("Train on CTRP (All); Test on CTRP (All)\n")
ret <- model.train.and.test(train.dss.arg = ctrp.dss.orig, train.expr.arg = ctrp.expr.filt[common.genes, ], train.genomic.arg = NULL, train.clinical.arg = NULL, 
                            train.common.drugs = ctrp.common.drugs, train.drugs = ctrp.query.drugs, train.drug.col = ctrp.drug.col,
                            train.patient.col = ctrp.patient.col, train.response.col = train.response.col, 
                            test.dss.arg = ctrp.dss.orig, test.expr.arg = ctrp.expr.filt[common.genes, ], test.genomic.arg = NULL, test.clinical.arg = NULL, 
                            test.common.drugs = ctrp.common.drugs, test.drugs = ctrp.query.drugs, test.drug.col = ctrp.drug.col,
                            test.patient.col = ctrp.patient.col, test.response.col = test.response.col, seed = 1234, use.rf = FALSE, use.svm = FALSE)  
rets[["ctrp.all.vs.ctrp.all"]] <- ret


cat("Train on CTRP (Heme); Test on CTRP (Heme)\n")
ret <- model.train.and.test(train.dss.arg = ctrp.dss.orig, train.expr.arg = ctrp.expr.filt[common.genes, heme.expr.cols], train.genomic.arg = NULL, train.clinical.arg = NULL, 
                            train.common.drugs = ctrp.common.drugs, train.drugs = ctrp.query.drugs, train.drug.col = ctrp.drug.col,
                            train.patient.col = ctrp.patient.col, train.response.col = train.response.col, 
                            test.dss.arg = ctrp.dss.orig, test.expr.arg = ctrp.expr.filt[common.genes, heme.expr.cols], test.genomic.arg = NULL, test.clinical.arg = NULL, 
                            test.common.drugs = ctrp.common.drugs, test.drugs = ctrp.query.drugs, test.drug.col = ctrp.drug.col,
                            test.patient.col = ctrp.patient.col, test.response.col = test.response.col, seed = 1234, use.rf = FALSE, use.svm = FALSE)  
rets[["ctrp.heme.vs.ctrp.heme"]] <- ret

cat("Train on CTRP (AML); Test on CTRP (AML)\n")
ret <- model.train.and.test(train.dss.arg = ctrp.dss.orig, train.expr.arg = ctrp.expr.filt[common.genes, aml.expr.cols], train.genomic.arg = NULL, train.clinical.arg = NULL, 
                            train.common.drugs = ctrp.common.drugs, train.drugs = ctrp.query.drugs, train.drug.col = ctrp.drug.col,
                            train.patient.col = ctrp.patient.col, train.response.col = train.response.col, 
                            test.dss.arg = ctrp.dss.orig, test.expr.arg = ctrp.expr.filt[common.genes, aml.expr.cols], test.genomic.arg = NULL, test.clinical.arg = NULL, 
                            test.common.drugs = ctrp.common.drugs, test.drugs = ctrp.query.drugs, test.drug.col = ctrp.drug.col,
                            test.patient.col = ctrp.patient.col, test.response.col = test.response.col, seed = 1234, use.rf = FALSE, use.svm = FALSE)  
rets[["ctrp.aml.vs.ctrp.aml"]] <- ret


cat("Train on CTRP (All); Test on FIMM\n")
ret <- model.train.and.test(train.dss.arg = ctrp.dss.orig, train.expr.arg = ctrp.expr.filt[common.genes, ], train.genomic.arg = NULL, train.clinical.arg = NULL, 
                            train.common.drugs = ctrp.common.drugs, train.drugs = ctrp.query.drugs, train.drug.col = "master_cpd_id",
                            train.patient.col = "ccl_name", train.response.col = train.response.col, 
                            test.dss.arg = fimm.dss.orig, test.expr.arg = fimm.expr.filt[common.genes, ], test.genomic.arg = NULL, test.clinical.arg = NULL, 
                            test.common.drugs = fimm.common.drugs, test.drugs = fimm.query.drugs, test.drug.col = "DRUG_ID", 
                            test.patient.col = "SCREEN_ID", test.response.col = test.response.col, seed = 1234, use.rf = FALSE, use.svm = FALSE)  
rets[["ctrp.all.vs.fimm"]] <- ret

cat("Train on CTRP (Heme); Test on FIMM\n")
ret <- model.train.and.test(train.dss.arg = ctrp.dss.orig, train.expr.arg = ctrp.expr.filt[common.genes, heme.expr.cols], train.genomic.arg = NULL, train.clinical.arg = NULL, 
                            train.common.drugs = ctrp.common.drugs, train.drugs = ctrp.query.drugs, train.drug.col = "master_cpd_id",
                            train.patient.col = "ccl_name", train.response.col = train.response.col, 
                            test.dss.arg = fimm.dss.orig, test.expr.arg = fimm.expr.filt[common.genes, ], test.genomic.arg = NULL, test.clinical.arg = NULL, 
                            test.common.drugs = fimm.common.drugs, test.drugs = fimm.query.drugs, test.drug.col = "DRUG_ID", 
                            test.patient.col = "SCREEN_ID", test.response.col = test.response.col, seed = 1234, use.rf = FALSE, use.svm = FALSE)  
rets[["ctrp.heme.vs.fimm"]] <- ret

cat("Train on CTRP (AML); Test on FIMM\n")
ret <- model.train.and.test(train.dss.arg = ctrp.dss.orig, train.expr.arg = ctrp.expr.filt[common.genes, aml.expr.cols], train.genomic.arg = NULL, train.clinical.arg = NULL, 
                            train.common.drugs = ctrp.common.drugs, train.drugs = ctrp.query.drugs, train.drug.col = "master_cpd_id",
                            train.patient.col = "ccl_name", train.response.col = train.response.col, 
                            test.dss.arg = fimm.dss.orig, test.expr.arg = fimm.expr.filt[common.genes, ], test.genomic.arg = NULL, test.clinical.arg = NULL, 
                            test.common.drugs = fimm.common.drugs, test.drugs = fimm.query.drugs, test.drug.col = "DRUG_ID", 
                            test.patient.col = "SCREEN_ID", test.response.col = test.response.col, seed = 1234, use.rf = FALSE, use.svm = FALSE)  
rets[["ctrp.aml.vs.fimm"]] <- ret

cat("Train on FIMM; Test on FIMM\n")

td <- extract.results(rets[["ohsu.vs.fimm"]])
ohsu.label  <- paste0("OHSU\n(n=", num.ohsu.samples, ")")
fimm.label <- paste0("FIMM\n(n=", num.fimm.samples, ")")
td$train.set <- ohsu.label
td$test.set <- fimm.label

te <- extract.results(rets[["ctrp.all.vs.fimm"]])
ctrp.all.label <- paste0("CTRP (All)\n(n=", num.ctrp.all.samples, ")")
te$train.set <- ctrp.all.label
te$test.set <- fimm.label

tf <- extract.results(rets[["ctrp.heme.vs.fimm"]])
ctrp.heme.label <- paste0("CTRP (Heme)\n(n=", num.ctrp.heme.samples, ")")
tf$train.set <- ctrp.heme.label
tf$test.set <- fimm.label

tg <- extract.results(rets[["ctrp.aml.vs.fimm"]])
ctrp.aml.label <- paste0("CTRP (AML)\n(n=", num.ctrp.aml.samples, ")")
tg$train.set <- ctrp.aml.label
tg$test.set <- fimm.label

ti <- extract.results(rets[["fimm.vs.fimm"]])
ti$train.set <- fimm.label
ti$test.set <- fimm.label

num.samples <- unique(tg[,c("test.drug", "n.train")])
num.samples$test.drug <- as.character(num.samples$test.drug)
num.samples.per.drug <- num.samples$n.train
names(num.samples.per.drug) <- num.samples$test.drug

ns2 <- unique(td[, c("train.drug", "test.drug", "n.train")])
ns2$test.drug <- as.character(ns2$test.drug)

num.samples <- merge(num.samples, ns2, by = c("test.drug"), suffixes = c(".ctrp.aml", ".ohsu"))

## We can only downsample OHSU if it has more samples than AML.  Hence, exclude the only case in which
## this isn't true
cat("Num drugs for which n.train.ctrp.aml < n.train.ohsu\n")
table(num.samples$n.train.ctrp.aml < num.samples$n.train.ohsu)

sub <- subset(num.samples, n.train.ctrp.aml < n.train.ohsu)
ohsu.query.drugs.down <- as.character(sub$train.drug)
fimm.query.drugs.down <- as.character(sub$test.drug)

cat("Train on OHSU (with downsampling); Test on FIMM\n")
if(TRUE) {
ret <- model.train.and.test.with.downsampling(train.dss.arg = ohsu.dss.orig, train.expr.arg = ohsu.expr.filt[common.genes, ], train.genomic.arg = NULL, train.clinical.arg = NULL, 
                                              train.common.drugs = ohsu.common.drugs, train.drugs = ohsu.query.drugs.down, train.drug.col = "inhibitor", 
                                              train.patient.col = "SeqID", train.response.col = train.response.col, 
                                              test.dss.arg = fimm.dss.orig, test.expr.arg = fimm.expr.filt[common.genes, ], test.genomic.arg = NULL, test.clinical.arg = NULL, 
                                              test.common.drugs = fimm.common.drugs, test.drugs = fimm.query.drugs.down, test.drug.col = "DRUG_ID", 
                                              test.patient.col = "SCREEN_ID", test.response.col = test.response.col, seed = 1234, use.rf = FALSE, use.svm = FALSE,
                                              num.samples.per.drug = num.samples.per.drug, num.iterations = 50, with.replacement = FALSE)  
rets[["ohsu.vs.fimm.down.aml"]] <- ret
th <- extract.nested.results(rets[["ohsu.vs.fimm.down.aml"]])
ohsu.label.down  <- paste0("OHSU\n(down sample)")
th$train.set <- ohsu.label.down
th$test.set <- fimm.label
}



stop("stop")

## END modeling

## Make violin plots showing correlation of FIMM predicted vs FIMM actual (validation) when trained against the following 4 data sets:
## (1) CTRP (All)
## (2) CTRP (Heme)
## (3) CTRP (AML)
## (4) OHSU

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
  pdf(paste0("ctrp-ohsu-fimm-boxplot-", met, ".pdf"))
  print(g)
  d <- dev.off()

  g <- ggplot(data = tmp, aes(x = train.set, y = val))
  g <- g + geom_violin()
  g <- g + geom_beeswarm()
  g <- g + facet_wrap(~ model)
  g <- g + xlab(paste0("Training Data Set [Test Data Set = FIMM (n=", num.fimm.samples, ")]"))
  g <- g + ylab(met)
  g <- g + theme(axis.text.x=element_text(angle = 45, hjust = 1))
  pdf(paste0("ctrp-ohsu-fimm-beeswarm-", met, ".pdf"))
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
                     tbl <- merge(train.err, val.err, by = c("train.drug"), suffixes = c(".train", ".test"))
                     df <- tbl[, c("val.train", "val.test")]
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
    pdf(paste0(mdl, "-", met, "-train-vs-test-correlations.pdf"))
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
pdf("ctrp-ohsu-fimm-beeswarm-mse.pdf")
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