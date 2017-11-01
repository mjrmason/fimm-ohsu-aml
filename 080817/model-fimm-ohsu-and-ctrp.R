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
common.drugs <- ohsu.dss.orig[, c("inhibitor", "FIMM_Batch_ID_Drug.ll4", "FIMM_DRUG_NAME.ll4", "Mechanism.Targets.ll4", "master_cpd_id.ll4")]
colnames(common.drugs) <- c("ohsu.inhibitor", "fimm.DRUG_ID", "DRUG_NAME", "Mechanism", "master_cpd_id")
common.drugs <- subset(common.drugs, fimm.DRUG_ID %in% fimm.dss.orig$DRUG_ID)
common.drugs <- subset(common.drugs, master_cpd_id %in% ctrp.dss.orig$master_cpd_id)

un <- unique(common.drugs)
fimm.common.drugs <- un$fimm.DRUG_ID
ohsu.common.drugs <- un$ohsu.inhibitor
ctrp.common.drugs <- un$master_cpd_id

fimm.drugs <- un$fimm.DRUG_ID
ohsu.drugs <- un$ohsu.inhibitor
ctrp.drugs <- un$master_cpd_id

ohsu.query.drugs <- ohsu.common.drugs
fimm.query.drugs <- fimm.common.drugs
ctrp.query.drugs <- ctrp.common.drugs

## END define drugs common to FIMM, OHSU, and CTRP

## BEGIN modeling setup

source("../common/data-preprocessing.R")
source("models.R")

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
ctrp.log2.expr.filt <- ctrp.log2.expr[iqrs > median(iqrs[iqrs > 0]),]

iqrs <- unlist(apply(ohsu.expr, 1, IQR))
ohsu.expr.filt <- ohsu.expr[iqrs > median(iqrs[iqrs > 0]),]

iqrs <- unlist(apply(fimm.expr, 1, IQR))
fimm.expr.filt <- fimm.expr[iqrs > median(iqrs[iqrs > 0]),]

## Restrict analysis to genes common to all 3 data sets
common.genes <- intersect(rownames(ctrp.log2.expr.filt), rownames(ohsu.expr.filt))
common.genes <- intersect(common.genes, rownames(fimm.expr.filt))

## plot(density(unlist(na.omit(ohsu.expr.filt[common.genes, ]))))
## lines(density(unlist(na.omit(ctrp.log2.expr.filt[common.genes, ]))))
## lines(density(unlist(na.omit(fimm.expr.filt[common.genes, ]))))

## plot(density(na.omit(ohsu.dss.orig[, "auc.ll4"])))
## lines(density(na.omit(fimm.dss.orig[, "auc.ll4"])))
## lines(density(na.omit(ctrp.dss.orig[, ctrp.response.col])))

num.ohsu.samples <- length(unique(intersect(ohsu.dss.orig$SeqID, colnames(ohsu.expr))))
num.fimm.samples <- length(unique(intersect(fimm.dss.orig$SCREEN_ID, colnames(fimm.expr))))
ctrp.samples <- unique(intersect(ctrp.dss.orig$ccl_name, colnames(ctrp.log2.expr)))
num.ctrp.all.samples <- length(ctrp.samples)
num.ctrp.heme.samples <- length(intersect(ctrp.samples, heme.expr.cols))
num.ctrp.aml.samples <- length(intersect(ctrp.samples, aml.expr.cols))

## END modeling setup

## BEGIN modeling

rets <- list()

cat("Train on OHSU; Test on CTRP (All)\n")
ret <- model.train.and.test(train.dss.arg = ohsu.dss.orig, train.expr.arg = ohsu.expr[common.genes, ], train.genomic.arg = NULL, train.clinical.arg = NULL, 
                            train.common.drugs = ohsu.common.drugs, train.drugs = ohsu.query.drugs, train.drug.col = "inhibitor", 
                            train.patient.col = "SeqID", train.response.col = train.response.col, 
                            test.dss.arg = ctrp.dss.orig, test.expr.arg = ctrp.log2.expr[common.genes, ], test.genomic.arg = NULL, test.clinical.arg = NULL, 
                            test.common.drugs = ctrp.common.drugs, test.drugs = ctrp.query.drugs, test.drug.col = "master_cpd_id", 
                            test.patient.col = "ccl_name", test.response.col = test.response.col, seed = 1234) 
rets[["ohsu.vs.ctrp.all"]] <- ret
ta <- extract.results(rets[["ohsu.vs.ctrp.all"]])
ta
ta$cmp <- "CTRP (All)"

cat("Train on OHSU; Test on CTRP (Heme)\n")
ret <- model.train.and.test(train.dss.arg = ohsu.dss.orig, train.expr.arg = ohsu.expr[common.genes, ], train.genomic.arg = NULL, train.clinical.arg = NULL, 
                            train.common.drugs = ohsu.common.drugs, train.drugs = ohsu.query.drugs, train.drug.col = "inhibitor", 
                            train.patient.col = "SeqID", train.response.col = train.response.col, 
                            test.dss.arg = ctrp.dss.orig, test.expr.arg = ctrp.log2.expr[common.genes, heme.expr.cols], test.genomic.arg = NULL, test.clinical.arg = NULL, 
                            test.common.drugs = ctrp.common.drugs, test.drugs = ctrp.query.drugs, test.drug.col = "master_cpd_id", 
                            test.patient.col = "ccl_name", test.response.col = test.response.col, seed = 1234) 
rets[["ohsu.vs.ctrp.heme"]] <- ret
tb <- extract.results(rets[["ohsu.vs.ctrp.heme"]])
tb$cmp <- "CTRP (Heme)"

cat("Train on OHSU; Test on CTRP (AML)\n")
ret <- model.train.and.test(train.dss.arg = ohsu.dss.orig, train.expr.arg = ohsu.expr[common.genes, ], train.genomic.arg = NULL, train.clinical.arg = NULL, 
                            train.common.drugs = ohsu.common.drugs, train.drugs = ohsu.query.drugs, train.drug.col = "inhibitor", 
                            train.patient.col = "SeqID", train.response.col = train.response.col, 
                            test.dss.arg = ctrp.dss.orig, test.expr.arg = ctrp.log2.expr[common.genes, aml.expr.cols], test.genomic.arg = NULL, test.clinical.arg = NULL, 
                            test.common.drugs = ctrp.common.drugs, test.drugs = ctrp.query.drugs, test.drug.col = "master_cpd_id", 
                            test.patient.col = "ccl_name", test.response.col = test.response.col, seed = 1234) 
rets[["ohsu.vs.ctrp.aml"]] <- ret
tc <- extract.results(rets[["ohsu.vs.ctrp.aml"]])
tc$cmp <- "CTRP (AML)"

cat("Train on OHSU; Test on FIMM\n")
ret <- model.train.and.test(train.dss.arg = ohsu.dss.orig, train.expr.arg = ohsu.expr[common.genes, ], train.genomic.arg = NULL, train.clinical.arg = NULL, 
                            train.common.drugs = ohsu.common.drugs, train.drugs = ohsu.query.drugs, train.drug.col = "inhibitor", 
                            train.patient.col = "SeqID", train.response.col = train.response.col, 
                            test.dss.arg = fimm.dss.orig, test.expr.arg = fimm.expr[common.genes, ], test.genomic.arg = NULL, test.clinical.arg = NULL, 
                            test.common.drugs = fimm.common.drugs, test.drugs = fimm.query.drugs, test.drug.col = "DRUG_ID", 
                            test.patient.col = "SCREEN_ID", test.response.col = test.response.col, seed = 1234) 
rets[["ohsu.vs.fimm"]] <- ret
td <- extract.results(rets[["ohsu.vs.fimm"]])
ohsu.label  <- paste0("OHSU\n(n=", num.ohsu.samples, ")")
td$cmp <- ohsu.label

cat("Train on CTRP (All); Test on FIMM\n")
ret <- model.train.and.test(train.dss.arg = ctrp.dss.orig, train.expr.arg = ctrp.expr[common.genes, ], train.genomic.arg = NULL, train.clinical.arg = NULL, 
                            train.common.drugs = ctrp.common.drugs, train.drugs = ctrp.query.drugs, train.drug.col = "master_cpd_id",
                            train.patient.col = "ccl_name", train.response.col = train.response.col, 
                            test.dss.arg = fimm.dss.orig, test.expr.arg = fimm.expr[common.genes, ], test.genomic.arg = NULL, test.clinical.arg = NULL, 
                            test.common.drugs = fimm.common.drugs, test.drugs = fimm.query.drugs, test.drug.col = "DRUG_ID", 
                            test.patient.col = "SCREEN_ID", test.response.col = test.response.col, seed = 1234) 
rets[["ctrp.all.vs.fimm"]] <- ret
te <- extract.results(rets[["ctrp.all.vs.fimm"]])
ctrp.all.label <- paste0("CTRP (All)\n(n=", num.ctrp.all.samples, ")")
te$cmp <- ctrp.all.label

cat("Train on CTRP (Heme); Test on FIMM\n")
ret <- model.train.and.test(train.dss.arg = ctrp.dss.orig, train.expr.arg = ctrp.expr[common.genes, heme.expr.cols], train.genomic.arg = NULL, train.clinical.arg = NULL, 
                            train.common.drugs = ctrp.common.drugs, train.drugs = ctrp.query.drugs, train.drug.col = "master_cpd_id",
                            train.patient.col = "ccl_name", train.response.col = train.response.col, 
                            test.dss.arg = fimm.dss.orig, test.expr.arg = fimm.expr[common.genes, ], test.genomic.arg = NULL, test.clinical.arg = NULL, 
                            test.common.drugs = fimm.common.drugs, test.drugs = fimm.query.drugs, test.drug.col = "DRUG_ID", 
                            test.patient.col = "SCREEN_ID", test.response.col = test.response.col, seed = 1234) 
rets[["ctrp.heme.vs.fimm"]] <- ret
tf <- extract.results(rets[["ctrp.heme.vs.fimm"]])
tf$cmp <- "CTRP (Heme)"
ctrp.heme.label <- paste0("CTRP (Heme)\n(n=", num.ctrp.heme.samples, ")")
tf$cmp <- ctrp.heme.label

cat("Train on CTRP (AML); Test on FIMM\n")
ret <- model.train.and.test(train.dss.arg = ctrp.dss.orig, train.expr.arg = ctrp.expr[common.genes, aml.expr.cols], train.genomic.arg = NULL, train.clinical.arg = NULL, 
                            train.common.drugs = ctrp.common.drugs, train.drugs = ctrp.query.drugs, train.drug.col = "master_cpd_id",
                            train.patient.col = "ccl_name", train.response.col = train.response.col, 
                            test.dss.arg = fimm.dss.orig, test.expr.arg = fimm.expr[common.genes, ], test.genomic.arg = NULL, test.clinical.arg = NULL, 
                            test.common.drugs = fimm.common.drugs, test.drugs = fimm.query.drugs, test.drug.col = "DRUG_ID", 
                            test.patient.col = "SCREEN_ID", test.response.col = test.response.col, seed = 1234) 
rets[["ctrp.aml.vs.fimm"]] <- ret
tg <- extract.results(rets[["ctrp.aml.vs.fimm"]])
tg$cmp <- "CTRP (AML)"
ctrp.aml.label <- paste0("CTRP (AML)\n(n=", num.ctrp.aml.samples, ")")
tg$cmp <- ctrp.aml.label

## END modeling

## Make violin plots showing correlation of FIMM predicted vs FIMM actual (validation) when trained against the following 4 data sets:
## (1) CTRP (All)
## (2) CTRP (Heme)
## (3) CTRP (AML)
## (4) OHSU

all.results <- rbind(td, te, tf, tg)
all.results$cmp <- factor(all.results$cmp, levels = c(ctrp.all.label, ctrp.heme.label, ctrp.aml.label, ohsu.label))
g <- ggplot(data = all.results, aes(x = cmp, y = cor))
g <- g + geom_violin()
g <- g + geom_boxplot(width = 0.1)
g <- g + facet_wrap(~ model)
g <- g + xlab(paste0("Training Data Set [Test Data Set = FIMM (n=", num.fimm.samples, ")]"))
g <- g + ylab("Correlation")
g <- g + theme(axis.text.x=element_text(angle = 45, hjust = 1))
pdf("ctrp-ohsu-fimm-boxplot.pdf")
print(g)
d <- dev.off()

all.results <- rbind(td, te, tf, tg)
all.results$cmp <- factor(all.results$cmp, levels = c(ctrp.all.label, ctrp.heme.label, ctrp.aml.label, ohsu.label))
g <- ggplot(data = all.results, aes(x = cmp, y = cor))
g <- g + geom_violin()
g <- g + geom_beeswarm()
g <- g + facet_wrap(~ model)
g <- g + xlab(paste0("Training Data Set [Test Data Set = FIMM (n=", num.fimm.samples, ")]"))
g <- g + ylab("Correlation")
g <- g + theme(axis.text.x=element_text(angle = 45, hjust = 1))
pdf("ctrp-ohsu-fimm-beeswarm.pdf")
print(g)
d <- dev.off()

ret <- model.train.and.test(train.dss.arg = ohsu.dss.orig, train.expr.arg = ohsu.expr[common.genes, ], train.genomic.arg = NULL, train.clinical.arg = NULL, 
                            train.common.drugs = ohsu.common.drugs, train.drugs = ohsu.query.drugs, train.drug.col = "inhibitor", 
                            train.patient.col = "SeqID", train.response.col = train.response.col, 
                            test.dss.arg = fimm.dss.orig, test.expr.arg = fimm.expr[common.genes, ], test.genomic.arg = NULL, test.clinical.arg = NULL, 
                            test.common.drugs = fimm.common.drugs, test.drugs = fimm.query.drugs, test.drug.col = "DRUG_ID", 
                            test.patient.col = "SCREEN_ID", test.response.col = test.response.col, seed = 1234) 
t2 <- extract.results(ret)
t2

ret <- model.train.and.test(train.dss.arg = ohsu.dss.orig, train.expr.arg = ohsu.expr[common.genes, ], train.genomic.arg = NULL, train.clinical.arg = NULL, 
                            train.common.drugs = ohsu.common.drugs, train.drugs = ohsu.query.drugs, train.drug.col = "inhibitor", 
                            train.patient.col = "SeqID", train.response.col = train.response.col, 
                            test.dss.arg = ctrp.dss.orig, test.expr.arg = ctrp.log2.expr[common.genes, aml.expr.cols], test.genomic.arg = NULL, test.clinical.arg = NULL, 
                            test.common.drugs = ctrp.common.drugs, test.drugs = ctrp.query.drugs, test.drug.col = "master_cpd_id", 
                            test.patient.col = "ccl_name", test.response.col = test.response.col, seed = 1234) 
t3 <- extract.results(ret)
t3


cols <- heme.expr.cols
cols <- all.expr.cols
half <- floor(length(cols)/2)
ctpr.train.samples <- cols[1:half]
ctpr.test.samples <- cols[!(cols %in% ctpr.train.samples)]
expr <- ctrp.log2.expr.filt
expr <- ctrp.log2.expr
ret <- model.train.and.test(train.dss.arg = ctrp.dss.orig, train.expr.arg = expr[common.genes, ctpr.train.samples], train.genomic.arg = NULL, train.clinical.arg = NULL, 
                            train.common.drugs = ctrp.common.drugs, train.drugs = ctrp.query.drugs, train.drug.col = "master_cpd_id", 
                            train.patient.col = "ccl_name", train.response.col = ctrp.response.col, 
                            test.dss.arg = ctrp.dss.orig, test.expr.arg = expr[common.genes, ctpr.test.samples], test.genomic.arg = NULL, test.clinical.arg = NULL, 
                            test.common.drugs = ctrp.common.drugs, test.drugs = ctrp.query.drugs, test.drug.col = "master_cpd_id", 
                            test.patient.col = "ccl_name", test.response.col = ctrp.response.col, seed = 1234) 
t4 <- extract.results(ret)
t4


