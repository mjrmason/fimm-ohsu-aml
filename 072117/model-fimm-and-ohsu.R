## Set the download and output directory
download.path <- "../input/"
if (!file.exists(download.path)) {
  dir.create(download.path)
}

output.path <- "output/"
if (!file.exists(output.path)) {
  dir.create(output.path)
}

beat.aml.cache <- "../input/beat-aml-files/"

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
suppressPackageStartupMessages(library(ReporteRs))
suppressPackageStartupMessages(library(rJava))

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

## Set the java heap to 8GB
##options(java.parameters = "-Xmx8000m")

## Define a java garbage collection function
jgc <- function()
{
  gc()
  .jcall("java/lang/System", method = "gc")
} 

cat("Simultaneously fit FIMM and OHSU DSS2 values using LASSO\n")

## parentId is the OHSU processed data folder
parentId <- "syn10083332"
tbl <- synQuery(paste0("SELECT name,id FROM file WHERE parentId==\'", parentId, "\'"))

## Read in the t=0 thresholded DSS values
file.name <- "ohsu.common.drugs.dss.t0.tsv"
synId <- tbl[tbl$file.name == file.name, "file.id"]
obj <- synGet(synId, downloadFile = TRUE)
file <- getFileLocation(obj)
ohsu.dss.orig <- as.data.frame(fread(file))

## parentId is the FIMM processed data folder
parentId <- "syn8270577"
tbl <- synQuery(paste0("SELECT name,id FROM file WHERE parentId==\'", parentId, "\'"))

## Read in the t=0 thresholded DSS values
file.name <- "fimm.common.drugs.dss.t0.tsv"
synId <- tbl[tbl$file.name == file.name, "file.id"]
obj <- synGet(synId, downloadFile = TRUE)
file <- getFileLocation(obj)
fimm.dss.orig <- as.data.frame(fread(file))

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


source("../common/data-preprocessing.R")

## Subset to analysis to drugs common to FIMM and OHSU
common.drugs <- ohsu.dss.orig[, c("inhibitor", "FIMM_Batch_ID_Drug.ll4", "FIMM_DRUG_NAME.ll4", "Mechanism.Targets.ll4")]
colnames(common.drugs) <- c("ohsu.inhibitor", "fimm.DRUG_ID", "DRUG_NAME", "Mechanism")
common.drugs <- subset(common.drugs, fimm.DRUG_ID %in% fimm.dss.orig$DRUG_ID)

un <- unique(common.drugs)
fimm.common.drugs <- un$fimm.DRUG_ID
ohsu.common.drugs <- un$ohsu.inhibitor

## Further restrict to MEK inhibitors
common.drugs <- subset(un, grepl(Mechanism, pattern="MEK"))
cat(paste0("Testing the following drugs: ", paste(unique(common.drugs$ohsu.inhibitor), collapse=", "), "\n"))

## Actually, just limit to a _single_ MEK inhibitor
## common.drugs <- common.drugs[1,,drop=F]

fimm.mek.drugs <- common.drugs$fimm.DRUG_ID
ohsu.mek.drugs <- common.drugs$ohsu.inhibitor

fimm.drugs <- un$fimm.DRUG_ID
ohsu.drugs <- un$ohsu.inhibitor

source("models.R")

## Model using DSS2 values from L.4 (logistic, not log logistic) fit.
## ret <- model.fimm.and.ohsu(ohsu.dss.orig, ohsu.expr, ohsu.drugs, fimm.dss.orig, fimm.expr, fimm.drugs, response.col = "dss2.l4")

ohsu.dss.arg <- ohsu.dss.orig
fimm.dss.arg <- fimm.dss.orig
rm(ohsu.dss.orig)
rm(fimm.dss.orig)
gc()
response.col <- "dss2.l4"

ohsu.query.drugs <- ohsu.common.drugs
fimm.query.drugs <- fimm.common.drugs

## ohsu.query.drugs <- ohsu.mek.drugs
## fimm.query.drugs <- fimm.mek.drugs

rets <- list()
responses <- c("dss2.l4", "auc.l4", "ic50.l4")
## responses <- c("dss2.l4")
data.types <- c("dna", "both", "rna")
## data.types <- c("dna")
for(response.col in responses) {
  rets[[response.col]] <- list()

  for(data.type in data.types) {
    ohsu.expr.arg <- NULL
    fimm.expr.arg <- NULL
    ohsu.genomic.arg <- NULL
    fimm.genomic.arg <- NULL

    if(data.type %in% c("rna", "both")) {
      ohsu.expr.arg <- ohsu.expr
      fimm.expr.arg <- fimm.expr
    }
    if(data.type %in% c("dna", "both")) {
      ohsu.genomic.arg <- ohsu.genomic
      fimm.genomic.arg <- fimm.genomic
    }
    cat(paste0("Fitting ", response.col, " using data ", data.type, "\n"))
    rets[[response.col]][[data.type]] <- model.fimm.and.ohsu(ohsu.dss.arg, ohsu.expr.arg, ohsu.genomic.arg, ohsu.common.drugs, ohsu.query.drugs, fimm.dss.arg, fimm.expr.arg, fimm.genomic.arg, fimm.common.drugs, fimm.query.drugs, response.col)
    gc()
    save.image(".Rdata")
  }
}

## response.col <- "dss2.l4"
## ret.dss2 <- model.fimm.and.ohsu(ohsu.dss.arg, ohsu.expr.arg, ohsu.common.drugs, ohsu.query.drugs, fimm.dss.arg, fimm.expr.arg, fimm.common.drugs, fimm.query.drugs, response.col)

## mek drugs
## PD184352
## Selumetinib
## Trametinib (this one is missing because the FIMM drug FIMM003751 has 2 drug sets, one of which only has a single concentration point)
##    PATIENT_ID                  SCREEN_ID  DRUG_SET    DRUG_ID min.conc.nM
##708  FHRB_1867 FHRB_1867_26102012_1030_BM  Fimm-01A FIMM003751        2.50
##760  FHRB_1867 FHRB_1867_26102012_1030_BM FO2AB1-aq FIMM003751        0.03
##    max.conc.nM
##708         2.5
##760       250.0

##    PATIENT_ID                  SCREEN_ID  DRUG_SET    DRUG_ID min.conc.nM
##586  FHRB_1690 FHRB_1690_12092012_1100_BM  Fimm-01A FIMM003751        2.50
##638  FHRB_1690 FHRB_1690_12092012_1100_BM FO2AB1-aq FIMM003751        0.03
##    max.conc.nM
##586         2.5
##638       250.0

cat("Finished modeling\n")

mydoc <- pptx()

doc.title <- paste0('Fit to OHSU; test on FIMM:\nIC50, AUC, and DSS (t = ', 0, ')')  
doc.file <- paste0("071817-fit-ohsu-test-fimm-t-", 0, ".pptx")
mydoc <- addSlide( mydoc, slide.layout = 'Title Slide' )
mydoc <- addTitle( mydoc, doc.title )
mydoc <- addSubtitle( mydoc , 'July 18, 2017')

tbls <- list()
for(response.col in names(rets)) {
  tbls[[response.col]] <- list()
  for(data.type in names(rets[[response.col]])) {
  tbls[[response.col]][[data.type]] <-
     ldply(rets[[response.col]][[data.type]], 
             .fun = function(foo) {
                ldply(foo,
             .fun = function(df) {
                ret <- NULL
                if(!is.na(df$predicted)) {
                  ct <- cor.test(df$predicted, df$actual)
                  ret <- data.frame(cor = ct$estimate, pval = ct$p.value, ohsu.drug = df$ohsu.drug, fimm.drug = df$fimm.drug, alpha = df$alpha, s = df$s, model = df$model)
                } else {
                  ret <- data.frame(cor = NA, pval = NA, ohsu.drug = df$ohsu.drug, fimm.drug = df$fimm.drug, alpha = df$alpha, s = df$s, model = df$model)
                }
                ret
             }) })
  metric <- response.col
  metric <- gsub(metric, pattern=".l4", replacement="")
  metric <- toupper(metric)

  tbl <- tbls[[response.col]][[data.type]]
  t2 <- subset(tbl, is.na(s) | (s == "lambda.min"))
  t2$model <- as.character(t2$model)
  t2[(t2$model == "glmnet") & (t2$alpha == 0),"model"] <- "ridge"
  t2[(t2$model == "glmnet") & (t2$alpha == 1),"model"] <- "lasso"
  t2 <- subset(t2, model != "glmnet")

  g <- plot.correlation.vs.drug(t2)
  g <- g + ggtitle(paste0(metric, " ", data.type))
  pdf(paste0(metric, "-", data.type, "-correlation-vs-drug.pdf"))
  print(g)
  d <- dev.off()

  mydoc <- addSlide( mydoc, slide.layout = "Title and Content" )
  title <- paste0(metric, " ", data.type, ": Correlation vs Drug")
  mydoc <- addTitle( mydoc, title )
  mydoc <- addPlot( doc = mydoc, fun = print, x = g)

  g <- plot.correlation.vs.model(t2)
  g <- g + ggtitle(paste0(metric, " ", data.type))
  pdf(paste0(metric, "-", data.type, "-correlation-vs-model.pdf"))
  print(g)
  d <- dev.off()

  mydoc <- addSlide( mydoc, slide.layout = "Title and Content" )
  title <- paste0(metric, " ", data.type, ": Correlation vs Model")
  mydoc <- addTitle( mydoc, title )
  mydoc <- addPlot( doc = mydoc, fun = print, x = g)
}
}

writeDoc(mydoc, doc.file)

cat("Finished successfully\n")

gc()
save.image(".Rdata")

cat("Finished writing .Rdata successfully\n")

## Only keep ridge, lasso, svm, and rf from table and rename with caps
sanitize.table <- function(tbl) {
  tbl.tmp <- subset(tbl, is.na(s) | (s == "lambda.min"))
  tbl.tmp$model <- as.character(tbl.tmp$model)
  tbl.tmp[(tbl.tmp$model == "glmnet") & (tbl.tmp$alpha == 0),"model"] <- "Ridge"
  tbl.tmp[(tbl.tmp$model == "glmnet") & (tbl.tmp$alpha == 1),"model"] <- "LASSO"
  tbl.tmp <- subset(tbl.tmp, model != "glmnet")
  tbl.tmp[(tbl.tmp$model == "rf"),"model"] <- "RF"
  tbl.tmp[(tbl.tmp$model == "svm"),"model"] <- "SVM"
  tbl.tmp$model <- factor(tbl.tmp$model, levels = c("LASSO", "Ridge", "RF", "SVM"))
  colnames(tbl.tmp)[colnames(tbl.tmp) == "model"] <- "Model"
  tbl.tmp$cor <- as.numeric(tbl.tmp$cor)
  tbl.tmp$ohsu.drug <- factor(tbl.tmp$ohsu.drug)
  tbl.tmp
}

## Plot IC50, DSS2, and AUC for Ridge RNA
t1 <- tbls[["ic50.l4"]][["rna"]]
t1$Metric <- "IC50"
t2 <- tbls[["dss2.l4"]][["rna"]]
t2$Metric <- "DSS2"
t3 <- tbls[["auc.l4"]][["rna"]]
t3$Metric <- "AUC"
tbl <- rbind(t1, t2, t3)
tbl <- sanitize.table(tbl)
tbl <- subset(tbl, Model == "Ridge")
tbl$Metric <- factor(tbl$Metric, levels = c("IC50", "AUC", "DSS2"))

g <- ggplot(data = tbl, aes(x = Metric, y = cor, colour = Metric))
g <- g + geom_violin()
g <- g + geom_beeswarm()
g <- g + xlab("Metric")
g <- g + ylab("Correlation")
g <- g + theme(axis.text.x=element_text(angle = 45, hjust = 1))
tiff("rna-ridge-vs-metric.tiff")
print(g)
d <- dev.off()

## Plot RNA, DNA, and both for Ridge DSS2
t1 <- tbls[["dss2.l4"]][["rna"]]
t1$Data <- "RNA"
t2 <- tbls[["dss2.l4"]][["dna"]]
t2$Data <- "DNA"
t3 <- tbls[["dss2.l4"]][["both"]]
t3$Data <- "Both"
tbl <- rbind(t1, t2, t3)
tbl <- sanitize.table(tbl)
tbl <- subset(tbl, Model == "Ridge")
tbl$Data <- factor(tbl$Data, levels = c("DNA", "RNA", "Both"))

g <- ggplot(data = tbl, aes(x = Data, y = cor, colour = Data))
g <- g + geom_violin()
g <- g + geom_beeswarm()
g <- g + xlab("Data")
g <- g + ylab("Correlation")
g <- g + theme(axis.text.x=element_text(angle = 45, hjust = 1))
tiff("dss2-ridge-vs-data.tiff")
print(g)
d <- dev.off()

## Plot LASSO, RF, Ridge, and SVM for DSS2 RNA
tbl <- tbls[["dss2.l4"]][["rna"]]
tbl <- sanitize.table(tbl)

g <- ggplot(data = tbl, aes(x = Model, y = cor, colour = Model))
g <- g + geom_violin()
g <- g + geom_beeswarm()
g <- g + xlab("Model")
g <- g + ylab("Correlation")
g <- g + theme(axis.text.x=element_text(angle = 45, hjust = 1))
tiff("dss2-rna-vs-model.tiff")
print(g)
d <- dev.off()

## Plot DSS2 RNA LASSO, RF, Ridge, and SVM vs drugs
tbl <- tbls[["dss2.l4"]][["rna"]]
tbl <- sanitize.table(tbl)
rftbl <- subset(tbl, Model == "RF")
rftbl <- rftbl[order(rftbl$cor, decreasing=TRUE),]
ordered.levels <- unique(rftbl$ohsu.drug)
tbl$ohsu.drug <- factor(tbl$ohsu.drug, levels = ordered.levels)
g <- ggplot(data = tbl, aes(x = ohsu.drug, y = cor, colour = Model, group = Model))
g <- g + geom_point() + geom_line()
g <- g + xlab("Drug")
g <- g + ylab("Correlation")
g <- g + theme(axis.text.x=element_text(angle = 45, hjust = 1, size = 6))
tiff("dss2-rna-vs-drugs.tiff")
print(g)
d <- dev.off()

save.image(".Rdata")

cat("Finished plotting successfully\n")

q(status=0)

