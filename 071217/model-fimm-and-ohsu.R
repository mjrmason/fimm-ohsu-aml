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

## Read in the FIMM expression data (RNASeq.CPM.log2.bc.csv)
## Load (batch corrected) FIMM expression data
obj <- synGet(id="syn8270602", downloadFile = TRUE)
file <- getFileLocation(obj)
fimm.expr <- read.table(file, header=TRUE, sep=",")
fimm.expr <- t(fimm.expr)

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
ohsu.expr.arg <- ohsu.expr
fimm.expr.arg <- fimm.expr
response.col <- "dss2.l4"

ohsu.query.drugs <- ohsu.common.drugs
fimm.query.drugs <- fimm.common.drugs

## ohsu.query.drugs <- ohsu.mek.drugs
## fimm.query.drugs <- fimm.mek.drugs

rets <- list()
for(response.col in c("dss2.l4", "auc.l4", "ic50.l4")) {
  cat(paste0("Fitting ", response.col, "\n"))
  rets[[response.col]] <- model.fimm.and.ohsu(ohsu.dss.arg, ohsu.expr.arg, ohsu.common.drugs, ohsu.query.drugs, fimm.dss.arg, fimm.expr.arg, fimm.common.drugs, fimm.query.drugs, response.col)
  save.image(".Rdata")
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
  tbls[[response.col]] <-
     ldply(rets[[response.col]], 
             .fun = function(foo) {
                ldply(foo,
             .fun = function(df) {
                ct <- cor.test(df$predicted, df$actual)
                data.frame(cor = ct$estimate, pval = ct$p.value, ohsu.drug = df$ohsu.drug, fimm.drug = df$fimm.drug, alpha = df$alpha, s = df$s, model = df$model)
             }) })
  metric <- response.col
  metric <- gsub(metric, pattern=".l4", replacement="")
  metric <- toupper(metric)

  tbl <- tbls[[response.col]]
  t2 <- subset(tbl, is.na(s) | (s == "lambda.min"))
  t2$model <- as.character(t2$model)
  t2[(t2$model == "glmnet") & (t2$alpha == 0),"model"] <- "ridge"
  t2[(t2$model == "glmnet") & (t2$alpha == 1),"model"] <- "lasso"
  t2 <- subset(t2, model != "glmnet")

  g <- plot.correlation.vs.drug(t2)
  g <- g + ggtitle(metric)
  pdf(paste0(metric, "-correlation-vs-drug.pdf"))
  print(g)
  d <- dev.off()

  mydoc <- addSlide( mydoc, slide.layout = "Title and Content" )
  title <- paste0(metric, ": Correlation vs Drug")
  mydoc <- addTitle( mydoc, title )
  mydoc <- addPlot( doc = mydoc, fun = print, x = g)

  g <- plot.correlation.vs.model(t2)
  g <- g + ggtitle(metric)
  pdf(paste0(metric, "-correlation-vs-model.pdf"))
  print(g)
  d <- dev.off()

  mydoc <- addSlide( mydoc, slide.layout = "Title and Content" )
  title <- paste0(metric, ": Correlation vs Model")
  mydoc <- addTitle( mydoc, title )
  mydoc <- addPlot( doc = mydoc, fun = print, x = g)

}

writeDoc(mydoc, doc.file)

cat("Finished successfully\n")
q(status=0)

