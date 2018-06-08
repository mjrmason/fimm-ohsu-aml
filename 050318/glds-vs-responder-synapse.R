rm(list = ls())

suppressPackageStartupMessages(library(glmnet))
suppressPackageStartupMessages(library("randomForest"))
suppressPackageStartupMessages(library("e1071"))
suppressPackageStartupMessages(library("foreach"))
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(tidyr))

suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("gridExtra"))
suppressPackageStartupMessages(library("caTools"))
suppressPackageStartupMessages(library("data.table"))

suppressPackageStartupMessages(library(ggbeeswarm))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(synapseClient))
suppressPackageStartupMessages(library(mygene))

## Log in to Synapse
synapseLogin()

## Bring in common code
source("../common/data-preprocessing.R")
source("../common/models.R")
source("../common/plotting.R")
source("../common/utils.R")

## Register parallelization
num.cores <- detectCores()
if(!is.na(num.cores) && (num.cores > 1)) {
  suppressPackageStartupMessages(library("doMC"))
  cat(paste("Registering ", num.cores-1, " cores.\n", sep=""))
  registerDoMC(cores=(num.cores-1))
}
num.processes <- num.cores - 1

## Begin setup (this will need to be changed as the data sets in the intersection change)

source("fimm-ohsu-setup-040218.R")

## Turn _off_ gene subsetting
use.subset <- FALSE

source("process-drc-and-expr.R")

synId <- "syn8149174"
obj <- synGet(id=synId, downloadFile = TRUE)
file <- getFileLocation(obj)
tmp <- read.xlsx(file, sheet = 1)
ohsu.metadata <- tmp

tmp <- read.xlsx(file, sheet = 3)
tmp <- tmp[grepl(pattern="FAB", ignore.case=TRUE, tmp$lab_type),]
tmp <- unique(tmp[, c("patient_id", "lab_result")])
colnames(tmp) <- c("patient_id", "fab")
ohsu.metadata <- merge(ohsu.metadata, tmp, by = "patient_id", all = TRUE)

patient.age.tbl <- read.xlsx(file, sheet = 5)
patient.age.tbl <- patient.age.tbl[!(!is.na(patient.age.tbl$age) & (as.numeric(patient.age.tbl$age) < 0)),]

## Map lab_id to seq_ids
synId <- "syn10083817"
obj <- synGet(synId, downloadFile = TRUE)
file <- getFileLocation(obj)
ohsu.rnaseq.sample.summary <- read.table(file, header=TRUE, sep="\t")
ohsu.metadata <- merge(ohsu.metadata, ohsu.rnaseq.sample.summary[, c("Original_LabID", "SeqID")], by.x = "lab_id", by.y = "Original_LabID", all = FALSE)
## patient.age.tbl <- merge(patient.age.tbl, ohsu.rnaseq.sample.summary[, c("Original_LabID", "SeqID")], by.x = "lab_id", by.y = "Original_LabID", all = FALSE)
ohsu.metadata <- merge(ohsu.metadata, patient.age.tbl, by = c("lab_id", "patient_id"))

write.table(file = "ohsu-patient-metadata.tsv", ohsu.metadata, sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(file = "ohsu-drug-response.tsv", orig.drcs[["ohsu"]], sep="\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
write.table(file = "ohsu-drug-metadata.tsv", drug.map, sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)



stop("stop")


## Find patients that have only one type of response across all samples (e.g., _not_ "Complete Response" and "Refractory")
patient.response.tbl <- patient.age.tbl
un <- unique(patient.response.tbl[, c("patient_id", "response")])
un <- un[!my.dup(un$patient_id),]
un <- subset(un, !is.na(response))
un <- subset(un, response %in% c("Complete Response", "Refractory"))

seq.response <- un
rownames(seq.response) <- un$SeqID
rownames(seq.response) <- un$patient_id

common.samples <- intersect(rownames(seq.response), colnames(orig.drcs[["ohsu"]]))
ohsu.glds <- colMeans(t(scale(t(orig.drcs[["ohsu"]][, common.samples]))), na.rm=TRUE)
seq.response <- seq.response[common.samples, ]
seq.response$glds <- as.numeric(ohsu.glds)

lm <- lm(glds ~ response, data = seq.response)
sm <- summary(lm)
pval <- round(1 - pf(sm$fstatistic[1], sm$fstatistic[2], sm$fstatistic[3]), digits = 3)
wt <- wilcox.test(glds ~ response, data = seq.response)
pval <- wt$p.value

g <- ggplot(seq.response, aes(x = response, y = glds))
g <- g + geom_violin()
g <- g + geom_beeswarm()
g <- g + ylab("GLDS")
g <- g + ggtitle(paste0("Response ~ GLDS: p = ", as.numeric(pval)))
pdf("ohsu-glds-vs-response.pdf")
g
d <- dev.off()

