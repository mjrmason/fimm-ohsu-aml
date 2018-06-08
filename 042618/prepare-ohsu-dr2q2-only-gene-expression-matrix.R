# Do this for read.xlsx2
options(java.parameters = "- Xmx1024m")


cat("Script to prepare the OHSU gene expression matrix for Data Release 2 Q2\n")

suppressPackageStartupMessages(library("synapseClient"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("foreach"))
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("gplots"))

num.cores <- detectCores()
if(!is.na(num.cores) && (num.cores > 1)) {
  suppressPackageStartupMessages(library("doMC"))
  cat(paste("Registering ", num.cores-1, " cores.\n", sep=""))
  registerDoMC(cores=(num.cores-1))
}

## Set the download directory
download.path <- "../input/"
if (!file.exists(download.path)) {
  dir.create(download.path)
}

synapseLogin()

suppressPackageStartupMessages(library("openxlsx"))

## This file has a map of patient_id to lab_id
cat("Downloading diagnosis file\n")
## Diagnosis_Labs_Treatments_Outcomes_2018_01_12.xlsx
synId <- "syn11766999"
obj <- synGet(id=synId, downloadFile = TRUE, downloadLocation = download.path)
diagnosis.tbl <- read.xlsx(getFileLocation(obj), sheet=1)
aml.diagnosis <- "ACUTE MYELOID LEUKAEMIA (AML) AND RELATED PRECURSOR NEOPLASMS"
aml.diagnosis.tbl <- subset(diagnosis.tbl, most_recent_diagnosis %in% aml.diagnosis)
aml.diagnosis.tbl <- subset(aml.diagnosis.tbl, diagnosis_at_time_of_specimen_acquisition %in% aml.diagnosis)
cat("AML patient samples\n")
patient.diagnosis.tbl <- diagnosis.tbl[,c("patient_id", "lab_id", "most_recent_diagnosis", "diagnosis_at_time_of_specimen_acquisition")]

outlier.sample.ids <- c("13-00330", "13-00636", "14-00065", "14-00494", "13-00537", "15-00309", "20-00063", 
                        "20-00071", "20-00343", "20-00090", "20-00096", "20-00314", "20-00321")


## Load OHSU expression data
load.ohsu.expr.data <- function() {
  ## Download the AML CPM expression data
  ## NB: linking to the original Beat AML project (syn2942337), since this has more updated files/releases
  ## (as well as the raw fastqs)
  ## Read in BeatAML_DR1_RNASeq_log2_cpm_2016_08_02.csv
  cat("Downloading and reading DR1 data\n")
  rna.dr1.obj <- synGet("syn7124177", downloadFile=TRUE)
  # Load the data
  rna.dr1.tbl <- as.data.frame(fread(getFileLocation(rna.dr1.obj), sep=",", header=TRUE))
  rna.dr1.seq.ids <- gsub("\\.", "-", gsub("^X", "", names(rna.dr1.tbl[-c(1:9)])))

  cat("Downloading and reading DR2 Q1 data\n")
  ## Read in BeatAML_DR2q1_RNASeq_log2_cpm_2016_12_16.csv
  rna.dr2q1.obj <- synGet("syn8149220", downloadFile=TRUE)
  rna.dr2q1.tbl <- as.data.frame(fread(getFileLocation(rna.dr2q1.obj), sep=",", header=TRUE))
  rna.dr2q1.seq.ids <- gsub("\\.", "-", gsub("^X", "", names(rna.dr2q1.tbl[-c(1:9)])))

  cat("Downloading and reading DR2 Q2 data\n")
  ## Read in BeatAML_DR2q2_RNASeq_log2_cpm_2017_01_11.csv
  rna.dr2q2.obj <- synGet("syn8149259", downloadFile=TRUE)
  rna.dr2q2.tbl <- as.data.frame(fread(getFileLocation(rna.dr2q2.obj), sep=",", header=TRUE))
  rna.dr2q2.seq.ids <- gsub("\\.", "-", gsub("^X", "", names(rna.dr2q2.tbl[-c(1:9)])))

  cat("Downloading and reading DR2 rnaseq data\n")
  ## Read in BeatAML_DR2_RNASeq_log2_cpm_2017_06_08.csv
  rna.dr2.obj <- synGet("syn10008812", downloadFile=TRUE)
  rna.dr2.tbl <- as.data.frame(fread(getFileLocation(rna.dr2.obj), sep=",", header=TRUE))
  rna.dr2.seq.ids <- gsub("\\.", "-", gsub("^X", "", names(rna.dr2.tbl[-c(1:9)])))

  cat("Downloading and reading DR3 rnaseq data\n")
  ## Read in BeatAML_DR3q1_RNASeq_log2_cpm_2018_01_26.csv
  rna.dr3q1.obj <- synGet("syn11767059", downloadFile=TRUE, downloadLocation = download.path)
  rna.dr3q1.tbl <- as.data.frame(fread(getFileLocation(rna.dr3q1.obj), sep=",", header=TRUE))
  rna.dr3q1.seq.ids <- gsub("\\.", "-", gsub("^X", "", names(rna.dr3q1.tbl[-c(1:9)])))

  ## Make sure we have a common set of genes across all 3 expression data sets
  lst <- list(rna.dr1.tbl$Gene, rna.dr2q1.tbl$Gene, rna.dr2q2.tbl$Gene, rna.dr2.tbl$Gene, rna.dr3q1.tbl$Gene)
  common.genes <- Reduce(intersect, lst)
  rownames(rna.dr1.tbl) <- rna.dr1.tbl$Gene
  rownames(rna.dr2q1.tbl) <- rna.dr2q1.tbl$Gene
  rownames(rna.dr2q2.tbl) <- rna.dr2q2.tbl$Gene
  rownames(rna.dr2.tbl) <- rna.dr2.tbl$Gene
  rownames(rna.dr3q1.tbl) <- rna.dr3q1.tbl$Gene
  rna.dr1.tbl <- rna.dr1.tbl[common.genes,]
  rna.dr2q1.tbl <- rna.dr2q1.tbl[common.genes,]
  rna.dr2q2.tbl <- rna.dr2q2.tbl[common.genes,]
  rna.dr2.tbl <- rna.dr2.tbl[common.genes,]
  rna.dr3q1.tbl <- rna.dr3q1.tbl[common.genes,]

  ## Drop the annotation columns from the expression matrices
  rna.dr1.tbl <- rna.dr1.tbl[,-c(1:9)]
  rna.dr2q1.tbl <- rna.dr2q1.tbl[,-c(1:9)]
  rna.dr2q2.tbl <- rna.dr2q2.tbl[,-c(1:9)]
  rna.dr2.tbl <- rna.dr2.tbl[,-c(1:9)]
  rna.dr3q1.tbl <- rna.dr3q1.tbl[,-c(1:9)]
  gc()

  ## Ensure that DR2 Q1 and Q2 are non-overlapping.
  dr2q1.and.2.samples <- c(colnames(rna.dr2q1.tbl), colnames(rna.dr2q2.tbl))
  if(any(duplicated(dr2q1.and.2.samples))) {
    stop("Duplicated samples between DR2 Q1 and Q2!\n")
  } else {
    cat("\nAs expected, no duplicated samples between DR2 Q1 and Q2\n")
  }

  ## DR2 is not a strict superset of DR2 Q1 and Q2:
  missing.samples <- sort(dr2q1.and.2.samples[!(dr2q1.and.2.samples %in% colnames(rna.dr2.tbl))])
  cat(paste0(length(missing.samples), " of ", length(dr2q1.and.2.samples), " samples from DR2 Q1 and DR2 Q2 are not in DR2:\n"))
  cat(paste0(paste(missing.samples, collapse=","), "\n"))

  ## Presumably those few missing from DR2 were excluded for a reason.  
  ## As described in the README (datawave-2-q3_rnaseq_README.pdf):
  ## "Note that the exact log2(cpm) values will differ if samples are added or removed 
  ## due to the scaling of the prior.count value."
  ## Hence, the values in DR2 Q1 and Q2 for sample X is slightly different than for that
  ## same sample in DR2.  Let's just use DR2.

  ## Also, note the following under section Genotyping Quality Control in datawave-2-q3_rnaseq_README.pdf
  ## "NOTE: Three samples released in earlier Data Wave 2 quarterly updates did not pair as expected and have
  ## been removed from Data Lock 2: SeqID 20-00063, 20-00071, and 20-00343. Four samples released in earlier
  ## Data Wave 2 quarterly updates have been found to be ‘outliers’: 20-00090, 20-00096, 20-00314, 20-00321.

  ## Three samples from Data Lock 1 did not pair as expected in subsequent quality control assessment: SeqID
  ## 13-00330, 13-00636, 14-00065, and 14-00494; we recommend not using these samples. Two samples were
  ## found to be ‘outliers’ during subsequent QC assessment: SeqID 13-00537 and 15-00309; use of these files in
  ## standard analyses in not recommended."
  outlier.sample.ids <- sort(outlier.sample.ids)
  cat(paste0("The following samples were excluded from Data Wave 2 because of QC issues:\n", paste(outlier.sample.ids, collapse=","), "\n"))

  ## Ensure that there is no overlap between DR1, DR2, and DR3.
  lst <- list(colnames(rna.dr1.tbl), colnames(rna.dr2.tbl), colnames(rna.dr3q1.tbl))
  all.samples <- Reduce(c, lst)
  if(any(duplicated(all.samples))) {
    warning("\nSome RNA-seq samples are duplicated across the 3 releases\n")
  } else {
    cat("\nAs expected, no samples are duplicated across the 3 releases\n")
  }
  ret <- list()
  ret[["rna.dr1"]] <- rna.dr1.tbl
  ret[["rna.dr2q1"]] <- rna.dr2q1.tbl
  ret[["rna.dr2q2"]] <- rna.dr2q2.tbl
  ret[["rna.dr2"]] <- rna.dr2.tbl
  ret[["rna.dr3q1"]] <- rna.dr3q1.tbl
  return(ret)

  if(FALSE) {
  rm(rna.dr2q1.tbl)
  rm(rna.dr2q2.tbl)
  gc()

  ## Since no samples are duplicated, we can just merge the columns
  ohsu.expr <- cbind(rna.dr1.tbl, rna.dr2.tbl, rna.dr3q1.tbl)
  rm(rna.dr1.tbl)
  rm(rna.dr2.tbl)
  rm(rna.dr3.q1.tbl)
  gc()
  ohsu.expr
  }
}

ohsu.exprs <- load.ohsu.expr.data()

ohsu.dr2q2.only.expr <- ohsu.exprs[["rna.dr2q2"]]
## Exclude any outliers
ohsu.dr2q2.only.expr <-	ohsu.dr2q2.only.expr[, !(colnames(ohsu.dr2q2.only.expr) %in% outlier.sample.ids)]

## Exclude any samples in DR1 or DR2 Q1
for(ds in c("rna.dr1", "rna.dr2q1")) {
  ohsu.dr2q2.only.expr <-	ohsu.dr2q2.only.expr[, !(colnames(ohsu.dr2q2.only.expr) %in% colnames(ohsu.exprs[[ds]]))]
}

## Exclude any _patients_ in DR1 or DR2 Q1
## Read in the latest rna-seq dashboard
## BeatAML_rnaseq_2018_01_29_public_dashboard.xlsx (syn11767026)
## As stated in the README (datawave-2-q3_rnaseq_README.pdf),
## this is a cumulative dashboard that covers all current releases (i.e., DR1, DR2, and DR3)
rna.obj <- synGet("syn11767026", downloadFile=TRUE)
file <- getFileLocation(rna.obj)

## The samples summary is the 4th sheet
## This sample.summary table provides a map from patient id to sequence id.
## Expect that since the sample sequence ids where column headers in our expression
## data, they were corrupted from AA-BBBBB (where A and B are digits) to XAA.BBBBB.
## So change that.
cat(paste0("Reading sample summary file ", file, "\n"))
rnaseq.sample.summary.tbl <- read.xlsx(file, sheet=4)
rnaseq.sample.summary.tbl$Original_LabID <- as.character(rnaseq.sample.summary.tbl$Original_LabID)
rnaseq.sample.summary.tbl$SeqID <- as.character(rnaseq.sample.summary.tbl$SeqID)

library(plyr)
library(dplyr)
patients <- llply(ohsu.exprs, 
                  .fun = function(mat) {
                           tmp <- subset(rnaseq.sample.summary.tbl, !is.na(SeqID) & (SeqID %in% colnames(mat)))
                           return(unique(tmp[, c("SeqID", "PatientID")]))
                         })

for(ds in c("rna.dr1", "rna.dr2q1")) {
  flag <- patients[["rna.dr2q2"]]$PatientID %in% patients[[ds]]$PatientID 
  if(any(flag)) {
    cat(paste0("Dropping patients: ", paste(patients[["rna.dr2q2"]]$PatientID[flag], collapse=","), "\n"))
    cat(paste0("Num RNA-seq patients before dropping: ", ncol(ohsu.dr2q2.only.expr), "\n"))
    ohsu.dr2q2.only.expr <- ohsu.dr2q2.only.expr[, !(colnames(ohsu.dr2q2.only.expr) %in% patients[["rna.dr2q2"]]$SeqID[flag])]
    cat(paste0("Num RNA-seq patients after dropping: ", ncol(ohsu.dr2q2.only.expr), "\n"))
  }
}

path <- "ohsu.dr2q2.only.expr.tsv"
parentId <- "syn10083332"
write.table(file=path, ohsu.dr2q2.only.expr, row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)
f <- File(path, parentId = parentId, synapseStore = TRUE)
synStore(f, executed = NULL, forceVersion = FALSE)

