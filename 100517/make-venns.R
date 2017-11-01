# Do this for read.xlsx2
options(java.parameters = "- Xmx1024m")


cat("Script to get the intersection of drug response and expression data for OHSU, FIMM, and CTRP\n")
cat("Individually and in aggregrate\n")

suppressPackageStartupMessages(library("synapseClient"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("foreach"))
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library("data.table"))

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

beat.aml.cache <- "../input/beat-aml-files/"

synapseLogin()

suppressPackageStartupMessages(library("openxlsx"))

## Download the 3 expression matrices.
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

## Download the 3 drug data sets.
## path <- "ohsu.dss.t0.tsv"
## synId <- "syn10083494"

## path <- "fimm.dss.t0.tsv"
## synId <- "syn10083488"

## This is ctrp.l4.and.ll4.fits.non.censored.pts.non.censored.curves.tsv
## synId <- "syn10307313"

drug.synIds <- c("syn10083494", "syn10083488", "syn10307313")
names(drug.synIds) <- c("ohsu", "fimm", "ctrp")

drug.matrices <- list()

for(i in 1:length(drug.synIds)) {
  synId <- drug.synIds[i]
  obj <- synGet(synId, downloadFile = TRUE)
  file <- getFileLocation(obj)
  drug.matrices[[names(drug.synIds)[i]]] <- as.data.frame(fread(file))
}

library(gplots)

## Create OHSU venn.
## Add the SeqID to the OHSU drug data.
## Read in the rna-seq summary table, which provides a map from the lab_id's used
## in the drug response data to the SeqID's used in the expression data.
synId <- "syn10083817"
obj <- synGet(synId, downloadFile = TRUE)
file <- getFileLocation(obj)
ohsu.rnaseq.sample.summary <- read.table(file, header=TRUE, sep="\t")
drug.matrices[["ohsu"]] <- merge(drug.matrices[["ohsu"]], ohsu.rnaseq.sample.summary, by.x = "lab_id", by.y = "Original_LabID")

drug <- unique(drug.matrices[["ohsu"]][, "SeqID"])
ohsu.drugs <- data.frame(ohsu.drug = unique(drug.matrices[["ohsu"]][, "inhibitor"]))
expr <- unique(colnames(expr.matrices[["ohsu"]]))
vennList = list("Drug"=drug, "RNA-seq" = expr)
pdf("ohsu-drug-expr-venn.pdf")
venn(vennList)
title(main="BEAT AML Patient Overlap")
d <- dev.off()

## Create FIMM venn
drug <- unique(drug.matrices[["fimm"]][, "SCREEN_ID"])
fimm.drugs <- data.frame(fimm.drug = unique(drug.matrices[["fimm"]][, "DRUG_ID"]))
expr <- unique(colnames(expr.matrices[["fimm"]]))
vennList = list("Drug"=drug, "RNA-seq" = expr)
pdf("fimm-drug-expr-venn.pdf")
venn(vennList)
title(main="FIMM Patient Overlap")
d <- dev.off()

## Create CTRP venn

synId <- "syn5632192"
obj <- synGet(synId, downloadFile = TRUE)
ctrp.cell.line.metadata <- read.table(getFileLocation(obj), sep="\t", header=TRUE, as.is=TRUE, comment.char="", quote="")

all.ccls <- unique(drug.matrices[["ctrp"]][, "ccl_name"])
heme.ccls <- all.ccls[all.ccls %in% ctrp.cell.line.metadata[ctrp.cell.line.metadata$ccle_primary_hist == "haematopoietic_neoplasm", "ccl_name"]]
aml.ccls <- all.ccls[all.ccls %in% ctrp.cell.line.metadata[ctrp.cell.line.metadata$ccle_hist_subtype_1 == "acute_myeloid_leukaemia", "ccl_name"]]
flag <- drug.matrices[["ctrp"]][, "ccl_name"] %in% aml.ccls
drug.matrices[["ctrp"]] <- drug.matrices[["ctrp"]][flag, ]

drug <- unique(drug.matrices[["ctrp"]][, "ccl_name"])
expr <- unique(colnames(expr.matrices[["ctrp"]]))
vennList = list("Drug"=drug, "RNA-seq" = expr)
pdf("ctrp-aml-drug-expr-venn.pdf")
venn(vennList)
title(main="CTRP AML Patient Overlap")
d <- dev.off()


## What samples are in common across the expression and drug response data in each data set.

## Load the map of drugs shared between OHSU and FIMM (FIMM_OHSU_Drug_Concentrations.xlsx)
synId <- "syn10163669"
obj <- synGet(synId, downloadFile = TRUE)
ohsu.fimm.drugs <- openxlsx:::read.xlsx(getFileLocation(obj), sheet=1)
ohsu.fimm.drugs <- ohsu.fimm.drugs[, c("OHSU_DRUG_NAME", "FIMM_DRUG_NAME", "Non.proprietary.name", "Trade.names", "Aliases", "FIMM_Batch_ID_Drug")]
ohsu.fimm.drugs <- merge(ohsu.fimm.drugs, fimm.drugs, by.x = "FIMM_Batch_ID_Drug", by.y = "fimm.drug", all = TRUE)
ohsu.fimm.drugs <- merge(ohsu.fimm.drugs, ohsu.drugs, by.x = "OHSU_DRUG_NAME", by.y = "ohsu.drug", all = TRUE)

## Load the CTRP drugs
## CTRP folder on Synapse: syn5622707
## CTRP drug info: "v20.meta.per_compound.txt" (syn5632193)
synId <- "syn5632193"
obj <- synGet(synId, downloadFile = TRUE)
ctrp.compounds <- read.table(getFileLocation(obj), sep="\t", header=TRUE, as.is=TRUE, comment.char="", quote="")

## CTRP experiment metadata (including DMSO): "v20.meta.per_experiment.txt" (syn5632194)
synId <- "syn5632194"
obj <- synGet(synId, downloadFile = TRUE)
ctrp.expt.metadata <- read.table(getFileLocation(obj), sep="\t", header=TRUE, as.is=TRUE, comment.char="", quote="")

## Match CTRP drugs to the drugs in common between FIMM and OHSU
ohsu.fimm.drugs$CTRP_DRUG_NAME <- NA
for(i in 1:nrow(ohsu.fimm.drugs)) {
  drugs <- unique(unname(unlist(lapply(ohsu.fimm.drugs[i,c("OHSU_DRUG_NAME", "FIMM_DRUG_NAME", "Non.proprietary.name", "Trade.names", "Aliases")], function(str) unlist(strsplit(unname(str), split=",[ ]*"))))))
  drugs <- drugs[!is.na(drugs)]
  flag <- grepl(pattern="\\(", drugs)
  if(any(flag)) {
    drugs <- c(drugs, unlist(lapply(drugs[flag], function(str) c(gsub(".*\\((.*)\\).*", "\\1", str), gsub("(.*) \\(.*\\).*", "\\1", str)))))
  }
  matches <- c()
  for(drug in tolower(drugs)) {
    if(any(tolower(ctrp.compounds$cpd_name) == drug)) {
      matches <- c(matches, which(tolower(ctrp.compounds$cpd_name) == drug))
    }
  }
  matches <- unique(matches)
  if(length(matches) > 1) { stop("Got multiple matches to CTRP!\n") }
  if(length(matches) > 0) {
    ohsu.fimm.drugs$CTRP_DRUG_NAME[i] <- ctrp.compounds$cpd_name[matches[1]]
  }
}

## What drugs overlap?

ohsu.fimm.drugs$id <- 1:nrow(ohsu.fimm.drugs)
uniq.ohsu.drugs <- ohsu.fimm.drugs$id[!is.na(ohsu.fimm.drugs$OHSU_DRUG_NAME)]
uniq.fimm.drugs <- ohsu.fimm.drugs$id[!is.na(ohsu.fimm.drugs$FIMM_Batch_ID_Drug)]
uniq.ctrp.drugs <- ohsu.fimm.drugs$id[!is.na(ohsu.fimm.drugs$CTRP_DRUG_NAME)]

vennList = list("OHSU"=uniq.ohsu.drugs, "FIMM"=uniq.fimm.drugs, "CTRP"=uniq.ctrp.drugs)
pdf("ohsu-fimm-ctrp-drug-venn.pdf")
venn(vennList)
title(main="OHSU, FIMM, CTRP AML Drug Overlap")
d <- dev.off()


## Heatmap of drug vs expression.

## This file has a map of patient_id to lab_id
cat("Downloading diagnosis file\n")
## Diagnosis_Labs_Treatments_Outcomes_2017_05_17.xlsx
obj <- synGet(id="syn10008777", downloadFile = TRUE, downloadLocation = download.path)
diagnosis.tbl <- read.xlsx(getFileLocation(obj), sheet=1)
aml.diagnosis <- "ACUTE MYELOID LEUKAEMIA (AML) AND RELATED PRECURSOR NEOPLASMS"
aml.diagnosis.tbl <- subset(diagnosis.tbl, most_recent_diagnosis %in% aml.diagnosis)
aml.diagnosis.tbl <- subset(aml.diagnosis.tbl, diagnosis_at_time_of_specimen_acquisition %in% aml.diagnosis)
cat("AML patient samples\n")
patient.diagnosis.tbl <- diagnosis.tbl[,c("patient_id", "lab_id", "most_recent_diagnosis", "diagnosis_at_time_of_specimen_acquisition")]

## Load OHSU expression data
load.ohsu.expr.data <- function() {
  ## Download the AML CPM expression data
  ## NB: linking to the original Beat AML project (syn2942337), since this has more updated files/releases
  ## (as well as the raw fastqs)
  ## Read in BeatAML_DR1_RNASeq_log2_cpm_2016_08_02.csv
  cat("Downloading and reading DR1 data\n")
  rna.dr1.obj <- synGet("syn7124177", downloadFile=TRUE, downloadLocation = download.path)
  # Load the data
  rna.dr1.tbl <- as.data.frame(fread(getFileLocation(rna.dr1.obj), sep=",", header=TRUE))
  rna.dr1.lab.ids <- gsub("\\.", "-", gsub("^X", "", names(rna.dr1.tbl[-c(1:9)])))

  cat("Downloading and reading DR2 Q1 data\n")
  ## Read in BeatAML_DR2q1_RNASeq_log2_cpm_2016_12_16.csv
  rna.dr2q1.obj <- synGet("syn8149220", downloadFile=TRUE, downloadLocation = download.path)
  rna.dr2q1.tbl <- as.data.frame(fread(getFileLocation(rna.dr2q1.obj), sep=",", header=TRUE))

  rna.dr2q1.lab.ids <- gsub("\\.", "-", gsub("^X", "", names(rna.dr2q1.tbl[-c(1:9)])))

  cat("Downloading and reading DR2 Q2 data\n")
  ## Read in BeatAML_DR2q2_RNASeq_log2_cpm_2017_01_11.csv
  rna.dr2q2.obj <- synGet("syn8149259", downloadFile=TRUE, downloadLocation = download.path)
  rna.dr2q2.tbl <- as.data.frame(fread(getFileLocation(rna.dr2q2.obj), sep=",", header=TRUE))
  rna.dr2q2.lab.ids <- gsub("\\.", "-", gsub("^X", "", names(rna.dr2q2.tbl[-c(1:9)])))

  cat("Downloading and reading DR2 rnaseq data\n")
  rna.dr2.obj <- synGet("syn10008812", downloadFile=TRUE, downloadLocation = download.path)
  rna.dr2.tbl <- as.data.frame(fread(getFileLocation(rna.dr2.obj), sep=",", header=TRUE))
  rna.dr2.lab.ids <- gsub("\\.", "-", gsub("^X", "", names(rna.dr2.tbl[-c(1:9)])))

  ## Make sure we have a common set of genes across all 3 expression data sets
  common.genes <- intersect(rna.dr1.tbl$Gene, intersect(rna.dr2q1.tbl$Gene, intersect(rna.dr2q2.tbl$Gene, rna.dr2.tbl$Gene)))
  rownames(rna.dr1.tbl) <- rna.dr1.tbl$Gene
  rownames(rna.dr2q1.tbl) <- rna.dr2q1.tbl$Gene
  rownames(rna.dr2q2.tbl) <- rna.dr2q2.tbl$Gene
  rownames(rna.dr2.tbl) <- rna.dr2.tbl$Gene
  rna.dr1.tbl <- rna.dr1.tbl[common.genes,]
  rna.dr2q1.tbl <- rna.dr2q1.tbl[common.genes,]
  rna.dr2q2.tbl <- rna.dr2q2.tbl[common.genes,]
  rna.dr2.tbl <- rna.dr2.tbl[common.genes,]

  ## Drop the annotation columns from the expression matrices
  rna.dr1.tbl <- rna.dr1.tbl[,-c(1:9)]
  rna.dr2q1.tbl <- rna.dr2q1.tbl[,-c(1:9)]
  rna.dr2q2.tbl <- rna.dr2q2.tbl[,-c(1:9)]
  rna.dr2.tbl <- rna.dr2.tbl[,-c(1:9)]
  gc()

  ## Ensure that DR2 Q1 and Q2 are non-overlapping.
  dr2q1.and.2.samples <- c(colnames(rna.dr2q1.tbl), colnames(rna.dr2q2.tbl))
  if(any(duplicated(dr2q1.and.2.samples))) {
    stop("Duplicated samples between DR2 Q1 and Q2!\n")
  } else {
    cat("As expected, no duplicated samples between DR2 Q1 and Q2\n")
  }

  ## DR2 is not a strict superset of DR2 Q1 and Q2:
  cat(paste0(length(which(!(dr2q1.and.2.samples %in% colnames(rna.dr2.tbl)))), " of ", length(dr2q1.and.2.samples), " samples from DR2 Q1 and DR2 Q2 are not in DR2\n"))
  ## Presumably those few missing from DR2 were excluded for a reason.  
  ## As described in the README (datawave-2-q3_rnaseq_README.pdf):
  ## "Note that the exact log2(cpm) values will differ if samples are added or removed 
  ## due to the scaling of the prior.count value."
  ## Hence, the values in DR2 Q1 and Q2 for sample X is slightly different than for that
  ## same sample in DR2.  Let's just use DR2.

  ## Ensure that there is no overlap between DR1 and DR2.
  rm(rna.dr2q1.tbl)
  rm(rna.dr2q2.tbl)
  gc()

  all.samples <- c(colnames(rna.dr1.tbl), colnames(rna.dr2.tbl))
  if(any(duplicated(all.samples))) {
    warning("Some RNA-seq samples are duplicated across the 2 releases\n")
  } else {
    cat("As expected, no samples are duplicated across the 2 releases\n")
  }
  ## Since no samples are duplicated, we can just merge the columns
  ohsu.expr <- cbind(rna.dr1.tbl, rna.dr2.tbl)
  rm(rna.dr1.tbl)
  rm(rna.dr2.tbl)
  gc()
  ohsu.expr
}

ohsu.expr <- load.ohsu.expr.data()

rna.lab.ids <- gsub("\\.", "-", gsub("^X", "", colnames(ohsu.expr)))

## Make sure that the lab ids we get from reading the expression data are the same we get from the fastqs
## Get the RNA-seq fastq files
parent.id <- "syn5522959"
## For the rna-seq files, lab_ids are embedded with 00-00002 being the lab_id
## RNA150113BD_00-00002_AACGTGAT_L001_R1_001.fastq.gz 
rna.query <- synQuery(paste0("SELECT name FROM file WHERE parentId==\"", parent.id, "\""), blockSize = 500)
rna.fastq.df <- rna.query$fetch()
rna.fastq.df <- rna.query$collectAll()

rna.fastq.df$lab_id <- unlist(lapply(rna.fastq.df$file.name, 
                                     function(file.name) {
                                       re <- regexpr(text=file.name, pattern="\\d+-\\d+")
                                       substr(file.name, re[1], re[1] + attr(re, "match.length")-1)
                                     }))
rna.fastq.lab.ids <- unique(rna.fastq.df$lab_id)
if(all(rna.lab.ids %in% rna.fastq.lab.ids)) {
  cat("As expected, all lab ids in RNA-seq expression data correspond to (one or more) fastq files\n")
} else {
  warning("UNEXPECTED: some lab ids in RNA-seq expression data have no corresponding fastq files\n")
}

# Pull in the drug response data

## this MD5 is broken
## Read in the latest rna-seq dashboard
## BeatAML_rnaseq_2017_06_09_public_dashboard.xlsx (syn10008793)
## As stated in the README (datawave-2-q3_rnaseq_README.pdf),
## this is a cumulative dashboard that covers all current releases (i.e., DR1 and DR2)
## rna.obj <- synGet("syn10008793", downloadFile=TRUE, downloadLocation = download.path)
## file <- getFileLocation(rna.obj)
file <- paste0(beat.aml.cache, "BeatAML_rnaseq_2017_06_09_public_dashboard.xlsx")

## The samples summary is the 4th sheet
## This sample.summary table provides a map from patient id to sequence id.
## Expect that since the sample sequence ids where column headers in our expression
## data, they were corrupted from AA-BBBBB (where A and B are digits) to XAA.BBBBB.
## So change that.
cat(paste0("Reading sample summary file ", file, "\n"))
rnaseq.sample.summary.tbl <- read.xlsx(file, sheet=4)
rnaseq.sample.summary.tbl$Original_LabID <- as.character(rnaseq.sample.summary.tbl$Original_LabID)
rnaseq.sample.summary.tbl$SeqID <- as.character(rnaseq.sample.summary.tbl$SeqID)

rna.run.details.tbl <- read.xlsx(file, sheet=2)
## Merge in the control samples on sheet 3
rna.run.details.tbl <- rbind(rna.run.details.tbl, read.xlsx(file, sheet=3))

## Read in the raw OHSU inhibitor data inhibitor_data_points_2017_05_17.txt
## This latest file seems to have data from the prior releases
inh.obj <- synGet("syn10008778", downloadFile=TRUE, downloadLocation = download.path)
inh.tbl <- fread(getFileLocation(inh.obj))
inh.tbl <- as.data.frame(inh.tbl)

## Old inhibitor data:
inh.obj <- synGet("syn8149180", downloadFile=TRUE, downloadLocation = download.path)
inh.tbl.old <- fread(getFileLocation(inh.obj))
inh.tbl.old <- as.data.frame(inh.tbl.old)

m <- merge(inh.tbl, inh.tbl.old, by = c("patient_id", "lab_id", "inhibitor", "replicant", "well_concentration", "time_of_read"), all=TRUE)

## Drop the data--we are only interested here in the specimens (i.e., lab_id) and patients (patient_id)
## that have drug response data.
inh.tbl <- unique(inh.tbl[,c("patient_id", "inhibitor", "lab_id")])
pts.with.drug.response <- unique(inh.tbl[,c("patient_id"), drop=F])
specimens.with.drug.response <- unique(inh.tbl[,c("patient_id", "lab_id"), drop=F])
## As expected, no specimens are duplicated across patients
if(any(duplicated(specimens.with.drug.response$lab_id))) {
  warning("UNEXPECTED: some lab ids were duplicated across patients\n")
} else {
  cat("As expected, no lab ids were duplicated across patients\n")
}

## What is the difference between SeqID and Original_LabID
## Original_LabID looks to be the same as the lab_id (specimen id) in the drug response data.
## But, confusingly, as labId in the rnaseq data.
## That Original_LabID corresponds to lab_id in the drug response data is indicated in the README (syn8241614: datawave-2-q2_rnaseq_README.pdf)
## "Note that an additional unique barcode, termed a ‘SeqID’, has been assigned to all Wave 2 sequencing samples
## at the FASTQ stage. Mappings between the clinical LabID (termed ‘Original_LabID’) and the SeqID are
## provided in the dashboards."

## Of the rnaseq cases in which the Original_LabID is not the same as the SeqID ...
flag <- rnaseq.sample.summary.tbl$Original_LabID != rnaseq.sample.summary.tbl$SeqID
tmp <- rnaseq.sample.summary.tbl[flag, ]
## ... none of the SeqID's match the lab_id of the drug response data
if(any(tmp$SeqID %in% specimens.with.drug.response$lab_id)) {
  warning("UNEXPECTED: in cases in which SeqID != Original_LabID, some of the SeqIDs match lab_id's in the original drug response\n")
} else {
  cat("As expected, in cases in which SeqID != Original_LabID, none of the SeqIDs match lab_id's in the original drug response\n")
}
## [1] FALSE
## ... whereas 22 of the 78 Original_LabIDs match lab_ids
cat("Overlap of Original_LabID and drug response lab_id: unlike above, some match\n")
table(tmp$Original_LabID %in% specimens.with.drug.response$lab_id)
## FALSE  TRUE 
## 22    56 
## This shows that seqID, not original_LabID, is the same as labId in the rna.fastq.df
cat("Overlap of SeqID and lab_id embedded in RNA-seq filename\n")
table(tmp$SeqID %in% rna.fastq.df$lab_id)
## TRUE 
## 78 
if(any(tmp$Original_LabID %in% specimens.with.drug.response$lab_id)) {
  cat("As expected, in cases in which SeqID != Original_LabID, some of the Original_LabIDs match lab_id's in the original drug response\n")
} else {
  warning("UNEXPECTED: in cases in which SeqID != Original_LabID, none of the Original_LabIDs match lab_id's in the original drug response\n")
}

## NB: merging with SeqID, not Original_LabID, here:
rnaseq.fastq.and.summary.tbl <- merge(rnaseq.sample.summary.tbl, rna.fastq.df, by.y = "lab_id", by.x = "SeqID")
tmp <- merge(rnaseq.sample.summary.tbl, rna.fastq.df, by.y = "lab_id", by.x = "Original_LabID")
if(nrow(rnaseq.fastq.and.summary.tbl) > nrow(tmp)) {
  cat("As expected, get more mappings when merge with SeqID\n")
} else {
  cat("UNEXPECTED: got more mappings merging with Original_LabID\n")
}
rnaseq.fastq.and.summary.tbl <- merge(rnaseq.sample.summary.tbl[,c("SeqID", "PatientID", "Diagnosis", "SpecimenHasSeqCap", "PatientHasSeqCap")], rna.fastq.df, by.y = "lab_id", by.x = "SeqID", all = TRUE)
if(any(duplicated(rna.run.details.tbl$SeqID))) {
  cat("UNEXPECTED: duplicated SeqIDs\n")
} else {
  cat("As expected, no duplicated SeqIDs\n")
}
rnaseq.fastq.and.summary.tbl <- merge(rnaseq.fastq.and.summary.tbl, rna.run.details.tbl, by.x = "SeqID", by.y = "SeqID", all = TRUE)

## Merge in diagnosis info (most recent and at specimen acquisition)
old.nrows <- nrow(rnaseq.fastq.and.summary.tbl)
rnaseq.fastq.and.summary.tbl <- merge(rnaseq.fastq.and.summary.tbl, patient.diagnosis.tbl, by.x = c("PatientID", "Original_LabID"), by.y = c("patient_id", "lab_id"), all.x=TRUE)
new.nrows <- nrow(rnaseq.fastq.and.summary.tbl)
if(old.nrows != new.nrows) {
  warning("UNEXPECTED: rna-seq summary num rows changed!\n")
} else {
  cat("As expected, rna-seq summary num rows did not change\n")
}

cat("\nOverlap of Diagnosis and most_recent_diagnosis\n")
table(rnaseq.fastq.and.summary.tbl[,c("Diagnosis", "most_recent_diagnosis")])

cat("\nOverlap of Diagnosis and diagnosis_at_time_of_specimen_acquisition\n")
table(rnaseq.fastq.and.summary.tbl[,c("Diagnosis", "diagnosis_at_time_of_specimen_acquisition")])


rnaseq.fastq.and.summary.tbl$PatientID <- as.character(rnaseq.fastq.and.summary.tbl$PatientID)
rnaseq.fastq.and.summary.tbl$SeqID <- as.character(rnaseq.fastq.and.summary.tbl$SeqID)
rnaseq.fastq.and.summary.tbl$Original_LabID <- as.character(rnaseq.fastq.and.summary.tbl$Original_LabID)

rnaseq.fastq.and.summary.tbl$SpecimenHasDrugResponse <- rep(FALSE, nrow(rnaseq.fastq.and.summary.tbl))
rnaseq.fastq.and.summary.tbl$PatientHasDrugResponse <- rep(FALSE, nrow(rnaseq.fastq.and.summary.tbl))
flag <- !is.na(rnaseq.fastq.and.summary.tbl$Original_LabID) & (rnaseq.fastq.and.summary.tbl$Original_LabID %in% specimens.with.drug.response$lab_id)
rnaseq.fastq.and.summary.tbl$SpecimenHasDrugResponse[flag] <- TRUE
flag <- !is.na(rnaseq.fastq.and.summary.tbl$PatientID) & (rnaseq.fastq.and.summary.tbl$PatientID %in% specimens.with.drug.response$patient_id)
rnaseq.fastq.and.summary.tbl$PatientHasDrugResponse[flag] <- TRUE

flag <- rnaseq.fastq.and.summary.tbl$SpecimenHasDrugResponse == TRUE
cat(paste0("Num specimens in RNA-seq: ", nrow(rnaseq.fastq.and.summary.tbl), "\n"))
cat(paste0("Num specimens in RNA-seq with drug response: ", nrow(rnaseq.fastq.and.summary.tbl[flag,]), "\n"))

flag <- rnaseq.fastq.and.summary.tbl$PatientHasDrugResponse == TRUE
cat(paste0("Num patients in RNA-seq: ", length(unique(rnaseq.fastq.and.summary.tbl$PatientID)), "\n"))
## 400
cat(paste0("Num patients in RNA-seq with drug response: ", length(unique(rnaseq.fastq.and.summary.tbl$PatientID[rnaseq.fastq.and.summary.tbl$PatientHasDrugResponse == TRUE])), "\n"))
## 319

## Save the table having RNA-seq fastq files and whether the specimen/patient has drug-response/dna-seq
write.table(file="ohsu-rnaseq-fastq-summary.tsv", rnaseq.fastq.and.summary.tbl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

## BEGIN DNA

## Get lab ids from both the VCFs and from the fastq files and ensure that those derived from
## the fastq files are a strict superset of those from the vcf files--i.e., we have all of the
## raw data we expect to have.

## Read in the VCFs from Files/seqcap/vcf/ and parse the ids from the file names
## Like with the rna-seq data above, the dna-seq labIds map to Original_LabIDs,
## which are the ids used in the drug response data.  This is explicitly stated in
## the README (datawave-2_seqcap_README.pdf: syn8241603):
## "Note that an additional unique barcode, termed a ‘SeqID’, has been assigned to all Wave 2 sequencing samples
## at the FASTQ stage. Mappings between the clinical LabID (termed ‘Original_LabID’) and the SeqID are
## provided in the dashboards."
vcfdf    <- synQuery("SELECT name FROM file WHERE parentId==\"syn5522833\"")
## The first id in the name (following Pair) seems to be the tumor
vcfdf$lab_id = gsub("^Pair_","",gsub("_AML_.*$","",gsub("\\..*$","",vcfdf$file.name)))

dna.lab.ids <- unique(vcfdf$lab_id)

## Make sure that the lab ids we get from reading the VCF data are the same we get from the fastqs
## Get the DNA-seq fastq files
parent.id <- "syn5522891"
## For the dna-seq files, lab_ids are embedded with 00-00002 being the lab_id
## 141209_SN632_0415_BC5RKPACXX_DNA140626JT_27_14-00331_TAAGGCGA_L001_R1_001.fastq.gz
dna.query <- synQuery(paste0("SELECT name FROM file WHERE parentId==\"", parent.id, "\""), blockSize = 500)
dna.fastq.df <- dna.query$fetch()
dna.fastq.df <- dna.query$collectAll()

dna.fastq.df$lab_id <- unlist(lapply(dna.fastq.df$file.name, 
                                     function(file.name) {
                                       re <- regexpr(text=file.name, pattern="\\d+-\\d+")
                                       substr(file.name, re[1], re[1] + attr(re, "match.length")-1)
                                     }))
dna.fastq.lab.ids <- unique(dna.fastq.df$lab_id)
if(all(dna.lab.ids %in% dna.fastq.lab.ids)) {
  cat("As expected, all lab ids in DNA-seq expression data correspond to (one or more) fastq files\n")
} else {
  warning("UNEXPECTED: some lab ids in DNA-seq expression data have no corresponding fastq files\n")
}

## Read in the OHSU seqcap dashboard file (so we can subset to those patients that have AML)
## For DNA sheet 4 is all samples (tumor and normal)
## Sheet 5 is tumor only.
## Tumor-only means that there was no normal comparator.
## I believe that sheet 5 is a strict subset of sheet 4, but to be be sure
## merge the two.
## Read in the latest dna-seq dashboard
## BeatAML_seqcap_2017_06_09_public_dashboard.xlsx (syn10008791)
## As stated in the README (datawave-2_seqcap_README.pdf), this dashboard is
## cumulative for all of the current releases (i.e., DR1 and DR2)
## working
## This MD5 is broken
## dna.obj <- synGet("syn10008791", downloadFile=TRUE, downloadLocation = download.path)
## file <- getFileLocation(dna.obj)
file <- paste0(beat.aml.cache, "BeatAML_seqcap_2017_06_09_public_dashboard.xlsx")


## For DNA sheet 4 is all samples (tumor and normal)
## Sheet 5 is tumor only.  
## Tumor-only means that there was no normal comparator.
## I believe that sheet 5 is a strict subset of sheet 4, but to be be sure
## merge the two.
cat(paste0("Reading sample summary file ", file, "\n"))
dnaseq.sample.summary.tbl <- read.xlsx(file, sheet=4)
dnaseq.sample.summary.tbl <- unique(rbind(dnaseq.sample.summary.tbl, read.xlsx(file, sheet=5)))
dnaseq.sample.summary.tbl$Original_LabID <- as.character(dnaseq.sample.summary.tbl$Original_LabID)
dnaseq.sample.summary.tbl$SeqID <- as.character(dnaseq.sample.summary.tbl$SeqID)

## Unlike with rna-seq, there are no control samples for dna-sesq
dna.run.details.tbl <- read.xlsx(file, sheet=3)

## Now that some fields in the DNA run details columns have a semicolon (;), which leads
## to multiple lines output by read.xlsx2 with the same SeqID.  Since we don't care about
## these fields (FlowCell, Lane, CoreRunID, CoreFlowCellID)
dna.run.details.tbl <- unique(dna.run.details.tbl[, !(colnames(dna.run.details.tbl) %in% c("FlowCell", "Lane", "CoreRunID", "CoreFlowCellID"))])

## As we did above with rna-seq, empirically confirm that the Original_LabID of the dna-seq
## files, rather than the SeqID, corresponds to labID in the drug response data.

## Of the dnaseq cases in which the Original_LabID is not the same as the SeqID ...
flag <- dnaseq.sample.summary.tbl$Original_LabID != dnaseq.sample.summary.tbl$SeqID
tmp <- dnaseq.sample.summary.tbl[flag, ]
## ... none of the SeqID's match the lab_id of the drug response data
if(any(tmp$SeqID %in% specimens.with.drug.response$lab_id)) {
  warning("UNEXPECTED: in dna cases in which SeqID != Original_LabID, some of the SeqIDs match lab_id's in the original drug response\n")
} else {
  cat("As expected, in dna cases in which SeqID != Original_LabID, none of the SeqIDs match lab_id's in the original drug response\n")
}
## ... whereas 133 of the 236 Original_LabIDs match lab_ids
cat("Overlap of Original_LabID and drug response lab_id: unlike above, some match\n")
table(tmp$Original_LabID %in% specimens.with.drug.response$lab_id)
## FALSE  TRUE 
## 103    133 
## This shows that seqID, not original_LabID, is the same as labId in the dna.fastq.df:
cat("Overlap of SeqID and lab_id embedded in DNA-seq filename\n")
table(tmp$SeqID %in% dna.fastq.df$lab_id)
## TRUE 
## 236 
if(any(tmp$Original_LabID %in% specimens.with.drug.response$lab_id)) {
  cat("As expected, in DNA cases in which SeqID != Original_LabID, some of the Original_LabIDs match lab_id's in the original drug response\n")
} else {
  warning("UNEXPECTED: in DNA cases in which SeqID != Original_LabID, none of the Original_LabIDs match lab_id's in the original drug response\n")
}

## NB: merging with SeqID, not Original_LabID, here:
dnaseq.fastq.and.summary.tbl <- merge(dnaseq.sample.summary.tbl, dna.fastq.df, by.y = "lab_id", by.x = "SeqID")
tmp <- merge(dnaseq.sample.summary.tbl, dna.fastq.df, by.y = "lab_id", by.x = "Original_LabID")
if(nrow(dnaseq.fastq.and.summary.tbl) > nrow(tmp)) {
  cat("As expected, get more mappings when merge with SeqID\n")
} else {
  warning("UNEXPECTED: got more mappings merging with Original_LabID\n")
}
dnaseq.fastq.and.summary.tbl <- merge(dnaseq.sample.summary.tbl[,c("SeqID", "PatientID", "Diagnosis", "Specimen_HasRNASeq", "Patient_HasRNASeq")], dna.fastq.df, by.y = "lab_id", by.x = "SeqID", all = TRUE)
if(any(duplicated(dna.run.details.tbl$SeqID))) {
  warning("UNEXPECTED: duplicated SeqIDs\n")
} else {
  cat("As expected, no duplicated SeqIDs\n")
}
dnaseq.fastq.and.summary.tbl <- merge(dnaseq.fastq.and.summary.tbl, dna.run.details.tbl, by.x = "SeqID", by.y = "SeqID", all = TRUE)

## Merge in diagnosis info (most recent and at specimen acquisition)
old.nrows <- nrow(dnaseq.fastq.and.summary.tbl)
dnaseq.fastq.and.summary.tbl <- merge(dnaseq.fastq.and.summary.tbl, patient.diagnosis.tbl, by.x = c("PatientID", "Original_LabID"), by.y = c("patient_id", "lab_id"), all.x=TRUE)
new.nrows <- nrow(dnaseq.fastq.and.summary.tbl)
if(old.nrows != new.nrows) {
  warning("UNEXPECTED: dna-seq summary num rows changed!\n")
} else {
  cat("As expected, dna-seq summary num rows did not change\n")
}

cat("\nOverlap of Diagnosis and most_recent_diagnosis\n")
table(dnaseq.fastq.and.summary.tbl[,c("Diagnosis", "most_recent_diagnosis")])

cat("\nOverlap of Diagnosis and diagnosis_at_time_of_specimen_acquisition\n")
table(dnaseq.fastq.and.summary.tbl[,c("Diagnosis", "diagnosis_at_time_of_specimen_acquisition")])

dnaseq.fastq.and.summary.tbl$PatientID <- as.character(dnaseq.fastq.and.summary.tbl$PatientID)
dnaseq.fastq.and.summary.tbl$SeqID <- as.character(dnaseq.fastq.and.summary.tbl$SeqID)
dnaseq.fastq.and.summary.tbl$Original_LabID <- as.character(dnaseq.fastq.and.summary.tbl$Original_LabID)

dnaseq.fastq.and.summary.tbl$SpecimenHasDrugResponse <- rep(FALSE, nrow(dnaseq.fastq.and.summary.tbl))
dnaseq.fastq.and.summary.tbl$PatientHasDrugResponse <- rep(FALSE, nrow(dnaseq.fastq.and.summary.tbl))
flag <- !is.na(dnaseq.fastq.and.summary.tbl$Original_LabID) & (dnaseq.fastq.and.summary.tbl$Original_LabID %in% specimens.with.drug.response$lab_id)
dnaseq.fastq.and.summary.tbl$SpecimenHasDrugResponse[flag] <- TRUE
flag <- !is.na(dnaseq.fastq.and.summary.tbl$PatientID) & (dnaseq.fastq.and.summary.tbl$PatientID %in% specimens.with.drug.response$patient_id)
dnaseq.fastq.and.summary.tbl$PatientHasDrugResponse[flag] <- TRUE

flag <- dnaseq.fastq.and.summary.tbl$SpecimenHasDrugResponse == TRUE
cat(paste0("Num specimens in DNA-seq: ", nrow(dnaseq.fastq.and.summary.tbl), "\n"))
cat(paste0("Num specimens in DNA-seq with drug response: ", nrow(dnaseq.fastq.and.summary.tbl[flag,]), "\n"))

flag <- dnaseq.fastq.and.summary.tbl$PatientHasDrugResponse == TRUE
cat(paste0("Num patients in DNA-seq: ", length(unique(dnaseq.fastq.and.summary.tbl$PatientID)), "\n"))
## 452
cat(paste0("Num patients in DNA-seq with drug response: ", length(unique(dnaseq.fastq.and.summary.tbl$PatientID[dnaseq.fastq.and.summary.tbl$PatientHasDrugResponse == TRUE])), "\n"))
## 308
## Save the table having DNA-seq fastq files and whether the specimen/patient has drug-response/rna-seq
write.table(file="ohsu-dnaseq-fastq-summary.tsv", dnaseq.fastq.and.summary.tbl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

## END DNA

## TODO:

library(gplots)
inh.patient.ids <- unique(inh.tbl$patient_id)
rna.patient.ids <- unique(rnaseq.fastq.and.summary.tbl$PatientID)
dna.patient.ids <- unique(dnaseq.fastq.and.summary.tbl$PatientID)
aml.patient.ids <- unique(aml.diagnosis.tbl$patient_id)
seq.aml.patient.ids <- union(rnaseq.fastq.and.summary.tbl$PatientID[!is.na(rnaseq.fastq.and.summary.tbl) & (rnaseq.fastq.and.summary.tbl$Diagnosis %in% aml.diagnosis)], dnaseq.fastq.and.summary.tbl$PatientID[!is.na(dnaseq.fastq.and.summary.tbl) & (dnaseq.fastq.and.summary.tbl$Diagnosis %in% aml.diagnosis)])
aml.patient.ids <- unique(intersect(aml.patient.ids, seq.aml.patient.ids))
vennList = list("DRC"=inh.patient.ids,"RNA-seq"= rna.patient.ids,"WES"=dna.patient.ids, "AML"=aml.patient.ids)

png("venn-drc-rna-dna.png")
venn(vennList)
title(main="BEAT AML Assay Overlap")
d <- dev.off()

save.image(".RData")

