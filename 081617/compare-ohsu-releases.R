
library(synapseClient)
library(data.table)
library(openxlsx)

synapseLogin()

## Load the OHSU diagnostics -- the latest appears to have all patient info for all data releases
## Diagnosis_Labs_Treatments_Outcomes_2017_01_12.xlsx (syn8149174)
synId <- "syn8149174"
obj <- synGet(id=synId, downloadFile = TRUE)
file <- getFileLocation(obj)
ohsu.diagnosis.tbl.1.12.17 <- read.xlsx(file, sheet = 1)

## Get the latest OHSU diagnostics as well.
## Diagnosis_Labs_Treatments_Outcomes_2017_05_17.xlsx (syn10008777)
synId <- "syn10008777"
obj <- synGet(id=synId, downloadFile = TRUE)
file <- getFileLocation(obj)
ohsu.diagnosis.tbl.5.17.17 <- read.xlsx(file, sheet = 1)

## Read in the diagnosis that are in the RNA-seq dashboard
## BeatAML_rnaseq_2017_01_19_public_dashboard.xlsx (syn8149218)
## As stated in the README (datawave-2-q3_rnaseq_README.pdf),
## this is a cumulative dashboard that covers all current releases (i.e., DR1 and DR2)
synId <- "syn8149218"
rna.obj <- synGet(synId, downloadFile=TRUE)
file <- getFileLocation(rna.obj)
## The samples summary is the 4th sheet
## This sample.summary table provides a map from patient id to sequence id.
## Except that since the sample sequence ids where column headers in our expression
## data, they were corrupted from AA-BBBBB (where A and B are digits) to XAA.BBBBB.
## So change that.
cat(paste0("Reading sample summary file ", file, "\n"))
rnaseq.sample.summary.tbl.1.19.17 <- openxlsx:::read.xlsx(file, sheet=4)

## Read in the latest RNA-seq dashboard
## BeatAML_rnaseq_2017_06_09_public_dashboard.xlsx (syn10008793)
synId <- "syn10008793"
rna.obj <- synGet(synId, downloadFile=TRUE)
file <- getFileLocation(rna.obj)
rnaseq.sample.summary.tbl.6.9.17 <- openxlsx:::read.xlsx(file, sheet=4)

all.diagnoses.1.17 <- merge(ohsu.diagnosis.tbl.1.12.17[, c("patient_id", "most_recent_diagnosis", "diagnosis_at_time_of_specimen_acquisition")], rnaseq.sample.summary.tbl.1.19.17[, c("Original_LabID", "Diagnosis")], all = TRUE, by.x = "patient_id", by.y = "Original_LabID")
all.diagnoses.5.17 <- merge(ohsu.diagnosis.tbl.5.17.17[, c("patient_id", "most_recent_diagnosis", "diagnosis_at_time_of_specimen_acquisition")], rnaseq.sample.summary.tbl.6.9.17[, c("Original_LabID", "Diagnosis")], all = TRUE, by.x = "patient_id", by.y = "Original_LabID")

# Read in the inhibitor_data_points_2017_01_12.txt inhibitor data
synId <- "syn8149180"
inh.obj <- synGet(synId, downloadFile=TRUE)
ohsu.inh.1.12.tbl <- fread(getFileLocation(inh.obj))
ohsu.inh.1.12.tbl <- as.data.frame(ohsu.inh.1.12.tbl)

## Read in the inhibitor_data_points_2017_05_17.txt inhibitor data
synId <- "syn10008778"
inh.obj <- synGet(synId, downloadFile=TRUE)
ohsu.inh.5.17.tbl <- fread(getFileLocation(inh.obj))
ohsu.inh.5.17.tbl <- as.data.frame(ohsu.inh.5.17.tbl)

## Are there any drugs unique 2017_01_12 data
df <- data.frame(inhibitor = unique(ohsu.inh.1.12.tbl$inhibitor)[!(unique(ohsu.inh.1.12.tbl$inhibitor) %in% unique(ohsu.inh.5.17.tbl$inhibitor))])
write.table(file="inhibitors-missing-from-inh-5-17.tsv", df, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

missing.patients <- data.frame(patient_id = unique(ohsu.inh.1.12.tbl$patient_id)[!(unique(ohsu.inh.1.12.tbl$patient_id) %in% unique(ohsu.inh.5.17.tbl$patient_id))])
df <- missing.patients
missing.patients.1.17 <- merge(missing.patients, all.diagnoses.1.17, all.x=TRUE)
missing.patients.5.17 <- merge(missing.patients, all.diagnoses.5.17, all.x=TRUE)

write.table(file="patient_ids-missing-from-inh-5-17-diagnosis-1-17.xls", missing.patients.1.17, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
write.table(file="patient_ids-missing-from-inh-5-17-diagnosis-5-17.xls", missing.patients.5.17, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

df <- data.frame(lab_id = unique(ohsu.inh.1.12.tbl$lab_id)[!(unique(ohsu.inh.1.12.tbl$lab_id) %in% unique(ohsu.inh.5.17.tbl$lab_id))])
write.table(file="lab_ids-missing-from-inh-5-17.tsv", df, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
