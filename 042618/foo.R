

rna.seq.ids <- gsub("\\.", "-", gsub("^X", "", colnames(ohsu.expr)))

## Make sure that the SeqIDs we get from reading the expression data are the same we get from the fastqs
## Get the RNA-seq fastq files
parent.id <- "syn5522959"
## For the rna-seq files, seqIDs are embedded with 00-00002 being the seqID
## RNA150113BD_00-00002_AACGTGAT_L001_R1_001.fastq.gz 
rna.query <- synQuery(paste0("SELECT name FROM file WHERE parentId==\"", parent.id, "\""), blockSize = 500)
rna.fastq.df <- rna.query$fetch()
rna.fastq.df <- rna.query$collectAll()

rna.fastq.df$SeqID <- unlist(lapply(rna.fastq.df$file.name, 
                                     function(file.name) {
                                       re <- regexpr(text=file.name, pattern="\\d+-\\d+")
                                       substr(file.name, re[1], re[1] + attr(re, "match.length")-1)
                                     }))
rna.fastq.seq.ids <- unique(rna.fastq.df$SeqID)
if(all(rna.seq.ids %in% rna.fastq.seq.ids)) {
  cat("\nAs expected, SeqIDs from RNA-seq expression data correspond to (one or more) fastq files\n")
} else {
  warning("\nUNEXPECTED: some SeqIDs in RNA-seq expression data have no corresponding fastq files\n")
}

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

rna.run.details.tbl <- read.xlsx(file, sheet=2)
## Merge in the control samples on sheet 3
rna.run.details.tbl <- rbind(rna.run.details.tbl, read.xlsx(file, sheet=3))

## Read in the raw OHSU inhibitor data inhibitor_data_points_2018_01_12.txt
## This latest file seems to have data from the prior releases
inh.id <- "syn11767000"
inh.obj <- synGet(inh.id, downloadFile=TRUE)
inh.tbl <- fread(getFileLocation(inh.obj))
inh.tbl <- as.data.frame(inh.tbl)

inh.and.rnaseq.samples <- unique(inh.tbl[, c("patient_id", "lab_id")])
inh.and.rnaseq.samples <- merge(inh.and.rnaseq.samples, unique(rnaseq.sample.summary.tbl[, c("Original_LabID", "SeqID")]),
                                by.x = "lab_id", by.y = "Original_LabID", all = FALSE)

rna.dr2q2.seq.ids.only <- rna.dr2q2.seq.ids[!(rna.dr2q2.seq.ids %in% rna.dr2q1.seq.ids)]
rna.dr2q2.seq.ids.only <- rna.dr2q2.seq.ids.only[!(rna.dr2q2.seq.ids.only %in% rna.dr1.seq.ids)]
inh.and.rnaseq.dr2q2.samples <- subset(inh.and.rnaseq.samples, SeqID %in% rna.dr2q2.seq.ids.only)

## ## Read in the raw OHSU inhibitor data inhibitor_data_points_2017_05_17.txt
## ## This latest file seems to have data from the prior releases
## inh.obj <- synGet("syn10008778", downloadFile=TRUE, downloadLocation = download.path)
## inh.tbl <- fread(getFileLocation(inh.obj))
## inh.tbl <- as.data.frame(inh.tbl)

## Old inhibitor data:
## inh.obj <- synGet("syn8149180", downloadFile=TRUE, downloadLocation = download.path)
## inh.tbl.old <- fread(getFileLocation(inh.obj))
## inh.tbl.old <- as.data.frame(inh.tbl.old)
## m <- merge(inh.tbl, inh.tbl.old, by = c("patient_id", "lab_id", "inhibitor", "replicant", "well_concentration", "time_of_read"), all=TRUE)

## Drop the data--we are only interested here in the specimens (i.e., lab_id) and patients (patient_id)
## that have drug response data.
inh.tbl <- unique(inh.tbl[,c("patient_id", "inhibitor", "lab_id")])
pts.with.drug.response <- unique(inh.tbl[,c("patient_id"), drop=F])
specimens.with.drug.response <- unique(inh.tbl[,c("patient_id", "lab_id"), drop=F])
## As expected, no specimens are duplicated across patients
if(any(duplicated(specimens.with.drug.response$lab_id))) {
  warning("\nUNEXPECTED: some lab ids were duplicated across patients\n")
} else {
  cat("\nAs expected, no lab ids were duplicated across patients\n")
}

## Original_LabID in the dashboard files is the lab_id (specimen id) in the drug response data.
## This is stated in syn8241614: datawave-2-q2_rnaseq_README.pdf:
## "Note that an additional unique barcode, termed a ‘SeqID’, has been assigned to all Wave 2 sequencing samples
## at the FASTQ stage. Mappings between the clinical LabID (termed ‘Original_LabID’) and the SeqID are
## provided in the dashboards."

## The following simply confirms that SeqID, rather than Original_LabID, is the ID used in the RNA-seq data.
## Of the rnaseq cases in which the Original_LabID is not the same as the SeqID ...
flag <- rnaseq.sample.summary.tbl$Original_LabID != rnaseq.sample.summary.tbl$SeqID
tmp <- rnaseq.sample.summary.tbl[flag, ]
## ... none of the SeqID's match the lab_id of the drug response data
if(any(tmp$SeqID %in% specimens.with.drug.response$lab_id)) {
  warning("\nUNEXPECTED: in cases in which SeqID != Original_LabID, some of the SeqIDs match lab_id's in the original drug response\n")
} else {
  cat("\nAs expected, in cases in which SeqID != Original_LabID, none of the SeqIDs match lab_id's in the original drug response\n")
}
## [1] FALSE
## ... whereas 22 of the 78 Original_LabIDs match lab_ids
cat("\nOverlap of Original_LabID and drug response lab_id: unlike above, some match\n")
table(tmp$Original_LabID %in% specimens.with.drug.response$lab_id)
## FALSE  TRUE 
## 22    56 

## The following shows that seqID, not original_LabID, is the same as SeqID in the rna.fastq.df
cat("\nOverlap of SeqID from dashboard and id embedded in RNA-seq filename\n")
table(tmp$SeqID %in% rna.fastq.df$SeqID)
## TRUE 
## 78 
if(any(tmp$Original_LabID %in% specimens.with.drug.response$lab_id)) {
  cat("\nAs expected, in cases in which SeqID != Original_LabID, some of the Original_LabIDs match lab_id's in the original drug response\n")
} else {
  warning("\nUNEXPECTED: in cases in which SeqID != Original_LabID, none of the Original_LabIDs match lab_id's in the original drug response\n")
}

## Compare merging id extraced from RNA-seq filename (i.e., rna.fastq.df$SeqID) with SeqID here:
rnaseq.fastq.and.summary.tbl <- merge(rnaseq.sample.summary.tbl, rna.fastq.df, by.y = "SeqID", by.x = "SeqID")
## With merging id extract from RNA-seq filename with Original_LabID here:
tmp <- merge(rnaseq.sample.summary.tbl, rna.fastq.df, by.y = "SeqID", by.x = "Original_LabID")
if(nrow(rnaseq.fastq.and.summary.tbl) > nrow(tmp)) {
  cat("\nAs expected, get more mappings when merge with SeqID\n")
} else {
  warning("\nUNEXPECTED: got more mappings merging with Original_LabID\n")
}

## Merge ids extracted from RNA-seq filename with SeqID, but with all = TRUE and keeping 
## only the subset of the columns we are interested in.
rnaseq.fastq.and.summary.tbl <- merge(rnaseq.sample.summary.tbl[,c("SeqID", "PatientID", "Diagnosis", "SpecimenHasSeqCap", "PatientHasSeqCap")], rna.fastq.df, by.y = "SeqID", by.x = "SeqID", all = TRUE)
if(any(duplicated(rna.run.details.tbl$SeqID))) {
  cat("\nUNEXPECTED: duplicated RNA-seq SeqIDs\n")
} else {
  cat("\nAs expected, no duplicated RNA-seq SeqIDs\n")
}
rnaseq.fastq.and.summary.tbl <- merge(rnaseq.fastq.and.summary.tbl, rna.run.details.tbl, by.x = "SeqID", by.y = "SeqID", all = TRUE)

## Merge in diagnosis info (most recent and at specimen acquisition)--here using Original_LabID, which was merged
## from dashboard, and lab_id from patient diagnosis file.
old.nrows <- nrow(rnaseq.fastq.and.summary.tbl)
rnaseq.fastq.and.summary.tbl <- merge(rnaseq.fastq.and.summary.tbl, patient.diagnosis.tbl, by.x = c("PatientID", "Original_LabID"), by.y = c("patient_id", "lab_id"), all.x=TRUE)
new.nrows <- nrow(rnaseq.fastq.and.summary.tbl)
if(old.nrows != new.nrows) {
  warning("\nUNEXPECTED: rna-seq summary num rows changed!\n")
} else {
  cat("\nAs expected, rna-seq summary num rows did not change\n")
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

## Get SeqIDs from both the VCFs and from the fastq files and ensure that those derived from
## the fastq files are a strict superset of those from the vcf files--i.e., we have all of the
## raw data we expect to have.

## Read in the VCFs from Files/seqcap/vcf/ and parse the ids from the file names
## Like with the rna-seq data above, the mapping between SeqIDs (used for the dna-seq data) and the
## Original_LabIDs is provided in the dashboard.  This is explicitly stated in
## the README (datawave-2_seqcap_README.pdf: syn8241603):
## "Note that an additional unique barcode, termed a ‘SeqID’, has been assigned to all Wave 2 sequencing samples
## at the FASTQ stage. Mappings between the clinical LabID (termed ‘Original_LabID’) and the SeqID are
## provided in the dashboards."
## Confusingly, datawave-2_seqcap_README.pdf lists Sample_DNA*_LabID (as opposed to SeqID as used in the rna-seq README)
## under the Folder Structure.  This is probably a typo--let's confirm though, as we did above, that the id 
## embedded in the filename is actually a SeqID.
vcfdf    <- synQuery("SELECT name FROM file WHERE parentId==\"syn5522833\"")
## The first id in the name (following Pair) seems to be the tumor
vcfdf$SeqID = gsub("^Pair_","",gsub("_AML_.*$","",gsub("\\..*$","",vcfdf$file.name)))

dna.seq.ids <- unique(vcfdf$SeqID)

## Make sure that the ids we are extracting from the VCF filenames are the same we get from the fastqs
## Get the DNA-seq fastq files
parent.id <- "syn5522891"
## For the dna-seq files, SeqIDs are embedded within the filename with 14-00331 being the SeqID:
## 141209_SN632_0415_BC5RKPACXX_DNA140626JT_27_14-00331_TAAGGCGA_L001_R1_001.fastq.gz
dna.query <- synQuery(paste0("SELECT name FROM file WHERE parentId==\"", parent.id, "\""), blockSize = 500)
dna.fastq.df <- dna.query$fetch()
dna.fastq.df <- dna.query$collectAll()

dna.fastq.df$SeqID <- unlist(lapply(dna.fastq.df$file.name, 
                                     function(file.name) {
                                       re <- regexpr(text=file.name, pattern="\\d+-\\d+")
                                       substr(file.name, re[1], re[1] + attr(re, "match.length")-1)
                                     }))
dna.fastq.seq.ids <- unique(dna.fastq.df$SeqID)
if(all(dna.seq.ids %in% dna.fastq.seq.ids)) {
  cat("\nAs expected, all ids extracted from DNA-seq VCF files correspond to (one or more) fastq files\n")
} else {
  warning("\nUNEXPECTED: some ids extracted from DNA-seq VCF files have no corresponding fastq files\n")
}

## Read in the OHSU seqcap dashboard file (so we can subset to those patients that have AML)
## For DNA sheet 4 is all samples (tumor and normal)
## Sheet 5 is tumor only.
## Tumor-only means that there was no normal comparator.
## I believe that sheet 5 is a strict subset of sheet 4, but to be be sure
## merge the two.
## Read in the latest dna-seq dashboard
## BeatAML_seqcap_2018_01_29_public_dashboard.xlsx (syn11767025)
## As stated in the README (datawave-2_seqcap_README.pdf), this dashboard is
## cumulative for all of the current releases (i.e., DR1 and DR2)
dna.obj <- synGet("syn11767025", downloadFile=TRUE, downloadLocation = download.path)
file <- getFileLocation(dna.obj)

## For DNA sheet 4 is all samples (tumor and normal)
## Sheet 5 is tumor only.  
## Tumor-only means that there was no normal comparator.
## I believe that sheet 6 is a strict subset of sheet 4, but to be be sure
## merge the two.
cat(paste0("Reading sample summary file ", file, "\n"))
dnaseq.sample.summary.tbl <- read.xlsx(file, sheet=4)
dnaseq.sample.summary.tbl <- unique(rbind(dnaseq.sample.summary.tbl, read.xlsx(file, sheet=6)))
dnaseq.sample.summary.tbl$Original_LabID <- as.character(dnaseq.sample.summary.tbl$Original_LabID)
dnaseq.sample.summary.tbl$SeqID <- as.character(dnaseq.sample.summary.tbl$SeqID)

## Unlike with rna-seq, there are no control samples for dna-sesq
dna.run.details.tbl <- read.xlsx(file, sheet=3)

## Note that some fields in the DNA run details columns have a semicolon (;), which leads
## to multiple lines output by read.xlsx2 with the same SeqID.  Since we don't care about
## these fields (FlowCell, Lane, CoreRunID, CoreFlowCellID), just drop them.
dna.run.details.tbl <- unique(dna.run.details.tbl[, !(colnames(dna.run.details.tbl) %in% c("FlowCell", "Lane", "CoreRunID", "CoreFlowCellID"))])

## As we did above with rna-seq, empirically confirm that the Original_LabID of the dna-seq
## files, rather than the SeqID, corresponds to labID in the drug response data.

## Of the dnaseq cases in which the Original_LabID is not the same as the SeqID ...
flag <- dnaseq.sample.summary.tbl$Original_LabID != dnaseq.sample.summary.tbl$SeqID
tmp <- dnaseq.sample.summary.tbl[flag, ]
## ... none of the SeqID's match the lab_id of the drug response data
if(any(tmp$SeqID %in% specimens.with.drug.response$lab_id)) {
  warning("\nUNEXPECTED: in DNA cases in which SeqID != Original_LabID, some of the SeqIDs match lab_id's in the original drug response\n")
} else {
  cat("\nAs expected, in DNA cases in which SeqID != Original_LabID, none of the SeqIDs match lab_id's in the original drug response\n")
}
## ... whereas 133 of the 236 Original_LabIDs match lab_ids
cat("\nOverlap of Original_LabID and drug response lab_id: unlike above, some match\n")
table(tmp$Original_LabID %in% specimens.with.drug.response$lab_id)
## FALSE  TRUE 
## 103    133 

## This shows that seqID, not original_LabID, is the same as the id extracted from the fastq filenames (in dna.fastq.df$SeqID)
cat("\nOverlap of SeqID and lab_id embedded in DNA-seq filename\n")
table(tmp$SeqID %in% dna.fastq.df$SeqID)
## TRUE 
## 236 
if(any(tmp$Original_LabID %in% specimens.with.drug.response$lab_id)) {
  cat("\nAs expected, in DNA cases in which SeqID != Original_LabID, some of the Original_LabIDs match lab_id's in the original drug response\n")
} else {
  warning("\nUNEXPECTED: in DNA cases in which SeqID != Original_LabID, none of the Original_LabIDs match lab_id's in the original drug response\n")
}


## Compare merging id extraced from DNA-seq fastq filename (i.e., dna.fastq.df$SeqID) with SeqID from dashboard here:
dnaseq.fastq.and.summary.tbl <- merge(dnaseq.sample.summary.tbl, dna.fastq.df, by.y = "SeqID", by.x = "SeqID")
## With merging id extract from RNA-seq filename with Original_LabID here:
tmp <- merge(dnaseq.sample.summary.tbl, dna.fastq.df, by.y = "SeqID", by.x = "Original_LabID")
if(nrow(dnaseq.fastq.and.summary.tbl) > nrow(tmp)) {
  cat("\nAs expected, get more mappings when merge with SeqID\n")
} else {
  warning("\nUNEXPECTED: got more mappings merging with Original_LabID\n")
}
dnaseq.fastq.and.summary.tbl <- merge(dnaseq.sample.summary.tbl[,c("SeqID", "PatientID", "Diagnosis", "Specimen_HasRNASeq", "Patient_HasRNASeq")], dna.fastq.df, by.y = "SeqID", by.x = "SeqID", all = TRUE)

## Actually found some a single case/pair of duplicated SeqIDs:
## SeqID Original_LabID SampleGroup CaptureGroup Sequential_CaptureGroup CoreProjectID CoreSampleID Outliers
## 14-00165 14-00165 9 1 1 DNA140624JT DNA140624JT_10_14-00165 NA
## 14-00165 14-00165 9 2 2 DNA140624JT DNA140624JT_10_14-00165 NA
## But these differ only in fields we don't care about: CaptureGroup and Sequential_CaptureGroup
dna.run.details.tbl <- unique(dna.run.details.tbl[, !(colnames(dna.run.details.tbl) %in% c("CaptureGroup", "Sequential_CaptureGroup"))])

if(any(duplicated(dna.run.details.tbl$SeqID))) {
  warning("\nUNEXPECTED: duplicated DNA-seq SeqIDs\n")
  tmp <- dna.run.details.tbl
  tmp <- tmp[order(tmp$SeqID),]
  write.table(tmp[duplicated(tmp$SeqID, fromLast = TRUE) | duplicated(tmp$SeqID, fromLast = FALSE), ], quote=FALSE, col.names=TRUE, row.names=FALSE)
} else {
  cat("\nAs expected, no duplicated DNA-seq SeqIDs\n")
}
dnaseq.fastq.and.summary.tbl <- merge(dnaseq.fastq.and.summary.tbl, dna.run.details.tbl, by.x = "SeqID", by.y = "SeqID", all = TRUE)

## Merge in diagnosis info (most recent and at specimen acquisition)
old.nrows <- nrow(dnaseq.fastq.and.summary.tbl)
dnaseq.fastq.and.summary.tbl <- merge(dnaseq.fastq.and.summary.tbl, patient.diagnosis.tbl, by.x = c("PatientID", "Original_LabID"), by.y = c("patient_id", "lab_id"), all.x=TRUE)
new.nrows <- nrow(dnaseq.fastq.and.summary.tbl)
if(old.nrows != new.nrows) {
  warning("\nUNEXPECTED: dna-seq summary num rows changed!\n")
} else {
  cat("\nAs expected, dna-seq summary num rows did not change\n")
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

inh.patient.ids <- unique(inh.tbl$patient_id)
rna.patient.ids <- unique(rnaseq.fastq.and.summary.tbl$PatientID)
dna.patient.ids <- unique(dnaseq.fastq.and.summary.tbl$PatientID)
aml.patient.ids <- unique(aml.diagnosis.tbl$patient_id)
seq.aml.patient.ids <- union(rnaseq.fastq.and.summary.tbl$PatientID[!is.na(rnaseq.fastq.and.summary.tbl) & (rnaseq.fastq.and.summary.tbl$Diagnosis %in% aml.diagnosis)], dnaseq.fastq.and.summary.tbl$PatientID[!is.na(dnaseq.fastq.and.summary.tbl) & (dnaseq.fastq.and.summary.tbl$Diagnosis %in% aml.diagnosis)])
aml.patient.ids <- unique(intersect(aml.patient.ids, seq.aml.patient.ids))
vennList = list("screens"=inh.patient.ids,"RNA-seq"= rna.patient.ids,"WES"=dna.patient.ids, "AML"=aml.patient.ids)

tiff("venn-drc-rna-dna-dr3.tiff")
venn(vennList)
title(main="OHSU Assay Overlap")
d <- dev.off()

q(status = 0)

## Write some summary statistics
cat("\n")
cat(paste0("Number DNA-seq fastqs: ", nrow(dna.fastq.df), "\n"))

## Get DNA-seq VCFs (break out by mutect, varscan2 indel, and varscan2 SNP)
parent.id <- "syn5522833"
query <- synQuery(paste0("SELECT name FROM file WHERE parentId==\"", parent.id, "\""), blockSize = 500)
df <- query$fetch()
df <- query$collectAll()
cat(paste0("Number DNA-seq mutect VCFs: ", length(df$file.name[grepl(df$file.name, pattern="mutect", ignore.case=TRUE)]), "\n"))
cat(paste0("Number DNA-seq varscan2 SNP VCFs: ", length(df$file.name[grepl(df$file.name, pattern="varscan2.*snp", ignore.case=TRUE)]), "\n"))
cat(paste0("Number DNA-seq varscan2 indel VCFs: ", length(df$file.name[grepl(df$file.name, pattern="varscan2.*indel", ignore.case=TRUE)]), "\n"))

## Get DNA-seq MAFs (break out by mutect, varscan2 indel, and varscan2 SNP)
parent.id <- "syn5522834"
query <- synQuery(paste0("SELECT name FROM file WHERE parentId==\"", parent.id, "\""), blockSize = 500)
df <- query$fetch()
df <- query$collectAll()
cat(paste0("Number DNA-seq mutect MAFs: ", length(df$file.name[grepl(df$file.name, pattern="mutect", ignore.case=TRUE)]), "\n"))
cat(paste0("Number DNA-seq varscan2 SNP MAFs: ", length(df$file.name[grepl(df$file.name, pattern="varscan2.*snp", ignore.case=TRUE)]), "\n"))
cat(paste0("Number DNA-seq varscan2 indel MAFs: ", length(df$file.name[grepl(df$file.name, pattern="varscan2.*indel", ignore.case=TRUE)]), "\n"))

## Get DNA-seq BAMs
parent.id <- "syn5522831"
query <- synQuery(paste0("SELECT name FROM file WHERE parentId==\"", parent.id, "\""), blockSize = 500)
df <- query$fetch()
df <- query$collectAll()
cat(paste0("Number DNA-seq BAMs: ", length(df$file.name[grepl(df$file.name, pattern="bam", ignore.case=TRUE)]), "\n"))

cat(paste0("\nNumber RNA-seq fastqs: ", nrow(rna.fastq.df), "\n"))
## Get RNA-seq BAMs
parent.id <- "syn5522955"
query <- synQuery(paste0("SELECT name FROM file WHERE parentId==\"", parent.id, "\""), blockSize = 500)
df <- query$fetch()
df <- query$collectAll()
cat(paste0("Number RNA-seq BAMs: ", length(df$file.name[grepl(df$file.name, pattern="bam", ignore.case=TRUE)]), "\n"))

## Get RNA-seq fusion files
parent.id <- "syn7765693"
query <- synQuery(paste0("SELECT name FROM file WHERE parentId==\"", parent.id, "\""), blockSize = 500)
df <- query$fetch()
df <- query$collectAll()
cat(paste0("Number RNA-seq fusion files: ", length(df$file.name[grepl(df$file.name, pattern="junctions.bed", ignore.case=TRUE)]), "\n"))

## Output the number of drugs
synId <- inh.id
inh.obj <- synGet(synId, downloadFile=TRUE)
ohsu.inh.tbl <- fread(getFileLocation(inh.obj))
ohsu.inh.tbl <- as.data.frame(ohsu.inh.tbl)
ohsu.inh.tbl <- subset(ohsu.inh.tbl, diagnosis == aml.diagnosis)
ohsu.inh.tbl <- unique(ohsu.inh.tbl[, c("patient_id", "lab_id", "inhibitor", "diagnosis")])
inhibitors <- unique(ohsu.inh.tbl$inhibitor)
rm(ohsu.inh.tbl)
gc()

combo.therapies <- inhibitors[grepl(inhibitors, pattern=" - ")]
mono.therapies <- inhibitors[!(inhibitors %in% combo.therapies)]
cat("Number mono therapies: ", length(mono.therapies), "\n")
cat("Number combination therapies: ", length(combo.therapies), "\n")

save.image(".Rdata.venn")

