library(synapseClient)

synapseLogin()

## Download the table of rna-seq files
synId <- "syn8739660"
obj <- synGet(synId, downloadFile=TRUE)
path <- getFileLocation(obj)
rnaseq.files <- read.table(path, sep="\t", header=TRUE)

aml.diagnosis <- "ACUTE MYELOID LEUKAEMIA (AML) AND RELATED PRECURSOR NEOPLASMS"
diagnosis.cols <- c("Diagnosis", "most_recent_diagnosis", "diagnosis_at_time_of_specimen_acquisition")

any.aml.diagnosis.rna.flag <- unlist(apply(rnaseq.files[, diagnosis.cols], 1, function(row) { aml.diagnosis %in% row }))
rnaseq.files.for.reprocessing <- rnaseq.files[any.aml.diagnosis.rna.flag,]

## We've been having problems downloading those from sample group 11.
## Restrict to those cases and try to download.

rnaseq.files.for.reprocessing <- subset(rnaseq.files.for.reprocessing, SampleGroup == 11)

download.path <- "download"
dir.create(download.path)

cat("Attempting to download ", length(rnaseq.files.for.reprocessing$file.id), " files\n")
for(synId in rnaseq.files.for.reprocessing$file.id) {
    cat(paste0("Downloading ", synId, "\n"))
    synGet(synId, downloadFile=TRUE, downloadLocation = download.path)
}
