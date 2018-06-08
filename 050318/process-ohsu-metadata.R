synId <- "syn8149174"
obj <- synGet(id=synId, downloadFile = TRUE)
file <- getFileLocation(obj)
tmp <- read.xlsx(file, sheet = 1)
ohsu.metadata <- tmp[, c("patient_id", "gender", "lab_id", "specimen_type")]

tmp <- read.xlsx(file, sheet = 3)
tmp <- tmp[grepl(pattern="FAB", ignore.case=TRUE, tmp$lab_type),]
tmp <- unique(tmp[, c("patient_id", "lab_result")])
colnames(tmp) <- c("patient_id", "fab")
ohsu.metadata <- merge(ohsu.metadata, tmp, by = "patient_id", all = TRUE)

patient.age.tbl <- read.xlsx(file, sheet = 5)
patient.age.tbl <- patient.age.tbl[as.numeric(patient.age.tbl$age) > 0,]

## Map lab_id to seq_ids
synId <- "syn10083817"
obj <- synGet(synId, downloadFile = TRUE)
file <- getFileLocation(obj)
ohsu.rnaseq.sample.summary <- read.table(file, header=TRUE, sep="\t")
ohsu.metadata <- merge(ohsu.metadata, ohsu.rnaseq.sample.summary[, c("Original_LabID", "SeqID")], by.x = "lab_id", by.y = "Original_LabID", all = FALSE)

if(FALSE) {
tmp <- patient.age.tbl
tmp <- tmp[, c("patient_id", "lab_id", "age")]
tmp <- ddply(unique(tmp[, c("patient_id", "lab_id", "age")]), c("patient_id", "lab_id"), .fun = function(df) { df$age <- mean(df$age); return(unique(df)) })
ohsu.metadata <- merge(ohsu.metadata, tmp[, c("patient_id", "lab_id", "age")], by = c("patient_id", "lab_id"), all = TRUE)

tmp <- patient.age.tbl
tmp <- tmp[as.numeric(tmp$age) > 0,]
tmp <- tmp[, c("patient_id", "lab_id", "response")]
ohsu.simplified.response.tbl <- ddply(tmp[!is.na(tmp$response),], c("patient_id", "lab_id"), 
                                      .fun = function(df) {
                                        if(any(grepl(x=df$response, pattern="Complete Response"))) { return("Complete Response")}
                                        if(any(df$response == "Hematologic CR")) { return("Hematologic CR")}
                                        if(any(df$response == "Refractory")) { return("Refractory")}
                                        if(any(df$response == "Supportive/Palliative Care")) { return("Supportive/Palliative Care")}
                                        if(any(df$response == "Unknown")) { return("Unknown")}
                                        warning(paste0("What is ", df$response, "\n"))
                                      })
colnames(ohsu.simplified.response.tbl) <- c("patient_id", "lab_id", "response")
ohsu.metadata <- merge(ohsu.metadata, ohsu.simplified.response.tbl, by = c("patient_id", "lab_id"), all = TRUE)

## Remove NAs
ohsu.metadata <- ohsu.metadata[!is.na(ohsu.metadata$patient_id),]

}

