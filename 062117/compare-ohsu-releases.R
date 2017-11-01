
library(synapseClient)
library(data.frame)

synapseLogin()

# Read in the inhibitor_data_points_2017_01_12.txt inhibitor data
synId <- "syn8149180"
inh.obj <- synGet(synId, downloadFile=TRUE, downloadLocation = download.path)
ohsu.inh.1.12.tbl <- fread(getFileLocation(inh.obj))
ohsu.inh.1.12.tbl <- as.data.frame(ohsu.inh.1.12.tbl)

## Read in the inhibitor_data_points_2017_05_17.txt inhibitor data
synId <- "syn10008778"
inh.obj <- synGet(synId, downloadFile=TRUE, downloadLocation = download.path)
ohsu.inh.5.17.tbl <- fread(getFileLocation(inh.obj))
ohsu.inh.5.17.tbl <- as.data.frame(ohsu.inh.5.17.tbl)

## Are there any drugs unique 2017_01_12 data
df <- data.frame(inhibitor = unique(ohsu.inh.1.12.tbl$inhibitor)[!(unique(ohsu.inh.1.12.tbl$inhibitor) %in% unique(ohsu.inh.5.17.tbl$inhibitor))])
write.table(file="inhibitors-missing-from-inh-5-17.tsv", df, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

df <- data.frame(patient_id = unique(ohsu.inh.1.12.tbl$patient_id)[!(unique(ohsu.inh.1.12.tbl$patient_id) %in% unique(ohsu.inh.5.17.tbl$patient_id))])
write.table(file="patient_ids-missing-from-inh-5-17.tsv", df, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

df <- data.frame(lab_id = unique(ohsu.inh.1.12.tbl$lab_id)[!(unique(ohsu.inh.1.12.tbl$lab_id) %in% unique(ohsu.inh.5.17.tbl$lab_id))])
write.table(file="lab_ids-missing-from-inh-5-17.tsv", df, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
