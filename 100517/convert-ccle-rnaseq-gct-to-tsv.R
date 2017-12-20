library(synapseClient)
library(data.table)
synapseLogin()

cat("Convert the CCLE RNA-seq data from a GCT file to a tsv\n")
cat("This includes:\n")
cat("1. Dropping the Ensembl version: i.e., ENSGXXX.Y -> ENSGXXX\n")
cat("2. Dropping the Description column\n")
cat("3. Making the Name column the rownames\n")
cat("4. Converting the sample names\n")

## CCLE_RNAseq_081117.reads.gct was downloaded from https://portals.broadinstitute.org/ccle/data on 9/15/17

## The first line in the GCT is commented/preceeded by a '#'
## The second seems to be the dimensions of the file.
## Skip both of these.
## NB: by using fread we avoid having samples that begin with numbers having an X prepended
## (e.g., 2313287_STOMACH -> X2313287_STOMACH) as would read.table
tbl <- as.data.frame(fread("CCLE_RNAseq_081117.reads.gct", sep="\t", header=TRUE, skip=2))

## Drop the Ensembl version: i.e., ENSGXXX.Y -> ENSGXXX
tbl$Name <- unlist(lapply(tbl$Name, function(str) unlist(strsplit(str, split="\\."))[1]))

## Make the Name column the rownames
rownames(tbl) <- tbl$Name

## Drop the Description column
tbl <- tbl[, !(colnames(tbl) %in% c("Name", "Description"))]

## Convert the sample names
## Do so using the cell line annotation file (CCLE_sample_info_file_2012-10-18.txt)
## synId <- "syn7112975"
## ccle.anno <- read.table(getFileLocation(synGet(synId)), header=TRUE, sep="\t", as.is=TRUE)
## No, it seems that the CTRP/drug sample names are just the prefix of the names used here.
## e.g., 2313287_STOMACH -> 2313287 (note that the cell.line.primary.name in the annotation file is 23132/87)
colnames(tbl) <- unlist(lapply(colnames(tbl), function(str) unlist(strsplit(str, split="_"))[1]))

## Store in synapse
file <- "CCLE_RNAseq_081117.reads.tsv"
write.table(file=file, tbl, sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

## Store in "Cancer Cell Line Encyclopedia (CCLE) RNA-Seq" folder
parentId <- "syn5612998"
fileObj <- File(file, parentId = parentId)
synStore(fileObj)

cat("Successfully stored reformatted CCLE RNA-seq count data\n")