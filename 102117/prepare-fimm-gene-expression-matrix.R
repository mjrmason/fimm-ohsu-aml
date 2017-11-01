library(synapseClient)
library(data.table)
synapseLogin()

cat("Prepare FIMM Ensembl expression matrix in log2 space")

## The only change we make is to transpose to make columns the samples

## Read in the FIMM expression data (RNASeq.CPM.log2.bc.csv)
## Load (batch corrected) FIMM expression data
obj <- synGet(id="syn8270602", downloadFile = TRUE)
file <- getFileLocation(obj)
expr.matrix <- read.table(file, header=TRUE, sep=",")
## Make the columns samples
expr.matrix <- t(expr.matrix)

## Store in
## FIMM_BEAT AML Collaboration/Files/Fimm Data/Processed DAta
parentId <- "syn8270577"

path <- paste0("fimm-log-cpm-gene-expr.tsv")
write.table(file=path, expr.matrix, row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)
f <- File(path, parentId = parentId, synapseStore = TRUE)
synStore(f, executed = NULL, forceVersion = FALSE)

cat("Successfully stored gene expression matrix\n")