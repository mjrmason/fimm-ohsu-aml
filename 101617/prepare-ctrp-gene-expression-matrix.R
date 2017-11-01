library(synapseClient)
library("biomaRt")
synapseLogin()

cat("Log-cpm transform the CTRP data\n")

## Read in the CTRP expression data
## Read in the CTRP expression data (this is counts)
## CCLE_RNAseq_081117.reads.tsv
synId <- "syn10808299"
obj <- synGet(synId, downloadFile = TRUE)
file <- getFileLocation(obj)
ctrp.cnt.expr <- read.table(file, header=TRUE, sep="\t")

library(edgeR)
## Convert expr to cpm using edger and log2-transform
log.cpm.expr <- cpm(ctrp.cnt.expr, log=TRUE)

## Samples in CTRP drug data do not have .'s and are all in caps.  So change the expression
## sample names to conform to this:
colnames(log.cpm.expr) <- gsub(x=colnames(log.cpm.expr), pattern="\\.", replacement="")
colnames(log.cpm.expr) <- toupper(colnames(log.cpm.expr))

## Store in
## FIMM_BEAT AML Collaboration/Files/CTRPv2 Data
parentId <- "syn10288724"

path <- paste0("ctrpv2-log-cpm-gene-expr.tsv")
write.table(file=path, log.cpm.expr, row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)
f <- File(path, parentId = parentId, synapseStore = TRUE)
synStore(f, executed = NULL, forceVersion = FALSE)

