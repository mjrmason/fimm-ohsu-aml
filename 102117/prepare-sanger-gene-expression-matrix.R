library(synapseClient)
library(data.table)
synapseLogin()

cat("Read in Sanger RMA-normalized (which is in log space) Ensembl expression matrix.")

## The only change we will make is using the ensembl_gene column as the rownames and
## droppping it from the matrix.

## Read in the estimated counts at the transcript level
synId <- "syn11275816"  ## File: sanger1018_brainarray_ensemblgene_rma.txt
obj <- synGet(synId, downloadFile = TRUE)
expr.matrix <- as.data.frame(fread(getFileLocation(obj), header=TRUE))
print(head(expr.matrix[, 1:5]))

rownames(expr.matrix) <- expr.matrix$ensembl_gene
expr.matrix <- expr.matrix[, !(colnames(expr.matrix) %in% "ensembl_gene")]

print(head(expr.matrix[, 1:5]))

## Store in
## FIMM_BEAT AML Collaboration/Files/Sanger Data
parentId <- "syn11287866"

path <- paste0("sanger-rma-gene-expr.tsv")
write.table(file=path, expr.matrix, row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)
f <- File(path, parentId = parentId, synapseStore = TRUE)
synStore(f, executed = NULL, forceVersion = FALSE)

cat("Successfully stored gene expression matrix\n")