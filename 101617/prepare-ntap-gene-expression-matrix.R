library(synapseClient)
library("biomaRt")
synapseLogin()

cat("Read in the kallisto-derived estimated transcript-level expression counts and sum to get gene-level results\n")

## Read in the estimated counts at the transcript level
synId <- "syn5562376"
obj <- synGet(synId, downloadFile = TRUE)
trnscript.expr.matrix <- read.table(getFileLocation(obj), sep="\t", header=TRUE)
trnscript.expr.matrix <- as.matrix(trnscript.expr.matrix)
transcripts <- rownames(trnscript.expr.matrix)
transcripts <- unlist(lapply(transcripts, function(tr) paste0("ENST", gsub(x=tr, pattern=".*ENST(.+)", replacement="\\1"))))
rownames(trnscript.expr.matrix) <- transcripts

## These are in grch38
ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
gen <- getBM(attributes=c('ensembl_gene_id', 'ensembl_transcript_id'), filters = 'ensembl_transcript_id', values = transcripts, mart = ensembl)

rownames(gen) <- gen$ensembl_transcript_id
flag <- rownames(trnscript.expr.matrix) %in% rownames(gen)
trnscript.expr.matrix <- trnscript.expr.matrix[flag, ]
rownames(trnscript.expr.matrix) <- gen[rownames(trnscript.expr.matrix),"ensembl_gene_id"]


gene.expr.matrix <- round(rowsum(trnscript.expr.matrix, rownames(trnscript.expr.matrix)))

library(edgeR)
log.cpm.expr <- cpm(gene.expr.matrix, log=TRUE)

## The sample names are synapse ids -- translate them to sampleIdentifiers using the table holding
## the mapping between sample names and file names
synId <- "syn8397154"
ntap.sample.key <- get_table_df(synId)
ntap.sample.key <- ntap.sample.key[!is.na(ntap.sample.key$`RNASeq Data`),]

flag <- colnames(log.cpm.expr) %in% ntap.sample.key$`RNASeq Data`
if(all(flag)) {
  cat("All samples in expression matrix are listed in table\n")
} else {
  cat("The following samples in the expression matrix are not listed in the table\n")
  print(colnames(log.cpm.expr)[!flag])
}

map <- ntap.sample.key$`Sample Name`
names(map) <- ntap.sample.key$`RNASeq Data`
colnames(log.cpm.expr) <- as.vector(map[colnames(log.cpm.expr)])

## Store in
## Files/Analysis/Drug Sensitivity/NCATS + CTRP Drug Response Modeling (syn11244429)
parentId <- "syn11244429"

path <- paste0("ntap-log-cpm-gene-expr.tsv")
write.table(file=path, log.cpm.expr, row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)
f <- File(path, parentId = parentId, synapseStore = TRUE)
synStore(f, executed = NULL, forceVersion = FALSE)

