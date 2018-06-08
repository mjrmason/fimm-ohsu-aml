library(dplyr)
library(plyr)

manifest.file <- "MANIFEST.txt"
clinical.file <- "clinical.tsv"
data.root.dir <- "tcga-aml"

manifest <- read.table(manifest.file, sep="\t", header=TRUE, as.is=TRUE)

manifest <- manifest[grepl(manifest$filename, pattern="count"),]

files <- manifest$filename
names(files) <- manifest$id
indices <- 1:nrow(manifest)
names(indices) <- manifest$id

cnts <- llply(indices, .parallel = FALSE,
              .fun = function(i) {
                       gz.file <- manifest$filename[i]
                       print(gz.file)
                       cnt.file <- paste0(data.root.dir, "/", gsub(gz.file, pattern=".gz", replacement=""))
                       anno.file <- gsub(x = cnt.file, pattern = "([^//]+$)", replacement = "annotations.txt")
                       if(!file.exists(cnt.file)) { system(paste0("gunzip ", data.root.dir, "/", gz.file)) }
                       cnts <- read.table(cnt.file, sep="\t", header=FALSE)
                       anno <- read.table(anno.file, sep="\t", header=TRUE, as.is = TRUE)
                       id <- anno$entity_id[1]
                       colnames(cnts) <- c("gene", id)
                       cnts
              })

cat("Combining cnts\n")
cnt.tbl <- Reduce(function(x, y) merge(x, y, by = "gene"), cnts)

## Drop the version suffixes from the ensembl ids
cnt.tbl <- cnt.tbl[grepl(cnt.tbl$gene, pattern="ENSG"),]
cnt.tbl$gene <- as.character(cnt.tbl$gene)
cnt.tbl$gene <- gsub(x = cnt.tbl$gene, pattern="(.*)\\.\\d+", replacement="\\1")
rownames(cnt.tbl) <- cnt.tbl$gene
cnt.tbl <- cnt.tbl[, !(colnames(cnt.tbl) == "gene")]

## Convert the column ids to TCGA ids
clinical.tbl <- read.table(clinical.file, sep="\t", header=TRUE, as.is=TRUE)

clinical.tbl$tcga_id <- gsub(x = clinical.tbl$submitter_id, pattern="TCGA-AB-", replacement="")

rownames(clinical.tbl) <- clinical.tbl$case_id
ids <- intersect(colnames(cnt.tbl), rownames(clinical.tbl))
clinical.tbl <- clinical.tbl[ids,]
cnt.tbl <- cnt.tbl[, ids]
colnames(cnt.tbl) <- clinical.tbl$submitter_id

write.table(file = "aml-tcga-counts.tsv", cnt.tbl, sep="\t", row.names = TRUE, col.names = TRUE, quote = FALSE)