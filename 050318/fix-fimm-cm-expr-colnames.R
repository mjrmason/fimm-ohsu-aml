suppressPackageStartupMessages(library(synapseClient))

## The drug data screens have a "_CM" suffix (e.g., FH_5232_14042016_9999_BM_CM) but in the expression data 
## they do not (FH_5232_14042016_9999_BM). More importantly, 3 of the screens in the expression data lack "_BM". In expression data:

## FH_4599_19042016_1030
## FH_5232_14042016_9999
## FH_5305_18052016_9999

## Presumably the same screens in the drug data:

## FH_4599_19042016_1030_BM_CM
## FH_5232_14042016_9999_BM_CM
## FH_5305_18052016_9999_BM_CM    

synapseLogin()

## "RNASeq.CPM.log2.bc.csv"
synId <- "syn12176271"
obj <- synGet(id=synId, downloadFile = TRUE)
file <- getFileLocation(obj)
## This takes a long time to read in
expr <- read.table(file, sep=",", header=TRUE, as.is=TRUE)

missing.bm.pb <- !grepl(rownames(expr), pattern="_BM") & !grepl(rownames(expr), pattern="_PB")

if(!(all(sort(rownames(expr)[missing.bm.pb]) == sort(c("FH_4599_19042016_1030", "FH_5232_14042016_9999", "FH_5305_18052016_9999"))))) {
  stop("Was expecting a different set of samples missing the _BM or _PB suffix\n")
}

## I know from looking ahead at the drug data that these should have a _BM suffix (as opposed to PB)
rownames(expr)[missing.bm.pb] <- unlist(lapply(rownames(expr)[missing.bm.pb], function(x) paste0(x, "_BM")))

## Now add "_CM"
rownames(expr) <- unlist(lapply(rownames(expr), function(x) paste0(x, "_CM")))

## Transpose the matrix and add gene as a named column so we can read in using read.table
texpr <- as.data.frame(t(expr))
cols <- colnames(texpr)
##texpr$gene <- rownames(texpr)
##texpr <- texpr[, c("gene", cols)]

ofile <- "RNASeq.CPM.log2.bc.fixed.names.csv"
ofile <- "RNASeq.CPM.log2.bc.fixed.names.tsv"
output.dir <- "output"
dir.create(output.dir)
write.table(file = paste0(output.dir, "/", ofile), texpr, sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

parentId <- "syn12176270"
f <- File(paste0(output.dir, "/", ofile), parentId = parentId, synapseStore = TRUE)
synStore(f, executed = NULL, forceVersion = FALSE)
