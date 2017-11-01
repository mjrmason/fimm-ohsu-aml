suppressPackageStartupMessages(library("synapseClient"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("foreach"))
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("gplots"))
suppressPackageStartupMessages(library("openxlsx"))

synapseLogin()

## Read in the FIMM expression data (RNASeq.CPM.log2.bc.csv)
## Load (batch corrected) FIMM expression data
obj <- synGet(id="syn8270602", downloadFile = TRUE)
file <- getFileLocation(obj)
fimm.expr <- read.table(file, header=TRUE, sep=",")
## Make the columns samples
fimm.expr <- t(fimm.expr)

## Read in the FIMM genomic data
## The "bin" file simply is 0/1 matrix indicating whether the sample (row) has a mutation
## in the gene (column).  These use gene symbols.
## NB: a 1 indicates that the gene had a non-synonymous or splice site mutation.
## path <- "fimm.genomic.ns.ss.tsv"
## Columns correspond to genes.
synId <- "syn10220160"
obj <- synGet(id=synId, downloadFile = TRUE)
file <- getFileLocation(obj)
fimm.genomic <- read.table(file, header=TRUE, sep="\t")
## Make the columns samples
fimm.genomic <- t(fimm.genomic)

## Read in the FIMM Metadata and ensure that we only analyze AML data
synId <- "syn8270594"
obj <- synGet(id=synId, downloadFile = TRUE)
file <- getFileLocation(obj)
fimm.metadata <- read.table(file, header=TRUE, sep=",")
fimm.metadata <- subset(fimm.metadata, grepl(diagnosis, pattern="AML"))

inh.patient.ids <- unique(fimm.metadata$person)
rna.patient.ids <- unique(unlist(lapply(colnames(fimm.expr), function(str) {
                                                               re <- regexpr(text=str, pattern="^[A-z]+_\\d+")
                                                               substr(str, re[1], re[1] + attr(re, "match.length")-1) })))
rna.patient.ids <- gsub(rna.patient.ids, pattern="_", replacement=".")
dna.patient.ids <- unique(unlist(lapply(colnames(fimm.genomic), function(str) {
                                                               re <- regexpr(text=str, pattern="^[A-z]+_\\d+")
                                                               substr(str, re[1], re[1] + attr(re, "match.length")-1) })))
dna.patient.ids <- gsub(dna.patient.ids, pattern="_", replacement=".")
aml.patient.ids <- fimm.metadata$person[grepl(fimm.metadata$diagnosis, pattern="AML")]
vennList = list("DRC"=inh.patient.ids,"RNA-seq"= rna.patient.ids,"WES"=dna.patient.ids, "AML"=aml.patient.ids)

## Read in FIMM raw drug response data (FIMM_DSRT_DATA.xls)

synId <- "syn8488824"
obj <- synGet(id=synId, downloadFile = TRUE, downloadLocation = ".")
file <- getFileLocation(obj)
fimm.raw.dr <- read.xlsx(file)
cat("\nNumber drugs tested: ", length(unique(fimm.raw.dr$DRUG_ID)), "\n")

tiff("fimm-venn-drc-rna-dna.tiff")
venn(vennList)
title(main="FIMM Assay Overlap")
d <- dev.off()
