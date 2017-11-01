bioc.packages <- c("nplr", "openxlsx", "drc", "minpack.lm", "maftools", "vcd", "binom", "DMwR", "WGCNA", "VennDiagram")

install.bioc.packages.auto <- function(x) { 
  x <- as.character(x)
  print(x)
##  if(isTRUE(x %in% .packages(all.available=TRUE))) { 
##    eval(parse(text = sprintf("require(\"%s\")", x)))
##  } else { 
##    #update.packages(ask= FALSE) #update installed packages.
##    eval(parse(text = sprintf("install.packages(\"%s\", dependencies = TRUE)", x)))
##  }

  if(isTRUE(x %in% .packages(all.available=TRUE))) { 
    eval(parse(text = sprintf("require(\"%s\")", x)))
  } else {
    source("http://bioconductor.org/biocLite.R")
    #biocLite(character(), ask=FALSE) #update installed packages.
    eval(parse(text = sprintf("biocLite(\"%s\")", x)))
    eval(parse(text = sprintf("require(\"%s\")", x)))
  }
}

d <- unlist(lapply(bioc.packages, install.bioc.packages.auto))

## Set the download and output directory
download.path <- "../input/"
if (!file.exists(download.path)) {
  dir.create(download.path)
}

output.path <- "output/"
if (!file.exists(output.path)) {
  dir.create(output.path)
}

beat.aml.cache <- "../input/beat-aml-files/"

suppressPackageStartupMessages(library("nplr"))
suppressPackageStartupMessages(library("foreach"))
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(tidyr))

suppressPackageStartupMessages(library("drc"))
suppressPackageStartupMessages(library("minpack.lm"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("gridExtra"))
suppressPackageStartupMessages(library("caTools"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("maftools"))
suppressPackageStartupMessages(library("pcaMethods"))
## suppressPackageStartupMessages(library("GGally"))
suppressPackageStartupMessages(library("VennDiagram"))
suppressPackageStartupMessages(library(ReporteRs))
suppressPackageStartupMessages(library(rJava))

library(ggbeeswarm)
library(plyr)
library(stringr)
library(synapseClient)

synapseLogin()

## Set the java heap to 8GB
options(java.parameters = "-Xmx8000m")

## Define a java garbage collection function
jgc <- function()
{
  gc()
  .jcall("java/lang/System", method = "gc")
} 

cat("Translate lab ids (used in drug data) to seq ids (used in sequence data) in OHSU\n")

## Read in the rna-seq summary table, which provides a map from the lab_id's used
## in the drug response data to the SeqID's used in the expression data.
synId <- "syn10083817"
obj <- synGet(synId, downloadFile = TRUE)
file <- getFileLocation(obj)
ohsu.rnaseq.sample.summary <- read.table(file, header=TRUE, sep="\t")
lab.id.to.seq.id.tbl <- ohsu.rnaseq.sample.summary[, c("Original_LabID", "SeqID")]
colnames(lab.id.to.seq.id.tbl) <- c("lab_id", "SeqID")

## parentId is the OHSU processed data folder
parentId <- "syn10083332"
tbl <- synQuery(paste0("SELECT name,id FROM file WHERE parentId==\'", parentId, "\'"))

## For each of the ohsu.common.drugs.dss.*.tsv files, merge the lab_id in that file
## with the Original_LabID in the summary file.  The SeqID in the summary file will be 
## the id used in the expression data.

## parentId is the FIMM processed data folder
parentId <- "syn10083332"
tbl <- synQuery(paste0("SELECT name,id FROM file WHERE parentId==\'", parentId, "\'"))
tbl <- tbl[grepl(tbl$file.name, pattern="ohsu.common.drugs.dss"),]
for(i in 1:nrow(tbl)) {
  file.name <- tbl$file.name[i]
  file.id <- tbl$file.id[i]
  cat(paste0("Adding SeqID column to ", file.name, "\n"))
  obj <- synGet(file.id, downloadFile = TRUE)
  file <- getFileLocation(obj)
  ohsu.dss <- as.data.frame(fread(file))

  ## Drop the SeqID column, if it already exists
  ohsu.dss <- ohsu.dss[, !(colnames(ohsu.dss) %in% "SeqID")]
  cat(paste0("   Original num columns: ", nrow(ohsu.dss), "\n"))
  ohsu.dss <- merge(ohsu.dss, lab.id.to.seq.id.tbl, all.x = TRUE, by = "lab_id")
  cat(paste0("   Final num columns: ", nrow(ohsu.dss), "\n"))
  write.table(file = file.name, ohsu.dss, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
  f <- File(file.name, parentId = parentId, synapseStore = TRUE)
  synStore(f, forceVersion = FALSE)
}


