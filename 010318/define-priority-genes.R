rm(list = ls())

suppressPackageStartupMessages(library(glmnet))
suppressPackageStartupMessages(library("randomForest"))
suppressPackageStartupMessages(library("e1071"))
suppressPackageStartupMessages(library("foreach"))
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(tidyr))

suppressPackageStartupMessages(library("minpack.lm"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("gridExtra"))
suppressPackageStartupMessages(library("caTools"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("pcaMethods"))
suppressPackageStartupMessages(library("VennDiagram"))
suppressPackageStartupMessages(library("caret"))

suppressPackageStartupMessages(library(ggbeeswarm))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(synapseClient))

suppressPackageStartupMessages(library("GSVA"))
# use GSA to read in a gmt file
suppressPackageStartupMessages(library("GSA"))

suppressPackageStartupMessages(library(mygene))
suppressPackageStartupMessages(library(gplots))

## Log in to Synapse
synapseLogin()

## Bring in common code
source("../common/data-preprocessing.R")
source("../common/models.R")
source("../common/plotting.R")
source("../common/utils.R")
source("../common/drug-mapping-functions/mapping_helpers.R")
source("../common/drug-mapping-functions/drug_mapping.R")

## Register parallelization
num.cores <- detectCores()
if(!is.na(num.cores) && (num.cores > 1)) {
  suppressPackageStartupMessages(library("doMC"))
  cat(paste("Registering ", num.cores-1, " cores.\n", sep=""))
  registerDoMC(cores=(num.cores-1))
}
num.processes <- num.cores - 1

## Begin setup (this will need to be changed as the data sets in the intersection change)

## Purpose:
## This file models 2 data sets (OHSU and FIMM)
source("fimm-ohsu-setup-120817.R")

if(TRUE) {
## Overwrite the expression synIds to exclude outliers (as before), but _not_ to restrict to highly variable/cancer genes
## Batch-corrected log2 cpm (or rma) expression files (NB: these exclude 2 outliers in OHSU)
data.set.expr.synIds <- list("ohsu" = "syn11362256", "fimm" = "syn11362257")
## names(data.set.expr.synIds) <- data.sets
expr.folder.synIds <- list("ohsu" = "syn11361089", "fimm" = "syn11361089")
data.set.expr.files <- list("ohsu" = "ohsu-expr-fimm-ohsu-outlier1-combat.tsv", "fimm" = "fimm-expr-fimm-ohsu-outlier1-combat.tsv")
}

file.prefix <- "fimm-ohsu-overlap-120817"
rdata.file <- ".Rdata.overlap.120817"

## End setup

cat(paste0("Model AUC ~ expr for data sets: ", paste(data.sets, collapse=","), "\n"))

## Load in the fit files
fits <- llply(data.set.dss.fit.synIds, .parallel = FALSE,
              .fun = function(synId) {
                       obj <- synGet(synId, downloadFile = TRUE)
                       file <- getFileLocation(obj)
                       read.table(file, sep="\t", header=TRUE, as.is=TRUE, comment.char="", quote="\"")
              })

## Merge metadata to fit files
cat("Merging metadata to fits\n")
for(ds in names(fits)) {
  switch(ds,
    "ohsu" = {
      ## Add the SeqID to the OHSU drug data.
      ## Read in the rna-seq summary table, which provides a map from the lab_id's used
      ## in the drug response data to the SeqID's used in the expression data.
      synId <- "syn10083817"
      obj <- synGet(synId, downloadFile = TRUE)
      file <- getFileLocation(obj)
      ohsu.rnaseq.sample.summary <- read.table(file, header=TRUE, sep="\t")
## Looks like the SeqID has already been included in the fit file.
##      fits[[ds]] <- merge(fits[[ds]], ohsu.rnaseq.sample.summary, by.x = "lab_id", by.y = "Original_LabID")
    },
    "fimm" = {
    },
    "ctrp" = {
      fits[[ds]] <- merge(fits[[ds]], ctrp.cell.line.metadata[, c("master_ccl_id", "ccl_name")])
    },
    "ntap" = {
    },
    "sanger" = {
    },
    {
      stop("Unrecognized data set\n")
    })
}

## Drop any fits that are flagged as outliers
fct.postfixes <- list("LL.4" = ".ll4", "L.4" = ".l4")
fct <- "LL.4"
for(ds in names(fits)) {
    exclude.cols <- colnames(fits[[ds]])[grepl(pattern=paste0("exclude", fct.postfixes[[fct]]), x=colnames(fits[[ds]]))]
    any.excluded <- unlist(apply(fits[[ds]][, exclude.cols], 1, function(row) any(row)))
    orig.nrow <- nrow(fits[[ds]])
    fits[[ds]] <- fits[[ds]][!any.excluded, ]
    new.nrow <- nrow(fits[[ds]])
    cat(paste0(ds, ": filtered ", orig.nrow - new.nrow, " of the original ", orig.nrow, " fits.\n"))
}

## Load the map of drugs between data sets
obj <- synGet(drug.mapping.synId, downloadFile = TRUE)
drug.map <- read.table(getFileLocation(obj), sep="\t", header=TRUE, quote="\"", comment.char="")

## Drop any drugs that show up in multiple rows of the drug map--which would indicate an ambiguous mapping.
my.dup <- function(x) duplicated(x, fromLast = TRUE) | duplicated(x, fromLast = FALSE)
for(col in drug.map.drug.id.cols) {
  flag <- !my.dup(drug.map[, col])
  drug.map <- drug.map[flag, ]
}

## Restrict the fits to those drugs in common across the data sets -- do so by merging.
cat("Restricting fits to those involving common drugs\n")
for(ds in names(fits)) {
  fits[[ds]] <- merge(fits[[ds]], drug.map, by.x = data.set.drug.id.cols[[ds]], by.y = drug.map.drug.id.cols[[ds]], all = FALSE)
}

common.drugs <- Reduce("intersect", lapply(fits, function(fit) unique(fit[, drug.name.col])))
for(ds in names(fits)) {
  flag <- fits[[ds]][, drug.name.col] %in% common.drugs
  fits[[ds]] <- fits[[ds]][flag, ]
}

## Limit the drug map to common drugs (which, by design of the drug map file, should already be the case)
drug.name.tbl <- drug.map[drug.map[, drug.name.col] %in% common.drugs,]

cat("Translating symbols\n")

drug.name.tbl$Ensg.Targets <- unlist(lapply(drug.name.tbl$Gene.Targets,
                                            function(str) {
                                              if(is.na(str) || (str == "")) { return(NA) }
                                              symbols <- unlist(strsplit(as.character(str), ",[ ]*"))
                                              ensg.genes <- symbols.to.ensg.mapping(symbols)
                                              ensg.genes <- ensg.genes$ensg[!is.na(ensg.genes$ensg)]
                                              if(length(ensg.genes) == 0) { return(NA) }
                                              paste(ensg.genes, collapse=", ")
                                            }))
cat("Done translating symbols\n")

load("clusters.of.consistent.drugs.Rd")

## Translate clusters of consistent drugs from OHSU names to FIMM names
tbl <- drug.name.tbl[, c("OHSU_DRUG_NAME", "FIMM_DRUG_NAME")]
tbl <- unique(tbl)
rownames(tbl) <- tbl$OHSU_DRUG_NAME

fimm.clusters.of.consistent.drugs <- llply(clusters.of.consistent.drugs, .fun = function(cluster) as.character(tbl[cluster,"FIMM_DRUG_NAME"]))

common.drugs.by.data.set <- list()
drugs.by.data.set <- list()
for(ds in names(fits)) {
  common.drugs.by.data.set[[ds]] <- drug.name.tbl[, drug.name.col]
  drugs.by.data.set[[ds]] <- fimm.clusters.of.consistent.drugs
}

## Load in the expr files
exprs <- llply(data.set.expr.synIds, .parallel = FALSE,
              .fun = function(synId) {
                       obj <- synGet(synId, downloadFile = TRUE)
                       file <- getFileLocation(obj)
                       read.table(file, sep="\t", header=TRUE, as.is=TRUE, comment.char="", quote="\"")
              })

## Sanitize some of the column headers 
for(nm in names(exprs)) {
  switch(nm,
    "ohsu" = {
      ## Convert the column names from X20.00347 to 20-00347
      colnames(exprs[[nm]]) <- gsub("X(.*)\\.(.*)", "\\1-\\2", colnames(exprs[[nm]]))
    },
    "fimm" = {
    },
    "ctrp" = {
    },
    "ntap" = {
    },
    "sanger" = {
      colnames(exprs[[nm]]) <- gsub("X(.*)", "\\1", colnames(exprs[[nm]]))
    },
    {
      stop("Unrecognized data set\n")
    })
}

common.genes <- Reduce("intersect", lapply(exprs, rownames))

use.subset <- TRUE
if(use.subset) {
## This is file supplementary data 4
aml.logsdon.gene.tbl <- openxlsx:::read.xlsx("../common-resources/41467_2017_2465_MOESM7_ESM.xlsx", sheet = 1, startRow = 3)
aml.logsdon.gene.tbl <- aml.logsdon.gene.tbl[order(as.numeric(aml.logsdon.gene.tbl$SCORE), decreasing=TRUE),]
## Take the top 2000, even though the plateau in scores is around 3000
aml.logsdon.genes <- unique(aml.logsdon.gene.tbl$GENE[1:500])
aml.logsdon.genes <- symbols.to.ensg.mapping(aml.logsdon.genes)
aml.logsdon.genes <- na.omit(aml.logsdon.genes$ensg)

## Take all Cancer genes, not just those with Cancer Type == LAML or confidence = A (High), B (Medium), or C (low)
cancer.genes <- openxlsx::read.xlsx("../common-resources/TableS2A.xlsx", sheet = 1, startRow = 3)
cancer.genes <- unique(cancer.genes$Gene)
cancer.genes <- symbols.to.ensg.mapping(cancer.genes)
cancer.genes <- na.omit(cancer.genes$ensg)

## From Suleiman and Ammad
genes <- read.table("../common-resources/gene_names.tsv", sep="\t", header=FALSE, stringsAsFactors = FALSE)$V1
gene.subset <- genes[grepl(genes, pattern="ENSG")]
genes.to.translate <- genes[!grepl(genes, pattern="ENSG")]
trns <- symbols.to.ensg.mapping(genes.to.translate)
gene.subset <- unique(c(gene.subset, as.character(trns$ensg)))

fimm.genes.tbl <- data.frame(ensg = gene.subset, FIMM = TRUE)

## Include target genes as annotated by FIMM 

## Include target genes as provided by Robert in his database.
## Load in the drug structure annotations provided by Robert
drug.structures.tbl <- read.table("../common/drug-mapping-functions/drug_structures.txt", comment.char="", sep="\t", header=TRUE, as.is=TRUE)
cat("Calling getGenesFromDrugs\n")
foo <- getGenesfromDrugs(drug.structures.tbl$drug, drug.structures.tbl$smiles, tanimoto_threshold = 0.95, parallelized = F)
cat("Done calling getGenesFromDrugs\n")
evo.ensg.targets <- symbols.to.ensg.mapping(foo$hugo_gene)
evo.ensg.targets <- evo.ensg.targets$ensg

save.image("evo.Rd")

old.nrow <- nrow(drug.name.tbl)
m <- merge(drug.name.tbl, drug.structures.tbl, by.x = "FIMM_DRUG_NAME", by.y = "drug")
print(nrow(m))
if(old.nrow != nrow(m)) {
  stop("Size of drug name table changed when merging in structure info\n")
}
drug.name.tbl <- m

print(length(aml.logsdon.genes))
print(length(cancer.genes))
print(length(gene.subset))

cat("Adding Suleiman/Ammad genes\n")
cat(paste0("Length of genes: ", length(gene.subset), "\n"))

cat("Adding cancer genes\n")
gene.subset <- unique(union(cancer.genes, gene.subset))

intogen.genes.tbl <- data.frame(ensg = cancer.genes, intOGen = TRUE)
## fimm.genes.tbl <- data.frame(ensg = gene.subset, FIMM = TRUE)
logsdon.genes.tbl <- data.frame(ensg = aml.logsdon.genes, logsdon = TRUE)

cat("Adding Logsdon genes\n")
gene.subset <- unique(union(aml.logsdon.genes, gene.subset))

## Include target genes as annotated by FIMM 
all.ensg.targets <- unlist(llply(drug.name.tbl$Ensg.Targets, .fun = function(str) unlist(strsplit(str, split=",[ ]*"))))
all.ensg.targets <- na.omit(all.ensg.targets)
all.ensg.targets <- unique(all.ensg.targets)

fimm.targets.tbl <- data.frame(ensg = all.ensg.targets, FIMM.targets = TRUE)

cat(paste0("Length of genes: ", length(gene.subset), "\n"))
cat("Adding FIMM-annotated gene targets\n")
gene.subset <- unique(union(gene.subset, all.ensg.targets))
cat(paste0("Length of genes: ", length(gene.subset), "\n"))

cat(paste0("Adding evo-annotated gene targets\n"))
gene.subset <- unique(union(gene.subset, evo.ensg.targets))
cat(paste0("Length of genes: ", length(gene.subset), "\n"))

evo.targets.tbl <- data.frame(ensg = evo.ensg.targets, evo.targets = TRUE)

m <- merge(unique(intogen.genes.tbl), unique(fimm.genes.tbl), all = TRUE)
m <- merge(m, unique(logsdon.genes.tbl), all = TRUE)
m <- merge(m, unique(fimm.targets.tbl), all = TRUE)
m <- merge(m, unique(evo.targets.tbl), all = TRUE)

write.table(file = "prioritized-genes.tsv", m, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")


