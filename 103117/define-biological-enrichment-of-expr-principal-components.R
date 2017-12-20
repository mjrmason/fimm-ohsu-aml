## rm(list = ls())

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
suppressPackageStartupMessages(library(goseq))

suppressPackageStartupMessages(library("GSVA"))
# use GSA to read in a gmt file
suppressPackageStartupMessages(library("GSA"))

## Bring in common code
source("../common/data-preprocessing.R")
source("../common/models.R")
source("../common/plotting.R")

## Register parallelization
num.cores <- detectCores()
if(!is.na(num.cores) && (num.cores > 1)) {
  suppressPackageStartupMessages(library("doMC"))
  cat(paste("Registering ", num.cores-1, " cores.\n", sep=""))
  registerDoMC(cores=(num.cores-1))
}
num.processes <- num.cores - 1

## Log in to Synapse
synapseLogin()

## Begin setup (this will need to be changed as the data sets in the intersection change)

## Purpose:
## This file models 2 data sets (OHSU and FIMM)

## Data sets to be analyzed:
data.sets <- c("ohsu", "fimm")
data.sets.to.analyze <- data.sets

## The drug mapping file lists the drugs in common across all data sets with the drug.map.drug.id.cols
## providing the identifier for each drug for each respective data set.
## The drug mapping file may be created by a script such as 
## harmonize-fimm-ohsu-compounds.R
drug.mapping.synId <- "syn11361110" ## fimm-ohsu-drug-map.tsv
drug.map.drug.id.cols <- c("OHSU_DRUG_NAME", "FIMM_Batch_ID_Drug")
names(drug.map.drug.id.cols) <- data.sets

## The column name within the drug map file that will be used to identify drugs (across all data sets)
drug.name.col <- "FIMM_DRUG_NAME"

## These fit files are assumed to be for AUCs/DSS values calculated across shared concentration ranges across data sets.
## These files may be created by a script such as
## calculate-dss-fimm-ohsu.R
data.set.dss.fit.synIds <- c("syn11362885", "syn11362891")
names(data.set.dss.fit.synIds) <- data.sets
data.set.dss.fit.files <- c("ohsu-fimm-ohsu-outlier1-dss-t0.tsv", "fimm-fimm-ohsu-outlier1-dss-t0.tsv")

## Batch-corrected log2 cpm (or rma) expression files (NB: these exclude 2 outliers in OHSU)
## These files may be created by a script such as
## normalize-fimm-ohsu-expression-distributions.R
## Batch-corrected log2 cpm (or rma) expression files (NB: these exclude 2 outliers in OHSU)
data.set.expr.synIds <- c("syn11362256", "syn11362257")
names(data.set.expr.synIds) <- data.sets
data.set.expr.files <- c("ohsu-expr-fimm-ohsu-outlier1-combat.tsv", "fimm-expr-fimm-ohsu-outlier1-combat.tsv")

## The columns within each respective data set that identifies the drug (in the name space of that data set)
data.set.drug.id.cols <- list("ohsu" = "inhibitor", "fimm" = "DRUG_ID")

## The columns within each respective data set that identifies a screen (i.e., a particular replicate of a patient/cell line x drug)
screen.id.cols <- list("ohsu" = c("patient_id", "inhibitor", "lab_id", "SeqID", "replicant", "time_of_read"),
                       "fimm" = c("SCREEN_ID", "PATIENT_ID", "DRUG_ID"))

## The columns within each respective data set that identifies the patient/cell line/specimen in the drug fit files.
## These patients should match the patient identifiers/columns of the expression data.
## NB: in some cases, these columns may not exist in the drug fit files, but need to be merged in.
## In particular, "SeqID" (the identifier used in sequencing experiments) may not be included in the OHSU drug fit files.
## Instead, patients are identified by "lab_id"'s in the drug fit file.  Below, we merge the OHSU metadata/sample
## summary information that maps lab_id's to SeqID's.
patient.id.cols <- list("ohsu" = "SeqID", "fimm" = "SCREEN_ID")

## The columns within the drug fit files that are effectively annotations and are not dependent on the fit [e.g., using
## a 4-parameter logistic L.4 or 4-parameter log logistic LL.4 fit].
merge.cols <- list(
  "ohsu" = c(screen.id.cols[["ohsu"]], "min.conc", "max.conc", "all.concs.nM", "uniq.concs.nM"),
  "fimm" = c(screen.id.cols[["fimm"]], "min.conc", "max.conc", "all.concs.nM", "uniq.concs.nM"))

## Pull in the CTRP cell line metadata so we can see the cancer/histology type of each sample
synId <- "syn5632192"
obj <- synGet(synId, downloadFile = TRUE)
ctrp.cell.line.metadata <- read.table(getFileLocation(obj), sep="\t", header=TRUE, as.is=TRUE, comment.char="", quote="")

## A string that will be included in any output files.
## Here "foc-aml-s-aml" = "fimm-ohsu-ctrp-aml-s-aml"
file.prefix <- "fimm-ohsu-outlier-ds"

## End setup

cat(paste0("Look at biological enrichment of first principal component of batch-corrected expression across: ", paste(data.sets, collapse=","), "\n"))

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

## Restrict the matrices to those genes in common across data sets
for(i in 1:length(exprs)) {
  ds <- names(exprs)[i]
  exprs[[ds]] <- exprs[[ds]][common.genes, ]
}

## Download biocarta, kegg, reactome, and hallmark data sets that have been translated to ensembl ids

## Hallmark gene sets: The list in this R object is called hallmark.gene.sets
synId <- "syn10644516"
load(getFileLocation(synGet(synId)))

## KEGG gene sets: The list in this R object is called kegg.gene.sets
synId <- "syn10644519"
load(getFileLocation(synGet(synId)))

## Biocarta gene sets: The list in this R object is called biocarta.gene.sets
synId <- "syn10644517"
load(getFileLocation(synGet(synId)))

## Reactome gene sets: The list in this R object is called reactome.gene.sets
synId <- "syn10644520"
load(getFileLocation(synGet(synId)))

all.expr <- do.call("cbind", exprs)

all.pca <- prcomp(t(all.expr))

batch.data.sets <- unlist(lapply(names(exprs), function(nm) rep(nm, ncol(exprs[[nm]]))))
cols <- as.vector(batch.data.sets)
cols <- as.factor(cols)
g1 <- varplot(all.pca, cols=cols, Main = "Uncorrected (w/ Outliers)")
g1 <- g1 + coord_fixed()
print(g1)

cos.angle <- function(x,y){
  dot.prod <- x%*%y 
  norm.x <- norm(x,type="2")
  norm.y <- norm(y,type="2")
  as.numeric(dot.prod / (norm.x * norm.y))
}

## Look at projection of a vector of reals on to a gene set
project.real.vector.onto.gene.set <- function(vec, gene.set) {
  genes <- intersect(names(vec), gene.set)
  vec <- as.vector(vec[genes])
  gene.vec <- rep(1, length(vec))
  cos.angle(vec, gene.vec)
}

calculate.pval.of.projection.of.real.vector.onto.gene.set <- function(vec, gene.set) {
  val <- project.real.vector.onto.gene.set(vec, gene.set)
  num.iters <- 10
  num.intersect <- length(intersect(names(vec), gene.set))
  vals <- unlist(llply(1:num.iters, .parallel = FALSE,
                       .fun = function(i) {
                                rand.genes <- names(vec)[sample.int(n = length(vec), size = num.intersect)]
                                project.real.vector.onto.gene.set(vec, rand.genes)
                       }))
  pval <- length(which(abs(val) < abs(vals))) / num.iters
  ret <- c(val, pval)
  names(ret) <- c("cos.angle", "pval")
  ret
}

plot.loading.components <- function(vec, gene.set) {
  genes <- intersect(names(vec), gene.set)
  ## vec <- as.vector(vec[genes])
  df <- data.frame(val = as.vector(vec), gene = names(vec), in.set = FALSE, x = 1:length(vec))
  flag <- df$gene %in% gene.set
  df$in.set[flag] <- TRUE
  df <- df[order(df$in.set, decreasing=FALSE),]
  g <- ggplot(data = df, aes(x = x, y = val, colour = in.set))
  g <- g + geom_point(alpha = 0.5)
  print(g)
}

all.gene.sets <- c(hallmark.gene.sets, kegg.gene.sets, biocarta.gene.sets, reactome.gene.sets)

projections <- ldply(all.gene.sets, .parallel = TRUE,
                     .fun = function(gene.set) {
                              ## project.real.vector.onto.gene.set(all.pca$rotation[,1], gene.set)
                              calculate.pval.of.projection.of.real.vector.onto.gene.set(all.pca$rotation[,1], gene.set)
                            })

names(projections) <- c("gene.set", "cos.angle", "pval")
projections <- projections[order(abs(projections$cos.angle),decreasing=TRUE),]
projections <- projections[order(abs(projections$pval),decreasing=FALSE),]

flag <- grepl(projections$gene.set, pattern="cycle", ignore.case=TRUE)
projections[flag,]

flag <- grepl(projections$gene.set, pattern="replication", ignore.case=TRUE)
projections[flag,]

flag <- projections$pval < 0.05
projections[flag,]

plot.loading.components(all.pca$rotation[,1], all.gene.sets[["KEGG_DNA_REPLICATION"]])

## Select the extreme outliers and see if they are enriched in any gene sets
pc <- all.pca$rotation[,1]
mean <- mean(pc)
sd <- sd(pc)
num.sds <- 2
flag <- (all.pca$rotation[,1] > (mean + num.sds * sd)) 
pos.outliers <- names(all.pca$rotation[,1])[flag]

flag <- (all.pca$rotation[,1] < (mean - num.sds * sd)) 
neg.outliers <- names(all.pca$rotation[,1])[flag]

plot.loading.components(all.pca$rotation[,1], pos.outliers)
plot.loading.components(all.pca$rotation[,1], neg.outliers)

invert.list <- function(lst) {
  tbl <- ldply(1:length(lst), .parallel = FALSE,
               .fun = function(i) {
                  df <- data.frame(gene = lst[[i]])
                  df$gene.set <- names(lst)[i]
                  df
               })
  tbl
  dlply(tbl, .variables = "gene", .fun = function(df) df$gene.set)
}

inverted.gene.sets <- invert.list(all.gene.sets)
inverted.kegg.gene.sets <- invert.list(kegg.gene.sets)

do.goseq <- function(all.genes, pos.set, gene.sets) {
  genes <- ifelse(all.genes %in% pos.set, 1, 0)
  names(genes) <- all.genes
  pwf <- nullp(genes, "hg19", "ensGene")
  res <- goseq(pwf, gene2cat = gene.sets)
  res
}

kegg.goseq <- do.goseq(names(pc), pos.outliers, inverted.kegg.gene.sets)
kegg.goseq[kegg.goseq$over_represented_pvalue < 0.05,]

cols <- c("category", "over_represented_pvalue", "numDEInCat", "numInCat")
pdf(paste0("kegg-", num.sds, "sds-in-pc1.pdf"))
mytheme <- gridExtra::ttheme_default(
    core = list(fg_params=list(cex = 0.2)))

plot.table(kegg.goseq[kegg.goseq$over_represented_pvalue < 0.05,cols], suppress.rows=TRUE, theme = mytheme)
d <- dev.off()

kegg.goseq <- do.goseq(names(pc), neg.outliers, inverted.kegg.gene.sets)
kegg.goseq[kegg.goseq$over_represented_pvalue < 0.05,]
