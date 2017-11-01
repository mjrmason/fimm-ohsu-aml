
library(glmnet)
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
## suppressPackageStartupMessages(library("GGally"))
suppressPackageStartupMessages(library("VennDiagram"))
suppressPackageStartupMessages(library("caret"))

library(ggbeeswarm)
library(plyr)
library(stringr)
library(synapseClient)

num.cores <- detectCores()
if(!is.na(num.cores) && (num.cores > 1)) {
  suppressPackageStartupMessages(library("doMC"))
  cat(paste("Registering ", num.cores-1, " cores.\n", sep=""))
  registerDoMC(cores=(num.cores-1))
}

synapseLogin()

cat("Translate entrez identifiers in GSEA data sets to ensembl\n")

suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(Biobase))

symbols.to.ensg.mapping <- function(symbols) {
  # Get a mapping from gene symbol to ensembl id
  ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
  bm <- getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id'), 
              filters = 'hgnc_symbol', 
              values = symbols, 
              mart = ensembl)
  names(bm) <- c("symbol", "ensg")
  bm <- bm[!(bm$ensg %in% c("")),]
  bm
}

symbols.to.ensg.mapping <- function(symbols, ensembl, curl = NULL) {
  # Get a mapping from gene symbol to ensembl id
  bm <- getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id'), 
              filters = 'hgnc_symbol', 
              values = symbols, 
              mart = ensembl,
              curl = curl)
  names(bm) <- c("symbol", "ensg")
  bm <- bm[!(bm$ensg %in% c("")),]
  bm
}

# use GSA to read in a gmt file
suppressPackageStartupMessages(library("GSA"))

## Download biocarta, kegg, reactome, and hallmark data sets and store in synapse

translate.gsea.gene.sets <- function(gsea.gene.sets, ensembl, curl) {
  translated.gene.sets <- list()
  for(set.name in gsea.gene.sets$geneset.names) {
      symbol.set <- gsea.gene.sets$genesets[[which(gsea.gene.sets$geneset.names == set.name)]]
      cat(paste0("Translating symbols to ensembl IDs for gene set ", set.name, "\n"))
      ensg.set <- na.omit(symbols.to.ensg.mapping(symbol.set, ensembl, curl))
      if(length(ensg.set) > 0) {
        translated.gene.sets[[set.name]] <- unique(ensg.set$ensg)
      }
  }
  translated.gene.sets
}

ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
curl <- RCurl::getCurlHandle()

## HALLMARK
## http://software.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/6.0/h.all.v6.0.symbols.gmt
synId <- "syn10507487"
obj <- synGet(synId, downloadFile = TRUE)
file <- getFileLocation(obj)
gsea.gene.sets <- GSA.read.gmt(file)
hallmark.gene.sets <- translate.gsea.gene.sets(gsea.gene.sets, ensembl, curl)

## BIOCARTA
## http://software.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/6.0/c2.cp.biocarta.v6.0.symbols.gmt
synId <- "syn10507483"
obj <- synGet(synId, downloadFile = TRUE)
file <- getFileLocation(obj)
gsea.gene.sets <- GSA.read.gmt(file)
biocarta.gene.sets <- translate.gsea.gene.sets(gsea.gene.sets, ensembl, curl)

## KEGG
## http://software.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/6.0/c2.cp.kegg.v6.0.symbols.gmt
synId <- "syn10507485"
obj <- synGet(synId, downloadFile = TRUE)
file <- getFileLocation(obj)
gsea.gene.sets <- GSA.read.gmt(file)
kegg.gene.sets <- translate.gsea.gene.sets(gsea.gene.sets, ensembl, curl)

## REACTOME
## http://software.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/6.0/c2.cp.reactome.v6.0.symbols.gmt
synId <- "syn10507486"
obj <- synGet(synId, downloadFile = TRUE)
file <- getFileLocation(obj)
gsea.gene.sets <- GSA.read.gmt(file)
## vec <- gsea.gene.sets$geneset.names
## names(vec) <- vec
## reactome.gene.sets <-
##  llply(vec,
##        .parallel = TRUE,
##        .fun = function(set.name) {
##                symbol.set <- gsea.gene.sets$genesets[[which(gsea.gene.sets$geneset.names == set.name)]]
##                ensg.set <- na.omit(symbols.to.ensg.mapping(symbol.set))
##                unique(ensg.set$ensg)
##        })

reactome.gene.sets <- translate.gsea.gene.sets(gsea.gene.sets, ensembl, curl)

## Save all of the translated gene sets as R objects and upload them to synapse (in FIMM_BEAT AML Collaboration/Files/Auxiliary Data/Gene Sets)
save(hallmark.gene.sets, file="hallmark.gene.sets.ensg.Rd")
save(biocarta.gene.sets, file="biocarta.gene.sets.ensg.Rd")
save(kegg.gene.sets, file="kegg.gene.sets.ensg.Rd")
save(reactome.gene.sets, file="reactome.gene.sets.ensg.Rd")

parentId <- "syn10507476"

for(file.name in c("hallmark.gene.sets.ensg.Rd", "biocarta.gene.sets.ensg.Rd", "kegg.gene.sets.ensg.Rd", "reactome.gene.sets.ensg.Rd")) {
  f <- File(file.name, parentId = parentId, synapseStore = TRUE)
  synStore(f, forceVersion = FALSE)
}


