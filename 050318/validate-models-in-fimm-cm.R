
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library(ggbeeswarm))
suppressPackageStartupMessages(library(synapseClient))
suppressPackageStartupMessages(library(mygene))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(cowplot))

## Log in to Synapse
synapseLogin()

source("fimm-cm-setup-050318.R")

data.set <- "fimm-cm"

source("../common/data-preprocessing.R")
source("../common/utils.R")
source("../common/plotting.R")

## Read in the drug response AUCs
synId <- processed.drug.fit.synIds[[data.set]]
file <- getFileLocation(synGet(synId))
fits <- read.table(file, sep = "\t", header=TRUE, as.is = TRUE)

exclude.outliers <- TRUE
if(exclude.outliers) {
  exclude.cols <- colnames(fits)[grepl(colnames(fits), pattern="exclude") & grepl(colnames(fits), pattern="ll4")]
##  exclude.cols <- exclude.cols[!grepl(exclude.cols, pattern="gof")]
##  exclude.cols <- exclude.cols[!grepl(exclude.cols, pattern="any")]
  flag <- unlist(apply(fits[, exclude.cols], 1, function(row) { r <- row[!is.na(row)]; return(any(r)) }))
  fits <- fits[!flag, ]
}

## Read in the gene expression
synId <- data.set.orig.expr.synIds[[data.set]]
file <- getFileLocation(synGet(synId))
expr <- read.table(file, sep = "\t", header = TRUE, as.is = TRUE)

orig.expr <- expr


ret <- collapse.drug.response.and.expr.matrices.samples.to.patient(fits, expr, drugs = NULL, drug.col = data.set.drug.id.cols[[data.set]], 
                                                                   sample.col = sample.id.cols[[data.set]], 
                                                                   patient.col = patient.id.cols[[data.set]], drug.response.col = drug.response.cols[[data.set]])

venetoclax.feature.syms <- c("BCL2", "CD14", "FLT3", "MAFB", "LRP1")
venetoclax.feature.gene.tbl <- symbols.to.ensg.mapping(venetoclax.feature.syms)
venetoclax.feature.gene.tbl$ensg <- as.character(venetoclax.feature.gene.tbl$ensg)

drugs.of.interest <- c("Venetoclax", "TG100-115", "XAV-939", "Lovastatin", "Panobinostat", "Neratinib", "PD184352", "Trametinib", "Palbociclib", "Pelitinib", "Navitoclax")

drug.feature.syms <- list()
## XAV-939 has a poor actual/predicted correlation
## drug.feature.syms[["XAV-939"]] <- NA
## MEK inhibitors
drug.feature.syms[["Trametinib"]] <- c("CD300E", "TBXAS1", "CEBPA", "TAL1", "RAMP1")
drug.feature.syms[["Selumetinib"]] <- c("TNNT1", "CD14", "CD300E", "CD163", "LGALS2")
drug.feature.syms[["PD184352"]] <- c("LRP1", "CEBPA", "VCAN", "CD300E", "ARHGEF10L")
drug.feature.syms[["Venetoclax"]] <- venetoclax.feature.syms
drug.feature.syms[["Navitoclax"]] <- venetoclax.feature.syms

drug.feature.syms[["SNS-032"]] <- c("VCAN", "MEIS1", "CTSG", "ROBO1", "SLITRK5")
drug.feature.syms[["BI 2536"]] <- c("LINC01475", "NKX2-3", "MYCN", "PF4", "MPO")
drug.feature.syms[["Panobinostat"]] <- c("LRP1", "ZNF711", "CPA3", "CD163", "MAFB")
drug.feature.syms[["Palbociclib"]] <- c("CD14", "LRP1", "MAML2", "SKIDA1", "TNFRSF14")

## drug.feature.syms[["Selinexor"]] <- NA
## drug.feature.syms[["Lovastatin"]] <- NA
## drug.feature.syms[["AGI-5198"]] <- NA
## drug.feature.syms[["MEK"]] <- unique(c(drug.feature.syms[["Selumetinib"]], drug.feature.syms[["Trametinib"]], drug.feature.syms[["PD184352"]]))

drug.feature.name.to.id.map <- llply(drug.feature.syms,
                                     .fun = function(gene.syms) {
                                              df <- symbols.to.ensg.mapping(gene.syms)
                                              df$symbol <- as.character(df$symbol)
                                              df$ensg <- as.character(df$ensg)
                                              df <- df[, c("symbol", "ensg")]
                                              colnames(df) <- c("name", "id")
                                              df <- na.omit(df)
                                              df
                                     })

## Load in the drug map
synId <- drug.mapping.synId
file <- getFileLocation(synGet(synId))
drug.map <- read.table(file, sep = "\t", header=TRUE, as.is = TRUE, comment.char = "")

overlapping.drugs <- drug.map[drug.map[, "FIMM_DRUG_NAME"] %in% names(drug.feature.syms), ]
drug.to.feature.map <- data.frame(drug.name = overlapping.drugs[, "FIMM_DRUG_NAME"], 
                                  drug = overlapping.drugs[, drug.map.drug.id.cols[[data.set]]], features = overlapping.drugs$FIMM_DRUG_NAME, stringsAsFactors = FALSE)
dm <- unique(fits[, c("DRUG_NAME", data.set.drug.id.cols[[data.set]])])
flag <- dm$DRUG_NAME == "Navitoclax"
drug.to.feature.map <- rbind(drug.to.feature.map, c("Navitoclax", dm[flag, data.set.drug.id.cols[[data.set]]], "Venetoclax"))
flag <- dm$DRUG_NAME == "Selumetinib"
drug.to.feature.map <- rbind(drug.to.feature.map, c("Selumetinib", dm[flag, data.set.drug.id.cols[[data.set]]], "Venetoclax"))

processed.exprs <- list()
processed.exprs[[data.set]] <- ret[["expr.df"]]
processed.drcs <- list()
processed.drcs[[data.set]] <- ret[["drc.df"]]


my.dup <- function(x) duplicated(x, fromLast = TRUE) | duplicated(x, fromLast = FALSE)
drug.ranges <- merge(unique(fits[, c("DRUG_NAME", min.conc.cols[[data.set]], max.conc.cols[[data.set]])]), overlapping.drugs, by.x = "DRUG_NAME", by.y = "FIMM_DRUG_NAME")
drug.ranges <- drug.ranges[, c("DRUG_NAME", min.conc.cols[[data.set]], max.conc.cols[[data.set]])]
if(any(my.dup(drug.ranges$DRUG_NAME))) {
  stop("Was not expecting multiple drug ranges per drug\n")
}
colnames(drug.ranges) <- c("drug.name", "min.conc.nM", "max.conc.nM")

output.dir <- "output"
dir.create(output.dir)
output.correlations(data.set, processed.exprs, processed.drcs, drug.to.feature.map, drug.feature.name.to.id.map = drug.feature.name.to.id.map, 
                    drug.ranges = drug.ranges, prefix = paste0(output.dir, "/", data.set, "-validation"), title = paste0(toupper(data.set), ":"))



