
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

source("lee-setup-050318.R")

data.set <- "lee"

source("../common/data-preprocessing.R")
source("../common/utils.R")
source("../common/plotting.R")

## Read in the drug response AUCs
synId <- data.set.drug.fit.synIds[[data.set]]
file <- getFileLocation(synGet(synId))
fits <- read.table(file, sep = "\t", header=TRUE, as.is = TRUE)

exclude.outliers <- TRUE
if(exclude.outliers) {
  exclude.cols <- colnames(fits)[grepl(colnames(fits), pattern="exclude") & grepl(colnames(fits), pattern="ll4")]
  exclude.cols <- exclude.cols[!grepl(exclude.cols, pattern="gof")]
  exclude.cols <- exclude.cols[!grepl(exclude.cols, pattern="any")]
  flag <- unlist(apply(fits[, exclude.cols], 1, function(row) { r <- row[!is.na(row)]; return(any(r)) }))
  fits <- fits[!flag, ]
}

drc.refit <- prepare.drug.response.matrix(fits, drug.col = data.set.drug.id.cols[[data.set]], patient.col = expr.patient.id.cols[[data.set]], response.col = "auc.ll4", drugs = NULL)

## Read in the gene expression
synId <- data.set.expr.synIds[[data.set]]
file <- getFileLocation(synGet(synId))
expr <- read.table(file, sep = "\t", header = TRUE, as.is = TRUE)

venetoclax.feature.syms <- c("BCL2", "CD14", "FLT3", "MAFB", "LRP1")
venetoclax.feature.gene.tbl <- symbols.to.ensg.mapping(venetoclax.feature.syms)
venetoclax.feature.gene.tbl$ensg <- as.character(venetoclax.feature.gene.tbl$ensg)

drugs.of.interest <- c("Venetoclax", "TG100-115", "XAV-939", "Lovastatin", "Panobinostat", "Neratinib", "PD184352", "Trametinib", "Palbociclib", "Pelitinib")

drug.feature.syms <- list()
drug.feature.syms[["XAV-939"]] <- NA
## MEK inhibitors
drug.feature.syms[["Trametinib"]] <- c("CD300E", "TBXAS1", "CEBPA", "TAL1", "RAMP1")
drug.feature.syms[["Selumetinib"]] <- c("TNNT1", "CD14", "CD300E", "CD163", "LGALS2")
drug.feature.syms[["PD184352"]] <- c("LRP1", "CEBPA", "VCAN", "CD300E", "ARHGEF10L")
drug.feature.syms[["Venetoclax"]] <- venetoclax.feature.syms

drug.feature.syms[["SNS-032"]] <- c("VCAN", "MEIS1", "CTSG", "ROBO1", "SLITRK5")
drug.feature.syms[["BI 2536"]] <- c("LINC01475", "NKX2-3", "MYCN", "PF4", "MPO")
drug.feature.syms[["Panobinostat"]] <- c("LRP1", "ZNF711", "CPA3", "CD163", "MAFB")
drug.feature.syms[["Palbociclib"]] <- c("CD14", "LRP1", "MAML2", "SKIDA1", "TNFRSF14")

drug.feature.syms[["Selinexor"]] <- NA
drug.feature.syms[["Lovastatin"]] <- NA
drug.feature.syms[["AGI-5198"]] <- NA
drug.feature.syms[["MEK"]] <- unique(c(drug.feature.syms[["Selumetinib"]], drug.feature.syms[["Trametinib"]], drug.feature.syms[["PD184352"]]))

## Load in the drug map between Lee et al, FIMM, and OHSU

synId <- drug.mapping.synId
file <- getFileLocation(synGet(synId))
drug.map <- read.table(file, sep = "\t", header=TRUE, as.is = TRUE, comment.char = "")

overlapping.drugs <- drug.map[drug.map[, "FIMM_DRUG_NAME"] %in% names(drug.feature.syms), ]
gene.to.feature.map <- data.frame(drug = overlapping.drugs[, drug.map.drug.id.cols[[ds]]], features = overlapping.drugs$FIMM_DRUG_NAME, stringsAsFactors = FALSE)
gene.to.feature.map <- rbind(gene.to.feature.map, c("ABT-263", "Venetoclax"))
gene.to.feature.map <- rbind(gene.to.feature.map, c("AZD-6244", "MEK"))

                             
processed.exprs <- list()
processed.exprs[[data.set]] <- expr
processed.drcs.refit <- list()
processed.drcs.refit[[data.set]] <- drc.refit

obj <- synGet("syn12177231", downloadFile=TRUE)
drc <- as.data.frame(read.table(getFileLocation(obj), sep="\t", header=TRUE))

obj <- synGet("syn12177264", downloadFile=TRUE)
drc.ic50 <- as.data.frame(read.table(getFileLocation(obj), sep="\t", header=TRUE))

processed.drcs.published <- list()
## Invert the sign
processed.drcs.published[[data.set]] <- - drc
## processed.drcs <- list("lee" = drc.ic50)


output.correlations <- function(ds, processed.exprs, processed.drcs, prefix) {
  samples <- intersect(colnames(processed.exprs[[ds]]), colnames(processed.drcs[[ds]]))
  l_ply(1:nrow(gene.to.feature.map),
        .fun = function(i) {
                 drug <- gene.to.feature.map$drug[i]
                 feature.name <- gene.to.feature.map$features[i]
                 features <- drug.feature.syms[[feature.name]]
                 features <- na.omit(features)
                 features <- features[features %in% rownames(processed.exprs[[ds]])]
                 if(length(features) == 0) { return() }
  print(features)
                 lst <- llply(1:length(features),
                              .fun = function(indx) {
                                       sym <- features[indx]
                                       flag <- !is.na(processed.exprs[[ds]][sym, samples]) & !is.na(processed.drcs[[ds]][drug, samples])
                                       if(length(which(flag)) == 0) { return(NULL) }
                                       gene.expr <- scale(as.numeric(processed.exprs[[ds]][sym,samples]))
                                       drug.resp <- scale(as.numeric(processed.drcs[[ds]][drug,samples]))
                                       g <- plot.correlation(gene.expr, drug.resp, display.pval = TRUE, size = 6, colour = "blue")
                                       g <- g + ylab(paste("Drug Response", collapse=" "))
                                       g <- g + xlab(paste(sym, "Expression", collapse=" "))
                                       g <- g + theme(text = element_text(size = 15))
                                       g
                              })
                 pg <- plot_grid(plotlist = lst, labels = LETTERS[1:length(lst)])
                 title <- ggdraw() + draw_label(paste(prefix, drug, "vs", feature.name, "biomarkers", sep = " "))
                 pg <- plot_grid(title, pg, ncol=1, rel_heights=c(0.1, 1)) # rel_heights values control title margin
                 pg
                 file <- make.names(paste0(prefix, "-", drug, "-vs-", feature.name, "-biomarkers.pdf"))
                 ggsave(file, width = 14)
         })
}      

output.correlations(ds, processed.exprs, processed.drcs.refit, prefix = paste0(ds, " (refitted AUCs)"))
output.correlations(ds, processed.exprs, processed.drcs.published, prefix = paste0(ds, " (published AUCs)"))

gene.tbl <- subset(venetoclax.feature.gene.tbl, symbol %in% rownames(processed.exprs[[ds]]))



first.batch.samples <- c("AML017", "AML030", "AML044", "AML056", "AML059", "AML060", "AML061", "AML067", "AML068", "AML069", "AML070", "AML076")

drug <- "ABT-263"
lst <- llply(1:nrow(gene.tbl), 
             .fun = function(i) {
                      ensg <- gene.tbl$ensg[i]
                      sym <- gene.tbl$symbol[i]
                      samples <- intersect(colnames(processed.exprs[[ds]]), colnames(processed.drcs[[ds]]))
##                      samples <- samples[(samples %in% first.batch.samples)]
                      gene.expr <- scale(as.numeric(processed.exprs[[ds]][sym,samples]))
                      drug.resp <- scale(as.numeric(processed.drcs[[ds]][drug,samples]))
                      g <- plot.correlation(gene.expr, drug.resp, display.pval = TRUE, size = 6, colour = "blue")
##                      g <- g + ylab(paste(drug, "Actual", "Drug Response\n(Standardized)", collapse=" "))
##                      g <- g + ylab(paste("Drug Response\n(Standardized)", collapse=" "))
##                      g <- g + xlab(paste(sym, "Expression", "\n(Standardized)", collapse=" "))
                      g <- g + ylab(paste("Drug Response", collapse=" "))
                      g <- g + xlab(paste(sym, "Expression", collapse=" "))
                      g <- g + theme(text = element_text(size = 15))
                      g
                    })
suppressPackageStartupMessages(library("cowplot"))
pg <- plot_grid(lst[[1]], lst[[2]], lst[[3]], lst[[4]], lst[[5]], labels = c("A", "B", "C", "D", "E"))

pg
ggsave("venetoclax-lee-correlations.pdf", width = 14)
## d <- dev.off()
cat("Done in plot-venetoclax-model-and-gene-correlations.R")
