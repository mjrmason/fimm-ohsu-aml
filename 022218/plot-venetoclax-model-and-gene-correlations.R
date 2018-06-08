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
suppressPackageStartupMessages(library(mygene))
suppressPackageStartupMessages(library(GGally))

suppressPackageStartupMessages(library("GSVA"))
# use GSA to read in a gmt file
suppressPackageStartupMessages(library("GSA"))
suppressPackageStartupMessages(library("cowplot"))


## Log in to Synapse
synapseLogin()

## Bring in common code
source("../common/data-preprocessing.R")
source("../common/models.R")
source("../common/plotting.R")
source("../common/utils.R")
source("heatmap-utils.R")

## Register parallelization
num.cores <- detectCores()
if(!is.na(num.cores) && (num.cores > 1)) {
  suppressPackageStartupMessages(library("doMC"))
  cat(paste("Registering ", num.cores-1, " cores.\n", sep=""))
  registerDoMC(cores=(num.cores-1))
}
num.processes <- num.cores - 1

## Begin setup (this will need to be changed as the data sets in the intersection change)

source("fimm-ohsu-setup-120817.R")

source("process-drc-and-expr.R")

rdata.file <- ".Rdata.overlap.120817"
load(rdata.file)

response.col <- "auc.ll4"

processed.drcs <- list()
processed.exprs <- list()

ds <- "ohsu"
l <- prepare.drug.response.and.expr.matrices(orig.fits[[ds]], orig.exprs[[ds]], drugs = orig.drug.name.tbl[, drug.name.col],
                                             drug.col = drug.name.col, patient.col = expr.patient.id.cols[[ds]],
                                             response.col = response.col)

processed.drcs[[ds]] <- l[["drc.df"]]
processed.exprs[[ds]] <- l[["expr.df"]]

ds <- "fimm"
l <- prepare.drug.response.and.expr.matrices(orig.fits[[ds]], orig.exprs[[ds]], drugs = orig.drug.name.tbl[, drug.name.col],
                                             drug.col = drug.name.col, patient.col = expr.patient.id.cols[[ds]],
                                             response.col = response.col)

processed.drcs[[ds]] <- l[["drc.df"]]
processed.exprs[[ds]] <- l[["expr.df"]]

drug <- "Venetoclax"
drug.pred.actual <- extract.predicted.actual(res.aucs[["unpenalized.mean.feature"]][["all.comparisons"]][["ohsu"]][["fimm"]][["gene"]], train.drug = drug, alpha = 0, s = NA, model = "glmnet")
## g <- plot.correlation(drug.pred.actual$actual, drug.pred.actual$predicted, display.pval = TRUE, size = 6)
samples <- intersect(drug.pred.actual$sample.names, colnames(processed.drcs[[ds]]))
drug.resp <- scale(as.numeric(processed.drcs[[ds]][drug,samples]))
g <- plot.correlation(drug.resp, drug.pred.actual$predicted, display.pval = FALSE, size = 6)
##g <- g + xlab(paste(drug, "Actual", "Drug Response\n(Standardized)", collapse=" "))
##g <- g + ylab(paste(drug, "Predicted Drug", "Response\n(Standardized)", collapse=" "))
##g <- g + xlab(paste("Drug Response\n(Standardized)", collapse=" "))
##g <- g + ylab(paste("Predicted Drug", "\nResponse (Standardized)", collapse=" "))
g <- g + xlab(paste("Drug Response", collapse=" "))
g <- g + ylab(paste("Predicted Drug", "\nResponse", collapse=" "))
g <- g + theme(text = element_text(size = 15))

venetoclax.feature.syms <- c("BCL2", "CD14", "FLT3", "MAFB", "LRP1")
venetoclax.feature.gene.tbl <- symbols.to.ensg.mapping(venetoclax.feature.syms)
venetoclax.feature.gene.tbl$ensg <- as.character(venetoclax.feature.gene.tbl$ensg)

gene.tbl <- subset(venetoclax.feature.gene.tbl, ensg %in% rownames(processed.exprs[["fimm"]]))

lst <- llply(1:nrow(gene.tbl), 
             .fun = function(i) {
                      ensg <- gene.tbl$ensg[i]
                      sym <- gene.tbl$symbol[i]
                      ds <- "fimm"
                      samples <- intersect(colnames(processed.exprs[[ds]]), colnames(processed.drcs[[ds]]))
                      gene.expr <- scale(as.numeric(processed.exprs[[ds]][ensg,samples]))
                      drug.resp <- scale(as.numeric(processed.drcs[[ds]][drug,samples]))
                      g <- plot.correlation(gene.expr, drug.resp, display.pval = FALSE, size = 6)
##                      g <- g + ylab(paste(drug, "Actual", "Drug Response\n(Standardized)", collapse=" "))
##                      g <- g + ylab(paste("Drug Response\n(Standardized)", collapse=" "))
##                      g <- g + xlab(paste(sym, "Expression", "\n(Standardized)", collapse=" "))
                      g <- g + ylab(paste("Drug Response", collapse=" "))
                      g <- g + xlab(paste(sym, "Expression", collapse=" "))
                      g <- g + theme(text = element_text(size = 15))
                      g
                    })
pg <- plot_grid(g, lst[[1]], lst[[2]], lst[[3]], lst[[4]], lst[[5]], labels = c("A", "B", "C", "D", "E", "F"))

## pdf("venetoclax-fimm-correlations.pdf")
pg
ggsave("venetoclax-fimm-correlations.pdf", width = 14)
## d <- dev.off()
stop("stop")
