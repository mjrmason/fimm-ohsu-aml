
initialize <- TRUE

if(initialize) {
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

} # if(initialize)


obj <- synGet("syn12177231", downloadFile=TRUE)
drc <- as.data.frame(read.table(getFileLocation(obj), sep="\t", header=TRUE))

obj <- synGet("syn12177264", downloadFile=TRUE)
drc.ic50 <- as.data.frame(read.table(getFileLocation(obj), sep="\t", header=TRUE))

obj <- synGet("syn12176399", downloadFile=TRUE)
expr <- as.data.frame(read.table(getFileLocation(obj), sep="\t", header=TRUE))

processed.drcs <- list("lee" = drc)
## processed.drcs <- list("lee" = drc.ic50)

processed.exprs <- list("lee" = expr)

venetoclax.feature.syms <- c("BCL2", "CD14", "FLT3", "MAFB", "LRP1")
venetoclax.feature.gene.tbl <- symbols.to.ensg.mapping(venetoclax.feature.syms)
venetoclax.feature.gene.tbl$ensg <- as.character(venetoclax.feature.gene.tbl$ensg)

drugs.of.interest <- c("Venetoclax", "TG100-115", "XAV-939", "Lovastatin", "Panobinostat", "Neratinib", "PD184352", "Trametinib", "Palbociclib", "Pelitinib")

drug.feature.syms <- list()
drug.feature.syms[["XAV-939"]] <- c()
drug.feature.syms[["Trametinib"]] <- c()
drug.feature.syms[["Venetoclax"]] <- c()
drug.feature.syms[["Selumetinib"]] <- c()
drug.feature.syms[["Selinexor"]] <- c()
drug.feature.syms[["SNS-032"]] <- c()
drug.feature.syms[["Panobinostat"]] <- c()
drug.feature.syms[["Palbociclib"]] <- c()
drug.feature.syms[["PD184352"]] <- c()
drug.feature.syms[["Lovastatin"]] <- c()
drug.feature.syms[["BI 2536"]] <- c()
drug.feature.syms[["AGI-5198"]] <- c()





gene.tbl <- subset(venetoclax.feature.gene.tbl, symbol %in% rownames(processed.exprs[["lee"]]))

first.batch.samples <- c("AML017", "AML030", "AML044", "AML056", "AML059", "AML060", "AML061", "AML067", "AML068", "AML069", "AML070", "AML076")

drug <- "ABT-263"
lst <- llply(1:nrow(gene.tbl), 
             .fun = function(i) {
                      ensg <- gene.tbl$ensg[i]
                      sym <- gene.tbl$symbol[i]
                      ds <- "lee"
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
