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
## source("heatmap-utils.R")

## Register parallelization
num.cores <- detectCores()
if(!is.na(num.cores) && (num.cores > 1)) {
  suppressPackageStartupMessages(library("doMC"))
  cat(paste("Registering ", num.cores-1, " cores.\n", sep=""))
  registerDoMC(cores=(num.cores-1))
}
num.processes <- num.cores - 1

## Begin setup (this will need to be changed as the data sets in the intersection change)

source("fimm-ohsu-setup-040218.R")

source("process-drc-and-expr.R")

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

file <- "Heikki_FAB_media_annot_Feb_2018.xlsx"
fimm.fab <- openxlsx::read.xlsx(file, sheet = 1)
fimm.fab <- unique(fimm.fab[, c("IDs", "medium", "FAB")])
fimm.fab <- fimm.fab[grepl(fimm.fab$medium, pattern="MCM"),]
rownames(fimm.fab) <- fimm.fab$IDs

samples <- intersect(colnames(processed.drcs[["fimm"]]), fimm.fab$IDs)

ff <- subset(fimm.fab, IDs %in% colnames(processed.drcs[["fimm"]]))
ff[order(ff$IDs),]

ven <- processed.drcs[["fimm"]]["Venetoclax",samples]
fab <- fimm.fab[samples, "FAB"]
df <- data.frame(ven = as.numeric(ven), fab = fab)

g <- ggplot(data = df, aes(x = fab, y = ven))
g <- g + geom_boxplot()
g <- g + geom_point()
g <- g + xlab("FAB")
g <- g + ylab("Venetoclax Response")
pdf("venetoclax-response-vs-fab.pdf")
print(g)
d <- dev.off()

lmo <- lm(ven ~ 0 + fab, data = df)
lmo <- lm(ven ~ fab, data = df)