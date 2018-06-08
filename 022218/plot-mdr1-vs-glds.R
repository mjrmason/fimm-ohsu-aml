rm(list = ls())

suppressPackageStartupMessages(library("foreach"))
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library(tidyr))

suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("gridExtra"))
suppressPackageStartupMessages(library("caTools"))
suppressPackageStartupMessages(library("data.table"))

suppressPackageStartupMessages(library(ggbeeswarm))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(synapseClient))
suppressPackageStartupMessages(library(mygene))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library("cowplot"))

library(rJava)
library(rcdk)
library(fingerprint)
library(plyr)
library(tidyverse)

library(dynamicTreeCut)
library(MCL)
library(clusteval)

library(WeightedCluster)
library(igraph)

library("dendextend")
## library("dendextendRcpp")

library(corrplot)

## Log in to Synapse
synapseLogin()

## Bring in common code
source("../common/data-preprocessing.R")
source("../common/models.R")
source("../common/plotting.R")
source("../common/utils.R")
source("heatmap-utils.R")

## These are defnied in heatmap-utils
orig.ohsu.metadata <- ohsu.metadata
orig.fimm.metadata <- fimm.metadata

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

rdata.file <- ".Rdata.model.mean.response.120817"

## Over-write some of the variables
## These fit files are assumed to be for AUCs/DSS values calculated across shared concentration ranges across data sets.
## These files may be created by a script such as
## calculate-dss-fimm-ohsu.R
## data.set.dss.fit.synIds <- c("syn11362885", "syn11362891")
## names(data.set.dss.fit.synIds) <- data.sets
## data.set.dss.fit.files <- c("ohsu-fimm-ohsu-outlier1-dss-t0.tsv", "fimm-fimm-ohsu-outlier1-dss-t0.tsv")

## A string that will be included in any output files.
file.prefix <- "fimm-ohsu-validate-trained-models-112217-one-ohsu-no-filtering-all-mean-response"

## End setup

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

orig.glds <- llply(data.sets, .parallel = FALSE,
              .fun = function(ds) {
                       m <- processed.drcs[[ds]]
                       m <- t(scale(t(m), scale = TRUE, center = TRUE))
                       colMeans(m, na.rm=TRUE)
                     })

glist <- llply(data.sets,
               .fun = function(ds) {
                        common.cols <- intersect(names(orig.glds[[ds]]), colnames(processed.exprs[[ds]]))
                        ## Plot correation (not r2) and include p-value
                        g <- plot.correlation(as.numeric(orig.glds[[ds]][common.cols]), 
                                              as.numeric(processed.exprs[[ds]]["ENSG00000085563",common.cols]), 
                                              display.pval = TRUE, method = "pearson", display.r2 = FALSE, size = 5, colour = "blue")
                        g <- g + ggtitle(toupper(ds))
                        g <- g + xlab("MDR1 Expression")
                        g <- g + ylab("GLDS")
                        g <- g + theme(text = element_text(size = 20))
                      })
pg <- plot_grid(glist[[1]], glist[[2]], labels = c("A", "B"), nrow=1, rel_widths = c(10,10), align="none")
pg
ggsave("mdr1-vs-glds.png", width = 14)

l_ply(data.sets,
      .fun = function(ds) {
               common.cols <- intersect(names(orig.glds[[ds]]), colnames(processed.exprs[[ds]]))
               sum <- summary(lm(as.numeric(processed.exprs[[ds]]["ENSG00000085563",common.cols]) ~ as.numeric(orig.glds[[ds]][common.cols])))
               cat(paste0(ds, ": ", "n = ", length(common.cols), "\n", sep=""))
               cat(paste0(ds, ": percent variance in GLDS explained by MDR1 = ", sum$r.squared * 100, "\n", sep=""))
               print(sum)
             })

## Calculate total explained in matrix by GLDS
l_ply(data.sets,
      .fun = function(ds) {
               glds <- orig.glds[[ds]]
               m <- processed.drcs[[ds]]
               m <- t(scale(t(m), scale = TRUE, center = TRUE))
               tots <- llply(1:nrow(m), .parallel = FALSE, .fun = function(i) { v <- as.numeric(m[i,]); sum((v - mean(v, na.rm=TRUE))^2, na.rm=TRUE) })
               rsses <- llply(1:nrow(m), .parallel = FALSE, .fun = function(i) { v <- as.numeric(m[i,]); sum((v - glds)^2, na.rm=TRUE) })
               cors <- unlist(llply(1:nrow(m), .parallel = FALSE, .fun = function(i) { v <- as.numeric(m[i,]); cor.test(v, glds, method = "pearson")$estimate }))
               frac.explained <- Reduce('+', rsses) / Reduce('+', tots)
               fracs.explained <- unlist(rsses)/unlist(tots)
##               cat(paste0(ds, ": percent variance explained by GDS = ", format(frac.explained * 100, digits = 3), "\n", sep=""))
               qs <- quantile(fracs.explained, probs = c(0.25, 0.5, 0.75), names = TRUE)
               cat(paste0(ds, ": percent variance explained by GDS across drugs: ", 
                   paste(unlist(lapply(1:length(qs), function(i) paste0(names(qs)[i], " = ", format(100*qs[i], digits=3)))), collapse=" "), "\n"))
               qs <- quantile(cors, probs = c(0.25, 0.5, 0.75), names = TRUE)
               cat(paste0(ds, ": correlation of drugs with glds: ", 
                   paste(unlist(lapply(1:length(qs), function(i) paste0(names(qs)[i], " = ", format(qs[i], digits=3)))), collapse=" "), "\n"))
             })

## ohsu: percent variance explained by GDS across drugs: 25% = 55.8 50% = 66.8 75% = 87.1
## ohsu: correlation of drugs with glds: 25% = 0.398 50% = 0.589 75% = 0.686
## fimm: percent variance explained by GDS across drugs: 25% = 53.5 50% = 66.4 75% = 85
## fimm: correlation of drugs with glds: 25% = 0.421 50% = 0.584 75% = 0.724


## Do univariate analysis
glds.uni <- 
llply(data.sets,
      .fun = function(ds) {
               common.cols <- intersect(names(orig.glds[[ds]]), colnames(processed.exprs[[ds]]))
               glds <- as.numeric(orig.glds[[ds]][common.cols])
               indices <- 1:nrow(processed.exprs[[ds]])
               names(indices) <- rownames(processed.exprs[[ds]])
               df <- ldply(indices, .parallel = FALSE,
                           .fun = function(i) {
                                    v <- as.numeric(processed.exprs[[ds]][i,common.cols])
                                    ct <- cor.test(glds, v, method = "pearson")
                                    data.frame(est = ct$estimate, pval = ct$p.value)
                           })
               names(df)[1] <- "gene"
               df <- df[order(df$pval),]
               up.df <- subset(df, est > 0)
               up.df$padj <- p.adjust(up.df$pval)
               down.df <- subset(df, est < 0)
               down.df$padj <- p.adjust(down.df$pval)
               write.table(file = paste0(ds, "-up-glds-univariate.tsv"), up.df, sep="\t", row.names=FALSE, col.names=TRUE, quote = FALSE)
               write.table(file = paste0(ds, "-down-glds-univariate.tsv"), down.df, sep="\t", row.names=FALSE, col.names=TRUE, quote = FALSE)
               df
             })

m <- merge(glds.uni[["fimm"]], glds.uni[["ohsu"]], by = "gene", suffixes = c(".fimm", ".ohsu"))
m$up.fimm <- m$est.fimm > 0
m$up.ohsu <- m$est.ohsu > 0

library(metap)

ohsu.up.tbl <- subset(m, est.ohsu > 0)
ohsu.up.tbl$fisher.pval <- unlist(apply(ohsu.up.tbl[, c("pval.fimm", "pval.ohsu")], 1, function(row) sumlog(row)$p))
ohsu.up.tbl <- ohsu.up.tbl[order(ohsu.up.tbl$fisher.pval),]
ohsu.up.tbl <- subset(ohsu.up.tbl, (up.fimm == TRUE) & (up.ohsu == TRUE)) 
write.table(file = "fimm-and-ohsu-up-glds-fisher.tsv", ohsu.up.tbl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

ohsu.down.tbl <- subset(m, est.ohsu < 0)
ohsu.down.tbl$fisher.pval <- unlist(apply(ohsu.down.tbl[, c("pval.fimm", "pval.ohsu")], 1, function(row) sumlog(row)$p))
ohsu.down.tbl <- ohsu.down.tbl[order(ohsu.down.tbl$fisher.pval),]
ohsu.down.tbl <- subset(ohsu.down.tbl, (up.fimm == FALSE) & (up.ohsu == FALSE)) 
write.table(file = "fimm-and-ohsu-down-glds-fisher.tsv", ohsu.down.tbl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

