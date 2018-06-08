rm(list = ls())

rdata.file <- ".Rdata.model.mean.response.040218"
load(rdata.file)

suppressPackageStartupMessages(library(glmnet))
suppressPackageStartupMessages(library("randomForest"))
suppressPackageStartupMessages(library("e1071"))
suppressPackageStartupMessages(library("foreach"))
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library(openxlsx))
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

source("fimm-ohsu-setup-040218.R")

## Turn _off_ gene subsetting
use.subset <- FALSE

source("process-drc-and-expr.R")

## A string that will be included in any output files.
file.prefix <- "model-glds"

for(train.set in names(res.aucs[["no.transform"]][["all.comparisons"]])) {
  for(test.set in names(res.aucs[["no.transform"]][["all.comparisons"]][[train.set]])) {
    for(gene.set in names(res.aucs[["no.transform"]][["all.comparisons"]][[train.set]][[test.set]])) {
      models <- c("glmnet", "glmnet", "rf")
      svec <- c("lambda.min", "lambda.min", NA)
      alphavec <- c(1, 0, NA)
      nms <- c("lasso", "ridge", "rf")
      for(i in 1:length(models)) {
        non.rand <- extract.predicted.actual(res.aucs[["no.transform"]][["all.comparisons"]][[train.set]][[test.set]][[gene.set]],
                                             train.drug = "mean.resp", alpha = alphavec[i], model = models[i], s = svec[i])
 
        g1 <- plot.correlation(non.rand$actual, non.rand$predicted, display.pval = TRUE, size = 5, colour = "blue")
        g1 <- g1 + xlab(paste("Observed GLDS (FIMM)"))
        g1 <- g1 + ylab(paste("Predicted GLDS (FIMM)"))
        g1 <- g1 + theme(text = element_text(size = 20))
##        g1 <- plot.r2(non.rand$predicted, non.rand$actual)
##        g1 <- g1 + ggtitle(paste0("Mean response ", train.set, " vs ", test.set, " (", nms[i], ": ", gene.set, ")"))
##        g1 <- g1 + ylab("Actual Mean Response")
##        g1 <- g1 + xlab("Predicted Mean Response")

        if(FALSE) {
        rand <- extract.predicted.actual(res.aucs[["random"]][["all.comparisons"]][[train.set]][[test.set]][[gene.set]],
                                             train.drug = "mean.resp", alpha = alphavec[i], model = models[i], s = svec[i])
        g2 <- plot.r2(rand$predicted, rand$actual)
        g2 <- g2 + ggtitle(paste0("Random mean response ", train.set, " vs ", test.set, " (", nms[i], ": ", gene.set, ")"))
        g2 <- g2 + ylab("Actual Mean Response")
        g2 <- g2 + xlab("Predicted Mean Response")
        }
        file <- paste(file.prefix, train.set, "vs", test.set, nms[i], gene.set, "fit.pdf", sep="-")
        cat(paste0("Writing ", file, "\n"))
        pdf(file)
        ## grid.arrange(g1, g2)
        print(g1)
        d <- dev.off()

        names(non.rand$predicted) <- non.rand$sample.names
        file <- paste(file.prefix, train.set, "vs", test.set, nms[i], gene.set, "prediction.tsv", sep="-")
        write.table(file = file, t(non.rand$predicted), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
      }
    }
  }
}
