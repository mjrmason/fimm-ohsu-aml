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

suppressPackageStartupMessages(library("GSVA"))
# use GSA to read in a gmt file
suppressPackageStartupMessages(library("GSA"))

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

rdata.file <- ".Rdata.model.mean.response.040218"

## A string that will be included in any output files.
file.prefix <- "model-glds"

## End setup

cat(paste0("Model AUC ~ expr for data sets: ", paste(data.sets, collapse=","), "\n"))

expr.filt.matrices <- exprs
if(FALSE) {
  expr.filt.matrices <- list()
  for(ds in names(exprs)) {
    iqrs <- unlist(apply(exprs[[ds]], 1, IQR))
    expr.filt.matrices[[ds]] <- exprs[[ds]][iqrs > median(iqrs[iqrs > 0]),]
  }
}

study.genes <- lapply(exprs, rownames)
study.filt.genes <- lapply(expr.filt.matrices, rownames)
common.genes <- Reduce(intersect, study.filt.genes)

## Restrict the matrices to those genes in common across data sets
for(ds in names(fits)) {
  cat(paste0("Subset ", ds, "\n"))
  expr.filt.matrices[[ds]] <- expr.filt.matrices[[ds]][common.genes, ]
  cat(paste0("End subset ", ds, "\n"))
}

venetoclax.feature.syms <- c("BCL2", "CD14", "FLT3", "MAFB", "LRP1")
venetoclax.feature.gene.tbl <- symbols.to.ensg.mapping(venetoclax.feature.syms)
venetoclax.feature.gene.tbl$ensg <- as.character(venetoclax.feature.gene.tbl$ensg)
ensg <- venetoclax.feature.gene.tbl$ensg

cat(paste0("Length of genes: ", length(common.genes), "\n"))

## Do analysis at pathway- (hallmark, biocarta, kegg, and reactome) as well as at gene-level
gene.sets <- list("gene" = "gene", "kegg" = kegg.gene.sets, "hallmark" = hallmark.gene.sets, "biocarta" = biocarta.gene.sets, "reactome" = reactome.gene.sets)
## Do analysis only at gene-level
## gene.sets <- list("gene" = "gene")

## BEGIN modeling setup


## END modeling setup

## BEGIN Define training and test data sets

## Based on whether the data sets are listed in data.sets.to.analyze (defined above), the following code will define
source("define-training-and-test-sets.R")

## END Define training and test data sets

## BEGIN modeling

cat("Training and testing with AUC\n")

## Set the response to use AUC based on 4-parameter log logistic fit (LL.4)
train.response.cols <- rep("dss.auc.ll4", length(train.response.cols))
test.response.cols <- rep("dss.auc.ll4", length(test.response.cols))

res.aucs <- list()
  cat("Randoming OHSU expr\n")

  dummy.col <- "col"
  mean.resp.name <- "mean.resp"
  dummy.table <- data.frame(mean.resp.name)
  colnames(dummy.table) <- dummy.col

  dummy.train.drugs <- list()
  dummy.train.drug.cols <- list()
  dummy.test.drug.cols <- list()

  for(nm in 1:length(train.drug.cols)) {
    dummy.train.drugs[[nm]] <- mean.resp.name
    dummy.train.drug.cols[[nm]] <- dummy.col
  }

  for(nm in 1:length(test.drug.cols)) {
    dummy.test.drug.cols[[nm]] <- dummy.col
  }

  cat("No transform\n")
  prep <- prepare.train.and.test.crossproduct(gene.sets = gene.sets, drug.name.tbl = drug.name.tbl, train.set.names = train.set.names, train.dss.args = train.dss.args, 
                                         train.expr.args = train.expr.args, train.genomic.args = train.genomic.args, train.clinical.args = train.clinical.args, 
                                         train.common.drugs = train.common.drugs, train.drugs = train.drugs, train.drug.cols = train.drug.cols, train.patient.cols = train.patient.cols, 
                                         train.response.cols = train.response.cols, test.set.names = test.set.names, test.dss.args = test.dss.args, test.expr.args = test.expr.args, 
                                         test.genomic.args = test.genomic.args, test.clinical.args = test.clinical.args, test.common.drugs = test.common.drugs,
                                         test.drug.cols = test.drug.cols, test.patient.cols = test.patient.cols, test.response.cols = test.response.cols, 
                                         seed = 1234, use.rf = TRUE, use.svm = FALSE, use.mean = FALSE, num.processes = num.processes, alphas = c(0, 1), include.mean.response.as.feature = FALSE,
                                         subtract.mean.response = FALSE, keep.forest = TRUE, z.score.expr = TRUE, z.score.drc = TRUE)

  ## NB: prep[["train.drcs"]] and prep[["test.drcs"]] have been z-scored.  So, as desired, we will be training to predict the means of the z-scored responses.
  do.z.score <- FALSE
  train.y <- list()
  for(train.set in names(prep[["train.drcs"]])) {
    mean.resp <- t(data.frame("mean.resp" = colMeans(prep[["train.drcs"]][[train.set]], na.rm = TRUE)))
    train.y[[train.set]] <- mean.resp
    if(do.z.score) {
      train.y[[train.set]] <- t(data.frame("mean.resp" = as.vector(t(scale(colMeans(prep[["train.drcs"]][[train.set]], na.rm = TRUE))))))
      colnames(train.y[[train.set]]) <- colnames(prep[["train.drcs"]][[train.set]])
    }
  }

  test.y <- list()
  for(test.set in names(prep[["test.drcs"]])) {
    test.y[[test.set]] <- t(data.frame("mean.resp" = colMeans(prep[["test.drcs"]][[test.set]], na.rm = TRUE)))
    if(do.z.score) {
      test.y[[test.set]] <- t(data.frame("mean.resp" = as.vector(t(scale(colMeans(prep[["test.drcs"]][[test.set]], na.rm = TRUE))))))
      colnames(test.y[[test.set]]) <- colnames(prep[["test.drcs"]][[test.set]])
    }
  }

  res.aucs[["no.transform"]] <- train.and.test.crossproduct_(gene.sets = gene.sets, drug.name.tbl = dummy.table, train.set.names = train.set.names, train.y = train.y, train.x = prep[["train.x"]],
                                         train.drugs = dummy.train.drugs, train.drug.cols = dummy.train.drug.cols, test.set.names = test.set.names, test.y = test.y, test.x = prep[["test.x"]],
                                         test.drug.cols = dummy.test.drug.cols, seed = 1234, use.rf = TRUE, use.svm = FALSE, use.mean = FALSE, num.processes = num.processes, alphas = c(0, 1), 
                                         keep.forest = TRUE, use.glmnet.intercept = TRUE)



cat("Done training and testing with AUC\n")

save.image(rdata.file)

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

cat("Done modeling GLDS\n")