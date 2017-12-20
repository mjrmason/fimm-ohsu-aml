## rm(list = ls())
file.prefix <- "fimm-ohsu-overlap-112017"
rdata.file <- ".Rdata.overlap.112017"

load(rdata.file)

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

cat(paste0("length of common.genes = ", length(common.genes), "\n"))

cat("Defining training table\n")

extract.table <- function(res.aucs, train.set.names, gene.sets, s = "lambda.min") {
 trn.tbl <- ldply(train.set.names,
             .fun = function(train.set) {
                      ldply(names(gene.sets),
                        .fun = function(gene.set) {
                                 ldply(res.aucs[[1]][["all.fits"]][[train.set]][[gene.set]],
                                   .fun = function(lst) {
                                     ldply(lst, 
                                       .fun = function(df) {
                                         ## Only pull out the regression and rf results
                                         if((df$model != "glmnet") && (df$model != "rf")) { return(NULL) }
                                         coeffs <- NULL
                                         if(is.null(df$fit) || is.na(df$fit)) { return(NULL) }
                                         switch(df$model,
                                           "glmnet" = {
                                             coeffs <- coefficients(df$fit, s = s)
                                           },
                                           "rf" = {
                                             col <- colnames(importance(df$fit))[1]
                                             ## NB: %IncMSE = [ mse(j) - mse0 ] / mse0 * 100  [ mse(j): permuted; mse0: original MSE ]
                                             ## i.e., larger/more positive %IncMSE is better
                                             if(!grepl(col, pattern="IncMSE")) { stop(paste0("Was expected ", col, " to be %IncMSE\n")) }
                                             coeffs <- importance(df$fit)[,1,drop=FALSE]
                                           },
                                           { stop(paste0("Unknown model ", df$model, "\n")) })
                                         ns <- rownames(coeffs)
                                         coeffs <- as.vector(coeffs)
                                         names(coeffs) <- ns
                                         coeffs <- coeffs[!grepl(pattern="Intercept", names(coeffs))]
                                         n.features <- length(coeffs) 
                                         coeffs <- coeffs[order(names(coeffs))]
                                         coeff.vals <- paste(coeffs, collapse=",")
                                         coeffs <- paste(names(coeffs), collapse=",")
                                         ret <- data.frame(train.set = train.set, gene.set = gene.set, train.drug = df$train.drug,
                                                           alpha = df$alpha, model = df$model, n.features = n.features, coeffs = coeffs, coeff.vals = coeff.vals)
                                         ret
                                       })
                                   })
                        })
             })
  trn.tbl$model <- as.character(trn.tbl$model)
  trn.tbl$train.drug <- as.character(trn.tbl$train.drug)
  trn.tbl$coeffs <- as.character(trn.tbl$coeffs)
  trn.tbl[(trn.tbl$model == "glmnet") & (trn.tbl$alpha == 0),"model"] <- "ridge"
  trn.tbl[(trn.tbl$model == "glmnet") & (trn.tbl$alpha == 1),"model"] <- "lasso"
  trn.tbl
}

trn.tbl <- extract.table(res.aucs, train.set.names, gene.sets, s = "lambda.min") 
trn.tbl.1se <- extract.table(res.aucs, train.set.names, gene.sets, s = "lambda.1se") 


trn.tbl.glmnet <- subset(trn.tbl, model %in% c("glmnet", "lasso", "ridge"))
trn.tbl <- subset(trn.tbl, model != "glmnet")

trn.tbl.lasso <- subset(trn.tbl, model == "lasso")
trn.tbl.ridge <- subset(trn.tbl, model == "ridge")
trn.tbl.rf <- subset(trn.tbl, model == "rf")

trn.tbl.1se.lasso <- subset(trn.tbl.1se, model == "lasso")
trn.tbl.1se.ridge <- subset(trn.tbl.1se, model == "ridge")
trn.tbl.1se.rf <- subset(trn.tbl.1se, model == "rf")

save.image(rdata.file)

## save.image(".Rdata.overlap.fimm.ohsu")

cat("Done subsetting table by model\n")

## Looking at correlation of coefficients
## method = c("pearson", "spearman", "kendall")
find.overlap.of.dense.predictors <- function(tbl.input, method) {
  ret.tbl <- c()
  for(gs in unique(tbl.input$gene.set)) {
    tmp <- subset(tbl.input, gene.set == gs)
    tr.sets <- as.character(unique(tmp$train.set))
    for(i in 1:(length(tr.sets)-1)) {
      tr1 <- tr.sets[i]
      tmp1 <- subset(tmp, train.set == tr1)
      merge.col <- drug.name.col
      tmp1 <- merge(drug.name.tbl, tmp1, by.x = merge.col, by.y = "train.drug")
      for(j in (i+1):length(tr.sets)) {
        tr2 <- tr.sets[j]
        if(((tr1 == "ohsu") || (tr2 == "ohsu")) && (grepl(pattern="ohsu.set1", paste0(tr1,tr2)) || grepl(pattern="ohsu.set2", paste0(tr1,tr2)))) { next }
        if(((tr1 == "fimm") || (tr2 == "fimm")) && (grepl(pattern="fimm.set1", paste0(tr1,tr2)) || grepl(pattern="fimm.set2", paste0(tr1,tr2)))) { next }
        tmp2 <- subset(tmp, train.set == tr2)
        merge.col <- drug.name.col
        tmp2 <- merge(drug.name.tbl, tmp2, by.x = merge.col, by.y = "train.drug")
  
        m <- merge(tmp1, tmp2, by = merge.col, suffixes = c(".1", ".2"))
        m <- m[, c(merge.col, "n.features.1", "coeffs.1", "coeffs.2", "coeff.vals.1", "coeff.vals.2", "train.set.1", "train.set.2", "gene.set.1")]
        m$id <- 1:nrow(m)
        hyp <- ddply(m, .variables = "id",
               .fun = function(r) {
                        coeff.vals.1 <- as.numeric(unlist(strsplit(as.character(r$coeff.vals.1), split=",")))
                        coeffs.1 <- (unlist(strsplit(as.character(r$coeffs.1), split=",")))
                        names(coeff.vals.1) <- coeffs.1
                        coeffs.1 <- coeff.vals.1
                        coeff.vals.2 <- as.numeric(unlist(strsplit(as.character(r$coeff.vals.2), split=",")))
                        coeffs.2 <- (unlist(strsplit(as.character(r$coeffs.2), split=",")))
                        names(coeff.vals.2) <- coeffs.2
                        coeffs.2 <- coeff.vals.2
                        inter <- intersect(names(coeffs.1), names(coeffs.2))
                        ct <- cor.test(coeffs.1[inter], coeffs.2[inter], method = method)
                        pval <- ct$p.value
                        cor <- ct$estimate
                        vec <- c(r$id[1], cor, pval)
                        names(vec) <- c("id", "cor", "pval")
                        vec
                      })
        hyp <- merge(hyp, m, by = "id")
        ret.tbl <- rbind(ret.tbl, hyp)
      }
    }
  }
  ret.tbl
}


source("sparse-overlap.R")

list.sizes <- seq(from=10,to=5000,by=10)
list.sizes <- c(10,50,100,1000)
list.sizes <- c(10,100,1000)
list.sizes <- unique(c(seq(from=2, to=10, by=1), seq(from=10,to=1000,by=10), seq(from=1000, to=5000, by=100)))
list.sizes <- seq(from=10,to=10000,by=100)
list.sizes <- unique(c(seq(from=10,to=2000,by=10), seq(from=2000, to=10000, by=100)))
names(list.sizes) <- list.sizes

num.rf.sigs <- NULL
if(use.rf) {
  cat("Summarizing rf results\n")
  ## all.rf <- find.overlap.of.dense.predictors(trn.tbl.rf, method = "spearman")

  ## these are just for debugging
  tmp.rf <- subset(trn.tbl.rf, train.set %in% c("ohsu", "fimm"))
  for(sz in c(10,100,1000,2000,5000,7000,7500,7600,length(common.genes)-100,length(common.genes)-10,length(common.genes)-1,length(common.genes),length(common.genes)+1,8000)) {
cat(paste0("sz = ", sz, "\n"))
                                tmp <- find.overlap.of.sparse.predictors(tmp.rf, universe = common.genes, transform = "none", top = sz, verbose = TRUE)
  }

  all.rf <- find.overlap.of.sparse.predictors(tmp.rf, universe = common.genes, transform = "none", top = 100)
  all.rf$pval <- as.numeric(all.rf$pval)
  cat("Done summarizing rf results\n")
  num.rf.sigs <- llply(list.sizes, .parallel = TRUE,
                       .fun = function(sz) {
                                tmp <- find.overlap.of.sparse.predictors(tmp.rf, universe = common.genes, transform = "none", top = sz, verbose = FALSE)
                                tmp$pval <- as.numeric(tmp$pval)
                                all.sig <- subset(tmp, pval < 0.05)
                                if(nrow(all.sig) == 0) { return(0) }
                                cross.sig <- subset(all.sig, train.set.1 == "ohsu" & train.set.2 == "fimm")
                                if(nrow(cross.sig) == 0) { return(0) }
                                cross.sig <- subset(cross.sig, gene.set.1 == "gene")
                                nrow(cross.sig)
                              })
  pdf(paste0(file.prefix, "-fimm-vs-ohsu-", "rf", "-feature-overlap-vs-num-features.pdf"))
  plot(as.numeric(names(num.rf.sigs)), num.rf.sigs)
  d <- dev.off()
  cat("Summarizing rf results as a function of list size\n")
  cat("Done summarizing rf results as a function of list size\n")
}
save.image(rdata.file)

cat("Summarizing sparse glmnet results as a function of list size\n")


for(glmnet.alpha in glmnet.alphas) {
  cat(paste0("Summarizing glmnet alpha = ", glmnet.alpha, "\n"))
  pdf(paste0(file.prefix, "-fimm-vs-ohsu-", "glmnet-alpha-", glmnet.alpha, "-feature-overlap-vs-num-features.pdf"))
  tmp.alpha <- subset(trn.tbl.glmnet, alpha == glmnet.alpha)
  tmp.alpha <- subset(tmp.alpha, train.set %in% c("ohsu", "fimm"))
  num.sigs <- llply(list.sizes, .parallel = TRUE,
                          .fun = function(sz) {
                                   tmp <- find.overlap.of.sparse.predictors(tmp.alpha, universe = common.genes, transform = "abs.drop.zero", top = sz)
                                   tmp$pval <- as.numeric(tmp$pval)
                                   all.sig <- subset(tmp, pval < 0.05)
                                   if(nrow(all.sig) == 0) { return(0) }
                                   cross.sig <- subset(all.sig, train.set.1 == "ohsu" & train.set.2 == "fimm")
                                   if(nrow(cross.sig) == 0) { return(0) }
                                   cross.sig <- subset(cross.sig, gene.set.1 == "gene")
                                   nrow(cross.sig)
                                 })
  plot(as.numeric(names(num.sigs)), num.sigs, xlab = "Number of Top Features", ylab = "Number of Drugs with Overlapping Features", main = paste0("GLMNET: alpha = ", glmnet.alpha))
  d <- dev.off()
}
cat("Done summarizing sparse glmnet results as a function of list size\n")

if(FALSE) {
cat("Summarizing sparse ridge results\n")
num.ridge.sigs <- llply(list.sizes, .parallel = TRUE,
                        .fun = function(sz) {
                                 tmp <- find.overlap.of.sparse.predictors(trn.tbl.ridge, universe = common.genes, transform = "abs.drop.zero", top = sz)
                                 tmp$pval <- as.numeric(tmp$pval)
                                 all.sig <- subset(tmp, pval < 0.05)
                                 cross.sig <- subset(all.sig, train.set.1 == "ohsu" & train.set.2 == "fimm")
                                 length(which(cross.sig$gene.set.1 == "gene"))
                               })
pdf(paste0(file.prefix, "-fimm-vs-ohsu-", "ridge", "-feature-overlap-vs-num-features.pdf"))
plot(as.numeric(names(num.ridge.sigs)), num.ridge.sigs)
d <- dev.off()
cat("Done summarizing sparse ridge results\n")
}

if(FALSE) {
  cat("Summarizing dense ridge results\n")
  all.ridge <- find.overlap.of.dense.predictors(trn.tbl.ridge, method = "spearman")
  all.ridge$pval <- as.numeric(all.ridge$pval)
  cat("Done summarizing dense ridge results\n")

lst <- list("ridge" = all.ridge)
all.cross.sigs <- list()
for(nm in names(lst)) {

  all.sig <- subset(lst[[nm]], pval < 0.05 & cor > 0)
  cross.sig <- subset(all.sig, train.set.1 == "ohsu" & train.set.2 == "fimm")
  ohsu.sig <- subset(all.sig, (train.set.1 == "ohsu.set1" & train.set.2 == "ohsu.set2") | (train.set.2 == "ohsu.set1" & train.set.1 == "ohsu.set2"))
  fimm.sig <- subset(all.sig, (train.set.1 == "fimm.set1" & train.set.2 == "fimm.set2") | (train.set.2 == "fimm.set1" & train.set.1 == "fimm.set2"))
  all.cross.sigs[[nm]] <- cross.sig

  tbl <- as.data.frame(table(cross.sig$gene.set.1))
  colnames(tbl) <- c("Feature Set", "Drugs w/ Correlated Features")
  pdf(paste0(file.prefix, "-fimm-vs-ohsu-", nm, "-feature-overlap.pdf"))
  plot.table(tbl, main = "OHSU (train) vs FIMM (train)")
  d <- dev.off()

  tbl <- as.data.frame(table(ohsu.sig$gene.set.1))
  colnames(tbl) <- c("Feature Set", "Drugs w/ Correlated Features")
  pdf(paste0(file.prefix, "-ohsu1-vs-ohsu2-", nm, "-feature-overlap.pdf"))
  plot.table(tbl, main = "OHSU 1 (train) vs OHSU 2 (train)")
  d <- dev.off()

  tbl <- as.data.frame(table(fimm.sig$gene.set.1))
  colnames(tbl) <- c("Feature Set", "Drugs w/ Correlated Features")
  pdf(paste0(file.prefix, "-fimm1-vs-fimm2-", nm, "-feature-overlap.pdf"))
  plot.table(tbl, main = "FIMM 1 (train) vs FIMM 2 (train)")
  d <- dev.off()
}

} ## FALSE 

## save.image(".Rdata.overlap.fimm.ohsu")

all <- find.overlap.of.sparse.predictors(trn.tbl.lasso, universe = common.genes, transform = "drop.zero", top = 0)
all$pval <- as.numeric(all$pval)
cat("Done summarizing lasso results\n")

all.1se <- find.overlap.of.sparse.predictors(trn.tbl.1se.lasso, universe = common.genes, transform = "drop.zero", top = 0)
all.1se$pval <- as.numeric(all.1se$pval)
cat("Done summarizing lasso results\n")

## save.image(".Rdata.overlap.fimm.ohsu")

lst <- list("lasso" = all, "lasso.1se" = all.1se)
if(use.rf) {
  lst <- list("lasso" = all, "lasso.1se" = all.1se, "rf" = all.rf)
}
all.cross.sigs <- list()
for(nm in names(lst)) {
  all.sig <- subset(lst[[nm]], pval < 0.05)
  cross.sig <- subset(all.sig, train.set.1 == "ohsu" & train.set.2 == "fimm")
  all.cross <- subset(lst[[nm]], train.set.1 == "ohsu" & train.set.2 == "fimm")
  ohsu.sig <- subset(all.sig, (train.set.1 == "ohsu.set1" & train.set.2 == "ohsu.set2") | (train.set.2 == "ohsu.set1" & train.set.1 == "ohsu.set2"))
  fimm.sig <- subset(all.sig, (train.set.1 == "fimm.set1" & train.set.2 == "fimm.set2") | (train.set.2 == "fimm.set1" & train.set.1 == "fimm.set2"))
  all.cross.sigs[[nm]] <- cross.sig

  table(cross.sig$gene.set.1)
  table(ohsu.sig$gene.set.1)
  table(fimm.sig$gene.set.1)

  tbl <- as.data.frame(table(cross.sig$gene.set.1))
  colnames(tbl) <- c("Feature Set", "Drugs w/ Overlapping Features")
  pdf(paste0(file.prefix, "-ohsu-vs-fimm-", nm , "-feature-overlap.pdf"))
  plot.table(tbl, main = "OHSU (train) vs FIMM (train)")
  d <- dev.off()

  if(grepl(x=nm, pattern="lasso")) {  
    s <- subset(all.cross, n.1 > 0 & n.2 > 0)
    tmp <- s[,c("FIMM_DRUG_NAME", "train.set.1", "n.1", "train.set.2", "n.2", "gene.set.1", "n.both")]
    write.table(file = paste0(file.prefix, "-ohsu-vs-fimm-", nm , "-feature-overlap.tsv"), tmp, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  }

  tbl <- as.data.frame(table(ohsu.sig$gene.set.1))
  colnames(tbl) <- c("Feature Set", "Drugs w/ Overlapping Features")
  pdf(paste0(file.prefix, "-ohsu1-vs-ohsu2-", nm , "-feature-overlap.pdf"))
  plot.table(tbl, main = "OHSU 1 (train) vs OHSU 2 (train)")
  d <- dev.off()

  tbl <- as.data.frame(table(fimm.sig$gene.set.1))
  colnames(tbl) <- c("Feature Set", "Drugs w/ Overlapping Features")
  pdf(paste0(file.prefix, "-fimm1-vs-fimm2-", nm , "-feature-overlap.pdf"))
  plot.table(tbl, main = "FIMM 1 (train) vs FIMM 2 (train)")
  d <- dev.off()
}

## save.image(rdata.file)
cat("successfully completed")
