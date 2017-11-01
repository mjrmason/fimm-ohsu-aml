library(synapseClient)
library(data.table)
library(ggplot2)
library(gridExtra)

suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(caret))
suppressPackageStartupMessages(library(glmnet))
library( ReporteRs )
library(scales)

suppressPackageStartupMessages(library("parallel"))
num.cores <- detectCores()
if(!is.na(num.cores) && (num.cores > 1)) {
  suppressPackageStartupMessages(library("doMC"))
  cat(paste("Registering ", num.cores-1, " cores.\n", sep=""))
  registerDoMC(cores=(num.cores-1))
}

synapseLogin()

## Load in the FIMM and OHSU fits.

## path <- "fimm.dss.t0.tsv"
synId <- "syn10083488"
obj <- synGet(synId, downloadFile = TRUE)
file <- getFileLocation(obj)
fimm.dss.t0 <- as.data.frame(fread(file))

## The e parameter in both L.4 and LL.4 is IC50 (not log IC50).
fimm.dss.t0$ic50.ll4 <- fimm.dss.t0$e.ll4
fimm.dss.t0$ic50.l4 <- fimm.dss.t0$e.l4

## path <- "fimm.dss.t10.tsv"
synId <- "syn10083490"
obj <- synGet(synId, downloadFile = TRUE)
file <- getFileLocation(obj)
fimm.dss.t10 <- as.data.frame(fread(file))
fimm.dss.t10$ic50.ll4 <- fimm.dss.t10$e.ll4
fimm.dss.t10$ic50.l4 <- fimm.dss.t10$e.l4

## path <- "ohsu.dss.t0.tsv"
synId <- "syn10083494"
obj <- synGet(synId, downloadFile = TRUE)
file <- getFileLocation(obj)
ohsu.dss.t0 <- as.data.frame(fread(file))
ohsu.dss.t0$ic50.ll4 <- ohsu.dss.t0$e.ll4
ohsu.dss.t0$ic50.l4 <- ohsu.dss.t0$e.l4

## path <- "ohsu.dss.t10.tsv"
synId <- "syn10083499"
obj <- synGet(synId, downloadFile = TRUE)
file <- getFileLocation(obj)
ohsu.dss.t10 <- as.data.frame(fread(file))
ohsu.dss.t10$ic50.ll4 <- ohsu.dss.t10$e.ll4
ohsu.dss.t10$ic50.l4 <- ohsu.dss.t10$e.l4

## path <- "ohsu.fimm.drugs.tsv"
synId <- "syn10083888"
obj <- synGet(synId, downloadFile = TRUE)
file <- getFileLocation(obj)
ohsu.fimm.drugs <- as.data.frame(fread(file))

## Read in the OHSU expression data
## path <- "ohsu.expr.tsv"
synId <- "syn10083723"
obj <- synGet(synId, downloadFile = TRUE)
file <- getFileLocation(obj)
ohsu.expr <- read.table(file, header=TRUE, sep="\t")
## Convert the column names from X20.00347 to 20-00347
colnames(ohsu.expr) <- gsub("X(.*)\\.(.*)", "\\1-\\2", colnames(ohsu.expr))

## Read in the rna-seq summary table, which provides a map from the lab_id's used
## in the drug response data to the SeqID's used in the expression data.
synId <- "syn10083817"
obj <- synGet(synId, downloadFile = TRUE)
file <- getFileLocation(obj)
ohsu.rnaseq.sample.summary <- read.table(file, header=TRUE, sep="\t")

## For each drug in common between OHSU and FIMM, fit elastic net to response (auc) and bound
## performance via bootstrapping.
## Specifically:
## 1. For each of 100 bootstraps
## 2.   Fit elastic net (or lasso or ridge) indepedently for each drug
## 3.   Predict on test samples
## 4. Get distribution of predicted drug response (AUCs) for each sample and each drug
## 5. For each drug
## 6.   Plot distribution of predicted drug response for each sample (x axis) vs actual,
##      include correlation line between median predicted and actual (Pearson and Spearman)
## 8. Make a single plot of Pearson (and Spearman) correlation vs drug (x axis)

## ohsu.fimm.drugs$`Mechanism/Targets`
## Selumetinib (AZD6244)
## Trametinib (GSK1120212)
## CI-1040 (PD184352)

source("analysis.R")

## mek.inhibitors <- ohsu.fimm.drugs$ID_Drug.ohsu[ohsu.fimm.drugs$`Mechanism/Targets` == "MEK1/2 inhibitor"]
## lst <- fit.elastic.net.ohsu(ohsu.dss.t0, ohsu.expr, ohsu.rnaseq.sample.summary, mek.inhibitors, "auc.l4", regression.method = "lasso", seed = 1234)

## 
all.common.drugs <- na.omit(ohsu.fimm.drugs$ID_Drug.ohsu)
mek.inhibitors <- ohsu.fimm.drugs$ID_Drug.ohsu[ohsu.fimm.drugs$`Mechanism/Targets` == "MEK1/2 inhibitor"]
drugs <- all.common.drugs
## drugs <- mek.inhibitors
## lst <- fit.elastic.net.ohsu(ohsu.dss.t0, ohsu.expr, ohsu.rnaseq.sample.summary, drugs, "auc.l4", regression.method = "lasso", seed = 1234)

## save.image(".Rdata")

## lst <- prepare.ohsu.drug.response.and.expr.matrices(ohsu.dss.t0, ohsu.expr, ohsu.rnaseq.sample.summary, drugs, response.col = "auc.l4")

## save.image(".Rdata")



cat("Done preparing\n")

## lst <- bootstrap.elastic.net.ohsu(ohsu.dss.t0, ohsu.expr, ohsu.rnaseq.sample.summary, drugs, response.col = "auc.l4", alphas = c(0, 0.25, 0.5, 0.75, 1), num.bootstraps = 10)

lst <- bootstrap.elastic.net.ohsu(ohsu.dss.t0, ohsu.expr, ohsu.rnaseq.sample.summary, drugs, response.col = "auc.l4", alphas = c(0, 0.25, 0.5, 0.75, 1), num.bootstraps = 50)

cat("Done with bootstrap\n")

save.image(".Rdata")

## Flatten the predicted responses
predicted.responses <- ldply(1:length(lst), .parallel = FALSE, 
                             .fun = function(drug.i) {
                                      drug.lst <- lst[[drug.i]]
                                      ret.drug <- ldply(drug.lst, .parallel = FALSE,
                                                        .fun = function(fit) {
                                                                 ret.alpha <- ldply(1:length(fit$alphas), .parallel = FALSE,
                                                                                    .fun = function(i) {
                                                                                             resps <- fit$predictions[[i]][,1]
                                                                                             samples <- rownames(fit$predictions[[i]])
if(length(resps) <= 1) { return(c(alpha = fit$alphas[i], predicted.response = NA, predicted.samples = NA)) }
                                                                                             response <- data.frame(sample = samples, response = resps)
                                                                                             resp.str <- paste0(resps, collapse=",")
                                                                                             sample.str <- paste0(samples, collapse=",")
                                                                                             vec <- c(fit$alphas[i], resp.str, sample.str)
                                                                                             names(vec) <- c("alpha", "predicted.response", "predicted.samples")
                                                                                             vec
                                                                                    })
                                                                 ret.alpha$drug <- fit$drug
                                                                 ret.alpha$boot <- fit$boot
                                                                 ret.alpha
                                                               })
                                      ret.drug
                                    })

## Flatten the actual responses
actual.responses <- ldply(1:length(lst), .parallel = FALSE, 
                             .fun = function(drug.i) {
                                      drug.lst <- lst[[drug.i]]
                                      ret.drug <- ldply(drug.lst, .parallel = FALSE,
                                                        .fun = function(fit) {
                                                                 resps <- as.vector(t(fit$test.response))
                                                                 samples <- colnames(fit$test.response)
if(length(resps) <= 1) { return(c(drug = fit$drug, actual.response = NA, actual.samples = NA)) }
                                                                 resp.str <- paste0(resps, collapse=",")
                                                                 sample.str <- paste0(samples, collapse=",")
                                                                 vec <- c(fit$drug, fit$boot, resp.str, sample.str)
                                                                 names(vec) <- c("drug", "boot", "actual.response", "actual.samples")
                                                                 vec
                                                               })
                                      ret.drug
                                    })

actual.and.predicted.responses <- merge(predicted.responses, actual.responses)
actual.and.predicted.responses <- subset(actual.and.predicted.responses, !is.na(predicted.response))
actual.and.predicted.responses <- subset(actual.and.predicted.responses, !is.na(actual.response))

save.image(".Rdata")
cat("Done with assembling responses\n")

## 5. For each drug
## 6.   Calculate Pearson and Spearman correlations for actual drug response (AUC) vs median (within each sample) predicted auc
## 7.   Plot each of these
## Create a short-form matrix having columns response.type (= "predicted", "actual"), sample, response, drug, bootstrap, and alpha
short.df <- 
ldply(1:nrow(actual.and.predicted.responses), .parallel = FALSE,
      .fun = function(i) {
               drug <- actual.and.predicted.responses$drug[i]
               bootstrap <- actual.and.predicted.responses$bootstrap[i]
               alpha <- actual.and.predicted.responses$alpha[i]
               vec <- unlist(strsplit(actual.and.predicted.responses$predicted.response[i], split=","))
               name.vec <- unlist(strsplit(actual.and.predicted.responses$predicted.samples[i], split=","))
               pred.ret <- ldply(1:length(vec), .parallel = FALSE,
                                .fun = function(vec.i) {
                                         res <- c(sample = name.vec[vec.i], response = vec[vec.i], response.type = "predicted")
                                         res
                                })
               vec <- unlist(strsplit(actual.and.predicted.responses$actual.response[i], split=","))
               name.vec <- unlist(strsplit(actual.and.predicted.responses$actual.samples[i], split=","))
               actual.ret <- ldply(1:length(vec), .parallel = FALSE,
                                 .fun = function(vec.i) {
                                          res <- c(sample = name.vec[vec.i], response = vec[vec.i], response.type = "actual")
                                          res
                                 })
               ret <- rbind(pred.ret, actual.ret)
               ret$drug <- drug
               ret$bootstrap <- bootstrap
               ret$alpha <- alpha
               ret
      })

do.plot <- FALSE
actual.predicted.correlations.df <- 
ddply(short.df, .variables = c("alpha", "drug"), .parallel = FALSE,
      .fun = function(df) {
                df$response <- as.numeric(df$response)
                predicted.flag <- df$response.type == "predicted"
                median.pred <- ddply(df[predicted.flag,], .variables = "sample",
                                     .fun = function(df.sample) {
                                              vec <- c(predicted.response = median(as.numeric(df.sample$response)), sample = df.sample$sample[1])
                                              vec
                                     })
                actual.df <- unique(df[!predicted.flag, c("sample", "response", "alpha", "drug")])
                colnames(actual.df) <- c("sample", "actual.response", "alpha", "drug")
                m <- merge(median.pred, actual.df)
                m$predicted.response <- as.numeric(m$predicted.response)
                m$actual.response <- as.numeric(m$actual.response)
                if(do.plot) {
                  g <- ggplot(data = m, aes(x = actual.response, y = predicted.response))
                  g <- g + geom_point()
                  g <- g + ggtitle(paste0("Drug: ", df$drug[1], " Alpha: ", df$alpha[1]))
                  print(g)
                }
                ct.spearman <- cor.test(m$actual.response, m$predicted.response, method="spearman")
                ct.pearson <- cor.test(m$actual.response, m$predicted.response, method="pearson")
                vec <- c(ct.spearman$p.value, unname(ct.spearman$estimate), ct.pearson$p.value, unname(ct.pearson$estimate))
                names(vec) <- c("spearman.p", paste0("spearman.", names(ct.spearman$estimate)), "pearson.p", paste0("pearson.", names(ct.pearson$estimate)))
                vec
      })

df2 <- 
  ddply(short.df, .variables = c("alpha", "drug"), .parallel = FALSE,
        .fun = function(df) {
          stop("This does not use median")
          df$response <- as.numeric(df$response)
          predicted.flag <- df$response.type == "predicted"
          median.pred <- ddply(df[predicted.flag,], .variables = "sample",
                               .fun = function(df.sample) {
                                 vec <- c(predicted.response = median(as.numeric(df.sample$response)), sample = df.sample$sample[1])
                                 vec
                               })
          actual.df <- unique(df[!predicted.flag, c("sample", "response", "alpha", "drug")])
          colnames(actual.df) <- c("sample", "actual.response", "alpha", "drug")
          predicted.df <- df[predicted.flag, c("sample", "response", "alpha", "drug")]
          colnames(predicted.df) <- c("sample", "predicted.response", "alpha", "drug")
          
          m <- merge(predicted.df, actual.df)
          m$predicted.response <- as.numeric(m$predicted.response)
          m$actual.response <- as.numeric(m$actual.response)
          if(do.plot) {
            g <- ggplot(data = m, aes(x = actual.response, y = predicted.response))
            g <- g + geom_point()
            g <- g + ggtitle(paste0("Drug: ", df$drug[1], " Alpha: ", df$alpha[1]))
            print(g)
          }
          ct.spearman <- cor.test(m$actual.response, m$predicted.response, method="spearman")
          ct.pearson <- cor.test(m$actual.response, m$predicted.response, method="pearson")
          vec <- c(ct.spearman$p.value, unname(ct.spearman$estimate), ct.pearson$p.value, unname(ct.pearson$estimate))
          names(vec) <- c("spearman.p", paste0("spearman.", names(ct.spearman$estimate)), "pearson.p", paste0("pearson.", names(ct.pearson$estimate)))
          vec
        })

actual.predicted.correlations.df$spearman.ap <- rep(NA, nrow(actual.predicted.correlations.df))
actual.predicted.correlations.df$pearson.ap <- rep(NA, nrow(actual.predicted.correlations.df))
for(alpha in unique(actual.predicted.correlations.df$alpha)) {
  alpha.flag <- actual.predicted.correlations.df$alpha == alpha
  actual.predicted.correlations.df$spearman.ap[alpha.flag] <- p.adjust(actual.predicted.correlations.df$spearman.p[alpha.flag], method = "bonferroni")
  actual.predicted.correlations.df$pearson.ap[alpha.flag] <- p.adjust(actual.predicted.correlations.df$pearson.p[alpha.flag], method = "bonferroni")
}

pearson.num.sig <- ddply(actual.predicted.correlations.df, .variables = "alpha",
                         .fun = function(df) {
                            alpha <- unique(df$alpha)
                            num.sig <- length(which(df$pearson.ap < 0.05))
                            tot <- nrow(df)
                            c(alpha = alpha, num.sig = num.sig, tot.tested = tot)
                         })
pearson.num.sig$test <- "Pearson"

spearman.num.sig <- ddply(actual.predicted.correlations.df, .variables = "alpha",
                         .fun = function(df) {
                           alpha <- unique(df$alpha)
                           num.sig <- length(which(df$spearman.ap < 0.05))
                           tot <- nrow(df)
                           c(alpha = alpha, num.sig = num.sig, tot.tested = tot)
                         })
spearman.num.sig$test <- "Spearman"
num.sig <- rbind(pearson.num.sig, spearman.num.sig)

g <- ggplot(data = num.sig, aes(x = alpha, y = num.sig))
g <- g + geom_point()
g <- g + facet_wrap(~ test)
g <- g + xlab("")
g <- g + ylab("Number Significant (Bonferonni-corrected p < 0.05)")
g <- g + ggtitle("Number of Drugs (n = 79) with\nSignificant Prediction vs Median Response")
pdf("num-significant-drugs.pdf")
print(g)
d <- dev.off()

##cat("Bootstrapping 50\n")
##lst2 <- bootstrap.elastic.net.ohsu(ohsu.dss.t0, ohsu.expr, ohsu.rnaseq.sample.summary, drugs, response.col = "auc.l4", alphas = c(0, 0.25, 0.5, 0.75, 1), num.bootstraps = 50)

##cat("Done with bootstrap\n")

plts <- 
dlply(actual.predicted.correlations.df, .variables = "alpha", .parallel = FALSE,
     .fun = function(tbl) {
              alpha <- unique(tbl$alpha)
              pval.col <- "pearson.p"
              tbl <- tbl[order(tbl[,pval.col]),]
              x.col <- "drug"
              tbl[,x.col] <- factor(tbl[,x.col], levels=unique(tbl[,x.col]))
              file <- paste0("pearson-alpha-", alpha, ".pdf")
              pdf(file, onefile=FALSE)
              g1 <- plot.stacked.pval.correlation.figures(tbl, "pearson.p", "pearson.cor", "drug", y.cor.lab = "Pearson's Correlation", main = paste0("Alpha = ", alpha), text.size = 5)
              plot(g1)
              d <- dev.off()
              file <- paste0("spearman-alpha-", alpha, ".pdf")
              pdf(file, onefile=FALSE)
              g2 <- plot.stacked.pval.correlation.figures(tbl, "spearman.p", "spearman.rho", "drug", y.cor.lab = "Spearman's Rho", main = paste0("Alpha = ", alpha), text.size = 5)
              plot(g2)
              d <- dev.off()
              return(list(g1, g2))
           })

## Plot the distribution of correlation values as a function of alpha
g1 <- ggplot(data = actual.predicted.correlations.df, aes(x = factor(alpha), y = spearman.rho))
g1 <- g1 + geom_violin()
g1 <- g1 + ggtitle("Spearman's rho")
g1 <- g1 + geom_boxplot(width = 0.5)
g1 <- g1 + xlab("alpha")
g1 <- g1 + ylab("Spearman's rho")

g2 <- ggplot(data = actual.predicted.correlations.df, aes(x = factor(alpha), y = pearson.cor))
g2 <- g2 + geom_violin()
g2 <- g2 + ggtitle("Pearson's correlation")
g2 <- g2 + geom_boxplot(width = 0.5)
g2 <- g2 + xlab("alpha")
g2 <- g2 + ylab("Pearson's r")

pdf("correlation-vs-alpha.pdf")
grid.arrange(g1, g2)
d <- dev.off()

corr.alpha.g <- arrangeGrob(grobs = list(g1, g2))

tmp <- subset(actual.predicted.correlations.df, alpha == 0)
tmp <- tmp[order(tmp$pearson.p,decreasing=FALSE),]

individual.drug.fits <- list()
for(i in 1:4) {
  target.alpha <- tmp$alpha[i]
  target.drug <- tmp$drug[i]
  pdf(paste0(target.drug, "-alpha-", target.alpha, ".pdf"))
  drug.g <- plot.predicted.vs.actual.response(short.df, target.alpha = target.alpha, target.drug = target.drug)
  individual.drug.fits[[target.drug]] <- drug.g
  print(drug.g)
  d <- dev.off()
}

mydoc <- pptx()

mydoc <- addSlide( mydoc, slide.layout = 'Title Slide' )
mydoc <- addTitle( mydoc, 'Elastic Net Analysis of OHSU (AUC vs expr)' )
mydoc <- addSubtitle( mydoc , 'June 28, 2017')

mydoc <- addSlide( mydoc, slide.layout = "Title and Content" )
## title <- "FIMM (t = 0)"
## mydoc <- addTitle( mydoc, title )
mydoc <- addPlot( doc = mydoc, fun = plot, x = corr.alpha.g)

for(i in 1:length(plts)) {
  mydoc <- addSlide( mydoc, slide.layout = "Two Content" )
  ## title <- paste0("OHSU (t = ", dss.t, "): ", mechanism)
  ## mydoc <- addTitle( mydoc, title )
  mydoc <- addPlot( doc = mydoc, fun = plot, x = plts[[i]][[1]] )
  mydoc <- addPlot( doc = mydoc, fun = plot, x = plts[[i]][[2]] )
}

for(drug in names(individual.drug.fits)) {
  mydoc <- addSlide( mydoc, slide.layout = "Title and Content" )
  title <- paste0("Drug: ", drug, " Alpha: ", 0)
  mydoc <- addTitle( mydoc, title )
  mydoc <- addPlot( doc = mydoc, fun = plot, x = individual.drug.fits[[drug]])
}


writeDoc(mydoc, "062817-elastic-net-ohsu.pptx")

save.image(".Rdata")
cat("Done with creating correlation table\n")
q(status=0)



translate.ensg.to.symbol <- function(ids) {
  library(biomaRt)
  # Get a mapping from Ensembl id to symbol
  ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
  bm <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol'),
              filters = 'ensembl_gene_id', 
              values = ids, 
              mart = ensembl)
  names(bm) <- c("ID", "TRANSLATION")
  bm <- bm[!(bm$TRANSLATION %in% c("")),]
  bm
}

ensg.to.sym.tbl <- translate.ensg.to.symbol(rownames(ohsu.expr))

## For each drug with significant correlation between prediction and response:
## NB: "significance" for glmnet seems not to be well defined.  e.g., 
## https://stats.stackexchange.com/questions/34859/how-to-present-results-of-a-lasso-using-glmnet
## https://stackoverflow.com/questions/12937331/why-is-it-inadvisable-to-get-statistical-summary-information-for-regression-coef/17725220#17725220
## 1.  Plot prediction vs response
## 2.  Plot gene expression vs response
library( ReporteRs )
mydoc <- pptx()

mydoc <- addSlide( mydoc, slide.layout = 'Title Slide' )
mydoc <- addTitle( mydoc, 'LASSO analysis of OHSU (expr vs AUC)' )
mydoc <- addSubtitle( mydoc , 'June 20, 2017')

for(indx in 1:length(lst)) {
  drug <- lst[[indx]]$drug
  test.response <- lst[[indx]]$test.response
  predictions <- lst[[indx]]$predictions
  predictions <- as.numeric(predictions[colnames(test.response),1])
  model <- lst[[indx]]$model
  mt <- as.matrix(coef(lst[[indx]]$model))
  nz.genes <- rownames(mt[mt[,1] > 0,,drop=FALSE])
  nz.genes <- nz.genes[!grepl(pattern="intercept", ignore.case=TRUE, x = nz.genes)]
  if(length(nz.genes) > 0) {
    cat(paste0("Genes associated with drug = ", drug, ": ", paste(nz.genes, collapse=", "), "\n"))
    df <- data.frame(predictions = predictions, test.response = as.numeric(test.response))
    g <- ggplot(data = df, aes(x = predictions, y = test.response))
    g <- g + geom_point()
    g <- g + xlab("Predicted drug response (z-scored AUC)")
    g <- g + ylab("Actual drug response (z-scored AUC)")
    g <- g + ggtitle(drug)

    ## Plot the correlation of each gene expression with actual drug response
    glist <- list(g)
    for(gene in nz.genes) {
      df <- data.frame(expr = as.numeric(ohsu.expr[gene, colnames(test.response)]), test.response = as.numeric(test.response))
      g <- ggplot(data = df, aes(x = expr, y = test.response))
      g <- g + geom_point()
      g <- g + xlab(paste0(gene, " expression (z-scored)"))
      g <- g + ylab("Drug response (z-scored AUC)")
      symbols <- gene
      if(gene %in% ensg.to.sym.tbl$ID) {
        symbols <- paste(ensg.to.sym.tbl$TRANSLATION[ensg.to.sym.tbl$ID == gene], collapse=",")
      }
  
      fit <- lm(data = df, test.response ~ expr)
      r2 <- round(summary(fit)$r.squared, 2)
##      eqn <- bquote(r^2 == .(r2) * "," ~~ RMSE == .(rmse))
      eqn <- r2
      eqn <- substitute(italic(r)^2~"="~r2, 
                       list(r2 = format(summary(fit)$r.squared, digits = 3)))
      eqn <- as.character(as.expression(eqn));
      d.eqn <- data.frame(label = eqn, x = min(df$expr,na.rm=TRUE) + 0.1 * (max(df$expr, na.rm=TRUE) - min(df$expr, na.rm=TRUE)), y = max(df$test.response,na.rm=TRUE) - 0.1 * (max(df$test.response, na.rm=TRUE) - min(df$test.response, na.rm=TRUE)))
      g <- g + geom_text(data = d.eqn, aes(x = x, y = y, label = label), parse=TRUE)

      g <- g + ggtitle(symbols)
      g <- g + geom_smooth(method = "lm", se = FALSE)
      glist[[length(glist)+1]] <- g
    }
    pdf(paste0(drug, "-prediction-vs-response.pdf"))
    ## do.call("grid.arrange", c(glist, top = drug))
    ag <- do.call("arrangeGrob", c(glist, top = drug, ncol = 2))
    plot(ag)
    d <- dev.off()
    mydoc <- addSlide( mydoc, slide.layout = "Title and Content" )
    title <- paste0(drug, ": ", ohsu.fimm.drugs$`Mechanism/Targets`[ohsu.fimm.drugs$ID_Drug.ohsu == drug])
    mydoc <- addTitle( mydoc, title )
    
##    mydoc <- addPlot( doc = mydoc, fun = print, x = g)
    mydoc <- addPlot( doc = mydoc, fun = plot, x = ag)
  }
}

## Make a heatmap of drug (row; annotated by drug class) vs participating genes (col)
mat <- ldply(1:length(lst), 
             .fun = function(indx) {
                      drug <- lst[[indx]]$drug
                      model <- lst[[indx]]$model
                      mt <- as.matrix(coef(lst[[indx]]$model))
                      nz.genes <- rownames(mt[mt[,1] > 0,,drop=FALSE])
                      nz.genes <- nz.genes[!grepl(pattern="intercept", ignore.case=TRUE, x = nz.genes)]
                      if(length(nz.genes > 0)) {
                        genes <- unlist(lapply(nz.genes, function(gene) { ifelse(gene %in% ensg.to.sym.tbl$ID, ensg.to.sym.tbl$TRANSLATION[ensg.to.sym.tbl$ID == gene], gene) }))
                        return(data.frame(drug = drug, gene = genes))
                      }
                    })
mat <- spread(mat, gene, gene)
rownames(mat) <- mat$drug
mat <- as.matrix(mat[, !(colnames(mat) %in% "drug")])
mat[!is.na(mat)] <- 1
mat[is.na(mat)] <- 0
class(mat) <- "numeric"

library(gplots)
## heatmap.2(mat, scale="none", trace="none")
tmp <- ohsu.fimm.drugs$`Mechanism/Targets`
names(tmp) <- ohsu.fimm.drugs$ID_Drug.ohsu
row.annos <- tmp[rownames(mat)]
old.rows <- rownames(mat)
new.rows <- unlist(apply(cbind(old.rows, row.annos), 1, function(row) paste0(substr(row[2], start=1, stop=max(length(row[2]), 20)), ": ", row[1])))
## rownames(mat) <- row.annos
rownames(mat) <- new.rows

## heatmap.2(mat, trace="none", Colv=FALSE, scale="none", cexRow=0.5, dendrogram="row", cexCol=0.7)
## heatmap.2(mat, trace="none", Colv=FALSE, scale="none", cexRow=0.5, dendrogram="none", cexCol=0.7)

## heatmap(mat, Colv=NA, scale="none", cexRow=0.5, cexCol=0.7, Rowv=NA)

## pdf("lasso-sig-genes-ohsu.pdf")
mydoc <- addSlide( mydoc, slide.layout = "Title and Content" )
##    mydoc <- addPlot( doc = mydoc, fun = print, x = g)
mydoc <- addTitle( mydoc, 'Gene features vs drug (class)' )
mydoc <- addPlot( doc = mydoc, fun = function() { heatmap(1-mat[order(rownames(mat)),rev(colnames(mat))], Colv=NA, scale="none", cexRow=1, cexCol=1, Rowv=NA, margins = c(5,9)) })

## heatmap(1-mat[order(rownames(mat)),rev(colnames(mat))], Colv=NA, scale="none", cexRow=1, cexCol=1, Rowv=NA, margins = c(5,9))
## d <- dev.off()

writeDoc(mydoc, "062017-lasso.pptx")

cat("Done fitting elastic net to OHSU\n")

save.image(".Rdata")

library(pdist)

rdist.w.na <- function(X,Y)
{
  if (!is.matrix(X)) 
    X = as.matrix(X)
  if (!is.matrix(Y)) 
    Y = as.matrix(Y)
  distances <- matrix(pdist(X,Y)@dist, ncol=nrow(X), byrow = TRUE)
  #count NAs
  na.count <- sapply(1:nrow(X),function(i){rowSums(is.na(Y) | is.na(X[i,]))})
  #scaling to number of cols
  distances * sqrt(ncol(X)/(ncol(X) - na.count))
}

## heatmap.2(scale(t2), distfun = function(x, ...) rdist.w.na(x, x), na.color="blue")
rd <- rdist.w.na(scale(t2[rownames(t2) != "R406",]),scale(t2[rownames(t2) != "R406",]))
heatmap.2(scale(t2[rownames(t2) != "R406",]), distfun = function(x, ...) suppressWarnings(rdist.w.na(x, x)), na.color="blue", scale="none")

heatmap.2(t(scale(t(t2[rownames(t2) != "R406",]))), distfun = function(x, ...) suppressWarnings(rdist.w.na(x, x)), na.color="blue", scale="none", trace = "none", Rowv = NA, Colv = NA)
