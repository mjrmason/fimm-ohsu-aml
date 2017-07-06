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

## Read in the FIMM expression data
## Load (batch corrected) FIMM expression data
obj <- synGet(id="syn8270602", downloadFile = TRUE)
file <- getFileLocation(obj)
fimm.expr <- read.table(file, header=TRUE, sep=",")
##fimm.expr <- fread.with.rownames(file, sep=",")
fimm.expr <- t(fimm.expr)

## Load FIMM metadata (including relapse/refractory/diagnosis)
synId <- "syn8270594"
obj <- synGet(id=synId, downloadFile = TRUE)
file <- getFileLocation(obj)
fimm.metadata <- read.table(file, header=TRUE, sep=",")



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

all.ohsu.common.drugs <- na.omit(ohsu.fimm.drugs$ID_Drug.ohsu)
## FIMM drug ids
all.fimm.common.drugs <- na.omit(ohsu.fimm.drugs$ID_Drug)

mek.inhibitors <- ohsu.fimm.drugs$ID_Drug.ohsu[ohsu.fimm.drugs$`Mechanism/Targets` == "MEK1/2 inhibitor"]
ohsu.drugs <- all.ohsu.common.drugs
fimm.drugs <- all.fimm.common.drugs

save.image(".Rdata")

if(FALSE) {
  fimm.prep <- prepare.fimm.drug.response.and.expr.matrices(fimm.dss.t0, fimm.expr, fimm.drugs, response.col = "auc.l4")
  target.drug <- "FIMM000143"
  target.alpha <- 0
  X <- fimm.prep$expr.df
  y <- fimm.prep$drc.df[target.drug,]
  Xy <- t(as.data.frame(rbind(X,y)))
  Xy <- na.omit(Xy)
  
}

cat("Processing FIMM results\n")
num.bootstraps <- 50
## Let's exclude negative areas--NAH
## fimm.lst <- bootstrap.elastic.net.fimm(fimm.dss.t0, fimm.expr, fimm.drugs, response.col = "auc.l4", alphas = c(0, 0.25, 0.5, 0.75, 1), num.bootstraps = num.bootstraps)
##fimm.drugs <- c("FIMM100365", "FIMM003751", "FIMM023795", "FIMM003793", "FIMM003782", "FIMM100382")
##fimm.drugs <- "FIMM100365,FIMM003751,FIMM003793,FIMM023795,FIMM003775,FIMM003774,FIMM023820,FIMM003760,FIMM003730,FIMM003782"
##fimm.drugs <- unlist(strsplit(fimm.drugs, split=","))
sub <- subset(fimm.dss.t0, (gof.l4 > 0.5) & (b.l4 < 0))
fimm.lst <- bootstrap.elastic.net.fimm(sub, fimm.expr, fimm.drugs, response.col = "auc.l4", alphas = c(0,1), num.bootstraps = num.bootstraps)
cat("Done processing FIMM results\n")

save.image(".Rdata.fimm")

## Only take negative slopes
##sub <- subset(fimm.dss.t0, (gof.l4 > 0.5) & (b.l4 < 0))
##fimm.out.lst <- bootstrap.elastic.net.fimm(sub, fimm.expr, fimm.drugs, response.col = "auc.l4", alphas = c(0), num.bootstraps = num.bootstraps, remove.outlier = FALSE)
##cat("Done processing FIMM results\n")

save.image(".Rdata.fimm")


fimm.pp <- post.process.fits(fimm.lst, "fimm-")
cat("Done post-processing FIMM\n")


##fimm.out.pp <- post.process.fits(fimm.out.lst, "fimm-out-")
##cat("Done post-processing FIMM\n")

save.image(".Rdata.fimm")

##q()

cat("Processing OHSU results\n")
## ohsu.lst <- bootstrap.elastic.net.ohsu(ohsu.dss.t0, ohsu.expr, ohsu.rnaseq.sample.summary, ohsu.drugs, response.col = "auc.l4", alphas = c(0, 0.25, 0.5, 0.75, 1), num.bootstraps = num.bootstraps)
sub <- subset(ohsu.dss.t0, (gof.l4 > 0.5) & (b.l4 < 0))
ohsu.lst <- bootstrap.elastic.net.ohsu(sub, ohsu.expr, ohsu.rnaseq.sample.summary, ohsu.drugs, response.col = "auc.l4", alphas = c(0,1), num.bootstraps = num.bootstraps)
cat("Done processing OHSU results\n")

save.image(".Rdata.ohsu")

post.process.fits(ohsu.lst, "ohsu-")
cat("Done post-processing FIMM\n")

save.image(".Rdata.ohsu")
cat("Exiting successfully\n")
q(status = 0)

## tmp.drugs <- c("FIMM136387", "FIMM136402")
## lst <- bootstrap.elastic.net.fimm(fimm.dss.t0, fimm.expr, tmp.drugs, response.col = "auc.l4", alphas = c(0, 0.25, 0.5, 0.75, 1), num.bootstraps = 2)
## mek.drugs <- c("FIMM003708", "FIMM003751", "FIMM133840")
## lst <- bootstrap.elastic.net.fimm(fimm.dss.t0, fimm.expr, mek.drugs, response.col = "auc.l4", alphas = c(0, 1), num.bootstraps = 20)


cat("Done with bootstrap\n")


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
