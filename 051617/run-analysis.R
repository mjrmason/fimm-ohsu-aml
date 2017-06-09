## Load in the OHSU fits
load(".RData.ohsu.fits")

bioc.packages <- c("nplr", "openxlsx", "drc", "minpack.lm", "maftools", "vcd", "binom", "DMwR", "WGCNA")

install.bioc.packages.auto <- function(x) { 
  x <- as.character(x)
  print(x)
##  if(isTRUE(x %in% .packages(all.available=TRUE))) { 
##    eval(parse(text = sprintf("require(\"%s\")", x)))
##  } else { 
##    #update.packages(ask= FALSE) #update installed packages.
##    eval(parse(text = sprintf("install.packages(\"%s\", dependencies = TRUE)", x)))
##  }

  if(isTRUE(x %in% .packages(all.available=TRUE))) { 
    eval(parse(text = sprintf("require(\"%s\")", x)))
  } else {
    source("http://bioconductor.org/biocLite.R")
    #biocLite(character(), ask=FALSE) #update installed packages.
    eval(parse(text = sprintf("biocLite(\"%s\")", x)))
    eval(parse(text = sprintf("require(\"%s\")", x)))
  }
}

d <- unlist(lapply(bioc.packages, install.bioc.packages.auto))

suppressPackageStartupMessages(library("nplr"))
suppressPackageStartupMessages(library("foreach"))
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(tidyr))

suppressPackageStartupMessages(library("drc"))
suppressPackageStartupMessages(library("minpack.lm"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("gridExtra"))
suppressPackageStartupMessages(library("caTools"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("maftools"))
suppressPackageStartupMessages(library("pcaMethods"))

library(ggbeeswarm)
library(plyr)
library(stringr)
library(synapseClient)
synapseLogin()

output.path <- "output"
if (!file.exists(output.path)) {
  dir.create(output.path)
}

my.dup <- function(x) {
  duplicated(x, fromLast=TRUE) | duplicated(x, fromLast=FALSE)
}

suppressPackageStartupMessages(library("scales")) ## for scientific

suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(Biobase))

Sys.setenv(R_ZIPCMD= "/usr/bin/zip")  

num.cores <- detectCores()
if(!is.na(num.cores) && (num.cores > 1)) {
  suppressPackageStartupMessages(library("doMC"))
  cat(paste("Registering ", num.cores-1, " cores.\n", sep=""))
  registerDoMC(cores=(num.cores-1))
}

## Filter fits
## min.gof:  exclude fits having a gof < min.gof
## remove.non.finite.ic50:  remove ic50 values that are NA or NaN.
## require.ic50.be.in.range:  require that ic50 be between min and max concentration tested
filter.drug.response.fits <- function(fits.df, min.gof = 0, remove.non.finite.ic50 = TRUE, require.ic50.be.in.range = TRUE) {
  if(remove.non.finite.ic50) {
    old.nrow <- nrow(fits.df)
    flag <- is.finite(fits.df$IC50)
    fits.df <- fits.df[flag,]
    new.nrow <- nrow(fits.df)
    cat(paste0("Filtering by finite IC50 reduces ", old.nrow, " fits to ", new.nrow, "\n"))
  }

  if(require.ic50.be.in.range) {
    old.nrow <- nrow(fits.df)
    flag <- is.finite(fits.df$IC50) & is.finite(fits.df$Max.Conc.tested) & is.finite(fits.df$Min.Conc.tested) & ((10^fits.df$IC50) >= fits.df$Min.Conc.tested) & ((10^fits.df$IC50) <= fits.df$Max.Conc.tested)
    fits.df <- fits.df[flag,]
    new.nrow <- nrow(fits.df)
    cat(paste0("Filtering by requiring IC50 be in range reduces ", old.nrow, " fits to ", new.nrow, "\n"))
  }
  
  if(min.gof > 0) {
    old.nrow <- nrow(fits.df)
    flag <- is.finite(fits.df$gof) & (fits.df$gof >= min.gof)
    fits.df <- fits.df[flag,]
    new.nrow <- nrow(fits.df)
    cat(paste0("Filtering by GOF < ", min.gof, " reduces ", old.nrow, " fits to ", new.nrow, "\n"))
  }
  
  fits.df
}

ohsu.fits.gof.0.7.df <- filter.drug.response.fits(ohsu.fits.df, min.gof = 0.7, remove.non.finite.ic50 = TRUE, require.ic50.be.in.range = TRUE)

## Get a single IC50 value for each patient.  Do so by taking the median within a patient.
tmp <- ddply(ohsu.fits.gof.0.7.df, c("inhibitor", "patient_id"), 
             .parallel = TRUE,
             .fun = function(df) {
               median(df$IC50)
             })
colnames(tmp) <- c("inhibitor", "patient_id", "IC50")
ohsu.ic50s <- spread(tmp, key = patient_id, value = IC50)
rownames(ohsu.ic50s) <- ohsu.ic50s$inhibitor

## Restrict to those for which we have RNA-seq
ohsu.ic50s <- ohsu.ic50s[, colnames(ohsu.ic50s) %in% rnaseq.sample.summary.tbl$PatientID]

mat <- t(ohsu.ic50s)
## When we scale, we will compute sd--so ensure we have at least 2 values in each column
flag <- unlist(apply(mat, 2, function(col) length(which(!is.na(col))) > 1))
mat <- mat[, flag]

smat <- scale(mat)

## 75% of the data are missing
## length(which(is.na(as.vector(mat))))/length(mat)

## Plot the fraction of missing data according to drug
monotherapy <- !grepl(colnames(mat), pattern=" - ")

ms <- unlist(apply(mat[,monotherapy], 2, function(col) length(which(is.na(col))) / length(col)))
## y axis is fraction less than or equal to x
pdf("ecdf-drug.pdf")
plot(ecdf(ms), xlab = "Fraction Missing", main = "CDF of Fraction Missing\nwithin each Drug (Monotherapies Only)")
d <- dev.off()

ms <- unlist(apply(mat[,monotherapy], 1, function(row) length(which(is.na(row))) / length(row)))
pdf("ecdf-patient.pdf")
plot(ecdf(ms), xlab = "Fraction Missing", main = "CDF of Fraction Missing\nwithin each Patient (Monotherapies Only)")
d <- dev.off()

stop("stop")


## Determine the potential of using different (non-linear) PCA approaches, including
## those that tolerate missing data.

## Two general approaches to doing these comparisons:
## (1) Simulate data, leave some out, impute/PCA and compare to complete data
## (2) Impute from read data to get complete data, leave some out, impute/PCA and compare
##     to complete/imputed data.

## Do (2) first.

## We'll impute using kNN.  Choose the k that gives us the best match to the IQR of the
## complete data.
iqr.complete <- IQR(as.vector(mat), na.rm=TRUE)
## Note that variables/"genes"/drugs are columns and samples are rows
l <- llply(1:10, .parallel = TRUE,
      .fun = function(k) {
        imat <- llsImpute(mat, center=FALSE, k=k)
        imat
      })

k.opt <- which.min(unlist(lapply(l, function(imat) ks.test(as.vector(completeObs(imat)), na.omit(as.vector(mat)))$statistic)))

k <- k.opt
k <- 4
co <- completeObs(l[[unname(k)]])
plot(density(co[abs(as.vector(co)) < 3]))
lines(density(na.omit(as.vector(mat))))

## I'm going to go with k = 4.  This isn't the best by KS test, but is much better than 1:3 by ks,
## looks qualitatively similar to the best k.opt = 8, and seems more reasonable than 8 (e.g., 8
## drugs in a drug class seems large)

k.opt <- 4
imat <- llsImpute(mat, center=FALSE, k=k)
## (2) Impute from real data to get complete data, leave some out, impute/PCA and compare
##     to complete/imputed data.

simat <- scale(completeObs(imat))

methods <- c("nipals", "bpca", "ppca", "svdImpute", "nlpca")
methods <- c("nipals", "bpca", "ppca", "svdImpute")

bootstrap.pca <- function(mat, pca.method, nboots = 100, num.eigenvalues = 10) {
  llply(1:nboots, .parallel = TRUE,
        .fun = function(ind) {
          rows <- sample.int(nrow(mat), nrow(mat), replace=TRUE)
          smat <- scale(mat[rows,])
          pca(smat, method=pca.method, center=FALSE, nPcs=num.eigenvalues)
        })
}

pcas <- list()
for(pca.method in c("nipals", "bpca", "ppca", "svdImpute")) {
  print(pca.method)
  pcas[[pca.method]] <- bootstrap.pca(na.mat, pca.method, nboots=200)
}

nipals.pcas <- bootstrap.pca(na.mat, "nipals", nboots=20, num.eigenvalues=ncol(na.mat))

nipals.sdevs <- lapply(pcas[["nipals"]], function(x) sDev(x))
nipals.sdevs <- do.call("rbind", nipals.sdevs)

df <- gather(as.data.frame(nipals.sdevs), PC, sDev)
df$PC <- factor(df$PC, levels = unlist(lapply(1:10, function(x) paste0("PC", x))))
df$var <- df$sDev^2

g <- ggplot()
g <- g + geom_beeswarm(data = df, aes(x = PC, y = var))
print(g)

g <- ggplot()
g <- g + geom_boxplot(data = df, aes(x = PC, y = var))
g <- g + xlab("PC")
g <- g + ylab("Variance of PC")
g <- g + ggtitle("Sensitivity of NIPALS PCA (200 Bootstraps)")
print(g)

nipals.vars <- lapply(nipals.pcas, function(x) (sDev(x)[1:10]^2)/sum(sDev(x)^2))
nipals.vars <- do.call("rbind", nipals.vars)

df <- gather(as.data.frame(nipals.vars), PC, var)
df$PC <- factor(df$PC, levels = unlist(lapply(1:ncol(na.mat), function(x) paste0("PC", x))))

pdf("nipals-sensitivity.pdf")
g <- ggplot()
g <- g + geom_boxplot(data = df, aes(x = PC, y = var))
g <- g + xlab("PC")
g <- g + ylab("Fraction of Variance Explained")
g <- g + ggtitle("Sensitivity of NIPALS PCA (20 Bootstraps)")
print(g)
d <- dev.off()

nipals.pca <- pca(na.mat, method="nipals", center=TRUE, nPcs=ncol(na.mat))

loadings <- nipals.pca@loadings
loadings <- loadings[rownames(loadings) %in% ohsu.fimm.drugs$ID_Drug.ohsu,]
loading.colors <- data.frame(name = ohsu.fimm.drugs$ID_Drug.ohsu, class = ohsu.fimm.drugs$Class.explained)
rownames(loading.colors) <- loading.colors$name
plot(loadings[,1], loadings[,2], col=loading.colors[rownames(loadings),"class"])


plot.eigenvalues.with.missing.data <- function(mat, pca.methods = c("nipals", "bpca", "ppca", "svdImpute"), prefix = "prefix", num.eigenvalues = 2) {
  
  ## Here, just using the various pca methods on the complete data
  smat <- scale(mat)
  l <- llply(pca.methods, .parallel = TRUE,
             .fun = function(method) {
               if(method == "nlpca") {
                 return(pca(smat, method=method, center=FALSE, nPcs=num.eigenvalues, maxSteps=10000))
               } else {
                 return(pca(smat, method=method, center=FALSE, nPcs=num.eigenvalues))
               }
             })
  names(l) <- c(pca.methods)
  pca.list <- list()
  pca.list[["0"]] <- l
  
  ## Plot the eigenvalues
  gs <- list()
  frac <- 0
  for(method in pca.methods) {
    component <- 1
    df <- data.frame(x = 1:num.eigenvalues, y = sDev(pca.list[[as.character(frac)]][[method]]))
    g <- ggplot(df)
    g <- g + geom_point(aes(x = x, y = y))
    g <- g + xlab("Eigenvalue number")
    g <- g + ylab(method)
    gs[[method]] <- g
  }
  pdf(paste0(prefix, "-eigenvalues-", frac, "-missing.pdf"))
  do.call("grid.arrange", gs)
  d <- dev.off()
}

apply.pca.with.missing.data <- function(mat, pca.methods = c("nipals", "bpca", "ppca", "svdImpute"), prefix = "prefix") {

  ## Here, just using the various pca methods on the complete data
  smat <- scale(mat)
  l <- llply(c(pca.methods, "svd"), .parallel = TRUE,
             .fun = function(method) {
               if(method == "nlpca") {
                 return(pca(smat, method=method, center=FALSE, nPcs=5, maxSteps=10000))
               } else {
                 return(pca(smat, method=method, center=FALSE, nPcs=5))
               }
            })
  names(l) <- c(pca.methods, "svd")
  pca.list <- list()
  pca.list[["0"]] <- l

  pdf(paste0(prefix, "-", "-svd-pcs1-and-2.pdf"))
  df <- data.frame(x = pca.list[["0"]][["svd"]]@scores[,1], 
                   y = pca.list[["0"]][["svd"]]@scores[,2])
  g <- ggplot(df)
  g <- g + geom_point(aes(x = x, y = y))
  g <- g + xlab("SVD PC1")
  g <- g + ylab("SVD PC2")
  print(g)
  d <- dev.off()
   
  ## Plot SVD vs all other PCA approaches
  gs <- list()
  frac <- 0
  for(method in pca.methods) {
    component <- 1
    df <- data.frame(x = pca.list[[as.character(frac)]][["svd"]]@scores[,component], y = pca.list[[as.character(frac)]][[method]]@scores[,component])
    g <- ggplot(df)
    g <- g + geom_point(aes(x = x, y = y))
    g <- g + xlab("SVD")
    g <- g + ylab(method)
    gs[[method]] <- g
  }
  pdf(paste0(prefix, "-", frac, "-missing.pdf"))
  do.call("grid.arrange", gs)
  d <- dev.off()

  ## Plot the eigenvalues
  gs <- list()
  frac <- 0
  for(method in pca.methods) {
    component <- 1
    df <- data.frame(x = sDev(pca.list[[as.character(frac)]][["svd"]]), y = sDev(pca.list[[as.character(frac)]][[method]]))
    g <- ggplot(df)
    g <- g + geom_point(aes(x = x, y = y))
    g <- g + xlab("SVD")
    g <- g + ylab(method)
    gs[[method]] <- g
  }
  pdf(paste0(prefix, "-eigenvalues-", frac, "-missing.pdf"))
  do.call("grid.arrange", gs)
  d <- dev.off()
  
  
  methods.with.svd <- c(pca.methods, "imputeSVD")
  
  ## Leave out various fractions
  leave.out.fracs <- c(0.01, seq(from=0.05, to=0.75, by = 0.05))
  l2 <- llply(leave.out.fracs,
              .parallel = TRUE,
              .fun = function(leave.out.frac) {
                rmat.lo <- mat
                rmat.lo[sample(1:length(rmat.lo), floor(leave.out.frac*length(rmat.lo)), replace = FALSE)] <- NA
                srmat.lo <- scale(rmat.lo)
                
                l <- llply(methods.with.svd, .parallel = FALSE,
                           .fun = function(method) {
                             if(method == "nlpca") {
                               maxSteps <- 10000
                               maxSteps <- 100
                               return(pca(srmat.lo, method=method, center=FALSE, nPcs=5, maxSteps=maxSteps))
                             } else if(method == "imputeSVD") {
                               k.opt <- 4
                               imat <- llsImpute(rmat.lo, center=FALSE, k=k)
                               ## (2) Impute from real data to get complete data, leave some out, impute/PCA and compare
                               ##     to complete/imputed data.
                               
                               srmat.lo <- scale(completeObs(imat))
                               return(pca(srmat.lo, method="svd", center=FALSE, nPcs=5))
                             } else {
                               return(pca(srmat.lo, method=method, center=FALSE, nPcs=5))
                             }
                          })
                names(l) <- methods.with.svd
                l
            })
  names(l2) <- leave.out.fracs
  pca.list2 <- l2

  for(frac in leave.out.fracs) {
    gs <- list()
    for(method in methods.with.svd) {
      component <- 1
      df <- data.frame(x = pca.list[["0"]][["svd"]]@scores[,component], y = pca.list2[[as.character(frac)]][[method]]@scores[,component])
      ## df <- data.frame(x = sDev(pca.list[["0"]][["svd"]]), y = sDev(pca.list2[[as.character(frac)]][[method]]))
      g <- ggplot(df)
      g <- g + geom_point(aes(x = x, y = y))
      g <- g + xlab("SVD (complete data)")
      g <- g + ylab(paste0(method, "\n(", frac * 100, "% missing)"))
      gs[[method]] <- g
    }
    pdf(paste0(prefix, "-", frac, "-missing.pdf"))
    do.call("grid.arrange", gs)
    d <- dev.off()
  }
  
  for(frac in leave.out.fracs) {
    gs <- list()
    for(method in methods.with.svd) {
      component <- 1
      ## df <- data.frame(x = pca.list[["0"]][["svd"]]@scores[,component], y = pca.list2[[as.character(frac)]][[method]]@scores[,component])
      df <- data.frame(x = sDev(pca.list[["0"]][["svd"]]), y = sDev(pca.list2[[as.character(frac)]][[method]]))
      g <- ggplot(df)
      g <- g + geom_point(aes(x = x, y = y))
      g <- g + xlab("SVD (complete data)")
      g <- g + ylab(paste0(method, "\n(", frac * 100, "% missing)"))
      gs[[method]] <- g
    }
    pdf(paste0(prefix, "-eigenvalues-", frac, "-missing.pdf"))
    do.call("grid.arrange", gs)
    d <- dev.off()
  }
  
}


## Different approach--take a subset of the data that have fewer missing.
## First, exclude the worst 50% of the patients
ms <- unlist(apply(mat, 1, function(row) length(which(is.na(row))) / length(row)))
plot(ecdf(ms), xlab = "Fraction Missing", main = "CDF of Fraction Missing\nwithin each Patient")

mat.sub <- mat[ms < sort(ms, decreasing=TRUE)[floor(length(ms)/2)],]

## Now keep any drugs that have < 40% missing
ms <- unlist(apply(mat.sub, 2, function(col) length(which(is.na(col))) / length(col)))
mat.sub <- mat.sub[, ms < 0.4]
mat.sub.before <- mat.sub
for(col in 1:ncol(mat.sub)) {
  mu <- mean(mat.sub[,col], na.rm=TRUE)
  std <- sd(mat.sub[,col], na.rm=TRUE)
  flag <- is.na(mat.sub[,col])
  n.na <- length(which(flag))
  mat.sub[flag,col] <- rnorm(n.na, mean = mu, sd = std)
}

tmp <- mat[,monotherapy]
tmp[!is.na(tmp)] <- 1
tmp[is.na(tmp)] <- 0
pdf("monotherapy-missing.pdf")
frac.missing <- length(which(is.na(as.vector(mat[,monotherapy]))))/length(as.vector(mat[,monotherapy]))
per.missing = round(100 * frac.missing, digits=0)

heatmap.2(tmp, scale="none", trace="none", main = paste0(per.missing, "% missing data (red)"))
d <- dev.off()

ms <- unlist(apply(mat.sub.before, 1, function(row) length(which(is.na(row))) / length(row)))
plot(ecdf(ms), xlab = "Fraction Missing", main = "CDF of Fraction Missing\nwithin each Patient")

tmp <- mat.sub.before[ms < 0.2,]

ms <- unlist(apply(tmp, 1, function(row) length(which(is.na(row))) / length(row)))
plot(ecdf(ms), xlab = "Fraction Missing", main = "CDF of Fraction Missing\nwithin each Patient")

tm <- tmp
tm[!is.na(tm)] <- 1
tm[is.na(tm)] <- 0
heatmap(tm, scale="none")

res <- blockwiseModules(mat, replaceMissingAdjacencies = FALSE,
                        power = 10, corType="bicor", deepSplit = 3,
                        TOMType = "signed", minModuleSize = 100,
                        numericLabels = TRUE,
                        pamStage = TRUE, pamRespectsDendro = TRUE,
                        saveTOMs = TRUE,
                        saveTOMFileBase = "WGCNA",
                        verbose = 5,
                        maxBlockSize = 10000)

heatmap(mat, dist = function(x) {
  dst <- suppressWarnings(as.dist((1-cor(t(x), use="pairwise", method="pearson"))/2))
  dst.mat <- as.matrix(dst)
  dst.mat[!is.finite(dst.mat)] <- 1
  as.dist(dst.mat)
})


corRaw <- corNA(mat, use="pairwise")

dissimilarity <- 1 - corNA(mat, use="pairwise")
distance <- as.dist(dissimilarity)

corDist <- function(x) {
  dst <- suppressWarnings(as.dist((1-cor(t(x), use="pairwise", method="spearman"))))
  dst.mat <- as.matrix(dst)
  dst.mat[!is.finite(dst.mat)] <- 1
  as.dist(dst.mat)
}

flag <- unlist(apply(mat, 2, function(col) length(which(!is.na(col))) >= 10))
na.mat <- mat[,flag]
compound <- grepl(colnames(na.mat), pattern=" - ")
na.mat <- na.mat[,!compound]

distance <- corDist(t(na.mat))
hc <- hclust(distance)
attr(distance, "Labels") <- unlist(lapply(attr(distance, "Labels"), function(str) substr(str, 1, max(5:length(str)))))
clusters <- cutree(hc, h=0.4)
plot(hclust(distance))
pairs(na.mat[,names(clusters[clusters==3])])

res <- blockwiseModules(na.mat, replaceMissingAdjacencies = FALSE,
                        power = 10, corType="bicor", deepSplit = 3,
                        TOMType = "signed", minModuleSize = 2,
                        numericLabels = TRUE,
                        pamStage = TRUE, pamRespectsDendro = TRUE,
                        saveTOMs = TRUE,
                        saveTOMFileBase = "WGCNA",
                        verbose = 5,
                        maxBlockSize = 10000)

me <- moduleEigengenes(na.mat, colors = rainbow(ncol(na.mat)), impute=FALSE)

## Define clusters
## Plot variance explained by top PC (if > 2 drugs)
## Plot a few examples

## Plot cutree
## Plot var explained by top module PC
## Plot median R2 
## Get code for upper diagonal
## Plot a few examples

## Plot bootstrap

## Do correlation of all genes vs eigendrugs
panel.rank.cor <- function(x, y, digits=2, prefix="", cex.cor, ...) 
{
  usr <- par("usr"); on.exit(par(usr)) 
  par(usr = c(0, 1, 0, 1)) 
  ok <- is.finite(x) & is.finite(y) 
  if(length(which(ok)) < 3) { return() }
  test <- cor.test(x[ok], y[ok], method="spearman")
  r <- test$estimate 
  txt <- format(c(r, 0.123456789), digits=digits)[1] 
  txt <- paste(prefix, txt, sep="") 
  if(missing(cex.cor)) cex <- 0.8/strwidth(txt) 
  
  ##test <- cor.test(x[ok],y[ok]) 
  # borrowed from printCoefmat
  Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                   symbols = c("***", "**", "*", ".", " ")) 
  
  text(0.5, 0.5, txt, cex = cex * r) 
  text(.8, .8, Signif, cex=cex, col=2) 
}

panel.regression <- function (x, y, col = par("col"), bg = NA, pch = par("pch"), 
                              cex = 1, col.regres = "red", ...) 
{ 
  points(x, y, pch = pch, col = col, bg = bg, cex = cex) 
  ok <- is.finite(x) & is.finite(y) 
  if (length(which(ok)) > 3) 
    abline(stats::lm(y[ok] ~ x[ok]), col = col.regres, ...) 
} 

panel.rank.regression <- function (x, y, col = par("col"), bg = NA, pch = par("pch"), 
                                   cex = 1, col.regres = "red", ...) 
{ 
  ok <- is.finite(x) & is.finite(y) 
  points(rank(x[ok]), rank(y[ok]), pch = pch, col = col, bg = bg, cex = cex) 
  if (any(ok)) 
    abline(stats::lm(rank(y[ok]) ~ rank(x[ok])), col = col.regres, ...) 
} 



flag <- unlist(apply(mat, 2, function(col) length(which(!is.na(col))) >= 10))
na.mat <- mat[,flag]
compound <- grepl(colnames(na.mat), pattern=" - ")
na.mat <- na.mat[,!compound]

distance <- corDist(t(na.mat))
hc <- hclust(distance)
attr(distance, "Labels") <- unlist(lapply(attr(distance, "Labels"), function(str) substr(str, 1, max(5:length(str)))))
cutoff <- 0.7
clusters <- cutree(hc, h=cutoff)
num.clusters <- length(unique(clusters))
pdf("ic50-spearman-clustering-monotherapies.pdf")
main <- paste0("Drug Clustering\n(", num.clusters, " clusters at cutoff = ", cutoff, ")")
plot(hclust(distance), ylab = "1 - Spearman Correlation", main = main)
abline(h=cutoff)
d <- dev.off()

distance <- corDist((na.mat))
hc <- hclust(distance)
attr(distance, "Labels") <- unlist(lapply(attr(distance, "Labels"), function(str) substr(str, 1, max(5:length(str)))))
cutoff <- 1
clusters <- cutree(hc, h=cutoff)
num.clusters <- length(unique(clusters))
pdf("patient-spearman-clustering-monotherapies.pdf")
main <- paste0("Drug Clustering\n(", num.clusters, " clusters at cutoff = ", cutoff, ")")
plot(hclust(distance), ylab = "1 - Spearman Correlation", main = main)
abline(h=cutoff)
d <- dev.off()


distance <- corDist((na.mat))
clusters <- cutree(hc, h=0.7)
plot(hclust(distance), ylab = "1 - Spearman Correlation", main = "Patient Clustering")
abline(h=1)
d <- dev.off()


library(pvclust)

## Krister says that Mek inhibitors should strongly correlated in AML.
## Let's look at original, processed data
synId <- "syn7440451"
obj <- synGet(synId, downloadFile=TRUE)
ohsu.interpreted.data <- read.table(getFileLocation(obj), sep="\t", header=TRUE)
mek.annotations <- ohsu.fimm.drugs[grepl(ohsu.fimm.drugs$`Mechanism/Targets`, pattern="MEK"),]
ohsu.mek.data <- subset(ohsu.interpreted.data, drug %in% mek.annotations$ID_Drug.ohsu)
## Limit to AML
ohsu.mek.data <- subset(ohsu.mek.data, diagnosis == "ACUTE MYELOID LEUKAEMIA (AML) AND RELATED PRECURSOR NEOPLASMS")

## Look at several metrics
metrics <- c("IC50", "Area_under_the_curve")
metrics <- c("IC50")
metrics <- c("IC10", "IC25", "IC50", "IC75", "IC90", "Area_under_the_curve")
for(metric in metrics) {
  print(metric)
  ## Take the median of the metric across replicants
  tmp <- ddply(ohsu.mek.data[,c("patient_id", "drug", metric)], c("drug", "patient_id"), 
               .fun = function(df) { median(df[,metric]) })
  colnames(tmp) <- c("drug", "patient_id", metric)
  tmp <- spread_(tmp, key = "patient_id", value = metric)
  rownames(tmp) <- tmp$drug
  tmp <- tmp[,-1]
  ## Remove any outliers
  tmp <- as.matrix(tmp)
  if(grepl(metric, pattern="IC", ignore.case=TRUE)) { 
    ## Remove any IC50s that were set to the max
    for(row in 1:nrow(tmp)) {
      flag <- tmp[row,] == max(tmp[row,], na.rm=TRUE)
      tmp[row,flag] <- NA
    }
  }
  for(row in 1:nrow(tmp)) {
    flag <- unname(abs( (tmp[row,] - mean(tmp[row,], na.rm=TRUE)) / sd(tmp[row,], na.rm=TRUE) ) > 3)
    tmp[row,flag] <- NA
  }
  if(grepl(metric, pattern="IC", ignore.case=TRUE)) { 
    tmp <- log(tmp) 
    for(row in 1:nrow(tmp)) {
      tmp[row,] <- scale(tmp[row,])
    }
  }
  pdf(paste0("ohsu-mek-", metric, ".pdf"))
  pairs(t(tmp))
  d <- dev.off()
}

plot.pairs <- function(data, metrics, drug.col = "drug", prefix = NULL, log.transform.IC = TRUE) {
  for(metric in metrics) {
    print(metric)
    ## Take the median of the metric across replicants
    tmp <- ddply(data[,c("patient_id", drug.col, metric)], c(drug.col, "patient_id"), 
                 .fun = function(df) { median(df[,metric]) })
    colnames(tmp) <- c(drug.col, "patient_id", metric)
    tmp <- spread_(tmp, key = "patient_id", value = metric)
    rownames(tmp) <- tmp[,drug.col]
    tmp <- tmp[,-1]
    ## Remove any outliers
    tmp <- as.matrix(tmp)
    if(grepl(metric, pattern="IC", ignore.case=TRUE)) { 
      ## Remove any IC50s that were set to the max
      for(row in 1:nrow(tmp)) {
        flag <- tmp[row,] == max(tmp[row,], na.rm=TRUE)
        tmp[row,flag] <- NA
      }
    }
    for(row in 1:nrow(tmp)) {
      flag <- !is.na(tmp[row,]) & unname(abs( (tmp[row,] - mean(tmp[row,], na.rm=TRUE)) / sd(tmp[row,], na.rm=TRUE) ) > 3)
      tmp[row,flag] <- NA
    }
    if(grepl(metric, pattern="IC", ignore.case=TRUE) && log.transform.IC) { 
      tmp <- log(tmp) 
      for(row in 1:nrow(tmp)) {
        tmp[row,] <- scale(tmp[row,])
      }
      for(row in 1:nrow(tmp)) {
        flag <- !is.na(tmp[row,]) & ( abs(tmp[row,]) > 3 )
        tmp[row,flag] <- NA
      }
      
    }
    if(!is.null(prefix)) {
      pdf(paste0("ohsu-mek-", metric, ".pdf"))
    }
    pairs(t(tmp))
    if(!is.null(prefix)) {
      d <- dev.off()
    }
  }
}

raw.ohsu.mek.data <- subset(ohsu.fits.df, inhibitor %in% mek.annotations$ID_Drug.ohsu)
raw.filtered.ohsu.mek.data <- subset(ohsu.fits.gof.0.7.df, inhibitor %in% mek.annotations$ID_Drug.ohsu)
## Limit to AML

venn(list("Raw" = unique(as.character(raw.ohsu.mek.data$patient_id)), "Interpreted" = unique(as.character(ohsu.mek.data$patient_id))))

raw.ohsu.mek.data <- subset(raw.ohsu.mek.data, patient_id %in% ohsu.mek.data$patient_id)
raw.filtered.ohsu.mek.data <- subset(raw.filtered.ohsu.mek.data, patient_id %in% ohsu.mek.data$patient_id)

plot.pairs(raw.ohsu.mek.data, "IC50", drug.col = "inhibitor", prefix = NULL, log.transform.IC = FALSE)
plot.pairs(raw.filtered.ohsu.mek.data, "IC50", drug.col = "inhibitor", prefix = NULL, log.transform.IC = FALSE)

## Concentration parameters (ic50, x.min, and x.max are in real, not log, space)
## LL.4 uses llogistic to fit this function:
## f(x) = c + \frac{d-c}{1+\exp(b(\log(x)-\log(e)))}
## i.e., the parameter e is the IC50 (not log IC50)
## and x is in real (not log) space
## b = slope
## c = min.asymptote
## d = max.asymptote
## e = IC50
ll.4.func <- function(x, b, c, d, e) {
  c + ( (d-c)/(1+exp(b*(log(x)-log(e)))) )
}

l.4.func <- function(x, b, c, d, e) {
  c + ( (d-c)/(1+exp(b*(x-e))) )
}


plot.LL.4.curve <- function(slope, min.asymptote, max.asymptote, ic50, x.min, x.max) {
  x.seq <- seq(from = x.min, to = x.max, by = 0.01)
  b <- slope
  c <- min.asymptote
  d <- max.asymptote
  e <- ic50
  y <- unlist(lapply(x.seq, function(x) {
    c + ( (d-c)/(1+exp(b*(log(x)-log(e)))) )
  }))
  plot(log10(x.seq), y, type="l")
}

## L.4 fits the function
## f(x) = c + \frac{d-c}{(1+\exp(b(x - e)))}
## i.e., x is on log, not real, scale
plot.L.4.curve <- function(slope, min.asymptote, max.asymptote, ic50, x.min, x.max, add = FALSE) {
  x.seq <- seq(from = x.min, to = x.max, by = 0.01)
  b <- slope
  c <- min.asymptote
  d <- max.asymptote
  e <- ic50
  y <- unlist(lapply(x.seq, function(x) {
    c + ( (d-c)/(1+exp(b*(x-e))) )
  }))
  if(add==FALSE) {
    plot(log10(x.seq), y, type="l")
  } else {
    lines(log10(x.seq), y, type="l")
  }
}

h <- head(ohsu.raw.drug.data,6)
drc.fit <- drm(normalized_viability ~ well_concentration, data = h, fct=LL.4(fixed = c(NA, NA, NA,NA), names = c("SLOPE","MIN","MAX","IC50")))
drc.fit2 <- drm(normalized_viability ~ well_concentration, data = h, fct=L.4(fixed = c(NA, NA, NA,NA), names = c("SLOPE","MIN","MAX","IC50")))

plot.LL.4.curve(coef(drc.fit)[1], coef(drc.fit)[2], coef(drc.fit)[3], coef(drc.fit)[4], x.min=0.02, x.max=10)
plot.L.4.curve(coef(drc.fit2)[1], coef(drc.fit2)[2], coef(drc.fit2)[3], coef(drc.fit2)[4], x.min=0.02, x.max=10, add=TRUE)
points(log10(h$well_concentration), h$normalized_viability)

## Get all of the mek data and fit with either log-logistic (LL.4) or logistic (L.4)
mek.ohsu.raw.drug.data <- subset(ohsu.raw.drug.data, inhibitor %in% mek.annotations$ID_Drug.ohsu)
mek.ohsu.raw.drug.data <- subset(mek.ohsu.raw.drug.data, diagnosis == "ACUTE MYELOID LEUKAEMIA (AML) AND RELATED PRECURSOR NEOPLASMS")

fit.fct.to.ohsu <- function(data, fct = "LL.4") {
  dlply(data[,c("inhibitor", "well_concentration", "normalized_viability", "time_of_read", "lab_id", "replicant", "patient_id")],
        c("patient_id", "inhibitor", "lab_id", "replicant", "time_of_read"),
        .parallel = TRUE,
        .fun = function(df) {
          fit <- NULL
          if(fct == "LL.4") {
            tryCatch({
              fit <- suppressWarnings(drm(normalized_viability ~ well_concentration, data = df, fct=LL.4(fixed = c(NA, NA, NA,NA), names = c("SLOPE","MIN","MAX","IC50"))))
            }, 
            error = function(e) {
              NULL
            })
          } else {
            tryCatch({
              fit <- suppressWarnings(drm(normalized_viability ~ well_concentration, data = df, fct=L.4(fixed = c(NA, NA, NA,NA), names = c("SLOPE","MIN","MAX","IC50"))))
            }, 
            error = function(e) {
              NULL
            })
          }
          fit
        })
  }

ll4.fits <- fit.fct.to.ohsu(mek.ohsu.raw.drug.data, fct = "LL.4")
l4.fits <- fit.fct.to.ohsu(mek.ohsu.raw.drug.data, fct = "L.4")

## All concentrations/IC50s are in real space
compute.ll.4.auc <- function(slope, min.asymptote, max.asymptote, ic50, min.conc, max.conc) {
  ## FIMM defines logistic (not log-logistic) function as:
  ## f(x) = d + ( a - d ) * [ 1 + 10^(b * c - b * x) ]^-1
  ## with (indefinite) integral
  ## Y(x) = a * x + { ( 1 / b ) * (a - d) * log10 [ 1 + 10^(b * c - b * x) ] }
  ## 
  ## LL.4 in DRC is:
  ## f(x) = c + \frac{d-c}{1+\exp(b(\log(x)-\log(e)))}
  ## f(x) = c + ( d - c ) / ( 1 + exp(b * log(x) - b * log(f)) )
  b <- slope
  c <- min.asymptote
  d <- max.asymptote
  e <- ic50
  fnc <- function(x) { ll.4.func(x, b, c, d, e) }

  val <- integrate(fnc, min.conc, max.conc)
  auc <- 0
  if(val$message == "OK") { auc <- val$value }
  auc
}

compute.l.4.auc <- function(slope, min.asymptote, max.asymptote, ic50, min.conc, max.conc, numerical = TRUE) {
  ## FIMM defines logistic (not log-logistic) function as:
  ## f(x) = d + ( a - d ) * [ 1 + 10^(b * c - b * x) ]^-1
  ## with (indefinite) integral
  ## Y(x) = a * x + { ( 1 / b ) * (a - d) * log10 [ 1 + 10^(b * c - b * x) ] }
  ## DRC L.4 is
  ## f(x) = c' + ( d' - c' ) * [ 1 + exp(b' * x - b' * e' ) ]^-1
  ## i.e., 
  ## a = d' = max.asymptote
  ## b = - ln(10) b' = -ln(10) slope
  ## c = e' = ic50
  ## d = c' = min.asymptote
  auc <- 0
  b <- slope
  c <- min.asymptote
  d <- max.asymptote
  e <- ic50
  if(numerical) {  
    fnc <- function(x) { l.4.func(x, b, c, d, e) }
    val <- integrate(fnc, min.conc, max.conc)
    if(val$message == "OK") { auc <- val$value }
  } else {
    a <-  max.asymptote
    b <- - log(10) * slope
    c <- ic50
    d <- min.asymptote
    y.int <- function(x, a, b, c, d) {
      ## This exponential sometimes blows up
      res <- 0
      exponent <- (b * c - b * x)
      pwr <- 10^exponent
      ## This exponential sometimes blows up
      ## res <- a * x + ( ( 1 / b ) * (a - d) * log10(1 + 10^exponent ) )
      ## If exponent >> 1, then approximate log10(1 + 10^exponent) ~ log10(10^exponent) = exponent
      if((pwr > 0) && is.infinite(pwr)) { 
        res <- a * x + ( ( 1 / b ) * (a - d) * exponent )
      } else {
        res <- a * x + ( ( 1 / b ) * (a - d) * log10(1 + pwr ) )
      }
      ## res <- a * x + ( ( 1 / b ) * (a - d) * log10(1 + 10^(b * c - b * x) ) )
      res
    }
    auc <- y.int(max.conc, a, b, c, d) - y.int(min.conc, a, b, c, d)
  }
  auc
}


ll4.fits.df <- ldply(ll4.fits,
                     .parallel = FALSE,
                     .fun = function(fit) {
                       gof <- 0
                       ic50 <- 0
                       auc <- 0
                       if(!is.null(fit) && (fit$fit$convergence || fit$convergence)) {
                         ssres <- sum(residuals(fit)^2)
                         resp <- na.omit(fit$dataList$resp)
                         sstot <- sum((mean(resp)-resp)^2)
                         gof <- 1 - (ssres/sstot)
                         slope <- coef(fit)[1]
                         min.asymptote <- coef(fit)[2]
                         max.asymptote <- coef(fit)[3]
                         ic50 <- coef(fit)[4]
                         print(c(slope, min.asymptote, max.asymptote, ic50, min(fit$origData$well_concentration), max(fit$origData$well_concentration)))
                         print(fit$origData)
                         auc <- compute.ll.4.auc(slope, min.asymptote, max.asymptote, ic50, min(fit$origData$well_concentration), max(fit$origData$well_concentration))
                       }
                       c(ic50, gof, auc)
                     })
colnames(ll4.fits.df)[6:8] <- c("ic50", "gof", "auc")

l4.fits.df <- ldply(l4.fits,
                     .parallel = TRUE,
                     .fun = function(fit) {
                       gof <- 0
                       ic50 <- 0
                       auc <- 0
                       slope <- 0
                       min.asymptote <- 0
                       max.asymptote <- 0
                       if(!is.null(fit) && (fit$fit$convergence || fit$convergence)) {
                         ssres <- sum(residuals(fit)^2)
                         resp <- na.omit(fit$dataList$resp)
                         sstot <- sum((mean(resp)-resp)^2)
                         gof <- 1 - (ssres/sstot)
                         slope <- coef(fit)[1]
                         min.asymptote <- coef(fit)[2]
                         max.asymptote <- coef(fit)[3]
                         ic50 <- coef(fit)[4]
                         auc <- compute.l.4.auc(slope, min.asymptote, max.asymptote, ic50, min(fit$origData$well_concentration), max(fit$origData$well_concentration))
                       }
                       c(ic50, slope, min.asymptote, max.asymptote, gof, auc)
                     })
colnames(l4.fits.df)[6:11] <- c("ic50", "slope", "min.asymptote", "max.asymptote", "gof", "auc")

compute.l.4.auc(40.6925418, 44.10752, 1543.92624, -0.06423722, 0.02, 10, numerical=FALSE)
compute.l.4.auc(40.6925418, 44.10752, 1543.92624, -0.06423722, 0.02, 10, numerical=TRUE)


l4.anal.fits.df <- ldply(l4.fits,
                    .parallel = TRUE,
                    .fun = function(fit) {
                      gof <- 0
                      ic50 <- 0
                      auc <- 0
                      if(!is.null(fit) && (fit$fit$convergence || fit$convergence)) {
                        ssres <- sum(residuals(fit)^2)
                        resp <- na.omit(fit$dataList$resp)
                        sstot <- sum((mean(resp)-resp)^2)
                        gof <- 1 - (ssres/sstot)
                        slope <- coef(fit)[1]
                        min.asymptote <- coef(fit)[2]
                        max.asymptote <- coef(fit)[3]
                        ic50 <- coef(fit)[4]
                        auc <- compute.l.4.auc(slope, min.asymptote, max.asymptote, ic50, min(fit$origData$well_concentration), max(fit$origData$well_concentration), numerical = FALSE)
                      }
                      c(ic50, gof, auc)
                    })
colnames(l4.anal.fits.df)[6:8] <- c("ic50", "gof", "auc")

l4s <- merge(l4.fits.df, l4.anal.fits.df, by = c("patient_id", "inhibitor", "lab_id", "replicant", "time_of_read"), suffixes = c(".numerical", ".analytical"))
plot(l4s$auc.numerical, l4s$auc.analytical)

plot.pairs(l4.fits.df[l4.fits.df$gof > 0.9,], "ic50", drug.col = "inhibitor", prefix = NULL, log.transform.IC = FALSE)
plot.pairs(l4.fits.df[l4.fits.df$gof > 0.9,], "auc", drug.col = "inhibitor", prefix = NULL, log.transform.IC = FALSE)

pdf("mek-inhibitor-ic50s.pdf")
plot.pairs(ll4.fits.df[ll4.fits.df$gof > 0.9,], "ic50", drug.col = "inhibitor", prefix = NULL, log.transform.IC = TRUE)
d <- dev.off()

pdf("mek-inhibitor-aucs.pdf")
plot.pairs(ll4.fits.df[ll4.fits.df$gof > 0.9,], "auc", drug.col = "inhibitor", prefix = NULL, log.transform.IC = TRUE)
d <- dev.off()

plot.pairs(ll4.fits.df[ll4.fits.df$gof > 0,], "ic50", drug.col = "inhibitor", prefix = NULL, log.transform.IC = TRUE)
plot.pairs(ll4.fits.df[ll4.fits.df$gof > 0,], "auc", drug.col = "inhibitor", prefix = NULL, log.transform.IC = TRUE)

## Let's correlate our IC50s with those from the interpreted data
m <- merge(ll4.fits.df, ohsu.mek.data, by.y = c("patient_id", "lab_id", "drug", "replicant"), by.x = c("patient_id", "lab_id", "inhibitor", "replicant"), suffixes = c(".raw", ".interp"))

m <- merge(l4.fits.df, ll4.fits.df, by = c("patient_id", "time_of_read", "inhibitor", "replicant", "lab_id"), suffixes = c(".l4", ".ll4"))
flag <- (abs(m$ic50.l4) < 1) & (abs(m$ic50.ll4) < 1)
plot(m$ic50.l4[flag], m$ic50.ll4[flag])

plot(m$auc.l4, m$auc.ll4)


## TODO:
## - look at IC50 of mek inhibitors before filtering
## - plot raw data vs curve
## - what does drc fit
## - calculate AUC of mek inhibitors
## - look at pairs
## - data dropped in mek example
## - origin of missingness
## - mike's plot (scatter of those with vs those without shared targets)

h <- head(ohsu.raw.drug.data,7)


spearman <- function(x, ...) {
  x <- as.matrix(x)
  res <- as.dist(1 - cor(x, method = "spearman", use = "pairwise"))
  res <- as.matrix(res)
  res[!is.finite(res)] <- 2
  res <- as.dist(res)
  attr(res, "method") <- "spearman"
  return(res)
}

result <- pvclust(t(na.mat.rename), method.dist=spearman, nboot=100)

na.mat.rename <- na.mat
colnames(na.mat.rename) <- unlist(lapply(colnames(na.mat.rename), function(str) substr(str, 1, max(5, length(str)))))
pv <- pvclust(t(na.mat.rename), method.dist = corDist, r = 1, nboot = 1000, parallel=TRUE)
pv.clusters <- pvpick(pv)

pv <- pvclust(na.mat, method.dist = spearman, nboot = 1000, parallel=TRUE)
cls <- pvpick(pv)
gc()

for(cluster in 1:length(cls$clusters)) {
  m <- na.omit(na.mat[,cls$clusters[[cluster]]])
  if(nrow(m) < 10) { next }
  if(ncol(m) < 2) { next }
  pdf(paste0("pvclust-", cluster, "-pairs.pdf"))
  pairs(m, lower.panel=panel.regression, upper.panel=panel.rank.cor, main = paste0("cluster ", cluster))
  d <- dev.off()
}


df <- as.data.frame(table(clusters))

var.list <- list()
spearman.list <- list()
for(cluster in unique(df$clusters[df$Freq >= 2])) {
  m <- na.omit(na.mat[,names(clusters[clusters==cluster])])
  if(nrow(m) < 10) { next }
  if(ncol(m) < 2) { next }
  print(cluster)
  pc <- pca(scale(m), method="svd", nPcs=2)
  slplot(pc)

  res <- svd(scale(m))
  var <- (res$d^2)/sum(res$d^2)
  var.list[[cluster]] <- var[1]
  print(var)
  pdf(paste0("cluster-", cluster, "-pairs.pdf"))
  pairs(m, lower.panel=panel.regression, upper.panel=panel.rank.cor, main = paste0("cluster ", cluster))
  d <- dev.off()
  cr.m <- cor(m, method="spearman")
  spearman.list[[cluster]] <- median(cr.m[lower.tri(cr.m)])
}

plt.df <- data.frame(cluster = names(spearman.list), cor = unlist(spearman.list), percent.var = unlist(var.list))
plt.df <- plt.df[order(plt.df$percent.var),]
plt.df$cluster <- factor(plt.df$cluster, levels = plt.df$cluster)

g1 <- ggplot(plt.df, aes(x = cluster, y = cor))
g1 <- g1 + geom_point()
g1 <- g1 + ylab("(Median) Spearman\nCorrelation")

g2 <- ggplot(plt.df, aes(x = cluster, y = percent.var))
g2 <- g2 + geom_point()
g2 <- g2 + ylab("Fraction Var Explained\nby Cluster Metadrug")

pdf("monotherapy-metadrug-var-explained.pdf")
grid.arrange(g2, g1)
d <- dev.off()


## Bootstrap one of the methods
## Do we get eigenvalues?
## Restrict to mono-therapies with at least 10 non-NAs.
## WGCNA
## cluster
## Justin's schedule
## Define eigenvalue from subset


## Plot pairwise correlations within cluster

                                          
                                          library(spatstat) # "im" function 
plot(im(corRaw[nrow(corRaw):1,]), main="Correlation Matrix Map")

## Now, drop a few that have many NA

tmp <- tmp[, !(colnames(tmp) %in% c("Erlotinib", "YM-155", "Ruxolitinib (INCB018424)", "Afatinib (BIBW-2992)", "PD173955", "GW-2580", "SNS-032 (BMS-387032)", "Dasatinib", "Dovitinib (CHIR-258)"))]
tm <- tmp
tm[!is.na(tm)] <- 1
tm[is.na(tm)] <- 0
heatmap(tm, scale="none")

sum(as.vector(tm))
length(tm)

plot.eigenvalues.with.missing.data(tmp, prefix="emp-restricted-NA", num.eigenvalues = 10)
k.opt <- 4
tmp.filled.in <- llsImpute(tmp, center=FALSE, k=k)

apply.pca.with.missing.data(completeObs(tmp.filled.in), prefix = "restricted-filled-in")
  
## Create a problem with structure and see if it is retain
## Make have the patients sensitive to one drug and have sensitive to another drug
mat.structured <- mat
half.patients <- sample(1:nrow(mat.structured), size = floor(nrow(mat.structured)/2), replace = FALSE)
sensitive.patient.flag <- (1:nrow(mat.structured)) %in% half.patients
num.sensitive.patients <- length(which(sensitive.patient.flag))
num.insensitive.patients <- length(which(!sensitive.patient.flag))

half.drugs <- sample(1:ncol(mat.structured), size = floor(ncol(mat.structured)/2), replace = FALSE)
sensitive.drug.flag <- (1:ncol(mat.structured)) %in% half.drugs
num.sensitive.drugs <- length(which(sensitive.drug.flag))
num.insensitive.drugs <- length(which(!sensitive.drug.flag))

sensitive.drug.mean <- -0.25
insensitive.drug.mean <- 0.25
drug.sd <- 1
mat.structured[sensitive.patient.flag, sensitive.drug.flag] <- rnorm(mean = sensitive.drug.mean, sd = drug.sd, n = num.sensitive.patients * num.sensitive.drugs)
mat.structured[!sensitive.patient.flag, sensitive.drug.flag] <- rnorm(mean = insensitive.drug.mean, sd = drug.sd, n = num.insensitive.patients * num.sensitive.drugs)
mat.structured[sensitive.patient.flag, !sensitive.drug.flag] <- rnorm(mean = insensitive.drug.mean, sd = drug.sd, n = num.sensitive.patients * num.insensitive.drugs)
mat.structured[!sensitive.patient.flag, !sensitive.drug.flag] <- rnorm(mean = sensitive.drug.mean, sd = drug.sd, n = num.insensitive.patients * num.insensitive.drugs)

heatmap(mat.structured, scale="none")

apply.pca.with.missing.data(mat.structured, prefix = "checkerboard")


## Empirically sample from the IC50 distribution for each drug
mat.sub.emp.dist <- mat.sub.before
for(col in 1:ncol(mat.sub.emp.dist)) {
  mu <- mean(mat.sub[,col], na.rm=TRUE)
  std <- sd(mat.sub[,col], na.rm=TRUE)
  flag <- is.na(mat.sub[,col])
  n.na <- length(which(flag))
  mat.sub.emp.dist[flag,col] <- rnorm(n.na, mean = mu, sd = std)
}


plot.eigenvalues.with.missing.data(mat.sub.before, prefix="emp-NA", num.eigenvalues = 20)

apply.pca.with.missing.data(mat.sub, prefix="empirical")

## Scale the matrix
smat.sub <- scale(mat.sub)

## Now, for each drug sample the NA from a distribution

## y axis is fraction less than or equal to x
plot(ecdf(ms), xlab = "Fraction Missing", main = "CDF of Fraction Missing\nwithin each Drug")


## Repeat above using simulated data
methods <- c("nipals", "bpca", "ppca", "svdImpute")
sim.pca.list <- list()

sim.data <- matrix(data = 0, nrow=nrow(mat), ncol=ncol(mat))
for(col in 1:ncol(sim.data)) {
  sim.data[,col] <- rnorm(n = nrow(sim.data), mean = mean(mat[,col], na.rm=TRUE), sd = sd(mat[,col], na.rm=TRUE))
}

apply.pca.with.missing.data(sim.data, prefix="sim-empirical")


## Here, just using the various pca methods on the complete (imputed) data
sim.data <- matrix(data = rnorm(n = prod(dim(mat))), nrow=nrow(mat))
ssim.data <- scale(sim.data)
sim.l <- llply(c(methods, "svd"), .parallel = TRUE,
           .fun = function(method) {
             if(method == "nlpca") {
               return(pca(ssim.data, method=method, center=FALSE, nPcs=5, maxSteps=10000))
             } else {
               return(pca(ssim.data, method=method, center=FALSE, nPcs=5))
             }
           })
names(sim.l) <- c(methods, "svd")
sim.pca.list <- list()
sim.pca.list[["0"]] <- sim.l

## Plot SVD vs all other PCA approaches
gs <- list()
frac <- 0
for(method in methods) {
  component <- 1
  df <- data.frame(x = sim.pca.list[[as.character(frac)]][["svd"]]@scores[,component], y = sim.pca.list[[as.character(frac)]][[method]]@scores[,component])
  g <- ggplot(df)
  g <- g + geom_point(aes(x = x, y = y))
  g <- g + xlab("SVD")
  g <- g + ylab(method)
  gs[[method]] <- g
}
do.call("grid.arrange", gs)

methods.with.svd <- c(methods, "imputeSVD")

## Leave out various fractions
leave.out.fracs <- c(0.01, seq(from=0.05, to=0.80, by = 0.05))
sim.l2 <- llply(leave.out.fracs,
            .parallel = TRUE,
            .fun = function(leave.out.frac) {
              rmat.lo <- sim.data
              rmat.lo[sample(1:length(rmat.lo), floor(leave.out.frac*length(rmat.lo)), replace = FALSE)] <- NA
              srmat.lo <- scale(rmat.lo)
              
              l <- llply(methods.with.svd, .parallel = FALSE,
                         .fun = function(method) {
                           if(method == "nlpca") {
                             maxSteps <- 10000
                             maxSteps <- 100
                             return(pca(srmat.lo, method=method, center=FALSE, nPcs=5, maxSteps=maxSteps))
                           } else if(method == "imputeSVD") {
                             k.opt <- 4
                             imat <- llsImpute(rmat.lo, center=FALSE, k=k.opt)
                             ## (2) Impute from real data to get complete data, leave some out, impute/PCA and compare
                             ##     to complete/imputed data.
                             
                             srmat.lo <- scale(completeObs(imat))
                             return(pca(srmat.lo, method="svd", center=FALSE, nPcs=5))
                           } else {
                             return(pca(srmat.lo, method=method, center=FALSE, nPcs=5))
                           }
                         })
              names(l) <- methods.with.svd
              l
            })
names(sim.l2) <- leave.out.fracs
sim.pca.list2 <- sim.l2

l_ply(leave.out.fracs,
      .parallel = TRUE,
      .fun = function(leave.out.frac) {
        gs <- list()
        frac <- leave.out.frac
        for(method in methods.with.svd) {
          component <- 1
          df <- data.frame(x = sim.pca.list[["0"]][["svd"]]@scores[,component], y = sim.pca.list2[[as.character(frac)]][[method]]@scores[,component])
          ## df <- data.frame(x = sDev(sim.pca.list[["0"]][["svd"]]), y = sDev(sim.pca.list2[[as.character(frac)]][[method]]))
          g <- ggplot(df)
          g <- g + geom_point(aes(x = x, y = y))
          g <- g + xlab("SVD (complete data)")
          g <- g + ylab(paste0(method, "\n(", frac * 100, "% missing)"))
          gs[[method]] <- g
        }
        pdf(paste0("sim-pca-", frac, "-missing.pdf"))
        do.call("grid.arrange", gs)
        d <- dev.off()
      })


lo.0.nlpca <- pca(rmat, method="nlpca", center=FALSE, nPcs=5, maxSteps=1000)
lo.0.svd <- pca(rmat, method="svd", center=FALSE, nPcs=5)
lo.0.nipals <- pca(rmat, method="nipals", center=FALSE, nPcs=5)
lo.0.bpca <- pca(rmat, method="bpca", center=FALSE, nPcs=5)

plot(lo.0.bpca@scores[,1], lo.0.svd@scores[,1])
plot(lo.0.bpca@scores[,2], lo.0.svd@scores[,2])
plot(sDev(lo.0.bpca), sDev(lo.0.svd))


plot(lo.0.nlpca@scores[,1], lo.0.svd@scores[,1])
plot(sDev(lo.0.nlpca), sDev(lo.0.svd))

plot(lo.0.nipals@scores[,1], lo.0.svd@scores[,1])
plot(sDev(lo.0.nipals), sDev(lo.0.svd))

leave.out.fracs <- c(0.01, seq(from=0.05, to=0.5, by = 0.05))
lo.nipals.res <- list()
lo.bpca.res <- list()
lo.nlpca.res <- list()

l_ply(leave.out.fracs,
      .parallel = TRUE,
      .fun = function(leave.out.frac) {
        cat(paste0("Setting ", leave.out.frac * 100, " percent of entries to NA\n"))
        rmat.lo <- rmat
        rmat.lo[sample(1:length(rmat.lo), floor(leave.out.frac*length(rmat.lo)), replace = FALSE)] <- NA
        lo.nlpca <- pca(rmat.lo, method="nlpca", center=FALSE, nPcs=5, maxSteps = 1000)
        lo.nlpca.res[[as.character(leave.out.frac)]] <- lo.nlpca
        lo.nipals <- pca(rmat.lo, method="nipals", center=FALSE, nPcs=5)
        lo.nipals.res[[as.character(leave.out.frac)]] <- lo.nipals
        lo.bpca <- pca(rmat.lo, method="bpca", center=FALSE, nPcs=5)
        lo.bpca.res[[as.character(leave.out.frac)]] <- lo.bpca
      })


## First, let's determine how much missing data we can handle.
## To do so: 
## (1) define a random matrix
## (2) randomly leave out 1%, 5%, 10%, ... 50%
## (3) compare first few PCs and eigenvalues to full matrix

## TODO:
## - match IQR to find best k for llsImpute
## - two-dimensional list [method and leave out]
## - also impute then do pca

## (1) define a random matrix
rmat <- matrix(data = rnorm(n = prod(dim(smat))), nrow=nrow(smat))
rmat <- scale(rmat)
lo.0.nlpca <- pca(rmat, method="nlpca", center=FALSE, nPcs=5, maxSteps=1000)
lo.0.svd <- pca(rmat, method="svd", center=FALSE, nPcs=5)
lo.0.nipals <- pca(rmat, method="nipals", center=FALSE, nPcs=5)
lo.0.bpca <- pca(rmat, method="bpca", center=FALSE, nPcs=5)

plot(lo.0.bpca@scores[,1], lo.0.svd@scores[,1])
plot(lo.0.bpca@scores[,2], lo.0.svd@scores[,2])
plot(sDev(lo.0.bpca), sDev(lo.0.svd))


plot(lo.0.nlpca@scores[,1], lo.0.svd@scores[,1])
plot(sDev(lo.0.nlpca), sDev(lo.0.svd))

plot(lo.0.nipals@scores[,1], lo.0.svd@scores[,1])
plot(sDev(lo.0.nipals), sDev(lo.0.svd))

leave.out.fracs <- c(0.01, seq(from=0.05, to=0.5, by = 0.05))
lo.nipals.res <- list()
lo.bpca.res <- list()
lo.nlpca.res <- list()

l_ply(leave.out.fracs,
      .parallel = TRUE,
      .fun = function(leave.out.frac) {
  cat(paste0("Setting ", leave.out.frac * 100, " percent of entries to NA\n"))
  rmat.lo <- rmat
  rmat.lo[sample(1:length(rmat.lo), floor(leave.out.frac*length(rmat.lo)), replace = FALSE)] <- NA
  lo.nlpca <- pca(rmat.lo, method="nlpca", center=FALSE, nPcs=5, maxSteps = 1000)
  lo.nlpca.res[[as.character(leave.out.frac)]] <- lo.nlpca
  lo.nipals <- pca(rmat.lo, method="nipals", center=FALSE, nPcs=5)
  lo.nipals.res[[as.character(leave.out.frac)]] <- lo.nipals
  lo.bpca <- pca(rmat.lo, method="bpca", center=FALSE, nPcs=5)
  lo.bpca.res[[as.character(leave.out.frac)]] <- lo.bpca
})

frac <- 0.50
plot(lo.nipals.res[[as.character(frac)]]@scores[,1], lo.0.svd@scores[,1])
plot(sDev(lo.nipals.res[[as.character(frac)]]), sDev(lo.0.svd))

frac <- 0.1
plot(lo.bpca.res[[as.character(frac)]]@scores[,1], lo.0.svd@scores[,1])
plot(sDev(lo.bpca.res[[as.character(frac)]]), sDev(lo.0.svd))


## Impute missing values
ismat <- knnImputation(smat, scale = FALSE)





nl <- nlpca(mat, maxSteps=100)
snl <- nlpca(smat, maxSteps=100)

resNLPCA <- pca(mat, method="nlpca", center=FALSE, nPcs=ncol(mat), maxSteps=300)

plot(nl@scores[,1], nl@scores[,2])

num.nas <- unlist(apply(ohsu.ic50s, 1, function(row) length(which(is.na(row)))))
num.nas <- unlist(apply(ohsu.ic50s, 2, function(col) length(which(is.na(col)))))

## Convert to a matrix of drug vs patient.

stop("stop")

## Add ensembl identifiers to the drug annotation table
## Translate gene symbols to ensg identifiers

symbols.to.ensg.mapping <- function(symbols) {
  # Get a mapping from gene symbol to ensembl id
  ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
  bm <- getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id'), 
              filters = 'hgnc_symbol', 
              values = symbols, 
              mart = ensembl)
  names(bm) <- c("SYMBOL", "ENSG")
  bm <- bm[!(bm$ENSG %in% c("")),]
  bm
}

all.targets <- unique(unlist(lapply(ohsu.fimm.drugs$Gene.Targets, function(x) unlist(strsplit(x, split=",[ ]*")))))
all.targets <- all.targets[!is.na(all.targets)]

sym.to.ensg <- symbols.to.ensg.mapping(all.targets)

ohsu.fimm.drugs$Ensg.Targets <- unlist(lapply(ohsu.fimm.drugs$Gene.Targets, 
                                              function(target.str) {
                                                targets <- unlist(strsplit(target.str, split=",[ ]*"))
                                                if(!(any(targets %in% sym.to.ensg$SYMBOL))) { 
                                                  return(NA)
                                                }
                                                targets <- targets[targets %in% sym.to.ensg$SYMBOL]
                                                targets <- unique(sym.to.ensg$ENSG[sym.to.ensg$SYMBOL %in% targets])
                                                paste(targets, collapse=", ")
                                              }))


## Merge the drug targets (ensembl identifiers and symbols) to the drug fits
old.nrow <- nrow(ohsu.fits.df)
ohsu.fits.df <- merge(ohsu.fits.df, unique(ohsu.fimm.drugs[,c("ID_Drug", "ID_Drug.ohsu", "Gene.Targets", "Ensg.Targets")]), by.x = "inhibitor", by.y = "ID_Drug.ohsu")
if(old.nrow != nrow(ohsu.fits.df)) {
  warning("UNEXPECTED size of ohsu.fits.df changed\n")
}

old.nrow <- nrow(fimm.fits.df)
fimm.fits.df <- merge(fimm.fits.df, unique(ohsu.fimm.drugs[,c("ID_Drug", "ID_Drug.ohsu", "Gene.Targets", "Ensg.Targets")]), by.x = "DRUG_ID", by.y = "ID_Drug")
if(old.nrow != nrow(fimm.fits.df)) {
  warning("UNEXPECTED size of fimm.fits.df changed\n")
}

## Add the translation from drug response labId to sequencing ID to the ohsu data
ohsu.fits.df <- merge(ohsu.fits.df, unique(rnaseq.sample.summary.tbl[,c("SeqID", "Original_LabID")]), by.x = "lab_id", by.y = "Original_LabID")

## When the sequence ids are used as column names, they will be converted from 12-00023 -> X12.00023.  Add this field so we can use to
## index columns so named.
ohsu.fits.df$ColSeqID <- gsub(pattern="-", x=ohsu.fits.df$SeqID, replacement=".")
ohsu.fits.df$ColSeqID <- unlist(lapply(ohsu.fits.df$ColSeqID, function(x) paste0("X", x)))

save.image(".RData.fits")
q()



common.drugs.tested <- intersect(unique(ohsu.fits.df$ID_Drug), unique(fimm.fits.df$DRUG_ID))
common.drugs.tested <- intersect(common.drugs.tested, ohsu.fimm.drugs$ID_Drug[!is.na(ohsu.fimm.drugs$Gene.Targets)])
common.drugs.tested <- intersect(common.drugs.tested, ohsu.fimm.drugs$ID_Drug[!is.na(ohsu.fimm.drugs$Ensg.Targets)])

## specimen ids are columns of expr.mat
## specimen.id.col is col of drug data frame holding the specimen id
## patient.id.col is col of drug data frame holding the patient id
correlate.expression.with.ic50 <- function(drug.df, drug.id.col, drugs.to.test, genes.to.test = NULL, expr.mat, specimen.id.col, patient.id.col) {
  df <- drug.df[drug.df[,drug.id.col] %in% drugs.to.test,]
  if(is.null(genes.to.test)) {
    df <- df[!is.na(df$Ensg.Targets),]
  }
  if(is.null(genes.to.test)) {
    flag <- unlist(lapply(df$Ensg.Targets, function(x) any(unlist(strsplit(x, split=",[ ]*")) %in% rownames(expr.mat))))
    df <- df[flag,]
  }
  if(!is.null(genes.to.test)) {
    if(!(any(genes.to.test %in% rownames(expr.mat)))) {
      warning("UNEXPECTED no genes to test are in data\n")
    }
  }
  ret <- ddply(df, drug.id.col, .parallel = TRUE,
               .fun = function(df.drug) {
                 cat(paste0("Correlating gene expr with IC50 for drug ", unique(df.drug[,drug.id.col]), "\n"))
                 drug.targets <- df.drug$Ensg.Targets[1]
                 drug.targets <- unlist(strsplit(drug.targets, split=",[ ]*"))
                 if(!is.null(genes.to.test)) {
                   drug.targets <- genes.to.test
                 }
                 drug.targets <- drug.targets[drug.targets %in% rownames(expr.mat)]
                 resp <- data.frame(IC50 = df.drug$IC50, PATIENT_ID = df.drug[,patient.id.col], id = df.drug[,specimen.id.col])
                 resp <- resp[!is.infinite(resp$IC50) & !is.nan(resp$IC50) & !is.na(resp$IC50),]
                 common.samples <- intersect(resp$id, colnames(expr.mat))
                 num.resp.samples <- length(resp$id)
                 num.unique.resp.samples <- length(unique(resp$id))
                 flag <- resp$id %in% common.samples
                 resp <- resp[flag,]
                 tmp <- as.data.frame(table(resp$PATIENT_ID))
                 colnames(tmp) <- c("PATIENT_ID", "weight")
                 tmp$weight <- 1/sqrt(tmp$weight)
                 resp <- merge(resp, tmp)
                 common.samples <- as.character(resp$id)
                 ldply(drug.targets, .parallel = FALSE,
                       .fun = function(gene) {
                         if(!(gene %in% rownames(expr.mat))) {
                           warning(paste0("UNEXPECTED: ", gene, " not in expr.mat\n"))
                         }
                         if(!(all(common.samples %in% colnames(expr.mat)))) {
                           warning("UNEXPECTED columns not in expr.mat\n")
                         }
                         gene.expr <- as.vector(unlist(expr.mat[gene, common.samples]))
                         lm.fit <- lm(resp$IC50 ~ gene.expr, weights = resp$weight)
                         sum <- summary(lm.fit)
                         f <- sum$fstatistic
                         r2 <- 0
                         p <- -1
                         if(is.numeric(f)) { 
                           r2 <- sum$r.squared
                           p <- pf(f[1],f[2],f[3],lower.tail=F)
                         }
                         num.unique.samples <- length(unique(common.samples))
                         num.samples <- length(common.samples)
                         res <- c(gene, p, r2, num.unique.samples, num.samples, num.unique.resp.samples, num.resp.samples)
                         names(res) <- c("gene", "pval", "r2", "num.uniq.samples", "num.samples", "num.uniq.resp.samples", "num.resp.samples")
                         res
                       })
                 })
  ret
}

common.genes <- intersect(rownames(ohsu.expr.data), rownames(fimm.expr.data))

ohsu.expr.corr <- correlate.expression.with.ic50(drug.df = ohsu.fits.gof.0.7.df, drug.id.col = "ID_Drug", genes.to.test = common.genes, drugs.to.test = common.drugs.tested, expr.mat = ohsu.expr.data, specimen.id.col = "SeqID", patient.id.col = "patient_id")
cat("Done correlating expression with IC50 for OHSU\n")
save.image(".RData")

## common.drugs.tested <- intersect(common.drugs.tested, c("FIMM136387", "FIMM133832", "FIMM133867", "FIMM133902"))

fimm.expr.corr <- correlate.expression.with.ic50(drug.df = fimm.fits.gof.0.7.df, drug.id.col = "DRUG_ID", genes.to.test = common.genes, drugs.to.test = common.drugs.tested, expr.mat = fimm.expr.data, specimen.id.col = "SCREEN_ID", patient.id.col = "PATIENT_ID")
cat("Done correlating expression with IC50 for FIMM\n")
cat("Done with expression\n")
save.image(".RData")
cat("Saved image\n")
stop("stop")

common.genes <- intersect(rownames(beat.mafs), rownames(fimm.genomic.bin))
ohsu.mut.corr <- correlate.genomic.with.ic50(drug.df = ohsu.fits.gof.0.7.df, drug.id.col = "ID_Drug", drugs.to.test = common.drugs.tested, genes.to.test = common.genes, mutation.mat = beat.mafs, specimen.id.col = "SeqID", patient.id.col = "patient_id")
fimm.mut.corr <- correlate.genomic.with.ic50(drug.df = fimm.fits.gof.0.7.df, drug.id.col = "DRUG_ID", drugs.to.test = common.drugs.tested, genes.to.test = common.genes, mutation.mat = fimm.genomic.bin, specimen.id.col = "SCREEN_ID", patient.id.col = "PATIENT_ID")

cat("Done with genomic\n")
save.image(".RData")
cat("Saved image\n")
stop("stop")

drug = "FIMM003774"
gene = "FLT3"
pdf(paste0(drug, "-", gene, "-ohsu-fimm-mut.pdf"))
g1 <- plot.genomic.vs.ic50(drug.df = ohsu.fits.gof.0.7.df, drug = drug, gene = gene, drug.id.col = "ID_Drug", mutation.mat = beat.mafs, specimen.id.col = "SeqID", patient.id.col = "patient_id")
g1 <- g1 + ggtitle("OHSU")
g2 <- plot.genomic.vs.ic50(drug.df = fimm.fits.gof.0.7.df, drug = drug, gene = gene, drug.id.col = "DRUG_ID", mutation.mat = fimm.genomic.bin, specimen.id.col = "SCREEN_ID", patient.id.col = "PATIENT_ID")
g2 <- g2 + ggtitle("FIMM")
grid.arrange(g1, g2)
d <- dev.off()

##for(drug in unique(fimm.mut.corr$DRUG_ID)) {
##  print(drug)
##  ohsu.mut.corr <- correlate.genomic.with.ic50(drug.df = ohsu.fits.gof.0.7.df, drug.id.col = "ID_Drug", drugs.to.test = drug, genes.to.test = common.genes, mutation.mat = beat.mafs, specimen.id.col = "SeqID", patient.id.col = "patient_id")
##}

ohsu.fimm.mut.corr <- merge(ohsu.mut.corr, fimm.mut.corr, by.x = c("ID_Drug", "gene"), by.y = c("DRUG_ID", "gene"), suffixes = c(".ohsu", ".fimm"))

library(gplots)
ohsu.fimm.mut.corr$drug.gene <- unlist(apply(ohsu.fimm.mut.corr[,c("ID_Drug", "gene")], 1, function(row) paste0(row[1], "-", row[2])))
vennList = list("FIMM"=ohsu.fimm.mut.corr$drug.gene[!is.na(ohsu.fimm.mut.corr$pval.fimm) & (ohsu.fimm.mut.corr$pval.fimm < 0.05)], "OHSU"=ohsu.fimm.mut.corr$drug.gene[!is.na(ohsu.fimm.mut.corr$pval.ohsu) & (ohsu.fimm.mut.corr$pval.ohsu < 0.05)])
pdf("ohsu-fimm-sig-mut-venn.pdf")
venn(vennList)
d <- dev.off()

ohsu.fimm.mut.corr$pval.fimm <- as.numeric(ohsu.fimm.mut.corr$pval.fimm)
ohsu.fimm.mut.corr$pval.ohsu <- as.numeric(ohsu.fimm.mut.corr$pval.ohsu)

## Plot pvalue distributions of FIMM and OHSU
pdf("ohsu-fimm-mut-pval-dist.pdf")
g1 <- ggplot(ohsu.fimm.mut.corr, aes(x=pval.fimm)) + geom_histogram() + xlab("FIMM p-value")
g2 <- ggplot(ohsu.fimm.mut.corr, aes(x=pval.ohsu)) + geom_histogram() + xlab("OHSU p-value")
grid.arrange(g1, g2)
d <- dev.off()


## specimen ids are columns of mutation.mat
## gene symbols are rows of mutation.mat
correlate.genomic.with.ic50 <- function(drug.df, drug.id.col, drugs.to.test, genes.to.test = NULL, mutation.mat, specimen.id.col, patient.id.col) {
  df <- drug.df[drug.df[,drug.id.col] %in% drugs.to.test,]
  if(is.null(genes.to.test)) {
    df <- df[!is.na(df$Gene.Targets),]
  }
  if(is.null(genes.to.test)) {
    flag <- unlist(lapply(df$Gene.Targets, function(x) any(unlist(strsplit(x, split=",[ ]*")) %in% rownames(mutation.mat))))
    df <- df[flag,]
  }
  if(!is.null(genes.to.test)) {
    if(!(any(genes.to.test %in% rownames(mutation.mat)))) {
      warning("UNEXPECTED no genes to test are in data\n")
    }
  }
  ret <- ddply(df, drug.id.col, .parallel = TRUE,
               .fun = function(df.drug) {
                 drug.targets <- df.drug$Gene.Targets[1]
                 drug.targets <- unlist(strsplit(drug.targets, split=",[ ]*"))
                 if(!is.null(genes.to.test)) {
                   drug.targets <- genes.to.test
                 }
                 drug.targets <- drug.targets[drug.targets %in% rownames(mutation.mat)]
                 resp <- data.frame(IC50 = df.drug$IC50, PATIENT_ID = df.drug[,patient.id.col], id = df.drug[,specimen.id.col])
                 resp <- resp[!is.infinite(resp$IC50) & !is.nan(resp$IC50) & !is.na(resp$IC50),]
                 common.samples <- intersect(resp$id, colnames(mutation.mat))
                 num.resp.samples <- length(resp$id)
                 num.unique.resp.samples <- length(unique(resp$id))
                 if(length(common.samples) == 0) {
                   return(ldply(drug.targets, .parallel = FALSE,
                                .fun = function(gene) {
                                  p <- NA
                                  num.unique.samples <- length(unique(common.samples))
                                  num.samples <- length(common.samples)
                                  res <- c(gene, p, num.unique.samples, num.samples, num.unique.resp.samples, num.resp.samples)
                                  names(res) <- c("gene", "pval", "num.uniq.samples", "num.samples", "num.uniq.resp.samples", "num.resp.samples")
                                  res
                                }))
                 }
                 flag <- resp$id %in% common.samples
                 resp <- resp[flag,]
                 tmp <- as.data.frame(table(resp$PATIENT_ID))
                 colnames(tmp) <- c("PATIENT_ID", "weight")
                 tmp$weight <- 1/sqrt(tmp$weight)
                 resp <- merge(resp, tmp)
                 common.samples <- as.character(resp$id)
                 ldply(drug.targets, .parallel = FALSE,
                       .fun = function(gene) {
                         if(!(gene %in% rownames(mutation.mat))) {
                           warning(paste0("UNEXPECTED: ", gene, " not in mutation.mat\n"))
                         }
                         if(!(all(common.samples %in% colnames(mutation.mat)))) {
                           warning("UNEXPECTED columns not in mutation.mat\n")
                         }
                         gene.mutation <- as.vector(unlist(mutation.mat[gene, common.samples]))
                         if(!(all(gene.mutation %in% c(0,1)))) {
                           warning(paste0("WARNING: non-binary gene mutations: ", paste(gene.mutation, collapse=", "), "\n"))
                         }
                         p <- NA
                         if(all((c(0,1) %in% gene.mutation))) {
                           test.df <- data.frame(x = factor(gene.mutation), y = resp$IC50)
                           wt <- wilcox.test(y ~ x, data = test.df)
                           p <- wt$p.value
                           print(wt)
                         }
                         num.unique.samples <- length(unique(common.samples))
                         num.samples <- length(common.samples)
                         res <- c(gene, p, num.unique.samples, num.samples, num.unique.resp.samples, num.resp.samples)
                         names(res) <- c("gene", "pval", "num.uniq.samples", "num.samples", "num.uniq.resp.samples", "num.resp.samples")
                         res
                       })
               })
  ret
}

plot.genomic.vs.ic50 <- function(drug.df, drug, gene, drug.id.col, mutation.mat, specimen.id.col, patient.id.col) {
  df.drug <- drug.df[drug.df[,drug.id.col] == drug,]
  resp <- data.frame(IC50 = df.drug$IC50, PATIENT_ID = df.drug[,patient.id.col], id = df.drug[,specimen.id.col])
  resp <- resp[!is.infinite(resp$IC50) & !is.nan(resp$IC50) & !is.na(resp$IC50),]
  common.samples <- intersect(resp$id, colnames(mutation.mat))
  flag <- resp$id %in% common.samples
  resp <- resp[flag,]
  common.samples <- as.character(resp$id)
  gene.mutation <- factor(as.vector(unlist(mutation.mat[gene, common.samples])))
  df <- data.frame(mutation = gene.mutation, resp = resp$IC50)  
  ggplot(df, aes(x = mutation, y = resp)) + 
    geom_beeswarm() +
    xlab(paste0(gene, "\nmutation status")) + ylab(paste0("Log10 ", drug, " IC50"))
}



## BEGIN load OHSU dna data


## beat.mafs <- read.mafs(maf.names = ohsu.beat.maf.names, maf.tbl = ohsu.beat.mafs, mutation.types.to.total = c("Missense_Mutation", "Nonsense_Mutation"))
beat.mafs <- read.mafs(maf.names = ohsu.beat.maf.names, maf.tbl = ohsu.beat.mafs, mutation.types.to.total = c("Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Splice_Site", "Translation_Start_Site"))
## Make the samples the columns
beat.mafs <- t(beat.mafs)

## END load OHSU dna data 


## HERE


cat("Fitting FIMM curves\n")
## This works with plyr 1.7.1 and 1.8.1, but not 1.8.4
fimm.fits <- dlply(fimm.raw.dr.subset[,c("DRUG_ID", "DRUG_NAME", "CONCENTRATION", "PERCENT_INHIBITION", "SCREEN_ID", "SAMPLE_DATE", "PATIENT_ID")],
                   c("SCREEN_ID", "DRUG_ID"),
                   .parallel = TRUE,
                   .fun = 
                     function(df) {
                       print(nrow(df))
                       id <- unique(df[, !(colnames(df) %in% c("CONCENTRATION", "PERCENT_INHIBITION")),drop=F])
                       if(nrow(df) < 2) {
                         return(list(id=id, tbl=NULL))
                       }
                       print(df)
                       dose <- df$CONCENTRATION
                       inhibition <- df$PERCENT_INHIBITION
                       patient.id <- df$PATIENT_ID[1]
                       inhibitor <- df$DRUG_NAME[1]
                       file <- paste0(patient.id, "-", inhibitor)
                       tbl <- NULL
                       fit <- NULL
                       res <- NULL
                       res <- suppressWarnings(drc.calc.ic50.robust(dose, inhibition, file, patient.id, inhibitor, return.fit = FALSE, return.residuals = FALSE, apply.constraints = FALSE, plot.fit = FALSE, robust = "mean", force.min.to.zero = FALSE))
                       if(!is.null(res) & !is.null(res$tbl)) {
                         tbl <- res$tbl
                       }
                       if(!is.null(res) & !is.null(res$res)) {
                         tbl <- res$tbl
                       }
                       list(id=id, tbl=tbl)
                    })


## WORKING

cat("Fitting OHSU curves\n")
## Fit curves to OHSU data (AML only, those drugs that overlap with FIMM drugs)
ohsu.inh.tbl.subset <- subset(ohsu.inh.tbl, (ohsu.inh.tbl$inhibitor %in% ohsu.fimm.drugs$ID_Drug.ohsu))
rm(ohsu.inh.tbl); gc()

ohsu.inh.tbl.subset <- ohsu.inh.tbl.subset[order(ohsu.inh.tbl.subset$inhibitor, ohsu.inh.tbl.subset$patient_id, ohsu.inh.tbl.subset$lab_id),]
ohsu.fits <- dlply(ohsu.inh.tbl.subset[,c("inhibitor", "well_concentration", "normalized_viability", "time_of_read", "lab_id", "replicant", "patient_id")],
                   c("patient_id", "inhibitor", "lab_id", "replicant", "time_of_read"),
                   .parallel = TRUE,
                   .fun = 
                     function(df) {
                       print(nrow(df))
                       id <- unique(df[, !(colnames(df) %in% c("well_concentration", "normalized_viability")),drop=F])
                       if(nrow(df) < 2) {
                         return(list(id=id, tbl=NULL))
                       }
                       print(df)
                       dose <- df$well_concentration
                       inhibition <- 100 - df$normalized_viability
                       patient.id <- df$patient_id[1]
                       inhibitor <- df$inhibitor[1]
                       file <- paste0(patient.id, "-", inhibitor)
                       tbl <- NULL
                       fit <- NULL
                       res <- NULL
                       res <- suppressWarnings(drc.calc.ic50.robust(dose, inhibition, file, patient.id, inhibitor, return.fit = FALSE, return.residuals = FALSE, apply.constraints = FALSE, plot.fit = FALSE, robust = "mean", force.min.to.zero = FALSE))
                       if(!is.null(res) & !is.null(res$tbl)) {
                         tbl <- res$tbl
                       }
                       if(!is.null(res) & !is.null(res$res)) {
                         tbl <- res$tbl
                       }
                       list(id=id, tbl=tbl)
                     })

cat("Done\n")
save.image(".RData")
stop("stop")

## Collect the results into a single table
tmp <- lapply(fimm.fits, function(x) x$tbl)
fimm.fits.tbl <- do.call("rbind", tmp)

tmp <- lapply(fimm.fits, function(x) x$id)
fimm.fits.id <- do.call("rbind", tmp)

fimm.fits.df <- cbind(fimm.fits.tbl, fimm.fits.id)
if(!(all(rownames(fimm.fits.tbl) == rownames(fimm.fits.id)))) {
  warning("UNEXPECTED fimm.fits.tbl rownames != fimm.fits.id rownames\n")
}
if(!(all(rownames(fimm.fits.tbl) == rownames(fimm.fits.df)))) {
  warning("UNEXPECTED fimm.fits.tbl rownames != fimm.fits.df rownames\n")
}
rm(fimm.fits.tbl)
rm(fimm.fits.id)
gc()

tmp <- lapply(ohsu.fits, function(x) x$tbl)
ohsu.fits.tbl <- do.call("rbind", tmp)

tmp <- lapply(ohsu.fits, function(x) x$id)
ohsu.fits.id <- do.call("rbind", tmp)

ohsu.fits.df <- cbind(ohsu.fits.tbl, ohsu.fits.id)
if(!(all(rownames(ohsu.fits.tbl) == rownames(ohsu.fits.id)))) {
  warning("UNEXPECTED ohsu.fits.tbl rownames != ohsu.fits.id rownames\n")
}
if(!(all(rownames(ohsu.fits.tbl) == rownames(ohsu.fits.df)))) {
  warning("UNEXPECTED ohsu.fits.tbl rownames != ohsu.fits.df rownames\n")
}
rm(ohsu.fits.tbl)
rm(ohsu.fits.id)
rm(tmp)
gc()

## Translate gene symbols to ensg identifiers
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(Biobase))

symbols.to.ensg.mapping <- function(symbols) {
  # Get a mapping from gene symbol to ensembl id
  ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
  bm <- getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id'), 
              filters = 'hgnc_symbol', 
              values = symbols, 
              mart = ensembl)
  names(bm) <- c("SYMBOL", "ENSG")
  bm <- bm[!(bm$ENSG %in% c("")),]
  bm
}

all.targets <- unique(unlist(lapply(ohsu.fimm.drugs$Gene.Targets, function(x) unlist(strsplit(x, split=",[ ]*")))))
all.targets <- all.targets[!is.na(all.targets)]

sym.to.ensg <- symbols.to.ensg.mapping(all.targets)

ohsu.fimm.drugs$Ensg.Targets <- unlist(lapply(ohsu.fimm.drugs$Gene.Targets, 
                                           function(target.str) {
                                             targets <- unlist(strsplit(target.str, split=",[ ]*"))
                                             if(!(any(targets %in% sym.to.ensg$SYMBOL))) { 
                                               return(NA)
                                             }
                                             targets <- targets[targets %in% sym.to.ensg$SYMBOL]
                                             targets <- unique(sym.to.ensg$ENSG[sym.to.ensg$SYMBOL %in% targets])
                                             paste(targets, collapse=", ")
                                           }))

## Look at correlation of drug IC50 with drug target gene expression for both
## FIMM and OHSU.  Compare p-values between the two.

old.nrow <- nrow(ohsu.fits.df)
ohsu.fits.df <- merge(ohsu.fits.df, unique(ohsu.fimm.drugs[,c("ID_Drug", "ID_Drug.ohsu", "Gene.Targets", "Ensg.Targets")]), by.x = "inhibitor", by.y = "ID_Drug.ohsu")
if(old.nrow != nrow(ohsu.fits.df)) {
  warning("UNEXPECTED size of ohsu.fits.df changed\n")
}

old.nrow <- nrow(fimm.fits.df)
fimm.fits.df <- merge(fimm.fits.df, unique(ohsu.fimm.drugs[,c("ID_Drug", "ID_Drug.ohsu", "Gene.Targets", "Ensg.Targets")]), by.x = "DRUG_ID", by.y = "ID_Drug")
if(old.nrow != nrow(fimm.fits.df)) {
  warning("UNEXPECTED size of fimm.fits.df changed\n")
}

g1 <- ggplot(data = ohsu.fits.df, aes(x = gof))
g1 <- g1 + geom_density() + xlab("OHSU GOF")
g2 <- ggplot(data = fimm.fits.df, aes(x = gof))
g2 <- g2 + geom_density() + xlab("FIMM GOF")
pdf("gof-density.pdf")
grid.arrange(g1, g2)
d <- dev.off()

all.ensg.targets <- unique(unlist(lapply(ohsu.fimm.drugs$Ensg.Targets, function(x) unlist(strsplit(x, split=",[ ]*")))))
all.ensg.targets <- all.ensg.targets[!is.na(all.ensg.targets)]

common.genes <- intersect(rownames(fimm.expr), rownames(ohsu.expr))
common.genes <- intersect(common.genes, all.ensg.targets)

fimm.expr.common <- fimm.expr[common.genes,]
ohsu.expr.common <- ohsu.expr[common.genes,]

## Correlate gene target expression with target IC50 for FIMM
df <- fimm.fits.df[fimm.fits.df$DRUG_ID %in% common.drugs.tested,]
flag <- unlist(lapply(df$Ensg.Targets, function(x) any(unlist(strsplit(x, split=",[ ]*")) %in% common.genes)))
df <- df[flag,]
fimm.expr.df <- merge(df, ohsu.fimm.drugs[,c("ID_Drug", "Gene.Targets")], by.x = "DRUG_ID", by.y = "ID_Drug")
fimm.expr.corr <- ddply(fimm.expr.df, c("DRUG_ID"), .parallel = TRUE,
                        .fun = function(df.drug) {
                          drug.targets <- df.drug$Ensg.Targets[1]
                          drug.targets <- unlist(strsplit(drug.targets, split=",[ ]*"))
                          drug.targets <- drug.targets[drug.targets %in% rownames(fimm.expr.common)]
                          resp <- data.frame(IC50 = df.drug$IC50, PATIENT_ID = df.drug$PATIENT_ID, id = df.drug$SCREEN_ID)
                          rownames(resp) <- resp$id
                          resp <- resp[!is.infinite(resp$IC50) & !is.nan(resp$IC50) & !is.na(resp$IC50),]
                          common.samples <- intersect(rownames(resp), colnames(fimm.expr.common))
                          resp <- resp[common.samples,]
                          tmp <- as.data.frame(table(resp$PATIENT_ID))
                          colnames(tmp) <- c("PATIENT_ID", "weight")
                          tmp$weight <- 1/sqrt(tmp$weight)
                          resp <- merge(resp, tmp)
                          rownames(resp) <- resp$id
                          if(!(all(common.samples %in% rownames(resp)))) {
                            warning("UNEXPECTED row names")
                          }
                          resp <- resp[common.samples,]
                          ldply(drug.targets, .parallel = FALSE,
                                .fun = function(gene) {
                                  if(!(gene %in% rownames(fimm.expr.common))) {
                                    warning(paste0("UNEXPECTED: ", gene, " not in fimm.expr\n"))
                                  }
                                  if(!(all(common.samples %in% colnames(fimm.expr.common)))) {
                                    warning("UNEXPECTED columns not in fimm expr\n")
                                  }
                                  gene.expr <- fimm.expr.common[gene, common.samples]
                                  lm.fit <- lm(resp$IC50 ~ gene.expr, weights = resp$weight)
                                  sum <- summary(lm.fit)
                                  f <- sum$fstatistic
                                  r2 <- sum$r.squared
                                  p <- pf(f[1],f[2],f[3],lower.tail=F)
                                  res <- c(gene, p, r2)
                                  num.unique.samples <- length(unique(common.samples))
                                  num.samples <- length(common.samples)
                                  res <- c(gene, p, r2, num.unique.samples, num.samples)
                                  names(res) <- c("gene", "pval", "r2", "num.uniq.samples", "num.samples")
                                  res
                                })
                        })

plot.fimm.expr.vs.ic50 <- function(df, drug, ensg.gene, expr.mat, xlab, ylab) {
  df.drug <- subset(df, DRUG_ID == drug)
  resp <- data.frame(IC50 = df.drug$IC50, PATIENT_ID = df.drug$PATIENT_ID, id = df.drug$SCREEN_ID)
  rownames(resp) <- resp$id
  resp <- resp[!is.infinite(resp$IC50) & !is.nan(resp$IC50) & !is.na(resp$IC50),]
  common.samples <- intersect(rownames(resp), colnames(expr.mat))
  resp <- resp[common.samples,]
  tmp <- as.data.frame(table(resp$PATIENT_ID))
  colnames(tmp) <- c("PATIENT_ID", "weight")
  tmp$weight <- 1/sqrt(tmp$weight)
  resp <- merge(resp, tmp)
  rownames(resp) <- resp$id
  if(!(all(common.samples %in% rownames(resp)))) {
    warning("UNEXPECTED row names")
  }
  if(!(ensg.gene %in% rownames(expr.mat))) {
    warning(paste0("UNEXPECTED: ", ensg.gene, " not in expr.mat\n"))
  }
  if(!(all(common.samples %in% colnames(expr.mat)))) {
    warning("UNEXPECTED columns not in expr.mat\n")
  }
  gene.expr <- fimm.expr.common[ensg.gene, common.samples]
  df <- data.frame(expr = gene.expr, resp = resp$IC50, weights=resp$weight)
  ggplot(df, aes(x = expr, y = resp, weight=weights)) + 
    geom_point() + 
    xlab(xlab) + ylab(ylab) +
    stat_smooth(method = "lm", col = "red")  
}


## Correlate gene target expression with target IC50 for OHSU
df <- ohsu.fits.df[ohsu.fits.df$ID_Drug %in% common.drugs.tested,]
flag <- unlist(lapply(df$Ensg.Targets, function(x) any(unlist(strsplit(x, split=",[ ]*")) %in% common.genes)))
df <- df[flag,]
df <- merge(df, ohsu.fimm.drugs[,c("ID_Drug", "Gene.Targets")], by.x = "ID_Drug", by.y = "ID_Drug")
## What is this lab_id
## Need to translate lab_id to seq_id
old.nrow <- nrow(df)
new.df <- merge(df, unique(aml.rnaseq.sample.summary.tbl[,c("SeqID", "Original_LabID")]), by.x = "lab_id", by.y = "Original_LabID")
ohsu.expr.df <- new.df
ohsu.expr.corr <- ddply(ohsu.expr.df, c("ID_Drug"), .parallel = TRUE,
                        .fun = function(df.drug) {
                          drug.targets <- df.drug$Ensg.Targets[1]
                          drug.targets <- unlist(strsplit(drug.targets, split=",[ ]*"))
                          drug.targets <- drug.targets[drug.targets %in% rownames(ohsu.expr.common)]
                          ids <- df.drug$SeqID
                          ids <- gsub(pattern="-", x=ids, replacement=".")
                          ids <- unlist(lapply(ids, function(x) paste0("X", x)))
                          resp <- data.frame(IC50 = df.drug$IC50, PATIENT_ID = df.drug$patient_id, id = ids)
                          resp <- resp[!is.infinite(resp$IC50) & !is.nan(resp$IC50) & !is.na(resp$IC50),]
                          common.samples <- intersect(resp$id, colnames(ohsu.expr.common))
                          flag <- resp$id %in% common.samples
                          resp <- resp[flag,]
                          tmp <- as.data.frame(table(resp$PATIENT_ID))
                          colnames(tmp) <- c("PATIENT_ID", "weight")
                          tmp$weight <- 1/sqrt(tmp$weight)
                          resp <- merge(resp, tmp)
                          common.samples <- as.character(resp$id)
                          ldply(drug.targets, .parallel = FALSE,
                                .fun = function(gene) {
                                  if(!(gene %in% rownames(ohsu.expr.common))) {
                                    warning(paste0("UNEXPECTED: ", gene, " not in ohsu.expr\n"))
                                  }
                                  if(!(all(common.samples %in% colnames(ohsu.expr.common)))) {
                                    warning("UNEXPECTED columns not in ohsu expr\n")
                                  }
                                  gene.expr <- as.vector(unlist(ohsu.expr.common[gene, common.samples]))
                                  lm.fit <- lm(resp$IC50 ~ gene.expr, weights = resp$weight)
                                  sum <- summary(lm.fit)
                                  f <- sum$fstatistic
                                  r2 <- sum$r.squared
                                  p <- pf(f[1],f[2],f[3],lower.tail=F)
                                  num.unique.samples <- length(unique(common.samples))
                                  num.samples <- length(common.samples)
                                  res <- c(gene, p, r2, num.unique.samples, num.samples)
                                  names(res) <- c("gene", "pval", "r2", "num.uniq.samples", "num.samples")
                                  res
                                })
                        })

plot.ohsu.expr.vs.ic50 <- function(df, drug, ensg.gene, expr.mat, xlab, ylab) {
  df.drug <- subset(df, ID_Drug == drug)
  ids <- df.drug$SeqID
  ids <- gsub(pattern="-", x=ids, replacement=".")
  ids <- unlist(lapply(ids, function(x) paste0("X", x)))
  resp <- data.frame(IC50 = df.drug$IC50, PATIENT_ID = df.drug$patient_id, id = ids)
  resp <- resp[!is.infinite(resp$IC50) & !is.nan(resp$IC50) & !is.na(resp$IC50),]
  common.samples <- intersect(resp$id, colnames(expr.mat))
  flag <- resp$id %in% common.samples
  resp <- resp[flag,]
  tmp <- as.data.frame(table(resp$PATIENT_ID))
  colnames(tmp) <- c("PATIENT_ID", "weight")
  tmp$weight <- 1/sqrt(tmp$weight)
  resp <- merge(resp, tmp)
  common.samples <- as.character(resp$id)
  if(!(ensg.gene %in% rownames(expr.mat))) {
    warning(paste0("UNEXPECTED: ", ensg.gene, " not in expr.mat\n"))
  }
  if(!(all(common.samples %in% colnames(expr.mat)))) {
    warning("UNEXPECTED columns not in expr.mat\n")
  }
  gene.expr <- as.vector(unlist(expr.mat[ensg.gene, common.samples]))
  df <- data.frame(expr = gene.expr, resp = resp$IC50, weights=resp$weight)
  ggplot(df, aes(x = expr, y = resp, weight=weights)) + 
    geom_point() + 
    xlab(xlab) + ylab(ylab) +
    stat_smooth(method = "lm", col = "red")  
}

ohsu.fimm.expr.corr <- merge(ohsu.expr.corr, fimm.expr.corr, by.x = c("ID_Drug", "gene"), by.y = c("DRUG_ID", "gene"), suffixes = c(".ohsu", ".fimm"))

subset(ohsu.fimm.expr.corr, (pval.ohsu < 0.05) & (pval.fimm < 0.05))
## ID_Drug            gene          pval.ohsu             r2.ohsu          pval.fimm           r2.fimm
## 100 FIMM003783 ENSG00000094631 0.0419410434254526  0.0494985601145971 0.0330721048250537 0.228393094777627
## 127 FIMM003794 ENSG00000122025 0.0140628306932017 0.00833458654791967 0.0449806673025468  0.24180451968676
ohsu.fimm.expr.corr$pval.fimm <- as.numeric(ohsu.fimm.expr.corr$pval.fimm)
ohsu.fimm.expr.corr$pval.ohsu <- as.numeric(ohsu.fimm.expr.corr$pval.ohsu)

## Plot pvalue distributions of FIMM and OHSU
pdf("ohsu-fimm-pval-dist.pdf")
g1 <- ggplot(ohsu.fimm.expr.corr, aes(x=pval.fimm)) + geom_histogram() + xlab("FIMM p-value")
g2 <- ggplot(ohsu.fimm.expr.corr, aes(x=pval.ohsu)) + geom_histogram() + xlab("OHSU p-value")
grid.arrange(g1, g2)
d <- dev.off()

ohsu.fimm.expr.corr$drug.gene <- unlist(apply(ohsu.fimm.expr.corr[,c("ID_Drug", "gene")], 1, function(row) paste0(row[1], "-", row[2])))
## Make Venn diagram of FIMM vs OHSU significant
library(gplots)
vennList = list("FIMM"=ohsu.fimm.expr.corr$drug.gene[ohsu.fimm.expr.corr$pval.fimm < 0.05], "OHSU"=ohsu.fimm.expr.corr$drug.gene[ohsu.fimm.expr.corr$pval.ohsu < 0.05])
pdf("ohsu-fimm-sig-venn.pdf")
venn(vennList)
d <- dev.off()

old.nrow <- nrow(ohsu.fimm.expr.corr)
ohsu.fimm.expr.corr <- merge(ohsu.fimm.expr.corr, unique(fimm.ohsu.drug.annotations[,c("ID_Drug", "DRUG_NAME")]))
if(old.nrow != nrow(ohsu.fimm.expr.corr)) {
  warning("UNEXPECTED drug correlation results grew\n")
}

old.nrow <- nrow(ohsu.fimm.expr.corr)
ohsu.fimm.expr.corr <- merge(ohsu.fimm.expr.corr, unique(sym.to.ensg), by.x = "gene", by.y = "ENSG")
if(old.nrow != nrow(ohsu.fimm.expr.corr)) {
  warning("UNEXPECTED drug correlation results grew\n")
}

## Plot few cases that have significant correlation in both
both.sig <- subset(ohsu.fimm.expr.corr, (pval.ohsu < 0.05) & (pval.fimm < 0.05))
for(i in 1:nrow(both.sig)) {
  ensg.gene <- both.sig$gene[i]
  drug.id <- both.sig$ID_Drug[i]
  pval.ohsu <- both.sig$pval.ohsu[i]
  pval.fimm <- both.sig$pval.fimm[i]
  
  gene.sym <- sym.to.ensg$SYMBOL[sym.to.ensg$ENSG == ensg.gene]
  drug.name <- fimm.ohsu.drug.annotations$DRUG_NAME[fimm.ohsu.drug.annotations$ID_Drug == drug.id]
  g1 <- plot.fimm.expr.vs.ic50(fimm.expr.df, drug.id, ensg.gene, fimm.expr.common, paste0("Log2 CPM (", gene.sym, ")"), paste0("Log10 IC50 (", drug.name, ")"))
  g1 <- g1 + ggtitle(paste0("FIMM (pval = ", pval.fimm, ")"))
  g2 <- plot.ohsu.expr.vs.ic50(ohsu.expr.df, drug.id, ensg.gene, ohsu.expr.common, paste0("Log2 CPM (", gene.sym, ")"), paste0("Log10 IC50 (", drug.name, ")"))
  g2 <- g2 + ggtitle(paste0("OHSU (pval = ", pval.ohsu, ")"))
  pdf(paste0("both-sig-", ensg.gene, "-", drug.id, ".pdf"))
  grid.arrange(g1, g2)
  d <- dev.off()
}

## Plot the top hits in each
sig <- subset(ohsu.fimm.expr.corr, pval.ohsu < 0.01)
for(i in 1:nrow(sig)) {
  ensg.gene <- sig$gene[i]
  drug.id <- sig$ID_Drug[i]
  pval.ohsu <- sig$pval.ohsu[i]
  pval.fimm <- sig$pval.fimm[i]
  
  gene.sym <- sym.to.ensg$SYMBOL[sym.to.ensg$ENSG == ensg.gene]
  drug.name <- fimm.ohsu.drug.annotations$DRUG_NAME[fimm.ohsu.drug.annotations$ID_Drug == drug.id]
  g1 <- plot.fimm.expr.vs.ic50(fimm.expr.df, drug.id, ensg.gene, fimm.expr.common, paste0("Log2 CPM (", gene.sym, ")"), paste0("Log10 IC50 (", drug.name, ")"))
  g1 <- g1 + ggtitle(paste0("FIMM (pval = ", pval.fimm, ")"))
  g2 <- plot.ohsu.expr.vs.ic50(ohsu.expr.df, drug.id, ensg.gene, ohsu.expr.common, paste0("Log2 CPM (", gene.sym, ")"), paste0("Log10 IC50 (", drug.name, ")"))
  g2 <- g2 + ggtitle(paste0("OHSU (pval = ", pval.ohsu, ")"))
  pdf(paste0("ohsu-sig-", ensg.gene, "-", drug.id, ".pdf"))
  grid.arrange(g1, g2)
  d <- dev.off()
}

sig <- subset(ohsu.fimm.expr.corr, pval.fimm < 0.01)
for(i in 1:nrow(sig)) {
  ensg.gene <- sig$gene[i]
  drug.id <- sig$ID_Drug[i]
  pval.ohsu <- sig$pval.ohsu[i]
  pval.fimm <- sig$pval.fimm[i]
  
  gene.sym <- sym.to.ensg$SYMBOL[sym.to.ensg$ENSG == ensg.gene]
  drug.name <- fimm.ohsu.drug.annotations$DRUG_NAME[fimm.ohsu.drug.annotations$ID_Drug == drug.id]
  g1 <- plot.fimm.expr.vs.ic50(fimm.expr.df, drug.id, ensg.gene, fimm.expr.common, paste0("Log2 CPM (", gene.sym, ")"), paste0("Log10 IC50 (", drug.name, ")"))
  g1 <- g1 + ggtitle(paste0("FIMM (pval = ", pval.fimm, ")"))
  g2 <- plot.ohsu.expr.vs.ic50(ohsu.expr.df, drug.id, ensg.gene, ohsu.expr.common, paste0("Log2 CPM (", gene.sym, ")"), paste0("Log10 IC50 (", drug.name, ")"))
  g2 <- g2 + ggtitle(paste0("OHSU (pval = ", pval.ohsu, ")"))
  pdf(paste0("fimm-sig-", ensg.gene, "-", drug.id, ".pdf"))
  grid.arrange(g1, g2)
  d <- dev.off()
}


## Plot the best cases in both

## - for intersected drugs, use symbols/ensg for genes


## - look at correlation of expr in gene targets with IC50
## - look at correlation of mutant in gene targets with IC50 (only those that have a lot of mutations)

## - 10-fold cross validation of ohsu using expr
## - 10-fold cross validation of ohsu using genomic

## Get some patients for Cristina
uniq.ohsu <- unique(ohsu.raw.drug.data[,c("inhibitor", "patient_id", "diagnosis", "specific_diagnosis", "lab_id", "specimen_type", "run_type", "replicant", "time_of_read")])

flag <- (uniq.ohsu$lab_id %in% colnames(ohsu.expr.data)) & (uniq.ohsu$lab_id %in% ohsu.beat.maf.names)

write.table(uniq.ohsu[flag, ], sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE, file="ohsu-drug-screens-with-rna-and-dna-seq.tsv")