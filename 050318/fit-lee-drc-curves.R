suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(synapseClient))
suppressPackageStartupMessages(library("foreach"))
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("ggplot2"))
##suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("MonoInc"))
suppressPackageStartupMessages(library("data.table"))

source("../common/dss.R")
source("../common/drc-fit.R")
source("../common/plotting.R")

num.cores <- detectCores()
if(!is.na(num.cores) && (num.cores > 1)) {
  suppressPackageStartupMessages(library("doMC"))
  cat(paste("Registering ", num.cores-1, " cores.\n", sep=""))
  registerDoMC(cores=(num.cores-1))
}

synapseLogin()

cat("Use DRC to fit L.4 (4-parameter logistic) and LL.4 (4-paramter log logistic) functions to Lee et al.\n")

## BEGIN setup

## plan
## what is concentration range?
## how often do we drop the first point? last point?
## make sure we drop these outside of range from auc calculation

source("lee-setup-050318.R")

prefix <- "lee"
data.set <- "lee"

conc.col <- conc.cols[[data.set]]
response.col <- response.cols[[data.set]]
response.is.viability <- responses.are.viabilities[[data.set]]
drug.col <- data.set.drug.id.cols[[data.set]]
drug.screen.cols <- screen.id.cols[[data.set]]
experimental.factors <- experimental.factor.cols[[data.set]]
parentId <- processed.folder.synIds[[data.set]]

output.file <- paste0("lee.l4.and.ll4.fits.050318.tsv")

## Read in raw drug response data
synId <- raw.drug.synIds[[data.set]]
obj <- synGet(id=synId, downloadFile = TRUE)
file <- getFileLocation(obj)
drug.data.long <- read.table(file, sep="\t", header=TRUE, as.is=TRUE)

## Read in the expression data
synId <- data.set.expr.synIds[[data.set]]
obj <- synGet(id=synId, downloadFile = TRUE)
file <- getFileLocation(obj)
expr <- read.table(file, sep="\t", header=TRUE, as.is=TRUE)

## Load FIMM metadata (including relapse/refractory/diagnosis)
drug.metadata <- NULL
if(FALSE) { ## Don't yet have metadata for CM
synId <- "syn8270594"
obj <- synGet(id=synId, downloadFile = TRUE, downloadLocation = ".")
file <- getFileLocation(obj)
drug.metadata <- read.table(file, header=TRUE, sep=",")
}

all.drugs <- unique(drug.data.long[, drug.col])
all.cols <- c(drug.screen.cols, conc.col, response.col)

## Limit fo AML
if(!is.null(drug.metadata)) {
  aml.drug.metadata <- subset(drug.metadata, grepl(diagnosis, pattern="AML"))
  drug.data.long <- subset(drug.data.long, (SCREEN_ID %in% rownames(aml.drug.metadata)))
}

drugs.to.exclude <- c()

## END setup

## Restrict to samples for which we have both expression and drug response
common.samples <- intersect(drug.data.long[, expr.patient.id.cols[[data.set]]], colnames(expr))
expr <- expr[, common.samples]
drug.data.long <- drug.data.long[drug.data.long[, expr.patient.id.cols[[data.set]]] %in% common.samples, ]

source("../common/annotate-drug-response-outliers.R")

## Define the acceptable response range at the minimum concentration (min.response.at.min.conc, max.response.at.min.conc)

## Define the distribution of the response at the minimum concentration
flag <- !(drug.data.long[, drug.col] %in% drugs.to.exclude)
drug.data.long <- drug.data.long[flag, ]
res <- define.response.at.min.conc.distribution(drug.data.long, drug.screen.cols, conc.col, response.col)
responses.at.min.conc <- res[["responses.at.min.conc"]][, response.col]
## min.conc.mean <- res[["mean"]]
## min.conc.sd <- res[["sd"]]
## min.response.at.min.conc <- min.conc.mean - 3 * min.conc.sd
## max.response.at.min.conc <- min.conc.mean + 3 * min.conc.sd

lst <- list("Response at min conc" = responses.at.min.conc)
indices <- 1:length(lst)
names(indices) <- names(lst)
file <- paste0("min-conc-distributions-", data.set, ".pdf")
pdf(file)
df <- Reduce("rbind", llply(indices, .fun = function(i) data.frame(color = names(lst)[i], x = lst[[i]])))
df <- data.frame(color = "foo", x = drug.data.long[, response.col])
g <- ggplot(data = df, aes(x = x, color = color))
g <- g + geom_density()
g <- g + xlab(response.col)
g <- g + ggtitle(toupper(data.set)) + theme(legend.title=element_blank())
g
d <- dev.off()
l_ply(indices, 
      .fun = function(i) {
         qnt <- quantile(lst[[i]], probs = c(0.05, 0.95))
         cat(paste0(names(lst)[i], " mean = ", mean(lst[[i]]), " sd = ", sd(lst[[i]]), " ",
             paste(unlist(lapply(1:length(qnt), function(i) paste0(names(qnt)[i], ": ", unname(qnt[i])))), collapse = " "), "\n"))

      })

res <- find.extremal.responses.qqplots_(drug.data.long, response.col)
min_response <- res$min.response
max_response <- res$max.response

drug.data.long <- annotate.drug.response.outliers_(drug.data.long, response.col, min_response, max_response)

## Plot responses vs each of the experimental factors.  Detect any outliers that we should exclude based on these experimental factors.
## In particular, for each level of each factor, calculate the fraction of samples that have response outside the range min.response to max.response.
## If that fraction is an outlying (and high) fraction across all such fractions for that factor, exclude the level. 
source("../common/drug-response-quality-control.R")
exclude <- correlate.drug.response.with.experimental.factors(drug.data.long, experimental.factors, response.col, min_response, max_response, prefix = prefix)

## Perform the fits, after excluding outliers
fits <- list()
fcts <- c("LL.4", "L.4")
outlier.col <- "outlier.response"
for(fct in fcts) {
  cat(paste0("Fitting drug response using ", fct, "\n"))
  flag <- !(drug.data.long[, drug.col] %in% drugs.to.exclude)
  fits[[fct]] <- fit.drc(drug.data.long[flag, c(all.cols, outlier.col)], fct = fct, drug.screen.cols = drug.screen.cols, conc.nM.col = conc.col, 
                         response.col = response.col, response.is.viability = response.is.viability, outlier.col = outlier.col, compute.auc = FALSE)

  cat(paste0("Done fitting drug response using ", fct, "\n"))

  ## Allow at most one outlier to be dropped per curve
  ## For some reason, some of these curves have many more responses -- possibly multiple replicates have
  ## been lumped together.  Filter these out. 
  num.responses <- unlist(lapply(fits[[fct]]$uniq.concs.nM, function(str) length(unlist(strsplit(str, split=",[ ]*")))))
  tbl <- as.data.frame(table(num.responses))
  expected.num.responses <- as.numeric(as.character(tbl$num.responses[tbl$Freq == max(tbl$Freq)]))
  accepted.responses <- (expected.num.responses-1):expected.num.responses

  ## Require at least 5 pts
  accepted.responses <- accepted.responses[accepted.responses >= 5]
  cat(paste0("Limiting curves to those with the following number of responses: ", paste(accepted.responses, collapse=", "), "\n"))
  exclusion.flag <- !(as.numeric(num.responses) %in% accepted.responses)
  exclusion.col <- "num.pts.exclude"
  exclusion.cols <- c(exclusion.col)
  factor.cols <- c(NA)
  fits[[fct]][, exclusion.col] <- FALSE
  fits[[fct]][exclusion.flag, exclusion.col] <- TRUE
 
  ## Exclude based on cutoffs established above
  for(nm in names(exclude)) {
    exclusion.col <- paste0(nm, ".exclude")
    exclusion.flag <- !is.na(fits[[fct]][, nm]) & (fits[[fct]][, nm] %in% exclude[[nm]])
    fits[[fct]][, exclusion.col] <- FALSE
    fits[[fct]][exclusion.flag, exclusion.col] <- TRUE
    exclusion.cols <- c(exclusion.cols, exclusion.col)
    factor.cols <- c(factor.cols, nm)
  }

  if(FALSE) {
    any.neg.response <- unlist(lapply(fits[[fct]]$all.responses, function(str) any(as.numeric(unlist(strsplit(str, split=",[ ]*"))) < -30)))
    monotone <- unlist(lapply(fits[[fct]]$all.responses, function(str) {
                                                            vec <- as.numeric(unlist(strsplit(str, split=",[ ]*")))
                                                            monotonic(vec, direction='inc') || monotonic(vec, direction='dec')
                                                          }))
  }

  ## Exclude based on gof
  ## Establish a gof cutoff by comparing approximately monotonic data vs non-approximately monotonic
  ## Consider two consecutive components of a vector v[i] and v[i+1] to be approximately monotonic increasing
  ## if (v[i] - v[i+1]) < 10 -- i.e., allowing a tolerance of 10%
  ## Consider two consecutive components of a vector v[i] and v[i+1] to be approximately _non_-monotonic increasing
  ## if (v[i] - v[i+1]) > 30 -- i.e., allowing a tolerance of 30%
  any.exclude <- unlist(apply(fits[[fct]][, grepl(pattern="exclude", x=colnames(fits[[fct]]))], 1, function(row) any(row)))
  tmp <- fits[[fct]][!any.exclude,]
  approx.monotone <- unlist(lapply(tmp$all.responses, function(str) {
                                                          vec <- as.numeric(unlist(strsplit(str, split=",[ ]*")))
                                                          approximately.monotonic(vec, direction='inc', tolerance = 10) || 
                                                             approximately.monotonic(vec, direction='dec', tolerance = 10)
                                                        }))
  approx.non.monotone <- unlist(lapply(tmp$all.responses, function(str) {
                                                          vec <- as.numeric(unlist(strsplit(str, split=",[ ]*")))
                                                          !approximately.monotonic(vec, direction='inc', tolerance = 30) &&
                                                             !approximately.monotonic(vec, direction='dec', tolerance = 30)
                                                        }))
  monotone <- rep(NA, nrow(tmp))
  flag <- approx.monotone == TRUE
  monotone[flag] <- TRUE
  flag <- approx.non.monotone == TRUE
  monotone[flag] <- FALSE
  flag <- !is.na(monotone)

  file <- paste0(prefix, "-", fct, "-gof-cutoff.pdf")
  pdf(file, onefile=FALSE)
  tbl <- data.frame(gof = as.numeric(tmp$gof[flag]), condition = monotone[flag])
##  gof.cutoff <- round(median(tbl$gof[tbl$condition == FALSE]), digits=2)
  suppressPackageStartupMessages(library(maxstat))
  formula <- as.formula("condition ~ gof")
  mt <- maxstat.test(formula, data=tbl, smethod="Wilcoxon", pmethod="none")
  gof.cutoff <- mt$estimate
  cat(paste0("Using GOF cutoff: ", gof.cutoff, "\n"))
  g1 <- ggplot(tbl)
  g1 <- g1 + geom_boxplot(aes(x = condition, y = gof))
  g1 <- g1 + xlab("Approximately Monotonic")

  g2 <- ggplot(tbl)
  g2 <- g2 + geom_density(aes(x = gof))
  g2 <- g2 + geom_vline(xintercept = gof.cutoff)
  grid.arrange(g1, g2, nrow=1, top = paste0(prefix, ": GOF cutoff = ", gof.cutoff))  
  d <- dev.off()

  exclusion.flag <- as.numeric(fits[[fct]]$converged) == 0 
  exclusion.flag <- exclusion.flag | ( !is.na(fits[[fct]]$gof) & (as.numeric(fits[[fct]]$gof) < gof.cutoff) )
  exclusion.col <- "gof.exclude"
  exclusion.cols <- c(exclusion.cols, exclusion.col)
  factor.cols <- c(factor.cols, NA)
  fits[[fct]][, exclusion.col] <- FALSE
  fits[[fct]][exclusion.flag, exclusion.col] <- TRUE

  any.exclude <- unlist(apply(fits[[fct]][, grepl(pattern="exclude", x=colnames(fits[[fct]]))], 1, function(row) any(row)))
  exclusion.flag <- any.exclude
  exclusion.col <- "any.exclude"
  exclusion.cols <- c(exclusion.cols, exclusion.col)
  factor.cols <- c(factor.cols, NA)
  fits[[fct]][, exclusion.col] <- FALSE
  fits[[fct]][exclusion.flag, exclusion.col] <- TRUE

##  for(exclusion.col in colnames(fits[[fct]])[grepl(pattern="exclude", x=colnames(fits[[fct]]))]) {
  for(indx in 1:length(exclusion.cols)) {
      exclusion.col <- exclusion.cols[indx]
      factor.col <- factor.cols[indx]
      cat(paste0("\n\nNum excluded: ", exclusion.col, "\n"))
      print(table(fits[[fct]][, exclusion.col]))
      exclude.flag <- fits[[fct]][, exclusion.col] == TRUE
      if(any(exclude.flag)) {
        if(!is.na(factor.col)) {
          excluded <- unique(fits[[fct]][exclude.flag, factor.col])
          cat(paste0("Num ", factor.col, " excluded: ", length(excluded), "\n"))
          cat(paste0(paste(excluded, collapse=", "), "\n"))
        }
      }
  }
}

common.cols <- c(drug.screen.cols, "min.conc.nM", "max.conc.nM", "all.concs.nM", "uniq.concs.nM", "all.responses",
                 "min.conc.w.outliers.nM", "max.conc.w.outliers.nM", "all.concs.w.outliers.nM", "uniq.concs.w.outliers.nM", "all.responses.w.outliers")
fits <- merge(fits[["LL.4"]], fits[["L.4"]], by = common.cols, suffixes = c(".ll4", ".l4"), all = TRUE)

my.dup <- function(x) duplicated(x, fromLast = TRUE) | duplicated(x, fromLast = FALSE)

## Make sure each drug is tested over a common range across patients

cols <- c(drug.col, "min.conc.w.outliers.nM", "max.conc.w.outliers.nM")
df <- fits[!fits$any.exclude.ll4, cols]
uniq.conc.orig <- ddply(df, .variables = cols, .fun = function(df) df$cnt <- nrow(df))
colnames(uniq.conc.orig) <- c(cols, "cnt")

cat("Note that the original ranges (before dropping any points) are the same for all of the drugs (modulo rounding errors)\n")
print(uniq.conc.orig[my.dup(uniq.conc.orig[, drug.col]),])

cols <- c(drug.col, "min.conc.nM", "max.conc.nM")
df <- fits[!fits$any.exclude.ll4, cols]
uniq.conc <- ddply(df, .variables = cols, .fun = function(df) df$cnt <- nrow(df))
colnames(uniq.conc) <- c(cols, "cnt")
uniq.conc <- uniq.conc[order(uniq.conc[, drug.col]),]
print(uniq.conc[my.dup(uniq.conc[, drug.col]),])


source("../common/restrict-concentration-to-common-range.R")

dummy.drug.map <- data.frame("inhibitor" = unique(fits[, drug.col]))
colnames(dummy.drug.map)[colnames(dummy.drug.map) == "inhibitor"] <- drug.col
dummy.drug.map.drug.id.cols <- list("data.set" = drug.col)
names(dummy.drug.map.drug.id.cols) <- c(data.set)
fit.list <- list("data.set" = fits)
names(fit.list) <- data.set
ret <- determine.drug.ranges(fit.list, screen.id.cols, data.set.drug.id.cols, 
                             dummy.drug.map, dummy.drug.map.drug.id.cols, shared.conc.postfix = ".auc.conc") 

print(uniq.conc[my.dup(uniq.conc[, drug.col]),])
flag <- ret$shared.drug.ranges[["LL.4"]][, drug.col]  %in% uniq.conc[my.dup(uniq.conc[, drug.col]),drug.col]
m <- merge(uniq.conc[my.dup(uniq.conc[, drug.col]), ], ret$shared.drug.ranges[["LL.4"]], by = drug.col)

## Compute the auc
tmp <- ret$fits.with.concs[[data.set]]
for(nm in names(tmp)) {
  ## Drop the columns with .ll4 prefix in the L.4 fits and vice versa
  pattern = "\\.l4"
  if(nm == "L.4") { pattern = "\\.ll4" }
  tmp[[nm]] <- tmp[[nm]][, !grepl(x = colnames(tmp[[nm]]), pattern = pattern)]
}

fits <- merge(tmp[["LL.4"]], tmp[["L.4"]], by = common.cols, all = TRUE)

source("../common/dss.R")

fcts <- list("LL.4" = "LL.4", "L.4" = "L.4")
postfixes <- list("LL.4" = ".ll4", "L.4" = ".l4")

for(nm in names(fcts)) {  
  min.conc.col <- paste0("min.auc.conc", postfixes[[nm]])
  max.conc.col <- paste0("max.auc.conc", postfixes[[nm]])
  b.param.col <- paste0("b", postfixes[[nm]])
  c.param.col <- paste0("c", postfixes[[nm]])
  d.param.col <- paste0("d", postfixes[[nm]])
  e.param.col <- paste0("e", postfixes[[nm]])
  converged.col <- paste0("converged", postfixes[[nm]])
  auc.output.col <- paste0("auc", postfixes[[nm]])
  shifted.auc.output.col <- paste0("shifted.auc", postfixes[[nm]])
  fits <- compute.all.auc(fits, fct = fcts[[nm]], min.conc.col = min.conc.col, max.conc.col = max.conc.col, 
                          b.param.col = b.param.col, c.param.col = c.param.col, d.param.col = d.param.col, e.param.col = e.param.col, 
                          converged.col = converged.col, auc.output.col = auc.output.col, shifted.auc.output.col = shifted.auc.output.col)
}


cat("Saving image\n")
save.image(rdata.file)

write.table(file=output.file, fits, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
f <- File(output.file, parentId = parentId, synapseStore = TRUE)
synStore(f, executed = NULL, forceVersion = FALSE)

cat("Done fitting drug response data using L.4 and LL.4\n")

q()
