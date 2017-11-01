## Set the download and output directory
download.path <- "../input/"
if (!file.exists(download.path)) {
  dir.create(download.path)
}

output.path <- "output/"
if (!file.exists(output.path)) {
  dir.create(output.path)
}

suppressPackageStartupMessages(library("foreach"))
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(tidyr))

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

num.cores <- detectCores()
if(!is.na(num.cores) && (num.cores > 1)) {
  suppressPackageStartupMessages(library("doMC"))
  cat(paste("Registering ", num.cores-1, " cores.\n", sep=""))
  registerDoMC(cores=(num.cores-1))
}

synapseLogin()

cat("Calculate AUC and DSS values based on 2-way intersection of concentrations\n")
cat("for each drug shared between FIMM and OHSU.\n")

## Load the map of drugs shared between OHSU and FIMM (FIMM_OHSU_Drug_Concentrations.xlsx)
synId <- "syn10163669"
obj <- synGet(synId, downloadFile = TRUE)
ohsu.fimm.drugs <- openxlsx:::read.xlsx(getFileLocation(obj), sheet=1)

ohsu.fimm.drugs <- ohsu.fimm.drugs[, c("FIMM_Batch_ID_Drug", "OHSU_DRUG_NAME", "FIMM_DRUG_NAME", "Mechanism.Targets", "Trade.names", "Class.explained", "Aliases")]

## path <- "fimm.dss.t0.tsv"
synId <- "syn10083488"
obj <- synGet(synId, downloadFile = TRUE)
file <- getFileLocation(obj)
fimm.dss.t0.orig <- as.data.frame(fread(file))

## Filter FIMM results to include only those drugs with at least 5 uniq concentration points
min.conc.pts <- 5
keep <- sapply(fimm.dss.t0.orig$uniq.concs, function(str) length(unique(unlist(strsplit(str, split=",")))) >= min.conc.pts)
cat(paste0("Keeping FIMM screens with at least ", min.conc.pts, " unique concentrations\n"))
print(table(keep))
fimm.dss.t0.orig <- fimm.dss.t0.orig[keep,]

## path <- "ohsu.dss.t0.tsv"
synId <- "syn10083494"
obj <- synGet(synId, downloadFile = TRUE)
file <- getFileLocation(obj)
ohsu.dss.t0.orig <- as.data.frame(fread(file))

min.conc.pts <- 5
max.conc.pts <- 7
keep <- unlist(apply(ohsu.dss.t0.orig[ ,c("all.concs.nM", "uniq.concs.nM")], 1, function(row) {
           ( length(unique(unlist(strsplit(row[1], split=",")))) <= max.conc.pts ) && ( length(unique(unlist(strsplit(row[2], split=",")))) >= min.conc.pts )
        }))
cat(paste0("Keeping OHSU screens with at least ", min.conc.pts, " unique concentrations and at most ", max.conc.pts, " total concentrations\n"))
print(table(keep))
ohsu.dss.t0.orig <- ohsu.dss.t0.orig[keep,]

## Confirm that all of the drugs in the shared drug annotation file are
## in the processed drug response file.
missing.drugs <- ohsu.fimm.drugs[!(ohsu.fimm.drugs$FIMM_Batch_ID_Drug %in% fimm.dss.t0.orig$DRUG_ID),]
if(nrow(missing.drugs) == 0) {
  cat("As expected, all drugs in shared file are in FIMM drug response data file\n")
} else {
  cat("Missing some drugs in shared file from FIMM drug response file\n")
  q(status = -1)
}

## Update as add/remove data sets
col.names <- list("OHSU_DRUG_NAME")
all.inhibitors <- list(unique(ohsu.dss.t0.orig$inhibitor))
data.set.names <- list("OHSU")

for(i in 1:length(col.names)) {
  missing.drugs <- ohsu.fimm.drugs[!(ohsu.fimm.drugs[, col.names[[i]]] %in% all.inhibitors[[i]]),]
  if(nrow(missing.drugs) == 0) {
    cat(paste0("As expected, all drugs in shared file are in ", data.set.names[[i]], " drug response data file\n"))
  } else {
    cat(paste0("Missing some drugs in shared file from ", data.set.names[[i]], " drug response file\n"))
    q(status = -1)
  }
}

my.dup <- function(x) {
  duplicated(x, fromLast=TRUE) | duplicated(x, fromLast=FALSE)
}

## Note (and drop) any drugs that have multiple entries (which would imply ambiguous mappings
## and/or inconsistent drug ranges)

## Update as add/remove data sets
drug.col.names <- list("OHSU_DRUG_NAME", "FIMM_Batch_ID_Drug")
dups <- rep(FALSE, nrow(ohsu.fimm.drugs))
for(i in 1:length(drug.col.names)) {
  dups <- dups | my.dup(ohsu.fimm.drugs[, drug.col.names[[i]]])
}
if(any(dups)) {
  cat("Dropping ambiguous mappings between drugs:\n")
  write.table(ohsu.fimm.drugs[dups,], sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
  ohsu.fimm.drugs <- ohsu.fimm.drugs[!dups,]
}

cat(paste0(nrow(ohsu.fimm.drugs), " common drugs\n"))
write.table(ohsu.fimm.drugs, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

## Update as add/remove data sets
## Restrict the processed results to those shared between both data sets
## And drop all of the pre-computed DSS columns.
fimm.dss.common <- subset(fimm.dss.t0.orig, DRUG_ID %in% ohsu.fimm.drugs$FIMM_Batch_ID_Drug)
fimm.dss.common <- fimm.dss.common[, !(grepl(x=colnames(fimm.dss.common), pattern="dss"))]
fimm.dss.common <- subset(fimm.dss.common, (converged.l4 == 1) | (converged.ll4 == 1))
rm(fimm.dss.t0.orig)

ohsu.dss.common <- subset(ohsu.dss.t0.orig, inhibitor %in% ohsu.fimm.drugs$OHSU_DRUG_NAME)
ohsu.dss.common <- ohsu.dss.common[, !(grepl(x=colnames(ohsu.dss.common), pattern="dss"))]
ohsu.dss.common <- subset(ohsu.dss.common, (converged.l4 == 1) | (converged.ll4 == 1))
rm(ohsu.dss.t0.orig)

## Update as add/remove data sets
## Note that ranges differ _within_ data sets.  It appears that the file
## FIMM_OHSU_Drug_Concentrations.xlsx was created by union'ing the concentrations of a drug
## across all patients and then finding the min and max.  That does not imply that all
## patients were tested for the drug in that range.  Instead, let's limit to the intersection,
## (i.e., max of mins and min of maxs), just as we do across data sets.

## First, note the occurrence of this in the FIMM raw data.
synId <- "syn8488824"
obj <- synGet(id=synId, downloadFile = TRUE, downloadLocation = download.path)
file <- getFileLocation(obj)
fimm.raw.dr <- read.xlsx(file)

## Drop any raw data that we did not fit.
fimm.raw.dr <- merge(fimm.raw.dr, fimm.dss.common[, c("SCREEN_ID", "PATIENT_ID", "DRUG_ID")])

## Compute the drug ranges from the _raw_ OHSU data.  i.e., before excluding outliers, which
## could conceivably affect the range.  However, limit to those screens for which we _attempted_ a fit.
## Do exclude those for which both L4 and LL4 failed
synId <- "syn8149180"
inh.obj <- synGet(synId, downloadFile=TRUE, downloadLocation = download.path)
ohsu.raw.dr <- fread(getFileLocation(inh.obj))
ohsu.raw.dr <- as.data.frame(ohsu.raw.dr)

## Drop any raw data that we did not fit.
ohsu.raw.dr <- merge(ohsu.raw.dr, ohsu.dss.common[, c("patient_id", "lab_id", "inhibitor", "replicant", "time_of_read")])

ohsu.raw.drug.ranges <- ddply(ohsu.raw.dr,
                              .variables = c("patient_id", "lab_id", "inhibitor", "replicant", "time_of_read"),
                              .fun = function(df) { 
                                       ## Adjust to be in nanomolar range
                                       v <- c(min(as.numeric(df$well_concentration)), max(as.numeric(df$well_concentration)))
                                       v <- v * 10^3
                                       names(v) <- c("min.conc.nM", "max.conc.nM")
                                       v 
                                     })

ohsu.raw.drug.ranges <- ddply(ohsu.raw.drug.ranges, .variables = c("inhibitor", "min.conc.nM", "max.conc.nM"),
                              .fun = function(df) {
                                       vec <- c(df$inhibitor[1], df$min.conc.nM[1], df$max.conc.nM[1], freq = nrow(df))
                                       names(vec) <- c("inhibitor", "min.conc.nM", "max.conc.nM", "freq")
                                       vec
                                 })

## In the OHSU data, we often see a drug that has several different lower concentrations, all of which
## are quite similar.  Generally, looks like a rounding issue.
cutoff <- 0.01
cutoff <- 0
cat(paste0("OHSU drugs with min or max concentrations differing by > ", cutoff, "\n"))
d_ply(ohsu.raw.drug.ranges, .variables = c("inhibitor"),
                             .fun = function(df) {
                                      df$min.conc.nM <- as.numeric(df$min.conc.nM); df$max.conc.nM <- as.numeric(df$max.conc.nM)
                                      if( ( (max(df$min.conc.nM) - min(df$min.conc.nM)) > cutoff ) || ( (max(df$max.conc.nM) - min(df$max.conc.nM)) > cutoff ) ) {
                                         cat("\n")
                                         write.table(df, col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
                                      }
                                    })

## Hence, take the most restrictive range (i.e., the intersection) for the OHSU drugs
ohsu.raw.drug.ranges <- ddply(ohsu.raw.drug.ranges, .variables = c("inhibitor"),
                              .fun = function(df) {
                                       v <- c(max(as.numeric(df$min.conc.nM)), min(as.numeric(df$max.conc.nM)))
                                       names(v) <- c("min.conc.nM", "max.conc.nM")
                                       v
                              })

## Do the same for FIMM
fimm.raw.drug.ranges <- ddply(fimm.raw.dr,
                              .variables = c("PATIENT_ID", "SCREEN_ID", "DRUG_SET", "DRUG_ID"),
                              .fun = function(df) { 
                                       v <- c(min(as.numeric(df$CONCENTRATION)), max(as.numeric(df$CONCENTRATION)))
                                       names(v) <- c("min.conc.nM", "max.conc.nM")
                                       v 
                                     })

## Now look at the range of one drug in particular, as generated from the raw data.  
## Note some go from 1.0 to 10,000 and some from 0.1 to 1,000
drug <- "FIMM023833"
cat(paste0("FIMM drug range for ", drug, " calculated from raw data on a per-patient basis\n"))
tail(subset(fimm.raw.drug.ranges, DRUG_ID == drug))

## Collapse duplicate ranges and drop patient ids
fimm.raw.drug.ranges <- ddply(fimm.raw.drug.ranges, .variables = c("DRUG_ID", "min.conc.nM", "max.conc.nM"),
                              .fun = function(df) {
                                       vec <- c(df$DRUG_ID[1], df$min.conc.nM[1], df$max.conc.nM[1], freq = nrow(df))
                                       names(vec) <- c("DRUG_ID", "min.conc.nM", "max.conc.nM", "freq")
                                       vec
                                 })


cat(paste0("FIMM drugs with min or max concentrations differing by >= ", cutoff, "\n"))
d_ply(fimm.raw.drug.ranges, .variables = c("DRUG_ID"),
                             .fun = function(df) {
                                      df$min.conc.nM <- as.numeric(df$min.conc.nM); df$max.conc.nM <- as.numeric(df$max.conc.nM)
                                      if(nrow(df) > 1) {
                                        if( ( (max(df$min.conc.nM) - min(df$min.conc.nM)) >= cutoff ) || ( (max(df$max.conc.nM) - min(df$max.conc.nM)) >= cutoff ) ) {
                                           cat("\n")
                                           write.table(df, col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
                                        }
                                      }
                                    })

## Discard up to max.screens.to.discard if it extends the drug concentration range
max.screens.to.discard <- 8
fimm.raw.drug.ranges <- ddply(fimm.raw.drug.ranges, .variables = c("DRUG_ID"),
                              .fun = function(df) {
                                       ## Consider the ranges in increasing order of their frequency--
                                       ## i.e., we will consider dropping the least frequent ranges first
                                       df <- df[order(as.numeric(df$freq), decreasing=FALSE),]
                                       flag <- as.numeric(df$freq) < max.screens.to.discard
                                       if((nrow(df) > 1) && any(flag)) {
                                         new.df <- df[!flag,]
                                         min.conc <- max(as.numeric(new.df$min.conc.nM))
                                         max.conc <- min(as.numeric(new.df$max.conc.nM))
                                         indices <- which(flag)
                                         tot.dropped <- 0
                                         for(indx in indices) {
                                           this.min <- as.numeric(df$min.conc.nM[indx])
                                           this.max <- as.numeric(df$max.conc.nM[indx])
                                           if( ( (tot.dropped + as.numeric(df$freq[indx])) < max.screens.to.discard ) && ( (this.min < min.conc) || (this.max > max.conc) ) ) {
                                             ## This would reduce the range, drop it.
                                             tot.dropped <- tot.dropped + as.numeric(df$freq[indx])
                                           } else {
                                             new.df <- rbind(new.df, df[indx,])
                                           }
                                         }
                                         df <- new.df
                                       }
                                       v <- c(max(as.numeric(df$min.conc.nM)), min(as.numeric(df$max.conc.nM)))
                                       names(v) <- c("min.conc.nM", "max.conc.nM")
                                       v
                              })

## Merge the ranges for all data sets together
shared.drug.anno <- ohsu.fimm.drugs
ohsu.raw.drug.ranges <- merge(shared.drug.anno, ohsu.raw.drug.ranges, by.y = "inhibitor", by.x = "OHSU_DRUG_NAME")
colnames(ohsu.raw.drug.ranges)[colnames(ohsu.raw.drug.ranges) == "min.conc.nM"] <- "ohsu.min.conc.nM"
colnames(ohsu.raw.drug.ranges)[colnames(ohsu.raw.drug.ranges) == "max.conc.nM"] <- "ohsu.max.conc.nM"

colnames(fimm.raw.drug.ranges)[colnames(fimm.raw.drug.ranges) == "min.conc.nM"] <- "fimm.min.conc.nM"
colnames(fimm.raw.drug.ranges)[colnames(fimm.raw.drug.ranges) == "max.conc.nM"] <- "fimm.max.conc.nM"

raw.drug.ranges <- merge(ohsu.raw.drug.ranges, fimm.raw.drug.ranges[, c("DRUG_ID", "fimm.min.conc.nM", "fimm.max.conc.nM")], by.x = "FIMM_Batch_ID_Drug", by.y = "DRUG_ID")

## Update as add/remove data sets
## DSS code will use min.conc and max.conc columns
cols <- c("ohsu.min.conc.nM", "fimm.min.conc.nM")
raw.drug.ranges$min.conc <- unlist(apply(raw.drug.ranges[, cols], 1, function(row) max(as.numeric(row))))

cols <- c("ohsu.max.conc.nM", "fimm.max.conc.nM")
raw.drug.ranges$max.conc <- unlist(apply(raw.drug.ranges[, cols], 1, function(row) min(as.numeric(row))))

## Store the drugs and concentration ranges common to OHSU and FIMM
## Update as add/remove data sets

## Merge the processed drug response tables with the range tables and exclude any screens that do 
## not completely cover the min/max range
fimm.dss.common <- merge(fimm.dss.common[, !(colnames(fimm.dss.common) %in% c("min.conc.nM", "max.conc.nM"))], raw.drug.ranges, by.x = "DRUG_ID", by.y = "FIMM_Batch_ID_Drug")
keep <- unlist(apply(fimm.dss.common[, c("uniq.concs.nM", "min.conc", "max.conc")], 1, 
                     function(row) {
                       concs <- as.numeric(unlist(strsplit(row[1], split=",")))
                       min.conc.nM <- as.numeric(row[2])
                       max.conc.nM <- as.numeric(row[3])
                       (min.conc.nM < max.conc.nM) & (min(concs) < max(concs)) & (min(concs) <= min.conc.nM) & (max(concs) >= max.conc.nM)
                     }))
fimm.dss.common <- fimm.dss.common[keep,]

ohsu.dss.common <- merge(ohsu.dss.common[, !(colnames(ohsu.dss.common) %in% c("min.conc.nM", "max.conc.nM"))], raw.drug.ranges, by.x = "inhibitor", by.y = "OHSU_DRUG_NAME")
keep <- unlist(apply(ohsu.dss.common[, c("uniq.concs.nM", "min.conc", "max.conc")], 1, 
                     function(row) {
                       concs <- as.numeric(unlist(strsplit(row[1], split=",")))
                       min.conc.nM <- as.numeric(row[2])
                       max.conc.nM <- as.numeric(row[3])
                       (min.conc.nM < max.conc.nM) & (min(concs) < max(concs)) & (min(concs) <= min.conc.nM) & (max(concs) >= max.conc.nM)
                     }))
ohsu.dss.common <- ohsu.dss.common[keep,]

save.image(".Rdata")

## Recalculate DSS using both L.4 and LL.4 fits (for FIMM and OHSU)
## Update as add/remove data sets

## Separate the tables into L4 and LL4 fits and strip off the .l4 and .ll4 prefixes,
## since the DSS code below expect parameters named b, c, d, e, not b.ll4, c.ll4, d.ll4, e.ll4
## Note that parameter e = ic50 (not log ic50) for both L.4 and LL.4 fits
l4.ohsu.common <- ohsu.dss.common[, !grepl(x = colnames(ohsu.dss.common), pattern="\\.ll4")]
ll4.ohsu.common <- ohsu.dss.common[, !grepl(x = colnames(ohsu.dss.common), pattern="\\.l4")]
colnames(l4.ohsu.common) <- gsub(x = colnames(l4.ohsu.common), pattern=".l4", replacement="")
colnames(ll4.ohsu.common) <- gsub(x = colnames(ll4.ohsu.common), pattern=".ll4", replacement="")
l4.ohsu.common$ic50 <- as.numeric(l4.ohsu.common$e)
ll4.ohsu.common$ic50 <- as.numeric(ll4.ohsu.common$e)

l4.fimm.common <- fimm.dss.common[, !grepl(x = colnames(fimm.dss.common), pattern="\\.ll4")]
ll4.fimm.common <- fimm.dss.common[, !grepl(x = colnames(fimm.dss.common), pattern="\\.l4")]
colnames(l4.fimm.common) <- gsub(x = colnames(l4.fimm.common), pattern=".l4", replacement="")
colnames(ll4.fimm.common) <- gsub(x = colnames(ll4.fimm.common), pattern=".ll4", replacement="")
l4.fimm.common$ic50 <- as.numeric(l4.fimm.common$e)
ll4.fimm.common$ic50 <- as.numeric(ll4.fimm.common$e)

source("../common/dss.R")

remove.outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

plot.l4.vs.ll4 <- function(data, col, col.name, use.log = FALSE, rm.outliers = FALSE) {
  x <- data[, paste0(col, ".l4")]
  y <- data[, paste0(col, ".ll4")]
  if(rm.outliers) {
    x <- remove.outliers(x)
    y <- remove.outliers(y)
  }
  df <- data.frame(x = x, y = y)
  xlab <- paste0(col.name, " (L.4)")
  ylab <- paste0(col.name, " (LL.4)")
  if(use.log) {
    df$x <- log(df$x)
    df$y <- log(df$y)
  }

  g <- ggplot(data = df, aes(x = x, y = y))
  g <- g + geom_point()
  g <- g + xlab(xlab)
  g <- g + ylab(ylab)
  g
}

plot.densities <- function(x, y, x.name, y.name, use.log = FALSE, rm.outliers = FALSE) {
  df <- data.frame(x = c(x, y), name = c(rep(x.name, length(x)), rep(y.name, length(y))))
  if(use.log) {
    df$x <- log(df$x)
  }
  if(rm.outliers) {
    df$x <- remove.outliers(df$x)
  }
  g <- ggplot(data = df, aes(x = x, colour = name))
  g <- g + geom_density()
  g  
}

thresholds <- c(0, 0.05, 0.1, 0.15, 0.2)

for(t in thresholds) {

  cat(paste0("Recalculating DSS for t = ", t, "\n"))

  ## Update as add/remove data sets

  ## Calculate and store DSS using a threshold of t for LL.4 and L.4 FIMM fits
  cat(paste0("Computing LL.4 FIMM DSS (t = ", t, ")\n"))

  ## DSS code will use min.conc and max.conc columns to establish integration bounds--these were set above.
  ll4.fimm.dss <- compute.all.dss(ll4.fimm.common, t = t, fct = "LL.4") 
  cat(paste0("Computing L.4 FIMM DSS (t = ", t, ")\n"))
  l4.fimm.dss <- compute.all.dss(l4.fimm.common, t = t, fct = "L.4") 

  common.fimm.cols <- c("SCREEN_ID", "PATIENT_ID", "DRUG_ID", "min.conc", "max.conc", "all.concs.nM", "uniq.concs.nM")
  fimm.dss <- merge(ll4.fimm.dss, l4.fimm.dss, by = common.fimm.cols, suffixes = c(".ll4", ".l4"), all = TRUE)
  colnames(fimm.dss)[colnames(fimm.dss) == "min.conc"] <- "min.conc.nM"
  colnames(fimm.dss)[colnames(fimm.dss) == "max.conc"] <- "max.conc.nM"

  ## Store the FIMM fits in FIMM Data/Processed Data (syn8270577)
  parentId <- "syn8270577"

  path <- paste0("fimm.fo.common.drugs.dss.t", t, ".tsv")
  write.table(file=path, fimm.dss, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
  f <- File(path, parentId = parentId, synapseStore = TRUE)
  synStore(f, executed = NULL, forceVersion = FALSE)

  ## Calculate and store DSS using a threshold of t for LL.4 and L.4 OHSU fits
  cat(paste0("Computing LL.4 OHSU DSS (t = ", t, ")\n"))
  ll4.ohsu.dss <- compute.all.dss(ll4.ohsu.common, t = t, fct = "LL.4") 

  cat(paste0("Computing L.4 OHSU DSS (t = ", t, ")\n"))
  l4.ohsu.dss <- compute.all.dss(l4.ohsu.common, t = t, fct = "L.4") 

  common.ohsu.cols <- c("patient_id", "inhibitor", "lab_id", "replicant", "time_of_read", "min.conc", "max.conc", "all.concs.nM", "uniq.concs.nM")
  ohsu.dss <- merge(ll4.ohsu.dss, l4.ohsu.dss, by = common.ohsu.cols, suffixes = c(".ll4", ".l4"), all = TRUE)
  colnames(ohsu.dss)[colnames(ohsu.dss) == "min.conc"] <- "min.conc.nM"
  colnames(ohsu.dss)[colnames(ohsu.dss) == "max.conc"] <- "max.conc.nM"

  ## Store the OHSU fits and expression in BEAT AML Data/Processed Data (syn10083332)
  parentId <- "syn10083332"

  path <- paste0("ohsu.fo.common.drugs.dss.t", t, ".tsv")
  write.table(file=path, ohsu.dss, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
  f <- File(path, parentId = parentId, synapseStore = TRUE)
  synStore(f, executed = NULL, forceVersion = FALSE)

  gc()
  save.image(".Rdata")

  datasets <- c("FIMM", "OHSU")

  model.suffices <- c("l4", "ll4")
  model.names <- c("L.4", "LL.4")

  if(t == thresholds[1]) {
    ## Note that IC50 is the e parameter/column
    
    ## Plot L.4 vs LL.4 for AUC and (log) IC50
    metrics <- c("auc", "ic50")
    metric.names <- c("AUC", "Log IC50")
    metric.use.log <- c(FALSE, TRUE)
    for(dataset in datasets) {
      dss <- fimm.dss
      if(dataset == "OHSU") { dss <- ohsu.dss }
      for(m in 1:length(metrics)) {
cat("Plotting l4 vs ll4\n")
        g <- plot.l4.vs.ll4(dss, metrics[m], metric.names[m], use.log = metric.use.log[m], rm.outliers = TRUE)
print(paste0(dataset, "-", metric.names[m], "-l4-vs-ll4.pdf"))
pdf(paste0(dataset, "-", metric.names[m], "-l4-vs-ll4.pdf"))
print(g)
d <- dev.off()
cat("Done plotting l4 vs ll4\n")
      }
    }

    ## For AUC and (log) IC50
    ##   For L.4 and LL.4
    ##      Plot FIMM and OHSU densities
    for(m in 1:length(metrics)) {
      for(o in 1:length(model.suffices)) {
        col <- paste0(metrics[m], ".", model.suffices[o])
cat("Plotting densities\n")
        g <- plot.densities(fimm.dss[,col,drop=T], ohsu.dss[,col,drop=T], "FIMM", "OHSU", use.log = metric.use.log[m], rm.outliers = TRUE)
pdf(paste0(col, "-fimm-vs-ohsu.pdf"))
print(g)
d <- dev.off()
cat("Done plotting densities\n")
      }
    }
  }

  ## Plot L.4 vs LL.4 for all DSS
  metrics <- c("dss1", "dss2", "dss3")
  metric.names <- c("DSS1", "DSS2", "DSS3")
  metric.use.log <- c(FALSE, FALSE, FALSE)
  for(dataset in datasets) {
    dss <- fimm.dss
    if(dataset == "OHSU") { dss <- ohsu.dss }
    for(m in 1:length(metrics)) {
cat("Plotting l4 vs l44 dss\n")
      g <- plot.l4.vs.ll4(dss, metrics[m], metric.names[m], use.log = metric.use.log[m], rm.outliers = TRUE)
cat("Done plotting l4 vs l44 dss\n")
    }
  }

  ## For all DSS
  ##   For L.4 and LL.4
  ##      Plot FIMM and OHSU densities
  for(m in 1:length(metrics)) {
    for(o in 1:length(model.suffices)) {
      col <- paste0(metrics[m], ".", model.suffices[o])
cat("Plotting dss densities\n")
      g <- plot.densities(fimm.dss[,col,drop=T], ohsu.dss[,col,drop=T], "FIMM", "OHSU", use.log = metric.use.log[m], rm.outliers = TRUE)
cat("Done plotting dss densities\n")
    }
  }

}


