suppressPackageStartupMessages(library("foreach"))
suppressPackageStartupMessages(library("parallel"))
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

## Read in the setup.
source("fimm-ohsu-setup-112017.R")

cat(paste0("Calculate AUC and DSS values based on ", length(data.sets), "-way intersection of concentrations\n"))
cat(paste0("for each drug shared between the data sets: ", paste(data.sets, collapse=","), "\n"))

## Load in the fit files
fits <- llply(data.set.drug.fit.synIds, .parallel = FALSE,
              .fun = function(synId) {
                       obj <- synGet(synId, downloadFile = TRUE)
                       file <- getFileLocation(obj)
                       read.table(file, sep="\t", header=TRUE, as.is=TRUE, comment.char="", quote="\"")
              })

for(ds in names(fits)) {
  switch(ds,
    "ohsu" = {
      ## Add the SeqID to the OHSU drug data.
      ## Read in the rna-seq summary table, which provides a map from the lab_id's used
      ## in the drug response data to the SeqID's used in the expression data.
      synId <- "syn10083817"
      obj <- synGet(synId, downloadFile = TRUE)
      file <- getFileLocation(obj)
      ohsu.rnaseq.sample.summary <- read.table(file, header=TRUE, sep="\t")
      fits[[ds]] <- merge(fits[[ds]], ohsu.rnaseq.sample.summary, by.x = "lab_id", by.y = "Original_LabID")
    },
    "fimm" = {
    },
    "ctrp" = {
      fits[[ds]] <- merge(fits[[ds]], ctrp.cell.line.metadata[, c("master_ccl_id", "ccl_name")])
    },
    "ntap" = {
    },
    "sanger" = {
    },
    {
      stop("Unrecognized data set\n")
    })
}


## Load in the expr files
exprs <- llply(data.set.expr.synIds, .parallel = FALSE,
              .fun = function(synId) {
                       obj <- synGet(synId, downloadFile = TRUE)
                       file <- getFileLocation(obj)
                       read.table(file, sep="\t", header=TRUE, as.is=TRUE, comment.char="", quote="\"")
              })

for(nm in names(exprs)) {
  switch(nm,
    "ohsu" = {
      ## Convert the column names from X20.00347 to 20-00347
      colnames(exprs[[nm]]) <- gsub("X(.*)\\.(.*)", "\\1-\\2", colnames(exprs[[nm]]))
    },
    "fimm" = {
    },
    "ctrp" = {
    },
    "ntap" = {
    },
    "sanger" = {
      colnames(exprs[[nm]]) <- gsub("X(.*)", "\\1", colnames(exprs[[nm]]))
    },
    {
      stop("Unrecognized data set\n")
    })
}

## Restrict to samples of interest (e.g., disease subset)
for(ds in names(fits)) {
  if(!is.null(patient.subsets[[ds]])) {
    flag <- fits[[ds]][, patient.id.cols[[ds]]] %in% patient.subsets[[ds]]
    orig.sz <- nrow(fits[[ds]])
    orig.num.pts <- length(unique(fits[[ds]][, patient.id.cols[[ds]]]))
    fits[[ds]] <- fits[[ds]][flag, ]
    num.pts <- length(unique(fits[[ds]][, patient.id.cols[[ds]]]))
    cat(paste0("Limiting ", ds, " data set to ", nrow(fits[[ds]]), " fits (across ", orig.num.pts, " patients) from ", orig.sz, " (across ", num.pts, " patients).\n"))
  }
}

## Restrict to samples with expression data
for(ds in names(fits)) {
  flag <- fits[[ds]][, expr.patient.id.cols[[ds]]] %in% colnames(exprs[[ds]])
  orig.sz <- nrow(fits[[ds]])
  orig.num.pts <- length(unique(fits[[ds]][, expr.patient.id.cols[[ds]]]))
  fits[[ds]] <- fits[[ds]][flag, ]
  num.pts <- length(unique(fits[[ds]][, expr.patient.id.cols[[ds]]]))
  cat(paste0("Limiting ", ds, " data set to ", nrow(fits[[ds]]), " fits with expression (across ", orig.num.pts, " patients) from ", orig.sz, " (across ", num.pts, " patients).\n"))
}

## Load the map of drugs between data sets
obj <- synGet(drug.mapping.synId, downloadFile = TRUE)
drug.map <- read.table(getFileLocation(obj), sep="\t", header=TRUE, quote="\"", comment.char="")

my.dup <- function(x) duplicated(x, fromLast = TRUE) | duplicated(x, fromLast = FALSE)
for(col in drug.map.drug.id.cols) {
  flag <- !my.dup(drug.map[, col])
  drug.map <- drug.map[flag, ]
}

## Restrict the fits to those drugs in common across the data sets
for(ds in names(fits)) {
  fits[[ds]] <- fits[[ds]][fits[[ds]][, data.set.drug.id.cols[[ds]]] %in% drug.map[, drug.map.drug.id.cols[[ds]]], ]
}

## Drop any precomputed DSS files (if they exist) in the the fit files
## and only keep fits if at least one of L.4 or LL.4 converged
for(ds in names(fits)) {
  flag <- grepl(x = colnames(fits[[ds]]), pattern = "dss")
  fits[[ds]] <- fits[[ds]][, !flag]
  fits[[ds]] <- subset(fits[[ds]], (converged.l4 == 1) | (converged.ll4 == 1))
}

decimalplaces <- function(x) {
    if ((x %% 1) != 0) {
##        nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed=TRUE)[[1]][[2]])
        nchar(strsplit(as.character(x), ".", fixed=TRUE)[[1]][[2]])
    } else {
        return(0)
    }
}

## Calculate drug ranges for each data set and each fit type (L.4 or LL.4)
## Drug ranges may differ between fits since we do different filtering (e.g., based on GOF)
drug.ranges <- list()
fcts <- c("LL.4", "L.4")
exclude.patterns <- list("LL.4" = "\\.l4", "L.4" = "\\.ll4")
include.patterns <- list("LL.4" = ".ll4", "L.4" = ".l4")
fct.postfixes <- list("LL.4" = ".ll4", "L.4" = ".l4")

for(ds in names(fits)) {
  drug.ranges[[ds]] <- list()
  for(fct in fcts) {
    ## Calculate the concentration range for each screen
    ## Do not consider filtered fits
    exclude.cols <- colnames(fits[[ds]])[grepl(pattern=paste0("exclude", fct.postfixes[[fct]]), x=colnames(fits[[ds]]))]
    any.excluded <- unlist(apply(fits[[ds]][, exclude.cols], 1, function(row) any(row)))
    screen.ranges <- ddply(fits[[ds]][!any.excluded, ], .variables = screen.id.cols[[ds]],
                                .fun = function(df) {
                                         v <- c(min(df$min.conc.nM), max(df$max.conc.nM), length(unique(unlist(strsplit(df$uniq.concs.nM[1], split=",")))))
                                         names(v) <- c("min.conc.nM", "max.conc.nM", "num.concs")
                                         v 
                                })

    options(scipen=999) ## disable scientific notation so we can compute number of decimals
    ## Round off insignificant differences between concentrations for a given drug.
    ## For any _min_ concentrations that are the same (to within rounding error), set to the max of the rounding-equivalent concentrations
    ## For any _max_ concentrations that are the same (to within rounding error), set to the min of the rounding-equivalent concentrations
    screen.ranges <- ddply(screen.ranges, .variables = c(data.set.drug.id.cols[[ds]]),
                                .fun = function(df) {
                                         df$min.conc.nM <- as.numeric(df$min.conc.nM)
                                         df$max.conc.nM <- as.numeric(df$max.conc.nM)
                                         for(j in 1:nrow(df)) {
                                           num.decimal <- decimalplaces(df$min.conc.nM[j])
                                           flag <- round(df$min.conc.nM[j], digits = num.decimal) == round(df$min.conc.nM, digits = num.decimal)
                                           df$min.conc.nM[flag] <- max(df$min.conc.nM[flag])
                                           num.decimal <- decimalplaces(df$max.conc.nM[j])
                                           flag <- round(df$max.conc.nM[j], digits = num.decimal) == round(df$max.conc.nM, digits = num.decimal)
                                           df$max.conc.nM[flag] <- min(df$max.conc.nM[flag])
                                         } 
                                         df        
                                })
    ## Calculate the concentration range for each drug
    drug.ranges[[ds]][[fct]] <- ddply(screen.ranges, .variables = c(data.set.drug.id.cols[[ds]], "min.conc.nM", "max.conc.nM"),
                              .fun = function(df) {
                                       vec <- c(df[1, data.set.drug.id.cols[[ds]]], df$num.concs[1], df$min.conc.nM[1], df$max.conc.nM[1], freq = nrow(df))
                                       names(vec) <- c(data.set.drug.id.cols[[ds]], "num.concs", "min.conc.nM", "max.conc.nM", "freq")
                                       vec
                                     })
  }
}

## Define a single range for each drug (in each data set) and each fit type (LL.4 or L.4)
## If any range occurs across more than 80% of the screens, use that range.
## Otherwise, drop the 20% most narrow ranges and use the most restrictive remaining range.                                 
## Define range by the number of concentration points (not the absolute range)
final.screen.ranges <- list()
min.frac.for.range <- 0.8
frac.to.discard <- 0.2
for(ds in names(fits)) {
  final.screen.ranges[[ds]] <- list()
  for(fct in fcts) {
    final.screen.ranges[[ds]][[fct]] <- ddply(drug.ranges[[ds]][[fct]], .variables = c(data.set.drug.id.cols[[ds]]),
                                .fun = function(df) {
                                         df$min.conc.nM <- as.numeric(df$min.conc.nM)
                                         df$max.conc.nM <- as.numeric(df$max.conc.nM)
                                         df$freq <- as.numeric(df$freq)
                                         df$rel.freq <- df$freq / sum(df$freq)
                                         flag <- ( df$rel.freq == max(df$rel.freq) ) & ( df$rel.freq >= min.frac.for.range )
                                         ret <- NULL
                                         if(any(flag)) {
                                           indx <- which(flag)[1]
                                           ret <- df[indx, c(data.set.drug.id.cols[[ds]], "min.conc.nM", "max.conc.nM")]
                                         } else {
                                           ## Discard the smallest frac.to.discard ranges and then define the restrictive range over what's left
  ##                                         df$range <- abs(df$max.conc.nM - df$min.conc.nM)
                                           df$range <- as.numeric(df$num.concs)
                                           df <- df[order(df$range, df$freq, decreasing=FALSE),]
                                           cum.freq <- 0
                                           first.indx.to.keep <- 1
                                           for(j in 1:(nrow(df)-1)) {
                                             if( ( cum.freq + df$rel.freq[j] ) > frac.to.discard ) { break }
                                             cum.freq <- cum.freq + df$rel.freq[j]
                                             first.indx.to.keep <- j + 1
                                           } 
                                           df <- df[first.indx.to.keep:nrow(df),,drop=F]
                                           ret <- df[1, data.set.drug.id.cols[[ds]], drop=F]
                                           ret$min.conc.nM <- max(df$min.conc.nM)
                                           ret$max.conc.nM <- min(df$max.conc.nM)
                                         }
                                         return(ret)
                                })
  }
}


## Merge the ranges for all data sets together
shared.drug.ranges <- list()
for(fct in fcts) {
  shared.drug.ranges[[fct]] <- drug.map
  for(ds in names(fits)) {
    tmp <- final.screen.ranges[[ds]][[fct]][, c(data.set.drug.id.cols[[ds]], "min.conc.nM", "max.conc.nM")]
    colnames(tmp) <- c(data.set.drug.id.cols[[ds]], paste0(data.sets[[ds]], ".min.conc.nM"), paste0(data.sets[[ds]], ".max.conc.nM"))
    cat(paste0("Merging by ", data.set.drug.id.cols[[ds]], " and ", drug.map.drug.id.cols[[ds]], "\n"))
    shared.drug.ranges[[fct]] <- merge(shared.drug.ranges[[fct]], tmp, by.y = data.set.drug.id.cols[[ds]], by.x = drug.map.drug.id.cols[[ds]], all = FALSE)
  }
}

min.conc.cols <- unlist(lapply(data.sets, function(ds) paste0(ds, ".min.conc.nM")))
max.conc.cols <- unlist(lapply(data.sets, function(ds) paste0(ds, ".max.conc.nM")))
## DSS code will use min.conc and max.conc columns
for(fct in fcts) {
  shared.drug.ranges[[fct]]$min.conc <- unlist(apply(shared.drug.ranges[[fct]][, min.conc.cols], 1, function(row) max(as.numeric(row))))
  shared.drug.ranges[[fct]]$max.conc <- unlist(apply(shared.drug.ranges[[fct]][, max.conc.cols], 1, function(row) min(as.numeric(row))))
}

## Merge the processed drug response tables with the range tables and exclude any screens that do 
## not completely cover the min/max range
fits.with.concs <- list()
for(ds in names(fits)) {
  fits.with.concs[[ds]] <- list()
  for(fct in fcts) {
    tmp <- shared.drug.ranges[[fct]][, c(drug.map.drug.id.cols[[ds]], "min.conc", "max.conc")]
    ## Exclude "min.conc.nM" and "max.conc.nM" from the fits, to avoid confusion
    fits.with.concs[[ds]][[fct]] <- fits[[ds]][, !(colnames(fits[[ds]]) %in% c("min.conc.nM", "max.conc.nM"))]
    cat(paste0("Merging by ", drug.map.drug.id.cols[[ds]], " and ", data.set.drug.id.cols[[ds]], "\n"))
    fits.with.concs[[ds]][[fct]] <- merge(fits.with.concs[[ds]][[fct]], tmp, by.y = drug.map.drug.id.cols[[ds]], by.x = data.set.drug.id.cols[[ds]], all = FALSE)
    keep <- unlist(apply(fits.with.concs[[ds]][[fct]][, c("uniq.concs.nM", "min.conc", "max.conc")], 1,   
                         function(row) {
                           concs <- as.numeric(unlist(strsplit(row[1], split=",")))
                           min.conc.nM <- as.numeric(row[2])
                           max.conc.nM <- as.numeric(row[3])
                           (min.conc.nM < max.conc.nM) & (min(concs) < max(concs)) & (min(concs) <= min.conc.nM) & (max(concs) >= max.conc.nM)
                         }))
    fits.with.concs[[ds]][[fct]] <- fits.with.concs[[ds]][[fct]][keep, ]
  }
}

## Recalculate DSS using both L.4 and LL.4 fits (for FIMM and OHSU) and using the one CTRP fit for varying thresholds t

source("../common/dss.R")

thresholds <- c(0, 0.05, 0.1, 0.15, 0.2)
thresholds <- c(0)

for(data.set in data.sets) {

  for(t in thresholds) {

    fits.dss <- list()
    for(fct in fcts) {

      cat(paste0("Recalculating DSS for t = ", t, " based on ", fct, " fits for data set = ", data.set, "\n"))
 
      fit <- fits.with.concs[[data.set]][[fct]]
      fit <- fit[, !grepl(x = colnames(fit), pattern=exclude.patterns[[fct]])]
      colnames(fit) <- gsub(x = colnames(fit), pattern=include.patterns[[fct]], replacement="")
      print(head(fit))
      ## DSS code expects that the fit data.frame has columns:
      ## converged, min.conc, max.conc, and fit parameters b, c, d, e.
      fits.dss[[fct]] <- compute.all.dss(fit, t = t, fct = fct) 

    } ## for(k in 1:length(fcts)) 

    common.cols <- merge.cols[[data.set]]
    fit.names <- names(fits.dss)
    cat(paste0("Merging by ", paste(common.cols, collapse=","), "\n"))
    dss.tbl <- merge(fits.dss[[fit.names[1]]], fits.dss[[fit.names[2]]], by = common.cols, suffixes = c(include.patterns[[fit.names[1]]], include.patterns[[fit.names[2]]]), all = TRUE)
    if(length(fit.names) > 2) {
      for(fit.name in fit.names[3:length(fit.names)]) {
        dss.tbl <- merge(dss.tbl, fits.dss[[fit.name]], by = common.cols, suffixes = c("", include.patterns[[fit.name]]), all = TRUE)
      }
    }
    colnames(dss.tbl)[colnames(dss.tbl) == "min.conc"] <- "min.conc.nM"
    colnames(dss.tbl)[colnames(dss.tbl) == "max.conc"] <- "max.conc.nM"

    ## Store the fits in synapse
    parentId <- fit.folder.synIds[[data.set]]

    path <- paste0(data.set, "-", fit.file.prefix, "-dss-t", t, ".tsv")
    write.table(file=path, dss.tbl, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
    f <- File(path, parentId = parentId, synapseStore = TRUE)
    synStore(f, executed = NULL, forceVersion = FALSE)

  } ## for(t in thresholds) 

} # for(data.set in data.sets) 

cat("Successfully calculated DSS.\n")