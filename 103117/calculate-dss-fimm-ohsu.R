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

## Begin setup (this will need to be changed as the data sets in the intersection change)

## Data sets to be analyzed:
data.sets <- c("ohsu", "fimm")
data.sets.to.analyze <- data.sets

## The drug mapping file lists the drugs in common across all data sets with the drug.map.drug.id.cols
## providing the identifier for each drug for each respective data set.
## The drug mapping file may be created by a script such as 
## harmonize-fimm-ohsu-compounds.R
drug.mapping.synId <- "syn11361110" ## fimm-ohsu-drug-map.tsv
drug.map.drug.id.cols <- c("OHSU_DRUG_NAME", "FIMM_Batch_ID_Drug")
names(drug.map.drug.id.cols) <- data.sets

## These were created by the script 070617/fit-ohsu-and-fimm-drug-response-data.R
data.set.drug.fit.synIds <- c("syn10083494", "syn10083488")
names(data.set.drug.fit.synIds) <- data.sets
data.set.drug.fit.files <- c("ohsu.dss.t0.tsv", "fimm.dss.t0.tsv")

## The columns within each respective data set that identifies the drug (in the name space of that data set)
data.set.drug.id.cols <- list("ohsu" = "inhibitor", "fimm" = "DRUG_ID")

## The columns within each respective data set that identifies a screen (i.e., a particular replicate of a patient/cell line x drug)
screen.id.cols <- list("ohsu" = c("patient_id", "inhibitor", "lab_id", "SeqID", "replicant", "time_of_read"),
                       "fimm" = c("SCREEN_ID", "PATIENT_ID", "DRUG_ID"))

## Batch-corrected log2 cpm (or rma) expression files (NB: these exclude 2 outliers in OHSU)
data.set.expr.synIds <- c("syn11362256", "syn11362257")
names(data.set.expr.synIds) <- data.sets
data.set.expr.files <- c("ohsu-expr-fimm-ohsu-outlier1-combat.tsv", "fimm-expr-fimm-ohsu-outlier1-combat.tsv")

## Columns whose values are not dependent on the (L.4 vs LL.4) fit.  
merge.cols <- list(
  "ohsu" = c(screen.id.cols[["ohsu"]], "min.conc", "max.conc", "all.concs.nM", "uniq.concs.nM"),
  "fimm" = c(screen.id.cols[["fimm"]], "min.conc", "max.conc", "all.concs.nM", "uniq.concs.nM"))

## The synapseIds of folders in which to store the fits for each respective data set.
## FIMM_BEAT AML Collaboration/Files/Analyses/fimm-ohsu/expr-outlier-1 ("syn11362287")
fit.folder.synIds <- list("ohsu" = "syn11362287", "fimm" = "syn11362287")
## A string that will be included in the fit file name.
## Here "foc-aml-s-aml" = "fimm-ohsu-ctrp-aml-s-aml"
fit.file.prefix <- "fimm-ohsu-outlier1"

patient.id.cols <- list("ohsu" = "patient_id", "fimm" = "PATIENT_ID")
patient.subsets <- list("ohsu" = NULL, "fimm" = NULL)

expr.patient.id.cols <- list("ohsu" = "SeqID", "fimm" = "SCREEN_ID")

## End setup

cat(paste0("Calculate AUC and DSS values based on ", length(data.sets), "-way intersection of concentrations\n"))
cat(paste0("for each drug shared between the data sets: ", paste(data.sets, collapse=","), "\n"))

## Load in the fit files
fits <- llply(data.set.drug.fit.synIds, .parallel = FALSE,
              .fun = function(synId) {
                       obj <- synGet(synId, downloadFile = TRUE)
                       file <- getFileLocation(obj)
                       read.table(file, sep="\t", header=TRUE, as.is=TRUE, comment.char="", quote="\"")
              })

## Filter some of the fits
## Filter FIMM results to include only those drugs with at least 5 uniq concentration points
for(ds in names(fits)) {
  switch(ds,
    "fimm" = {
      min.conc.pts <- 5
      keep <- sapply(fits[[ds]][, c("uniq.concs.nM")], function(str) length(unique(unlist(strsplit(str, split=",")))) >= min.conc.pts)
      cat(paste0("Keeping FIMM screens with at least ", min.conc.pts, " unique concentrations\n"))
      print(table(keep))
      fits[[ds]] <- fits[[ds]][keep,]
    }, 
    "ohsu" = {
      min.conc.pts <- 5
      max.conc.pts <- 7
      keep <- unlist(apply(fits[[ds]][ ,c("all.concs.nM", "uniq.concs.nM")], 1, function(row) {
                 ( length(unique(unlist(strsplit(row[1], split=",")))) <= max.conc.pts ) && ( length(unique(unlist(strsplit(row[2], split=",")))) >= min.conc.pts )
              }))
      cat(paste0("Keeping OHSU screens with at least ", min.conc.pts, " unique concentrations and at most ", max.conc.pts, " total concentrations\n"))
      print(table(keep))
      fits[[ds]] <- fits[[ds]][keep,]
    },
    {
    })
}

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
for(i in 1:length(fits)) {
  fits[[i]] <- fits[[i]][fits[[i]][, data.set.drug.id.cols[[i]]] %in% drug.map[, drug.map.drug.id.cols[i]], ]
}

## Drop any precomputed DSS files (if they exist) in the the fit files
## and only keep fits if at least one of L.4 or LL.4 converged
for(i in 1:length(fits)) {
  flag <- grepl(x = colnames(fits[[i]]), pattern = "dss")
  fits[[i]] <- fits[[i]][, !flag]
  fits[[i]] <- subset(fits[[i]], (converged.l4 == 1) | (converged.ll4 == 1))
}

decimalplaces <- function(x) {
    if ((x %% 1) != 0) {
##        nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed=TRUE)[[1]][[2]])
        nchar(strsplit(as.character(x), ".", fixed=TRUE)[[1]][[2]])
    } else {
        return(0)
    }
}

screen.ranges <- list()
drug.ranges <- list()
for(i in 1:length(fits)) {
  ## Calculate the concentration range for each screen
  screen.ranges[[i]] <- ddply(fits[[i]], .variables = screen.id.cols[[i]],
                              .fun = function(df) {
                                       v <- c(min(df$min.conc.nM), max(df$max.conc.nM), length(unique(unlist(strsplit(df$uniq.concs.nM[1], split=",")))))
                                       names(v) <- c("min.conc.nM", "max.conc.nM", "num.concs")
                                       v 
                              })
  names(screen.ranges)[i] <- names(fits)[i]
  options(scipen=999) ## disable scientific notation so we can compute number of decimals
  ## Round off insignificant differences between concentrations for a given drug.
  ## For any _min_ concentrations that are the same (to within rounding error), set to the max of the rounding-equivalent concentrations
  ## For any _max_ concentrations that are the same (to within rounding error), set to the min of the rounding-equivalent concentrations
  screen.ranges[[i]] <- ddply(screen.ranges[[i]], .variables = c(data.set.drug.id.cols[[i]]),
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
  drug.ranges[[i]] <- ddply(screen.ranges[[i]], .variables = c(data.set.drug.id.cols[[i]], "min.conc.nM", "max.conc.nM"),
                            .fun = function(df) {
                                     vec <- c(df[1, data.set.drug.id.cols[[i]]], df$num.concs[1], df$min.conc.nM[1], df$max.conc.nM[1], freq = nrow(df))
                                     names(vec) <- c(data.set.drug.id.cols[[i]], "num.concs", "min.conc.nM", "max.conc.nM", "freq")
                                     vec
                                   })
  names(drug.ranges)[i] <- names(fits)[i]
}

## Define a single range for each drug (in each data set)
## If any range occurs across more than 80% of the screens, use that range.
## Otherwise, drop the 20% most narrow ranges and use the most restrictive remaining range.                                 
## Define range by the number of concentration points (not the absolute range)
final.screen.ranges <- list()
min.frac.for.range <- 0.8
frac.to.discard <- 0.2
for(i in 1:length(fits)) {
  final.screen.ranges[[i]] <- ddply(drug.ranges[[i]], .variables = c(data.set.drug.id.cols[[i]]),
                              .fun = function(df) {
                                       df$min.conc.nM <- as.numeric(df$min.conc.nM)
                                       df$max.conc.nM <- as.numeric(df$max.conc.nM)
                                       df$freq <- as.numeric(df$freq)
                                       df$rel.freq <- df$freq / sum(df$freq)
                                       flag <- ( df$rel.freq == max(df$rel.freq) ) & ( df$rel.freq >= min.frac.for.range )
                                       ret <- NULL
                                       if(any(flag)) {
                                         indx <- which(flag)[1]
                                         ret <- df[indx, c(data.set.drug.id.cols[[i]], "min.conc.nM", "max.conc.nM")]
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
                                         ret <- df[1, data.set.drug.id.cols[[i]], drop=F]
                                         ret$min.conc.nM <- max(df$min.conc.nM)
                                         ret$max.conc.nM <- min(df$max.conc.nM)
                                       }
                                       return(ret)
                              })
}
names(final.screen.ranges) <- names(fits)

## Merge the ranges for all data sets together
shared.drug.ranges <- drug.map
for(i in 1:length(fits)) {
  tmp <- final.screen.ranges[[i]][, c(data.set.drug.id.cols[[i]], "min.conc.nM", "max.conc.nM")]
  colnames(tmp) <- c(data.set.drug.id.cols[[i]], paste0(data.sets[i], ".min.conc.nM"), paste0(data.sets[i], ".max.conc.nM"))
  shared.drug.ranges <- merge(shared.drug.ranges, tmp, by.y = data.set.drug.id.cols[[i]], by.x = drug.map.drug.id.cols[i], all = FALSE)
}

min.conc.cols <- unlist(lapply(data.sets, function(ds) paste0(ds, ".min.conc.nM")))
max.conc.cols <- unlist(lapply(data.sets, function(ds) paste0(ds, ".max.conc.nM")))
## DSS code will use min.conc and max.conc columns
shared.drug.ranges$min.conc <- unlist(apply(shared.drug.ranges[, min.conc.cols], 1, function(row) max(as.numeric(row))))
shared.drug.ranges$max.conc <- unlist(apply(shared.drug.ranges[, max.conc.cols], 1, function(row) min(as.numeric(row))))

## Merge the processed drug response tables with the range tables and exclude any screens that do 
## not completely cover the min/max range
for(i in 1:length(fits)) {
  tmp <- shared.drug.ranges[, c(drug.map.drug.id.cols[i], "min.conc", "max.conc")]
  ## Exclude "min.conc.nM" and "max.conc.nM" from the fits, to avoid confusion
  fits[[i]] <- fits[[i]][, !(colnames(fits[[i]]) %in% c("min.conc.nM", "max.conc.nM"))]
  fits[[i]] <- merge(fits[[i]], tmp, by.y = drug.map.drug.id.cols[i], by.x = data.set.drug.id.cols[[i]], all = FALSE)
  keep <- unlist(apply(fits[[i]][, c("uniq.concs.nM", "min.conc", "max.conc")], 1,   
                       function(row) {
                         concs <- as.numeric(unlist(strsplit(row[1], split=",")))
                         min.conc.nM <- as.numeric(row[2])
                         max.conc.nM <- as.numeric(row[3])
                         (min.conc.nM < max.conc.nM) & (min(concs) < max(concs)) & (min(concs) <= min.conc.nM) & (max(concs) >= max.conc.nM)
                       }))
  fits[[i]] <- fits[[i]][keep, ]
}

## Recalculate DSS using both L.4 and LL.4 fits (for FIMM and OHSU) and using the one CTRP fit for varying thresholds t

source("../common/dss.R")

thresholds <- c(0, 0.05, 0.1, 0.15, 0.2)
thresholds <- c(0)

for(data.set in data.sets) {

  for(t in thresholds) {

    fcts <- c("LL.4", "L.4")
    exclude.patterns <- c("\\.l4", "\\.ll4")
    include.patterns <- c(".ll4", ".l4")
    fits.dss <- list()
    for(k in 1:length(fcts)) {

      fct <- fcts[k]
      cat(paste0("Recalculating DSS for t = ", t, " based on ", fct, " fits for data set = ", data.set, "\n"))
 
      fit <- fits[[data.set]]
      fit <- fit[, !grepl(x = colnames(fit), pattern=exclude.patterns[k])]
      colnames(fit) <- gsub(x = colnames(fit), pattern=include.patterns[k], replacement="")
      print(head(fit))
      ## DSS code expects that the fit data.frame has columns:
      ## converged, min.conc, max.conc, and fit parameters b, c, d, e.
      fits.dss[[k]] <- compute.all.dss(fit, t = t, fct = fct) 

    } ## for(k in 1:length(fcts)) 

    common.cols <- merge.cols[[data.set]]
    dss.tbl <- merge(fits.dss[[1]], fits.dss[[2]], by = common.cols, suffixes = include.patterns[1:2], all = TRUE)
    if(length(fits.dss) > 2) {
      for(j in 3:length(fits.dss)) {
        dss.tbl <- merge(dss.tbl, fits.dss[[j]], by = common.cols, suffixes = c("", include.patterns[j]), all = TRUE)
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