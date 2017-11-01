bioc.packages <- c("nplr", "openxlsx", "drc", "minpack.lm", "maftools", "vcd", "binom", "DMwR", "WGCNA", "VennDiagram")

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

## Set the download and output directory
download.path <- "../input/"
if (!file.exists(download.path)) {
  dir.create(download.path)
}

output.path <- "output/"
if (!file.exists(output.path)) {
  dir.create(output.path)
}

beat.aml.cache <- "../input/beat-aml-files/"

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
## suppressPackageStartupMessages(library("GGally"))
suppressPackageStartupMessages(library("VennDiagram"))
suppressPackageStartupMessages(library(ReporteRs))
suppressPackageStartupMessages(library(rJava))

library(ggbeeswarm)
library(plyr)
library(stringr)
library(synapseClient)

synapseLogin()

## Set the java heap to 8GB
options(java.parameters = "-Xmx8000m")

## Define a java garbage collection function
jgc <- function()
{
  gc()
  .jcall("java/lang/System", method = "gc")
} 

## Load the FIMM and OHSU fits and recalculate the DSS values.
## Note that t0 and t10 have the same L.4 and LL.4 fits and include the same
## rows.  So, just load in the t0 and recompute for all values of t.

## path <- "fimm.dss.t0.tsv"
synId <- "syn10083488"
obj <- synGet(synId, downloadFile = TRUE)
file <- getFileLocation(obj)
fimm.dss.t0.orig <- as.data.frame(fread(file))

## Filter FMM results to include only those drugs with at least 5 uniq concentration points
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

## Load the map of drugs shared between OHSU and FIMM (FIMM_OHSU_Drug_Concentrations.xlsx)
synId <- "syn10163669"
obj <- synGet(synId, downloadFile = TRUE)
ohsu.fimm.drugs <- openxlsx:::read.xlsx(getFileLocation(obj), sheet=1)

## Folder for FIMM Data/Processed Data (syn8270577)
parentId <- "syn8270577"

## Folder for BEAT AML Data/Processed Data (syn10083332)
parentId <- "syn10083332"

## Confirm that all of the drugs in the shared drug annotation file are
## in the processed drug response file.
missing.drugs <- ohsu.fimm.drugs[!(ohsu.fimm.drugs$FIMM_Batch_ID_Drug %in% fimm.dss.t0.orig$DRUG_ID),]
if(nrow(missing.drugs) == 0) {
  cat("As expected, all drugs in shared file are in FIMM drug response data file\n")
} else {
  cat("Missing some drugs in shared file from FIMM drug response file\n")
  q(status = -1)
}

missing.drugs <- ohsu.fimm.drugs[!(ohsu.fimm.drugs$OHSU_DRUG_NAME %in% ohsu.dss.t0.orig$inhibitor),]
if(nrow(missing.drugs) == 0) {
  cat("As expected, all drugs in shared file are in OHSU drug response data file\n")
} else {
  cat("Missing some drugs in shared file from OHSU drug response file\n")
  q(status = -1)
}

my.dup <- function(x) {
  duplicated(x, fromLast=TRUE) | duplicated(x, fromLast=FALSE)
}

## Note (and drop) any drugs that have multiple entries (which would imply ambiguous mappings
## and/or inconsistent drug ranges)
dups <- my.dup(ohsu.fimm.drugs$OHSU_DRUG_NAME) | my.dup(ohsu.fimm.drugs$FIMM_Batch_ID_Drug)
if(any(dups)) {
  cat("Dropping ambiguous mappings between OHSU and FIMM drugs:\n")
  write.table(ohsu.fimm.drugs[dups,], sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
  ohsu.fimm.drugs <- ohsu.fimm.drugs[!dups,]
}

## Restrict the processed results to those shared between both data sets
## And drop all of the pre-computed DSS columns.
fimm.dss.common <- subset(fimm.dss.t0.orig, DRUG_ID %in% ohsu.fimm.drugs$FIMM_Batch_ID_Drug)
fimm.dss.common <- fimm.dss.common[, !(grepl(x=colnames(fimm.dss.common), pattern="dss"))]
rm(fimm.dss.t0.orig)

ohsu.dss.common <- subset(ohsu.dss.t0.orig, inhibitor %in% ohsu.fimm.drugs$OHSU_DRUG_NAME)
ohsu.dss.common <- ohsu.dss.common[, !(grepl(x=colnames(ohsu.dss.common), pattern="dss"))]
rm(ohsu.dss.t0.orig)

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

## Compute the drug ranges from the _raw_ OHSU data.  i.e., before excluding outliers, which
## could conceivably affect the range.
synId <- "syn8149180"
inh.obj <- synGet(synId, downloadFile=TRUE, downloadLocation = download.path)
ohsu.raw.dr <- fread(getFileLocation(inh.obj))
ohsu.raw.dr <- as.data.frame(ohsu.raw.dr)

## But limit to the patients and drugs we are analyzing
ohsu.raw.drug.ranges <- ddply(subset(ohsu.raw.dr, (patient_id %in% ohsu.dss.common$patient_id) & (inhibitor %in% ohsu.dss.common$inhibitor)), 
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
## are quite similar.
cutoff <- 0.01
cat(paste0("OHSU drugs with min or max concentrations differing by >= ", cutoff, "\n"))
d_ply(ohsu.raw.drug.ranges, .variables = c("inhibitor"),
                             .fun = function(df) {
                                      df$min.conc.nM <- as.numeric(df$min.conc.nM); df$max.conc.nM <- as.numeric(df$max.conc.nM)
                                      if( ( (max(df$min.conc.nM) - min(df$min.conc.nM)) >= cutoff ) || ( (max(df$max.conc.nM) - min(df$max.conc.nM)) >= cutoff ) ) {
                                         write.table(df, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
                                      }
                                    })

## Hence, take the most restrictive range (i.e., the intersection) for the OHSU drugs
ohsu.raw.drug.ranges <- ddply(ohsu.raw.drug.ranges, .variables = c("inhibitor"),
                              .fun = function(df) {
                                       v <- c(max(df$min.conc.nM), min(df$max.conc.nM))
                                       names(v) <- c("min.conc.nM", "max.conc.nM")
                                       v
                              })

## Do the same for FIMM
fimm.raw.drug.ranges <- ddply(subset(fimm.raw.dr, (SCREEN_ID %in% fimm.dss.common$SCREEN_ID) & (DRUG_ID %in% fimm.dss.common$DRUG_ID)), 
                                           .variables = c("PATIENT_ID", "SCREEN_ID", "DRUG_SET", "DRUG_ID"),
                                           .fun = function(df) { 
                                                    v <- c(min(df$CONCENTRATION), max(df$CONCENTRATION))
                                                    names(v) <- c("min.conc.nM", "max.conc.nM")
                                                    v 
                                                  })

## Now look at the range of one drug in particular, as generated from the raw data.  
## Note some go from 1.0 to 10,000 and some from 0.1 to 1,000
drug <- "FIMM023833"
cat(paste0("FIMM drug range for ", drug, " calculated from raw data on a per-patient basis\n"))
tail(subset(fimm.raw.drug.ranges, DRUG_ID == drug))

## Now look at the range of that same drug in the annotation table.
## Note that the min and max range for FIMM is 0.1 to 10,000--i.e., the union of the ranges
## across all patients tested for this drug.
cat(paste0("FIMM drug range for ", drug, " as annotated in FIMM_OHSU_Drug_Concentrations.xlsx\n"))
subset(ohsu.fimm.drugs, FIMM_Batch_ID_Drug == drug)

## Collapse duplicate ranges and drop patient ids
fimm.raw.drug.ranges <- ddply(fimm.raw.drug.ranges, .variables = c("DRUG_ID", "min.conc.nM", "max.conc.nM"),
                              .fun = function(df) {
                                       vec <- c(df$DRUG_ID[1], df$min.conc.nM[1], df$max.conc.nM[1], freq = nrow(df))
                                       names(vec) <- c("DRUG_ID", "min.conc.nM", "max.conc.nM", "freq")
                                       vec
                                 })


## Hence, take the most restrictive range (i.e., the intersection) for the OHSU drugs
fimm.raw.drug.ranges <- ddply(fimm.raw.drug.ranges, .variables = c("DRUG_ID"),
                              .fun = function(df) {
                                       v <- c(max(df$min.conc.nM), min(df$max.conc.nM))
                                       names(v) <- c("min.conc.nM", "max.conc.nM")
                                       v
                              })

if(FALSE) {
## In FIMM, there are a lot of instances in which the ranges differ substantially.
## We would be taking a big hit to limit to the most restrictive range of a few cases.
cat("FIMM drugs with min or max concentrations differing by >= 0.01\n")
d_ply(fimm.raw.drug.ranges, .variables = c("DRUG_ID"),
                             .fun = function(df) {
                                      df$min.conc.nM <- as.numeric(df$min.conc.nM); df$max.conc.nM <- as.numeric(df$max.conc.nM)
                                      if( ( (max(df$min.conc.nM) - min(df$min.conc.nM)) >= 0.01 ) || ( (max(df$max.conc.nM) - min(df$max.conc.nM)) >= 0.01 ) ) {
                                         write.table(df, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
                                      }
                                    })

## Hence, define the most frequent range.  Then discard any cases that do not include that range.
fimm.raw.drug.ranges <- ddply(fimm.raw.drug.ranges, .variables = c("DRUG_ID"),
                                          .fun = function(df) {
                                             df[as.numeric(df$freq) == max(as.numeric(df$freq)),,drop=F]
                                          })

}

## Merge the two ranges together
shared.drug.anno <- ohsu.fimm.drugs[, !(grepl(x=colnames(ohsu.fimm.drugs), pattern="conc", ignore.case=TRUE) | (colnames(ohsu.fimm.drugs) == "X1"))]
ohsu.raw.drug.ranges <- merge(shared.drug.anno, ohsu.raw.drug.ranges, by.y = "inhibitor", by.x = "OHSU_DRUG_NAME")
colnames(ohsu.raw.drug.ranges)[colnames(ohsu.raw.drug.ranges) == "min.conc.nM"] <- "ohsu.min.conc.nM"
colnames(ohsu.raw.drug.ranges)[colnames(ohsu.raw.drug.ranges) == "max.conc.nM"] <- "ohsu.max.conc.nM"

colnames(fimm.raw.drug.ranges)[colnames(fimm.raw.drug.ranges) == "min.conc.nM"] <- "fimm.min.conc.nM"
colnames(fimm.raw.drug.ranges)[colnames(fimm.raw.drug.ranges) == "max.conc.nM"] <- "fimm.max.conc.nM"
raw.drug.ranges <- merge(ohsu.raw.drug.ranges, fimm.raw.drug.ranges[, c("DRUG_ID", "fimm.min.conc.nM", "fimm.max.conc.nM")], by.x = "FIMM_Batch_ID_Drug", by.y = "DRUG_ID")

## DSS code will use min.conc and max.conc columns
raw.drug.ranges$min.conc <- pmax(raw.drug.ranges$ohsu.min.conc.nM, raw.drug.ranges$fimm.min.conc.nM)
raw.drug.ranges$max.conc <- pmin(raw.drug.ranges$ohsu.max.conc.nM, raw.drug.ranges$fimm.max.conc.nM)

write.table(file="raw-drug-ranges.tsv", raw.drug.ranges, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)

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

## Drop a redundant column
ohsu.dss.common <- ohsu.dss.common[, !(colnames(ohsu.dss.common) == "Mechanism.Targets.1")]
fimm.dss.common <- fimm.dss.common[, !(colnames(fimm.dss.common) == "Mechanism.Targets.1")]

write.table(file="ohsu.dss.common.tsv", ohsu.dss.common, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
write.table(file="fimm.dss.common.tsv", fimm.dss.common, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
save.image(".Rdata")

## Recalculate DSS using both L.4 and LL.4 fits for varying thresholds t


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

executed.url <- "https://github.com/bswhite/fimm-ohsu-aml/blob/master/070617/fit-ohsu-and-fimm-drug-response-data.R"

thresholds <- c(0, 0.05, 0.1, 0.15, 0.2)
for(t in thresholds) {

  mydoc <- pptx()

  doc.title <- paste0('QC of IC50, AUC, and DSS (t = ', t, ')')  
  doc.file <- paste0("070717-metric-qc-t-", t, ".pptx")
  mydoc <- addSlide( mydoc, slide.layout = 'Title Slide' )
  mydoc <- addTitle( mydoc, doc.title )
  mydoc <- addSubtitle( mydoc , 'July 7, 2017')

  cat(paste0("Recalculating DSS for t = ", t, "\n"))

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

  path <- paste0("fimm.common.drugs.dss.t", t, ".tsv")
  write.table(file=path, fimm.dss, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
  f <- File(path, parentId = parentId, synapseStore = TRUE)
  synStore(f, executed = list(list(url=executed.url)), forceVersion = FALSE)

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

  path <- paste0("ohsu.common.drugs.dss.t", t, ".tsv")
  write.table(file=path, ohsu.dss, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
  f <- File(path, parentId = parentId, synapseStore = TRUE)
  synStore(f, executed = list(list(url=executed.url)), forceVersion = FALSE)

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
        jgc()
        mydoc <- addSlide( mydoc, slide.layout = "Title and Content" )
        title <- paste0(dataset, " (", metric.names[m], "): L.4 vs LL.4")
        mydoc <- addTitle( mydoc, title )
        mydoc <- addPlot( doc = mydoc, fun = print, x = g)
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
        jgc()
        mydoc <- addSlide( mydoc, slide.layout = "Title and Content" )
        title <- paste0("FIMM and OHSU ", metric.names[m], " ", model.names[o], " Densities")
        mydoc <- addTitle( mydoc, title )
        mydoc <- addPlot( doc = mydoc, fun = print, x = g)
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
      jgc()
      mydoc <- addSlide( mydoc, slide.layout = "Title and Content" )
      title <- paste0(dataset, " (", metric.names[m], "): L.4 vs LL.4")
      mydoc <- addTitle( mydoc, title )
      mydoc <- addPlot( doc = mydoc, fun = print, x = g)
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
      jgc()
      mydoc <- addSlide( mydoc, slide.layout = "Title and Content" )
      title <- paste0("FIMM and OHSU ", metric.names[m], " ", model.names[o], " Densities")
      mydoc <- addTitle( mydoc, title )
      mydoc <- addPlot( doc = mydoc, fun = print, x = g)
    }
  }

  writeDoc(mydoc, doc.file)
  jgc()

  ## Store the powerpoint in the Presentations folder
  parentId <- "syn8202511"
  f <- File(doc.file, parentId = parentId, synapseStore = TRUE)
  synStore(f, executed = list(list(url=executed.url)), forceVersion = FALSE)

}


