suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(synapseClient))
suppressPackageStartupMessages(library("foreach"))
suppressPackageStartupMessages(library("parallel"))

num.cores <- detectCores()
if(!is.na(num.cores) && (num.cores > 1)) {
  suppressPackageStartupMessages(library("doMC"))
  cat(paste("Registering ", num.cores-1, " cores.\n", sep=""))
  registerDoMC(cores=(num.cores-1))
}

synapseLogin()

cat("Use L.4 and LL.4 fits for CTRPv2.\n")
cat("Exclude both censored points and censored curves.\n")

## Load the CTRP drugs
## CTRP folder on Synapse: syn5622707
## CTRP drug info: "v20.meta.per_compound.txt" (syn5632193)
synId <- "syn5632193"
obj <- synGet(synId, downloadFile = TRUE)
ctrp.compounds <- read.table(getFileLocation(obj), sep="\t", header=TRUE, as.is=TRUE, comment.char="", quote="")

## CTRP drug data: "v20.data.per_cpd_well.txt" (syn5632189)
synId <- "syn5632189"
obj <- synGet(synId, downloadFile = TRUE)
ctrp.raw.data <- read.table(getFileLocation(obj), sep="\t", header=TRUE, as.is=TRUE, comment.char="", quote="")

## CTRP assay metadata (including DMSO): "v20.meta.per_assay_plate.txt" (syn5632191)
synId <- "syn5632191"
obj <- synGet(synId, downloadFile = TRUE)
ctrp.assay.metadata <- read.table(getFileLocation(obj), sep="\t", header=TRUE, as.is=TRUE, comment.char="", quote="")

## CTRP experiment metadata (including DMSO): "v20.meta.per_experiment.txt" (syn5632194)
synId <- "syn5632194"
obj <- synGet(synId, downloadFile = TRUE)
ctrp.expt.metadata <- read.table(getFileLocation(obj), sep="\t", header=TRUE, as.is=TRUE, comment.char="", quote="")

## CTRP cell line metadata: "v20.meta.per_cell_line.txt" (syn5632192)
synId <- "syn5632192"
obj <- synGet(synId, downloadFile = TRUE)
ctrp.cell.line.metadata <- read.table(getFileLocation(obj), sep="\t", header=TRUE, as.is=TRUE, comment.char="", quote="")

## CTRP post qc "v20.data.per_cpd_post_qc.txt" (syn5623539)
synId <- "syn5623539"
obj <- synGet(synId, downloadFile = TRUE)
ctrp.post.qc <- read.table(getFileLocation(obj), sep="\t", header=TRUE, as.is=TRUE, comment.char="", quote="")

## CTRP post qc "v20.data.per_cpd_pre_qc.txt" (syn5623903)
synId <- "syn5623903"
obj <- synGet(synId, downloadFile = TRUE)
ctrp.pre.qc <- read.table(getFileLocation(obj), sep="\t", header=TRUE, as.is=TRUE, comment.char="", quote="")

## CTRP drug data: "v20.data.curves_post_qc.txt" (syn5622711)
synId <- "syn5622711"
obj <- synGet(synId, downloadFile = TRUE)
ctrp.curves <- read.table(getFileLocation(obj), sep="\t", header=TRUE, as.is=TRUE, comment.char="", quote="")

## Perform some sanity checks to ensure we understand how individual concentration points and entire curves were censored.
## Page 4 of the supp lists the number of curves censored
num.censored.curves <- 1151 + 423 + 359 + 78 + 6 + 401

## The main text says that 397,681 fits were successfully computed (though these were subsequently censored)
num.successful.fits <- 397681

## The main text says that 395,263 curves were kept for final analysis.
num.final.curves <- 395263

## The main text says that there were 6,239,426 average PV response-point measurements (i.e., concentrations)
num.avg.resp.pts <- 6239426

if(num.final.curves == dim(ctrp.curves)[1]) {
  cat("As expected, ctrp.curves holds the final, post-censoring set of curves\n")
} else {
  stop(paste0("ctrp.curves (n = ", dim(ctrp.curves)[1], ") is not the final, post-censoring set of curves (n = ", num.final.curves, ")\n"))
}

if( ( num.successful.fits - num.censored.curves ) == dim(ctrp.curves)[1] ) {
  cat(paste0("# fits (", num.successful.fits, ") - # censored curves (", num.censored.curves, ") = final # curves (", dim(ctrp.curves)[1], "\n"))
} else {
  stop(paste0("# fits (", num.successful.fits, ") - # censored curves (",num.censored.curves,") != final # curves (", dim(ctrp.curves)[1], "\n"))
}

if(dim(ctrp.pre.qc)[1] == num.avg.resp.pts) {
  cat(paste0("As expected, num pre-qc avg PV responst-pt measurements = ", num.avg.resp.pts, "\n"))
} else {
  stop(paste0("Num pre-qc avg PV responst-pt measurements (n = ", dim(ctrp.pre.qc)[1], ") != that reported (n = ", num.avg.resp.pts, ")\n"))
}

if(dim(ctrp.pre.qc)[1] > dim(ctrp.post.qc)[1]) {
  cat(paste0("As expected, # of pre-qc data pts (n = ", dim(ctrp.pre.qc)[1], ") > # of post-qc data pts (n = ", dim(ctrp.post.qc)[1], ")\n"))
} else {
  stop(paste0("# of pre-qc data pts (n = ", dim(ctrp.pre.qc)[1], ") <= # of post-qc data pts (n = ", dim(ctrp.post.qc)[1], ")\n"))
}

old.sz <- nrow(ctrp.raw.data)
ctrp.data <- merge(ctrp.raw.data, ctrp.compounds[, c("master_cpd_id", "cpd_name")], by = "master_cpd_id")
new.sz <- nrow(ctrp.data)
if(old.sz != new.sz) {
  stop(paste0("Size changed from ", old.sz, " to ", new.sz, "\n"))
}

## Only fit curves to points that were not excluded/censored as outliers in original manuscript
ctrp.data <- merge(ctrp.data, ctrp.post.qc[, c("experiment_id", "master_cpd_id", "cpd_conc_umol")])

## Merge in the cancer cell line identifier
old.sz <- nrow(ctrp.data)
## ctrp.data <- merge(ctrp.data, unique(ctrp.expt.metadata[, c("experiment_id", "baseline_signal")]))
ctrp.data <- merge(ctrp.data, unique(ctrp.expt.metadata[, c("experiment_id", "master_ccl_id")]))
new.sz <- nrow(ctrp.data)
if(old.sz != new.sz) {
  stop(paste0("Size changed from ", old.sz, " to ", new.sz, "\n"))
}

sanity.check <- FALSE
## Sanity check that bsub_value_log2 = raw_value_log2 - dmso_plate_avg_log2
if(sanity.check) {

  ## Merge in the DMSO controls
  old.sz <- nrow(ctrp.data)
  ctrp.data <- merge(ctrp.data, ctrp.assay.metadata[, c("experiment_id", "assay_plate_barcode", "dmso_plate_avg_log2", "dmso_expt_std_log2")])
  new.sz <- nrow(ctrp.data)
  if(old.sz != new.sz) {
    stop(paste0("Size changed from ", old.sz, " to ", new.sz, "\n"))
  }

  tmp <- ctrp.data$raw_value_log2 - ctrp.data$dmso_plate_avg_log2 - ctrp.data$bsub_value_log2
  if(max(abs(tmp)) <= 0.02) {
    cat(paste0("Confirmed that bsub_value_log2 = raw_value_log2 - dmso_plate_avg_log2 (max diff = ", max(abs(tmp)), ")\n"))
  } else {
    stop(paste0("bsub_value_log2 != raw_value_log2 - dmso_plate_avg_log2 (max diff = ", max(abs(tmp)), ")\n"))
  }

  ## Merge in the predicted percent-viability from curve fit
  old.sz <- nrow(ctrp.data)
  ctrp.data <- merge(ctrp.data, ctrp.post.qc)
  new.sz <- nrow(ctrp.data)
  if(old.sz != new.sz) {
    stop(paste0("Size changed from ", old.sz, " to ", new.sz, "\n"))  
  }

  ## Ensure that predicted percent-viability is similar to 2^bsub
  tmp <- 2^ctrp.data$bsub_value_log2 - ctrp.data$cpd_pred_pv
  cat(paste0("Median 2^bsub - cpd_prev_pv = ", median(tmp), "\n"))
  cat(paste0("Std Dev 2^bsub - cpd_prev_pv = ", sd(tmp), "\n"))
}

## Scale concentrations from micro- to nano-molar
ctrp.data$conc.nM <- ctrp.data$cpd_conc_umol * 10^3 

## Note that viabilities generally max out at 1. 
## i.e., they are fractional viabilities.
## Convert to percent viability.
ctrp.data$percent_viability <- 100 * (2^ctrp.data$bsub_value_log2)

source("../common/dss.R")
source("../common/drc-fit.R")

drug.screen.cols <- c("experiment_id", "master_cpd_id", "master_ccl_id")
stop("Need to redo this to account for replicates--by including assay_plate_barcode above.  This change will need to be made elsewhere as well\n")
conc.col <- "conc.nM"
response.col <- "percent_viability"
all.cols <- c(drug.screen.cols, conc.col, response.col)

cat("Fitting CTRPv2 using LL4\n")
ll4.fits <- fit.drc(ctrp.data[,all.cols], fct = "LL.4", drug.screen.cols = drug.screen.cols, conc.nM.col = conc.col, response.col = response.col, response.is.viability = TRUE)
cat("Done fitting CTRPv2 using LL4\n")
save.image(".Rdata")

cat("Fitting CTRPv2 using L4\n")
l4.fits <- fit.drc(ctrp.data[,all.cols], fct = "L.4", drug.screen.cols = drug.screen.cols, conc.nM.col = conc.col, response.col = response.col, response.is.viability = TRUE)
cat("Done fitting CTRPv2 using L4\n")
save.image(".Rdata")

common.cols <- c(drug.screen.cols, "min.conc.nM", "max.conc.nM", "all.concs.nM", "uniq.concs.nM")
fits <- merge(ll4.fits, l4.fits, by = common.cols, suffixes = c(".ll4", ".l4"), all = TRUE)

## Merge in the cancer cell line metadata
old.sz <- nrow(fits)
fits <- merge(fits, ctrp.cell.line.metadata)
new.sz <- nrow(fits)
if(old.sz != new.sz) {
  stop(paste0("Size changed from ", old.sz, " to ", new.sz, "\n"))
}

save.image(".Rdata")

## Store the CTRP fits in CTRP Data/Processed Data (syn10288724)
parentId <- "syn10288724"

path <- paste0("ctrp.l4.and.ll4.fits.non.censored.pts.tsv")
write.table(file=path, fits, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
f <- File(path, parentId = parentId, synapseStore = TRUE)
synStore(f, executed = NULL, forceVersion = FALSE)

## Finally, censor the curves that were filtered in the published manuscript
non.censored.fits <- merge(fits, ctrp.curves[, c("experiment_id", "master_cpd_id")])
path <- paste0("ctrp.l4.and.ll4.fits.non.censored.pts.non.censored.curves.tsv")
write.table(file=path, non.censored.fits, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
f <- File(path, parentId = parentId, synapseStore = TRUE)
synStore(f, executed = NULL, forceVersion = FALSE)

cat("Fit CTRPv2 data using L.4 and LL.4\n")