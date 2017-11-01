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

cat("Calculate AUC and DSS values based on 3-way intersection of concentrations\n")
cat("for each drug shared between FIMM, OHSU, and CTRPv2\n")

## Load the map of drugs shared between OHSU and FIMM (FIMM_OHSU_Drug_Concentrations.xlsx)
synId <- "syn10163669"
obj <- synGet(synId, downloadFile = TRUE)
ohsu.fimm.drugs <- openxlsx:::read.xlsx(getFileLocation(obj), sheet=1)

## Load the CTRP drugs
## CTRP folder on Synapse: syn5622707
## CTRP drug info: "v20.meta.per_compound.txt" (syn5632193)
synId <- "syn5632193"
obj <- synGet(synId, downloadFile = TRUE)
ctrp.compounds <- read.table(getFileLocation(obj), sep="\t", header=TRUE, as.is=TRUE, comment.char="", quote="")

## Match CTRP drugs to the drugs in common between FIMM and OHSU
ohsu.fimm.drugs$CTRP_DRUG_NAME <- NA
for(i in 1:nrow(ohsu.fimm.drugs)) {
  drugs <- unique(unname(unlist(lapply(ohsu.fimm.drugs[i,c("OHSU_DRUG_NAME", "FIMM_DRUG_NAME", "Non.proprietary.name", "Trade.names", "Aliases")], function(str) unlist(strsplit(unname(str), split=",[ ]*"))))))
  drugs <- drugs[!is.na(drugs)]
  flag <- grepl(pattern="\\(", drugs)
  if(any(flag)) {
    drugs <- c(drugs, unlist(lapply(drugs[flag], function(str) c(gsub(".*\\((.*)\\).*", "\\1", str), gsub("(.*) \\(.*\\).*", "\\1", str)))))
  }
  matches <- c()
  for(drug in tolower(drugs)) {
    if(any(tolower(ctrp.compounds$cpd_name) == drug)) {
      matches <- c(matches, which(tolower(ctrp.compounds$cpd_name) == drug))
    }
  }
  matches <- unique(matches)
  if(length(matches) > 1) { stop("Got multiple matches to CTRP!\n") }
  if(length(matches) > 0) {
    ohsu.fimm.drugs$CTRP_DRUG_NAME[i] <- ctrp.compounds$cpd_name[matches[1]]
  }
}

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

old.sz <- nrow(ctrp.raw.data)
ctrp.data <- merge(ctrp.raw.data, ctrp.compounds[, c("master_cpd_id", "cpd_name")], by = "master_cpd_id")
new.sz <- nrow(ctrp.data)
if(old.sz != new.sz) {
  stop(paste0("Size changed from ", old.sz, " to ", new.sz, "\n"))
}

## Limit the data to those that have OHSU and FIMM drugs
ctrp.data <- subset(ctrp.data, cpd_name %in% ohsu.fimm.drugs$CTRP_DRUG_NAME)

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
conc.col <- "conc.nM"
response.col <- "percent_viability"
all.cols <- c(drug.screen.cols, conc.col, response.col)

## CTRP drug data: "v20.data.curves_post_qc.txt" (syn5622711)
synId <- "syn5622711"
obj <- synGet(synId, downloadFile = TRUE)
ctrp.curves <- read.table(getFileLocation(obj), sep="\t", header=TRUE, as.is=TRUE, comment.char="", quote="")

min.conc <- ddply(ctrp.data, c("experiment_id", "master_cpd_id"), .fun = function(df) min(df$cpd_conc_umol))
max.conc <- ddply(ctrp.data, c("experiment_id", "master_cpd_id"), .fun = function(df) max(df$cpd_conc_umol))
colnames(min.conc) <- c("experiment_id", "master_cpd_id", "min_cpd_conc_umol")
colnames(max.conc) <- c("experiment_id", "master_cpd_id", "max_cpd_conc_umol")
ctrp.curves <- merge(ctrp.curves, min.conc)
ctrp.curves <- merge(ctrp.curves, max.conc)

aucs <- 
  ddply(ctrp.curves,
        c("experiment_id", "master_cpd_id"),
        .parallel = TRUE,
        .fun = function(df) {
          ## compute.ctrp.auc(b = df$p4_baseline, h = df$p3_total_decline, alpha = df$p1_center, beta = df$p2_slope, min.conc = log2(df$min_cpd_conc_umol), max.conc = log2(df$max_cpd_conc_umol))
          compute.ctrp.auc(b = df$p4_baseline, h = df$p3_total_decline, alpha = df$p1_center, beta = df$p2_slope, min.conc = df$min_cpd_conc_umol, max.conc = df$max_cpd_conc_umol)
        })
colnames(aucs) <- c("experiment_id", "master_cpd_id", "numerical_auc")    
m <- merge(ctrp.curves, aucs)
## plot(m$area_under_curve, m$numerical_auc)
## The AUCs we have computed are different than the published AUCs.  The published ones appear to
## integrate a logistic function from log2(min.conc) to log2(max.conc), but without a Jacobian of
## transformation.  Instead, we integrate a log-logistic function from min.conc to max.conc, which 
## properly handles the transformation from real to log space.
smoothScatter(m$area_under_curve, m$numerical_auc)
##lmf <- lm(area_under_curve ~ numerical_auc, data = m)
##coefficients(lmf)[2]

save.image(".Rdata")
