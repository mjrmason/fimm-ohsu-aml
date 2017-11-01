suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(synapseClient))
suppressPackageStartupMessages(library("foreach"))
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library("plyr"))

num.cores <- detectCores()
if(!is.na(num.cores) && (num.cores > 1)) {
  suppressPackageStartupMessages(library("doMC"))
  cat(paste("Registering ", num.cores-1, " cores.\n", sep=""))
  registerDoMC(cores=(num.cores-1))
}

synapseLogin()

cat("Use DRC to fit L.4 (4-parameter logistic) and LL.4 (4-paramter log logistic) functions to Sanger/GDSC data.\n")

## Sanger/GDSC data is located in the following synapse folder
sanger.syn.folder.id <- "syn7466648"

## Read in the drug sensitivity data.
synId <- "syn11275807"
drug.data <- openxlsx:::read.xlsx(getFileLocation(synGet(synId)), sheet=1)

## Reformat the sanger data such that there is one concentration/response per row.
## NB: concentrations are in uM (as described in drug_sensitivity_raw_data_definitions.xlsx -- syn11275820)
## NB: drug_sensitivity_raw_data_definitions.xlsx describes the data as:
## MAX_CONC (in micromolar)
## FOLD_DILUTION: fold change in dilution between drugged wells from the MAX_CONC
## raw_max: intensity reading for well treated with the maximum concnetration of drug.
## raw2: intensity reading for well drug treated with a concentration a fold dilution less than that of raw_max
## raw3: intensity reading for well drug treated with a concentration a fold dilution less than that of raw2
## etc.
## controlx: intensity reading for negative control well (cell but no drugs)
## blankx: intensity reading for positive control well (no cells and no drugs)
## The supplement indicates the the fits are to raw / control (A Landscape of Pharmacogenomic Interactions in Cancer; Iorio et al)
## The assumed dose-response model is a two-parameter sigmoid that models the relative viability. The latter is obtained by scaling the observed intensities to the mean intensities of the control wells. Since the assumed response is strictly between 0 and 100% relative viability, the values are capped at 0 and 100.
## Hence, these are viaibilty data.
colnames(drug.data)[colnames(drug.data) == "raw_max"] <- "raw1"
## And explicit concentration rows
data.cols <- colnames(drug.data)[grepl(pattern="raw", colnames(drug.data))]
control.cols <- colnames(drug.data)[grepl(pattern="control", colnames(drug.data))]
max.raw.well <- max(as.numeric(gsub(data.cols, pattern="raw", replacement="")))
if(max.raw.well != length(data.cols)) {
  stop(paste0("Unexpected non-contiguous number of raw wells: ", paste0(data.cols, collapse=","), "\n"))
}
conc.cols <- gsub(data.cols, pattern="raw", replacement="conc")
drug.data$conc1 <- drug.data$MAX_CONC
for(i in 2:max.raw.well) {
  prev.conc.col <- paste0("conc", i-1)
  conc.col <- paste0("conc", i)
  drug.data[, conc.col] <- as.numeric(drug.data[, prev.conc.col]) / as.numeric(drug.data$FOLD_DILUTION)
}

drug.screen.cols <- c("IC_RESULTS_ID", "COSMIC_ID", "CELL_LINE_NAME", "DRUG_ID")
if(any(duplicated(drug.data$IC_RESULTS_ID))) {
  stop("Unexpected: thought IC_RESULTS_ID was a unique identifier!\n")
}
tmp <- drug.data[, c(data.cols, conc.cols, control.cols, drug.screen.cols)]
drug.data.long <- ldply(1:nrow(tmp), .parallel = FALSE,
                             .fun = function(i) {
                               mean.cntrl <- mean(as.numeric(tmp[i, control.cols]), na.rm=TRUE)
                               dfs <- lapply(data.cols,
                                                    function(data.col) {
                                                      conc.col <- gsub(data.col, pattern="raw", replacement="conc")
                                                      ret <- tmp[i, c(drug.screen.cols, data.col, conc.col)]
                                                      colnames(ret) <- c(drug.screen.cols, "percent_viability", "conc.uM")
                                                      ret$percent_viability <- 100 * as.numeric(ret$percent_viability) / as.numeric(mean.cntrl)
                                                      ret
                                                    })
                               do.call("rbind", dfs)
                             })

save.image(".Rdata")

## Scale concentrations from micro- to nano-molar
drug.data.long$conc.nM <- drug.data.long$conc.uM * 10^3 

source("../common/dss.R")
source("../common/drc-fit.R")

conc.col <- "conc.nM"
response.col <- "percent_viability"
all.cols <- c(drug.screen.cols, conc.col, response.col)

cat("Fitting Sanger using LL4\n")
ll4.fits <- fit.drc(drug.data.long[,all.cols], fct = "LL.4", drug.screen.cols = drug.screen.cols, conc.nM.col = conc.col, response.col = response.col, 
                    response.is.viability = TRUE)
cat("Done fitting Sanger using LL4\n")
save.image(".Rdata")

if(nrow(ll4.fits) != nrow(drug.data)) {
  stop(paste0("Number of LL4 fits (", nrow(ll4.fits), ") != Number of original entries (", nrow(drug.data), ")\n"))
}

cat("Fitting Sanger using L4\n")
l4.fits <- fit.drc(drug.data.long[,all.cols], fct = "L.4", drug.screen.cols = drug.screen.cols, conc.nM.col = conc.col, response.col = response.col, response.is.viability = TRUE)
cat("Done fitting Sanger using L4\n")
save.image(".Rdata")

if(nrow(l4.fits) != nrow(drug.data)) {
  stop(paste0("Number of L4 fits (", nrow(l4.fits), ") != Number of original entries (", nrow(drug.data), ")\n"))
}

common.cols <- c(drug.screen.cols, "min.conc.nM", "max.conc.nM", "all.concs.nM", "uniq.concs.nM")
fits <- merge(ll4.fits, l4.fits, by = common.cols, suffixes = c(".ll4", ".l4"), all = TRUE)

## Store in
## FIMM_BEAT AML Collaboration/Files/Sanger Data
parentId <- "syn11287866"

path <- paste0("sanger.l4.and.ll4.fits.tsv")
write.table(file=path, fits, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
f <- File(path, parentId = parentId, synapseStore = TRUE)
synStore(f, executed = NULL, forceVersion = FALSE)

cat("Fit Sanger data using L.4 and LL.4\n")

q()
