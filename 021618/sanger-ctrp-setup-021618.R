## Begin setup (this will need to be changed as the data sets in the intersection change)

## Data sets to be analyzed:
data.sets <- list("ctrp" = "ctrp", "sanger" = "sanger")
data.sets.to.analyze <- unlist(data.sets)

## Check these variables that were just changed to lists
## data.set.drug.fit.synIds
## data.set.expr.synIds
## drug.map.drug.id.cols 
## data.set.drug.id.cols
## data.sets

use.subset <- FALSE

## The drug mapping file lists the drugs in common across all data sets with the drug.map.drug.id.cols
## providing the identifier for each drug for each respective data set.
## The drug mapping file may be created by a script such as 
## harmonize-ctrp-sanger-compounds.R
drug.mapping.synId <- "syn11297396" ## fimm-ohsu-drug-map.tsv
drug.map.drug.id.cols <- list("ctrp" = "master_cpd_id", "sanger" = "Drug.ID")

## These were created by the scripts in this directory:
## fit-sanger-drc-curves.R
## fit-ctrp-drc-curves-non-censored.R
data.set.drug.fit.synIds <- c("syn10307313", "syn11289087")
names(data.set.drug.fit.synIds) <- data.sets
data.set.drug.fit.files <- c("ctrp.l4.and.ll4.fits.non.censored.pts.non.censored.curves.tsv", "sanger.l4.and.ll4.fits.tsv")

## The column name within the drug map file that will be used to identify drugs (across all data sets)
drug.name.col <- "cpd_name"

## These fit files are assumed to be for AUCs/DSS values calculated across shared concentration ranges across data sets.
## These files may be created by a script such as
## calculate-dss.R
data.set.dss.fit.synIds <- data.set.drug.fit.synIds
data.set.dss.fit.files <- data.set.drug.fit.files

data.set.drug.id.cols <- list("ctrp" = "master_cpd_id", "sanger" = "DRUG_ID")

screen.id.cols <- list(
                       "ctrp" = c("experiment_id", "master_cpd_id", "master_ccl_id", "ccl_name"), 
                       "sanger" = c("IC_RESULTS_ID", "COSMIC_ID", "CELL_LINE_NAME", "DRUG_ID"))

## Columns whose values are not dependent on the (L.4 vs LL.4) fit.  
merge.cols <- list(
  "ctrp" = c(screen.id.cols[["ctrp"]], "min.conc", "max.conc", "all.concs.nM", "uniq.concs.nM", "all.responses",
             "min.conc.w.outliers.nM", "max.conc.w.outliers.nM", "all.concs.w.outliers.nM", "uniq.concs.w.outliers.nM", "all.responses.w.outliers"),
  "sanger" = c(screen.id.cols[["sanger"]], "min.conc", "max.conc", "all.concs.nM", "uniq.concs.nM", "all.responses",
             "min.conc.w.outliers.nM", "max.conc.w.outliers.nM", "all.concs.w.outliers.nM", "uniq.concs.w.outliers.nM", "all.responses.w.outliers"))


## Begin setup (this will need to be changed as the data sets in the intersection change)
data.sets <- c("ctrp", "sanger")
data.set.expr.synIds <- c("syn11261484", "syn11288436")
names(data.set.expr.synIds) <- data.sets
data.set.expr.files <- c("ctrpv2-log-cpm-gene-expr.tsv", "sanger-rma-gene-expr.tsv")

## Pull in the CTRP cell line metadata so we can see the cancer/histology type of each sample
synId <- "syn5632192"
obj <- synGet(synId, downloadFile = TRUE)
ctrp.cell.line.metadata <- read.table(getFileLocation(obj), sep="\t", header=TRUE, as.is=TRUE, comment.char="", quote="")

## Subset the cell lines to those of interest
ccle.subset <- ctrp.cell.line.metadata[!is.na(ctrp.cell.line.metadata$ccle_hist_subtype_1),]
ccle.subset <- subset(ccle.subset, ccle_hist_subtype_1 == "acute_myeloid_leukaemia")

## Read in the Sanger cell line metadata so we can see the cancer type of each sample (Cell_Lines_Details.xlsx)
synId <- "syn11275809"
sanger.cell.line.metadata <- openxlsx:::read.xlsx(getFileLocation(synGet(synId)), sheet=1)
sanger.subset <- sanger.cell.line.metadata[!is.na(sanger.cell.line.metadata$`Cancer.Type.(matching.TCGA.label)`),]
sanger.subset <- sanger.subset[sanger.subset$`Cancer.Type.(matching.TCGA.label)` == "LAML", ]

patient.subsets <- list("ctrp" = ccle.subset$ccl_name, "sanger" = sanger.subset$COSMIC.identifier)
patient.subsets <- list("ctrp" = ccle.subset$master_ccl_id, "sanger" = sanger.subset$COSMIC.identifier)

patient.id.cols <- list("ctrp" = "master_ccl_id", "sanger" = "COSMIC_ID")
expr.patient.id.cols <- list("ctrp" = "ccl_name", "sanger" = "COSMIC_ID")

## End setup
