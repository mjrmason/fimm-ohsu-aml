## Begin setup (this will need to be changed as the data sets in the intersection change)

## Data sets to be analyzed:
data.sets <- list("ohsu" = "ohsu", "fimm" = "fimm", "ctrp" = "ctrp", "sanger" = "sanger")
data.sets.to.analyze <- unlist(data.sets)

## data.set.drug.fit.synIds
## data.set.expr.synIds
## drug.map.drug.id.cols 
## data.set.drug.id.cols
## data.sets

use.subset <- FALSE

## The drug mapping file lists the drugs in common across all data sets with the drug.map.drug.id.cols
## providing the identifier for each drug for each respective data set.
## The drug mapping file may be created by a script such as 
## harmonize-fimm-ohsu-ctrp-sanger-compounds.R
drug.mapping.synId <- "syn11290968"
drug.map.drug.id.cols <- list("ohsu" = "OHSU_DRUG_NAME", "fimm" = "FIMM_Batch_ID_Drug", "ctrp" = "master_cpd_id", "sanger" = "Drug.ID")

## These were created by the scripts:
## fit-fimm-drc-curves.R
## fit-ohsu-drc-curves.R
## ./102117/fit-sanger-drc-curves.R
## ./080817/fit-ctrp-drc-curves-all.R
## ./080817/fit-ctrp-drc-curves-non-censored.R

data.set.drug.fit.synIds <- list("ohsu" = "syn11470803", "fimm" = "syn11471586", "ctrp" = "syn10307313", "sanger" = "syn11289087")
data.set.drug.fit.files <- c("ohsu.l4.and.ll4.fits.112017.tsv", "fimm.l4.and.ll4.fits.112017.tsv", "ctrp.l4.and.ll4.fits.non.censored.pts.non.censored.curves.tsv", "sanger.l4.and.ll4.fits.tsv")

## The column name within the drug map file that will be used to identify drugs (across all data sets)
drug.name.col <- "FIMM_DRUG_NAME"

## These fit files are assumed to be for AUCs/DSS values calculated across shared concentration ranges across data sets.
## These files may be created by a script such as
## calculate-dss.R
data.set.dss.fit.synIds <- NULL
data.set.dss.fit.files <- NULL

## The columns within each respective data set that identifies the drug (in the name space of that data set)
data.set.drug.id.cols <- list("ohsu" = "inhibitor", "fimm" = "DRUG_ID", "ctrp" = "master_cpd_id", "sanger" = "DRUG_ID")

## The columns within each respective data set that identifies a screen (i.e., a particular replicate of a patient/cell line x drug)
screen.id.cols <- list("ohsu" = c("patient_id", "inhibitor", "lab_id", "SeqID", "replicant", "time_of_read"),
                       "fimm" = c("SCREEN_ID", "PATIENT_ID", "DRUG_ID", "DRUG_NAME", "DRUG_SET", "SAMPLE_DATE"),
                       "ctrp" = c("experiment_id", "master_cpd_id", "master_ccl_id", "ccl_name"), 
                       "sanger" = c("IC_RESULTS_ID", "COSMIC_ID", "CELL_LINE_NAME", "DRUG_ID"))

## Batch-corrected log2 cpm (or rma) expression files (NB: these exclude 2 outliers in OHSU)
expr.folder.synIds <- list("ohsu" = "syn11288886", "fimm" = "syn11288886", "ctrp" = "syn11288886", "sanger" = "syn11288886")
data.set.expr.synIds <- list("ohsu" = "syn11307673", "fimm" = "syn11307680", "ctrp" = "syn11307674", "sanger" = "syn11307678")
data.set.expr.files <- list("ohsu" = "ohsu-expr-foc-aml-s-aml-outlier1-combat.tsv", "fimm" = "fimm-expr-foc-aml-s-aml-outlier1-combat.tsv",
                            "ctrp" = "ctrp-expr-foc-aml-s-aml-outlier1-combat.tsv", "sanger" = "sanger-expr-foc-aml-s-aml-outlier1-combat.tsv")

## Columns whose values are not dependent on the (L.4 vs LL.4) fit.  
merge.cols <- list(
  "ohsu" = c(screen.id.cols[["ohsu"]], "min.conc", "max.conc", "all.concs.nM", "uniq.concs.nM", "all.responses",
             "min.conc.w.outliers.nM", "max.conc.w.outliers.nM", "all.concs.w.outliers.nM", "uniq.concs.w.outliers.nM", "all.responses.w.outliers"),
  "fimm" = c(screen.id.cols[["fimm"]], "min.conc", "max.conc", "all.concs.nM", "uniq.concs.nM", "all.responses",
             "min.conc.w.outliers.nM", "max.conc.w.outliers.nM", "all.concs.w.outliers.nM", "uniq.concs.w.outliers.nM", "all.responses.w.outliers"),
  "ctrp" = c(screen.id.cols[["ctrp"]], "min.conc", "max.conc", "all.concs.nM", "uniq.concs.nM", "all.responses",
             "min.conc.w.outliers.nM", "max.conc.w.outliers.nM", "all.concs.w.outliers.nM", "uniq.concs.w.outliers.nM", "all.responses.w.outliers"),
  "sanger" = c(screen.id.cols[["sanger"]], "min.conc", "max.conc", "all.concs.nM", "uniq.concs.nM", "all.responses",
             "min.conc.w.outliers.nM", "max.conc.w.outliers.nM", "all.concs.w.outliers.nM", "uniq.concs.w.outliers.nM", "all.responses.w.outliers"))

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

patient.subsets <- list("ohsu" = NULL, "fimm" = NULL, "ctrp" = ccle.subset$ccl_name, "sanger" = sanger.subset$COSMIC.identifier)
patient.subsets <- list("ohsu" = NULL, "fimm" = NULL, "ctrp" = ccle.subset$master_ccl_id, "sanger" = sanger.subset$COSMIC.identifier)

patient.id.cols <- list("ohsu" = "patient_id", "fimm" = "PATIENT_ID", "ctrp" = "master_ccl_id", "sanger" = "COSMIC_ID")
expr.patient.id.cols <- list("ohsu" = "SeqID", "fimm" = "SCREEN_ID", "ctrp" = "ccl_name", "sanger" = "COSMIC_ID")

## End setup
