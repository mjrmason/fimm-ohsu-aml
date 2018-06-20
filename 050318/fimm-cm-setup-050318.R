## Begin setup (this will need to be changed as the data sets in the intersection change)

conc.cols <- list("fimm-cm" = "conc.nM", "fimm" = "conc.nM")
response.cols <- list("fimm-cm" = "PERCENT_INHIBITION", "fimm" = "PERCENT_INHIBITION")
responses.are.viabilities <- list("fimm-cm" = FALSE, "fimm" = FALSE)

## Data sets to be analyzed:
data.sets <- list("fimm-cm" = "fimm-cm", "fimm" = "fimm")
data.sets.to.analyze <- unlist(data.sets)

## Check these variables that were just changed to lists
## data.set.drug.fit.synIds
## data.set.expr.synIds
## drug.map.drug.id.cols 
## data.set.drug.id.cols
## data.sets

## The drug mapping file lists the drugs in common across the FIMM and OHSU data sets with the drug.map.drug.id.cols
## providing the identifier for each drug for each respective data set.
## The drug mapping file may be created by a script such as 
## harmonize-fimm-ohsu-compounds.R
drug.mapping.synId <- "syn11361110" ## fimm-ohsu-drug-map.tsv
drug.map.drug.id.cols <- list("ohsu" = "OHSU_DRUG_NAME", "fimm" = "FIMM_Batch_ID_Drug", "fimm-cm" = "FIMM_Batch_ID_Drug")


## These were created by the scripts in this directory:
## fit-fimm-drc-curves.R
## fit-ohsu-drc-curves.R
data.set.drug.fit.synIds <- list("ohsu" = "syn11470803", "fimm" = "syn11471586", "fimm-cm" = "syn12213740")
data.set.drug.fit.files <- c("ohsu.l4.and.ll4.fits.112017.tsv", "fimm.l4.and.ll4.fits.112017.tsv", "fimm.cm.l4.and.ll4.fits.050318.tsv")

## The column name within the drug map file that will be used to identify drugs (across all data sets)
drug.name.col <- "FIMM_DRUG_NAME"

## These fit files are assumed to be for AUCs/DSS values calculated across shared concentration ranges across data sets.
## These files may be created by a script such as
## calculate-dss.R
data.set.dss.fit.synIds <- list("ohsu" = "syn11472590", "fimm" = "syn11472591")
data.set.dss.fit.files <- c("ohsu-fimm-ohsu-112017-dss-t0.tsv", "fimm-fimm-ohsu-112017-dss-t0.tsv")

processed.drug.fit.synIds <- list("ohsu" = "syn11472590", "fimm" = "syn11472591", "fimm-cm" = "syn12213740")
processed.drug.fit.files <- c("ohsu-fimm-ohsu-112017-dss-t0.tsv", "fimm-fimm-ohsu-112017-dss-t0.tsv", "fimm.cm.l4.and.ll4.fits.050318.tsv")

## The columns within each respective data set that identifies the drug (in the name space of that data set)
data.set.drug.id.cols <- list("ohsu" = "inhibitor", "fimm" = "DRUG_ID", "fimm-cm" = "DRUG_ID")

## The columns within each respective data set that identifies a screen (i.e., a particular replicate of a patient/cell line x drug)
screen.id.cols <- list("ohsu" = c("patient_id", "inhibitor", "lab_id", "SeqID", "replicant", "time_of_read"),
                       "fimm" = c("SCREEN_ID", "PATIENT_ID", "DRUG_ID", "DRUG_NAME", "DRUG_SET", "SAMPLE_DATE"),
                       "fimm-cm" = c("SCREEN_ID", "PATIENT_ID", "DRUG_ID", "DRUG_NAME", "DRUG_SET", "SAMPLE_DATE"))

## Batch-corrected log2 cpm (or rma) expression files (NB: these exclude 2 outliers in OHSU and are highly variable/cancer genes)
data.set.expr.synIds <- list("ohsu" = "syn11498061", "fimm" = "syn11498062")
expr.folder.synIds <- list("ohsu" = "syn11361089", "fimm" = "syn11361089")
data.set.expr.files <- list("ohsu" = "ohsu-expr-fimm-ohsu-outlier1-combat-variable-and-cancer.tsv", "fimm" = "fimm-expr-fimm-ohsu-outlier1-combat-variable-and-cancer.tsv")

drug.response.cols <- list("ohsu" = "dss.auc.ll4", "fimm" = "dss.auc.ll4", "fimm-cm" = "auc.ll4")

min.conc.cols <- list("ohsu" = "min.conc.nM", "fimm" = "min.conc.nM", "fimm-cm" = "min.auc.conc.ll4")
max.conc.cols <- list("ohsu" = "max.conc.nM", "fimm" = "max.conc.nM", "fimm-cm" = "max.auc.conc.ll4")

if(TRUE) {
  ## Overwrite the expression synIds to exclude outliers (as before), but _not_ to restrict to highly variable/cancer genes
  ## Batch-corrected log2 cpm (or rma) expression files (NB: these exclude 2 outliers in OHSU)
  data.set.expr.synIds <- list("ohsu" = "syn11362256", "fimm" = "syn11362257")
  ## names(data.set.expr.synIds) <- data.sets
  expr.folder.synIds <- list("ohsu" = "syn11361089", "fimm" = "syn11361089")
  data.set.expr.files <- list("ohsu" = "ohsu-expr-fimm-ohsu-outlier1-combat.tsv", "fimm" = "fimm-expr-fimm-ohsu-outlier1-combat.tsv")
}

## These are the original expression matrices before combat
## data.set.orig.expr.synIds <- list("ohsu" = "syn10083723", "fimm" = "syn11288849", "lee" = "syn12176399", "fimm-cm" = "syn12498087")
## data.set.orig.expr.files <- list("ohsu" = "ohsu.expr.tsv", "fimm" = "fimm-log-cpm-gene-expr.tsv", "lee" = "lee-patient-expr.tsv", "fimm-cm" = "RNASeq.CPM.log2.bc.fixed.names.csv")

data.set.orig.expr.synIds <- list("ohsu" = "syn10083723", "fimm" = "syn11288849", "lee" = "syn12176399", "fimm-cm" = "syn12498090")
data.set.orig.expr.files <- list("ohsu" = "ohsu.expr.tsv", "fimm" = "fimm-log-cpm-gene-expr.tsv", "lee" = "lee-patient-expr.tsv", "fimm-cm" = "RNASeq.CPM.log2.bc.fixed.names.tsv")

## Columns whose values are not dependent on the (L.4 vs LL.4) fit.  
merge.cols <- list(
  "ohsu" = c(screen.id.cols[["ohsu"]], "min.conc", "max.conc", "all.concs.nM", "uniq.concs.nM", "all.responses",
             "min.conc.w.outliers.nM", "max.conc.w.outliers.nM", "all.concs.w.outliers.nM", "uniq.concs.w.outliers.nM", "all.responses.w.outliers"),
  "fimm" = c(screen.id.cols[["fimm"]], "min.conc", "max.conc", "all.concs.nM", "uniq.concs.nM", "all.responses",
             "min.conc.w.outliers.nM", "max.conc.w.outliers.nM", "all.concs.w.outliers.nM", "uniq.concs.w.outliers.nM", "all.responses.w.outliers"),
  "fimm-cm" = c(screen.id.cols[["fimm-cm"]], "min.conc", "max.conc", "all.concs.nM", "uniq.concs.nM", "all.responses",
             "min.conc.w.outliers.nM", "max.conc.w.outliers.nM", "all.concs.w.outliers.nM", "uniq.concs.w.outliers.nM", "all.responses.w.outliers"))

## The synapseIds of folders in which to store the fits for each respective data set.  FIMM_BEAT AML Collaboration/Files/Analyses/fimm-ohsu/112217 ("syn11498066")

fit.folder.synIds <- list("ohsu" = "syn11498066", "fimm" = "syn11498066")
## A string that will be included in the fit file name.
fit.file.prefix <- "fimm-ohsu-112217"
output.dir <- "output"

patient.id.cols <- list("ohsu" = "patient_id", "fimm" = "PATIENT_ID", "fimm-cm" = "PATIENT_ID")
patient.subsets <- list("ohsu" = NULL, "fimm" = NULL)

expr.patient.id.cols <- list("ohsu" = "SeqID", "fimm" = "SCREEN_ID", "fimm-cm" = "SCREEN_ID")
sample.id.cols <- list("ohsu" = "SeqID", "fimm" = "SCREEN_ID", "fimm-cm" = "SCREEN_ID")

use.subset <- TRUE

## HERE

metadata.synIds <- list("fimm-cm" = "syn12493054", "fimm" = "syn8270594")

raw.drug.synIds <- list("fimm-cm" = "syn12180768", "fimm" = "syn8488824")

experimental.factor.cols <- list("fimm-cm" = c("DRUG_NAME", "DRUG_SET", "SCREEN_ID", "PATIENT_ID"),
                                 "fimm" = c("DRUG_NAME", "DRUG_SET", "SCREEN_ID", "PATIENT_ID"))

## Store in
## FIMM_BEAT AML Collaboration/Files/FIMM Data/Processed Data
processed.folder.synIds <- list("fimm-cm" = "syn8270577", "fimm" = "syn8270577")


## End setup
