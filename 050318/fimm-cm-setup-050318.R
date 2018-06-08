## Begin setup (this will need to be changed as the data sets in the intersection change)

conc.cols <- list("fimm-cm" = "conc.nM")
response.cols <- list("fimm-cm" = "PERCENT_INHIBITION")
responses.are.viabilities <- list("fimm-cm" = FALSE)

## Data sets to be analyzed:
data.sets <- list("fimm-cm" = "fimm-cm")
data.sets.to.analyze <- unlist(data.sets)

## Check these variables that were just changed to lists
## data.set.drug.fit.synIds
## data.set.expr.synIds
## drug.map.drug.id.cols 
## data.set.drug.id.cols
## data.sets

raw.drug.synIds <- list("fimm-cm" = "syn12180768")

## The columns within each respective data set that identifies the drug (in the name space of that data set)
data.set.drug.id.cols <- list("fimm-cm" = "DRUG_ID")

## The columns within each respective data set that identifies a screen (i.e., a particular replicate of a patient/cell line x drug)
screen.id.cols <- list("fimm-cm" = c("SCREEN_ID", "PATIENT_ID", "DRUG_ID", "DRUG_NAME", "DRUG_SET", "SAMPLE_DATE"))
experimental.factor.cols <- list("fimm-cm" = c("DRUG_NAME", "DRUG_SET", "SCREEN_ID", "PATIENT_ID"))

## Store in
## FIMM_BEAT AML Collaboration/Files/FIMM Data/Processed Data
processed.folder.synIds <- list("fimm-cm" = "syn8270577")

## Columns whose values are not dependent on the (L.4 vs LL.4) fit.  
merge.cols <- list(
  "fimm-cm" = c(screen.id.cols[["fimm-cm"]], "min.conc", "max.conc", "all.concs.nM", "uniq.concs.nM", "all.responses",
             "min.conc.w.outliers.nM", "max.conc.w.outliers.nM", "all.concs.w.outliers.nM", "uniq.concs.w.outliers.nM", "all.responses.w.outliers"))

## End setup
