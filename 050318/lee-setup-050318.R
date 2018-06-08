## Begin setup (this will need to be changed as the data sets in the intersection change)

conc.cols <- list("lee" = "conc.nM")
response.cols <- list("lee" = "inhibition")
responses.are.viabilities <- list("lee" = FALSE)

## Data sets to be analyzed:
data.sets <- list("lee" = "lee")
data.sets.to.analyze <- unlist(data.sets)

## Check these variables that were just changed to lists
## data.set.drug.fit.synIds
## data.set.expr.synIds
## drug.map.drug.id.cols 
## data.set.drug.id.cols
## data.sets

raw.drug.synIds <- list("lee" = "syn12225291")

data.set.drug.fit.synIds <- list("lee" = "syn12230369")
data.set.drug.fit.files <- c("lee.l4.and.ll4.fits.050318.tsv")

drug.mapping.synId <- "syn12180493" ## ohsu-fimm-lee-drug-map.tsv
drug.map.drug.id.cols <- list("ohsu" = "OHSU_DRUG_NAME", "fimm" = "FIMM_Batch_ID_Drug", "lee" = "Name")


## The columns within each respective data set that identifies the drug (in the name space of that data set)
data.set.drug.id.cols <- list("lee" = "Compound")

## Expression data
data.set.expr.synIds <- list("lee" = "syn12176399")
expr.patient.id.cols <- list("lee" = "sample")

## The columns within each respective data set that identifies a screen (i.e., a particular replicate of a patient/cell line x drug)
screen.id.cols <- list("lee" = c("Compound", "sample"))
experimental.factor.cols <- list("lee" = c("Compound", "sample"))

## Store in
## FIMM_BEAT AML Collaboration/Files/FIMM Data/Processed Data
processed.folder.synIds <- list("lee" = "syn12176379")

## Columns whose values are not dependent on the (L.4 vs LL.4) fit.  
merge.cols <- list(
  "lee" = c(screen.id.cols[["lee"]], "min.conc", "max.conc", "all.concs.nM", "uniq.concs.nM", "all.responses",
             "min.conc.w.outliers.nM", "max.conc.w.outliers.nM", "all.concs.w.outliers.nM", "uniq.concs.w.outliers.nM", "all.responses.w.outliers"))

## End setup
