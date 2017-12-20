## Fit drug response curves.
source("fit-ohsu-drc-curves.R")
source("fit-fimm-drc-curves.R")

## Compute DSS based on above fits (we also we have already harmonized the drugs).
source("calculate-dss.R")

## Check the normality of the DSS/AUC/IC50 responses.
source("check-normality-of-drug-response.R")

## Fit FIMM and OHSU independently using ridge, lasso, and RF.
## Check to see what kind of overlap we have.
source("find-model-overlap.R")