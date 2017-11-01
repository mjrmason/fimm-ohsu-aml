library(plyr)
library(openxlsx)
library(synapseClient)

synapseLogin()

## Load the map of drugs shared between OHSU and FIMM (FIMM_OHSU_Drug_Concentrations.xlsx)
synId <- "syn10163669"
obj <- synGet(synId, downloadFile = TRUE)
ohsu.fimm.drugs <- openxlsx:::read.xlsx(getFileLocation(obj), sheet=1)

## Note that ranges differ _within_ data sets.  It appears that the file
## FIMM_OHSU_Drug_Concentrations.xlsx was created by union'ing the concentrations of a drug
## across all patients and then finding the min and max.  That does not imply that all
## patients were tested for the drug in that range.  Instead, let's limit to the intersection,
## (i.e., max of mins and min of maxs), just as we do across data sets.

## First, note the occurrence of this in the FIMM raw data.
synId <- "syn8488824"
obj <- synGet(id=synId, downloadFile = TRUE)
file <- getFileLocation(obj)
fimm.raw.dr <- read.xlsx(file)
fimm.raw.drug.ranges <- ddply(fimm.raw.dr, .variables = c("PATIENT_ID", "DRUG_ID"),
                                           .fun = function(df) { 
                                                    v <- c(min(df$CONCENTRATION), max(df$CONCENTRATION))
                                                    names(v) <- c("min", "max")
                                                    v 
                                                  })

## Now look at the range of one drug in particular, as generated from the raw data.  
## Note some go from 1.0 to 10,000 and some from 0.1 to 1,000
drug <- "FIMM023833"
cat(paste0("\n\n(Subset of) FIMM drug range for ", drug, " calculated from raw data on a per-patient basis\n"))
cat("Note that for some patients drug is in range 1.0 to 10,000 and for others it is 0.1 to 1,000\n\n")
tail(subset(fimm.raw.drug.ranges, DRUG_ID == drug))

## Now look at the range of that same drug in the annotation table.
## Note that the min and max range for FIMM is 0.1 to 10,000--i.e., the union of the ranges
## across all patients tested for this drug.
cat(paste0("\n\nFIMM drug range for ", drug, " as annotated in FIMM_OHSU_Drug_Concentrations.xlsx\n"))
cat("Note that the range is 0.1 to 10,000--i.e., the union, not the intersection, of the per-patient ranges\n\n")
subset(ohsu.fimm.drugs, FIMM_Batch_ID_Drug == drug)
