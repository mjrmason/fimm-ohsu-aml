suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(synapseClient))

synapseLogin()

## Extract data from supplemental file of
## A machine learning approach to integrate big data for precision medicine in acute myeloid leukemia
## Lee, Celik, Logsdon, et al.
## Nat Communications 2018

## Synapse folder holding Lee et al data
lee.syn.id <- "syn12176376"

## Synapse id of supplemental data 1 file
lee.sd1.file.syn.id <- "syn12176380"
obj <- synGet(id=lee.sd1.file.syn.id, downloadFile = TRUE)
lee.sd1.file <- getFileLocation(obj)

## Clinical and Drug Response Data folder
lee.drug.syn.id <- "syn12176377"

## Expression Data folder
lee.expr.syn.id <- "syn12176378"

## Processed Data foler
lee.processed.data.syn.id <- "syn12176379"

## Extract the _patient_ expression data
expr <- read.xlsx(lee.sd1.file, sheet = 1, startRow = 3)
expr <- expr[!is.na(expr$UID),]
rownames(expr) <- expr$UID
## Drop the sparrow column
expr <- expr[, !(colnames(expr) %in% c("sparrow.hub", "UID"))]

file <- "lee-patient-expr.tsv"
write.table(expr, file = file, row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)

obj <- File(file, parentId = lee.expr.syn.id)
synStore(obj)

## Extract the _patient_ drug response AUCs
aucs <- read.xlsx(lee.sd1.file, sheet = 2)
rownames(aucs) <- aucs$UID
aucs <- aucs[, !(colnames(aucs) %in% c("53.drugs", "UID"))]
colnames(aucs) <- gsub(colnames(aucs), pattern="\\.", replacement="")

file <- "lee-patient-processed-aucs.tsv"
write.table(aucs, file = file, row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)

obj <- File(file, parentId = lee.processed.data.syn.id)
synStore(obj)

## Extract the _patient_ drug response IC50s
## NB: I added a UID column name for the first column
ic50s <- read.xlsx(lee.sd1.file, sheet = 3)
ic50s <- ic50s[, grepl(colnames(ic50s), pattern="AML") | colnames(ic50s) == "UID"]
rownames(ic50s) <- ic50s$UID
ic50s <- ic50s[, !(colnames(ic50s) %in% c("53.drugs", "UID"))]
colnames(ic50s) <- gsub(colnames(ic50s), pattern="\\.", replacement="")

file <- "lee-patient-processed-ic50s.tsv"
ic50s[ic50s == " CCNU"] <- NA
write.table(ic50s, file = file, row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)

obj <- File(file, parentId = lee.processed.data.syn.id)
synStore(obj)

## Extract the _patient_ drug response raw data
dr <- read.xlsx(lee.sd1.file, sheet = 6)

dr <- dr[, !(colnames(dr) %in% c("54.drugs"))]

## The organization of this file is very odd.  Columns are _mostly_ patients.
## However, there are 3 "Compound" and "[Compound].M" columns.  
## Only the first "Compound" column has all non-blank entries.  e.g., the first 8 rows are
## for 5-Iodotubercidin.  The second two "Compound" columns only list 5-Iodotubercidin in the
## first of these rows, and leave the rest of the rows blank.
## Confirm this.
compound.col.indices <- which(colnames(dr) %in% c("Compound", "Compound.1", "Compound.2"))
for(indx in compound.col.indices[compound.col.indices != compound.col.indices[1]]) {
  flag <- !is.na(dr[, indx])
  if(!all(dr[flag, compound.col.indices[1]] == dr[flag, indx])) {
    stop("Mismatch compound names")
  }
}
## Drop all but the first compound column
dr <- dr[, -compound.col.indices[compound.col.indices != compound.col.indices[1]]]

## Note that concentrations are _not_ consistent.
conc.col.indices <- which(grepl(colnames(dr), pattern="\\[Compound\\].M"))
flag <- unlist(apply(dr[, conc.col.indices], 1, 
               function(row) {
                 r <- row[!is.na(row)]
                 if(length(r) < 2) { return(TRUE) }
                 r <- r - r[1]
                 eps <- 10^-14
                 all(abs(r) < eps)
               }))
all.compounds <- unique(dr$Compound)
compounds.with.same.conc.across.pts <- all.compounds[!(all.compounds %in% dr[!flag, "Compound"])]
for(indx in conc.col.indices[conc.col.indices != conc.col.indices[1]]) {
  flag <- !is.na(dr[, indx])
  if(!all(dr[flag, conc.col.indices[1]] == dr[flag, indx])) {
    cat("Mismatch concentrations\n")
    print(head(dr[flag & dr[flag, conc.col.indices[1]] != dr[flag, indx], conc.col.indices]))
  }
}

## Drop all but the first compound concentration column
## dr <- dr[, -conc.col.indices[conc.col.indices != conc.col.indices[1]]]

## Melt the data
library(dplyr)
library(plyr)
conc.col.indices <- which(grepl(colnames(dr), pattern="\\[Compound\\].M"))
indices <- c(conc.col.indices, ncol(dr)+1)
ranges <- lapply(1:(length(indices)-1),
                        function(i) unique(c(1, indices[i]:(indices[i+1]-1))))


library(reshape2)
drs <- llply(ranges, .fun = function(range) {
                              tmp <- dr[, range]
                              colnames(tmp)[2] <- "[Compound].M"
                              melt(tmp, id.vars = c("Compound", "[Compound].M"))
                            })
dr.all <- do.call("rbind", drs)

## Convert to nanomolar
dr.all[, 2] <- dr.all[,2] * 10^9
colnames(dr.all) <- c("Compound", "conc.nM", "sample", "viability")

## Drop "."s from sample names, since they do not occur in expr data
dr.all$sample <- gsub(dr.all$sample, pattern="\\.", replacement="")
samples <- unique(dr.all$sample)

missing <- samples[!(samples %in% colnames(expr))]
cat(paste0("The following patients have drug response data but not expr: ", paste(missing, collapse = ", "), "\n"))

dr.all$inhibition <- 100 - dr.all$viability

file <- "lee-drug-response.tsv"
write.table(dr.all, file = file, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
obj <- File(file, parentId = lee.drug.syn.id)
synStore(obj)


## File Name: Supplementary Data 2
## Description: The 160 drugs in our customized drug panel and their functional classes.

## File Name: Supplementary Data 3
## Description: Clinical information of the 30 AML patient samples.