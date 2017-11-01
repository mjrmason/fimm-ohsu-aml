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

cat("Use DRC to fit L.4 (4-parameter logistic) and LL.4 (4-paramter log logistic) functions to NTAP data.\n")

## Read in the drug sensitivity data, which are organized as one file per cell line in the 
## Single Agent Screens folder (syn5522627) of the NTAP Drug Screening of pNF Cell Lines project (syn4939906).
parentId <- "syn5522627"
tbl <- synQuery(paste0("SELECT name,id FROM file WHERE parentId==\'", parentId, "\'"))

## The csv files are the data files
tbl <- subset(tbl, grepl(pattern="csv", x=file.name))

## The columns are not consistently spelled or named across the files.  But the raw data
## columns (DATAn and Cn for n=1 to 10) seem to be OK.  Let's define all of the columns that
## are in common.  Also, define common annotations.
annos <- llply(tbl$file.id, .parallel = FALSE,
                        .fun = function(synId) {
                          obj <- synGet(synId, downloadFile = TRUE)
                          names(synGetAnnotations(obj))
                        })
common.annos <- Reduce("intersect", annos)

ntap.drug.data <- lapply(1:nrow(tbl), 
                        function(i) {
                          synId <- tbl$file.id[i]
                          file.name <- tbl$file.name[i]
                          obj <- synGet(synId, downloadFile = TRUE)
                          annos <- synGetAnnotations(obj)
                          csv <- read.table(getFileLocation(obj), sep=",", header=TRUE, as.is=TRUE, comment.char="", quote="\"")
                          ## Cell.line is sometimes spelled like Cell.Line (cap L).  Capitalize all cols so we can merge.
                          colnames(csv) <- toupper(colnames(csv))
                          ## Sometimes the cell line is not even included in the file.  Get it from the annotations.
                          if("CELL.LINE" %in% colnames(csv)) {
                            cell.line.data <- as.character(csv$CELL.LINE[1])
                            cell.line.anno <- as.character(annos["sampleIdentifier"])
                            if(cell.line.data == cell.line.anno) {
                              cat(paste0("Data file and annotations agree on cell line: ", cell.line.data, "\n"))
                            } else {
                              cat(paste0("Data file (", cell.line.data, ") and annotations (", cell.line.anno, ") do not agree.\n"))
                            }
                          }
                          ## Sometimes concentration columns are called CONCn, sometimes Cn.
                          colnames(csv) <- gsub(pattern="CONC", replacement="C", colnames(csv))
                          colnames(csv) <- gsub(pattern="SMILES", replacement="SMI", colnames(csv))
                          colnames(csv) <- gsub(pattern="NCGC.SID", replacement="SID", colnames(csv))
                          flag <- colnames(csv) == "NCGCID"
                          if(any(flag)) { colnames(csv)[flag] <- "SID" }
                          for(anno in names(annos)) {
                            csv[, anno] <- unname(annos[anno])
                          }
                          csv$file.name <- file.name         
                          csv
                        })

cols <- llply(ntap.drug.data, .fun = function(tbl) colnames(tbl))
common.cols <- Reduce("intersect", cols)
ntap.drug.data <- llply(ntap.drug.data, .fun = function(tbl) tbl[, common.cols])
ntap.drug.data <- do.call("rbind", ntap.drug.data)
ntap.drug.data$NAME <- gsub(ntap.drug.data$NAME, pattern="\n", replacement=" ")

get_table_df <- function(table_id) {
    syn_table_data <- synTableQuery(sprintf("select * from %s", table_id))
    return(syn_table_data@values)
}

## Get the map of all NTAP samples to ensure we have the proper cell line identifiers
synId <- "syn8397154"
ntap.sample.key <- get_table_df(synId)

## Check manually/visually that the sampleIdentifier seem sensible (not necessarily the same)
## as the file name
library(knitr)
cat("Check visually that the sampleIdentifier seem sensible (not necessarily the same) given the file name\n")
print(knitr::kable(unique(ntap.drug.data[, c("sampleIdentifier", "file.name")]), row.names=FALSE))
## write.table(unique(ntap.drug.data[, c("sampleIdentifier", "file.name")]), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

## Merge the sampleIdentifier annotation extracted from a drug file to the corresponding Sample Name in the table.
## Ensure that the drug file from which sampleIdentifier was extracted is the same as that listed in the table.
sample.file.tbl <- unique(ntap.drug.data[, c("sampleIdentifier", "file.name")])
sample.file.tbl <- merge(sample.file.tbl, ntap.sample.key[!is.na(ntap.sample.key$`Drug Sensitivity Data`), c("Sample Name", "Drug Sensitivity Data")], 
                         by.x = "sampleIdentifier", by.y = "Sample Name", all.x=TRUE)

sample.file.tbl$syn.file.name <- unlist(lapply(sample.file.tbl$`Drug Sensitivity Data`,
                                        function(synId) {
                                                 if(is.na(synId)) { return(NA) }
                                                 obj <- synGet(synId, downloadFile = FALSE)
                                                 properties(obj)$name
                                        }))

flag <- sample.file.tbl$syn.file.name != sample.file.tbl$file.name
if(any(flag)) {
  cat("The following are discrepancies between the file name from which a sampleIdentifier was extracted and\n")
  cat("the file listed in the table with the matching Sample Name\n")
  write.table(sample.file.tbl[flag,,drop=F])
} else {
  cat("No discrepancies between the file name from which a sampleIdentifier was extracted and\n")
  cat("the file listed in the table with the matching Sample Name\n")
}

data.cols <- common.cols[grepl(pattern="DATA", common.cols)]
conc.cols <- gsub(data.cols, pattern="DATA", replacement="C")
drug.screen.cols <- c("SID", "NAME", "SMI", "sampleIdentifier")

## Reformat the ntap data such that there is one concentration/response per row.
## NB: concentrations are in uM (as described in qhts-protocol-dump-headers.txt -- syn5522652)
## NB: qhts-protocol-dump-headers.txt also saves that the response in normalized relative to DMSO control (i.e., 100% == DMSO)
##     hence, response = viability
tmp <- ntap.drug.data[, c(data.cols, conc.cols, drug.screen.cols)]
ntap.drug.data.long <- ldply(1:nrow(tmp), .parallel = FALSE,
                             .fun = function(i) {
                               dfs <- lapply(data.cols,
                                                    function(data.col) {
                                                      conc.col <- gsub(data.col, pattern="DATA", replacement="C")
                                                      ret <- tmp[i, c(drug.screen.cols, data.col, conc.col)]
                                                      colnames(ret) <- c(drug.screen.cols, "percent_viability", "conc.uM")
                                                      ret
                                                    })
                               do.call("rbind", dfs)
                             })

##indx <- sample.int(nrow(tmp), size = 1)
##cl <- tmp$sampleIdentifier[indx]
##smi <- tmp$SMI[indx]
##dat <- subset(ntap.drug.data.long, ( sampleIdentifier == cl ) & ( SMI == smi ))
##plot(log(dat$conc.uM), dat$percent_viability)

## Scale concentrations from micro- to nano-molar
ntap.drug.data.long$conc.nM <- ntap.drug.data.long$conc.uM * 10^3 

source("../common/dss.R")
source("../common/drc-fit.R")

conc.col <- "conc.nM"
response.col <- "percent_viability"
all.cols <- c(drug.screen.cols, conc.col, response.col)

cat("Fitting NTAP using LL4\n")
ll4.fits <- fit.drc(ntap.drug.data.long[,all.cols], fct = "LL.4", drug.screen.cols = drug.screen.cols, conc.nM.col = conc.col, response.col = response.col, response.is.viability = TRUE)
cat("Done fitting NTAP using LL4\n")
save.image(".Rdata")

if(nrow(ll4.fits) != nrow(ntap.drug.data)) {
  stop(paste0("Number of LL4 fits (", nrow(ll4.fits), ") != Number of original entries (", nrow(ntap.drug.data), ")\n"))
}

cat("Fitting NTAP using L4\n")
l4.fits <- fit.drc(ntap.drug.data.long[,all.cols], fct = "L.4", drug.screen.cols = drug.screen.cols, conc.nM.col = conc.col, response.col = response.col, response.is.viability = TRUE)
cat("Done fitting NTAP using L4\n")
save.image(".Rdata")

if(nrow(l4.fits) != nrow(ntap.drug.data)) {
  stop(paste0("Number of L4 fits (", nrow(l4.fits), ") != Number of original entries (", nrow(ntap.drug.data), ")\n"))
}

common.cols <- c(drug.screen.cols, "min.conc.nM", "max.conc.nM", "all.concs.nM", "uniq.concs.nM")
fits <- merge(ll4.fits, l4.fits, by = common.cols, suffixes = c(".ll4", ".l4"), all = TRUE)

## Store in
## Files/Analysis/Drug Sensitivity/NCATS + CTRP Drug Response Modeling (syn11244429)
parentId <- "syn11244429"

path <- paste0("ntap.l4.and.ll4.fits.tsv")
write.table(file=path, fits, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
f <- File(path, parentId = parentId, synapseStore = TRUE)
synStore(f, executed = NULL, forceVersion = FALSE)

cat("Fit NTAP data using L.4 and LL.4\n")

q()

ntap.cell.lines <- unique(ntap.drug.data$sampleIdentifier)

missing.samples <- ntap.cell.lines[!(ntap.cell.lines %in% ntap.sample.key$`Sample Name`)]
if(!is.null(missing.samples) && (length(missing.samples) > 0)) {
  stop(paste0("The following samples have drug data, but are not listed in the sample table: ", paste(missing.samples, collapse=","), "\n"))
}
expected.missing.samples <- c("ipNF95.11bC_T", "ipNF95.11C")
unexpected.missing.samples <- missing.samples[!(missing.samples %in% expected.missing.samples)]
if(length(unexpected.missing.samples) > 0) {
  stop(paste0("The following samples have drug data, but are not listed in the sample table but are _expected_ to be in that table: ", 
     paste(unexpected.missing.samples, collapse=","), "\n"))
}


