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

cat("Harmonize NTAP and CTRP compounds based on their SMILES.\n")

## Robert Allaway's code to compare SMILES

## To run rcdk, needed to update Java JDK as described here:
## http://www.webupd8.org/2012/09/install-oracle-java-8-in-ubuntu-via-ppa.html
## Installed rJava from source as described here:
## install.packages("rJava", type = "source")
## Installed rcdk as described here:
## https://github.com/rajarshi/cdkr
## git clone https://github.com/rajarshi/cdkr.git
##   R CMD build rcdklibs
##   R CMD INSTALL rcdklibs_*gz
##   cd rcdkjar
##   ant clean jar
##   cd ../
##   R CMD build rcdk
##   R CMD INSTALL rcdk_*gz
library(dplyr)
library(rJava)
library(rcdk)
library(fingerprint)
library(parallel)

## input SMILES to fingerprint function
parseInputFingerprint <- function(input) {
  input.mol <- parse.smiles(input)
  fp.inp <- lapply(input.mol, get.fingerprint, type = "extended")
}

parseSingleInputFingerprint <- function(input) {
  input.mol <- NULL
  tryCatch({
    input.mol <- parse.smiles(input)
  }, error = function(e) { return(NULL) })
  fp.inp <- NULL
  tryCatch({
    fp.inp <- get.fingerprint(input.mol[[1]], type = "extended")
  }, error = function(e) { return(NULL) })
  fp.inp
}

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

ntap.drug.data <- llply(tbl$file.id, .parallel = FALSE,
                        .fun = function(synId) {
                          obj <- synGet(synId, downloadFile = TRUE)
                          csv <- read.table(getFileLocation(obj), sep=",", header=TRUE, as.is=TRUE, comment.char="", quote="\"")
                          ## Cell.line is sometimes spelled like Cell.Line (cap L).  Capitalize all cols so we can merge.
                          colnames(csv) <- toupper(colnames(csv))
                          ## Sometimes concentration columns are called CONCn, sometimes Cn.
                          colnames(csv) <- gsub(pattern="CONC", replacement="C", colnames(csv))
                          colnames(csv) <- gsub(pattern="SMILES", replacement="SMI", colnames(csv))
                          colnames(csv) <- gsub(pattern="NCGC.SID", replacement="SID", colnames(csv))
                          flag <- colnames(csv) == "NCGCID"
                          if(any(flag)) { colnames(csv)[flag] <- "SID" }
                          ## Sometimes the cell line is not even included in the file.  Get it from the annotations.
                          annos <- synGetAnnotations(obj)
                          for(anno in names(annos)) {
                            csv[, anno] <- unname(annos[anno])
                          }         
                          csv
                        })

cols <- llply(ntap.drug.data, .fun = function(tbl) colnames(tbl))
common.cols <- Reduce("intersect", cols)
ntap.drug.data <- llply(ntap.drug.data, .fun = function(tbl) tbl[, common.cols])
ntap.drug.data <- do.call("rbind", ntap.drug.data)
ntap.drug.data$NAME <- gsub(ntap.drug.data$NAME, pattern="\n", replacement=" ")

## Load the CTRP drugs
## CTRP folder on Synapse: syn5622707
## CTRP drug info: "v20.meta.per_compound.txt" (syn5632193)
synId <- "syn5632193"
obj <- synGet(synId, downloadFile = TRUE)
ctrp.compounds <- read.table(getFileLocation(obj), sep="\t", header=TRUE, as.is=TRUE, comment.char="", quote="")
ctrp.compounds <- ctrp.compounds[, c("cpd_name", "cpd_smiles", "master_cpd_id")]

## NTAP does not even have unique names for the same SMILE, so don't keep the name around
## ntap.compounds <- unique(ntap.drug.data[, c("NAME", "SMI")])
ntap.compounds <- unique(ntap.drug.data[, c("SMI"), drop=FALSE])

if(any(duplicated(ctrp.compounds$cpd_smiles)) | any(duplicated(ctrp.compounds$cpd_name)) | any(duplicated(ctrp.compounds$master_cpd_id))) {
  stop("CTRP has some duplicated drug identifiers\n")
}

##overlap <- intersect(ctrp.compounds$cpd_name, ntap.compounds$NAME)

##ctrp.overlap.compounds <- subset(ctrp.compounds, cpd_name %in% overlap)
##ntap.overlap.compounds <- subset(ntap.compounds, NAME %in% overlap)
##ctrp.overlap.compounds <- ctrp.overlap.compounds[order(ctrp.overlap.compounds$cpd_name),]
##ntap.overlap.compounds <- ntap.overlap.compounds[order(ntap.overlap.compounds$NAME),]

##fp.ctrp <- parseInputFingerprint(as.character(ctrp.overlap.compounds$cpd_smiles))
##fp.ntap <- parseInputFingerprint(as.character(ntap.overlap.compounds$SMI))

cat("Calculating drug fingerprints for CTRP\n")
fp.ctrp <- parseInputFingerprint(as.character(ctrp.compounds$cpd_smiles))

## Several of the NTAP smiles gives us problems and returns NULL
vec <- as.character(ntap.compounds$SMI)
names(vec) <- vec
cat("Calculating drug fingerprints for NTAP\n")
fp.ntap <- llply(vec, .parallel = FALSE, .fun = parseSingleInputFingerprint)
cat("Filtering drug fingerprints for NTAP\n")
flag <- unlist(lapply(fp.ntap, is.null))
fp.ntap <- fp.ntap[!flag]

## compares every smiles on list one across list 2
cat("Calculating all pairwise distances between CTRP drug fingerprints and NTAP drug fingerprints\n")
sims <- mclapply(fp.ctrp, function(i) {
  sim <- sapply(fp.ntap, function(j) {
    distance(i, j)
  })
  bar <- data.frame(sim=sim, match=names(sim))
  bar
}, mc.cores = detectCores())

cutoff <- 1.0
sims.cutoff <- llply(sims, .parallel = FALSE, 
                     .fun = function(entry) {
                       flag <- entry$sim >= cutoff
                       if(any(flag)) { return(entry$match[flag]) }
                       return(NULL)
                     })
flag <- unlist(lapply(sims.cutoff, is.null))
sims.cutoff <- sims.cutoff[!flag]
flag <- unlist(lapply(sims.cutoff, function(lst) length(lst) != 1))
sims.cutoff <- sims.cutoff[!flag]

ctrp.ntap.mapping <- ldply(1:length(sims.cutoff),
                           .fun = function(i) {
                             df <- data.frame(ntap.smiles = unlist(sims.cutoff[[i]]))
                             df$ctrp.smiles <- names(sims.cutoff)[i]
                             df
                           })

colnames(ctrp.compounds) <- c("ctrp.name", "ctrp.smiles", "ctrp.id")
colnames(ntap.compounds) <- c("ntap.smiles")
ctrp.ntap.mapping <- merge(ctrp.ntap.mapping, ctrp.compounds, by.x = "ctrp.smiles", by.y = "ctrp.smiles", all.x = TRUE)
ctrp.ntap.mapping <- merge(ctrp.ntap.mapping, ntap.compounds, by.x = "ntap.smiles", by.y = "ntap.smiles", all.x = TRUE)

## Save the mapping of CTRP to NTAP compounds in the folder
## Files/Analysis/Drug Sensitivity/NCATS + CTRP Drug Response Modeling (syn11244429)
parentId <- "syn11244429"
path <- paste0("ctrp-ntap-compound-mapping.tsv")
write.table(file=path, ctrp.ntap.mapping, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
f <- File(path, parentId = parentId, synapseStore = TRUE)
synStore(f, executed = NULL, forceVersion = FALSE)
