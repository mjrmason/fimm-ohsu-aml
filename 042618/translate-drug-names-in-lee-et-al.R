
## Read in Robert Allaway's table of compounds
compound.file.name <- "compound_names.rds"
githubURL <- paste0("https://github.com/Sage-Bionetworks/polypharmacology-db/raw/develop/Data/", compound.file.name)
githubURL <- paste0("https://raw.githubusercontent.com/Sage-Bionetworks/polypharmacology-db/develop/Data/", compound.file.name)
download.file(githubURL, compound.file.name, method="curl")
compound.tbl <- readRDS(compound.file.name)

## Read in the compounds from the Lee et al paper
## A machine learning approach to integrate big data for precision medicine in acute myeloid leukemia
## Lee, Celik, Logsdon, et al.
## Nat Communications 2018
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(synapseClient))
synapseLogin()

## Synapse id of supplemental data 2 file
lee.sd2.file.syn.id <- "syn12179873"
obj <- synGet(id=lee.sd2.file.syn.id, downloadFile = TRUE)
lee.sd2.file <- getFileLocation(obj)
lee.drug.anno.tbl <- read.xlsx(lee.sd2.file, sheet = 1, startRow = 2)
lee.drug.anno.tbl$internal_id <- unlist(lapply(lee.drug.anno.tbl$Name,
                                               function(name) {
                                                 flag <- grepl(compound.tbl$external_id, pattern=name, ignore.case=TRUE)
                                                 if(any(flag)) { return(compound.tbl$internal_id[flag][1]) }
                                                 flag <- grepl(compound.tbl$common_name, pattern=name, ignore.case=TRUE)
                                                 if(any(flag)) { return(compound.tbl$internal_id[flag][1]) }
                                                 return(NA)
                                               }))
lee.drug.anno.tbl$internal_id <- unlist(apply(lee.drug.anno.tbl[, c("Alternate.Name", "internal_id")], 1,
                                               function(row) {
                                                 names <- unlist(strsplit(row[1], split=",[ ]*"))
                                                 id <- row[2]
                                                 if(is.na(row[1])) { return(id) }
                                                 if(!is.na(id)) { return(id) }
                                                 for(name in names) {
                                                   new.id <- NA
                                                   flag <- grepl(compound.tbl$external_id, pattern=name, ignore.case=TRUE)
                                                   if(any(flag)) { new.id <- compound.tbl$internal_id[flag][1] }
                                                   flag <- grepl(compound.tbl$common_name, pattern=name, ignore.case=TRUE)
                                                   if(any(flag)) { new.id <- compound.tbl$internal_id[flag][1] }
                                                   if(!is.na(new.id)) { return(new.id) }
                                                 }
                                                 return(NA)
                                                 if(is.na(new.id)) new.id <- id
                                                 if(!is.na(id) && (id != new.id)) { stop(paste0("Inconsistent ids: ", id, " and ", new.id, " for ", name, "\n")) }
                                                 if(is.na(id)) { id <- new.id }
                                                 id
                                               }))

## Get the OHSU/FIMM drug map
drug.mapping.synId <- "syn11361110" ## fimm-ohsu-drug-map.tsv
obj <- synGet(drug.mapping.synId, downloadFile = TRUE)
drug.map <- read.table(getFileLocation(obj), sep="\t", header=TRUE, quote="\"", comment.char="")
fimm.ohsu.drug.map <- drug.map

cols <- c("OHSU_DRUG_NAME", "FIMM_DRUG_NAME", "Trade.names")
drug.map$internal_id <- unlist(apply(drug.map[, cols], 1, 
                                      function(row) {
                                        for(names in row) {
                                          if(is.na(names)) { next }
                                          names <- unlist(strsplit(names, split=",[ ]*"))
                                          for(name in names) {
                                            new.id <- NA
                                            flag <- grepl(compound.tbl$external_id, pattern=name, ignore.case=TRUE)
                                            if(any(flag)) { new.id <- compound.tbl$internal_id[flag][1] }
                                            flag <- grepl(compound.tbl$common_name, pattern=name, ignore.case=TRUE)
                                            if(any(flag)) { new.id <- compound.tbl$internal_id[flag][1] }
                                            if(!is.na(new.id)) { return(new.id) }
                                          }
                                        }
                                        return(NA)
                                      }))

m <- merge(subset(drug.map, !is.na(internal_id)), subset(lee.drug.anno.tbl, !is.na(internal_id)), by = "internal_id")
write.table(file="ohsu-fimm-lee-drug-map.tsv", m, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)

## Store in Files/Analyses/fimm-ohsu-lee
file <- "ohsu-fimm-lee-drug-map.tsv"
parentId <- "syn12180492"
obj <- File(file, parentId = parentId)
synStore(obj)


fimm-ohsu-sanger-drug-map.tsv

drugs.of.interest <- c("Venetoclax", "TG100-115", "XAV-939", "Lovastatin", "Panobinostat", "Neratinib", "PD184352", "Trametinib", "Palbociclib", "Pelitinib")
drugs.of.interest <- c("XAV-939", "Trametinib", "Venetoclax", "Selumetinib", "Selinexor", "SNS-032", "Panobinostat", "Palbociclib", "PD184352", "Lovastatin", "BI 2536", "AGI-5198",
                       "Navitoclax")
subset(drug.map, FIMM_DRUG_NAME %in% drugs.of.interest)[,1:5]

subset(m, FIMM_DRUG_NAME %in% drugs.of.interest)[,1:5]
subset(m, FIMM_DRUG_NAME %in% drugs.of.interest)[,c("OHSU_DRUG_NAME", "FIMM_DRUG_NAME", "Name", "Mechanism.Targets", "Class.explained")]
          OHSU_DRUG_NAME FIMM_DRUG_NAME         Name      Mechanism.Targets                         Class.explained
13  SNS-032 (BMS-387032)        SNS-032      SNS-032          CDK inhibitor                     B. Kinase inhibitor
26               BI-2536        BI 2536      BI-2536         PLK1 inhibitor                     B. Kinase inhibitor
35          Panobinostat   Panobinostat Panobinostat         HDAC inhibitor E. Differentiating/ epigenetic modifier
37 Selumetinib (AZD6244)    Selumetinib     AZD-6244       MEK1/2 inhibitor                     B. Kinase inhibitor
43           Palbociclib    Palbociclib    PD0332991 CDK inhibitor (Cdk4/6)                     B. Kinase inhibitor

## Read in the fimm-ohsu-drug-map.tsv
drug.mapping.synId <- "syn11361110" ## fimm-ohsu-drug-map.tsv
obj <- synGet(drug.mapping.synId, downloadFile = TRUE)
fimm.ohsu.drug.map <- read.table(getFileLocation(obj), sep="\t", header=TRUE, quote="\"", comment.char="")

## Read in the fimm-ohsu-ctrp-drug-map.tsv (syn12180273)
synId <- "syn12180273"
obj <- synGet(id=synId, downloadFile = TRUE)
file <- getFileLocation(obj)
fimm.ohsu.ctrp.drug.map <- read.table(file, sep="\t", header=TRUE, comment.char="")

## Read in the fimm-ohsu-sanger-drug-map.tsv (syn12180490)
synId <- "syn12180490"
obj <- synGet(id=synId, downloadFile = TRUE)
file <- getFileLocation(obj)
fimm.ohsu.sanger.drug.map <- read.table(file, sep="\t", header=TRUE, comment.char="")

## Read in the ohsu-fimm-lee-drug-map.tsv (syn12180493)
synId <- "syn12180493"
obj <- synGet(id=synId, downloadFile = TRUE)
file <- getFileLocation(obj)
fimm.ohsu.lee.drug.map <- read.table(file, sep="\t", header=TRUE, comment.char="")

fimm.ohsu.drug.map$OHSU <- TRUE
fimm.ohsu.drug.map$FIMM <- TRUE
fimm.ohsu.ctrp.drug.map$CTRP <- TRUE
fimm.ohsu.sanger.drug.map$Sanger <- TRUE
fimm.ohsu.lee.drug.map$Lee <- TRUE
fimm.ohsu.drug.map <- merge(fimm.ohsu.drug.map, fimm.ohsu.ctrp.drug.map, all.x = TRUE)
fimm.ohsu.drug.map <- merge(fimm.ohsu.drug.map, fimm.ohsu.sanger.drug.map, all.x = TRUE)
fimm.ohsu.drug.map <- merge(fimm.ohsu.drug.map, fimm.ohsu.lee.drug.map, all.x = TRUE)

cols <- c("OHSU_DRUG_NAME", "FIMM_DRUG_NAME", "Mechanism.Targets", "Gene.Targets", "Class.explained", "OHSU", "FIMM", "CTRP", "Sanger", "Lee")
fimm.ohsu.drug.map <- fimm.ohsu.drug.map[, c(cols, colnames(fimm.ohsu.drug.map)[!(colnames(fimm.ohsu.drug.map) %in% cols)])]
file <- "fimm-ohsu-ctrp-sanger-lee-drug-map.tsv"
write.table(file = file, fimm.ohsu.drug.map, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

flag <- fimm.ohsu.drug.map$FIMM_DRUG_NAME %in% drugs.of.interest
file <- "fimm-ohsu-ctrp-sanger-lee-drugs-of-interest-map.tsv"
tmp <- fimm.ohsu.drug.map[flag,]
tmp <- tmp[order(tmp$FIMM_DRUG_NAME),]
write.table(file = file, tmp, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

