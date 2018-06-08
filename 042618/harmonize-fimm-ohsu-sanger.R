suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(synapseClient))
synapseLogin()

drug.mapping.synId <- "syn11361110" ## fimm-ohsu-drug-map.tsv
obj <- synGet(drug.mapping.synId, downloadFile = TRUE)
drug.map <- read.table(getFileLocation(obj), sep="\t", header=TRUE, quote="\"", comment.char="")
fimm.ohsu.drug.map <- drug.map

if(FALSE) {
library(PharmacoGX)
gsdc.2013.pset <- downloadPSet("GDSC_2013")
gsdc.pset <- downloadPSet("GDSC")
gsdc1000.pset <- downloadPSet("GDSC1000")
fimm.pset <- downloadPSet("FIMM")

## Most FIMM/OHSU drugs are not in PharmacoGX
fimm.ohsu.drug.map$FIMM_DRUG_NAME %in% drugNames(fimm.pset)
}


## Read in Robert Allaway's table of compounds
compound.file.name <- "compound_names.rds"
githubURL <- paste0("https://github.com/Sage-Bionetworks/polypharmacology-db/raw/develop/Data/", compound.file.name)
githubURL <- paste0("https://raw.githubusercontent.com/Sage-Bionetworks/polypharmacology-db/develop/Data/", compound.file.name)
download.file(githubURL, compound.file.name, method="curl")
compound.tbl <- readRDS(compound.file.name)

synId <- "syn11275819"
obj <- synGet(synId, downloadFile = TRUE)
sanger.compounds <- openxlsx:::read.xlsx(getFileLocation(obj), sheet=1)

cols <- c("OHSU_DRUG_NAME", "FIMM_DRUG_NAME", "Trade.names")
fimm.ohsu.drug.map$internal_id <- unlist(apply(fimm.ohsu.drug.map[, cols], 1, 
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

cols <- c("Drug.Name", "Synonyms")
flag <- sanger.compounds$Synonyms == "-"
sanger.compounds$Synonyms[flag] <- NA

sanger.compounds$internal_id <- unlist(apply(sanger.compounds[, cols], 1, 
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

stop("stop")
m <- merge(subset(fimm.ohsu.drug.map, !is.na(internal_id)), subset(sanger.compounds, !is.na(internal_id)), by = "internal_id")
write.table(file="ohsu-fimm-sanger-drug-map.tsv", m, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)

## Save the mapping of Sanger to FIMM/OHSU compounds in the folder
## FIMM_BEAT AML Collaboration/Files/Analyses/fimm-ohsu-sanger (syn12180269)
## NB: there are multiple entries for some of these drugs -- it appears that Sanger uses unique IDs for different concentration ranges
parentId <- "syn12180269"
path <- paste0("fimm-ohsu-sanger-drug-map.tsv")
write.table(file=path, m, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
f <- File(path, parentId = parentId, synapseStore = TRUE)
synStore(f, executed = NULL, forceVersion = FALSE)


