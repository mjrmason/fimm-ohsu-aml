suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("maftools"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("synapseClient"))

synapseLogin()

## Read in the FIMM detailed (as opposed to binarized) genomic data
## to see what kind of mutations are being captured:
## Load in the FIMM detailed (as opposed to binarized) genomic data
synId <- "syn8270591"
obj <- synGet(id=synId, downloadFile = TRUE, downloadLocation = ".")
file <- getFileLocation(obj)
fimm.genomic.detail <- read.table(file, header=TRUE, sep=",", comment.char="", quote="\"")

## Mutation types reported for FIMM
## unique(fimm.genomic.detail$variant.effect)
## [1] non_synonymous_coding             codon_change_plus_codon_deletion 
## [3] codon_insertion                   frame_shift                      
## [5] stop_gained                       codon_change_plus_codon_insertion
## [7] codon_deletion                    splice_site_donor                
## [9] ITD                               splice_site_acceptor             
## [11] start_lost                       

## Load OHSU genomic data (just mutect for now)
origBeatMafs     <- synQuery("SELECT id,name FROM file WHERE parentId=='syn5522834'")
beatMafs     <- origBeatMafs[grep("mutect",origBeatMafs$file.name),]
beatMafNames <- gsub("^.*_(.*?)_AML_.*maf$","\\1",beatMafs$file.name)

## mutation.types.to.total: list of mutation types that will be considered a mutation (i.e., a "1")
##                          in the final 0/1 matrix.
## e.g., c("Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Splice_Site", "Translation_Start_Site")
read.mafs <- function(maf.names, maf.tbl, mutation.types.to.total = c("Missense_Mutation", "Nonsense_Mutation")) {
  ## Read in the MAFs
  maf.list  <- list()
  for(i in 1:nrow(maf.tbl)){ maf.list[[maf.names[i]]] <-  read.maf(getFileLocation(synGet(maf.tbl$file.id[i])), verbose = F); print(i)}
  
  maf.tbl <- ldply(maf.list, getGeneSummary)
  rm(maf.list); gc()
  
  ## Tally mutation.types.to.total (e.g., Missense_Mutation and Nonsense_Mutation)
  print(colnames(maf.tbl))
  maf.tbl$total <- unlist(apply(maf.tbl[,mutation.types.to.total], 1, function(row) sum(row, na.rm=TRUE)))
  print(colnames(maf.tbl))
  maf.tbl <- maf.tbl[,c(".id","Hugo_Symbol","total")]
  
  maf.tbl                  <- spread(maf.tbl, Hugo_Symbol,total); 
  rownames(maf.tbl)        <- maf.tbl$.id; 
  maf.tbl                  <- maf.tbl[,-1]; 
  maf.tbl[is.na(maf.tbl)] <- 0
  maf.tbl[maf.tbl > 0] <- 1
  maf.tbl
}

## The mutation types in the OHSU Mutect MAFs are:
## "Missense_Mutation"     
## "Nonsense_Mutation"     
## "Splice_Site"
## "Nonstop_Mutation"      
## "Translation_Start_Site"

## The mutation types in the FIMM detailed genomic data are:
## [1] non_synonymous_coding             codon_change_plus_codon_deletion 
## [3] codon_insertion                   frame_shift                      
## [5] stop_gained                       codon_change_plus_codon_insertion
## [7] codon_deletion                    splice_site_donor                
## [9] ITD                               splice_site_acceptor             
## [11] start_lost                       

## Let's consider the following types of mutations in OHSU
## "Missense_Mutation", "Nonsense_Mutation", "Splice_Site"

## And these in FIMM
## "non_synonymous_coding", "splice_site_donor", "splice_site_acceptor"             

beat.genomic <- read.mafs(maf.names = beatMafNames, maf.tbl = beatMafs, mutation.types.to.total = c("Missense_Mutation", "Nonsense_Mutation", "Splice_Site"))

## Make the columns correspond to patients
beat.genomic <- t(beat.genomic)
path <- "ohsu.genomic.ns.ss.tsv"
write.table(file = path, beat.genomic, sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

## Store the OHSU data on Synapse under the Beat AML Data/Processed Data folder (syn10083332)
parentId <- "syn10083332"
f <- File(path, parentId = parentId, synapseStore = TRUE)
synStore(f)

## Re-summarize the FIMM data so that it only calls a gene mutated if it has the following mutation types
## "non_synonymous_coding", "splice_site_donor", "splice_site_acceptor"             
fimm.genomic.detail <- subset(fimm.genomic.detail, variant.effect %in% c("non_synonymous_coding", "splice_site_donor", "splice_site_acceptor"))
fimm.genomic.detail <- unique(fimm.genomic.detail[, c("sample", "gene")])
fimm.genomic.detail <- subset(fimm.genomic.detail, gene != "")
fimm.genomic.detail$mutated <- 1

fimm.genomic.bin <- spread(fimm.genomic.detail, gene, mutated, fill = 0)
rownames(fimm.genomic.bin) <- fimm.genomic.bin$sample
fimm.genomic.bin <- fimm.genomic.bin[, !(colnames(fimm.genomic.bin) == "sample")]

## Here the columns correspond to genes.  We'll leave that as is to be consistent with the other FIMM files.
path <- "fimm.genomic.ns.ss.tsv"
write.table(file = path, fimm.genomic.bin, sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

## Store the FIMM data on Synapse under the FIMM Data/Processed Data folder (syn8270577)
parentId <- "syn8270577"
f <- File(path, parentId = parentId, synapseStore = TRUE)
synStore(f)
