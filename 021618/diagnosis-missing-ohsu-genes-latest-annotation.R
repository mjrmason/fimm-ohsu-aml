## Download the genome build against which the OHSU data were aligned
## wget -c ftp://ftp.ensembl.org/pub/release-91/gtf/homo_sapiens/Homo_sapiens.GRCh38.91.gtf.gz
## gunzip Homo_sapiens.GRCh38.91.gtf.gz

## Read in the GTF file of the genome file
library(refGenome)
ens <- ensemblGenome()
build <- "GRCh38.91"
read.gtf(ens, paste0("Homo_sapiens.", build, ".gtf"))
my_gene <- getGenePositions(ens)

## Read in the prioritized genes
prioritized.genes <- read.table("prioritized-genes.tsv", sep="\t", header=TRUE)

library(synapseClient)
synapseLogin()

## Read in the FIMM expression data (RNASeq.CPM.log2.bc.csv)
## Load (batch corrected) FIMM expression data
obj <- synGet(id="syn8270602", downloadFile = TRUE)
file <- getFileLocation(obj)
expr.matrix <- read.table(file, header=TRUE, sep=",")
## Make the columns samples
fimm.expr <- t(expr.matrix)

## Read in the OHSU expression data
synId <- "syn10083723"
obj <- synGet(synId, downloadFile = TRUE)
file <- getFileLocation(obj)
ohsu.expr <- read.table(file, header=TRUE, sep="\t")
## Convert the column names from X20.00347 to 20-00347
colnames(ohsu.expr) <- gsub("X(.*)\\.(.*)", "\\1-\\2", colnames(ohsu.expr))

## Find the prioritized genes that are missing from the OHSU expression data
## and report how the number that are or are not in the genome build
missing.flag <- !(prioritized.genes$ensg %in% rownames(ohsu.expr))

write.table(prioritized.genes[missing.flag,], file=paste0("missing-prioritized-genes", build, ".tsv"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

num.found <- length(which(!missing.flag))
num.missing <- length(which(missing.flag))
cat(paste("OHSU: Of ", nrow(prioritized.genes), " prioritized genes, ", num.found, " are included in the data set and ", num.missing, " are missing\n"))
missing.genes <- prioritized.genes$ensg[missing.flag]
included.genes <- prioritized.genes$ensg[!missing.flag]
missing.genes.in.genome.build <- missing.genes %in% my_gene$gene_id
num.missing.genes.in.genome.build <- length(which(missing.genes.in.genome.build))
num.missing.genes.not.in.genome.build <- length(which(!missing.genes.in.genome.build))
cat(paste("OHSU: Of the ", length(missing.genes), " missing prioritized genes, ", num.missing.genes.in.genome.build, " genes are in the genome build and ",
          num.missing.genes.not.in.genome.build, " genes are not in the genome build\n"))
included.genes.in.genome.build <- included.genes %in% my_gene$gene_id
num.included.genes.in.genome.build <- length(which(included.genes.in.genome.build))
num.included.genes.not.in.genome.build <- length(which(!included.genes.in.genome.build))
cat(paste("OHSU: Of the ", length(included.genes), " included prioritized genes, ", num.included.genes.in.genome.build, " genes are in the genome build and ",
          num.included.genes.not.in.genome.build, " genes are not in the genome build\n"))