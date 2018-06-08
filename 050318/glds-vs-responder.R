my.dup <- function(x) duplicated(x, fromLast = TRUE) | duplicated(x, fromLast = FALSE)

## Read in the drug response, patient metadata, and drug metadata.
## Note that these drugs are just those that also occur in FIMM.
ohsu.metadata <- read.table("ohsu-patient-metadata.tsv", header = TRUE, sep = "\t")
ohsu.drug.response <- read.table("ohsu-drug-response.tsv", header = TRUE, sep = "\t")
ohsu.drug.metadata <- read.table("ohsu-drug-metadata.tsv", header = TRUE, sep = "\t")

stop("stop")

## Find patients that have only one type of response across all samples (e.g., _not_ "Complete Response" and "Refractory")
patient.response.tbl <- ohsu.metadata
un <- unique(patient.response.tbl[, c("patient_id", "response")])
un <- un[!my.dup(un$patient_id),]
un <- subset(un, !is.na(response))
un <- subset(un, response %in% c("Complete Response", "Refractory"))

seq.response <- un
rownames(seq.response) <- un$patient_id

## Restrict to common patients in patient metadata and drug response
colnames(ohsu.drug.response) <- gsub(colnames(ohsu.drug.response), pattern="X", replacement="")
common.samples <- intersect(rownames(seq.response), colnames(ohsu.drug.response))

## Calculate GLDS as the mean of z-scored drug response
ohsu.glds <- colMeans(t(scale(t(ohsu.drug.response[, common.samples]))), na.rm=TRUE)
seq.response <- seq.response[common.samples, ]
seq.response$glds <- as.numeric(ohsu.glds)

## Look for association between response and glds.
wt <- wilcox.test(glds ~ response, data = seq.response)
pval <- wt$p.value

## Plot GLDS vs response.
library(ggplot2)
library(ggbeeswarm)
g <- ggplot(seq.response, aes(x = response, y = glds))
g <- g + geom_violin()
g <- g + geom_beeswarm()
g <- g + ylab("GLDS")
g <- g + ggtitle(paste0("Response ~ GLDS: p = ", as.numeric(pval)))
pdf("ohsu-glds-vs-response.pdf")
g
d <- dev.off()

