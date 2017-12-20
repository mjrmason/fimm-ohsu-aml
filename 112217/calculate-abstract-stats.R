
tbl <- read.table("fimm-ohsu-validate-trained-models-112217-one-ohsu-AUC-ridge-pearson-vs-gene-sets-all.xls", sep="\t", header=TRUE, as.is=TRUE)
fdr.cutoff <- 0.2
fdr.col <- "fimm.fdr"
corr.col <- "fimm.val"

tot <- nrow(tbl)
sig.flag <- !is.na(tbl[, fdr.col]) & (as.numeric(tbl[,fdr.col]) < fdr.cutoff) 
pos.flag <- !is.na(tbl[, corr.col]) & (as.numeric(tbl[,corr.col] > 0))
sig.and.pos.flag <- sig.flag & pos.flag
neither.flag <- !sig.flag & !pos.flag
sig.only.flag <- sig.flag & !pos.flag
pos.only.flag <- !sig.flag & pos.flag

mat <- rbind(c(length(which(sig.and.pos.flag)), length(which(sig.only.flag))), c(length(which(pos.only.flag)), length(which(neither.flag))))
fisher.test(mat, alternative = "greater")


