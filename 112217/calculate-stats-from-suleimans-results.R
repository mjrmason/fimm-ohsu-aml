tbl <- read.table("TR-Ridge.tsv", sep="\t", header=TRUE)

val.col <- "Ridge.Pearson.Correlation"
fdr.cutoff <- 0.2
fdr.col <- "Ridge.fdr"

pos.flag <- !is.na(tbl[, val.col]) & (as.numeric(tbl[, val.col] > 0))

sig.flag <- !is.na(tbl[, fdr.col]) & (as.numeric(tbl[, fdr.col]) < fdr.cutoff)
sig.and.pos.flag <- pos.flag & sig.flag
sig.only.flag <- sig.flag & !pos.flag
pos.only.flag <- !sig.flag & pos.flag
neither.flag <- !sig.flag & !pos.flag

mat <- rbind(c(length(which(sig.and.pos.flag)), length(which(sig.only.flag))), c(length(which(pos.only.flag)), length(which(neither.flag))))
print(sum(mat))

ft <- fisher.test(mat, alternative = "greater")
print(ft)

cat(paste0("Num with positive correlation for ridge model: ", length(which(pos.flag)), "\n"))
cat(paste0("Num with significant correlation for ridge model: ", length(which(sig.flag)), "\n"))
cat(paste0("Num with positive correlation for ridge model: ", length(which(sig.and.pos.flag)), "\n"))
ci <- quantile(tbl[, val.col], na.rm=TRUE, probs = c(0.025, 1 - 0.025))
cat(paste0("Mean ridge correlation = ", mean(tbl[, val.col], na.rm=TRUE), " sd: ", sd(tbl[, val.col], na.rm=TRUE), " 95% CI: ", ci[1], " - ", ci[2], "\n"))


## Do TR
val.col <- "TR.Pearson.Correlation"
fdr.cutoff <- 0.2
fdr.col <- "TR.fdr"

pos.flag <- !is.na(tbl[, val.col]) & (as.numeric(tbl[, val.col] > 0))

sig.flag <- !is.na(tbl[, fdr.col]) & (as.numeric(tbl[, fdr.col]) < fdr.cutoff)
sig.and.pos.flag <- pos.flag & sig.flag
sig.only.flag <- sig.flag & !pos.flag
pos.only.flag <- !sig.flag & pos.flag
neither.flag <- !sig.flag & !pos.flag

mat <- rbind(c(length(which(sig.and.pos.flag)), length(which(sig.only.flag))), c(length(which(pos.only.flag)), length(which(neither.flag))))
print(sum(mat))

ft <- fisher.test(mat, alternative = "greater")
print(ft)

cat(paste0("Num with positive correlation for TR model: ", length(which(pos.flag)), "\n"))
cat(paste0("Num with significant correlation for TR model: ", length(which(sig.flag)), "\n"))
cat(paste0("Num with positive correlation for TR model: ", length(which(sig.and.pos.flag)), "\n"))
ci <- quantile(tbl[, val.col], na.rm=TRUE, probs = c(0.025, 1 - 0.025))
cat(paste0("Mean TR correlation = ", mean(tbl[, val.col], na.rm=TRUE), " sd: ", sd(tbl[, val.col], na.rm=TRUE), " 95% CI: ", ci[1], " - ", ci[2], "\n"))
