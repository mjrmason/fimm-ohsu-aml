library(plyr)
library(ggplot2)
library(reshape2)
library(ggbeeswarm)

base.dir <- "top-feature-plots"
methods <- c("ridge", "lasso", "rf")
gsets <- c("path-go.BP", "path-kegg", "kegg", "gene", "biocarta", "hallmark", "reactome")
responses <- c("random", "mean.feature", "model.mean.feature", "no.transform", "subtract.mean", "unpenalized.mean.feature")

methods <- c("rf", "lasso")
gsets <- c("gene")

methods <- c("rf", "lasso", "ridge")
gsets <- c("path-go.BP", "gene")
gsets <- c("path-go.BP", "path-kegg", "kegg", "gene", "biocarta", "hallmark", "reactome")
## gsets <- c("path-go.BP")

for(method in methods) {
  for(gset in gsets) {
    df.all <- ldply(1:length(responses), .fun = function(i) {
      respi <- responses[i]
      diri <- paste0(base.dir, "/", method, "-", gset, "-", respi, "/")
      df <- ldply(list.files(diri, full.names = TRUE, pattern = "*.tsv"), 
                    .fun = function(filei) {
                             tbli <- read.table(filei, sep="\t", header=TRUE, as.is=TRUE)
                             if(nrow(tbli) < 1) { return(NULL) }
                             if("n.both" %in% colnames(tbli)) {
                               tbli$frac <- as.numeric(tbli$n.both) / as.numeric(tbli$top)
                               tbli$frac.expected <- as.numeric(tbli$n.expected) / as.numeric(tbli$top)
                             } else {
                               tbli$frac <- as.numeric(tbli$frac.overlap)
                               tbli$frac.expected <- as.numeric(tbli$exp.frac.overlap)
                             }
                             tbli$frac.diff <- tbli$frac - tbli$frac.expected
                             if(all(is.na(tbli$frac.diff))) { return(NULL) }
                             data.frame(area = sum(tbli$frac.diff, na.rm=TRUE), resp = respi)
                           }) })
     if(nrow(df.all) == 0) { next }
        
        g <- ggplot(data = df.all, aes(x = resp, y = area))
        g <- g + geom_violin()
##        g <- g + geom_beeswarm()
        g <- g + geom_boxplot(width = 0.1)
        g <- g + ylab("Sum observed - expected fraction overlapping features")
        g <- g + xlab("Response")
        g <- g + ggtitle(paste(method, gset, sep = " "))
        pdf.file <- paste(method, gset, "all-diff-frac-overlap-features.pdf", sep="-")
        pdf(pdf.file)
        print(g)
        d <- dev.off()
  }
}

cat("Done with all\n")

for(method in methods) {
  for(gset in gsets) {
    for(i in 1:(length(responses)-1)) {
      respi <- responses[i]
      diri <- paste0(base.dir, "/", method, "-", gset, "-", respi, "/")
      for(j in (i+1):length(responses)) {
        respj <- responses[j]
        print(c(method, gset, respi, respj))
        df <- ldply(list.files(diri, full.names = TRUE, pattern = "*.tsv"), .parallel = FALSE, 
                    .fun = function(filei) {
                             tbli <- read.table(filei, sep="\t", header=TRUE, as.is=TRUE)
                             if(nrow(tbli) < 1) { return(NULL) }
                             filej <- gsub(filei, pattern=respi, replacement=respj)
                             if(!file.exists(filej)) { return(NULL) }
                             tblj <- read.table(filej, sep="\t", header=TRUE, as.is=TRUE)
                             if(nrow(tblj) < 1) { return(NULL) }
                             if("n.both" %in% colnames(tbli)) {
                               tbli$frac <- as.numeric(tbli$n.both) / as.numeric(tbli$top)
                               tbli$frac.expected <- as.numeric(tbli$n.expected) / as.numeric(tbli$top)
                             } else {
                               tbli$frac <- as.numeric(tbli$frac.overlap)
                               tbli$frac.expected <- as.numeric(tbli$exp.frac.overlap)
                             }
                             tbli$frac.diff <- tbli$frac - tbli$frac.expected
                             if("n.both" %in% colnames(tblj)) {
                               tblj$frac <- as.numeric(tblj$n.both) / as.numeric(tblj$top)
                               tblj$frac.expected <- as.numeric(tblj$n.expected) / as.numeric(tblj$top)
                             } else {
                               tblj$frac <- as.numeric(tblj$frac.overlap)
                               tblj$frac.expected <- as.numeric(tblj$exp.frac.overlap)
                             }
                             tblj$frac.diff <- tblj$frac - tblj$frac.expected
                             tbl <- tryCatch({merge(tbli, tblj, by = "top", suffixes = c(paste0(".", respi), paste0(".", respj)))}, 
                                              error = function(e) { cat(paste0("error with\n", filei, "\n", filej, "\n")) })
                             if(nrow(tbl) == 0) { return(NULL) }
                             coli <- paste0("frac.diff.", respi)
                             colj <- paste0("frac.diff.", respj)
                             ret <- data.frame(areai = sum(tbl[, coli], na.rm=TRUE), areaj = sum(tbl[, colj], na.rm=TRUE), diffij = sum(tbl[,coli] - tbl[,colj], na.rm=TRUE))
                             if(all(is.na(tbl[, coli])) || all(is.na(tbl[, colj]))) { return(NULL) }
                             ret
                           })
        if(nrow(df) == 0) { next }
        pval <- tryCatch({t.test(df$areai, df$areaj, alternative = "two.sided", paired = TRUE)$p.value}, error = function(e) { return(NA) })
        dfm <- df[, c("areai", "areaj")]
        colnames(dfm) <- c(respi, respj)
        dfm <- melt(dfm)
        g <- ggplot(data = dfm, aes(x = variable, y = value))
        g <- g + geom_violin()
##        g <- g + geom_beeswarm()
        g <- g + geom_boxplot(width = 0.1)
        g <- g + ylab("Sum observed - expected fraction overlapping features")
        g <- g + xlab("Response")
        g <- g + ggtitle(paste(method, gset, "pval:", pval, sep = " "))
        pdf.file <- paste(method, gset, respi, "vs", respj, "diff-frac-overlap-features.pdf", sep="-")
        pdf(pdf.file)
        print(g)
        d <- dev.off()
      }
    }
  }
}

