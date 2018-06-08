## Make violin plots of correlations (model prediction of drug response vs actual drug response) faceted by
## training/test data set and type of feature (gene, biocart, hallmark, kegg, and reactome).
## Make violin plots of correlations (model prediction of drug response vs actual drug response) faceted by
## training/test data set and type of feature (gene, biocart, hallmark, kegg, and reactome).
tbl$gene.set <- factor(tbl$gene.set, levels = c("gene", "biocarta", "hallmark", "kegg", "reactome"))
tbl$train.set <- factor(tbl$train.set, levels = c("ctrp", "ctrp.all", "ctrp.heme", "ctrp.aml", "ohsu", "ntap", "fimm"))
tbl$pt.response <- factor(tbl$pt.response, levels = c("random", "random.class", "no.transform", "no.transform.class", "subtract.mean", "mean.feature"))
tbl$model <- factor(tbl$model, levels = c("lasso", "ridge", "rf"))
for(resp in unique(tbl$response)) {
  for(met in c("pearson", "spearman")) {
##    for(gset in unique(tbl$gene.set)) { }
##    for(mdl in c("ridge", "lasso", "rf")) { }
      sub <- subset(tbl, metric == met & response == resp & train.set == "ohsu" & test.set == "fimm")
      g <- ggplot(data = sub, aes(x = model, y = val))
      g <- g + facet_grid(gene.set ~ pt.response)
      g <- g + geom_violin()
##    g <- g + geom_beeswarm()
      g <- g + geom_boxplot(width = 0.5)
      g <- g + geom_hline(yintercept = 0, linetype = "dashed")
      g <- g + ylab(paste0(capwords(met), " Correlation"))
      g <- g + ggtitle(paste0(resp, " ", capwords(met), " "))
      g <- g + ylim(c(-1,1))
      file <- paste0(file.prefix, "-", resp, "-", met, "-vs-gene-sets.pdf")
      pdf(file)
      print(g)
      d <- dev.off()
##    glist[[length(glist)+1]] <- g
##    { }
  }
##  do.call("grid.arrange", glist)
}

## Make volcano plots of correlations (model prediction of drug response vs actual drug response) faceted by
## training/test data set and type of feature (gene, biocart, hallmark, kegg, and reactome).
for(resp in unique(tbl$response)) {
  for(pt.resp in unique(tbl$pt.response)) {
  for(met in c("pearson", "spearman")) {
    for(mdl in c("ridge", "lasso")) {
      file <- paste0(file.prefix, "-", pt.resp, "-", resp, "-", mdl, "-", met, "-vs-gene-sets-volcano.pdf")
      pdf(file)
      for(gs in unique(tbl$gene.set)) {
        sub <- subset(tbl, metric == met & model == mdl & pt.response == pt.resp & response == resp & gene.set == gs)
        tmp <- sub
        pvalue <- rep(">= 0.05", length(tmp$pval))
        pvalue[!is.na(tmp$pval) & (tmp$pval < 0.05)] <- "< 0.05"
        tmp$pvalue <- pvalue
        tmp$pval <- -log10(tmp$pval)
        ## Training on OHSU and testing on OHSU leads to very low pvalues--artificially
        ## adjust them so that they do not set the scale for the other comparisons.
##        flag <- sub$train.set == "ohsu.train" & sub$test.set == "ohsu.test"
##        if(any(flag)) { tmp$pval[flag] <- tmp$pval[flag] / max(tmp$pval[flag], na.rm=TRUE) }
##        if(any(flag)) { tmp$pval[flag] <- tmp$pval[flag] / 10 }
        g <- ggplot(data = tmp, aes(x = val, y = pval, colour = pvalue))
        g <- g + facet_grid(train.set ~ test.set)
        g <- g + geom_point()
        g <- g + xlab(paste0(capwords(met), " Correlation"))
        g <- g + geom_hline(yintercept = -log10(0.05), linetype = "dashed")
        g <- g + geom_vline(xintercept = 0, linetype = "dashed")
        g <- g + ylab("-Log10 p-value")
        flag <- !is.na(tmp$pval) & (tmp$pval > -log10(0.05))
        add.labels <- FALSE
        if(add.labels && any(flag)) {
          labels <- tmp[flag,]
          g <- g + geom_text(data = labels, aes(x = val, y = pval, label = train.drug))
        }
        g <- g + ggtitle(paste0(resp, " ", capwords(met), " ", mdl, " ", gs))
        print(g)
      }
      d <- dev.off()
##    glist[[length(glist)+1]] <- g
    }
  }
  }
##  do.call("grid.arrange", glist)
}

## Make a scatter plot of no.transform vs no.transform.class correlations.
## i.e., if a drug x is in drug class y compare the correlation of drug class y to that of x alone
drug.class.cor <- subset(tbl, pt.response == "no.transform.class")
drug.cor <- subset(tbl, pt.response == "no.transform")

indices <- 1:nrow(drug.class.cor)
drug.class.cor <- 
  ldply(indices, .parallel = FALSE,
      .fun = function(i) {
               df <- drug.class.cor[i, ,drop=FALSE]
               drugs <- unlist(strsplit(df$train.drug, split=",[ ]*"))
               df <- df[rep(1, length(drugs)),]
               df$train.drug <- drugs
               df
      }) 

m <- merge(drug.class.cor, drug.cor, by = c("train.drug", "alpha", "model", "s", "metric", "train.set", "test.set", "response"), suffixes = c(".class", ".drug"))
tmp <- subset(m, model == "ridge" & metric == "pearson")
plot(tmp$val.class, tmp$val.drug)
abline(a=0, b=1)


for(mdl in unique(m$model)) {
  for(met in c("pearson", "spearman")) {
    tmp <- subset(m, model == mdl & metric == met)
    plot(tmp$val.class, tmp$val.drug)
    abline(a=0, b=1)
    print(c(mdl, met))
    print(table(tmp$val.class > tmp$val.drug))
  }
}

## Make scatter plot of actual vs predicted for each class, coloring by drug
res <- res.aucs[["no.transform.class"]][["all.comparisons"]][["ohsu"]][["fimm"]][["gene"]]
for(i in 1:length(res)) {
  for(j in 1:length(res[[i]])) {
    sub <- res[[i]][[j]]
    train.drug <- sub$train.drug
    model <- sub$model
    alpha <- sub$alpha
    s <- sub$s
    tmp <- extract.predicted.actual(res.aucs[["no.transform.class"]][["all.comparisons"]][["ohsu"]][["fimm"]][["gene"]],
                                    train.drug = train.drug, alpha = alpha, model = model, s = s)

    if((model == "glmnet") && !is.na(alpha) && (alpha == 0)) {
      model <- "ridge"
    }
    if((model == "glmnet") && !is.na(alpha) && (alpha == 1)) {
      model <- "lasso"
    }
    predicted <- tmp$predicted
    actual <- tmp$actual
    drug.names <- sub$drug.names
    drug.indicators <- sub$drug.indicators
    str <- gsub(pattern=",", x=train.drug, replacement="-")
    str <- gsub(pattern=",", x=str, replacement="-")
if(FALSE) {
    pdf(paste0("cluster-model-", str, "-actual-vs-predicted-", "gene", "-", model, ".pdf"))
    g1 <- plot.r2(predicted, actual, colors = drug.names[drug.indicators])
    g1 <- g1 + ggtitle(paste0(train.drug, " (", model, ")"))
    g1 <- g1 + ylab("Actual Response")
    g1 <- g1 + xlab("Predicted Response")
    print(g1)
    d <- dev.off()
}
  }
}

res <- res.aucs[["no.transform.class"]][["all.comparisons"]][["ohsu"]][["fimm"]][["gene"]]

for(train.set in names(res.aucs[["no.transform.class"]][["all.comparisons"]])) {
  for(test.set in names(res.aucs[["no.transform.class"]][["all.comparisons"]][[train.set]])) {
    for(gene.set in names(res.aucs[["no.transform.class"]][["all.comparisons"]][[train.set]][[test.set]])) {

      res <- res.aucs[["no.transform.class"]][["all.comparisons"]][[train.set]][[test.set]][[gene.set]]
      for(i in 1:length(res)) {
        for(j in 1:length(res[[i]])) {
          sub <- res[[i]][[j]]
          train.drug <- sub$train.drug
          model <- sub$model
          alpha <- sub$alpha
          s <- sub$s
          tmp <- extract.predicted.actual(res, train.drug = train.drug, alpha = alpha, model = model, s = s)

          if((model == "glmnet") && !is.na(alpha) && (alpha == 0)) {
            model <- "ridge"
          }
          if((model == "glmnet") && !is.na(alpha) && (alpha == 1)) {
            model <- "lasso"
          }

          str <- gsub(pattern=",", x=train.drug, replacement="-")
          str <- gsub(pattern=",", x=str, replacement="-")

          predicted <- tmp$predicted
          actual <- tmp$actual
          drug.names <- sub$drug.names
          drug.indicators <- sub$drug.indicators

          file <- paste("cluster-model", train.set, "vs", test.set, str, model, gene.set, "actual-vs-predicted.pdf", sep="-")
          pdf(file)
          g1 <- plot.r2(predicted, actual, colors = drug.names[drug.indicators])
          g1 <- g1 + ggtitle(paste0(train.drug, " (", paste(model, gene.set, paste0(train.set, " vs ", test.set), sep=", "), ")"))
          g1 <- g1 + ylab(paste0(test.set, " Actual Response"))
          g1 <- g1 + xlab(paste0(test.set, " Predicted Response"))
          print(g1)
          d <- dev.off()

          names(tmp$predicted) <- tmp$sample.names
          file <- paste(file.prefix, train.set, "vs", test.set, str, model, gene.set, "prediction.tsv", sep="-")
          write.table(file = file, t(tmp$predicted), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
        }

      }
    }
  }
}

cat("Done\n")

