gsets <- unique(trn.tbl$gene.set)
## gsets <- gsets[gsets != "gene"]
gsets <- c("biocarta")
gsets <- c("hallmark")
gsets <- c("reactome")
gsets <- c("gene")
gsets <- c("kegg")
gsets <- unique(as.character(trn.tbl$gene.set))
gsets <- c("gene", gsets[gsets != "gene"])
max.deviance.df <- c()
for(gset in gsets) {
  tbl.gset <- subset(trn.tbl, gene.set == gset)
  models <-  unique(as.character(tbl.gset$model))
##  models <- c("rf")
  models <- c("rf", models[models != "rf"])
  for(mdl in models) {
    tbl.model <- subset(tbl.gset, model == mdl)
    responses <- unique(as.character(tbl.model$pt.response))
##    responses <- "no.transform"
    responses <- c("random", "mean.feature", "no.transform", responses[!(responses %in% c("mean.feature", "random", "no.transform"))])
    for(resp in responses) {
      tbl.resp <- subset(tbl.model, pt.response == resp)
      universe <- universes[[gset]]
      if(resp == "mean.feature") { universe <- c(universe, "mean.resp") }

## We will not consider "mean.resp" in the feature overlap.
      mx <- length(universe[!grepl(pattern="mean", universe)])
      by <- 10
      from <- 10
      if(gset != "gene") { by <- 5; from <- 5 }
      list.sizes <- c(seq(from=from,to=1000, by=by), 10,100,1000,2000,5000,7000,7500,7600,mx-100,
                      mx-10,mx-1, mx,mx+1,8000) 
      list.sizes <- list.sizes[list.sizes <= mx]
      list.sizes <- list.sizes[list.sizes > 0]
      list.sizes <- sort(list.sizes)
      names(list.sizes) <- list.sizes



      print(c(gset, mdl, resp))

      ## Look at overlap between OHSU and FIMM
      tbl.sub <- subset(tbl.full, model == mdl & gene.set == gset & pt.response == resp & train.set == "ohsu" & test.set == "fimm" & metric == "pearson")
      colnames(tbl.sub)[colnames(tbl.sub) == "pval"] <- "corr.pval"
      colnames(tbl.sub)[colnames(tbl.sub) == "val"] <- "corr"

      cutoff <- 0.05
      if(any(tbl.sub$corr.pval < cutoff)) {
        drugs <- tbl.sub$train.drug[!is.na(tbl.sub$corr.pval) & (tbl.sub$corr.pval < cutoff) & !is.na(tbl.sub$corr) & (tbl.sub$corr > 0)]
        names(drugs) <- drugs

        overlapping.markers <- NULL
        switch(mdl,
          "lasso" = { overlapping.markers <- llply(drugs, .fun = function(drug) common.markers(res.aucs[[resp]], gset, "fimm", drug, "ohsu", drug, alpha = 1, model = "glmnet", top = 25)) },
          "ridge" = { overlapping.markers <- llply(drugs, .fun = function(drug) common.markers(res.aucs[[resp]], gset, "fimm", drug, "ohsu", drug, alpha = 0, model = "glmnet", top = 25)) },
          "rf" = { overlapping.markers <- llply(drugs, .fun = function(drug) common.markers(res.aucs[[resp]], gset, "fimm", drug, "ohsu", drug, alpha = NA, model = "rf", top = 25)) },
          { stop(paste0("Unexpected model: ", mdl, "\n")) })

        tmp <- ldply(overlapping.markers, .fun = function(lst) data.frame(gene = lst))

        baz <- NULL
        if(nrow(tmp) > 0) {
          colnames(tmp) <- c("drug", "gene")
          if(gset == "gene") {
            baz <- queryMany(tmp$gene, scopes="ensembl.gene", fields=c("symbol"), species="human")
            if(nrow(baz) > 0) {
              baz <- as.data.frame(baz[, c("query", "symbol")])
##              baz <- rbind(baz, c("mean.resp", "mean.resp"))
            }
          }
        }
        mat <- NULL
        if((nrow(tmp) > 0) && !is.null(baz) && (nrow(baz) > 0)) {
          overlapping.markers <- merge(tmp, unique(as.data.frame(baz[,c("query","symbol")])), by.x = "gene", by.y = "query")
          overlapping.markers <- overlapping.markers[!is.na(overlapping.markers$drug) & !is.na(overlapping.markers$gene) & !is.na(overlapping.markers$symbol), ]
          mat <- spread(overlapping.markers[, c("drug", "gene", "symbol")], key = drug, value = gene)
        }
        if((nrow(tmp) > 0) && (gset != "gene")) {
          mat <- spread(tmp, key = drug, value = drug)
        }
        if(!is.null(mat)) {
          rownames(mat) <- mat[, 1]
          mat <- as.matrix(mat[, -1])
          mat[!is.na(mat)] <- 1
          mat[is.na(mat)] <- 0
          class(mat) <- "numeric"
          flag <- rowSums(mat) > 1
          if(length(which(flag)) < 2) {
            flag <- rowSums(mat) > 0
          } 
          if((nrow(mat[flag,,drop=F]) >= 2) && (ncol(mat[flag,,drop=F]) >= 2)) {
            pdf(paste0(file.prefix, "-fimm-vs-ohsu-", paste(mdl, gset, resp, sep="-"), "-pearson-sig-drugs-feature-heatmap.pdf"))
            ## heatmap(mat[flag,,drop=F], Rowv = NA, Colv = NA, scale = "none", main = paste0(paste(mdl, gset, resp, collapse=" "), ":\nsignificant (p < 0.05), positively correlated drugs"))
            heatmap(mat[flag,,drop=F], Rowv = NA, Colv = NA, scale = "none", main = paste0(paste(mdl, gset, resp, collapse=" "), ": pos cor; p < ", cutoff))
            d <- dev.off()
          }
        }
      }

      num.sigs <- ldply(list.sizes, .parallel = TRUE,
                           .fun = function(sz) {
                                    tmp <- NULL
                                    switch(mdl,
                                      "rf" = { tmp <- find.overlap.of.sparse.predictors(tbl.resp, universe = universe, transform = "none", top = sz, verbose = FALSE) },
                                      "lasso" = { tmp <- find.overlap.of.sparse.predictors(tbl.resp, universe = universe, transform = "abs.drop.zero", top = sz, verbose = FALSE) },
                                      "ridge" = { tmp <- find.overlap.of.sparse.predictors(tbl.resp, universe = universe, transform = "abs.drop.zero", top = sz, verbose = FALSE) },
                                      { stop(paste0("Wasn't expecting ", mdl, "\n")) })
                                    tmp <- tmp[, c("FIMM_DRUG_NAME", "train.set.1", "train.set.2", "n.1", "n.2", "n.both", "n.neither", "n.expected", "pval", "top", "gene.set.1")]
                                    tmp
                                  })

      tmp <- merge(subset(num.sigs, (train.set.1 == "ohsu") & (train.set.2 == "fimm")), tbl.sub, by.x = c("FIMM_DRUG_NAME", "gene.set.1"), by.y = c("test.drug", "gene.set"))
      tmp$top <- as.numeric(tmp$top)
      tmp$n.both <- as.numeric(tmp$n.both)
      tmp$n.expected <- as.numeric(tmp$n.expected)
      tmp$pval <- as.numeric(tmp$pval)
      tmp$corr.pval <- as.numeric(tmp$corr.pval)
      tmp <- tmp[order(tmp$corr.pval), ]
      tmp$FIMM_DRUG_NAME <- factor(as.character(tmp$FIMM_DRUG_NAME), levels = unique(tmp$FIMM_DRUG_NAME))
      path <- paste(mdl, gset, resp, sep="-")
      dir.create(path)
      gs <- dlply(tmp, .variables = c("FIMM_DRUG_NAME"), .parallel = FALSE,
            .fun = function(df) {
                     pdf(paste0(path, "/", file.prefix, "-", paste(mdl, gset, resp, df$FIMM_DRUG_NAME[1], sep="-"), "-ohsu-vs-fimm-top-feature-overlap.pdf"))
                     flag <- !is.na(df$n.both)
                     corr <- df$corr[1]; corr.pval <- df$corr.pval[1]
    ##                 if((corr < 0) || (corr.pval > 0.05)) { return() }
    ##                 if(!any(df$pval < 0.05)) { return() }
    ##                 plot(df$top[flag], df$n.both[flag], main = paste0(df$FIMM_DRUG_NAME[1], ": corr = ", signif(df$corr[1], digits=2), " pval = ", 
    ##                      signif(df$corr.pval[1], digits=2)))
##                     g1 <- ggplot(data = df[flag,]) + geom_point(aes(x = top, y = n.both)) + xlab("Num Top Features") + ylab("Num Overlapping")
                     df$frac <- df$n.both / df$top
                     df$frac.expected <- df$n.expected / df$top
                     g1 <- ggplot(data = df[flag,]) + geom_point(aes(x = top, y = frac)) + xlab("Num Top Features") + ylab("Fraction Overlapping")
                     g1 <- g1 + ggtitle(paste0(paste(df$FIMM_DRUG_NAME[1], mdl, gset, resp, sep=" "), ": corr = ", signif(df$corr[1], digits=2), " pval = ", signif(df$corr.pval[1], digits=2)))
##                     g1 <- g1 + geom_abline(intercept = 0, slope = 1)
                     ## g1 <- g1 + geom_hline(yintercept = frac.expected)
                     g1 <- g1 + geom_line(aes(x = top, y = frac.expected), linetype = "dashed")
              
                     g2 <- ggplot(data = df[flag,]) + geom_point(aes(x = top, y = -log10(pval))) + xlab("Num Top Features") + ylab("Overlap -log10 p-value")
                     g2 <- g2 + geom_hline(yintercept = -log10(0.05))
                     grid.arrange(g1, g2)
                     d <- dev.off()
                     max.deviance <- max(df$frac - df$frac.expected, na.rm=TRUE)
                     min.pval <- NA
                     if(any(!is.na(df$pval))) { min.pval <- min(df$pval, na.rm=TRUE) }
                     list("g1" = g1, "g2" = g2, "max.deviance" = max.deviance, "min.pval" = min.pval)
            })
      md <- unlist(lapply(gs, function(l) l$max.deviance))
      mp <- unlist(lapply(gs, function(l) l$min.pval))
      md.df <- data.frame(max.deviance = md, min.pval = mp, gene.set = gset, model = mdl, response = resp, cmp = "ohsu-vs-fimm", stringsAsFactors = FALSE)
      max.deviance.df <- rbind(max.deviance.df, md.df)
      pdf(paste0(file.prefix, "-", paste(mdl, gset, resp, "all", sep="-"), "-ohsu-vs-fimm-top-feature-overlap.pdf"))
      for(i in 1:length(gs)) {
        grid.arrange(gs[[i]]$g1, gs[[i]]$g2)
      }
      d <- dev.off()

      ## Look at overlap between OHSU1 and OHSU2
      tbl.sub <- subset(tbl.full, model == mdl & gene.set == gset & pt.response == resp & train.set == "ohsu.set1" & test.set == "ohsu.set2" & metric == "pearson")
      colnames(tbl.sub)[colnames(tbl.sub) == "pval"] <- "corr.pval"
      colnames(tbl.sub)[colnames(tbl.sub) == "val"] <- "corr"

      if(any(tbl.sub$corr.pval < cutoff)) {
        drugs <- tbl.sub$train.drug[!is.na(tbl.sub$corr.pval) & (tbl.sub$corr.pval < cutoff) & !is.na(tbl.sub$corr) & (tbl.sub$corr > 0)]
        names(drugs) <- drugs

        overlapping.markers <- NULL
        switch(mdl,
          "lasso" = { overlapping.markers <- llply(drugs, .fun = function(drug) common.markers(res.aucs[[resp]], gset, "fimm", drug, "ohsu", drug, alpha = 1, model = "glmnet", top = 25)) },
          "ridge" = { overlapping.markers <- llply(drugs, .fun = function(drug) common.markers(res.aucs[[resp]], gset, "fimm", drug, "ohsu", drug, alpha = 0, model = "glmnet", top = 25)) },
          "rf" = { overlapping.markers <- llply(drugs, .fun = function(drug) common.markers(res.aucs[[resp]], gset, "fimm", drug, "ohsu", drug, alpha = NA, model = "rf", top = 25)) },
          { stop(paste0("Unexpected model: ", mdl, "\n")) })

        tmp <- ldply(overlapping.markers, .fun = function(lst) data.frame(gene = lst))

        baz <- NULL
        if(nrow(tmp) > 0) {
          colnames(tmp) <- c("drug", "gene")
          if(gset == "gene") {
            baz <- queryMany(tmp$gene, scopes="ensembl.gene", fields=c("symbol"), species="human")
            if(nrow(baz) > 0) {
              baz <- as.data.frame(baz[, c("query", "symbol")])
##              baz <- rbind(baz, c("mean.resp", "mean.resp"))
            }
          }
        }
        mat <- NULL
        if((nrow(tmp) > 0) && !is.null(baz) && (nrow(baz) > 0)) {
          overlapping.markers <- merge(tmp, unique(as.data.frame(baz[,c("query","symbol")])), by.x = "gene", by.y = "query")
          overlapping.markers <- overlapping.markers[!is.na(overlapping.markers$drug) & !is.na(overlapping.markers$gene) & !is.na(overlapping.markers$symbol), ]
          mat <- spread(overlapping.markers[, c("drug", "gene", "symbol")], key = drug, value = gene)
        }
        if((nrow(tmp) > 0) && (gset != "gene")) {
          mat <- spread(tmp, key = drug, value = drug)
        }
        if(!is.null(mat)) {
          rownames(mat) <- mat[, 1]
          mat <- as.matrix(mat[, -1])
          mat[!is.na(mat)] <- 1
          mat[is.na(mat)] <- 0
          class(mat) <- "numeric"
          flag <- rowSums(mat) > 1
          if(length(which(flag)) < 2) {
            flag <- rowSums(mat) > 0
          } 
          if((nrow(mat[flag,,drop=F]) >= 2) && (ncol(mat[flag,,drop=F]) >= 2)) {
            pdf(paste0(file.prefix, "-ohsu-vs-ohsu-", paste(mdl, gset, resp, sep="-"), "-pearson-sig-drugs-feature-heatmap.pdf"))
            ## heatmap(mat[flag,,drop=F], Rowv = NA, Colv = NA, scale = "none", main = paste0(paste(mdl, gset, resp, collapse=" "), ":\nsignificant (p < 0.05), positively correlated drugs"))
            heatmap(mat[flag,,drop=F], Rowv = NA, Colv = NA, scale = "none", main = paste0(paste(mdl, gset, resp, collapse=" "), ": pos cor; p < ", cutoff))
            d <- dev.off()
          }
        }
      }

      tmp <- merge(subset(num.sigs, (train.set.1 == "ohsu.set1") & (train.set.2 == "ohsu.set2")), tbl.sub, by.x = c("FIMM_DRUG_NAME", "gene.set.1"), by.y = c("test.drug", "gene.set"))
      tmp$top <- as.numeric(tmp$top)
      tmp$n.both <- as.numeric(tmp$n.both)
      tmp$n.expected <- as.numeric(tmp$n.expected)
      tmp$pval <- as.numeric(tmp$pval)
      tmp$corr.pval <- as.numeric(tmp$corr.pval)
      tmp <- tmp[order(tmp$corr.pval), ]
      tmp$FIMM_DRUG_NAME <- factor(as.character(tmp$FIMM_DRUG_NAME), levels = unique(tmp$FIMM_DRUG_NAME))
      gs <- dlply(tmp, .variables = c("FIMM_DRUG_NAME"), .parallel = FALSE,
            .fun = function(df) {
                     pdf(paste0(path, "/", file.prefix, "-", paste(mdl, gset, resp, df$FIMM_DRUG_NAME[1], sep="-"), "-ohsu-vs-ohsu-top-feature-overlap.pdf"))
                     flag <- !is.na(df$n.both)
                     corr <- df$corr[1]; corr.pval <- df$corr.pval[1]
    ##                 if((corr < 0) || (corr.pval > 0.05)) { return() }
    ##                 if(!any(df$pval < 0.05)) { return() }
    ##                 plot(df$top[flag], df$n.both[flag], main = paste0(df$FIMM_DRUG_NAME[1], ": corr = ", signif(df$corr[1], digits=2), " pval = ", 
    ##                      signif(df$corr.pval[1], digits=2)))
##                     g1 <- ggplot(data = df[flag,]) + geom_point(aes(x = top, y = n.both)) + xlab("Num Top Features") + ylab("Num Overlapping")
                     df$frac <- df$n.both / df$top
                     df$frac.expected <- df$n.expected / df$top
                     g1 <- ggplot(data = df[flag,]) + geom_point(aes(x = top, y = frac)) + xlab("Num Top Features") + ylab("Fraction Overlapping")
                     g1 <- g1 + ggtitle(paste0(paste(df$FIMM_DRUG_NAME[1], mdl, gset, resp, sep=" "), ": corr = ", signif(df$corr[1], digits=2), " pval = ", signif(df$corr.pval[1], digits=2)))
##                     g1 <- g1 + geom_abline(intercept = 0, slope = 1)
                     ## g1 <- g1 + geom_hline(yintercept = frac.expected)
                     g1 <- g1 + geom_line(aes(x = top, y = frac.expected), linetype = "dashed")
                
                     g2 <- ggplot(data = df[flag,]) + geom_point(aes(x = top, y = -log10(pval))) + xlab("Num Top Features") + ylab("Overlap -log10 p-value")
                     g2 <- g2 + geom_hline(yintercept = -log10(0.05))
                     grid.arrange(g1, g2)
                     list("g1" = g1, "g2" = g2)
                     d <- dev.off()
                     max.deviance <- max(df$frac - df$frac.expected, na.rm=TRUE)
                     min.pval <- NA
                     if(any(!is.na(df$pval))) { min.pval <- min(df$pval, na.rm=TRUE) }
                     list("g1" = g1, "g2" = g2, "max.deviance" = max.deviance, "min.pval" = min.pval)
            })
      md <- unlist(lapply(gs, function(l) l$max.deviance))
      mp <- unlist(lapply(gs, function(l) l$min.pval))
      md.df <- data.frame(max.deviance = md, min.pval = mp, gene.set = gset, model = mdl, response = resp, cmp = "ohsu-vs-ohsu", stringsAsFactors = FALSE)
      max.deviance.df <- rbind(max.deviance.df, md.df)

      pdf(paste0(file.prefix, "-", paste(mdl, gset, resp, "all", sep="-"), "-ohsu-vs-ohsu-top-feature-overlap.pdf"))
      for(i in 1:length(gs)) {
        grid.arrange(gs[[i]]$g1, gs[[i]]$g2)
      }
      d <- dev.off()

    }

    cat("Doing max deviance\n")
    save(max.deviance.df, file = "md.df.Rd")
    max.deviance.df$max.deviance <- as.numeric(max.deviance.df$max.deviance)
    max.deviance.df$min.pval <- as.numeric(max.deviance.df$min.pval)
    max.deviance.df$response <- as.character(max.deviance.df$response)
    max.deviance.df$gene.set <- as.character(max.deviance.df$gene.set)
    max.deviance.df$model <- as.character(max.deviance.df$model)
    max.deviance.df$cmp <- as.character(max.deviance.df$cmp)
    print(unique(max.deviance.df[, c("gene.set", "model", "cmp")]))
    df <- subset(max.deviance.df, gene.set == as.character(gset) & model == as.character(mdl) & cmp == "ohsu-vs-fimm")
    g1 <- ggplot(data = df, aes(x = response, y = max.deviance))
    g1 <- g1 + geom_violin()
    g1 <- g1 + geom_boxplot(width = 0.5)
    g1 <- g1 + geom_beeswarm()
    g1 <- g1 + xlab("response")
    g1 <- g1 + ylab("max deviance from expected frac")
    g1 <- g1 + ggtitle(paste(gset, mdl, "ohsu-vs-fimm", sep=" "))
    pdf(paste0(file.prefix, "-", paste(mdl, gset, "ohsu-vs-fimm", "max-deviance.pdf", sep="-")))
    print(g1)
    d <- dev.off()

    df <- subset(max.deviance.df, gene.set == as.character(gset) & model == as.character(mdl) & cmp == "ohsu-vs-ohsu")
    g1 <- ggplot(data = df, aes(x = response, y = max.deviance))
    g1 <- g1 + geom_violin()
    g1 <- g1 + geom_boxplot(width = 0.5)
    g1 <- g1 + geom_beeswarm()
    g1 <- g1 + xlab("response")
    g1 <- g1 + ylab("max deviance from expected frac")
    g1 <- g1 + ggtitle(paste(gset, mdl, "ohsu-vs-ohsu", sep=" "))
    pdf(paste0(file.prefix, "-", paste(mdl, gset, "ohsu-vs-ohsu", "max-deviance.pdf", sep="-")))
    print(g1)
    d <- dev.off()


  }
}

save.image(".Rdata")

gsets <- unique(as.character(max.deviance.df$gene.set))
for(gset in gsets) {
  models <-  unique(as.character(max.deviance.df$model))
  for(mdl in models) {
    df <- subset(max.deviance.df, gene.set == gset & model == mdl & cmp == "ohsu-vs-fimm")
    g1 <- ggplot(data = df, aes(x = response, y = max.deviance))
    g1 <- g1 + geom_violin()
    g1 <- g1 + geom_boxplot(width = 0.5)
    g1 <- g1 + geom_beeswarm()
    g1 <- g1 + xlab("response")
    g1 <- g1 + ylab("max deviance from expected frac")
    g1 <- g1 + ggtitle(paste(gset, mdl, "ohsu-vs-fimm", sep=" "))
    pdf(paste0(file.prefix, "-", paste(mdl, gset, "ohsu-vs-fimm", "max-deviance.pdf", sep="-")))
    print(g1)
    d <- dev.off()

    pval.cutoff <- 0.05
    df.pval <- df[!is.na(df$min.pval) & (df$min.pval < pval.cutoff),]
    g1 <- ggplot(data = df.pval, aes(x = response))
    g1 <- g1 + geom_bar()
    g1 <- g1 + xlab("response")
    g1 <- g1 + ylab(paste0("Num drugs with sig overlap (p < ", pval.cutoff, ")"))
    g1 <- g1 + ggtitle(paste(gset, mdl, "ohsu-vs-fimm", sep=" "))
    pdf(paste0(file.prefix, "-", paste(mdl, gset, "ohsu-vs-fimm", "min-pval.pdf", sep="-")))
    print(g1)
    d <- dev.off()

    df <- subset(max.deviance.df, gene.set == gset & model == mdl & cmp == "ohsu-vs-ohsu")
    g1 <- ggplot(data = df, aes(x = response, y = max.deviance))
    g1 <- g1 + geom_violin()
    g1 <- g1 + geom_boxplot(width = 0.5)
    g1 <- g1 + geom_beeswarm()
    g1 <- g1 + xlab("response")
    g1 <- g1 + ylab("max deviance from expected frac")
    g1 <- g1 + ggtitle(paste(gset, mdl, "ohsu-vs-ohsu", sep=" "))
    pdf(paste0(file.prefix, "-", paste(mdl, gset, "ohsu-vs-ohsu", "max-deviance.pdf", sep="-")))
    print(g1)
    d <- dev.off()

    pval.cutoff <- 0.05
    df.pval <- df[!is.na(df$min.pval) & (df$min.pval < pval.cutoff),]
    g1 <- ggplot(data = df.pval, aes(x = response))
    g1 <- g1 + geom_bar()
    g1 <- g1 + xlab("response")
    g1 <- g1 + ylab(paste0("Num drugs with sig overlap (p < ", pval.cutoff, ")"))
    g1 <- g1 + ggtitle(paste(gset, mdl, "ohsu-vs-ohsu", sep=" "))
    pdf(paste0(file.prefix, "-", paste(mdl, gset, "ohsu-vs-ohsu", "min-pval.pdf", sep="-")))
    print(g1)
    d <- dev.off()

  }
}