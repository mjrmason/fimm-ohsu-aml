initialize <- TRUE

load(".Rdata.overlap.120817")

if(initialize) {
suppressPackageStartupMessages(library(glmnet))
suppressPackageStartupMessages(library("randomForest"))
suppressPackageStartupMessages(library("e1071"))
suppressPackageStartupMessages(library("foreach"))
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(tidyr))

suppressPackageStartupMessages(library("minpack.lm"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("gridExtra"))
suppressPackageStartupMessages(library("caTools"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("pcaMethods"))
suppressPackageStartupMessages(library("VennDiagram"))
suppressPackageStartupMessages(library("caret"))

suppressPackageStartupMessages(library(ggbeeswarm))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(synapseClient))

suppressPackageStartupMessages(library("GSVA"))
# use GSA to read in a gmt file
suppressPackageStartupMessages(library("GSA"))

suppressPackageStartupMessages(library(mygene))
suppressPackageStartupMessages(library(gplots))

## Bring in common code
source("../common/data-preprocessing.R")
source("../common/models.R")
source("../common/plotting.R")
source("../common/utils.R")

## Register parallelization
num.cores <- detectCores()
if(!is.na(num.cores) && (num.cores > 1)) {
  suppressPackageStartupMessages(library("doMC"))
  cat(paste("Registering ", num.cores-1, " cores.\n", sep=""))
  registerDoMC(cores=(num.cores-1))
}
num.processes <- num.cores - 1

## Log in to Synapse
synapseLogin()

## Begin setup (this will need to be changed as the data sets in the intersection change)

## Purpose:
## This file models 2 data sets (OHSU and FIMM)
source("fimm-ohsu-setup-120817.R")
} # end initialize


## ensg.to.pathway.list <- list()

cat("Beginning\n")
print(names(ensg.to.pathway.list))
gsets <- unique(trn.tbl$gene.set)
## gsets <- gsets[gsets != "gene"]
gsets <- c("biocarta")
gsets <- c("hallmark")
gsets <- c("reactome")
gsets <- c("gene")
gsets <- c("kegg")
gsets <- unique(as.character(trn.tbl$gene.set))
gsets <- c("gene", gsets[gsets != "gene"])
## gsets <- c("gene")
max.deviance.df <- c()
for(gset in gsets) {
  cat(paste0("Doing gset: ", gset, "\n"))
  tbl.gset <- subset(trn.tbl, gene.set == gset)
  models <-  unique(as.character(tbl.gset$model))
##  models <- c("rf")
  models <- c("ridge", "rf", models[!(models %in% c("ridge", "rf"))])
##  models <- c("ridge", models[!(models %in% c("ridge", "rf"))])
##  models <- c("ridge")
##  models <- c("rf")
##  models <- c("rf", "crf")
  models <- models[models %in% tbl.gset$model]
  for(mdl in models) {
    cat(paste0("Doing model: ", mdl, "\n"))
    tbl.model <- subset(tbl.gset, model == mdl)
    responses <- unique(as.character(tbl.model$pt.response))
##    responses <- "no.transform"
    responses <- c("random", "mean.feature", "no.transform", "model.mean.feature", responses[!(responses %in% c("mean.feature", "random", "no.transform", "model.mean.feature"))])
##    responses <- c("random", "mean.feature")
    responses <- responses[responses %in% tbl.model$pt.response]

    expected.fracs <- list()
    if(gset == "gene") {
      pthwys <- names(ensg.to.pathway.list)
      if(length(pthwys) > 0) {
      for(gs.nm in pthwys[pthwys %in% c("go.BP", "kegg")]) {
        cat(paste0("Doing gs.nm: ", gs.nm, "\n"))
        universe <- universes[[gset]]
print(names(ensg.to.pathway.list))
print(head(universe))
print(length(universe[!grepl(pattern="mean", universe)]))
print(head(names(ensg.to.pathway.list[[gs.nm]])))
        mx <- max(length(universe[!grepl(pattern="mean", universe)]), length(names(ensg.to.pathway.list[[gs.nm]])))
print(mx)
        list.sizes <- c(seq(from=5,to=100, by=5), seq(from=100, to=500, by=10), seq(from=500, to=1000, by=25))
        list.sizes <- list.sizes[list.sizes <= mx]
        list.sizes <- list.sizes[list.sizes > 0]
        list.sizes <- unique(sort(list.sizes))
        names(list.sizes) <- list.sizes
print(head(list.sizes))

cat(paste0("Calculating expected overlap for ", gs.nm, " and ", mdl, "\n"))
      expected.fracs[[gs.nm]] <- llply(list.sizes, .parallel = TRUE,
                           .fun = function(sz) {
        switch(mdl,
          "rf" = { calc.expected.overlap.of.pathways(universe = universe, transform = "none", top = sz, pathway.set = ensg.to.pathway.list[[gs.nm]]) },
          "crf" = { calc.expected.overlap.of.pathways(universe = universe, transform = "none", top = sz, pathway.set = ensg.to.pathway.list[[gs.nm]]) },
          "lasso" = { calc.expected.overlap.of.pathways(universe = universe, transform = "abs.drop.zero", top = sz, pathway.set = ensg.to.pathway.list[[gs.nm]]) },
          "ridge" = { calc.expected.overlap.of.pathways(universe = universe, transform = "abs.drop.zero", top = sz, pathway.set = ensg.to.pathway.list[[gs.nm]]) },
          { stop(paste0("Wasn't expecting ", mdl, "\n")) }) })
      }
cat(paste0("Done calculating expected overlap for ", gs.nm, " and ", mdl, "\n"))
print(head(names(expected.fracs[[gs.nm]])))
print(head(expected.fracs[[gs.nm]][[list.sizes[1]]]))
      }
    }

    for(resp in responses) {
      tbl.resp <- subset(tbl.model, pt.response == resp)
      universe <- universes[[gset]]
      if(grepl(resp, pattern= "mean.feature")) { universe <- c(universe, "mean.resp") }

      print(c(gset, mdl, resp))

      ## Look at overlap between OHSU and FIMM
      tbl.sub <- subset(tbl.full, model == mdl & gene.set == gset & pt.response == resp & train.set == "ohsu" & test.set == "fimm" & metric == "pearson")
      colnames(tbl.sub)[colnames(tbl.sub) == "pval"] <- "corr.pval"
      colnames(tbl.sub)[colnames(tbl.sub) == "val"] <- "corr"

      cutoff <- 0.05
if(FALSE) {
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
}

      ppath <- "top-feature-plots"
      dir.create(ppath)
      if(gset == "gene") {
        pthwys <- names(ensg.to.pathway.list)
        if(length(pthwys) > 0) { 
        for(gs.nm in pthwys[pthwys %in% c("go.BP", "kegg")]) {

          path <- paste0(ppath, "/", paste(mdl, "path", gs.nm, resp, sep="-"))
          dir.create(path)
          ofile <- paste0(path, "/", file.prefix, "-", paste(mdl, "path", gs.nm, resp, sep="-"), "-ohsu-vs-fimm-top-feature-overlap")
          main <- paste(mdl, "path", gs.nm, resp, sep=" ")

## We will not consider "mean.resp" in the feature overlap.
      mx <- max(length(universe[!grepl(pattern="mean", universe)]), length(names(ensg.to.pathway.list[[gs.nm]])))
##      list.sizes <- c(seq(from=5,to=100, by=5), seq(from=100, to=500, by=10), seq(from=500, to=1000, by=25), 2000,5000,7000,7500,7600,mx-100, mx-10,mx-1, mx,mx+1,8000) 
      list.sizes <- c(seq(from=5,to=100, by=5), seq(from=100, to=500, by=10), seq(from=500, to=1000, by=25))
      list.sizes <- list.sizes[list.sizes <= mx]
      list.sizes <- list.sizes[list.sizes > 0]
      list.sizes <- unique(sort(list.sizes))
      names(list.sizes) <- list.sizes

          save(tbl.resp, file="tbl.resp.Rd")
          save(tbl.sub, file="tbl.sub.Rd")
          save(universe, file="universe.Rd")
          save(list.sizes, file="list.sizes.Rd") 
          save(ensg.to.pathway.list, file="ensg.to.pathway.list.Rd")
          save(gs.nm, file="gs.nm.Rd")
          save(ofile, file="ofile.Rd")
          save(main, file="main.Rd")

cat(paste0("Calling plot.overlap.of.pathways: ", gs.nm, "\n"))

          switch(mdl,
            "rf" = { 
##                      plot.overlap.of.pathways(tbl.resp, tbl.sub, universe = universe, transform = "none", list.sizes = list.sizes, ofile.prefix = ofile, main = main, pathway.set = ensg.to.pathway.list[[gs.nm]]) 
                      plot.overlap.of.pathways(tbl.resp, tbl.sub, universe = universe, transform = "none", list.sizes = list.sizes, ofile.prefix = ofile, main = main, pathway.set = ensg.to.pathway.list[[gs.nm]], expected.fracs = expected.fracs[[gs.nm]]) 
            },
            "crf" = { 
##                      plot.overlap.of.pathways(tbl.resp, tbl.sub, universe = universe, transform = "none", list.sizes = list.sizes, ofile.prefix = ofile, main = main, pathway.set = ensg.to.pathway.list[[gs.nm]]) 
                      plot.overlap.of.pathways(tbl.resp, tbl.sub, universe = universe, transform = "none", list.sizes = list.sizes, ofile.prefix = ofile, main = main, pathway.set = ensg.to.pathway.list[[gs.nm]], expected.fracs = expected.fracs[[gs.nm]]) 
            },
            "lasso" = { 
##                        plot.overlap.of.pathways(tbl.resp, tbl.sub, universe = universe, transform = "abs.drop.zero", list.sizes = list.sizes, ofile.prefix = ofile, main = main, pathway.set = ensg.to.pathway.list[[gs.nm]]) 
                        plot.overlap.of.pathways(tbl.resp, tbl.sub, universe = universe, transform = "abs.drop.zero", list.sizes = list.sizes, ofile.prefix = ofile, main = main, pathway.set = ensg.to.pathway.list[[gs.nm]], expected.fracs = expected.fracs[[gs.nm]]) 
            },
            "ridge" = { 
##                        plot.overlap.of.pathways(tbl.resp, tbl.sub, universe = universe, transform = "abs.drop.zero", list.sizes = list.sizes, ofile.prefix = ofile, main = main, pathway.set = ensg.to.pathway.list[[gs.nm]]) 
                        plot.overlap.of.pathways(tbl.resp, tbl.sub, universe = universe, transform = "abs.drop.zero", list.sizes = list.sizes, ofile.prefix = ofile, main = main, pathway.set = ensg.to.pathway.list[[gs.nm]], expected.fracs = expected.fracs[[gs.nm]]) 
            },
            { stop(paste0("Wasn't expecting ", mdl, "\n")) })


##          plot.overlap.of.pathways(tbl.resp, tbl.sub, universe = universe, transform = "abs.drop.zero", list.sizes = list.sizes, ofile.prefix = ofile, main = main, pathway.set = ensg.to.pathway.list[[gs.nm]])
cat("Done calling plot.overlap.of.pathways\n")
        }
        }
      }


## We will not consider "mean.resp" in the feature overlap.
      mx <- length(universe[!grepl(pattern="mean", universe)])
      list.sizes <- c(seq(from=5,to=100, by=5), seq(from=100, to=500, by=10), seq(from=500, to=1000, by=25))
      list.sizes <- list.sizes[list.sizes <= mx]
      list.sizes <- list.sizes[list.sizes > 0]
      list.sizes <- unique(sort(list.sizes))
      names(list.sizes) <- list.sizes

##          save(tbl.resp, file="tbl.resp.Rd")
##          save(universe, file="universe.Rd")
##          save(list.sizes, file="list.sizes.Rd") 

      cat(paste0("About to do feature-level: ", gset, "\n"))
      num.sigs <- ldply(list.sizes, .parallel = TRUE,
                           .fun = function(sz) {
                                    tmp <- NULL
                                    switch(mdl,

                                      "rf" = { tmp <- find.overlap.of.sparse.predictors(tbl.resp, universe = universe, transform = "none", top = sz, verbose = FALSE, account.for.sign.of.coeff = FALSE) },
                                      "crf" = { tmp <- find.overlap.of.sparse.predictors(tbl.resp, universe = universe, transform = "none", top = sz, verbose = FALSE, account.for.sign.of.coeff = FALSE) },
                                      "lasso" = { tmp <- find.overlap.of.sparse.predictors(tbl.resp, universe = universe, transform = "abs.drop.zero", top = sz, verbose = FALSE, account.for.sign.of.coeff = TRUE) },
                                      "ridge" = { tmp <- find.overlap.of.sparse.predictors(tbl.resp, universe = universe, transform = "abs.drop.zero", top = sz, verbose = FALSE, account.for.sign.of.coeff = TRUE) },
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
##      flag <- !is.na(tmp$corr.pval) & !is.na(tmp$corr) & (tmp$corr.pval < 0.05) & (tmp$corr > 0)
##      if(any(flag)) { }
##      tmp <- tmp[flag,,drop=F]
      tmp$FIMM_DRUG_NAME <- factor(as.character(tmp$FIMM_DRUG_NAME), levels = unique(tmp$FIMM_DRUG_NAME))
      ppath <- "top-feature-plots"
      dir.create(ppath)
      path <- paste0(ppath, "/", paste(mdl, gset, resp, sep="-"))
      dir.create(path)
      gs <- dlply(tmp, .variables = c("FIMM_DRUG_NAME"), .parallel = FALSE,
            .fun = function(df) {
                     pdf(paste0(path, "/", file.prefix, "-", paste(mdl, gset, resp, df$FIMM_DRUG_NAME[1], sep="-"), "-ohsu-vs-fimm-top-feature-overlap.pdf"))
write.table(file = paste0(path, "/", file.prefix, "-", paste(mdl, gset, resp, df$FIMM_DRUG_NAME[1], sep="-"), "-ohsu-vs-fimm-top-feature-overlap.tsv"), df, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
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
                     list("g1" = g1, "g2" = g2, "max.deviance" = max.deviance, "min.pval" = min.pval, "corr" = df$corr[1], "corr.pval" = df$corr.pval[1])
            })
      md <- unlist(lapply(gs, function(l) l$max.deviance))
      mp <- unlist(lapply(gs, function(l) l$min.pval))
      md.df <- data.frame(max.deviance = md, min.pval = mp, gene.set = gset, model = mdl, response = resp, cmp = "ohsu-vs-fimm", stringsAsFactors = FALSE)
      max.deviance.df <- rbind(max.deviance.df, md.df)
      pdf(paste0(path, "/", file.prefix, "-", paste(mdl, gset, resp, "all", sep="-"), "-ohsu-vs-fimm-top-feature-overlap.pdf"))
      any.sig <- FALSE
      for(i in 1:length(gs)) {
        grid.arrange(gs[[i]]$g1, gs[[i]]$g2)
        if( !is.na(gs[[i]]$corr.pval) && !is.na(gs[[i]]$corr) && (as.numeric(gs[[i]]$corr.pval) < 0.05) && (as.numeric(gs[[i]]$corr) > 0) ) { any.sig <- TRUE }
      }
      d <- dev.off()
      if(any.sig) {
        pdf(paste0(path, "/", file.prefix, "-", paste(mdl, gset, resp, "all-sig", sep="-"), "-ohsu-vs-fimm-top-feature-overlap.pdf"))
        for(i in 1:length(gs)) {
          if( !is.na(gs[[i]]$corr.pval) && !is.na(gs[[i]]$corr) && (as.numeric(gs[[i]]$corr.pval) < 0.05) && (as.numeric(gs[[i]]$corr) > 0) ) { 
            grid.arrange(gs[[i]]$g1, gs[[i]]$g2)
          }
        }
        d <- dev.off()
      }

##      } ##       if((any(flag)))
if(FALSE) {
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
      flag <- tmp$corr.pval < 0.05 & tmp$corr > 0

      tmp <- tmp[flag,,drop=F]
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
    }
if(FALSE) {
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
}
if(FALSE) {
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

cat("Done with foo.R\n")