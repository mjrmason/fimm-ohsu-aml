
initialize <- TRUE

if(initialize) {
rm(list = ls())

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
suppressPackageStartupMessages(library(mygene))
suppressPackageStartupMessages(library(GGally))

suppressPackageStartupMessages(library("GSVA"))
# use GSA to read in a gmt file
suppressPackageStartupMessages(library("GSA"))
suppressPackageStartupMessages(library("cowplot"))


## Log in to Synapse
synapseLogin()

## Bring in common code
source("../common/data-preprocessing.R")
source("../common/models.R")
source("../common/plotting.R")
source("../common/utils.R")
source("heatmap-utils.R")

## Register parallelization
num.cores <- detectCores()
if(!is.na(num.cores) && (num.cores > 1)) {
  suppressPackageStartupMessages(library("doMC"))
  cat(paste("Registering ", num.cores-1, " cores.\n", sep=""))
  registerDoMC(cores=(num.cores-1))
}
num.processes <- num.cores - 1

## Begin setup (this will need to be changed as the data sets in the intersection change)

source("fimm-ohsu-setup-040218.R")

source("process-drc-and-expr.R")

file.prefix <- "fimm-ohsu-overlap-120817"

rdata.file <- ".Rdata.overlap.120817"
load(rdata.file)

} # if(initialize)

source("plot-pvals.R")

tbl <- tbl.full

renamed.tbl <- tbl
model.renamings <- list("lasso" = "LASSO", "ridge" = "Ridge", "rf" = "Random Forest")
renamed.tbl$model <- as.character(renamed.tbl$model)
for(i in 1:length(model.renamings)) {
  flag <- renamed.tbl$model == names(model.renamings)[i]
  if(any(flag)) { renamed.tbl$model[flag] <- model.renamings[[i]] }
}
response.renamings <- list("mean.feature.only" = "GLDS Only", 
                        "random" = "Shuffled Gene Expression",
                        "model.mean.feature" = "Gene Expression + Penalized Modeled GLDS",
                        "unpenalized.model.mean.feature" = "Gene Expression + Modeled GLDS",
                        "model.mean.feature.only" = "Modeled GLDS Only",
                        "mean.feature" = "Gene Expression + Penalized GLDS",
                        "unpenalized.mean.feature" = "Gene Expression + GLDS",
                        "no.transform" = "Gene Expression")
renamed.tbl$pt.response <- as.character(renamed.tbl$pt.response)
for(i in 1:length(response.renamings)) {
  flag <- renamed.tbl$pt.response == names(response.renamings)[i]
  if(any(flag)) { renamed.tbl$pt.response[flag] <- response.renamings[[i]] }
}

use.ggplot.2.2.1.limit.code <- TRUE

## Make a stripped down version of the above for poster/publication
## Poster publication
tmp.tbl <- renamed.tbl
tmp.tbl <- subset(tmp.tbl, pt.response %in% c("Shuffled Gene Expression", "Gene Expression", "Gene Expression + Modeled GLDS"))
tmp.tbl$pt.response <- factor(tmp.tbl$pt.response, levels = c("Shuffled Gene Expression", "Gene Expression", "Gene Expression + Modeled GLDS"))
tmp.tbl <- subset(tmp.tbl, model %in% c("LASSO", "Ridge"))
tmp.tbl$model <- factor(tmp.tbl$model, levels = c("LASSO", "Ridge"))
## Limit to genes
tmp.tbl <- subset(tmp.tbl, gene.set == "gene")
for(resp in unique(tmp.tbl$response)) {
  for(met in c("pearson", "spearman")) {
##    for(gset in unique(tmp.tbl$gene.set)) { }
##    for(mdl in c("ridge", "lasso", "rf")) { }
      tr.set <- "ohsu"
      tt.set <- "fimm"
      sub <- subset(tmp.tbl, metric == met & response == resp & train.set == tr.set & test.set == tt.set)
      ## Report statistics
      d_ply(sub, .variables = c("pt.response", "gene.set", "model"), 
            .fun = function(df) {
                     m <- df$model[1]
                     p <- df$pt.response[1]
                     g <- df$gene.set[1]
                     med <- median(df$val, na.rm=TRUE)
                     avg <- mean(df$val, na.rm=TRUE)
                     bs <- boxplot.stats(df$val)
                     cat(paste0("stats: metric = ", met, " response = ", resp, " set = ", g, 
                               " model = ", m, " features = ", p, " n = ", length(!is.na(df$val)), 
                               " median = ", med, " mean = ", avg, 
                               " boxplot stats (n = ", bs$n, ") = ", paste(bs$stats, collapse = ", "), "\n"))
            })
      ## Compare groups
      for(m in unique(sub$model)) {
        for(g in unique(sub$gene.set)) {
          df <- subset(sub, (model == m) & (gene.set == g))
          features <- unique(df$pt.response)
          names(features) <- features
          feature.vals <- dlply(df, .variables = c("pt.response"), .fun = function(df2) df2)
          f1 <- "Gene Expression"
          f2 <- "Gene Expression + Modeled GLDS"
          mer <- merge(feature.vals[[f1]], feature.vals[[f2]], by = "train.drug", suffixes = c(".x", ".y"))
          wt <- wilcox.test(mer$val.x, mer$val.y, alternative = "less") 
          cat(paste0("stats: metric = ", met, " response = ", resp, " set = ", g, 
                     " model = ", m, " feature comparison = ", f1, " vs ", f2, 
                     " wilcox.test(", f1, ", ", f2, ", ", wt$alternative, ") p = ", wt$p.value, "\n"))

          if(FALSE) {
          for(i in 1:(length(features)-1)) {
            for(j in (i+1):length(features)) {
              mer <- merge(feature.vals[[features[i]]], feature.vals[[features[j]]], by = "train.drug") 
              if(nrow(mer) != nrow(feature.values[[features[i]]])) {
                stop("Num rows changed after merge\n")
              }
            }
          }
          }
        }
      }
      g <- ggplot(data = sub, aes(x = model, y = val))
      if(length(unique(sub$gene.set)) == 1) {
        g <- g + facet_grid(. ~ pt.response, labeller = labeller(pt.response = label_wrap_gen(15)))
      } else {
        g <- g + facet_grid(gene.set ~ pt.response)
      }
##      g <- g + geom_violin()
      g <- g + geom_boxplot(width = 0.5, outlier.shape = NA)
      g <- g + geom_beeswarm()
      g <- g + geom_hline(yintercept = 0, linetype = "dashed")
      g <- g + ylab(paste0(capwords(met), " Correlation (", toupper(tt.set), " Actual vs Predicted)"))
      g <- g + xlab("Model")
      g <- g + theme(text = element_text(size = 20))
##      g <- g + ggtitle(paste0(resp, " ", capwords(met), " "))
      g <- g + ylim(c(-1,1))


      ## Get the length of the y axis
      ##        ylimits <- ggplot_build(g)$panel$ranges[[1]]$y.range
      ylimits <- NULL
      if(use.ggplot.2.2.1.limit.code) { 
          ylimits <- ggplot_build(g)$layout$panel_ranges[[1]]$y.range
      } else {
          ylimits <- ggplot_build(g)$layout$panel_params[[1]]$y.range
      }
        
      ## When we add error bars, the spacing between them will be in units of yoffset
      yoffset <- 0.01 * ( max(ylimits, na.rm=TRUE) - min(ylimits, na.rm=TRUE) )
      ## We will add 2 error bars
      num.error.bars <- 2
      ## ylimits[2] <- ylimits[2] + 2 * 7 * num.error.bars * yoffset
      ylimits[2] <- ylimits[2] + 2 * 6 * num.error.bars * yoffset
      g <- g + ylim(ylimits)

      gb <- ggplot_build(g)
      ggt <- ggplot_gtable(gb)

      ## ggplot2 doesn't use native units in data space
      ## instead, the data is rescaled to npc, i.e., from 0 to 1
      ## so we need to use the build info to convert from data to [0,1]
      ## ranges <- gb$panel$ranges
      ranges <- NULL
      if(use.ggplot.2.2.1.limit.code) { 
          ranges <- gb$layout$panel_ranges
      } else {
          ranges <- gb$layout$panel_params
      }
      tmp <- ddply(sub, .variables = c("pt.response", "model"), .fun = function(df) max(df$val, na.rm=TRUE))
      panel.maxs <- tmp$V1
      names(panel.maxs) <- unlist(apply(tmp[, 1:2], 1, function(row) paste0(row[1], "-", row[2])))

      facet1 <- "Gene Expression"
      facet2 <- "Gene Expression + Modeled GLDS"
      facetIndx1 <- 2
      facetIndx2 <- 3
      xis <- c("LASSO", "Ridge")

      x1 <- "LASSO"
      x2 <- "LASSO"
      df <- subset(sub, (model == x1) & (gene.set == "gene"))
      features <- unique(df$pt.response)
      names(features) <- features
      feature.vals <- dlply(df, .variables = c("pt.response"), .fun = function(df2) df2)
      mer <- merge(feature.vals[[facet1]], feature.vals[[facet2]], by = "train.drug", suffixes = c(".x", ".y"))
      wt <- wilcox.test(mer$val.x, mer$val.y, alternative = "less") 
      pval <- wt$p.value
      tmp <- draw.err.bar(pval, ranges, panel.maxs, ggt, facet1, x1, facet2, x2, facetIndx1, facetIndx2, xis, yoffset, 
                          plot.pvals.as.stars = TRUE, stat.name = "p") 
      ggt <- tmp[["g"]]
      panel.maxs <- tmp[["panel.maxs"]]

      x1 <- "Ridge"
      x2 <- "Ridge"
      df <- subset(sub, (model == x1) & (gene.set == "gene"))
      features <- unique(df$pt.response)
      names(features) <- features
      feature.vals <- dlply(df, .variables = c("pt.response"), .fun = function(df2) df2)
      mer <- merge(feature.vals[[facet1]], feature.vals[[facet2]], by = "train.drug", suffixes = c(".x", ".y"))
      wt <- wilcox.test(mer$val.x, mer$val.y, alternative = "less") 
      pval <- wt$p.value
      tmp <- draw.err.bar(pval, ranges, panel.maxs, ggt, facet1, x1, facet2, x2, facetIndx1, facetIndx2, xis, yoffset, 
                          plot.pvals.as.stars = TRUE, stat.name = "p") 
      ggt <- tmp[["g"]]
      panel.maxs <- tmp[["panel.maxs"]]

      ## Turn clipping off to see the line across the panels
      ggt$layout$clip <- "off"

      file <- paste0(file.prefix, "-simple-", resp, "-", met, "-", tr.set, "-vs-", tt.set, "-vs-gene-sets.pdf")
      pdf(file)
      grid.draw(ggt)
      d <- dev.off()
##    glist[[length(glist)+1]] <- g
##    { }
  }
##  do.call("grid.arrange", glist)
}

## Poster publication
## Note that this is using "empirical" GLDS (not modeled GLDS)
pt.responses1 <- c("GLDS Only")
pt.responses2 <- c("Gene Expression + GLDS")

## for(mdl in unique(tbl$model)) {
for(mdl in c("LASSO", "Ridge")) {
  for(resp in unique(renamed.tbl$response)) {
    for(gset in c("gene")) {
      for(met in c("pearson", "spearman")) {
        train.ds <- "ohsu"; test.ds <- "fimm"
        for(indx in 1:length(pt.responses1)) {
          pt.response1 <- pt.responses1[indx]; pt.response2 <- pt.responses2[indx]
          sub1 <- subset(renamed.tbl, metric == met & response == resp & gene.set == gset & train.set == "ohsu" & test.set == "fimm" & model == mdl & pt.response == pt.response1)
          sub2 <- subset(renamed.tbl, metric == met & response == resp & gene.set == gset & train.set == "ohsu" & test.set == "fimm" & model == mdl & pt.response == pt.response2)
          sub <- merge(sub1, sub2, by = "train.drug", suffixes = c(".1", ".2"))
          sub$diff <- sub$val.1 - sub$val.2

            write.table(file = paste(file.prefix, mdl, resp, met, gset, make.names(pt.response1), "vs", make.names(pt.response2), "all-points.tsv", sep="-"), sub, sep="\t", row.names=FALSE, col.names=TRUE)

##          flag <- !is.na(sub$val.1) & !is.na(sub$val.2) & ( abs(sub$val.1 - sub$val.2) > 0.1 )
          flag <- !is.na(sub$val.1) & !is.na(sub$val.2) & ( (sub$val.1 - sub$val.2) < -0.1 )
          flag <- flag | ( sub$train.drug == "Venetoclax" )
          sub$color <- "black"
          sub$color[flag] <- "blue"
          g <- ggplot(data = sub)
          g <- g + geom_abline(slope = 1, intercept = 0, linetype = "dashed")
          if(any(flag)) {
            labels <- sub[flag,]
            extremal.drugs <- merge(labels, drug.map[, c("OHSU_DRUG_NAME", "FIMM_DRUG_NAME")], by.x = "train.drug", by.y = "FIMM_DRUG_NAME")
            if(nrow(labels) != nrow(extremal.drugs)) {
              stop("Unexpected change in extremal drug table")
            }
            write.table(file = paste(file.prefix, mdl, resp, met, gset, make.names(pt.response1), "vs", make.names(pt.response2), "extremal-points.tsv", sep="-"), extremal.drugs, sep="\t", row.names=FALSE, col.names=TRUE)
          sub.large <- subset(sub, train.drug %in% c("Venetoclax"))
          sub.small <- subset(sub, !(train.drug %in% c("Venetoclax")))
          g <- g + geom_point(data = sub.small, aes(x = val.1, y = val.2, col = color), show.legend = FALSE)
          if(nrow(sub.large) >= 1) {
            g <- g + geom_point(data = sub.large, aes(x = val.1, y = val.2, col = color, size = 10), show.legend = FALSE)
          }
            nudge <- 0.025
set.seed(1234)
##            g <- g + geom_text(data = labels, aes(x = val.1, y = val.2, label = train.drug), position = position_jitter(width=nudge, height=nudge), hjust = "right", size = 8)
##            g <- g + geom_text(data = labels, aes(x = val.1, y = val.2, label = train.drug), nudge_x = - 0.02, hjust = "right", size = 6)
              left.labels <- labels
              right.labels <- NULL
              large.labels <- NULL
              if((mdl == "Ridge") && (met == "pearson")) {
                right.labels <- subset(labels, train.drug %in% c("XAV-939", "Selumetinib"))
              }
              large.labels <- subset(labels, train.drug %in% c("Venetoclax"))
              if(!is.null(right.labels)) {
                left.labels <- subset(left.labels, !(train.drug %in% right.labels$train.drug))
              }
              if(!is.null(large.labels)) {
                left.labels <- subset(left.labels, !(train.drug %in% large.labels$train.drug))
              }
              if(!is.null(left.labels) && (nrow(left.labels) >= 1)) { g <- g + geom_text(data = left.labels, aes(x = val.1, y = val.2, label = train.drug), nudge_x = - 0.02, hjust = "right", size = 6) }
              if(!is.null(large.labels) && (nrow(large.labels) >= 1)) { g <- g + geom_text(data = large.labels, aes(x = val.1, y = val.2, label = train.drug), nudge_x = - 0.02, hjust = "right", size = 10) }
              if(!is.null(right.labels) && (nrow(right.labels) >= 1)) { g <- g + geom_text(data = right.labels, aes(x = val.1, y = val.2, label = train.drug), nudge_x = 0.02, hjust = "left", size = 6) }
          }
##          g <- g + geom_point(data = sub, aes(x = val.1, y = val.2, col = color), show.legend = FALSE)
          lb <- min(min(c(sub$val.1, sub$val.2)) - 0.1, -0.2)
          g <- g + xlim(c(lb,1)) + ylim(c(lb,1))
          g <- g + xlab(capwords(paste(pt.response1, met, "correlation", sep=" ")))
          g <- g + ylab(capwords(paste(pt.response2, met, "correlation", sep=" ")))
          g <- g + theme(text = element_text(size = 20))
          ## NB: this will cutoff a point.  Do it to focus on points of interest.
          g <- g + xlim(c(-0.5, 1)) + ylim(c(-0.5, 1))
##          g <- g + ggtitle(paste0(resp, " ", capwords(met), " ", mdl, " ", gset))
          pdf(paste(file.prefix, mdl, resp, met, gset, make.names(pt.response1), "vs", make.names(pt.response2), "comparison.pdf", sep="-"))
          print(g)
          d <- dev.off()
        }
      }
    }
  } 
}

plot.observed.vs.expected.overlap <- function(trn.tbl, tbl.full, drugs, gene.set = "gene", model = "ridge", pt.response = "unpenalized.mean.feature", train.set = "ohsu", test.set = "fimm", metric = "pearson", universe, prefix) {
  gset <- gene.set
  mdl <- model
  resp <- pt.response
  training.set <- train.set
  validation.set <- test.set
  met <- metric

  tbl.tmp <- subset(trn.tbl, (gene.set == gset) & (model == mdl) & (pt.response == resp) & (train.drug %in% drugs))
  tbl.sub <- subset(tbl.full, (model == mdl) & (gene.set == gset) & (pt.response == resp) & (train.set == training.set) & (test.set == validation.set) & (metric == met) & (train.drug %in% drugs))
  colnames(tbl.sub)[colnames(tbl.sub) == "pval"] <- "corr.pval"
  colnames(tbl.sub)[colnames(tbl.sub) == "val"] <- "corr"

  mx <- length(universe[!grepl(pattern="mean", universe)])
  list.sizes <- c(seq(from=5,to=100, by=5), seq(from=100, to=500, by=10), seq(from=500, to=1000, by=25))
  list.sizes <- list.sizes[list.sizes <= mx]
  list.sizes <- list.sizes[list.sizes > 0]
  list.sizes <- unique(sort(list.sizes))
  names(list.sizes) <- list.sizes

  num.sigs <- ldply(list.sizes, .parallel = TRUE,
                    .fun = function(sz) {
                             tmp <- NULL
                             switch(mdl,
                                    "rf" = { tmp <- find.overlap.of.sparse.predictors(tbl.tmp, universe = universe, transform = "none", top = sz, verbose = FALSE, account.for.sign.of.coeff = FALSE) },
                                    "crf" = { tmp <- find.overlap.of.sparse.predictors(tbl.tmp, universe = universe, transform = "none", top = sz, verbose = FALSE, account.for.sign.of.coeff = FALSE) },
                                    "lasso" = { tmp <- find.overlap.of.sparse.predictors(tbl.tmp, universe = universe, transform = "abs.drop.zero", top = sz, verbose = FALSE, account.for.sign.of.coeff = TRUE) },
                                    "ridge" = { tmp <- find.overlap.of.sparse.predictors(tbl.tmp, universe = universe, transform = "abs.drop.zero", top = sz, verbose = FALSE, account.for.sign.of.coeff = TRUE) },
                                    { stop(paste0("Wasn't expecting ", mdl, "\n")) })
                             tmp <- tmp[, c("FIMM_DRUG_NAME", "train.set.1", "train.set.2", "n.1", "n.2", "n.both", "n.neither", "n.expected", "pval", "top", "gene.set.1"), drop = FALSE]
                             tmp
                           })
  tmp <- merge(subset(num.sigs, (train.set.1 == training.set) & (train.set.2 == validation.set)), tbl.sub, by.x = c("FIMM_DRUG_NAME", "gene.set.1"), by.y = c("test.drug", "gene.set"))
  tmp$top <- as.numeric(tmp$top)
  tmp$n.both <- as.numeric(tmp$n.both)
  tmp$n.expected <- as.numeric(tmp$n.expected)
  tmp$pval <- as.numeric(tmp$pval)
  tmp$corr.pval <- as.numeric(tmp$corr.pval)
  tmp <- tmp[order(tmp$corr.pval), , drop=FALSE]
  tmp$FIMM_DRUG_NAME <- factor(as.character(tmp$FIMM_DRUG_NAME), levels = unique(tmp$FIMM_DRUG_NAME))
  gs <- 
    dlply(tmp, .variables = c("FIMM_DRUG_NAME"), .parallel = FALSE,
            .fun = function(df) {
                     pdf(paste0(prefix, "-", df$FIMM_DRUG_NAME[1], ".pdf"))
                     flag <- !is.na(df$n.both)
                     corr <- df$corr[1]; corr.pval <- df$corr.pval[1]
                     df$frac <- df$n.both / df$top
                     df$frac.expected <- df$n.expected / df$top
                     df$Overlap <- "Observed"; df$frac.expected.color = "Expected"
                     g1 <- ggplot(data = df[flag,]) + geom_point(aes(x = top, y = frac, colour = Overlap))
                     g1 <- g1 + xlab("Number of Genes in Venetoclax Ridge Model\n(Ordered by Absolute Value)")
                     g1 <- g1 + ylab("Fraction of Overlapping Genes")
                     g1 <- g1 + geom_line(aes(x = top, y = frac.expected, colour = frac.expected.color), linetype = "dashed")
                     g1 <- g1 + theme(text = element_text(size = 20), legend.position = "top", legend.title.align = 0.5, legend.justification = "center")
                     g1 <- g1 + scale_colour_manual(values = c("black", "black"),labels = c("Observed", "Expected"),
                                                    guide = guide_legend(override.aes = list(linetype = c("blank", "dashed"), shape = c(16, NA))))
                     print(g1)
                     d <- dev.off()
                     list("g1" = g1)
            })
  list(gs = gs)
}

plot.observed.vs.expected.overlap(trn.tbl = trn.tbl, tbl.full = tbl.full, drugs = c("Venetoclax"), gene.set = "gene", model = "ridge", 
                                  pt.response = "unpenalized.mean.feature", train.set = "ohsu", test.set = "fimm",
                                  metric = "pearson", universe = universes[["gene"]], prefix = "fimm-vs-ohsu-ridge-gene-overlap")

