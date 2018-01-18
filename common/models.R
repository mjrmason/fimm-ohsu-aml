suppressPackageStartupMessages(library("hdi"))

get.fit.coeff <- function(fit, model) {

  if(model == "glmnet") {
    tmp <- coefficients(fit, s = "lambda.min")
    c1 <- as.vector(tmp)
    names(c1) <- rownames(tmp)
  } else {
    col <- colnames(importance(fit))[1]
    ## NB: %IncMSE = [ mse(j) - mse0 ] / mse0 * 100  [ mse(j): permuted; mse0: original MSE ]
    ## i.e., larger/more positive %IncMSE is better
    if(!grepl(col, pattern="IncMSE")) { stop(paste0("Was expected ", col, " to be %IncMSE\n")) }
    coeffs <- importance(fit)[,1,drop=FALSE]
    c1 <- coeffs[,1]
    names(c1) <- rownames(coeffs)
  }
  return(c1)
}

res.list.to.table <- function(res.list) {
  tbl.full <- ldply(names(res.list[["all.comparisons"]]),
               .fun = function(train.set) {
                        ldply(names(res.list[["all.comparisons"]][[train.set]]), 
                              .fun = function(test.set) {
                                       ldply(names(res.list[["all.comparisons"]][[train.set]][[test.set]]),
                                             .fun = function(gene.set) {
                                                      tmp.auc <- extract.results(res.list[["all.comparisons"]][[train.set]][[test.set]][[gene.set]])
                                                      tmp.auc$train.set <- train.set
                                                      tmp.auc$test.set <- test.set
                                                      tmp.auc$gene.set <- gene.set
                                                      tmp.auc$response <- "AUC"

                                                      tmp.features <- extract.features(res.list[["all.fits"]][[train.set]][[gene.set]], s = "lambda.min")
                                                      tmp.auc <- merge(tmp.auc, tmp.features, all.x = TRUE, all.y = FALSE)
                                                      return(tmp.auc)                         
                                                    }) 
                                     }) 
                      })
  tbl.full
}


common.markers <- function(res.auc, gene.set, ds1, drug1, ds2, drug2, alpha = 0, model = "glmnet", top = 20) {
  foo <- extract.fit(res.auc[["all.fits"]][[ds1]][[gene.set]], train.drug = drug1, alpha = alpha, model = model)
  bar <- extract.fit(res.auc[["all.fits"]][[ds2]][[gene.set]], train.drug = drug2, alpha = alpha, model = model)
  c1 <- get.fit.coeff(foo, model)
  c2 <- get.fit.coeff(bar, model)
  c1 <- c1[!grepl(pattern="Intercept", names(c1))]
  c2 <- c2[!grepl(pattern="Intercept", names(c2))]
  c2 <- c2[names(c1)]
  
  flag <- which(names(c1) %in% intersect(names(c1),names(c2)))
  if(top != 0) { 
    if(FALSE) {
    r1 <- rank(c1)
    r2 <- rank(c2)
    flag <- (r1 < top) & (r2 < top)

    r1 <- rank(-c1)
    r2 <- rank(-c2)
    flag <- flag | ( (r1 < top) & (r2 < top) )
    }
    switch(model,
      "glmnet" = {
        r1 <- rank(abs(c1))
        r2 <- rank(abs(c2))
      },
      "rf" = {
        r1 <- rank(-1 * c1)
        r2 <- rank(-1 * c2)
      },
      { stop(paste0("Unknown model: ", model, "\n")) })
    flag <- (r1 < top) & (r2 < top) & sign(c1) == sign(c2)
  }
  return(names(c1)[flag])
}

extract.table <- function(res.aucs, train.set.names, gene.sets, s = "lambda.min") {
 trn.tbl <- ldply(train.set.names,
             .fun = function(train.set) {
                      ldply(names(gene.sets),
                        .fun = function(gene.set) {
                                 ldply(res.aucs[[1]][["all.fits"]][[train.set]][[gene.set]],
                                   .fun = function(lst) {
                                     ldply(lst, 
                                       .fun = function(df) {
                                         ## Only pull out the regression and rf results
                                         if((df$model != "glmnet") && (df$model != "rf")) { return(NULL) }
                                         coeffs <- NULL
                                         if(is.null(df$fit) || is.na(df$fit)) { return(NULL) }
                                         switch(df$model,
                                           "glmnet" = {
                                             coeffs <- coefficients(df$fit, s = s)
                                           },
                                           "rf" = {
                                             col <- colnames(importance(df$fit))[1]
                                             ## NB: %IncMSE = [ mse(j) - mse0 ] / mse0 * 100  [ mse(j): permuted; mse0: original MSE ]
                                             ## i.e., larger/more positive %IncMSE is better
                                             if(!grepl(col, pattern="IncMSE")) { stop(paste0("Was expected ", col, " to be %IncMSE\n")) }
                                             coeffs <- importance(df$fit)[,1,drop=FALSE]
                                           },
                                           { stop(paste0("Unknown model ", df$model, "\n")) })
                                         ns <- rownames(coeffs)
                                         coeffs <- as.vector(coeffs)
                                         names(coeffs) <- ns
                                         coeffs <- coeffs[!grepl(pattern="Intercept", names(coeffs))]
##                                         n.features <- length(coeffs) 
n.features <- df$n.features
                                         coeffs <- coeffs[order(names(coeffs))]
                                         coeff.vals <- paste(coeffs, collapse=",")
                                         coeffs <- paste(names(coeffs), collapse=",")
                                         ret <- data.frame(train.set = train.set, gene.set = gene.set, train.drug = df$train.drug,
                                                           alpha = df$alpha, model = df$model, n.features = n.features, coeffs = coeffs, coeff.vals = coeff.vals)
                                         ret
                                       })
                                   })
                        })
             })
  trn.tbl$model <- as.character(trn.tbl$model)
  trn.tbl$train.drug <- as.character(trn.tbl$train.drug)
  trn.tbl$coeffs <- as.character(trn.tbl$coeffs)
  trn.tbl[(trn.tbl$model == "glmnet") & (trn.tbl$alpha == 0),"model"] <- "ridge"
  trn.tbl[(trn.tbl$model == "glmnet") & (trn.tbl$alpha == 1),"model"] <- "lasso"
  trn.tbl
}

extract.table.2 <- function(res.auc, s = "lambda.min") {
 trn.tbl <- ldply(names(res.auc[["all.fits"]]),
             .fun = function(train.set) {
                      ldply(names(res.auc[["all.fits"]][[train.set]]),
                        .fun = function(gene.set) {
                                 ldply(res.auc[["all.fits"]][[train.set]][[gene.set]],
                                   .fun = function(lst) {
                                     ldply(lst, 
                                       .fun = function(df) {
                                         ## Only pull out the regression and rf results
                                         if((df$model != "glmnet") && (df$model != "rf")) { return(NULL) }
                                         coeffs <- NULL
                                         if(is.null(df$fit) || is.na(df$fit)) { return(NULL) }
                                         switch(df$model,
                                           "glmnet" = {
                                             coeffs <- coefficients(df$fit, s = s)
                                           },
                                           "rf" = {
                                             col <- colnames(importance(df$fit))[1]
                                             ## NB: %IncMSE = [ mse(j) - mse0 ] / mse0 * 100  [ mse(j): permuted; mse0: original MSE ]
                                             ## i.e., larger/more positive %IncMSE is better
                                             if(!grepl(col, pattern="IncMSE")) { stop(paste0("Was expected ", col, " to be %IncMSE\n")) }
                                             coeffs <- importance(df$fit)[,1,drop=FALSE]
                                           },
                                           { stop(paste0("Unknown model ", df$model, "\n")) })
                                         ns <- rownames(coeffs)
                                         coeffs <- as.vector(coeffs)
                                         names(coeffs) <- ns
                                         coeffs <- coeffs[!grepl(pattern="Intercept", names(coeffs))]
##                                         n.features <- length(coeffs) 
n.features <- df$n.features
                                         coeffs <- coeffs[order(names(coeffs))]
                                         coeff.vals <- paste(coeffs, collapse=",")
                                         coeffs <- paste(names(coeffs), collapse=",")
                                         ret <- data.frame(train.set = train.set, gene.set = gene.set, train.drug = df$train.drug,
                                                           alpha = df$alpha, model = df$model, n.features = n.features, coeffs = coeffs, coeff.vals = coeff.vals)
                                         ret
                                       })
                                   })
                        })
             })
  trn.tbl$model <- as.character(trn.tbl$model)
  trn.tbl$train.drug <- as.character(trn.tbl$train.drug)
  trn.tbl$coeffs <- as.character(trn.tbl$coeffs)
  trn.tbl[(trn.tbl$model == "glmnet") & (trn.tbl$alpha == 0),"model"] <- "ridge"
  trn.tbl[(trn.tbl$model == "glmnet") & (trn.tbl$alpha == 1),"model"] <- "lasso"
  trn.tbl
}

plot.overlap.of.pathways <- function(tbl.resp, tbl.sub, universe, transform = "abs.drop.zero", list.sizes, ofile.prefix, main, pathway.set, expected.fracs = NULL) {

      cat("Calling find.overlap.of.sparse.predictors.simulated\n")
      flag <- !is.na(tbl.sub$corr) & !is.na(tbl.sub$corr.pval) & (tbl.sub$corr > 0) & (tbl.sub$corr.pval < 0.05)
      if(length(which(flag)) < 1) { return }
      keep <- unlist(apply(tbl.sub[flag,c("train.drug", "gene.set"),drop=F], 1, function(row) gsub(x=paste(row,collapse="-"), pattern=" ", replacement="")))
      tbl.sub <- tbl.sub[flag,,drop=F]
      keep.flag <- unlist(apply(tbl.resp[,c("train.drug", "gene.set")], 1, function(row) gsub(x=paste(row,collapse="-"), pattern=" ", replacement="") %in% keep))
      tbl.resp <- tbl.resp[keep.flag,,drop=F]
      if(nrow(tbl.resp) < 1) { return }
      num.sigs <- ldply(list.sizes, .parallel = TRUE,
                           .fun = function(sz) {
if(!(as.character(sz) %in% names(expected.fracs))) {
  save(expected.fracs, file="expected.fracs.prob.Rd")
  stop(paste0("Could not find ", sz, " in expected.fracs\n"))
}
                                    tmp <- find.overlap.of.sparse.predictors.simulated(tbl.resp, universe = universe, transform = transform, top = sz, pathway.set = pathway.set, verbose = FALSE, expected.fracs = expected.fracs[[as.character(sz)]])
                                    tmp <- tmp[, c("FIMM_DRUG_NAME", "train.set.1", "train.set.2", "frac.overlap", "exp.frac.overlap", "lb", "ub", "pval", "top", "gene.set.1")]
                                    tmp
                                  })
      cat("Done calling find.overlap.of.sparse.predictors.simulated\n")
      tmp <- merge(subset(num.sigs, (train.set.1 == "ohsu") & (train.set.2 == "fimm")), tbl.sub, by.x = c("FIMM_DRUG_NAME", "gene.set.1"), by.y = c("test.drug", "gene.set"))
      if(nrow(tmp) < 1) { return }
      tmp$top <- as.numeric(tmp$top)
      tmp$frac.overlap <- as.numeric(tmp$frac.overlap)
      tmp$exp.frac.overlap <- as.numeric(tmp$exp.frac.overlap)
##print(tmp[tmp$exp.frac.overlap > 0,])
      tmp$lb <- as.numeric(tmp$lb)
      tmp$ub <- as.numeric(tmp$ub)
      tmp$pval <- as.numeric(tmp$pval)
      tmp$corr.pval <- as.numeric(tmp$corr.pval)
      tmp <- tmp[order(tmp$corr.pval), ]
      flag <- !is.na(tmp$corr.pval) & !is.na(tmp$corr) & (tmp$corr.pval < 0.05) & (tmp$corr > 0)
##      flag <- !is.na(tmp$corr.pval) & !is.na(tmp$corr)
##      flag <- tmp$corr.pval < 0.05 & tmp$corr > 0
      if(!(any(flag))) { return }
      tmp <- tmp[flag,,drop=F]
      tmp$FIMM_DRUG_NAME <- factor(as.character(tmp$FIMM_DRUG_NAME), levels = unique(tmp$FIMM_DRUG_NAME))
      gs <- dlply(tmp, .variables = c("FIMM_DRUG_NAME"), .parallel = FALSE,
            .fun = function(df) {
                     flag <- !is.na(df$frac.overlap)
##if(TRUE) {
                     if(length(which(flag)) >= 2) { 
                       write.table(file = paste0(paste(ofile.prefix, df$FIMM_DRUG_NAME[1], sep="-"), "-ohsu-vs-fimm-top-feature-overlap.tsv"), df, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
                       pdf(paste0(paste(ofile.prefix, df$FIMM_DRUG_NAME[1], sep="-"), "-ohsu-vs-fimm-top-feature-overlap.pdf"))
                       corr <- df$corr[1]; corr.pval <- df$corr.pval[1]
                       df2 <- df[flag,,drop=F]
                       df2$value <- df2$frac.overlap
                       df2$Fraction <- rep("Observed", nrow(df2))
                       df2$colour = "black"
                       df3 <- df2
                       df3$value <- df3$exp.frac.overlap
                       df3$Fraction <- rep("Expected", nrow(df3))
                       df3$colour = "red"
                       df3 <- rbind(df2, df3)
##                       df3$colour <- factor(df3$colour, levels=c("red", "black"))
                       df3$colour <- factor(df3$colour, levels=c("black", "red"))
##                       g1 <- ggplot(data = df[flag,,drop=F])
                       g1 <- ggplot()
##                       g1 <- g1 + geom_errorbar(aes(x = top, ymin = lb, ymax = ub, colour = "red"))
##                       g1 <- g1 + geom_point(aes(x = top, y = exp.frac.overlap, colour = "red"))
##                       g1 <- g1 + geom_point(aes(x = top, y = frac.overlap)) + xlab("Num Top Features") + ylab("Fraction Overlapping")
##                       g1 <- g1 + geom_errorbar(data = df3[df3$Fraction == "Expected",,drop=F], aes(x = top, ymin = lb, ymax = ub, colour = colour))
                       g1 <- g1 + geom_errorbar(data = df3[df3$Fraction == "Expected",,drop=F], aes(x = top, ymin = lb, ymax = ub, colour = Fraction))
##                       g1 <- g1 + geom_point(data = df3, aes(x = top, y = value, colour = colour))
                       g1 <- g1 + geom_point(data = df3[df3$Fraction == "Expected",], aes(x = top, y = value, colour = Fraction))
                       g1 <- g1 + geom_point(data = df3[df3$Fraction == "Observed",], aes(x = top, y = value, colour = Fraction))
                       g1 <- g1 + xlab("Num Top Features") + ylab("Fraction Overlapping")
                       g1 <- g1 + ggtitle(paste0(paste(df$FIMM_DRUG_NAME[1], main, sep=" "), ": corr = ", signif(df$corr[1], digits=2), " pval = ", signif(df$corr.pval[1], digits=2)))
                       g1 <- g1 + theme(legend.position = "top")
                       g1 <- g1 + scale_color_manual(breaks = c("Expected", "Observed"), values=c("red", "black"))
##                       g1 <- g1 + scale_colour_manual(values=c("black", "red"), name="Fraction", labels=c("Observed", "Expected"))
                       g2 <- ggplot(data = df[flag,]) + geom_point(aes(x = top, y = -log10(pval))) + xlab("Num Top Features") + ylab("Overlap -log10 p-value")
                       g2 <- g2 + geom_hline(yintercept = -log10(0.05))
                       grid.arrange(g1, g2)
                       d <- dev.off()
                       return(list("g1" = g1, "g2" = g2))
                     }
                     return(list("g1" = NULL, "g2" = NULL))
            })
      pdf(paste0(paste(ofile.prefix, "all", sep="-"), "-ohsu-vs-fimm-top-feature-overlap.pdf"))
      for(i in 1:length(gs)) {
        if(!is.null(gs[[i]]$g1) && !is.null(gs[[i]]$g2)) {
          grid.arrange(gs[[i]]$g1, gs[[i]]$g2)
        }
      }
      d <- dev.off()

}

## transform = c("none", "abs", "negate")
## top = # of features to consider in overall
## top = 0 -> use all
## account.for.sign.of.coeff = TRUE -> only consider coeffs overlapping if they have the same sign
find.overlap.of.sparse.predictors <- function(tbl.input, universe, transform, top, gene.sets = NULL, verbose = FALSE, account.for.sign.of.coeff = TRUE) {
  ret.tbl <- c()
  for(gs in unique(tbl.input$gene.set)) {
    tmp <- subset(tbl.input, gene.set == gs)
    tr.sets <- as.character(unique(tmp$train.set))
    for(i in 1:(length(tr.sets)-1)) {
      tr1 <- tr.sets[i]
      tmp1 <- subset(tmp, train.set == tr1)
      merge.col <- drug.name.col
      tmp1 <- merge(drug.name.tbl, tmp1, by.x = merge.col, by.y = "train.drug")
      for(j in (i+1):length(tr.sets)) {
        tr2 <- tr.sets[j]
        if(((tr1 == "ohsu") || (tr2 == "ohsu")) && (grepl(pattern="ohsu.set1", paste0(tr1,tr2)) || grepl(pattern="ohsu.set2", paste0(tr1,tr2)))) { next }
        if(((tr1 == "fimm") || (tr2 == "fimm")) && (grepl(pattern="fimm.set1", paste0(tr1,tr2)) || grepl(pattern="fimm.set2", paste0(tr1,tr2)))) { next }
        tmp2 <- subset(tmp, train.set == tr2)
        merge.col <- drug.name.col
        tmp2 <- merge(drug.name.tbl, tmp2, by.x = merge.col, by.y = "train.drug")
        tmp1$id.1 <- 1:nrow(tmp1)
        m <- merge(tmp1, tmp2, by = merge.col, suffixes = c(".1", ".2"))
        m <- m[, c(merge.col, "id.1", "n.features.1", "n.features.2", "coeffs.1", "coeffs.2", "coeff.vals.1", "coeff.vals.2", "train.set.1", "train.set.2", "gene.set.1")]

        m$id <- 1:nrow(m)
        hyp <- ddply(m, .variables = "id.1", .parallel = FALSE,
                     .fun = function(r1) {
drug.name <- r1[1,drug.name.col]
                        n.f1 <- as.numeric(r1$n.features.1)
                        coeff.vals.1 <- as.numeric(unlist(strsplit(as.character(r1$coeff.vals.1), split=",")))
                        coeffs.1 <- (unlist(strsplit(as.character(r1$coeffs.1), split=",")))
                        names(coeff.vals.1) <- coeffs.1
                        coeffs.1 <- coeff.vals.1
                        ## append "-" to the name if coeff is negative
                        if(account.for.sign.of.coeff) {
                          names(coeffs.1) <- unlist(lapply(1:length(coeffs.1), function(ind) { ifelse(is.na(coeffs.1[ind]) || (coeffs.1[ind] >= 0), names(coeffs.1)[ind], paste0("-", names(coeffs.1)[ind])) }))
                        }
##save(coeffs.1, file=paste0(drug.name, "-coeffs.1.Rd"))
                        switch(transform,
                          "drop.zero" = { 
                             coeffs.1 <- coeffs.1[coeffs.1 > 0]
                          },
                          "none" = { },
                          "abs" = {
                             coeffs.1 <- abs(coeffs.1)
                          },
                          "abs.drop.zero" = {
                             coeffs.1 <- abs(coeffs.1)
                             coeffs.1 <- coeffs.1[coeffs.1 != 0]
                          },
                          "negate" = {
                             coeffs.1 <- - coeffs.1
                          })
                        coeffs.1 <- sort(coeffs.1, decreasing=TRUE)
##save(coeffs.1, file=paste0(drug.name, "-sorted-coeffs.1.Rd"))
                        ddply(r1, .variables = "id", .parallel = FALSE,
                              .fun = function(r) {
                        n.f2 <- as.numeric(r$n.features.2)
                                 coeff.vals.2 <- as.numeric(unlist(strsplit(as.character(r$coeff.vals.2), split=",")))
                                 coeffs.2 <- (unlist(strsplit(as.character(r$coeffs.2), split=",")))
                                 names(coeff.vals.2) <- coeffs.2
                                 coeffs.2 <- coeff.vals.2
                                 if(account.for.sign.of.coeff) {
                                   names(coeffs.2) <- unlist(lapply(1:length(coeffs.2), function(ind) { ifelse(is.na(coeffs.2[ind]) || (coeffs.2[ind] >= 0), names(coeffs.2)[ind], paste0("-", names(coeffs.2)[ind])) }))
                                 }
##save(coeffs.2, file=paste0(drug.name, "-coeffs.2.Rd"))
                                 switch(transform,
                                   "drop.zero" = { 
                                      coeffs.2 <- coeffs.2[coeffs.2 > 0]
                                   },
                                   "none" = { },
                                   "abs" = {
                                      coeffs.2 <- abs(coeffs.2)
                                   },
                                   "abs.drop.zero" = {
                                      coeffs.2 <- abs(coeffs.2)
                                      coeffs.2 <- coeffs.2[coeffs.2 != 0]
                                   },
                                   "negate" = {
                                      coeffs.2 <- - coeffs.2
                                   })
                                 coeffs.2 <- sort(coeffs.2, decreasing=TRUE)
##save(coeffs.2, file=paste0(drug.name, "-sorted-coeffs.2.Rd"))
                                 ldply(top, .parallel = FALSE,
                                       .fun = function(sz) {
                                                c1 <- coeffs.1
                                                c2 <- coeffs.2
                                                if(length(c1) > 0) {
                                                  c1 <- c1[!duplicated(names(c1))]
                                                }
                                                if(length(c2) > 0) {
                                                  c2 <- c2[!duplicated(names(c2))]
                                                }
c1 <- c1[is.finite(c1)]
c2 <- c2[is.finite(c2)]
## Drop any mean response
c1 <- c1[!grepl(pattern="mean", names(c1))]
c2 <- c2[!grepl(pattern="mean", names(c2))]
if(any(grepl(pattern="mean", universe))) {
n.f1 <- n.f1 - 1
n.f2 <- n.f2 - 1
}
## Because sign matters, need to double the number of features in the universe
if(account.for.sign.of.coeff) {
  n.f1 <- 2 * n.f1
  n.f2 <- 2 * n.f2
}

if((sz > length(c1)) || (sz > length(c2))) {
                                                vec <- c(r$id[1], length(c1), length(c2), NA, NA, NA, NA, NA, sz)
                                                names(vec) <- c("id", "n.1", "n.2", "n.both", "n.neither", "n.expected", "both", "pval", "top")
                                                return(vec)
}

                                                c1.len <- min(sz, length(c1))
                                                if((sz != 0) && (length(c1) > 0)) {
                                                  c1 <- c1[1:c1.len]
                                                }

                                                
                                                c2.len <- min(sz, length(c2))
                                                if((sz != 0) && (length(c2) > 0)) {
                                                  c2 <- c2[1:c2.len]
                                                }
                                                n.1 <- length(c1)
                                                n.2 <- length(c2)
                                                both <- intersect(names(c1), names(c2))
                                                n.both <- length(both)
                                                both <- paste(both, collapse=",")
                                                ## pval <- phyper(q = n.both - 1, m = n.1, n = n.f1 - n.1, k = n.2, lower.tail = FALSE) 
## if(n.f1 != n.f2) { stop(paste0(r[,drug.name.col], ": n.f1 (", n.f1, ") != n.f2 (", n.f2, ")\n")) }
n.expected <- c1.len * c2.len / max(n.f1, n.f2)

##if(any(grepl(pattern="mean", names(c1)))) {
##  stop(paste0(paste(names(c1)[grepl(pattern="mean", names(c1))], collapse=", "), "\n"))
##}
##if(any(grepl(pattern="mean", names(c2)))) {
##  stop(paste0(paste(names(c2)[grepl(pattern="mean", names(c2))], collapse=", "), "\n"))
##}
n.1 <- length(setdiff(names(c1),names(c2)))
n.2 <- length(setdiff(names(c2),names(c1)))
## n.neither <- length(setdiff(universe, union(names(c1),names(c2))))
n.neither <- max(n.f1, n.f2) - n.1 - n.2 - n.both
if(n.neither < 0) { stop(paste0("n.neither < 0: ", paste(n.neither, sz, c1.len, c2.len, n.f1, n.f2, n.1, n.2, n.both, sep=" "), "\n")) }
if(n.neither < 0) {
stop(paste0(paste(setdiff(names(c2),names(c1)), collapse=","), "\n"))
stop(paste0(paste(setdiff(names(c1),names(c2)), collapse=","), "\n"))
stop(paste0(both, "\n"))
}
mat <- cbind(c(n.both, n.1), c(n.2, n.neither))
##                                                mat <- cbind(c(n.both, n.1 - n.both), c(n.2 - n.both, n.f1 - n.1))
                                                fet <- fisher.test(mat, alternative = "greater")
                                                pval <- fet$p.value
if(verbose) {
  print(mat)
  print(fet)
  print(pval)
}
                                                vec <- c(r$id[1], n.1, n.2, n.both, n.neither, n.expected, both, pval, sz)
                                                names(vec) <- c("id", "n.1", "n.2", "n.both", "n.neither", "n.expected", "both", "pval", "top")
                                                vec
                                  })
                               })
                             })
        hyp <- hyp[, !grepl(pattern="id.1", x=colnames(hyp))]
        m <- m[, !grepl(pattern="id.1", x=colnames(m))]
        hyp <- merge(hyp, m, by = "id")
        ret.tbl <- rbind(ret.tbl, hyp)
      }
    }
  }
  ret.tbl
}

## transform = c("none", "abs", "negate")
## top = # of features to consider in overall
## top = 0 -> use all
calc.expected.overlap.of.pathways.old <- function(tbl.input, universe, transform, top, pathway.set, verbose = FALSE) {
  expected.fracs <- list()
  for(gs in unique(tbl.input$gene.set)) {
    tmp <- subset(tbl.input, gene.set == gs)
    tr.sets <- as.character(unique(tmp$train.set))
    for(i in 1:(length(tr.sets)-1)) {
      tr1 <- tr.sets[i]
      tmp1 <- subset(tmp, train.set == tr1)
      merge.col <- drug.name.col
      tmp1 <- merge(drug.name.tbl, tmp1, by.x = merge.col, by.y = "train.drug")
      for(j in (i+1):length(tr.sets)) {
        tr2 <- tr.sets[j]
        if(((tr1 == "ohsu") || (tr2 == "ohsu")) && (grepl(pattern="ohsu.set1", paste0(tr1,tr2)) || grepl(pattern="ohsu.set2", paste0(tr1,tr2)))) { next }
        if(((tr1 == "fimm") || (tr2 == "fimm")) && (grepl(pattern="fimm.set1", paste0(tr1,tr2)) || grepl(pattern="fimm.set2", paste0(tr1,tr2)))) { next }
        tmp2 <- subset(tmp, train.set == tr2)
        merge.col <- drug.name.col
        tmp2 <- merge(drug.name.tbl, tmp2, by.x = merge.col, by.y = "train.drug")
        tmp1$id.1 <- 1:nrow(tmp1)
        m <- merge(tmp1, tmp2, by = merge.col, suffixes = c(".1", ".2"))
        m <- m[, c(merge.col, "id.1", "n.features.1", "n.features.2", "coeffs.1", "coeffs.2", "coeff.vals.1", "coeff.vals.2", "train.set.1", "train.set.2", "gene.set.1")]
        m$id <- 1:nrow(m)
        hyp <- ddply(m, .variables = "id.1", .parallel = FALSE,
                     .fun = function(r1) {
                        n.f1 <- as.numeric(r1$n.features.1)
                        coeff.vals.1 <- as.numeric(unlist(strsplit(as.character(r1$coeff.vals.1), split=",")))
                        coeffs.1 <- (unlist(strsplit(as.character(r1$coeffs.1), split=",")))
                        names(coeff.vals.1) <- coeffs.1
                        coeffs.1 <- coeff.vals.1
                        ## append "-" to the name if coeff is negative
## not clear how to handle this since we map to pathways
##                        names(coeffs.1) <- unlist(lapply(1:length(coeffs.1), function(ind) { ifelse(is.na(coeffs.1[ind]) || (coeffs.1[ind] >= 0), names(coeffs.1)[ind], paste0("-", names(coeffs.1)[ind])) }))
                        switch(transform,
                          "drop.zero" = { 
                             coeffs.1 <- coeffs.1[coeffs.1 > 0]
                          },
                          "none" = { },
                          "abs" = {
                             coeffs.1 <- abs(coeffs.1)
                          },
                          "abs.drop.zero" = {
                             coeffs.1 <- abs(coeffs.1)
                             coeffs.1 <- coeffs.1[coeffs.1 != 0]
                          },
                          "negate" = {
                             coeffs.1 <- - coeffs.1
                          })
                          coeffs.1 <- sort(coeffs.1, decreasing=TRUE)
drug.name <- r1[1,drug.name.col]
                        ddply(r1, .variables = "id", .parallel = FALSE,
                              .fun = function(r) {
                        n.f2 <- as.numeric(r$n.features.2)
                                 coeff.vals.2 <- as.numeric(unlist(strsplit(as.character(r$coeff.vals.2), split=",")))
                                 coeffs.2 <- (unlist(strsplit(as.character(r$coeffs.2), split=",")))
                                 names(coeff.vals.2) <- coeffs.2
                                 coeffs.2 <- coeff.vals.2
##                                 names(coeffs.2) <- unlist(lapply(1:length(coeffs.2), function(ind) { ifelse(is.na(coeffs.2[ind]) || (coeffs.2[ind] >= 0), names(coeffs.2)[ind], paste0("-", names(coeffs.2)[ind])) }))
                                 switch(transform,
                                   "drop.zero" = { 
                                      coeffs.2 <- coeffs.2[coeffs.2 > 0]
                                   },
                                   "none" = { },
                                   "abs" = {
                                      coeffs.2 <- abs(coeffs.2)
                                   },
                                   "abs.drop.zero" = {
                                      coeffs.2 <- abs(coeffs.2)
                                      coeffs.2 <- coeffs.2[coeffs.2 != 0]
                                   },
                                   "negate" = {
                                      coeffs.2 <- - coeffs.2
                                   })
                                 coeffs.2 <- sort(coeffs.2, decreasing=TRUE)
                                 ldply(top, .parallel = FALSE,
                                       .fun = function(sz) {
                                                c1 <- coeffs.1
                                                c2 <- coeffs.2
                                                if(length(c1) > 0) {
                                                  c1 <- c1[!duplicated(names(c1))]
                                                }
                                                if(length(c2) > 0) {
                                                  c2 <- c2[!duplicated(names(c2))]
                                                }
c1 <- c1[is.finite(c1)]
c2 <- c2[is.finite(c2)]
n1 <- min(length(c1), sz)
n2 <- min(length(c2), sz)
## Drop any mean response
c1 <- c1[!grepl(pattern="mean", names(c1))]
c2 <- c2[!grepl(pattern="mean", names(c2))]

if((sz > length(c1)) || (sz > length(c2))) {
                                                vec <- c(r$id[1], NA, NA, NA, NA, NA, sz) 
                                                names(vec) <- c("id", "frac.overlap", "exp.frac.overlap", "lb", "ub", "pval", "top")

                                                return(vec)
}


                                                c1.len <- min(sz, length(c1))
                                                if((sz != 0) && (length(c1) > 0)) {
                                                  c1 <- c1[1:c1.len]
                                                }

                                                
                                                c2.len <- min(sz, length(c2))
                                                if((sz != 0) && (length(c2) > 0)) {
                                                  c2 <- c2[1:c2.len]
                                                }
                                                c1 <- names(c1)
                                                c2 <- names(c2)
                                                if(any(c1 %in% names(pathway.set))) {
                                                  c1 <- unique(unlist(lapply(c1[c1 %in% names(pathway.set)], function(gene) pathway.set[[gene]])))
##                                                  c1 <- unique(unlist(pathway.set[[c1[c1 %in% names(pathway.set)]]]))
                                                } else { c1 <- "" }
                                                if(any(c2 %in% names(pathway.set))) {
                                                  c2 <- unique(unlist(lapply(c2[c2 %in% names(pathway.set)], function(gene) pathway.set[[gene]])))
##                                                  c2 <- unique(unlist(pathway.set[[c2[c2 %in% names(pathway.set)]]]))
                                                } else { c2 <- "" }
                                                c1 <- c1[!is.null(c1) & !is.na(c1) & (c1 != "")]
                                                c2 <- c2[!is.null(c2) & !is.na(c2) & (c2 != "")]
                                                obs.frac.overlap <- 0
                                                both <- intersect(c1, c2)
                                                if((length(c1) > 0) && (length(c2) > 0)) {
                                                  obs.frac.overlap <- length(both)/min(length(c1), length(c2))
                                                }

                                                n.iters <- 1000
set.seed(1234)
cat("Calculating frac overlaps\n")
                                                  expected.fracs[[drug.name]] <- unlist(llply(1:n.iters, .parallel = FALSE,
                                                                              .fun = function(i) {
                                                                                       i1 <- sample.int(length(universe), size = n1, replace = FALSE) 
                                                                                       i2 <- sample.int(length(universe), size = n2, replace = FALSE) 
                                                                                       c1 <- universe[i1]
                                                                                       c2 <- universe[i2]
                                                if(any(c1 %in% names(pathway.set))) {
                                                  c1 <- unique(unlist(lapply(c1[c1 %in% names(pathway.set)], function(gene) pathway.set[[gene]])))
##                                                  c1 <- unique(unlist(pathway.set[[c1[c1 %in% names(pathway.set)]]]))
                                                } else { c1 <- "" }
                                                if(any(c2 %in% names(pathway.set))) {
                                                  c2 <- unique(unlist(lapply(c2[c2 %in% names(pathway.set)], function(gene) pathway.set[[gene]])))
##                                                  c2 <- unique(unlist(pathway.set[[c2[c2 %in% names(pathway.set)]]]))
                                                } else { c2 <- "" }
                                                c1 <- c1[!is.null(c1) & !is.na(c1) & (c1 != "")]
                                                c2 <- c2[!is.null(c2) & !is.na(c2) & (c2 != "")]

                                                                                       frac <- 0
                                                                                       if((length(c1) > 0) && (length(c2) > 0)) {

                                                                                         frac <- length(intersect(c1,c2))/min(length(c1), length(c2))
                                                                                       }
                                                                                       frac
                                                                                     }))
                                                
                                  })
                               })
                             })
      }
    }
  }
  expected.fracs
}

calc.expected.overlap.of.pathways <- function(universe, transform, top, pathway.set, verbose = FALSE) {
  n.iters <- 10000
  set.seed(1234)
  cat("Calculating frac overlaps\n")
  iter.vec <- 1:n.iters
  names(iter.vec) <- iter.vec
  expected.fracs <- unlist(llply(iter.vec, .parallel = FALSE,
                           .fun = function(i) {
                                    i1 <- sample.int(length(universe), size = top, replace = FALSE) 
                                    i2 <- sample.int(length(universe), size = top, replace = FALSE) 
                                    c1 <- universe[i1]
                                    c2 <- universe[i2]
                                    if(any(c1 %in% names(pathway.set))) {
                                      c1 <- unique(unlist(lapply(c1[c1 %in% names(pathway.set)], function(gene) pathway.set[[gene]])))
##                                    c1 <- unique(unlist(pathway.set[[c1[c1 %in% names(pathway.set)]]]))
                                    } else { c1 <- "" }
                                    if(any(c2 %in% names(pathway.set))) {
                                      c2 <- unique(unlist(lapply(c2[c2 %in% names(pathway.set)], function(gene) pathway.set[[gene]])))
##                                    c2 <- unique(unlist(pathway.set[[c2[c2 %in% names(pathway.set)]]]))
                                    } else { c2 <- "" }
                                    c1 <- c1[!is.null(c1) & !is.na(c1) & (c1 != "")]
                                    c2 <- c2[!is.null(c2) & !is.na(c2) & (c2 != "")]

                                    frac <- 0
                                    if((length(c1) > 0) && (length(c2) > 0)) {
                                      frac <- length(intersect(c1,c2))/min(length(c1), length(c2))
                                    }
                                    frac
                                  }))
  expected.fracs
}

find.overlap.of.sparse.predictors.simulated <- function(tbl.input, universe, transform, top, pathway.set, verbose = FALSE, expected.fracs = NULL) {
  ret.tbl <- c()
  for(gs in unique(tbl.input$gene.set)) {
    tmp <- subset(tbl.input, gene.set == gs)
    tr.sets <- as.character(unique(tmp$train.set))
    for(i in 1:(length(tr.sets)-1)) {
      tr1 <- tr.sets[i]
      tmp1 <- subset(tmp, train.set == tr1)
      merge.col <- drug.name.col
      tmp1 <- merge(drug.name.tbl, tmp1, by.x = merge.col, by.y = "train.drug")
      for(j in (i+1):length(tr.sets)) {
        tr2 <- tr.sets[j]
        if(((tr1 == "ohsu") || (tr2 == "ohsu")) && (grepl(pattern="ohsu.set1", paste0(tr1,tr2)) || grepl(pattern="ohsu.set2", paste0(tr1,tr2)))) { next }
        if(((tr1 == "fimm") || (tr2 == "fimm")) && (grepl(pattern="fimm.set1", paste0(tr1,tr2)) || grepl(pattern="fimm.set2", paste0(tr1,tr2)))) { next }
        tmp2 <- subset(tmp, train.set == tr2)
        merge.col <- drug.name.col
        tmp2 <- merge(drug.name.tbl, tmp2, by.x = merge.col, by.y = "train.drug")
        tmp1$id.1 <- 1:nrow(tmp1)
        m <- merge(tmp1, tmp2, by = merge.col, suffixes = c(".1", ".2"))
        m <- m[, c(merge.col, "id.1", "n.features.1", "n.features.2", "coeffs.1", "coeffs.2", "coeff.vals.1", "coeff.vals.2", "train.set.1", "train.set.2", "gene.set.1")]
        m$id <- 1:nrow(m)
        hyp <- ddply(m, .variables = "id.1", .parallel = FALSE,
                     .fun = function(r1) {
                        n.f1 <- as.numeric(r1$n.features.1)
                        coeff.vals.1 <- as.numeric(unlist(strsplit(as.character(r1$coeff.vals.1), split=",")))
                        coeffs.1 <- (unlist(strsplit(as.character(r1$coeffs.1), split=",")))
                        names(coeff.vals.1) <- coeffs.1
                        coeffs.1 <- coeff.vals.1
                        ## append "-" to the name if coeff is negative
## not clear how to handle this since we map to pathways
##                        names(coeffs.1) <- unlist(lapply(1:length(coeffs.1), function(ind) { ifelse(is.na(coeffs.1[ind]) || (coeffs.1[ind] >= 0), names(coeffs.1)[ind], paste0("-", names(coeffs.1)[ind])) }))
                        switch(transform,
                          "drop.zero" = { 
                             coeffs.1 <- coeffs.1[coeffs.1 > 0]
                          },
                          "none" = { },
                          "abs" = {
                             coeffs.1 <- abs(coeffs.1)
                          },
                          "abs.drop.zero" = {
                             coeffs.1 <- abs(coeffs.1)
                             coeffs.1 <- coeffs.1[coeffs.1 != 0]
                          },
                          "negate" = {
                             coeffs.1 <- - coeffs.1
                          })
                          coeffs.1 <- sort(coeffs.1, decreasing=TRUE)
drug.name <- r1[1,drug.name.col]
                        ddply(r1, .variables = "id", .parallel = FALSE,
                              .fun = function(r) {
                        n.f2 <- as.numeric(r$n.features.2)
                                 coeff.vals.2 <- as.numeric(unlist(strsplit(as.character(r$coeff.vals.2), split=",")))
                                 coeffs.2 <- (unlist(strsplit(as.character(r$coeffs.2), split=",")))
                                 names(coeff.vals.2) <- coeffs.2
                                 coeffs.2 <- coeff.vals.2
##                                 names(coeffs.2) <- unlist(lapply(1:length(coeffs.2), function(ind) { ifelse(is.na(coeffs.2[ind]) || (coeffs.2[ind] >= 0), names(coeffs.2)[ind], paste0("-", names(coeffs.2)[ind])) }))
                                 switch(transform,
                                   "drop.zero" = { 
                                      coeffs.2 <- coeffs.2[coeffs.2 > 0]
                                   },
                                   "none" = { },
                                   "abs" = {
                                      coeffs.2 <- abs(coeffs.2)
                                   },
                                   "abs.drop.zero" = {
                                      coeffs.2 <- abs(coeffs.2)
                                      coeffs.2 <- coeffs.2[coeffs.2 != 0]
                                   },
                                   "negate" = {
                                      coeffs.2 <- - coeffs.2
                                   })
                                 coeffs.2 <- sort(coeffs.2, decreasing=TRUE)
                                 ldply(top, .parallel = FALSE,
                                       .fun = function(sz) {
                                                c1 <- coeffs.1
                                                c2 <- coeffs.2
                                                if(length(c1) > 0) {
                                                  c1 <- c1[!duplicated(names(c1))]
                                                }
                                                if(length(c2) > 0) {
                                                  c2 <- c2[!duplicated(names(c2))]
                                                }
c1 <- c1[is.finite(c1)]
c2 <- c2[is.finite(c2)]
n1 <- min(length(c1), sz)
n2 <- min(length(c2), sz)
## Drop any mean response
c1 <- c1[!grepl(pattern="mean", names(c1))]
c2 <- c2[!grepl(pattern="mean", names(c2))]

if((sz > length(c1)) || (sz > length(c2))) {
                                                vec <- c(r$id[1], NA, NA, NA, NA, NA, sz) 
                                                names(vec) <- c("id", "frac.overlap", "exp.frac.overlap", "lb", "ub", "pval", "top")

                                                return(vec)
}


                                                c1.len <- min(sz, length(c1))
                                                if((sz != 0) && (length(c1) > 0)) {
                                                  c1 <- c1[1:c1.len]
                                                }

                                                
                                                c2.len <- min(sz, length(c2))
                                                if((sz != 0) && (length(c2) > 0)) {
                                                  c2 <- c2[1:c2.len]
                                                }
                                                c1 <- names(c1)
                                                c2 <- names(c2)
                                                if(any(c1 %in% names(pathway.set))) {
                                                  c1 <- unique(unlist(lapply(c1[c1 %in% names(pathway.set)], function(gene) pathway.set[[gene]])))
##                                                  c1 <- unique(unlist(pathway.set[[c1[c1 %in% names(pathway.set)]]]))
                                                } else { c1 <- "" }
                                                if(any(c2 %in% names(pathway.set))) {
                                                  c2 <- unique(unlist(lapply(c2[c2 %in% names(pathway.set)], function(gene) pathway.set[[gene]])))
##                                                  c2 <- unique(unlist(pathway.set[[c2[c2 %in% names(pathway.set)]]]))
                                                } else { c2 <- "" }
                                                c1 <- c1[!is.null(c1) & !is.na(c1) & (c1 != "")]
                                                c2 <- c2[!is.null(c2) & !is.na(c2) & (c2 != "")]
                                                obs.frac.overlap <- 0
                                                both <- intersect(c1, c2)
                                                if((length(c1) > 0) && (length(c2) > 0)) {
                                                  obs.frac.overlap <- length(both)/min(length(c1), length(c2))
                                                }

                                                n.iters <- 1000
set.seed(1234)
                                                frac.overlaps <- NULL
                                                if(!is.null(expected.fracs)) {
cat("Using memoized frac overlaps\n")
                                                  frac.overlaps <- expected.fracs
                                                } else {
cat("Calculating frac overlaps\n")
                                                  frac.overlaps <- unlist(llply(1:n.iters, .parallel = TRUE,
                                                                              .fun = function(i) {
                                                                                       i1 <- sample.int(length(universe), size = n1, replace = FALSE) 
                                                                                       i2 <- sample.int(length(universe), size = n2, replace = FALSE) 
                                                                                       c1 <- universe[i1]
                                                                                       c2 <- universe[i2]
                                                if(any(c1 %in% names(pathway.set))) {
                                                  c1 <- unique(unlist(lapply(c1[c1 %in% names(pathway.set)], function(gene) pathway.set[[gene]])))
##                                                  c1 <- unique(unlist(pathway.set[[c1[c1 %in% names(pathway.set)]]]))
                                                } else { c1 <- "" }
                                                if(any(c2 %in% names(pathway.set))) {
                                                  c2 <- unique(unlist(lapply(c2[c2 %in% names(pathway.set)], function(gene) pathway.set[[gene]])))
##                                                  c2 <- unique(unlist(pathway.set[[c2[c2 %in% names(pathway.set)]]]))
                                                } else { c2 <- "" }
                                                c1 <- c1[!is.null(c1) & !is.na(c1) & (c1 != "")]
                                                c2 <- c2[!is.null(c2) & !is.na(c2) & (c2 != "")]

                                                                                       frac <- 0
                                                                                       if((length(c1) > 0) && (length(c2) > 0)) {

                                                                                         frac <- length(intersect(c1,c2))/min(length(c1), length(c2))
                                                                                       }
                                                                                       frac
                                                                                     }))
                                                }
                                                exp.frac.overlap <- mean(frac.overlaps)
                                                qs <- quantile(frac.overlaps, probs=c(0.05, 0.95))
##print(c(obs.frac.overlap, exp.frac.overlap, frac.overlaps))
                                                lb.frac.overlap <- as.numeric(qs[1])
                                                ub.frac.overlap <- as.numeric(qs[2])
                                                pval <- length(which(frac.overlaps >= obs.frac.overlap))/length(frac.overlaps)
                                                vec <- c(r$id[1], obs.frac.overlap, exp.frac.overlap, lb.frac.overlap, ub.frac.overlap, pval, sz) 
                                                names(vec) <- c("id", "frac.overlap", "exp.frac.overlap", "lb", "ub", "pval", "top")
                                                vec
                                  })
                               })
                             })
        hyp <- hyp[, !grepl(pattern="id.1", x=colnames(hyp))]
        m <- m[, !grepl(pattern="id.1", x=colnames(m))]
        hyp <- merge(hyp, m, by = "id")
        ret.tbl <- rbind(ret.tbl, hyp)
      }
    }
  }
  ret.tbl
}


## transform = c("none", "abs", "negate")
## top = # of features to consider in overall
## top = 0 -> use all
old.find.overlap.of.sparse.predictors <- function(tbl.input, transform, top) {
  ret.tbl <- c()
  for(gs in unique(tbl.input$gene.set)) {
    tmp <- subset(tbl.input, gene.set == gs)
    tr.sets <- as.character(unique(tmp$train.set))
    for(i in 1:(length(tr.sets)-1)) {
      tr1 <- tr.sets[i]
      tmp1 <- subset(tmp, train.set == tr1)
      merge.col <- drug.name.col
      tmp1 <- merge(drug.name.tbl, tmp1, by.x = merge.col, by.y = "train.drug")
      for(j in (i+1):length(tr.sets)) {
        tr2 <- tr.sets[j]
        if(((tr1 == "ohsu") || (tr2 == "ohsu")) && (grepl(pattern="ohsu.set1", paste0(tr1,tr2)) || grepl(pattern="ohsu.set2", paste0(tr1,tr2)))) { next }
        if(((tr1 == "fimm") || (tr2 == "fimm")) && (grepl(pattern="fimm.set1", paste0(tr1,tr2)) || grepl(pattern="fimm.set2", paste0(tr1,tr2)))) { next }
        tmp2 <- subset(tmp, train.set == tr2)
        merge.col <- drug.name.col
        tmp2 <- merge(drug.name.tbl, tmp2, by.x = merge.col, by.y = "train.drug")

        m <- merge(tmp1, tmp2, by = merge.col, suffixes = c(".1", ".2"))
        m <- m[, c(merge.col, "n.features.1", "coeffs.1", "coeffs.2", "coeff.vals.1", "coeff.vals.2", "train.set.1", "train.set.2", "gene.set.1")]
        m$id <- 1:nrow(m)
        hyp <- ddply(m, .variables = "id", 
               .fun = function(r) {
                        n.f <- max(as.numeric(r$n.features.1), as.numeric(r$n.features.2))
                        coeff.vals.1 <- as.numeric(unlist(strsplit(as.character(r$coeff.vals.1), split=",")))
                        coeffs.1 <- (unlist(strsplit(as.character(r$coeffs.1), split=",")))
                        names(coeff.vals.1) <- coeffs.1
                        coeffs.1 <- coeff.vals.1
                        coeff.vals.2 <- as.numeric(unlist(strsplit(as.character(r$coeff.vals.2), split=",")))
                        coeffs.2 <- (unlist(strsplit(as.character(r$coeffs.2), split=",")))
                        names(coeff.vals.2) <- coeffs.2
                        coeffs.2 <- coeff.vals.2
                        switch(transform,
                          "drop.zero" = { 
                             coeffs.1 <- coeffs.1[coeffs.1 != 0]
                             coeffs.2 <- coeffs.2[coeffs.2 != 0]
                          },
                          "none" = { },
                          "abs.drop.zero" = {
                             coeffs.1 <- abs(coeffs.1)
                             coeffs.2 <- abs(coeffs.2)
                             coeffs.1 <- coeffs.1[coeffs.1 != 0]
                             coeffs.2 <- coeffs.2[coeffs.2 != 0]
                          },
                          "abs" = {
                             coeffs.1 <- abs(coeffs.1)
                             coeffs.2 <- abs(coeffs.2)
                          },
                          "negate" = {
                             coeffs.1 <- - coeffs.1
                             coeffs.2 <- - coeffs.2
                          })
                        if(top != 0) {
                          coeffs.1 <- sort(coeffs.1, decreasing=TRUE)
                          coeffs.2 <- sort(coeffs.2, decreasing=TRUE)
                          coeffs.1 <- coeffs.1[1:min(top, length(coeffs.1))]
                          coeffs.2 <- coeffs.2[1:min(top, length(coeffs.2))]
                        }
                        n.1 <- length(setdiff(names(coeffs.1), names(coeffs.2)))
                        n.2 <- length(setdiff(names(coeffs.2), names(coeffs.1)))
                        both <- intersect(names(coeffs.1), names(coeffs.2))
                        n.both <- length(both)
                        n.neither <- n.f - n.1 - n.2 - n.both
                        both <- paste(both, collapse=",")
                        ## pval <- phyper(q = n.both - 1, m = n.1, n = n.f - n.1, k = n.2, lower.tail = FALSE) 
                        mat <- cbind(c(n.both, n.1), c(n.2, n.neither))
                        fet <- fisher.test(mat, alternative = "greater")
                        pval <- fet$p.value
                        vec <- c(r$id[1], n.1, n.2, n.both, both, pval)
                        names(vec) <- c("id", "n.1", "n.2", "n.both", "both", "pval")
                        vec
                      })
        hyp <- merge(hyp, m, by = "id")
        ret.tbl <- rbind(ret.tbl, hyp)
      }
    }
  }
  ret.tbl
}

## Looking at correlation of coefficients
## method = c("pearson", "spearman", "kendall")
find.overlap.of.dense.predictors <- function(tbl.input, method) {
  ret.tbl <- c()
  for(gs in unique(tbl.input$gene.set)) {
    tmp <- subset(tbl.input, gene.set == gs)
    tr.sets <- as.character(unique(tmp$train.set))
    for(i in 1:(length(tr.sets)-1)) {
      tr1 <- tr.sets[i]
      tmp1 <- subset(tmp, train.set == tr1)
      merge.col <- drug.name.col
      tmp1 <- merge(drug.name.tbl, tmp1, by.x = merge.col, by.y = "train.drug")
      for(j in (i+1):length(tr.sets)) {
        tr2 <- tr.sets[j]
        if(((tr1 == "ohsu") || (tr2 == "ohsu")) && (grepl(pattern="ohsu.set1", paste0(tr1,tr2)) || grepl(pattern="ohsu.set2", paste0(tr1,tr2)))) { next }
        if(((tr1 == "fimm") || (tr2 == "fimm")) && (grepl(pattern="fimm.set1", paste0(tr1,tr2)) || grepl(pattern="fimm.set2", paste0(tr1,tr2)))) { next }
        tmp2 <- subset(tmp, train.set == tr2)
        merge.col <- drug.name.col
        tmp2 <- merge(drug.name.tbl, tmp2, by.x = merge.col, by.y = "train.drug")
  
        m <- merge(tmp1, tmp2, by = merge.col, suffixes = c(".1", ".2"))
        m <- m[, c(merge.col, "n.features.1", "coeffs.1", "coeffs.2", "coeff.vals.1", "coeff.vals.2", "train.set.1", "train.set.2", "gene.set.1")]
        m$id <- 1:nrow(m)
        hyp <- ddply(m, .variables = "id",
               .fun = function(r) {
                        coeff.vals.1 <- as.numeric(unlist(strsplit(as.character(r$coeff.vals.1), split=",")))
                        coeffs.1 <- (unlist(strsplit(as.character(r$coeffs.1), split=",")))
                        names(coeff.vals.1) <- coeffs.1
                        coeffs.1 <- coeff.vals.1
                        coeff.vals.2 <- as.numeric(unlist(strsplit(as.character(r$coeff.vals.2), split=",")))
                        coeffs.2 <- (unlist(strsplit(as.character(r$coeffs.2), split=",")))
                        names(coeff.vals.2) <- coeffs.2
                        coeffs.2 <- coeff.vals.2
                        inter <- intersect(names(coeffs.1), names(coeffs.2))
                        ct <- cor.test(coeffs.1[inter], coeffs.2[inter], method = method)
                        pval <- ct$p.value
                        cor <- ct$estimate
                        vec <- c(r$id[1], cor, pval)
                        names(vec) <- c("id", "cor", "pval")
                        vec
                      })
        hyp <- merge(hyp, m, by = "id")
        ret.tbl <- rbind(ret.tbl, hyp)
      }
    }
  }
  ret.tbl
}

model.fimm.and.ohsu <- function(ohsu.dss.arg, ohsu.expr.arg, ohsu.genomic.arg, ohsu.clinical.arg, ohsu.common.drugs, ohsu.drugs, fimm.dss.arg, fimm.expr.arg, fimm.genomic.arg, fimm.clinical.arg, fimm.common.drugs, fimm.drugs, response.col, seed = 1234, train.on.ohsu = TRUE, alphas = c(0, 1)) {

   l <- prepare.ohsu.drug.response.and.expr.matrices(ohsu.dss.arg, ohsu.expr.arg, rnaseq.summary.df = NULL, drugs = ohsu.common.drugs, response.col = response.col)
   ohsu.drc <- l[["drc.df"]]
   ohsu.expr <- l[["expr.df"]]

   l <- prepare.fimm.drug.response.and.expr.matrices(fimm.dss.arg, fimm.expr.arg, drugs = fimm.common.drugs, response.col = response.col)
   fimm.drc <- l[["drc.df"]]
   fimm.expr <- l[["expr.df"]]
 
   ## For the moment, just fit to a single drug
##   if(length(ohsu.drugs) > 1) { 
##     cat("Code only fits a single drug at the moment!\n")
##     stop("stop")
##   }

##   if(length(fimm.drugs) > 1) { 
##     cat("Code only fits a single drug at the moment!\n")
##     stop("stop")
##   }

   if(is.null(ohsu.expr) && !is.null(fimm.expr)) {
     stop("OHSU expr is null, but FIMM expr is not null!\n")
   }
   if(!is.null(ohsu.expr) && is.null(fimm.expr)) {
     stop("OHSU expr is not null, but FIMM expr is null!\n")
   }
   if(!is.null(ohsu.expr) && !is.null(fimm.expr)) {
     common.genes <- intersect(rownames(ohsu.expr), rownames(fimm.expr))
     ohsu.expr <- ohsu.expr[common.genes,]
     fimm.expr <- fimm.expr[common.genes,]

     ## Z-score the genes/rows (individually)
     z.score <- TRUE
     if(z.score) {
       ohsu.expr <- t(scale(t(ohsu.expr), center = TRUE, scale = TRUE))
       fimm.expr <- t(scale(t(fimm.expr), center = TRUE, scale = TRUE))
     }
   }

   ohsu.samples <- colnames(ohsu.drc)
   fimm.samples <- colnames(fimm.drc)

   ohsu.genomic <- ohsu.genomic.arg
   fimm.genomic <- fimm.genomic.arg
   if(is.null(ohsu.genomic) && !is.null(fimm.genomic)) {
     stop("OHSU genomic is null, but FIMM genomic is not null!\n")
   }
   if(!is.null(ohsu.genomic) && is.null(fimm.genomic)) {
     stop("OHSU genomic is not null, but FIMM genomic is null!\n")
   }

   if(!is.null(ohsu.genomic) && !is.null(fimm.genomic)) {
     common.genes <- intersect(rownames(ohsu.genomic), rownames(fimm.genomic))

     ohsu.genomic <- ohsu.genomic[common.genes,]
     fimm.genomic <- fimm.genomic[common.genes,]

     ohsu.samples <- intersect(ohsu.samples, colnames(ohsu.genomic))
     fimm.samples <- intersect(fimm.samples, colnames(fimm.genomic))
   }

   ohsu.clinical <- ohsu.clinical.arg
   fimm.clinical <- fimm.clinical.arg
   if(is.null(ohsu.clinical) && !is.null(fimm.clinical)) {
     stop("OHSU clinical is null, but FIMM clinical is not null!\n")
   }
   if(!is.null(ohsu.clinical) && is.null(fimm.clinical)) {
     stop("OHSU clinical is not null, but FIMM clinical is null!\n")
   }

   if(!is.null(ohsu.clinical) && !is.null(fimm.clinical)) {
     common.annotations <- intersect(rownames(ohsu.clinical), rownames(fimm.clinical))

     ohsu.clinical <- ohsu.clinical[common.annotations,]
     fimm.clinical <- fimm.clinical[common.annotations,]

     ohsu.samples <- intersect(ohsu.samples, colnames(ohsu.clinical))
     fimm.samples <- intersect(fimm.samples, colnames(fimm.clinical))
   }

   if(!is.null(ohsu.expr) && !is.null(fimm.expr)) {
     ohsu.samples <- intersect(ohsu.samples, colnames(ohsu.expr))
     fimm.samples <- intersect(fimm.samples, colnames(fimm.expr))
   }

   ohsu.drc <- ohsu.drc[, ohsu.samples]
   fimm.drc <- fimm.drc[, fimm.samples]

   if(!is.null(ohsu.genomic) && !is.null(fimm.genomic)) {
     ohsu.genomic <- ohsu.genomic[, ohsu.samples]
     fimm.genomic <- fimm.genomic[, fimm.samples]
   }
   if(!is.null(ohsu.clinical) && !is.null(fimm.clinical)) {
     ohsu.clinical <- ohsu.clinical[, ohsu.samples]
     fimm.clinical <- fimm.clinical[, fimm.samples]
   }
   if(!is.null(ohsu.expr) && !is.null(fimm.expr)) {
     ohsu.expr <- ohsu.expr[, ohsu.samples]
     fimm.expr <- fimm.expr[, fimm.samples]
   }

   ## Adjust the location (mean or median) and scale (mad or sd) of the two drug-response data sets to be the same--
   ## based on all drugs in common, even if we don't test them here.
##   fimm.location <- median(unlist(fimm.drc), na.rm=TRUE)
##   ohsu.location <- median(unlist(ohsu.drc), na.rm=TRUE)
##   fimm.scale <- mad(unlist(fimm.drc), na.rm=TRUE)
##   ohsu.scale <- mad(unlist(ohsu.drc), na.rm=TRUE)

   ## Adjust by mean and sd
   fimm.location <- mean(unlist(fimm.drc), na.rm=TRUE)
   ohsu.location <- mean(unlist(ohsu.drc), na.rm=TRUE)
   fimm.scale <- sd(unlist(fimm.drc), na.rm=TRUE)
   ohsu.scale <- sd(unlist(ohsu.drc), na.rm=TRUE)

   fimm.drc <- ( ohsu.scale * ( fimm.drc - fimm.location ) / fimm.scale ) + ohsu.location

   ## By default, train on OHSU data and validate on FIMM
   train.drc <- ohsu.drc
   train.drugs <- ohsu.drugs
   train.genomic <- ohsu.genomic
   train.clinical <- ohsu.clinical
   train.expr <- ohsu.expr
   test.drc <- fimm.drc
   test.drugs <- fimm.drugs
   test.genomic <- fimm.genomic
   test.clinical <- fimm.clinical
   test.expr <- fimm.expr

   if(!train.on.ohsu) {
     ## Instead, train on FIMM data and validate on OHSU
     test.drc <- ohsu.drc
     test.drugs <- ohsu.drugs
     test.genomic <- ohsu.genomic
     test.clinical <- ohsu.clinical
     test.expr <- ohsu.expr
     train.drc <- fimm.drc
     train.drugs <- fimm.drugs
     train.genomic <- fimm.genomic
     train.clinical <- fimm.clinical
     train.expr <- fimm.expr
   }
   rm(ohsu.drc)
   rm(ohsu.drugs)
   rm(ohsu.genomic)
   rm(ohsu.clinical)
   rm(ohsu.expr)
   rm(fimm.drc)
   rm(fimm.drugs)
   rm(fimm.genomic)
   rm(fimm.clinical)
   rm(fimm.expr)
   gc()

   ## Fit ridge to training data and apply to validation for each drug independently.
   ret <- llply(1:length(train.drugs),
                .parallel = TRUE,
                .fun = function(drug.i) {
                         train.drug <- train.drugs[drug.i]
                         test.drug <- test.drugs[drug.i]
                         cat(paste0("Modeling ", train.drug, " using ", response.col, "\n"))

                         x <- rbind(train.expr, train.genomic, train.clinical)
                         y <- as.numeric(train.drc[train.drug,])
                         flag <- is.na(y)
                         x <- x[, !flag]
                         y <- y[!flag]

                         newx <- rbind(test.expr, test.genomic, test.clinical)
                         newy <- as.numeric(test.drc[test.drug,])
                         flag <- is.na(newy)
                         newx <- newx[, !flag]
                         newy <- newy[!flag]

                         n.train <- length(y)
                         n.test <- length(newy)

                         ## Exclude any genes that have no variation in either data set (or were set to NA/NaN by
                         ## z-scoring above)
                         constant.genes <- rownames(x)[unlist(apply(x, 1, function(row) any(!is.finite(row)) | all(row == row[1])))]
                         constant.genes <- unique(c(constant.genes, rownames(newx)[unlist(apply(newx, 1, function(row) any(!is.finite(row)) | all(row == row[1])))]))
                         x <- x[!(rownames(x) %in% constant.genes),]
                         newx <- newx[!(rownames(newx) %in% constant.genes),]
##print(dim(x))
##print(dim(newx))
##print(head(x[,1]))
##print(head(x[1,]))
##print(head(newx[,1]))
##print(head(newx[1,]))


                         nfolds <- 5 
                         ret.list <- list()
                         ## Apply fit to test/validation data

                         ## Use glmnet
                         ## alphas <- c(0, 0.25, 0.5, 0.7, 1)
                         ## alphas <- c(0, 1)
                         N <- length(y)
                         foldid <- sample(rep(seq(nfolds), length = N))

                         svalues <- c("lambda.min", "lambda.1se")
                         svalues <- c("lambda.min")
                         for(alpha in alphas) {
                           for(s in svalues) {
                             set.seed(seed)
                             fit <- tryCatch({cv.glmnet(x = t(x), y, family = "gaussian", type.measure = "mse", nfolds = nfolds, foldid = foldid, alpha = alpha, standardize = FALSE, parallel = TRUE)}, error = function(e) { cat("glmnet err\n"); return(NA) })
                             pred <- tryCatch({predict(fit, newx=t(newx), s = s)}, error = function(e) { cat("glmnet pred err\n"); return(NA) })
                             ret.list[[length(ret.list)+1]] <- list(train.drug = train.drug, test.drug = test.drug, predicted = pred, actual = newy, fit = fit, alpha = alpha, s = s, model = "glmnet", n.train = n.train, n.test = n.test)
                           }
                         }

                         ## Use RF
                         set.seed(seed)
                         fit <- tryCatch({randomForest(x = t(x), y, importance = TRUE, keep.forest = TRUE)}, error = function(e) { cat("rf error\n"); return(NA) })
                         pred <- tryCatch({predict(fit, t(newx))}, error = function(e) { cat("rf pred error\n"); return(NA) })
                         ret.list[[length(ret.list)+1]] <- list(train.drug = train.drug, test.drug = test.drug, predicted = pred, actual = newy, fit = fit, alpha = NA, s = NA, model = "rf", n.train = n.train, n.test = n.test)

                         ## Use SVM
                         set.seed(seed)
                         fit <- tryCatch({svm(x = t(x), y, scale = FALSE)}, error = function(e) { cat("svm error\n"); return(NA) })
                         pred <- tryCatch({predict(fit, t(newx))}, error = function(e) { cat("svm pred error\n"); return(NA) })
                         ret.list[[length(ret.list)+1]] <- list(train.drug = train.drug, test.drug = test.drug, predicted = pred, actual = newy, fit = fit, alpha = NA, s = NA, model = "svm", n.train = n.train, n.test = n.test)
                         ret.list
                })
   ret
} 

model.train.and.test <- function(train.dss.arg, train.expr.arg, train.genomic.arg, train.clinical.arg, train.common.drugs, train.drugs, train.drug.col, train.patient.col, train.response.col, test.dss.arg, test.expr.arg, test.genomic.arg, test.clinical.arg, test.common.drugs, test.drugs, test.drug.col, test.patient.col, test.response.col, seed = 1234, use.rf = FALSE, use.svm = FALSE, use.ridge = TRUE, use.mean = TRUE, alphas = c(0, 1)) {

   l <- prepare.drug.response.and.expr.matrices(train.dss.arg, train.expr.arg, drugs = train.common.drugs, drug.col = train.drug.col, patient.col = train.patient.col, response.col = train.response.col)
   train.drc <- l[["drc.df"]]
   train.expr <- l[["expr.df"]]

   l <- prepare.drug.response.and.expr.matrices(test.dss.arg, test.expr.arg, drugs = test.common.drugs, drug.col = test.drug.col, patient.col = test.patient.col, response.col = test.response.col)
   test.drc <- l[["drc.df"]]
   test.expr <- l[["expr.df"]]
 
   if(is.null(train.expr) && !is.null(test.expr)) {
     stop("train expr is null, but FIMM expr is not null!\n")
   }
   if(!is.null(train.expr) && is.null(test.expr)) {
     stop("train expr is not null, but FIMM expr is null!\n")
   }
   if(!is.null(train.expr) && !is.null(test.expr)) {
     common.genes <- intersect(rownames(train.expr), rownames(test.expr))
     train.expr <- train.expr[common.genes,]
     test.expr <- test.expr[common.genes,]

     ## Z-score the genes/rows (individually)
     z.score <- TRUE
     if(z.score) {
       train.expr <- t(scale(t(train.expr), center = TRUE, scale = TRUE))
       test.expr <- t(scale(t(test.expr), center = TRUE, scale = TRUE))
     }
   }

   train.samples <- colnames(train.drc)
   test.samples <- colnames(test.drc)

   train.genomic <- train.genomic.arg
   test.genomic <- test.genomic.arg
   if(is.null(train.genomic) && !is.null(test.genomic)) {
     stop("train genomic is null, but FIMM genomic is not null!\n")
   }
   if(!is.null(train.genomic) && is.null(test.genomic)) {
     stop("train genomic is not null, but FIMM genomic is null!\n")
   }

   if(!is.null(train.genomic) && !is.null(test.genomic)) {
     common.genes <- intersect(rownames(train.genomic), rownames(test.genomic))

     train.genomic <- train.genomic[common.genes,]
     test.genomic <- test.genomic[common.genes,]

     train.samples <- intersect(train.samples, colnames(train.genomic))
     test.samples <- intersect(test.samples, colnames(test.genomic))
   }

   train.clinical <- train.clinical.arg
   test.clinical <- test.clinical.arg
   if(is.null(train.clinical) && !is.null(test.clinical)) {
     stop("train clinical is null, but FIMM clinical is not null!\n")
   }
   if(!is.null(train.clinical) && is.null(test.clinical)) {
     stop("train clinical is not null, but FIMM clinical is null!\n")
   }

   if(!is.null(train.clinical) && !is.null(test.clinical)) {
     common.annotations <- intersect(rownames(train.clinical), rownames(test.clinical))

     train.clinical <- train.clinical[common.annotations,]
     test.clinical <- test.clinical[common.annotations,]

     train.samples <- intersect(train.samples, colnames(train.clinical))
     test.samples <- intersect(test.samples, colnames(test.clinical))
   }

   if(!is.null(train.expr) && !is.null(test.expr)) {
     train.samples <- intersect(train.samples, colnames(train.expr))
     test.samples <- intersect(test.samples, colnames(test.expr))
   }

   train.drc <- train.drc[, train.samples]
   test.drc <- test.drc[, test.samples]

   if(!is.null(train.genomic) && !is.null(test.genomic)) {
     train.genomic <- train.genomic[, train.samples]
     test.genomic <- test.genomic[, test.samples]
   }
   if(!is.null(train.clinical) && !is.null(test.clinical)) {
     train.clinical <- train.clinical[, train.samples]
     test.clinical <- test.clinical[, test.samples]
   }
   if(!is.null(train.expr) && !is.null(test.expr)) {
     train.expr <- train.expr[, train.samples]
     test.expr <- test.expr[, test.samples]
   }

   ## Adjust the location (mean or median) and scale (mad or sd) of the two drug-response data sets to be the same--
   ## based on all drugs in common, even if we don't test them here.
##   test.location <- median(unlist(test.drc), na.rm=TRUE)
##   train.location <- median(unlist(train.drc), na.rm=TRUE)
##   test.scale <- mad(unlist(test.drc), na.rm=TRUE)
##   train.scale <- mad(unlist(train.drc), na.rm=TRUE)

   ## Adjust by mean and sd
   test.location <- mean(unlist(test.drc), na.rm=TRUE)
   train.location <- mean(unlist(train.drc), na.rm=TRUE)
   test.scale <- sd(unlist(test.drc), na.rm=TRUE)
   train.scale <- sd(unlist(train.drc), na.rm=TRUE)

   ## test.drc <- ( train.scale * ( test.drc - test.location ) / test.scale ) + train.location
   z.score.drc <- TRUE
   if(z.score.drc) {
       train.drc <- t(scale(t(train.drc), center = TRUE, scale = TRUE))
       test.drc <- t(scale(t(test.drc), center = TRUE, scale = TRUE))
   }

   ## Fit ridge to training data and apply to validation for each drug independently.
   ret <- llply(1:length(train.drugs),
                .parallel = TRUE,
                .fun = function(drug.i) {
                         train.drug <- train.drugs[drug.i]
                         test.drug <- test.drugs[drug.i]
                         cat(paste0("Modeling ", train.drug, " using ", train.response.col, "\n"))

                         x <- rbind(train.expr, train.genomic, train.clinical)
                         y <- as.numeric(train.drc[as.character(train.drug),])
                         flag <- is.na(y)
                         x <- x[, !flag]
                         y <- y[!flag]

                         newx <- rbind(test.expr, test.genomic, test.clinical)
                         newy <- as.numeric(test.drc[as.character(test.drug),])
                         flag <- is.na(newy)
                         newx <- newx[, !flag]
                         newy <- newy[!flag]

                         n.train <- length(y)
                         n.test <- length(newy)

                         ## Exclude any genes that have no variation in either data set (or were set to NA/NaN by
                         ## z-scoring above)
                         constant.genes <- rownames(x)[unlist(apply(x, 1, function(row) any(!is.finite(row)) | all(row == row[1])))]
                         constant.genes <- unique(c(constant.genes, rownames(newx)[unlist(apply(newx, 1, function(row) any(!is.finite(row)) | all(row == row[1])))]))
                         x <- x[!(rownames(x) %in% constant.genes),]
                         newx <- newx[!(rownames(newx) %in% constant.genes),]

                         nfolds <- 5 
                         ret.list <- list()
                         ## Apply fit to test/validation data
                         N <- length(y)
                         foldid <- sample(rep(seq(nfolds), length = N))

                         ## Use glmnet
                         ## alphas <- c(0, 0.25, 0.5, 0.7, 1)
                         ## alphas <- c(0, 1)
                         if(!use.ridge) { alphas <- alphas[alphas != 0] }
                         svalues <- c("lambda.min", "lambda.1se")
                         svalues <- c("lambda.min")
                         for(alpha in alphas) {
                           for(s in svalues) {
                             set.seed(seed)
                             fit <- tryCatch({cv.glmnet(x = t(x), y, family = "gaussian", type.measure = "mse", nfolds = nfolds, foldid = foldid, alpha = alpha, standardize = FALSE, parallel = TRUE)}, error = function(e) { cat("glmnet err\n"); return(NA) })
                             pred <- tryCatch({predict(fit, newx=t(newx), s = s)}, error = function(e) { cat("glmnet pred err\n"); return(NA) })
                             ret.list[[length(ret.list)+1]] <- list(train.drug = train.drug, test.drug = test.drug, predicted = pred, actual = newy, fit = fit, alpha = alpha, s = s, model = "glmnet", n.train = n.train, n.test = n.test)
                           }
                         }

                         ## Use RF
                         if(use.rf) {
                           set.seed(seed)
                           fit <- tryCatch({randomForest(x = t(x), y, importance = TRUE, keep.forest = TRUE)}, error = function(e) { cat("rf error\n"); return(NA) })
                           pred <- tryCatch({predict(fit, t(newx))}, error = function(e) { cat("rf pred error\n"); return(NA) })
                           ret.list[[length(ret.list)+1]] <- list(train.drug = train.drug, test.drug = test.drug, predicted = pred, actual = newy, fit = fit, alpha = NA, s = NA, model = "rf", n.train = n.train, n.test = n.test)
                         }

                         ## Use SVM
                         if(use.svm) {
                           set.seed(seed)
                           fit <- tryCatch({svm(x = t(x), y, scale = FALSE)}, error = function(e) { cat("svm error\n"); return(NA) })
                           pred <- tryCatch({predict(fit, t(newx))}, error = function(e) { cat("svm pred error\n"); return(NA) })
                           ret.list[[length(ret.list)+1]] <- list(train.drug = train.drug, test.drug = test.drug, predicted = pred, actual = newy, fit = fit, alpha = NA, s = NA, model = "svm", n.train = n.train, n.test = n.test)
                         }

                         ## Return a baseline in which the prediction is just the mean of training data set
                         if(!use.mean) { 
                           len <- length(newy) 
                           pred <- rep(mean(y, na.rm = TRUE), len)
                           ret.list[[length(ret.list)+1]] <- list(train.drug = train.drug, test.drug = test.drug, predicted = pred, actual = newy, fit = NA, alpha = NA, s = NA, model = "mean", n.train = n.train, n.test = n.test)
                         }

                         ret.list
                })
   ret
} 

## fits is as returned from train.model
## drug.name.tbl provides a mapping between the drug name in the training set (train.drug.col) and in the test set (test.drug.col)
test.model <- function(fits, drug.name.tbl, train.drug.col, test.dss.arg, test.expr.arg, test.genomic.arg, test.clinical.arg, test.common.drugs, test.drugs, test.drug.col, test.patient.col, test.response.col, seed = 1234, use.rf = FALSE, use.svm = FALSE, flatten.downsampled.fits = FALSE) {

   l <- prepare.drug.response.and.expr.matrices(test.dss.arg, test.expr.arg, drugs = test.common.drugs, drug.col = test.drug.col, patient.col = test.patient.col, response.col = test.response.col)
   test.drc <- l[["drc.df"]]
   test.expr <- l[["expr.df"]]

   if(!is.null(test.expr)) {
     ## Z-score the genes/rows (individually)
     z.score <- TRUE
     if(z.score) {
       test.expr <- t(scale(t(test.expr), center = TRUE, scale = TRUE))
     }
   }

   test.samples <- colnames(test.drc)

   test.genomic <- test.genomic.arg

   if(!is.null(test.genomic)) {
     test.samples <- intersect(test.samples, colnames(test.genomic))
   }

   test.clinical <- test.clinical.arg

   if(!is.null(test.clinical)) {
     test.samples <- intersect(test.samples, colnames(test.clinical))
   }

   if(!is.null(test.expr)) {
     test.samples <- intersect(test.samples, colnames(test.expr))
   }

   test.drc <- test.drc[, test.samples]

   if(!is.null(test.genomic)) {
     test.genomic <- test.genomic[, test.samples, drop = FALSE]
   }
   if(!is.null(test.clinical)) {
     test.clinical <- test.clinical[, test.samples, drop = FALSE]
   }
   if(!is.null(test.expr)) {
     test.expr <- test.expr[, test.samples, drop = FALSE]
   }

   z.score.drc <- TRUE
   if(z.score.drc) {
       test.drc <- t(scale(t(test.drc), center = TRUE, scale = TRUE))
   }

   ## Iterate through the fits (on training data) and apply them to the test data.
   ## fits is a list (of lists).  The ith entry is a list of fits for drug i.
   ## Each entry of _that_ list is a fit using a different method.
   ret <- llply(fits,
                .parallel = FALSE,
                .fun = function(drug.fits) {

                         newx.orig <- rbind(test.expr, test.genomic, test.clinical)

                         ret.list <- list()
                         if(flatten.downsampled.fits == TRUE) { drug.fits <- unlist(drug.fits, recursive = FALSE) }
                         for(i in 1:length(drug.fits)) {
                             train.drug <- drug.fits[[i]]$train.drug
                             test.drug <- drug.name.tbl[drug.name.tbl[, train.drug.col] == train.drug, test.drug.col]

                             newy <- as.numeric(test.drc[as.character(test.drug),])
                             flag <- is.na(newy)
                             newx <- newx.orig[, !flag]
                             newy <- newy[!flag]
                             n.test <- length(newy)

                             if(length(test.drug) != 1) { next }
                             fit <- drug.fits[[i]]$fit
                             if(is.na(fit)) { next }
                             alpha <- drug.fits[[i]]$alpha
                             s <- NA
                             model <- drug.fits[[i]]$model

                             if(model == "glmnet") {
                               flag <- rownames(newx) %in% rownames(coefficients(fit))
                               newx <- newx[flag, ]
                               svalues <- c("lambda.min")

                               for(s in svalues) {
                                 set.seed(seed)
                                 pred <- tryCatch({predict(fit, newx=t(newx), s = s)}, error = function(e) { cat("glmnet pred err\n"); return(NA) })
                                 ret.list[[length(ret.list)+1]] <- list(train.drug = train.drug, test.drug = test.drug, alpha = alpha, s = s, model = model, predicted = pred, actual = newy, n.test = n.test)
                               }
                             } else if(model == "rf") {
                               if(use.rf) {
                                 pred <- tryCatch({predict(fit, t(newx))}, error = function(e) { cat("glmnet pred err\n"); return(NA) })
                                 ret.list[[length(ret.list)+1]] <- list(train.drug = train.drug, test.drug = test.drug, alpha = alpha, s = s, model = model, predicted = pred, actual = newy, n.test = n.test)
                               }
                             } else if(model == "svm") {
                               if(use.svm) {
                                 pred <- tryCatch({predict(fit, newx=t(newx))}, error = function(e) { cat("glmnet pred err\n"); return(NA) })
                                 ret.list[[length(ret.list)+1]] <- list(train.drug = train.drug, test.drug = test.drug, alpha = alpha, s = s, model = model, predicted = pred, actual = newy, n.test = n.test)
                               }
                             } else if(model == "mean") {
                               pred <- rep(fit, length(newy))
                               ret.list[[length(ret.list)+1]] <- list(train.drug = train.drug, test.drug = test.drug, alpha = alpha, s = s, model = model, predicted = pred, actual = newy, n.test = n.test)
                             } else {
                               stop(paste0("Unknown model: ", model, "\n"))
                             }
                         }
                         ret.list
                })
   ret
} 

test.model.with.downsampling_ <- function(fits, drug.name.tbl, train.drug.col, test.expr, test.genomic, test.clinical, test.drc, test.drug.col, seed = 1234, use.rf = FALSE, use.svm = FALSE, use.ridge = TRUE, use.mean = TRUE, flatten.downsampled.fits = FALSE, num.samples.per.drug = NULL, with.replacement = FALSE, include.mean.response.as.feature = FALSE, subtract.mean.response = FALSE) {
   z.score.drc <- TRUE
   z.score <- TRUE
   ## Iterate through the fits (on training data) and apply them to the test data.
   ## fits is a list (of lists).  The ith entry is a list of fits for drug i.
   ## Each entry of _that_ list is a fit using a different method.
   ret <- llply(fits,
                .parallel = TRUE,
                .fun = function(drug.fits) {

                         ret.list <- list()
                         if(flatten.downsampled.fits == TRUE) { drug.fits <- unlist(drug.fits, recursive = FALSE) }
                         for(i in 1:length(drug.fits)) {
                             train.drug <- drug.fits[[i]]$train.drug
                             test.drug <- drug.name.tbl[drug.name.tbl[, train.drug.col] == train.drug, test.drug.col]
                             cat(paste0("Testing drug ", test.drug, "\n"))
                             y <- as.numeric(test.drc[as.character(test.drug),])
                             flag <- is.na(y)
                             y <- y[!flag]

                             t.expr <- test.expr
                             t.genomic <- test.genomic
                             t.clinical <- test.clinical
                             if(!is.null(t.expr)) { t.expr <- t.expr[, !flag, drop = FALSE] }
                             if(!is.null(t.genomic)) { t.genomic <- t.genomic[, !flag, drop = FALSE] }
                             if(!is.null(t.clinical)) { t.clinical <- t.clinical[, !flag, drop = FALSE] }
                             t.drc <- test.drc[, !flag]

                             if(!is.null(num.samples.per.drug)) {
                               set.seed(i)
                               n.samples <- num.samples.per.drug[num.samples.per.drug[, train.drug.col] == as.character(train.drug), "num.samples"]
                               cat(paste0("Testing ", test.drug, " using num samples = ", n.samples, "\n"))
                               if(n.samples > length(y)) {
                                 stop("Attempting to test ", test.drug, " using sampling ", ifelse(with.replacement, "with", "without"), " replacement of ", n.samples, " samples, but length of response is only ", length(y), "\n") 
                               }
                               indices <- sample.int(length(y), size = n.samples, replace = with.replacement)
                               y <- y[indices]
                               if(!is.null(t.expr)) { t.expr <- t.expr[, indices, drop = FALSE] }
                               if(!is.null(t.genomic)) { t.genomic <- t.genomic[, indices, drop = FALSE] }
                               if(!is.null(t.clinical)) { t.clinical <- t.clinical[, indices, drop = FALSE] }

                               t.drc <- t.drc[, indices]
                             }

                             ## Subtract off the mean response _before_ z-scoring
                             if(subtract.mean.response) {
                               t.drc.mean <- unname(colMeans(t.drc, na.rm=TRUE))
                               t.drc <- t.drc - matrix(t.drc.mean, nrow=nrow(t.drc), ncol=ncol(t.drc), byrow=TRUE)
                             }

                             if(z.score.drc) {
                               y <- (y - mean(y)) / sd(y)
                               t.drc <- t(scale(t(t.drc), center = TRUE, scale = TRUE))
                             }

                             if(z.score) {
                               t.expr <- t(scale(t(t.expr), center = TRUE, scale = TRUE))
                             }

                             ## Include mean _z-scored_ response
                             if(include.mean.response.as.feature && !is.null(t.expr)) {
                               y.mean <- colMeans(t.drc, na.rm=TRUE)
                               t.expr <- rbind(t.expr, "mean.resp" = unname(y.mean[colnames(t.expr)]))
                             }

                             newx <- rbind(t.expr, t.genomic, t.clinical)
                             newy <- y
                             n.test <- length(newy)

                             if(length(test.drug) != 1) { next }
                             fit <- drug.fits[[i]]$fit
                             if(is.na(fit)) { next }
                             alpha <- drug.fits[[i]]$alpha
                             s <- NA
                             model <- drug.fits[[i]]$model

                             if(model == "glmnet") {
                               flag <- rownames(newx) %in% rownames(coefficients(fit))
                               newx <- newx[flag, ]
                               svalues <- c("lambda.min")

                               if((alpha != 0) || (use.ridge == TRUE)) {
                                 for(s in svalues) {
                                   set.seed(seed)
                                   pred <- tryCatch({predict(fit, newx=t(newx), s = s)}, error = function(e) { cat("glmnet pred err\n"); return(NA) })
                                   ret.list[[length(ret.list)+1]] <- list(train.drug = train.drug, test.drug = test.drug, alpha = alpha, s = s, model = model, predicted = pred, actual = newy, n.test = n.test)
                                 }
                               }
                             } else if(model == "rf") {
                               if(use.rf) {
                                 pred <- tryCatch({predict(fit, t(newx))}, error = function(e) { cat("glmnet pred err\n"); return(NA) })
                                 ret.list[[length(ret.list)+1]] <- list(train.drug = train.drug, test.drug = test.drug, alpha = alpha, s = s, model = model, predicted = pred, actual = newy, n.test = n.test)
                               }
                             } else if(model == "svm") {
                               if(use.svm) {
                                 pred <- tryCatch({predict(fit, newx=t(newx))}, error = function(e) { cat("glmnet pred err\n"); return(NA) })
                                 ret.list[[length(ret.list)+1]] <- list(train.drug = train.drug, test.drug = test.drug, alpha = alpha, s = s, model = model, predicted = pred, actual = newy, n.test = n.test)
                               }
                             } else if(model == "mean") {
                               if(use.mean) {
                                 pred <- rep(fit, length(newy))
                                 ret.list[[length(ret.list)+1]] <- list(train.drug = train.drug, test.drug = test.drug, alpha = alpha, s = s, model = model, predicted = pred, actual = newy, n.test = n.test)
                               }
                             } else {
                               stop(paste0("Unknown model: ", model, "\n"))
                             }
                         }
                         ret.list
                })
   ret
}


test.model_ <- function(fits, drug.name.tbl, train.drug.col, x.arg, test.drc, test.drug.col, seed = 1234, use.rf = FALSE, use.svm = FALSE, use.ridge = TRUE, use.mean = TRUE, flatten.downsampled.fits = FALSE) {
   ## Iterate through the fits (on training data) and apply them to the test data.
   ## fits is a list (of lists).  The ith entry is a list of fits for drug i.
   ## Each entry of _that_ list is a fit using a different method.
   ret <- llply(fits,
                .parallel = FALSE,
                .fun = function(drug.fits) {

                         newx.orig <- x.arg

                         ret.list <- list()
                         if(flatten.downsampled.fits == TRUE) { drug.fits <- unlist(drug.fits, recursive = FALSE) }
                         for(i in 1:length(drug.fits)) {
                             train.drug <- drug.fits[[i]]$train.drug
                             test.drug <- drug.name.tbl[drug.name.tbl[, train.drug.col] == train.drug, test.drug.col]
                             newy <- as.numeric(test.drc[as.character(test.drug),])
                             flag <- is.na(newy)
                             newx <- newx.orig[, !flag, drop = F]
                             newy <- newy[!flag]

                             n.test <- length(newy)

                             if(length(test.drug) != 1) { next }
                             fit <- drug.fits[[i]]$fit
                             if(is.na(fit)) { next }
                             alpha <- drug.fits[[i]]$alpha
                             s <- NA
                             model <- drug.fits[[i]]$model

                             cat(paste0("Testing drug ", test.drug, " (n = ", n.test, ")\n"))

                             if(model == "glmnet") {
                               flag <- rownames(newx) %in% rownames(coefficients(fit))
                               newx <- newx[flag, , drop = F]
                               svalues <- c("lambda.min")

                               if((alpha != 0) || (use.ridge == TRUE)) {
                                 for(s in svalues) {
                                   set.seed(seed)
                                   pred <- tryCatch({predict(fit, newx=t(newx), s = s)}, error = function(e) { cat("glmnet pred err\n"); return(NA) })
                                   ret.list[[length(ret.list)+1]] <- list(train.drug = train.drug, test.drug = test.drug, alpha = alpha, s = s, model = model, predicted = pred, actual = newy, n.test = n.test)
                                 }
                               }
                             } else if(model == "rf") {
                               if(use.rf) {
                                 pred <- tryCatch({predict(fit, t(newx))}, error = function(e) { cat("glmnet pred err\n"); return(NA) })
                                 ret.list[[length(ret.list)+1]] <- list(train.drug = train.drug, test.drug = test.drug, alpha = alpha, s = s, model = model, predicted = pred, actual = newy, n.test = n.test)
                               }
                             } else if(model == "svm") {
                               if(use.svm) {
                                 pred <- tryCatch({predict(fit, newx=t(newx))}, error = function(e) { cat("glmnet pred err\n"); return(NA) })
                                 ret.list[[length(ret.list)+1]] <- list(train.drug = train.drug, test.drug = test.drug, alpha = alpha, s = s, model = model, predicted = pred, actual = newy, n.test = n.test)
                               }
                             } else if(model == "mean") {
                               if(use.mean) {
                                 pred <- rep(fit, length(newy))
                                 ret.list[[length(ret.list)+1]] <- list(train.drug = train.drug, test.drug = test.drug, alpha = alpha, s = s, model = model, predicted = pred, actual = newy, n.test = n.test)
                               }
                             } else {
                               stop(paste0("Unknown model: ", model, "\n"))
                             }
                         }
                         ret.list
                })
   ret
} 


## This function does not format the drug or expression data
train.model_ <- function(x.arg, train.drc, train.drugs, seed = 1234, use.rf = FALSE, use.svm = FALSE, use.ridge = TRUE, use.mean = TRUE, calc.coefficient.pvals = FALSE, alphas = c(0, 1), keep.forest = TRUE, use.glmnet.intercept = TRUE) {

   ## Fit ridge to training data and apply to validation for each drug independently.
   ret <- llply(1:length(train.drugs),
                .parallel = TRUE,
                .fun = function(drug.i) {
                         train.drug <- train.drugs[drug.i]
                         x <- x.arg
                         y <- as.numeric(train.drc[as.character(train.drug),])
                         flag <- is.na(y)
                         x <- x[, !flag, drop=F]
                         y <- y[!flag]
                         n.train <- length(y)
                         cat(paste0("Modeling ", train.drug, " (n = ", n.train, ")\n"))

                         ## Exclude any genes that have no variation in training data set (or were set to NA/NaN by
                         ## z-scoring above)
                         constant.genes <- rownames(x)[unlist(apply(x, 1, function(row) any(!is.finite(row))))]
##                         constant.genes <- rownames(x)[unlist(apply(x, 1, function(row) any(!is.finite(row)) || all(row == row[1])))]
##                         x <- x[!(rownames(x) %in% constant.genes), , drop=F]

                         n.features <- nrow(x)

                         nfolds <- 5 
                         ret.list <- list()
                         ## Apply fit to test/validation data
                         N <- length(y)
                         set.seed(seed)
                         foldid <- sample(rep(seq(nfolds), length = N))

                         ## Use glmnet
                         ## alphas <- c(0, 0.25, 0.5, 0.7, 1)
                         ## alphas <- c(0, 1)
                         if(!use.ridge) { alphas <- alphas[alphas != 0] }
                         for(alpha in alphas) {
                           pvals <- NA
                           if( (alpha == 1) && (calc.coefficient.pvals == TRUE) ) {
                             ## lasso
                             ## Calculate the p-values for the individual components
                             fit.lasso <- lasso.proj(t(x), y)
                             pvals <- fit.lasso$pval
                           }
                           set.seed(seed)
                           fit <- tryCatch({cv.glmnet(x = t(x), y, family = "gaussian", type.measure = "mse", nfolds = nfolds, foldid = foldid, alpha = alpha, standardize = FALSE, parallel = TRUE, intercept = use.glmnet.intercept)}, error = function(e) { print(e); cat("glmnet err\n"); return(NA) })
                           ret.list[[length(ret.list)+1]] <- list(train.drug = train.drug, fit = fit, alpha = alpha, model = "glmnet", n.train = n.train, pvals = pvals, n.features = n.features)
                         }

                         ## Use RF
                         if(use.rf) {
                           set.seed(seed)
                           fit <- tryCatch({randomForest(x = t(x), y, importance = TRUE, keep.forest = keep.forest)}, error = function(e) { cat("rf error\n"); return(NA) })
                           ret.list[[length(ret.list)+1]] <- list(train.drug = train.drug, fit = fit, alpha = NA, model = "rf", n.train = n.train, pvals = NA, n.features = n.features)
                         }

                         ## Use SVM
                         if(use.svm) {
                           set.seed(seed)
                           fit <- tryCatch({svm(x = t(x), y, scale = FALSE)}, error = function(e) { cat("svm error\n"); return(NA) })
                           ret.list[[length(ret.list)+1]] <- list(train.drug = train.drug, fit = fit, alpha = NA, model = "svm", n.train = n.train, pvals = NA, n.features = n.features)
                         }

                         ## Return a baseline in which the prediction is just the mean of training data set
                         if(use.mean) { 
                           pred <- mean(y, na.rm = TRUE)
                           ret.list[[length(ret.list)+1]] <- list(train.drug = train.drug, fit = pred, alpha = NA, model = "mean", n.train = n.train, pvals = NA, n.features = n.features)
                         }
                         ret.list
                })
   ret
} 

train.model.with.downsampling_ <- function(train.expr, train.genomic, train.clinical, train.drc, train.drugs, train.drug.col, train.response.col, seed = 1234, use.rf = FALSE, use.svm = FALSE, use.ridge = TRUE, use.mean = TRUE, calc.coefficient.pvals = FALSE, num.samples.per.drug, num.iterations = 100, with.replacement = FALSE, alphas = c(0, 1), include.mean.response.as.feature = FALSE, subtract.mean.response = FALSE, keep.forest = TRUE) {
   z.score.drc <- TRUE
   z.score <- TRUE
   ret <- llply(1:length(train.drugs),
                .parallel = FALSE,
                .fun = function(drug.i) {
                         llply(1:num.iterations, 
                               .parallel = FALSE,
                               .fun = function(iter) {
                                        train.drug <- train.drugs[drug.i]

                                        ## Z-score training _after_ sub-sampling
                                        y <- as.numeric(train.drc[as.character(train.drug),])
                                        flag <- is.na(y)
                                        y <- y[!flag]
                                        t.drc <- train.drc[, !flag]
                                        t.expr <- train.expr
                                        t.genomic <- train.genomic
                                        t.clinical <- train.clinical
                                        if(!is.null(t.expr)) { t.expr <- t.expr[, !flag, drop = FALSE] }
                                        if(!is.null(t.genomic)) { t.genomic <- t.genomic[, !flag, drop = FALSE] }
                                        if(!is.null(t.clinical)) { t.clinical <- t.clinical[, !flag, drop = FALSE] }

                                        set.seed(iter)

                                        n.samples <- num.samples.per.drug[num.samples.per.drug[, train.drug.col] == as.character(train.drug), "num.samples"]
                                        cat(paste0("Modeling ", train.drug, " using ", train.response.col, " and num samples = ", n.samples, "\n"))
                                        if(n.samples > length(y)) {
                                          stop("Attempting to model ", train.drug, " using sampling ", ifelse(with.replacement, "with", "without"), " replacement of ", n.samples, " samples, but length of response is only ", length(y), "\n") 
                                        }
                                        indices <- sample.int(length(y), size = n.samples, replace = with.replacement)
                                        y <- y[indices]
                                        if(!is.null(t.expr)) { t.expr <- t.expr[, indices, drop = FALSE] }
                                        if(!is.null(t.genomic)) { t.genomic <- t.genomic[, indices, drop = FALSE] }
                                        if(!is.null(t.clinical)) { t.clinical <- t.clinical[, indices, drop = FALSE] }

                                        t.drc <- t.drc[, indices]

                                        ## Subtract off the mean response _before_ z-scoring
                                        if(subtract.mean.response) {
                                          t.drc.mean <- unname(colMeans(t.drc, na.rm=TRUE))
                                          t.drc <- t.drc - matrix(t.drc.mean, nrow=nrow(t.drc), ncol=ncol(t.drc), byrow=TRUE)
                                        }

                                        if(z.score.drc) {
                                          y <- (y - mean(y)) / sd(y)
                                          t.drc <- t(scale(t(t.drc), center = TRUE, scale = TRUE))
                                        }

                                        if(z.score) {
                                          t.expr <- t(scale(t(t.expr), center = TRUE, scale = TRUE))
                                        }

                                        ## Include mean _z-scored_ response
                                        if(include.mean.response.as.feature && !is.null(t.expr)) {
                                          y.mean <- colMeans(t.drc, na.rm=TRUE)
                                          t.expr <- rbind(t.expr, "mean.resp" = unname(y.mean[colnames(t.expr)]))
                                        }

                                        x <- rbind(t.expr, t.genomic, t.clinical)
                                        n.train <- length(y)

                                        ## Exclude any genes that have no variation in training data set (or were set to NA/NaN by
                                        ## z-scoring above)
                                        constant.genes <- rownames(x)[unlist(apply(x, 1, function(row) any(!is.finite(row)) | all(row == row[1])))]
                                        x <- x[!(rownames(x) %in% constant.genes),]

                                        nfolds <- 5 
                                        ret.list <- list()
                                        ## Apply fit to test/validation data
                                        N <- length(y)
                           set.seed(seed)
                                        foldid <- sample(rep(seq(nfolds), length = N))

                                        ## Use glmnet
                                        ## alphas <- c(0, 0.25, 0.5, 0.7, 1)
                                        ## alphas <- c(0, 1)
                                        for(alpha in alphas) {
                                          pvals <- NA
                                          if( (alpha == 1) && (calc.coefficient.pvals == TRUE) ) {
                                            ## lasso
                                            ## Calculate the p-values for the individual components
                                            fit.lasso <- lasso.proj(t(x), y)
                                            pvals <- fit.lasso$pval
                                          }
                                          set.seed(seed)
                                          fit <- tryCatch({cv.glmnet(x = t(x), y, family = "gaussian", type.measure = "mse", nfolds = nfolds, foldid = foldid, alpha = alpha, standardize = FALSE, parallel = TRUE)}, error = function(e) { cat("glmnet err\n"); return(NA) })
                                          ret.list[[length(ret.list)+1]] <- list(train.drug = train.drug, fit = fit, alpha = alpha, model = "glmnet", n.train = n.train, pvals = pvals, iter = iter)
                                        }

                                        ## Use RF
                                        if(use.rf) {
                                          set.seed(seed)
                                          fit <- tryCatch({randomForest(x = t(x), y, importance = TRUE, keep.forest = keep.forest)}, error = function(e) { cat("rf error\n"); return(NA) })
                                          ret.list[[length(ret.list)+1]] <- list(train.drug = train.drug, fit = fit, alpha = NA, model = "rf", n.train = n.train, pvals = NA, iter = iter)
                                        }

                                        ## Use SVM
                                        if(use.svm) {
                                          set.seed(seed)
                                          fit <- tryCatch({svm(x = t(x), y, scale = FALSE)}, error = function(e) { cat("svm error\n"); return(NA) })
                                          ret.list[[length(ret.list)+1]] <- list(train.drug = train.drug, fit = fit, alpha = NA, model = "svm", n.train = n.train, pvals = NA, iter = iter)
                                        }

                                        ## Return a baseline in which the prediction is just the mean of training data set
                                        if(use.mean) {
                                          pred <- mean(y, na.rm = TRUE)
                                          ret.list[[length(ret.list)+1]] <- list(train.drug = train.drug, fit = pred, alpha = NA, model = "mean", n.train = n.train, pvals = NA, iter = iter)
                                        }
                                        ret.list
                                      })  
                })
   ret
}

train.model <- function(train.dss.arg, train.expr.arg, train.genomic.arg, train.clinical.arg, train.common.drugs, train.drugs, train.drug.col, train.patient.col, train.response.col, seed = 1234, use.rf = FALSE, use.svm = FALSE, calc.coefficient.pvals = FALSE, alphas = c(0, 1)) {

   l <- prepare.drug.response.and.expr.matrices(train.dss.arg, train.expr.arg, drugs = train.common.drugs, drug.col = train.drug.col, patient.col = train.patient.col, response.col = train.response.col)
   train.drc <- l[["drc.df"]]
   train.expr <- l[["expr.df"]]

   if(!is.null(train.expr)) {
     ## Z-score the genes/rows (individually)
     z.score <- TRUE
     if(z.score) {
       train.expr <- t(scale(t(train.expr), center = TRUE, scale = TRUE))
     }
   }

   train.samples <- colnames(train.drc)

   train.genomic <- train.genomic.arg

   if(!is.null(train.genomic)) {
     train.samples <- intersect(train.samples, colnames(train.genomic))
   }

   train.clinical <- train.clinical.arg

   if(!is.null(train.clinical)) {
     train.samples <- intersect(train.samples, colnames(train.clinical))
   }

   if(!is.null(train.expr)) {
     train.samples <- intersect(train.samples, colnames(train.expr))
   }

   train.drc <- train.drc[, train.samples]

   if(!is.null(train.genomic)) {
     train.genomic <- train.genomic[, train.samples, drop = FALSE]
   }
   if(!is.null(train.clinical)) {
     train.clinical <- train.clinical[, train.samples, drop = FALSE]
   }
   if(!is.null(train.expr)) {
     train.expr <- train.expr[, train.samples, drop = FALSE]
   }

   z.score.drc <- TRUE
   if(z.score.drc) {
       train.drc <- t(scale(t(train.drc), center = TRUE, scale = TRUE))
   }

   ## Fit ridge to training data and apply to validation for each drug independently.
   ret <- llply(1:length(train.drugs),
                .parallel = TRUE,
                .fun = function(drug.i) {
                         train.drug <- train.drugs[drug.i]
                         cat(paste0("Modeling ", train.drug, " using ", train.response.col, "\n"))

                         x <- rbind(train.expr, train.genomic, train.clinical)
                         y <- as.numeric(train.drc[as.character(train.drug),])
                         flag <- is.na(y)
                         x <- x[, !flag]
                         y <- y[!flag]

                         n.train <- length(y)

                         ## Exclude any genes that have no variation in training data set (or were set to NA/NaN by
                         ## z-scoring above)
                         constant.genes <- rownames(x)[unlist(apply(x, 1, function(row) any(!is.finite(row)) | all(row == row[1])))]
                         x <- x[!(rownames(x) %in% constant.genes),]

                         nfolds <- 5 
                         ret.list <- list()
                         ## Apply fit to test/validation data
                         N <- length(y)
                           set.seed(seed)
                         foldid <- sample(rep(seq(nfolds), length = N))

                         ## Use glmnet
                         ## alphas <- c(0, 0.25, 0.5, 0.7, 1)
                         ## alphas <- c(0, 1)
                         for(alpha in alphas) {
                           pvals <- NA
                           if( (alpha == 1) && (calc.coefficient.pvals == TRUE) ) {
                             ## lasso
                             ## Calculate the p-values for the individual components
                             fit.lasso <- lasso.proj(t(x), y)
                             pvals <- fit.lasso$pval
                           }
                           set.seed(seed)
                           fit <- tryCatch({cv.glmnet(x = t(x), y, family = "gaussian", type.measure = "mse", nfolds = nfolds, foldid = foldid, alpha = alpha, standardize = FALSE, parallel = TRUE)}, error = function(e) { cat("glmnet err\n"); return(NA) })
                           ret.list[[length(ret.list)+1]] <- list(train.drug = train.drug, fit = fit, alpha = alpha, model = "glmnet", n.train = n.train, pvals = pvals)
                         }

                         ## Use RF
                         if(use.rf) {
                           set.seed(seed)
                           fit <- tryCatch({randomForest(x = t(x), y, importance = TRUE, keep.forest = TRUE)}, error = function(e) { cat("rf error\n"); return(NA) })
                           ret.list[[length(ret.list)+1]] <- list(train.drug = train.drug, fit = fit, alpha = NA, model = "rf", n.train = n.train, pvals = NA)
                         }

                         ## Use SVM
                         if(use.svm) {
                           set.seed(seed)
                           fit <- tryCatch({svm(x = t(x), y, scale = FALSE)}, error = function(e) { cat("svm error\n"); return(NA) })
                           ret.list[[length(ret.list)+1]] <- list(train.drug = train.drug, fit = fit, alpha = NA, model = "svm", n.train = n.train, pvals = NA)
                         }

                         ## Return a baseline in which the prediction is just the mean of training data set
                         pred <- mean(y, na.rm = TRUE)
                         ret.list[[length(ret.list)+1]] <- list(train.drug = train.drug, fit = pred, alpha = NA, model = "mean", n.train = n.train, pvals = NA)

                         ret.list
                })
   ret
} 


fit.elastic.net.alpha <- function(alpha, x, y, foldid, nfolds, seed = 1234, standardize = TRUE) {
  N <- nrow(x)
  cv <- cv.glmnet(x, y, family = "gaussian", type.measure = "mse", foldid = foldid, nfolds = nfolds, alpha = alpha, standardize = standardize, parallel = FALSE)
  cv
}



bootstrap.model <- function(train.dss.arg, train.expr.arg, train.common.drugs, train.drugs, train.drug.col, train.patient.col, train.response.col, seed = 1234, num.bootstraps = 500, alphas = c(0, 1)) {

   l <- prepare.drug.response.and.expr.matrices(train.dss.arg, train.expr.arg, drugs = train.common.drugs, drug.col = train.drug.col, patient.col = train.patient.col, response.col = train.response.col)
   train.drc <- l[["drc.df"]]
   train.expr <- l[["expr.df"]]

   train.samples <- colnames(train.drc)

   if(!is.null(train.expr)) {
     train.samples <- intersect(train.samples, colnames(train.expr))
   }

   train.drc <- train.drc[, train.samples]

   if(!is.null(train.expr)) {
     train.expr <- train.expr[, train.samples]
   }

   z.score.drc <- TRUE
   if(z.score.drc) {
       train.drc <- t(scale(t(train.drc), center = TRUE, scale = TRUE))
   }

   ## Fit ridge to training data and apply to validation for each drug independently.
   ret <- llply(1:length(train.drugs),
                .parallel = TRUE,
                .fun = function(drug.i) {
                         train.drug <- train.drugs[drug.i]
                         cat(paste0("Modeling ", train.drug, " using ", train.response.col, "\n"))

                         x.orig <- train.expr
                         y.orig <- as.numeric(train.drc[as.character(train.drug),])
                         flag <- is.na(y.orig)
                         x.orig <- x.orig[, !flag]
                         y.orig <- y.orig[!flag]

                         ## Exclude any genes that have no variation in training data set (or were set to NA/NaN by
                         ## z-scoring above)
                         constant.genes <- rownames(x.orig)[unlist(apply(x.orig, 1, function(row) any(!is.finite(row)) | all(row == row[1])))]
                         x.orig <- x.orig[!(rownames(x.orig) %in% constant.genes),]

                         nfolds <- 5
                         ret.list <- list() 
                         ## alphas <- seq(from = 0, to = 1, by = 0.05)
                         ## alphas <- c(0, 0.25, 0.5, 0.7, 1)
                         N <- length(y.orig)
                           set.seed(seed)
                         foldid <- sample(rep(seq(nfolds), length = N))

                         ret.list <- llply(1:num.bootstraps,
                               .parallel = FALSE,
                               .fun = function(boot.i) {

                                        set.seed(boot.i)
                                        indices <- sample.int(n = N, size = N, replace=TRUE)
                                        x.boot <- x.orig[, indices]
                                        y.boot <- y.orig[indices]

                                        train.models <- llply(alphas, .parallel = FALSE,
                                                              .fun = function(alpha) { 
                                                                       fit <- tryCatch({fit.elastic.net.alpha(alpha, x = t(x.boot), y = as.numeric(y.boot), foldid = foldid, nfolds = nfolds, standardize = TRUE, seed = boot.i)}, 
                                                                                       error = function(e) { return(NA) })
                                                                       fit
                                                              })
                                        ## lambda - lambda.1se below is used because rounding my prevent any lambda from equaling lambda.1se
                                        mse <- unlist(lapply(train.models, function(x) ifelse(class(x) == "cv.glmnet", x$cvm[which(abs(x$lambda - x$lambda.1se) == min(abs(x$lambda - x$lambda.1se)))[1]], NA)))
                                        ## Choose the alpha parameter that gives the lowest MSE
                                        best.model <- NULL
                                        best.alpha <- NULL
                                        coeffs <- NULL
                                        if(any(is.finite(mse)) && any(!is.na(mse))) {
                                          best.alpha <- alphas[which(mse == min(mse, na.rm=TRUE))]
                                          best.model <- train.models[[which(alphas == best.alpha)]]
                                          coeffs <- coef(best.model)
                                        }
##                                        list(train.drug = train.drug, fit = best.model, alpha = best.alpha, model = "glmnet", n.train = N, boot.i = boot.i)
                                        list(train.drug = train.drug, coeffs = coeffs, alpha = best.alpha, model = "glmnet", n.train = N, boot.i = boot.i)
                                      })
                         ret.list
                       })
   ret
} 


## Do all train x test data set x metric (e.g., pearson, spearman) x response (e.g., auc, dss2) x data type (e.g., gene-level, pathway-level expression) comparisons
train.and.test.crossproduct <- function(gene.sets, drug.name.tbl, train.set.names, train.dss.args, train.expr.args, train.genomic.args, train.clinical.args, train.common.drugs, train.drugs, train.drug.cols, train.patient.cols, train.response.cols, test.set.names, test.dss.args, test.expr.args, test.genomic.args, test.clinical.args, test.common.drugs, test.drug.cols, test.patient.cols, test.response.cols, seed = 1234, use.rf = FALSE, use.svm = FALSE, use.ridge = TRUE, use.mean = TRUE, num.processes = 1, alphas = c(0, 1), include.mean.response.as.feature = FALSE, subtract.mean.response = FALSE, keep.forest = TRUE) {

   z.score.expr <- TRUE
   z.score.drc <- TRUE

   ## Define the feature spaces for all of the training and test data.  Make sure they are common across all of these data sets. 
   ## Also sensible might be just that each pairwise train vs test set has the same feature space--but that would not enable
   ## fair comparisons across training sets.
   train.drcs <- list()
   train.exprs <- list()
   train.genomics <- list()
   train.clinicals <- list()

   for(train.indx in 1:length(train.set.names)) {
     train.set <- train.set.names[[train.indx]]
     cat(paste0("Preparing train set: ", train.set, "\n"))
     l <- prepare.drug.response.and.expr.matrices(train.dss.args[[train.indx]], train.expr.args[[train.indx]], drugs = train.common.drugs[[train.indx]], 
                                                  drug.col = train.drug.cols[[train.indx]], patient.col = train.patient.cols[[train.indx]], response.col = train.response.cols[[train.indx]])

     train.drc <- l[["drc.df"]]
     train.expr <- l[["expr.df"]]
     train.samples <- colnames(train.drc)

     train.genomic <- train.genomic.args[[train.indx]]

     if(!is.null(train.genomic)) {
       train.samples <- intersect(train.samples, colnames(train.genomic))
     }

     train.clinical <- train.clinical.args[[train.indx]]

     if(!is.null(train.clinical)) {
       train.samples <- intersect(train.samples, colnames(train.clinical))
     }

     if(!is.na(train.expr) && !is.null(train.expr)) {
       train.samples <- intersect(train.samples, colnames(train.expr))
     }

     train.drc <- train.drc[, train.samples]

     ## Subtract off the mean response _before_ z-scoring
     if(subtract.mean.response) {
       cat("Subtracting mean response\n")
       train.drc.mean <- unname(colMeans(train.drc, na.rm=TRUE))
       train.drc <- train.drc - matrix(train.drc.mean, nrow=nrow(train.drc), ncol=ncol(train.drc), byrow=TRUE)
     }

     if(z.score.drc) {
         # z-score
         train.drc <- t(scale(t(train.drc), center = TRUE, scale = TRUE))
     }

     if(!is.null(train.genomic)) {
       train.genomic <- train.genomic[, train.samples, drop=F]
     }
     if(!is.null(train.clinical)) {
       train.clinical <- train.clinical[, train.samples, drop=F]
     }
     if(!is.na(train.expr) && !is.null(train.expr)) {
       train.expr <- train.expr[, train.samples]
       flag <- unlist(apply(train.expr, 1, function(row) any(!is.finite(row)) || (any(is.na(row))) || (all(row == row[1]))))
       train.expr <- train.expr[!flag, , drop=F]
     }

     train.drcs[[train.set]] <- train.drc
     train.exprs[[train.set]] <- train.expr
     train.genomics[[train.set]] <- train.genomic
     train.clinicals[[train.set]] <- train.clinical
     rm(train.drc)
     rm(train.expr)
     rm(train.genomic)
     rm(train.clinical)
     gc()
   }

   test.drcs <- list()
   test.exprs <- list()
   test.genomics <- list()
   test.clinicals <- list()
   if(length(test.set.names) > 0) {
     for(test.indx in 1:length(test.set.names)) {
       test.set <- test.set.names[[test.indx]]
       cat(paste0("Preparing test set: ", test.set, "\n"))
       l <- prepare.drug.response.and.expr.matrices(test.dss.args[[test.indx]], test.expr.args[[test.indx]], drugs = test.common.drugs[[test.indx]], 
                                                    drug.col = test.drug.cols[[test.indx]], patient.col = test.patient.cols[[test.indx]], response.col = test.response.cols[[test.indx]])

       test.drc <- l[["drc.df"]]
       test.expr <- l[["expr.df"]]

       test.samples <- colnames(test.drc)

       test.genomic <- test.genomic.args[[test.indx]]
  
       if(!is.null(test.genomic)) {
         test.samples <- intersect(test.samples, colnames(test.genomic))
       }
  
       test.clinical <- test.clinical.args[[test.indx]]
  
       if(!is.null(test.clinical)) {
         test.samples <- intersect(test.samples, colnames(test.clinical))
       }
  
       if(!is.na(test.expr) && !is.null(test.expr)) {
         test.samples <- intersect(test.samples, colnames(test.expr))
       }
  
       test.drc <- test.drc[, test.samples]
       ## Subtract off the mean response _before_ z-scoring
       if(subtract.mean.response) {
         test.drc.mean <- unname(colMeans(test.drc, na.rm=TRUE))
         test.drc <- test.drc - matrix(test.drc.mean, nrow=nrow(test.drc), ncol=ncol(test.drc), byrow=TRUE)
       }

       if(z.score.drc) {
           # z-score
           test.drc <- t(scale(t(test.drc), center = TRUE, scale = TRUE))
       }
       if(!is.null(test.genomic)) {
         test.genomic <- test.genomic[, test.samples, drop=F]
       }
       if(!is.null(test.clinical)) {
         test.clinical <- test.clinical[, test.samples, drop=F]
       }
       if(!is.na(test.expr) && !is.null(test.expr)) {
         test.expr <- test.expr[, test.samples]
         flag <- unlist(apply(test.expr, 1, function(row) any(!is.finite(row)) || (any(is.na(row))) || (all(row == row[1]))))
         test.expr <- test.expr[!flag, ,drop=F]
       }
  
       test.drcs[[test.set]] <- test.drc
       test.exprs[[test.set]] <- test.expr
       test.genomics[[test.set]] <- test.genomic
       test.clinicals[[test.set]] <- test.clinical
       rm(test.drc)
       rm(test.expr)
       rm(test.genomic)
       rm(test.clinical)
       gc()
     }
   }

   expr.train.features <- Reduce(intersect, lapply(train.set.names, function(set.name) rownames(train.exprs[[set.name]])))
   expr.features <- expr.train.features

   genomic.train.features <- Reduce(intersect, lapply(train.set.names, function(set.name) rownames(train.genomics[[set.name]])))
   genomic.features <- genomic.train.features

   clinical.train.features <- Reduce(intersect, lapply(train.set.names, function(set.name) rownames(train.clinicals[[set.name]])))
   clinical.features <- clinical.train.features

   if(length(test.set.names) > 0) {
     expr.test.features <- Reduce(intersect, lapply(test.set.names, function(set.name) rownames(test.exprs[[set.name]])))
     expr.features <- intersect(expr.train.features, expr.test.features)

     genomic.test.features <- Reduce(intersect, lapply(test.set.names, function(set.name) rownames(test.genomics[[set.name]])))
     genomic.features <- intersect(genomic.train.features, genomic.test.features)

     clinical.test.features <- Reduce(intersect, lapply(test.set.names, function(set.name) rownames(test.clinicals[[set.name]])))
     clinical.features <- intersect(clinical.train.features, clinical.test.features)
   }

   ## Assemble the design matrices (combination of expression, genomic, and clinical) using only common features.
   for(train.indx in 1:length(train.set.names)) {
     train.set <- train.set.names[[train.indx]]

     if(!is.na(train.exprs[[train.set]]) && !is.null(train.exprs[[train.set]])) {
       train.exprs[[train.set]] <- train.exprs[[train.set]][expr.features, ]
     }

     if(!is.null(train.genomics[[train.set]])) {
       train.genomics[[train.set]] <- train.genomics[[train.set]][genomic.features, ,drop=F]
     }

     if(!is.null(train.clinicals[[train.set]])) {
       train.clinicals[[train.set]] <- train.clinicals[[train.set]][clinical.features, ,drop=F]
     }

   }

   if(length(test.set.names) > 0) {
     for(test.indx in 1:length(test.set.names)) {
       test.set <- test.set.names[[test.indx]]
  
       if(!is.na(test.exprs[[test.set]]) && !is.null(test.exprs[[test.set]])) {
         test.exprs[[test.set]] <- test.exprs[[test.set]][expr.features, ]
       }
  
       if(!is.null(test.genomics[[test.set]])) {
         test.genomics[[test.set]] <- test.genomics[[test.set]][genomic.features, ,drop=F]
       }
  
       if(!is.null(test.clinicals[[test.set]])) {
         test.clinicals[[test.set]] <- test.clinicals[[test.set]][clinical.features, ,drop=F]
       }
  
     }
   }

   train.expr.sets <- list()
   for(train.indx in 1:length(train.set.names)) {
     train.set <- train.set.names[[train.indx]]
     train.expr.sets[[train.set]] <- list()
     for(gene.set in names(gene.sets)) {
       train.expr.sets[[train.set]][[gene.set]] <- train.exprs[[train.set]]
       if(!(gene.set == "gene") && !is.na(train.exprs[[train.set]]) && !is.null(train.exprs[[train.set]])) { 
         cat(paste0("Computing GSVA for training set ", train.set, " and gene set ", gene.set, "\n"))
         suppressPackageStartupMessages(library("GSVA"))
         train.expr.sets[[train.set]][[gene.set]] <- gsva(as.matrix(train.exprs[[train.set]]), gene.sets[[gene.set]], parallel.sz=num.processes, verbose=TRUE)$es.obs
       }
       ## z-score
       if(z.score.expr && !is.na(train.exprs[[train.set]]) && !is.null(train.exprs[[train.set]])) {
         train.expr.sets[[train.set]][[gene.set]] <- t(scale(t(train.expr.sets[[train.set]][[gene.set]]), center = TRUE, scale = TRUE))
       }
       ## Add the mean of the _z-scored_ responses as a feature (but do not zscore them)
       if(include.mean.response.as.feature) {
         y.mean <- (colMeans(train.drcs[[train.set]], na.rm=TRUE))
         if(!is.na(train.expr.sets[[train.set]][[gene.set]]) && !is.null(train.expr.sets[[train.set]][[gene.set]])) {
           train.expr.sets[[train.set]][[gene.set]] <- rbind(train.expr.sets[[train.set]][[gene.set]], "mean.resp" = unname(y.mean[colnames(train.expr.sets[[train.set]][[gene.set]])]))
         } else {
           train.expr.sets[[train.set]][[gene.set]] <- t(data.frame("dummy" = rep(0,length(y.mean)), "mean.resp" = y.mean))
           colnames(train.expr.sets[[train.set]][[gene.set]]) <- names(y.mean)
cat("train\n")
print(head(y.mean))
print(head(train.expr.sets[[train.set]][[gene.set]]))
print(class(train.expr.sets[[train.set]][[gene.set]]))
         }
##cat(paste0("mean response train:\n"))
##print(train.expr.sets[[train.set]][[gene.set]]["mean.resp",])
       }
     }
   }

   test.expr.sets <- list()
   if(length(test.set.names) > 0) {
     for(test.indx in 1:length(test.set.names)) {
       test.set <- test.set.names[[test.indx]]
       test.expr.sets[[test.set]] <- list()
       for(gene.set in names(gene.sets)) {
         test.expr.sets[[test.set]][[gene.set]] <- test.exprs[[test.set]]
         if(!(gene.set == "gene") && !is.na(test.exprs[[test.set]]) && !is.null(test.exprs[[test.set]])) { 
           cat(paste0("Computing GSVA for testing set ", test.set, " and gene set ", gene.set, "\n"))
           test.expr.sets[[test.set]][[gene.set]] <- gsva(as.matrix(test.exprs[[test.set]]), gene.sets[[gene.set]], parallel.sz=num.processes, verbose=TRUE)$es.obs
         }
         ## z-score
         if(z.score.expr && !is.na(test.exprs[[test.set]]) && !is.null(test.exprs[[test.set]])) {
           test.expr.sets[[test.set]][[gene.set]] <- t(scale(t(test.expr.sets[[test.set]][[gene.set]]), center = TRUE, scale = TRUE))
         }
         ## Add the mean of the _z-scored_ responses as a feature (but do not zscore them)
         if(include.mean.response.as.feature) {
           y.mean <- (colMeans(test.drcs[[test.set]], na.rm=TRUE))
           if(!is.na(test.expr.sets[[test.set]][[gene.set]]) && !is.null(test.expr.sets[[test.set]][[gene.set]])) {
             test.expr.sets[[test.set]][[gene.set]] <- rbind(test.expr.sets[[test.set]][[gene.set]], "mean.resp" = unname(y.mean[colnames(test.expr.sets[[test.set]][[gene.set]])]))
           } else {
             test.expr.sets[[test.set]][[gene.set]] <- t(data.frame("dummy" = rep(0,length(y.mean)), "mean.resp" = y.mean))
             colnames(test.expr.sets[[test.set]][[gene.set]]) <- names(y.mean)
cat("test\n")
print(head(y.mean))
print(head(test.expr.sets[[test.set]][[gene.set]]))
           }
cat(paste0("mean response test:\n"))
print(test.expr.sets[[test.set]][[gene.set]]["mean.resp",])
         }
       }
     }
   }

   all.comparisons <- list()
   all.fits <- list()
   for(train.indx in 1:length(train.set.names)) {

     train.set <- train.set.names[[train.indx]]
     all.comparisons[[train.set]] <- list()
     all.fits[[train.set]] <- list()

     if(length(test.set.names) > 0) {
       for(test.indx in 1:length(test.set.names)) {
         test.set <- test.set.names[[test.indx]]
         all.comparisons[[train.set]][[test.set]] <- list()
       }
     }

     for(gene.set in names(gene.sets)) {

       x <- rbind(train.expr.sets[[train.set]][[gene.set]], train.genomics[[train.set]], train.clinicals[[train.set]])
       cat(paste0("Fitting models for train set ", train.set, " and gene.set ", gene.set, "\n"))
       fits <- train.model_(x, train.drcs[[train.set]], train.drugs[[train.indx]], seed = seed, use.rf = use.rf, use.svm = use.svm, use.ridge = use.ridge, use.mean = use.mean, alphas = alphas, keep.forest = keep.forest)
       all.fits[[train.set]][[gene.set]] <- fits
       rm(x)
       gc()

       if(length(test.set.names) > 0) {
         for(test.indx in 1:length(test.set.names)) {
  
           test.set <- test.set.names[[test.indx]]
           cat(paste0("Training on ", train.set, " and testing on ", test.set, " using gene.set ", gene.set, "\n"))
           newx <- rbind(test.expr.sets[[test.set]][[gene.set]], test.genomics[[test.set]], test.clinicals[[test.set]])
           all.comparisons[[train.set]][[test.set]][[gene.set]] <- test.model_(fits, drug.name.tbl, train.drug.cols[[train.indx]], x.arg = newx, test.drc = test.drcs[[test.set]], test.drug.cols[[test.indx]], seed = seed, use.rf = use.rf, use.svm = use.svm, use.ridge = use.ridge, use.mean = use.mean, flatten.downsampled.fits = FALSE)

         } ## for(test.indx in 1:length(test.set.names))
       }
##       rm(fits)
##       gc()
     } ## for(gene.set in names(gene.sets))

   } ## for(train.indx in 1:length(train.set.names)) 
   return(list("all.comparisons" = all.comparisons, "all.fits" = all.fits))
}

prepare.train.and.test.crossproduct <- function(gene.sets, drug.name.tbl, train.set.names, train.dss.args, train.expr.args, train.genomic.args, train.clinical.args, train.common.drugs, train.drugs, train.drug.cols, train.patient.cols, train.response.cols, test.set.names, test.dss.args, test.expr.args, test.genomic.args, test.clinical.args, test.common.drugs, test.drug.cols, test.patient.cols, test.response.cols, seed = 1234, use.rf = FALSE, use.svm = FALSE, use.ridge = TRUE, use.mean = TRUE, num.processes = 1, alphas = c(0, 1), include.mean.response.as.feature = FALSE, subtract.mean.response = FALSE, keep.forest = TRUE, z.score.expr = TRUE, z.score.drc = TRUE) {

   train.x <- list()
   test.x <- list()
   
   ## Define the feature spaces for all of the training and test data.  Make sure they are common across all of these data sets. 
   ## Also sensible might be just that each pairwise train vs test set has the same feature space--but that would not enable
   ## fair comparisons across training sets.
   train.drcs <- list()
   train.exprs <- list()
   train.genomics <- list()
   train.clinicals <- list()

   for(train.indx in 1:length(train.set.names)) {
     train.set <- train.set.names[[train.indx]]
     cat(paste0("Preparing train set: ", train.set, "\n"))
     l <- prepare.drug.response.and.expr.matrices(train.dss.args[[train.indx]], train.expr.args[[train.indx]], drugs = train.common.drugs[[train.indx]], 
                                                  drug.col = train.drug.cols[[train.indx]], patient.col = train.patient.cols[[train.indx]], response.col = train.response.cols[[train.indx]])

     train.drc <- l[["drc.df"]]
     train.expr <- l[["expr.df"]]
     train.samples <- colnames(train.drc)

     train.genomic <- train.genomic.args[[train.indx]]

     if(!is.null(train.genomic)) {
       train.samples <- intersect(train.samples, colnames(train.genomic))
     }

     train.clinical <- train.clinical.args[[train.indx]]

     if(!is.null(train.clinical)) {
       train.samples <- intersect(train.samples, colnames(train.clinical))
     }

     if(!is.null(train.expr)) {
       train.samples <- intersect(train.samples, colnames(train.expr))
     }

     train.drc <- train.drc[, train.samples]

     ## Subtract off the mean response _before_ z-scoring
     if(subtract.mean.response) {
       cat("Subtracting mean response\n")
       train.drc.mean <- unname(colMeans(train.drc, na.rm=TRUE))
       train.drc <- train.drc - matrix(train.drc.mean, nrow=nrow(train.drc), ncol=ncol(train.drc), byrow=TRUE)
     }

     if(z.score.drc) {
         # z-score
         train.drc <- t(scale(t(train.drc), center = TRUE, scale = TRUE))
     }

     if(!is.null(train.genomic)) {
       train.genomic <- train.genomic[, train.samples]
     }
     if(!is.null(train.clinical)) {
       train.clinical <- train.clinical[, train.samples]
     }
     if(!is.null(train.expr)) {
       train.expr <- train.expr[, train.samples]
       flag <- unlist(apply(train.expr, 1, function(row) any(!is.finite(row)) || (any(is.na(row))) || (all(row == row[1]))))
       train.expr <- train.expr[!flag, ]
     }

     train.drcs[[train.set]] <- train.drc
     train.exprs[[train.set]] <- train.expr
     train.genomics[[train.set]] <- train.genomic
     train.clinicals[[train.set]] <- train.clinical
     rm(train.drc)
     rm(train.expr)
     rm(train.genomic)
     rm(train.clinical)
     gc()
   }

   test.drcs <- list()
   test.exprs <- list()
   test.genomics <- list()
   test.clinicals <- list()
   if(length(test.set.names) > 0) {
     for(test.indx in 1:length(test.set.names)) {
       test.set <- test.set.names[[test.indx]]
       cat(paste0("Preparing test set: ", test.set, "\n"))
       l <- prepare.drug.response.and.expr.matrices(test.dss.args[[test.indx]], test.expr.args[[test.indx]], drugs = test.common.drugs[[test.indx]], 
                                                    drug.col = test.drug.cols[[test.indx]], patient.col = test.patient.cols[[test.indx]], response.col = test.response.cols[[test.indx]])

       test.drc <- l[["drc.df"]]
       test.expr <- l[["expr.df"]]

       test.samples <- colnames(test.drc)

       test.genomic <- test.genomic.args[[test.indx]]
  
       if(!is.null(test.genomic)) {
         test.samples <- intersect(test.samples, colnames(test.genomic))
       }
  
       test.clinical <- test.clinical.args[[test.indx]]
  
       if(!is.null(test.clinical)) {
         test.samples <- intersect(test.samples, colnames(test.clinical))
       }
  
       if(!is.null(test.expr)) {
         test.samples <- intersect(test.samples, colnames(test.expr))
       }
  
       test.drc <- test.drc[, test.samples]
       ## Subtract off the mean response _before_ z-scoring
       if(subtract.mean.response) {
         test.drc.mean <- unname(colMeans(test.drc, na.rm=TRUE))
         test.drc <- test.drc - matrix(test.drc.mean, nrow=nrow(test.drc), ncol=ncol(test.drc), byrow=TRUE)
       }

       if(z.score.drc) {
           # z-score
           test.drc <- t(scale(t(test.drc), center = TRUE, scale = TRUE))
       }
       if(!is.null(test.genomic)) {
         test.genomic <- test.genomic[, test.samples]
       }
       if(!is.null(test.clinical)) {
         test.clinical <- test.clinical[, test.samples]
       }
       if(!is.null(test.expr)) {
         test.expr <- test.expr[, test.samples]
         flag <- unlist(apply(test.expr, 1, function(row) any(!is.finite(row)) || (any(is.na(row))) || (all(row == row[1]))))
         test.expr <- test.expr[!flag, ]
       }
  
       test.drcs[[test.set]] <- test.drc
       test.exprs[[test.set]] <- test.expr
       test.genomics[[test.set]] <- test.genomic
       test.clinicals[[test.set]] <- test.clinical
       rm(test.drc)
       rm(test.expr)
       rm(test.genomic)
       rm(test.clinical)
       gc()
     }
   }

   expr.train.features <- Reduce(intersect, lapply(train.set.names, function(set.name) rownames(train.exprs[[set.name]])))
   expr.features <- expr.train.features

   genomic.train.features <- Reduce(intersect, lapply(train.set.names, function(set.name) rownames(train.genomics[[set.name]])))
   genomic.features <- genomic.train.features

   clinical.train.features <- Reduce(intersect, lapply(train.set.names, function(set.name) rownames(train.clinicals[[set.name]])))
   clinical.features <- clinical.train.features

   if(length(test.set.names) > 0) {
     expr.test.features <- Reduce(intersect, lapply(test.set.names, function(set.name) rownames(test.exprs[[set.name]])))
     expr.features <- intersect(expr.train.features, expr.test.features)

     genomic.test.features <- Reduce(intersect, lapply(test.set.names, function(set.name) rownames(test.genomics[[set.name]])))
     genomic.features <- intersect(genomic.train.features, genomic.test.features)

     clinical.test.features <- Reduce(intersect, lapply(test.set.names, function(set.name) rownames(test.clinicals[[set.name]])))
     clinical.features <- intersect(clinical.train.features, clinical.test.features)
   }

   ## Assemble the design matrices (combination of expression, genomic, and clinical) using only common features.
   for(train.indx in 1:length(train.set.names)) {
     train.set <- train.set.names[[train.indx]]

     if(!is.null(train.exprs[[train.set]])) {
       train.exprs[[train.set]] <- train.exprs[[train.set]][expr.features, ]
     }

     if(!is.null(train.genomics[[train.set]])) {
       train.genomics[[train.set]] <- train.genomics[[train.set]][genomic.features, ]
     }

     if(!is.null(train.clinicals[[train.set]])) {
       train.clinicals[[train.set]] <- train.clinicals[[train.set]][clinical.features, ]
     }

   }

   if(length(test.set.names) > 0) {
     for(test.indx in 1:length(test.set.names)) {
       test.set <- test.set.names[[test.indx]]
  
       if(!is.null(test.exprs[[test.set]])) {
         test.exprs[[test.set]] <- test.exprs[[test.set]][expr.features, ]
       }
  
       if(!is.null(test.genomics[[test.set]])) {
         test.genomics[[test.set]] <- test.genomics[[test.set]][genomic.features, ]
       }
  
       if(!is.null(test.clinicals[[test.set]])) {
         test.clinicals[[test.set]] <- test.clinicals[[test.set]][clinical.features, ]
       }
  
     }
   }

   train.expr.sets <- list()
   for(train.indx in 1:length(train.set.names)) {
     train.set <- train.set.names[[train.indx]]
     train.expr.sets[[train.set]] <- list()
     for(gene.set in names(gene.sets)) {
       if(gene.set == "gene") { 
         train.expr.sets[[train.set]][[gene.set]] <- train.exprs[[train.set]]
       } else {
         cat(paste0("Computing GSVA for training set ", train.set, " and gene set ", gene.set, "\n"))
         suppressPackageStartupMessages(library("GSVA"))
         train.expr.sets[[train.set]][[gene.set]] <- gsva(as.matrix(train.exprs[[train.set]]), gene.sets[[gene.set]], parallel.sz=num.processes, verbose=TRUE)$es.obs
       }
       ## z-score
       if(z.score.expr) {
         train.expr.sets[[train.set]][[gene.set]] <- t(scale(t(train.expr.sets[[train.set]][[gene.set]]), center = TRUE, scale = TRUE))
       }
       ## Add the mean of the _z-scored_ responses as a feature (but do not zscore them)
       if(include.mean.response.as.feature) {
         y.mean <- (colMeans(train.drcs[[train.set]], na.rm=TRUE))
         train.expr.sets[[train.set]][[gene.set]] <- rbind(train.expr.sets[[train.set]][[gene.set]], "mean.resp" = unname(y.mean[colnames(train.expr.sets[[train.set]][[gene.set]])]))
##cat(paste0("mean response train:\n"))
##print(train.expr.sets[[train.set]][[gene.set]]["mean.resp",])
       }
     }
   }

   test.expr.sets <- list()
   if(length(test.set.names) > 0) {
     for(test.indx in 1:length(test.set.names)) {
       test.set <- test.set.names[[test.indx]]
       test.expr.sets[[test.set]] <- list()
       for(gene.set in names(gene.sets)) {
         if(gene.set == "gene") { 
           test.expr.sets[[test.set]][[gene.set]] <- test.exprs[[test.set]]
         } else {
           cat(paste0("Computing GSVA for testing set ", test.set, " and gene set ", gene.set, "\n"))
           test.expr.sets[[test.set]][[gene.set]] <- gsva(as.matrix(test.exprs[[test.set]]), gene.sets[[gene.set]], parallel.sz=num.processes, verbose=TRUE)$es.obs
         }
         ## z-score
         if(z.score.expr) {
           test.expr.sets[[test.set]][[gene.set]] <- t(scale(t(test.expr.sets[[test.set]][[gene.set]]), center = TRUE, scale = TRUE))
         }
         ## Add the mean of the _z-scored_ responses as a feature (but do not zscore them)
         if(include.mean.response.as.feature) {
           y.mean <- (colMeans(test.drcs[[test.set]], na.rm=TRUE))
           test.expr.sets[[test.set]][[gene.set]] <- rbind(test.expr.sets[[test.set]][[gene.set]], "mean.resp" = unname(y.mean[colnames(test.expr.sets[[test.set]][[gene.set]])]))
cat(paste0("mean response test:\n"))
print(test.expr.sets[[test.set]][[gene.set]]["mean.resp",])
         }
       }
     }
   }

   all.comparisons <- list()
   all.fits <- list()
   for(train.indx in 1:length(train.set.names)) {

     train.set <- train.set.names[[train.indx]]
     train.x[[train.set]] <- list()

     for(gene.set in names(gene.sets)) {

       x <- rbind(train.expr.sets[[train.set]][[gene.set]], train.genomics[[train.set]], train.clinicals[[train.set]])
       train.x[[train.set]][[gene.set]] <- x

     }
   }

   if(length(test.set.names) > 0) {
     for(test.indx in 1:length(test.set.names)) {
       test.set <- test.set.names[[test.indx]]
       test.x[[test.set]] <- list()

       for(gene.set in names(gene.sets)) {
         newx <- rbind(test.expr.sets[[test.set]][[gene.set]], test.genomics[[test.set]], test.clinicals[[test.set]])
         test.x[[test.set]][[gene.set]] <- newx
       }
     }
   }

   return(list("train.x" = train.x, "train.drcs" = train.drcs, "test.x" = test.x, "test.drcs" = test.drcs))
}

train.and.test.crossproduct_ <- function(gene.sets, drug.name.tbl, train.set.names, train.y, train.x, train.drugs, train.drug.cols, test.set.names, test.y, test.x, test.drug.cols, seed = 1234, use.rf = FALSE, use.svm = FALSE, use.ridge = TRUE, use.mean = TRUE, num.processes = 1, alphas = c(0, 1), keep.forest = TRUE, use.glmnet.intercept = TRUE) {

   all.comparisons <- list()
   all.fits <- list()
   for(train.indx in 1:length(train.set.names)) {

     train.set <- train.set.names[[train.indx]]
     all.comparisons[[train.set]] <- list()
     all.fits[[train.set]] <- list()

     if(length(test.set.names) > 0) {
       for(test.indx in 1:length(test.set.names)) {
         test.set <- test.set.names[[test.indx]]
         all.comparisons[[train.set]][[test.set]] <- list()
       }
     }

     for(gene.set in names(gene.sets)) {

       x <- train.x[[train.set]][[gene.set]]
       cat(paste0("Fitting models for train set ", train.set, " and gene.set ", gene.set, "\n"))
       fits <- train.model_(x, train.y[[train.set]], train.drugs[[train.indx]], seed = seed, use.rf = use.rf, use.svm = use.svm, use.ridge = use.ridge, use.mean = use.mean, alphas = alphas, keep.forest = keep.forest, use.glmnet.intercept = use.glmnet.intercept)
       all.fits[[train.set]][[gene.set]] <- fits
       rm(x)
       gc()

       if(length(test.set.names) > 0) {
         for(test.indx in 1:length(test.set.names)) {
  
           test.set <- test.set.names[[test.indx]]
           cat(paste0("Training on ", train.set, " and testing on ", test.set, " using gene.set ", gene.set, "\n"))
           newx <- test.x[[test.set]][[gene.set]]
           all.comparisons[[train.set]][[test.set]][[gene.set]] <- test.model_(fits, drug.name.tbl, train.drug.cols[[train.indx]], x.arg = newx, test.drc = test.y[[test.set]], test.drug.cols[[test.indx]], seed = seed, use.rf = use.rf, use.svm = use.svm, use.ridge = use.ridge, use.mean = use.mean, flatten.downsampled.fits = FALSE)

         } ## for(test.indx in 1:length(test.set.names))
       }
##       rm(fits)
##       gc()
     } ## for(gene.set in names(gene.sets))

   } ## for(train.indx in 1:length(train.set.names)) 
   return(list("all.comparisons" = all.comparisons, "all.fits" = all.fits))
}

train.and.test.crossproduct.with.downsampling <- function(gene.sets, drug.name.tbl, train.set.names, train.dss.args, train.expr.args, train.genomic.args, train.clinical.args, train.common.drugs, train.drugs, train.drug.cols, train.patient.cols, train.response.cols, test.set.names, test.dss.args, test.expr.args, test.genomic.args, test.clinical.args, test.common.drugs, test.drug.cols, test.patient.cols, test.response.cols, seed = 1234, use.rf = FALSE, use.svm = FALSE, use.ridge = TRUE, use.mean = TRUE, num.processes = 1, num.train.samples.per.drug, num.test.samples.per.drug, num.iterations = 100, with.replacement = FALSE, return.fits = TRUE, include.mean.response.as.feature = FALSE, subtract.mean.response = FALSE, keep.forest = TRUE) {
   z.score.expr <- TRUE
   z.score.drc <- TRUE

   ## Define the feature spaces for all of the training and test data.  Make sure they are common across all of these data sets. 
   ## Also sensible might be just that each pairwise train vs test set has the same feature space--but that would not enable
   ## fair comparisons across training sets.
   train.drcs <- list()
   train.exprs <- list()
   train.genomics <- list()
   train.clinicals <- list()

   for(train.indx in 1:length(train.set.names)) {
     train.set <- train.set.names[[train.indx]]
     cat(paste0("Preparing train set: ", train.set, "\n"))
     l <- prepare.drug.response.and.expr.matrices(train.dss.args[[train.indx]], train.expr.args[[train.indx]], drugs = train.common.drugs[[train.indx]], 
                                                  drug.col = train.drug.cols[[train.indx]], patient.col = train.patient.cols[[train.indx]], response.col = train.response.cols[[train.indx]])

     train.drc <- l[["drc.df"]]
     train.expr <- l[["expr.df"]]
     train.samples <- colnames(train.drc)

     train.genomic <- train.genomic.args[[train.indx]]

     if(!is.null(train.genomic)) {
       train.samples <- intersect(train.samples, colnames(train.genomic))
     }

     train.clinical <- train.clinical.args[[train.indx]]

     if(!is.null(train.clinical)) {
       train.samples <- intersect(train.samples, colnames(train.clinical))
     }

     if(!is.null(train.expr)) {
       train.samples <- intersect(train.samples, colnames(train.expr))
     }

     train.drc <- train.drc[, train.samples]

## We will zscore below after downsampling
##     if(z.score.drc) {
##         # z-score
##         train.drc <- t(scale(t(train.drc), center = TRUE, scale = TRUE))
##     }

     if(!is.null(train.genomic)) {
       train.genomic <- train.genomic[, train.samples]
     }
     if(!is.null(train.clinical)) {
       train.clinical <- train.clinical[, train.samples]
     }
     if(!is.null(train.expr)) {
       train.expr <- train.expr[, train.samples]
       flag <- unlist(apply(train.expr, 1, function(row) any(!is.finite(row)) || (any(is.na(row))) || (all(row == row[1]))))
       train.expr <- train.expr[!flag, ]
     }

     train.drcs[[train.set]] <- train.drc
     train.exprs[[train.set]] <- train.expr
     train.genomics[[train.set]] <- train.genomic
     train.clinicals[[train.set]] <- train.clinical
     rm(train.drc)
     rm(train.expr)
     rm(train.genomic)
     rm(train.clinical)
     gc()
   }

   test.drcs <- list()
   test.exprs <- list()
   test.genomics <- list()
   test.clinicals <- list()
   if(length(test.set.names) > 0) {
     for(test.indx in 1:length(test.set.names)) {
       test.set <- test.set.names[[test.indx]]
       cat(paste0("Preparing test set: ", test.set, "\n"))
       l <- prepare.drug.response.and.expr.matrices(test.dss.args[[test.indx]], test.expr.args[[test.indx]], drugs = test.common.drugs[[test.indx]], 
                                                    drug.col = test.drug.cols[[test.indx]], patient.col = test.patient.cols[[test.indx]], response.col = test.response.cols[[test.indx]])

       test.drc <- l[["drc.df"]]
       test.expr <- l[["expr.df"]]

       test.samples <- colnames(test.drc)

       test.genomic <- test.genomic.args[[test.indx]]
  
       if(!is.null(test.genomic)) {
         test.samples <- intersect(test.samples, colnames(test.genomic))
       }
  
       test.clinical <- test.clinical.args[[test.indx]]
  
       if(!is.null(test.clinical)) {
         test.samples <- intersect(test.samples, colnames(test.clinical))
       }
  
       if(!is.null(test.expr)) {
         test.samples <- intersect(test.samples, colnames(test.expr))
       }
  
       test.drc <- test.drc[, test.samples]
## We will zscore below after downsampling
##       if(z.score.drc) {
##           # z-score
##           test.drc <- t(scale(t(test.drc), center = TRUE, scale = TRUE))
##       }
  
       if(!is.null(test.genomic)) {
         test.genomic <- test.genomic[, test.samples]
       }
       if(!is.null(test.clinical)) {
         test.clinical <- test.clinical[, test.samples]
       }
       if(!is.null(test.expr)) {
         test.expr <- test.expr[, test.samples]
         flag <- unlist(apply(test.expr, 1, function(row) any(!is.finite(row)) || (any(is.na(row))) || (all(row == row[1]))))
         test.expr <- test.expr[!flag, ]
       }
  
       test.drcs[[test.set]] <- test.drc
       test.exprs[[test.set]] <- test.expr
       test.genomics[[test.set]] <- test.genomic
       test.clinicals[[test.set]] <- test.clinical
       rm(test.drc)
       rm(test.expr)
       rm(test.genomic)
       rm(test.clinical)
       gc()
     }
   }

   expr.train.features <- Reduce(intersect, lapply(train.set.names, function(set.name) rownames(train.exprs[[set.name]])))
   expr.features <- expr.train.features

   genomic.train.features <- Reduce(intersect, lapply(train.set.names, function(set.name) rownames(train.genomics[[set.name]])))
   genomic.features <- genomic.train.features

   clinical.train.features <- Reduce(intersect, lapply(train.set.names, function(set.name) rownames(train.clinicals[[set.name]])))
   clinical.features <- clinical.train.features

   if(length(test.set.names) > 0) {
     expr.test.features <- Reduce(intersect, lapply(test.set.names, function(set.name) rownames(test.exprs[[set.name]])))
     expr.features <- intersect(expr.train.features, expr.test.features)

     genomic.test.features <- Reduce(intersect, lapply(test.set.names, function(set.name) rownames(test.genomics[[set.name]])))
     genomic.features <- intersect(genomic.train.features, genomic.test.features)

     clinical.test.features <- Reduce(intersect, lapply(test.set.names, function(set.name) rownames(test.clinicals[[set.name]])))
     clinical.features <- intersect(clinical.train.features, clinical.test.features)
   }

   ## Assemble the design matrices (combination of expression, genomic, and clinical) using only common features.
   for(train.indx in 1:length(train.set.names)) {
     train.set <- train.set.names[[train.indx]]

     if(!is.null(train.exprs[[train.set]])) {
       train.exprs[[train.set]] <- train.exprs[[train.set]][expr.features, ]
     }

     if(!is.null(train.genomics[[train.set]])) {
       train.genomics[[train.set]] <- train.genomics[[train.set]][genomic.features, ]
     }

     if(!is.null(train.clinicals[[train.set]])) {
       train.clinicals[[train.set]] <- train.clinicals[[train.set]][clinical.features, ]
     }

   }

   if(length(test.set.names) > 0) {
     for(test.indx in 1:length(test.set.names)) {
       test.set <- test.set.names[[test.indx]]
  
       if(!is.null(test.exprs[[test.set]])) {
         test.exprs[[test.set]] <- test.exprs[[test.set]][expr.features, ]
       }
  
       if(!is.null(test.genomics[[test.set]])) {
         test.genomics[[test.set]] <- test.genomics[[test.set]][genomic.features, ]
       }
  
       if(!is.null(test.clinicals[[test.set]])) {
         test.clinicals[[test.set]] <- test.clinicals[[test.set]][clinical.features, ]
       }
  
     }
   }

   train.expr.sets <- list()
   for(train.indx in 1:length(train.set.names)) {
     train.set <- train.set.names[[train.indx]]
     train.expr.sets[[train.set]] <- list()
     for(gene.set in names(gene.sets)) {
       if(gene.set == "gene") { 
         train.expr.sets[[train.set]][[gene.set]] <- train.exprs[[train.set]]
       } else {
         cat(paste0("Computing GSVA for training set ", train.set, " and gene set ", gene.set, "\n"))
         suppressPackageStartupMessages(library("GSVA"))
         train.expr.sets[[train.set]][[gene.set]] <- gsva(as.matrix(train.exprs[[train.set]]), gene.sets[[gene.set]], parallel.sz=num.processes, verbose=TRUE)$es.obs
       }
## We will zscore below after downsampling
##       ## z-score
##       if(z.score.expr) {
##         train.expr.sets[[train.set]][[gene.set]] <- t(scale(t(train.expr.sets[[train.set]][[gene.set]]), center = TRUE, scale = TRUE))
##       }
     }
   }

   test.expr.sets <- list()
   if(length(test.set.names) > 0) {
     for(test.indx in 1:length(test.set.names)) {
       test.set <- test.set.names[[test.indx]]
       test.expr.sets[[test.set]] <- list()
       for(gene.set in names(gene.sets)) {
         if(gene.set == "gene") { 
           test.expr.sets[[test.set]][[gene.set]] <- test.exprs[[test.set]]
         } else {
           cat(paste0("Computing GSVA for testing set ", test.set, " and gene set ", gene.set, "\n"))
           test.expr.sets[[test.set]][[gene.set]] <- gsva(as.matrix(test.exprs[[test.set]]), gene.sets[[gene.set]], parallel.sz=num.processes, verbose=TRUE)$es.obs
         }
## We will zscore below after downsampling
##         ## z-score
##         if(z.score.expr) {
##           test.expr.sets[[test.set]][[gene.set]] <- t(scale(t(test.expr.sets[[test.set]][[gene.set]]), center = TRUE, scale = TRUE))
##         }
       }
     }
   }

   all.comparisons <- list()
   all.fits <- list()
   for(train.indx in 1:length(train.set.names)) {

     train.set <- train.set.names[[train.indx]]
     all.comparisons[[train.set]] <- list()
     if(return.fits) {
       all.fits[[train.set]] <- list()
     }
     if(length(test.set.names) > 0) {
       for(test.indx in 1:length(test.set.names)) {
         test.set <- test.set.names[[test.indx]]
         all.comparisons[[train.set]][[test.set]] <- list()
       }
     }

     for(gene.set in names(gene.sets)) {

       cat(paste0("Fitting models for train set ", train.set, " and gene.set ", gene.set, "\n"))
       fits <- train.model.with.downsampling_(train.expr = train.expr.sets[[train.set]][[gene.set]], train.genomic = train.genomics[[train.set]][[gene.set]],
                                              train.clinical = train.clinicals[[train.set]][[gene.set]], train.drcs[[train.set]], train.drugs[[train.indx]], 
                                              train.drug.col = train.drug.cols[[train.indx]], train.response.col = train.response.cols[[train.indx]],
                                              seed = seed, use.rf = use.rf, use.svm = use.svm, use.ridge = use.ridge, use.mean = use.mean,
                                              num.samples.per.drug = num.train.samples.per.drug, num.iterations = num.iterations, with.replacement = with.replacement,
                                              include.mean.response.as.feature = include.mean.response.as.feature, subtract.mean.response = subtract.mean.response, keep.forest = keep.forest)
       if(return.fits) { all.fits[[train.set]][[gene.set]] <- fits }

       if(length(test.set.names) > 0) {
         for(test.indx in 1:length(test.set.names)) {
  
           test.set <- test.set.names[[test.indx]]
           cat(paste0("Training on ", train.set, " and testing on ", test.set, " using gene.set ", gene.set, "\n"))
           all.comparisons[[train.set]][[test.set]][[gene.set]] <-  
                test.model.with.downsampling_(fits, drug.name.tbl, train.drug.cols[[train.indx]], test.expr = test.expr.sets[[test.set]][[gene.set]], 
                                              test.genomic = test.genomics[[test.set]], test.clinical = test.clinicals[[test.set]],
                                              test.drc = test.drcs[[test.set]], test.drug.cols[[test.indx]], 
                                              seed = seed, use.rf = use.rf, use.svm = use.svm, use.ridge = use.ridge, use.mean = use.mean, flatten.downsampled.fits = TRUE, 
                                              num.samples.per.drug = num.test.samples.per.drug, with.replacement = with.replacement,
                                              include.mean.response.as.feature = include.mean.response.as.feature, subtract.mean.response = subtract.mean.response)

         } ## for(test.indx in 1:length(test.set.names))
       }
       
       rm(fits)
       gc()
     } ## for(gene.set in names(gene.sets))

   } ## for(train.indx in 1:length(train.set.names)) 
   return(list("all.comparisons" = all.comparisons, "all.fits" = all.fits))
}

train.model.with.downsampling <- function(train.dss.arg, train.expr.arg, train.genomic.arg, train.clinical.arg, train.common.drugs, train.drugs, train.drug.col, train.patient.col, train.response.col, seed = 1234, use.rf = FALSE, use.svm = FALSE, calc.coefficient.pvals = FALSE, num.samples.per.drug, num.iterations = 100, with.replacement = FALSE, alphas = c(0, 1), keep.forest = TRUE) {  

   l <- prepare.drug.response.and.expr.matrices(train.dss.arg, train.expr.arg, drugs = train.common.drugs, drug.col = train.drug.col, patient.col = train.patient.col, response.col = train.response.col)
   train.drc <- l[["drc.df"]]
   train.expr <- l[["expr.df"]]

   if(!is.null(train.expr)) {
     ## Z-score the genes/rows (individually)
     z.score <- TRUE
     if(z.score) {
       ## z-score training expression data below after sub-sampling
##       train.expr <- t(scale(t(train.expr), center = TRUE, scale = TRUE))
     }
   }

   train.samples <- colnames(train.drc)

   train.genomic <- train.genomic.arg

   if(!is.null(train.genomic)) {
     train.samples <- intersect(train.samples, colnames(train.genomic))
   }

   train.clinical <- train.clinical.arg

   if(!is.null(train.clinical)) {
     train.samples <- intersect(train.samples, colnames(train.clinical))
   }

   if(!is.null(train.expr)) {
     train.samples <- intersect(train.samples, colnames(train.expr))
   }

   train.drc <- train.drc[, train.samples]

   if(!is.null(train.genomic)) {
     train.genomic <- train.genomic[, train.samples, drop = FALSE]
   }
   if(!is.null(train.clinical)) {
     train.clinical <- train.clinical[, train.samples, drop = FALSE]
   }
   if(!is.null(train.expr)) {
     train.expr <- train.expr[, train.samples, drop = FALSE]
   }

   z.score.drc <- TRUE
   if(z.score.drc) {
       ## z-score drug data below after sub-sampling
##       train.drc <- t(scale(t(train.drc), center = TRUE, scale = TRUE))
   }

   ## Fit ridge to training data and apply to validation for each drug independently.
   ret <- llply(1:length(train.drugs),
                .parallel = TRUE,
                .fun = function(drug.i) {
                         llply(1:num.iterations, 
                               .parallel = FALSE,
                               .fun = function(iter) {
                                        train.drug <- train.drugs[drug.i]

                                        ## Z-score training _after_ sub-sampling
                                        y <- as.numeric(train.drc[as.character(train.drug),])
                                        flag <- is.na(y)
                                        y <- y[!flag]

                                        t.expr <- train.expr
                                        t.genomic <- train.genomic
                                        t.clinical <- train.clinical
                                        if(!is.null(t.expr)) { t.expr <- t.expr[, !flag, drop = FALSE] }
                                        if(!is.null(t.genomic)) { t.genomic <- t.genomic[, !flag, drop = FALSE] }
                                        if(!is.null(t.clinical)) { t.clinical <- t.clinical[, !flag, drop = FALSE] }

                                        set.seed(iter)
                                        n.samples <- num.samples.per.drug[num.samples.per.drug[, train.drug.col] == as.character(train.drug), "num.samples"]
                                        cat(paste0("Modeling ", train.drug, " using ", train.response.col, " and num samples = ", n.samples, "\n"))
                                        if(n.samples > length(y)) {
                                          stop("Attempting to model ", train.drug, " using sampling ", ifelse(with.replacement, "with", "without"), " replacement of ", n.samples, " samples, but length of response is only ", length(y), "\n") 
                                        }
                                        indices <- sample.int(length(y), size = n.samples, replace = with.replacement)
                                        y <- y[indices]
                                        if(!is.null(t.expr)) { t.expr <- t.expr[, indices, drop = FALSE] }
                                        if(!is.null(t.genomic)) { t.genomic <- t.genomic[, indices, drop = FALSE] }
                                        if(!is.null(t.clinical)) { t.clinical <- t.clinical[, indices, drop = FALSE] }

                                        if(z.score.drc) {
                                          y <- (y - mean(y)) / sd(y)
                                        }

                                        if(z.score) {
                                          t.expr <- t(scale(t(t.expr), center = TRUE, scale = TRUE))
                                        }

                                        x <- rbind(t.expr, t.genomic, t.clinical)
                                        n.train <- length(y)

                                        ## Exclude any genes that have no variation in training data set (or were set to NA/NaN by
                                        ## z-scoring above)
                                        constant.genes <- rownames(x)[unlist(apply(x, 1, function(row) any(!is.finite(row)) | all(row == row[1])))]
                                        x <- x[!(rownames(x) %in% constant.genes),]

                                        nfolds <- 5 
                                        ret.list <- list()
                                        ## Apply fit to test/validation data

                                        ## Use glmnet
                                        ## alphas <- c(0, 0.25, 0.5, 0.7, 1)
                                        ## alphas <- c(0, 1)
                                        N <- length(y)
                           set.seed(seed)
                                        foldid <- sample(rep(seq(nfolds), length = N))
                                        for(alpha in alphas) {
                                          pvals <- NA
                                          if( (alpha == 1) && (calc.coefficient.pvals == TRUE) ) {
                                            ## lasso
                                            ## Calculate the p-values for the individual components
                                            fit.lasso <- lasso.proj(t(x), y)
                                            pvals <- fit.lasso$pval
                                          }
                                          set.seed(seed)
                                          fit <- tryCatch({cv.glmnet(x = t(x), y, family = "gaussian", type.measure = "mse", nfolds = nfolds, foldid = foldid, alpha = alpha, standardize = FALSE, parallel = TRUE)}, error = function(e) { cat("glmnet err\n"); return(NA) })
                                          ret.list[[length(ret.list)+1]] <- list(train.drug = train.drug, fit = fit, alpha = alpha, model = "glmnet", n.train = n.train, pvals = pvals, iter = iter)
                                        }

                                        ## Use RF
                                        if(use.rf) {
                                          set.seed(seed)
                                          fit <- tryCatch({randomForest(x = t(x), y, importance = TRUE, keep.forest = keep.forest)}, error = function(e) { cat("rf error\n"); return(NA) })
                                          ret.list[[length(ret.list)+1]] <- list(train.drug = train.drug, fit = fit, alpha = NA, model = "rf", n.train = n.train, pvals = NA, iter = iter)
                                        }

                                        ## Use SVM
                                        if(use.svm) {
                                          set.seed(seed)
                                          fit <- tryCatch({svm(x = t(x), y, scale = FALSE)}, error = function(e) { cat("svm error\n"); return(NA) })
                                          ret.list[[length(ret.list)+1]] <- list(train.drug = train.drug, fit = fit, alpha = NA, model = "svm", n.train = n.train, pvals = NA, iter = iter)
                                        }

                                        ## Return a baseline in which the prediction is just the mean of training data set
                                        pred <- mean(y, na.rm = TRUE)
                                        ret.list[[length(ret.list)+1]] <- list(train.drug = train.drug, fit = pred, alpha = NA, model = "mean", n.train = n.train, pvals = NA, iter = iter)


                                        ret.list
                               })
                })
   ret
} 

                                              

model.train.and.test.with.downsampling <- function(train.dss.arg, train.expr.arg, train.genomic.arg, train.clinical.arg, train.common.drugs, train.drugs, train.drug.col, train.patient.col, train.response.col, test.dss.arg, test.expr.arg, test.genomic.arg, test.clinical.arg, test.common.drugs, test.drugs, test.drug.col, test.patient.col, test.response.col, seed = 1234, use.rf = FALSE, use.svm = FALSE, num.samples.per.drug, num.iterations = 100, with.replacement = FALSE, alphas = c(0, 1)) {  

   l <- prepare.drug.response.and.expr.matrices(train.dss.arg, train.expr.arg, drugs = train.common.drugs, drug.col = train.drug.col, patient.col = train.patient.col, response.col = train.response.col)
   train.drc <- l[["drc.df"]]
   train.expr <- l[["expr.df"]]

   l <- prepare.drug.response.and.expr.matrices(test.dss.arg, test.expr.arg, drugs = test.common.drugs, drug.col = test.drug.col, patient.col = test.patient.col, response.col = test.response.col)
   test.drc <- l[["drc.df"]]
   test.expr <- l[["expr.df"]]
 
   if(is.null(train.expr) && !is.null(test.expr)) {
     stop("train expr is null, but FIMM expr is not null!\n")
   }
   if(!is.null(train.expr) && is.null(test.expr)) {
     stop("train expr is not null, but FIMM expr is null!\n")
   }
   if(!is.null(train.expr) && !is.null(test.expr)) {
     common.genes <- intersect(rownames(train.expr), rownames(test.expr))
     train.expr <- train.expr[common.genes,]
     test.expr <- test.expr[common.genes,]

   }

   train.samples <- colnames(train.drc)
   test.samples <- colnames(test.drc)

   train.genomic <- train.genomic.arg
   test.genomic <- test.genomic.arg
   if(is.null(train.genomic) && !is.null(test.genomic)) {
     stop("train genomic is null, but FIMM genomic is not null!\n")
   }
   if(!is.null(train.genomic) && is.null(test.genomic)) {
     stop("train genomic is not null, but FIMM genomic is null!\n")
   }

   if(!is.null(train.genomic) && !is.null(test.genomic)) {
     common.genes <- intersect(rownames(train.genomic), rownames(test.genomic))

     train.genomic <- train.genomic[common.genes,]
     test.genomic <- test.genomic[common.genes,]

     train.samples <- intersect(train.samples, colnames(train.genomic))
     test.samples <- intersect(test.samples, colnames(test.genomic))
   }

   train.clinical <- train.clinical.arg
   test.clinical <- test.clinical.arg
   if(is.null(train.clinical) && !is.null(test.clinical)) {
     stop("train clinical is null, but FIMM clinical is not null!\n")
   }
   if(!is.null(train.clinical) && is.null(test.clinical)) {
     stop("train clinical is not null, but FIMM clinical is null!\n")
   }

   if(!is.null(train.clinical) && !is.null(test.clinical)) {
     common.annotations <- intersect(rownames(train.clinical), rownames(test.clinical))

     train.clinical <- train.clinical[common.annotations,]
     test.clinical <- test.clinical[common.annotations,]

     train.samples <- intersect(train.samples, colnames(train.clinical))
     test.samples <- intersect(test.samples, colnames(test.clinical))
   }

   if(!is.null(train.expr) && !is.null(test.expr)) {
     train.samples <- intersect(train.samples, colnames(train.expr))
     test.samples <- intersect(test.samples, colnames(test.expr))
   }

   train.drc <- train.drc[, train.samples]
   test.drc <- test.drc[, test.samples]

   if(!is.null(train.genomic) && !is.null(test.genomic)) {
     train.genomic <- train.genomic[, train.samples]
     test.genomic <- test.genomic[, test.samples]
   }
   if(!is.null(train.clinical) && !is.null(test.clinical)) {
     train.clinical <- train.clinical[, train.samples]
     test.clinical <- test.clinical[, test.samples]
   }

   z.score <- TRUE
   if(!is.null(train.expr) && !is.null(test.expr)) {
     train.expr <- train.expr[, train.samples]
     test.expr <- test.expr[, test.samples]

     ## Z-score the genes/rows (individually)
     if(z.score) {
       ## z-score training expression data below after sub-sampling
       ## train.expr <- t(scale(t(train.expr), center = TRUE, scale = TRUE))
       test.expr <- t(scale(t(test.expr), center = TRUE, scale = TRUE))
     }


   }

   ## Adjust the location (mean or median) and scale (mad or sd) of the two drug-response data sets to be the same--
   ## based on all drugs in common, even if we don't test them here.
##   test.location <- median(unlist(test.drc), na.rm=TRUE)
##   train.location <- median(unlist(train.drc), na.rm=TRUE)
##   test.scale <- mad(unlist(test.drc), na.rm=TRUE)
##   train.scale <- mad(unlist(train.drc), na.rm=TRUE)

   ## Adjust by mean and sd
##   test.location <- mean(unlist(test.drc), na.rm=TRUE)
##   train.location <- mean(unlist(train.drc), na.rm=TRUE)
##   test.scale <- sd(unlist(test.drc), na.rm=TRUE)
##   train.scale <- sd(unlist(train.drc), na.rm=TRUE)

   ## test.drc <- ( train.scale * ( test.drc - test.location ) / test.scale ) + train.location

   z.score.drc <- TRUE
   if(z.score.drc) {
##       train.drc <- t(scale(t(train.drc), center = TRUE, scale = TRUE))
       test.drc <- t(scale(t(test.drc), center = TRUE, scale = TRUE))
   }

   ## Fit ridge to training data and apply to validation for each drug independently.
   ret <- llply(1:length(train.drugs),
                .parallel = TRUE,
                .fun = function(drug.i) {
                         llply(1:num.iterations, 
                               .parallel = FALSE,
                               .fun = function(iter) {
                                        train.drug <- train.drugs[drug.i]
                                        test.drug <- test.drugs[drug.i]
                                        cat(paste0("Modeling ", train.drug, " using ", train.response.col, "\n"))

                                        ## Z-score training _after_ sub-sampling
                                        y <- as.numeric(train.drc[as.character(train.drug),])
                                        flag <- is.na(y)
                                        y <- y[!flag]

                                        t.expr <- train.expr
                                        t.genomic <- train.genomic
                                        t.clinical <- train.clinical
                                        if(!is.null(t.expr)) { t.expr <- t.expr[, !flag] }
                                        if(!is.null(t.genomic)) { t.genomic <- t.genomic[, !flag] }
                                        if(!is.null(t.clinical)) { t.clinical <- t.clinical[, !flag] }

                                        set.seed(iter)
                                        n.samples <- num.samples.per.drug[as.character(test.drug)]
                                        indices <- sample.int(length(y), size = n.samples, replace = with.replacement)
                                        y <- y[indices]
                                        if(!is.null(t.expr)) { t.expr <- t.expr[, indices] }
                                        if(!is.null(t.genomic)) { t.genomic <- t.genomic[, indices] }
                                        if(!is.null(t.clinical)) { t.clinical <- t.clinical[, indices] }

                                        if(z.score.drc) {
                                          y <- (y - mean(y)) / sd(y)
                                        }

                                        if(z.score) {
                                          t.expr <- t(scale(t(t.expr), center = TRUE, scale = TRUE))
                                        }

                                        x <- rbind(t.expr, t.genomic, t.clinical)

                                        newx <- rbind(test.expr, test.genomic, test.clinical)
                                        newy <- as.numeric(test.drc[as.character(test.drug),])
                                        flag <- is.na(newy)
                                        newx <- newx[, !flag]
                                        newy <- newy[!flag]

                                        n.train <- length(y)
                                        n.test <- length(newy)

                                        ## Exclude any genes that have no variation in either data set (or were set to NA/NaN by
                                        ## z-scoring above)
                                        constant.genes <- rownames(x)[unlist(apply(x, 1, function(row) any(!is.finite(row)) | all(row == row[1])))]
                                        constant.genes <- unique(c(constant.genes, rownames(newx)[unlist(apply(newx, 1, function(row) any(!is.finite(row)) | all(row == row[1])))]))
                                        x <- x[!(rownames(x) %in% constant.genes),]
                                        newx <- newx[!(rownames(newx) %in% constant.genes),]

                                        nfolds <- 5 
                                        ret.list <- list()
                                        ## Apply fit to test/validation data

                                        ## Use glmnet
                                        ## alphas <- c(0, 0.25, 0.5, 0.7, 1)
                                        ## alphas <- c(0, 1)
                                        N <- length(y)
                           set.seed(seed)
                                        foldid <- sample(rep(seq(nfolds), length = N))

                                        svalues <- c("lambda.min", "lambda.1se")
                                        svalues <- c("lambda.min")
                                        for(alpha in alphas) {
                                          for(s in svalues) {
                                            set.seed(seed)
                                            fit <- tryCatch({cv.glmnet(x = t(x), y, family = "gaussian", type.measure = "mse", nfolds = nfolds, foldid = foldid, alpha = alpha, standardize = FALSE, parallel = TRUE)}, error = function(e) { cat("glmnet err\n"); return(NA) })
                                            pred <- tryCatch({predict(fit, newx=t(newx), s = s)}, error = function(e) { cat("glmnet pred err\n"); return(NA) })
                                            ret.list[[length(ret.list)+1]] <- list(train.drug = train.drug, test.drug = test.drug, predicted = pred, actual = newy, fit = fit, alpha = alpha, s = s, model = "glmnet", n.train = n.train, n.test = n.test, iter = iter)
                                          }
                                        }

                                        ## Use RF
                                        if(use.rf) {
                                          set.seed(seed)
                                          fit <- tryCatch({randomForest(x = t(x), y, importance = TRUE, keep.forest = TRUE)}, error = function(e) { cat("rf error\n"); return(NA) })
                                          pred <- tryCatch({predict(fit, t(newx))}, error = function(e) { cat("rf pred error\n"); return(NA) })
                                          ret.list[[length(ret.list)+1]] <- list(train.drug = train.drug, test.drug = test.drug, predicted = pred, actual = newy, fit = fit, alpha = NA, s = NA, model = "rf", n.train = n.train, n.test = n.test, iter = iter)
                                        }

                                        ## Use SVM
                                        if(use.svm) {
                                          set.seed(seed)
                                          fit <- tryCatch({svm(x = t(x), y, scale = FALSE)}, error = function(e) { cat("svm error\n"); return(NA) })
                                          pred <- tryCatch({predict(fit, t(newx))}, error = function(e) { cat("svm pred error\n"); return(NA) })
                                          ret.list[[length(ret.list)+1]] <- list(train.drug = train.drug, test.drug = test.drug, predicted = pred, actual = newy, fit = fit, alpha = NA, s = NA, model = "svm", n.train = n.train, n.test = n.test, iter = iter)
                                        }

                                        ## Return a baseline in which the prediction is just the mean of training data set
                                        len <- length(newy) 
                                        pred <- rep(mean(y, na.rm = TRUE), len)
                                        ret.list[[length(ret.list)+1]] <- list(train.drug = train.drug, test.drug = test.drug, predicted = pred, actual = newy, fit = NA, alpha = NA, s = NA, model = "mean", n.train = n.train, n.test = n.test, iter = iter)

                                        ret.list
                               })
                })
   ret
}

annotate.ranked.list <- function(lst) {
  foo <- llply(lst, .fun = function(df) {
    df.10 <- df[!grepl(pattern="Intercept", df$gene),]
    df.10 <- df.10[order(as.numeric(df.10$grand.rank), decreasing=FALSE),]
    df.10 <- df.10[1:10,]
    top.10.sym <- ensg.to.sym.mapping(df.10$gene)
    df.10 <- merge(df.10, top.10.sym, by.x = "gene", by.y = "ensg", all = TRUE)
    df.10 <- df.10[order(as.numeric(df.10$grand.rank), decreasing=FALSE),]
    ## top.10.rank.str <- paste(df.10$grand.rank, collapse=", ")
    top.10.rank.str <- paste(unlist(apply(df.10[!is.na(df.10$symbol),c("symbol", "grand.rank")], 1, function(row) paste0(row[1], " (", as.numeric(row[2]), ")"))), collapse=", ")
    top.10.ensg <- sort(as.character(df.10$gene))
    top.10.sym <- top.10.sym$symbol
    top.10.sym <- top.10.sym[!is.na(top.10.sym)]
    top.10.sym.str <- paste(sort(as.character(top.10.sym)), collapse=", ")
    top.10.ensg.str <- paste(top.10.ensg, collapse=", ")
    vec <- c(top.10.rank.str, top.10.ensg.str, top.10.sym.str)
    names(vec) <- c("ranks", "genes.ensg", "genes.symbol")
    vec
  })
  tbl <- as.data.frame(do.call("rbind", foo))
  tbl$drug <- names(foo)
  tbl <- tbl[, c("drug", "genes.symbol", "ranks", "genes.ensg")]
  tbl
}

library(grid)
library(gridExtra)
library(gtable)

plot.table <- function(tbl, main = "", suppress.rows = TRUE) {
  t1 <- NULL
  if(suppress.rows == TRUE) {
    t1 <- tableGrob(tbl, rows=NULL)
  } else {
    t1 <- tableGrob(tbl)
  }
  title <- textGrob(main,gp=gpar(fontsize=15))
  padding <- unit(5,"mm")

  table <- gtable_add_rows(
       t1, 
       heights = grobHeight(title) + padding,
       pos = 0)
  table <- gtable_add_grob(
      table, 
      title, 
      1, 1, 1, ncol(table))

  grid.newpage()
  grid.draw(table)
}


annotate.with.overlapping.gene.sets <- function(tbl, gene.sets, gene.set.col.name, tbl.gene.col, tbl.target.gene.col) {
  tbl[,gene.set.col.name] <- ""
  for(i in 1:nrow(tbl)) {
    features <- as.character(tbl[i, tbl.gene.col])
    target.genes <- as.character(tbl[i, tbl.target.gene.col])
    if(is.na(target.genes) || is.null(target.genes)) { next }
    target.genes <- unlist(strsplit(target.genes, split=",[ ]*"))
    genes.in.overlapping.pathways <- c()
    for(indx in 1:length(gene.sets$geneset.names)) {
      genes.in.pathway <- gene.sets$genesets[[indx]]
      genes.in.pathway <- na.omit(genes.in.pathway)
      if(any(genes.in.pathway %in% target.genes)) {
        genes.in.overlapping.pathways <- unique(c(genes.in.overlapping.pathways, genes.in.pathway))
      }
    }
    if(length(genes.in.overlapping.pathways) > 0) {
      gs <- genes.in.overlapping.pathways[genes.in.overlapping.pathways %in% features]
      tbl[i, gene.set.col.name] <- paste(gs, collapse=", ")
    }
  }
  tbl
}

plot.correlation.vs.drug <- function(tbl, model = "rf", metric = "spearman") {
  rftbl <- subset(tbl, model == model)
  rftbl <- subset(tbl, metric == metric)
  rftbl <- rftbl[order(rftbl$cor, decreasing=TRUE),]
  ordered.levels <- unique(rftbl$train.drug)
  tbl$model <- factor(tbl$model)
  tbl$train.drug <- factor(tbl$train.drug, levels = ordered.levels)
  tbl$cor <- as.numeric(tbl$cor)
  g <- ggplot(data = tbl, aes(x = train.drug, y = cor, colour = model, group = model))
  g <- g + geom_point() + geom_line()
  g <- g + xlab("Drug")
  g <- g + ylab("Correlation")
  g <- g + theme(axis.text.x=element_text(angle = 45, hjust = 1))
  g
} 

plot.correlation.vs.model <- function(tbl, metric = "spearman") {
  tbl <- subset(tbl, metric == metric)
  tbl$model <- factor(tbl$model)
  tbl$train.drug <- factor(tbl$train.drug)
  tbl$cor <- as.numeric(tbl$cor)
  g <- ggplot(data = tbl, aes(x = model, y = cor, colour = model))
  g <- g + geom_violin()
library(ggbeeswarm)
##  g <- g + geom_point() + geom_boxplot()
g <- g + geom_beeswarm()
  g <- g + xlab("Model")
  g <- g + ylab("Correlation")
  g <- g + theme(axis.text.x=element_text(angle = 45, hjust = 1))
  g
} 

extract.fit <- function(ret, train.drug, alpha, model) {
  for(i in 1:length(ret)) {
    foo <- ret[[i]]
    for(j in 1:length(foo)) {
        df <- foo[[j]]
        flag <- (df$train.drug == train.drug) && (df$model == model)
        if(!is.na(alpha)) { flag <- flag && (df$alpha == alpha) }
        if(flag) { return(df$fit) }
      }
  }
  return(NA)
}

extract.features <- function(ret, s = "lambda.min") {
  tbl <-
     ldply(ret,
             .fun = function(foo) {
                ldply(foo,
             .fun = function(df) {
                      ret <- data.frame(train.drug = df$train.drug, alpha = df$alpha, model = df$model, n.train = df$n.train, features = "")
                      if((df$model != "glmnet")) { return(ret) }
                      if(as.numeric(df$alpha) != 1) { return(ret) }
                      if(is.null(df$fit) || is.na(df$fit)) { return(ret) }
                      coeffs <- NULL
                      switch(df$model,
                        "glmnet" = {
                          coeffs <- coefficients(df$fit, s = s)
                        },
                        "rf" = {
                          col <- colnames(importance(df$fit))[1]
                          ## NB: %IncMSE = [ mse(j) - mse0 ] / mse0 * 100  [ mse(j): permuted; mse0: original MSE ]
                          ## i.e., larger/more positive %IncMSE is better
                          if(!grepl(col, pattern="IncMSE")) { stop(paste0("Was expected ", col, " to be %IncMSE\n")) }
                          coeffs <- importance(df$fit)[,1,drop=FALSE]
                        },
                        { stop(paste0("Unknown model ", df$model, "\n")) })
                                         ns <- rownames(coeffs)
                                         coeffs <- as.vector(coeffs)
                                         names(coeffs) <- ns
                                         coeffs <- coeffs[!grepl(pattern="Intercept", names(coeffs))]
                                         coeffs <- coeffs[coeffs != 0]
                                         n.features <- length(coeffs) 
                                         coeffs <- coeffs[order(names(coeffs))]
                                         coeff.vals <- paste(coeffs, collapse=",")
                                         coeffs <- paste(names(coeffs), collapse=",")

                      ret$features <- coeffs
                      ret
                    }) })
  tbl$model <- as.character(tbl$model)
  tbl[(tbl$model == "glmnet") & (tbl$alpha == 0),"model"] <- "ridge"
  tbl[(tbl$model == "glmnet") & (tbl$alpha == 1),"model"] <- "lasso"
  tbl
}

extract.all.features <- function(ret, s = "lambda.min") {
  tbl <-
     ldply(ret,
             .fun = function(foo) {
                ldply(foo,
             .fun = function(df) {
                      ret <- data.frame(train.drug = df$train.drug, alpha = df$alpha, model = df$model, n.train = df$n.train, features = "")
                      if(is.null(df$fit) || is.na(df$fit)) { return(ret) }
                      coeffs <- NULL
                      switch(df$model,
                        "glmnet" = {
                          coeffs <- coefficients(df$fit, s = s)
                        },
                        "rf" = {
                          col <- colnames(importance(df$fit))[1]
                          ## NB: %IncMSE = [ mse(j) - mse0 ] / mse0 * 100  [ mse(j): permuted; mse0: original MSE ]
                          ## i.e., larger/more positive %IncMSE is better
                          if(!grepl(col, pattern="IncMSE")) { stop(paste0("Was expected ", col, " to be %IncMSE\n")) }
                          coeffs <- importance(df$fit)[,1,drop=FALSE]
                        },
                        { stop(paste0("Unknown model ", df$model, "\n")) })
                                         ns <- rownames(coeffs)
                                         coeffs <- as.vector(coeffs)
                                         names(coeffs) <- ns
                                         coeffs <- coeffs[!grepl(pattern="Intercept", names(coeffs))]
                                         coeffs <- coeffs[coeffs != 0]
                                         n.features <- length(coeffs) 
                                         coeffs <- coeffs[order(names(coeffs))]
                                         coeff.vals <- paste(coeffs, collapse=",")
                                         coeffs <- paste(names(coeffs), collapse=",")

                      ret$features <- coeffs
                      ret
                    }) })
  tbl$model <- as.character(tbl$model)
  tbl[(tbl$model == "glmnet") & (tbl$alpha == 0),"model"] <- "ridge"
  tbl[(tbl$model == "glmnet") & (tbl$alpha == 1),"model"] <- "lasso"
  tbl
}

extract.predicted.actual <- function(ret, train.drug, alpha, s, model) {
  tbl <-
     ldply(ret,
             .fun = function(foo) {
                ldply(foo,
             .fun = function(df) {
                ret <- NULL
                if(!is.na(df$predicted) && (class(df$predicted) != "logical") && (df$n.test > 0)) {
                  if(df$model != "rf") {
                    ret <- data.frame(train.drug = df$train.drug, test.drug = df$test.drug, alpha = df$alpha, s = df$s, model = df$model, actual = df$actual, predicted = df$predicted[,1], sample.names = rownames(df$predicted))
                  } else {
                    ret <- data.frame(train.drug = df$train.drug, test.drug = df$test.drug, alpha = NA, s = NA, model = df$model, actual = as.vector(df$actual), predicted = as.vector(unname(df$predicted)), sample.names = as.vector(names(df$predicted)))
                  }
                }
                ret
             }) })
  ## tbl <- tbl[(tbl$train.drug == train.drug) & (tbl$alpha == alpha) & (tbl$s == s) & (tbl$model == model),]
  tbl <- tbl[(tbl$train.drug == train.drug) & (tbl$model == model),]
  if(!is.na(alpha)) {
    tbl <- tbl[tbl$alpha == alpha, ]
  }
  if(!is.na(s)) {
    tbl <- tbl[tbl$s == s,]
  }
  list(actual = tbl$actual, predicted = tbl$predicted, sample.names = as.character(tbl$sample.names))
}

extract.results <- function(ret) {
  tbl <-
     ldply(ret,
             .fun = function(foo) {
                ldply(foo,
             .fun = function(df) {
                ret <- NULL
                if(!is.na(df$predicted) && (class(df$predicted) != "logical") && (df$n.test > 0)) {
                  method <- "pearson"
##                  ct <- cor.test(df$predicted, df$actual, method = method)
                  ret <- tryCatch({ ct <- cor.test(df$predicted, df$actual, method = method); data.frame(val = ct$estimate, pval = ct$p.value, train.drug = df$train.drug, test.drug = df$test.drug, alpha = df$alpha, s = df$s, model = df$model, metric = method, n.test = df$n.test) }, error = function(e) { data.frame(val = NA, pval = NA, train.drug = df$train.drug, test.drug = df$test.drug, alpha = df$alpha, s = df$s, model = df$model, metric = "pearson", n.test = df$n.test) })
                  method <- "spearman"
##                  ct <- cor.test(df$predicted, df$actual, method = method)
##                  ret2 <- data.frame(val = ct$estimate, pval = ct$p.value, train.drug = df$train.drug, test.drug = df$test.drug, alpha = df$alpha, s = df$s, model = df$model, metric = method, n.test = df$n.test)
                  ret2 <- tryCatch({ ct <- cor.test(df$predicted, df$actual, method = method); data.frame(val = ct$estimate, pval = ct$p.value, train.drug = df$train.drug, test.drug = df$test.drug, alpha = df$alpha, s = df$s, model = df$model, metric = method, n.test = df$n.test) }, error = function(e) { data.frame(val = NA, pval = NA, train.drug = df$train.drug, test.drug = df$test.drug, alpha = df$alpha, s = df$s, model = df$model, metric = "spearman", n.test = df$n.test) })
                  mse <- mean((df$predicted-df$actual)^2)
                  ret3 <- data.frame(val = mse, pval = NA, train.drug = df$train.drug, test.drug = df$test.drug, alpha = df$alpha, s = df$s, model = df$model, metric = "mse", n.test = df$n.test)
                  ret <- rbind(ret, ret2, ret3)
                  
                } else {
                  ret <- data.frame(val = NA, pval = NA, train.drug = df$train.drug, test.drug = df$test.drug, alpha = df$alpha, s = df$s, model = df$model, metric = "pearson", n.test = df$n.test)
                  ret2 <- data.frame(val = NA, pval = NA, train.drug = df$train.drug, test.drug = df$test.drug, alpha = df$alpha, s = df$s, model = df$model, metric = "spearman", n.test = df$n.test)
                  ret3 <- data.frame(val = NA, pval = NA, train.drug = df$train.drug, test.drug = df$test.drug, alpha = df$alpha, s = df$s, model = df$model, metric = "mse", n.test = df$n.test)
                  ret <- rbind(ret, ret2, ret3)
                }
                ret
             }) })
  t2 <- subset(tbl, is.na(s) | (s == "lambda.min"))
  t2$model <- as.character(t2$model)
  t2[(t2$model == "glmnet") & (t2$alpha == 0),"model"] <- "ridge"
  t2[(t2$model == "glmnet") & (t2$alpha == 1),"model"] <- "lasso"
  t2 <- subset(t2, model != "glmnet")
  t2$train.drug <- as.character(t2$train.drug)
  t2$test.drug <- as.character(t2$test.drug)
  t2$s <- as.character(t2$s)
  t2$model <- as.character(t2$model)
  t2$metric <- as.character(t2$metric)

  t2
}

extract.training.info <- function(ret, flatten.downsampled.fits = FALSE) {
  tbl <-
     ldply(ret,
             .fun = function(foo) {
                if(flatten.downsampled.fits == TRUE) { foo <- unlist(foo, recursive = FALSE) }
                ldply(foo,
             .fun = function(df) {
                  ret <- data.frame(train.drug = df$train.drug, alpha = df$alpha, model = df$model, n.train = df$n.train)
                ret
             })  })
  t2 <- tbl
  t2$model <- as.character(t2$model)
  t2[(t2$model == "glmnet") & (t2$alpha == 0),"model"] <- "ridge"
  t2[(t2$model == "glmnet") & (t2$alpha == 1),"model"] <- "lasso"
  t2 <- subset(t2, model != "glmnet")
  t2$train.drug <- as.character(t2$train.drug)
  t2$model <- as.character(t2$model)

  t2
}

extract.nested.results <- function(ret) {
  tbl <-
     ldply(ret,
             .fun = function(bar) {
                ldply(bar,
             .fun = function(foo) {
                ldply(foo,
             .fun = function(df) {
                ret <- NULL

                if((class(df$predicted) != "logical") && !is.na(df$predicted)) {
                  method <- "pearson"
                  ## ct <- cor.test(df$predicted, df$actual, method = method)
                  ## ret <- data.frame(val = ct$estimate, pval = ct$p.value, train.drug = df$train.drug, test.drug = df$test.drug, alpha = df$alpha, s = df$s, model = df$model, metric = method, n.train = df$n.train, n.test = df$n.test)
                  ret <- tryCatch({ ct <- cor.test(df$predicted, df$actual, method = method); data.frame(val = ct$estimate, pval = ct$p.value, train.drug = df$train.drug, test.drug = df$test.drug, alpha = df$alpha, s = df$s, model = df$model, metric = method, n.test = df$n.test) }, error = function(e) { data.frame(val = NA, pval = NA, train.drug = df$train.drug, test.drug = df$test.drug, alpha = df$alpha, s = df$s, model = df$model, metric = "pearson", n.train = df$n.train, n.test = df$n.test) })
                  method <- "spearman"
                  ## ct <- cor.test(df$predicted, df$actual, method = method)
                  ## ret2 <- data.frame(val = ct$estimate, pval = ct$p.value, train.drug = df$train.drug, test.drug = df$test.drug, alpha = df$alpha, s = df$s, model = df$model, metric = method, n.train = df$n.train, n.test = df$n.test)
                  ret2 <- tryCatch({ ct <- cor.test(df$predicted, df$actual, method = method); data.frame(val = ct$estimate, pval = ct$p.value, train.drug = df$train.drug, test.drug = df$test.drug, alpha = df$alpha, s = df$s, model = df$model, metric = method, n.test = df$n.test) }, error = function(e) { data.frame(val = NA, pval = NA, train.drug = df$train.drug, test.drug = df$test.drug, alpha = df$alpha, s = df$s, model = df$model, metric = "spearman", n.train = df$n.train, n.test = df$n.test) })
                  mse <- mean((df$predicted-df$actual)^2)
                  ret3 <- data.frame(val = mse, pval = NA, train.drug = df$train.drug, test.drug = df$test.drug, alpha = df$alpha, s = df$s, model = df$model, metric = "mse", n.train = df$n.train, n.test = df$n.test)
                  ret <- rbind(ret, ret2, ret3)
                  
                } else {
                  ret <- data.frame(val = NA, pval = NA, train.drug = df$train.drug, test.drug = df$test.drug, alpha = df$alpha, s = df$s, model = df$model, metric = "pearson", n.train = df$n.train, n.test = df$n.test)
                  ret2 <- data.frame(val = NA, pval = NA, train.drug = df$train.drug, test.drug = df$test.drug, alpha = df$alpha, s = df$s, model = df$model, metric = "spearman", n.train = df$n.train, n.test = df$n.test)
                  ret3 <- data.frame(val = NA, pval = NA, train.drug = df$train.drug, test.drug = df$test.drug, alpha = df$alpha, s = df$s, model = df$model, metric = "mse", n.train = df$n.train, n.test = df$n.test)
                  ret <- rbind(ret, ret2, ret3)
                }
                ret
             }) }) })
  t2 <- subset(tbl, is.na(s) | (s == "lambda.min"))
  t2$model <- as.character(t2$model)
  t2[(t2$model == "glmnet") & (t2$alpha == 0),"model"] <- "ridge"
  t2[(t2$model == "glmnet") & (t2$alpha == 1),"model"] <- "lasso"
  t2 <- subset(t2, model != "glmnet")
  t2$train.drug <- as.character(t2$train.drug)
  t2$test.drug <- as.character(t2$test.drug)
  t2$s <- as.character(t2$s)
  t2$model <- as.character(t2$model)
  t2$metric <- as.character(t2$metric)

  t2
}

plot.trend <- function(df) {
  g <- ggplot(df, aes(x = train.set, y = fisher))
  g <- g + geom_point()
  g <- g + geom_errorbar(aes(ymin = fisher.lb, ymax = fisher.ub))
  g
}

plot.all.trends <- function(results) {
  trend.res <- list()
  library(psych)
  for(test.set in test.sets) {
    tmp <- results[[test.set]]
    tmp <- tmp[tmp$model != "mean",] 
    tmp <- tmp[tmp$metric != "mse",]
    tmp$fisher <- fisherz(tmp$val)
    tmp$fisher.lb <- unlist(apply(tmp[, c("val", "n.test")], 1, function(row) row[1] - 1/(sqrt(row[2] - 3)))) 
    tmp$fisher.ub <- unlist(apply(tmp[, c("val", "n.test")], 1, function(row) row[1] + 1/(sqrt(row[2] - 3)))) 
    ## Create an ordered train.set factor so we can sort by it below
    tmp$train.set <- factor(tmp$train.set, levels = c("ctrp.all", "ctrp.heme", "ctrp.aml", "ohsu.train"))
    tbl <- ddply(tmp, .variables = c("alpha", "s", "model"),
                 .fun = function(df1) {
                          ddply(df1, .variables = c("metric"),
                                .fun = function(df2) {
                                         ddply(df2, .variables = c("test.drug"),
                                               .fun = function(df3) {
                                                        df3$trend <- NA
                                                        df4 <- df3[df3$train.set != "ohsu.train",]
                                                        df4 <- df4[, c("fisher", "fisher.lb", "fisher.ub")]
                                                        df4 <- na.omit(df4)
                                                        if(nrow(df4) > 1) {
                                                          all.overlap <- TRUE
                                                          mono.increasing <- TRUE
                                                          for(i in 1:(nrow(df4)-1)) {
                                                            if(df4$fisher.ub[i] >= df4$fisher.lb[i+1]) {
                                                              mono.increasing <- FALSE
                                                            }
                                                            for(j in (i+1):nrow(df4)) { 
                                                              if(df4$fisher.lb[i] > df4$fisher.ub[j]) {
                                                                all.overlap <- FALSE
                                                              }
                                                              if(df4$fisher.ub[i] < df4$fisher.lb[j]) {
                                                                all.overlap <- FALSE
                                                              }
                                                            }
                                                          }
                                                          if(all.overlap && mono.increasing) {
                                                            stop("Was not expecting all overlapping and increasing\n")
                                                          }
                                                          if(mono.increasing) { df3$trend <- "mono" }
                                                          if(all.overlap) { df3$trend <- "overlap" }
                                                          p <- ggplot(df3, aes(x = train.set, y = fisher))
                                                          p <- p + geom_point()
                                                          p <- p + geom_errorbar(aes(ymin = fisher.lb, ymax = fisher.ub))
                                                          pdf(paste0(df3$test.drug[1], "-", test.set, "-", df3$model, "-", df3$metric, "-trend.pdf"))
                                                          print(p)
                                                          dev.off()
                                                        }
                                                        df3
                                                      })
                                })
                        })
  
    trend.res[[test.set]] <- tbl
  }
  cat("Done plotting trends\n")
  trend.res
}

## Even in cases in which we repeatedly downsample the same size as the data set (i.e., the training data set
## does not change), we might have multiple, different test results because of the randomness of glmnet
## in choosing lambda.1se/lambda.min.  Yes, we use a seed, but that seed will be iteration dependent--
## since each downsample has a different iteration, it will give a different lambda.min, even if the data
## set does not change.
plot.all.downsampled.trends <- function(results) {
  trend.res <- list()
  library(psych)
  sets <- test.sets
##  sets <- c("fimm")
  for(test.set in sets) {
    tmp <- results[[test.set]]
    tmp <- tmp[tmp$model != "mean",] 
    tmp <- tmp[tmp$metric != "mse",]
    tmp$fisher <- fisherz(tmp$val)
## tmp <- subset(tmp, test.drug == tmp$test.drug[1])
## tmp <- subset(tmp, test.drug == "FIMM003792")
    ## Create an ordered train.set factor so we can sort by it below
    tmp$train.set <- factor(tmp$train.set, levels = c("ctrp.all", "ctrp.heme", "ctrp.aml", "ohsu.train"))
    tbl <- ddply(tmp, .variables = c("alpha", "s", "model"),
                 .fun = function(df1) {
                          ddply(df1, .variables = c("metric"), .parallel = FALSE,
                                .fun = function(df2) {
                                         ddply(df2, .variables = c("inhibitor"), .parallel = FALSE,
                                               .fun = function(df3) {
                                                        df3$trend <- NA
                                                        df4 <- df3[df3$train.set != "ohsu.train",]

                                                        erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
                                                        min.prob <- 0.5 - erf(1/sqrt(2))/2
                                                        max.prob <- 0.5 + erf(1/sqrt(2))/2
                                                        ## df4 <- ddply(df4, .variables = c("train.set"), .fun = function(df) { c(min(df$fisher, na.rm=TRUE), max(df$fisher, na.rm=TRUE)) })
                                                        df4 <- ddply(df4, .variables = c("train.set"), .fun = function(df) { c(quantile(df$fisher, probs = min.prob, na.rm=TRUE), quantile(df$fisher, probs = max.prob, na.rm=TRUE)) })
                                                        colnames(df4) <- c("train.set", "fisher.lb", "fisher.ub")
                                                        df4 <- df4[, c("fisher.lb", "fisher.ub")]
                                                        df4 <- na.omit(df4)
                                                        if(nrow(df4) > 1) {
                                                          all.overlap <- TRUE
                                                          mono.increasing <- TRUE
                                                          for(i in 1:(nrow(df4)-1)) {
                                                            if(df4$fisher.ub[i] >= df4$fisher.lb[i+1]) {
                                                              mono.increasing <- FALSE
                                                            }
                                                            for(j in (i+1):nrow(df4)) { 
                                                              if(df4$fisher.lb[i] > df4$fisher.ub[j]) {
                                                                all.overlap <- FALSE
                                                              }
                                                              if(df4$fisher.ub[i] < df4$fisher.lb[j]) {
                                                                all.overlap <- FALSE
                                                              }
                                                            }
                                                          }
                                                          if(all.overlap && mono.increasing) {
                                                            stop("Was not expecting all overlapping and increasing\n")
                                                          }
                                                          if(mono.increasing) { df3$trend <- "mono" }
                                                          if(all.overlap) { df3$trend <- "overlap" }
                                                          if(is.na(df3$trend)) { df3$trend <- "none" }
                                                          if(!is.na(df3$trend)) {
                                                            p <- ggplot(df3, aes(x = train.set, y = fisher))
                                                            p <- p + geom_violin()
                                                            p <- p + geom_point()
                                                            p <- p + ggtitle(paste0(df3$inhibitor[1], " test set:", test.set, " ", df3$model[1], " ", df3$metric[1]))
                                                            ## p <- p + geom_errorbar(aes(ymin = fisher.lb, ymax = fisher.ub))
                                                            file <- paste0(df3$inhibitor[1], "-", test.set, "-", df3$model[1], "-", df3$metric[1], "-trend-", df3$trend[1], ".pdf")
                                                            pdf(file)
                                                            print(p)
                                                            dev.off()
                                                          }
                                                        }
                                                        df3
                                                      })
                                })
                        })
  
    trend.res[[test.set]] <- tbl
  }
  cat("Done plotting trends\n")
  trend.res
}

