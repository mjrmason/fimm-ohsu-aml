## For each drug in common between OHSU and FIMM, fit elastic net to response (auc)
## Specifically:
## 1. Split into 70%/30% training/test
## 2. Standardize training/test sets separately for X (gene expr) and y (drug response)
## 3. Optimize elastic net parameters on standardized training set with cv.glmnet
## 4. Run prediction on test

library(synapseClient)
library(data.table)
library(ggplot2)
library(gridExtra)

suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(caret))
suppressPackageStartupMessages(library(glmnet))
library( ReporteRs )
library(scales)


library("scales")
reverselog_trans <- function(base = exp(1)) {
    trans <- function(x) -log(x, base)
    inv <- function(x) base^(-x)
    trans_new(paste0("reverselog-", format(base)), trans, inv, 
              log_breaks(base = base), 
              domain = c(1e-100, Inf))
}

plot.stacked.pval.correlation.figures <- function(tbl, pval.col, correlation.col, x.col, y.cor.lab, text.size = NULL, main = NULL) {
    library(gridExtra)
    library(grid)
    g1 <- ggplot(data = tbl)
    ## g1 <- g1 + geom_point(aes(x = gene.set, y = -log10(pval), col = Direction))
    g1 <- g1 + geom_point(aes_string(x = x.col, y = pval.col))
##    g1 <- g1 + geom_violin(aes_string(x = x.col, y = pval.col))
    g1 <- g1 + scale_y_continuous(trans=reverselog_trans(10))
    g1 <- g1 + theme(axis.title.x=element_blank(), axis.text.x=element_text(angle = 45, hjust = 1), axis.ticks.x=element_blank())
    if(!is.null(text.size)) {
        g1 <- g1 + theme(axis.text.x=element_text(angle = 45, hjust = 1, size = text.size))
    }
    g1 <- g1 + geom_hline(yintercept = 0.01)
    g1 <- g1 + ylab(bquote(italic('p')-value))
    if(!is.null(main)) {
        g1 <- g1 + ggtitle(main) + theme(plot.title = element_text(hjust = 0.5))
    }
    ## Grab the axis labels
    gxlabels <- gtable::gtable_filter(ggplotGrob(g1), "axis-b")
    ## Then remove it
    g1 <- g1 + theme(axis.title.x=element_blank(), axis.text.x=element_blank())
    ## gr1 <- ggplotGrob(g1)
    
    g2 <- ggplot(data = tbl)
    ## g2 <- g2 + geom_point(aes(x = gene.set, y = -log10(pval), col = Direction))
    g2 <- g2 + geom_point(aes_string(x = x.col, y = correlation.col))
##    g2 <- g2 + geom_violin(aes_string(x = x.col, y = correlation.col))
    g2 <- g2 + theme(axis.title.x=element_blank(), axis.text.x=element_text(angle = 45, hjust = 1), axis.ticks.x=element_blank())
    if(!is.null(text.size)) {
        g2 <- g2 + theme(axis.text.x=element_text(angle = 45, hjust = 1, size = text.size))
    }
##    g2 <- g2 + geom_hline(yintercept = 0.01)
    g2 <- g2 + ylab(y.cor.lab)
##    g2 <- g2 + ggtitle(main2) + theme(plot.title = element_text(hjust = 0.5))

    gr1 <- ggplot_gtable(ggplot_build(g1))
    gr2 <- ggplot_gtable(ggplot_build(g2))
  
    ## Set the heights to be the same for both plots
    gr1$heights[which(gr1$layout$name == "panel")] <- unit(5, units="cm")
    gr2$heights[which(gr2$layout$name == "panel")] <- unit(5, units="cm")    

##    gr1$widths[which(gr1$layout$name == "axis-l")] <- unit(2, units="cm")
##    gr2$widths[which(gr2$layout$name == "axis-l")] <- unit(2, units="cm")    

    maxWidth = unit.pmax(gr1$widths[2:3], gr2$widths[2:3])

    gr1$widths[2:3] <- maxWidth
    gr2$widths[2:3] <- maxWidth

    ## grid.arrange(gr1, gr2)
    arrangeGrob(grobs = list(gr1, gr2))
}
    
fit.elastic.net <- function(X, y) {
  N <- nrow(X)
  nfolds <- 5
  foldid <- sample(rep(seq(nfolds), length = N))
  cv.list <- list()
  alphas <- seq(from = 0, to = 1, by = 0.05)
  for(i in 1:length(alphas)) {
     cv <- cv.glmnet(X, y, family = "gaussian", type.measure = "mse", foldid = foldid, nfolds = nfolds, alpha = alphas[i], standardize = FALSE, parallel = TRUE)
    cv.list[[i]] <- cv
  }
  cv.list
}

fit.elastic.net.alpha <- function(alpha, x, y, foldid, nfolds, seed = 1234) {
  set.seed(seed)
  N <- nrow(x)
  cv <- cv.glmnet(x, y, family = "gaussian", type.measure = "mse", foldid = foldid, nfolds = nfolds, alpha = alpha, standardize = FALSE, parallel = FALSE)
  cv
}

prepare.ohsu.drug.response.and.expr.matrices <- function(drc.df, expr.df, rnaseq.summary.df, drugs, response.col) {

  drc.df <- subset(drc.df, inhibitor %in% drugs)

  ## Convert lab_id as used in drug response to SeqIDs used in expr matrix.
  drc.df <- merge(drc.df, rnaseq.summary.df, by.x = "lab_id", by.y = "Original_LabID")

  ## Exclude any NAs in response
  flag <- !is.na(drc.df[, response.col])
  drc.df <- drc.df[flag, ]

  ## Take the median of replicates in the X and y data
  drc.df <- ddply(drc.df, c("inhibitor", "SeqID"), .fun = function(df) median(df[, response.col]))
  colnames(drc.df) <- c("inhibitor", "SeqID", response.col)

  ## Convert drug response table into a matrix in which columns are samples and
  ## rows are drug response
  drc.df <- spread_(drc.df, "SeqID", response.col)  
  rownames(drc.df) <- drc.df$inhibitor
  drc.df <- drc.df[, !(colnames(drc.df) == "inhibitor")]

  ## Restrict X and response to be over the same samples
  common.samples <- intersect(colnames(drc.df), colnames(expr.df))
  drc.df <- drc.df[, common.samples]
  expr.df <- expr.df[, common.samples]

  ## Restrict to drugs of interest
  drc.df <- drc.df[rownames(drc.df) %in% drugs, ]

  return(list(drc.df = drc.df, expr.df = expr.df))
}

prepare.fimm.drug.response.and.expr.matrices <- function(drc.df, expr.df, drugs, response.col) {

  drug.col <- "DRUG_ID"
  patient.col <- "SCREEN_ID"

  flag <- drc.df[, drug.col] %in% drugs
  drc.df <- drc.df[flag, ]

  ## Exclude any NAs in response
  flag <- !is.na(drc.df[, response.col])
  drc.df <- drc.df[flag, ]

  ## Take the median of replicates in the X and y data
  drc.df <- drc.df[, c(drug.col, patient.col, response.col)]

  ## Convert drug response table into a matrix in which columns are samples and
  ## rows are drug response
  drc.df <- spread_(drc.df, patient.col, response.col)  
  rownames(drc.df) <- drc.df[, drug.col]
  drc.df <- drc.df[, !(colnames(drc.df) == drug.col)]

  ## Restrict X and response to be over the same samples
  common.samples <- intersect(colnames(drc.df), colnames(expr.df))
  drc.df <- drc.df[, common.samples]
  expr.df <- expr.df[, common.samples]

  ## Restrict to drugs of interest
  drc.df <- drc.df[rownames(drc.df) %in% drugs, ]

  return(list(drc.df = drc.df, expr.df = expr.df))
}

plot.predicted.vs.actual.response <- function(df, target.alpha, target.drug) {
  df <- subset(df, (alpha == target.alpha) & (drug == target.drug))
  df$response <- as.numeric(df$response)
  predicted.flag <- df$response.type == "predicted"
  predicted.df <- df[predicted.flag, c("sample", "response", "alpha", "drug")]
  colnames(predicted.df) <- c("sample", "predicted.response", "alpha", "drug")
  actual.df <- unique(df[!predicted.flag, c("sample", "response", "alpha", "drug")])
  colnames(actual.df) <- c("sample", "actual.response", "alpha", "drug")
  m <- merge(predicted.df, actual.df)
  m$predicted.response <- as.numeric(m$predicted.response)
  m$actual.response <- as.numeric(m$actual.response)
  
  ct.spearman <- cor.test(m$actual.response, m$predicted.response, method="spearman")
  ct.pearson <- cor.test(m$actual.response, m$predicted.response, method="pearson")
  cat("Actual vs all predicted: spearman corr = ", unname(ct.spearman$estimate), " pval = ", unname(ct.spearman$p.value), "\n")
  cat("Actual vs all predicted: pearson corr = ", unname(ct.pearson$estimate), " pval = ", unname(ct.pearson$p.value), "\n")

  median.pred <- ddply(df[predicted.flag,], .variables = "sample",
                       .fun = function(df.sample) {
                         vec <- c(predicted.response = median(as.numeric(df.sample$response)), sample = df.sample$sample[1])
                         vec
                       })
  actual.df <- unique(df[!predicted.flag, c("sample", "response", "alpha", "drug")])
  colnames(actual.df) <- c("sample", "actual.response", "alpha", "drug")
  m <- merge(median.pred, actual.df)
  m$predicted.response <- as.numeric(m$predicted.response)
  m$actual.response <- as.numeric(m$actual.response)
  ct.spearman <- cor.test(m$actual.response, m$predicted.response, method="spearman")
  ct.pearson <- cor.test(m$actual.response, m$predicted.response, method="pearson")
  cat("Actual vs median predicted: spearman corr = ", unname(ct.spearman$estimate), " pval = ", unname(ct.spearman$p.value), "\n")
  cat("Actual vs median predicted: pearson corr = ", unname(ct.pearson$estimate), " pval = ", unname(ct.pearson$p.value), "\n")
  
  print(head(m))
  g <- ggplot(m, aes(x = actual.response, y = predicted.response))
  g <- g + geom_point()
  
  ## g <- g + geom_violin(data = m, aes(x = factor(actual.response), y = predicted.response))
  fit <- lm(data = m, predicted.response ~ actual.response)
  r2 <- round(summary(fit)$r.squared, 2)
  rc <- cor(x = m$actual.response, y = m$predicted.response, method = "pearson")
  ##      eqn <- bquote(r^2 == .(r2) * "," ~~ RMSE == .(rmse))
  eqn <- r2
  eqn <- substitute(italic(r)^2~"="~r2, 
                    list(r2 = format(summary(fit)$r.squared, digits = 3)))
  eqn <- substitute(italic(r)~"="~rc, 
                    list(rc = format(rc, digits = 3)))
  eqn <- as.character(as.expression(eqn));
  d.eqn <- data.frame(label = eqn, x = min(m$actual.response,na.rm=TRUE) + 0.1 * (max(m$actual.response, na.rm=TRUE) - min(m$actual.response, na.rm=TRUE)), y = max(m$predicted.response,na.rm=TRUE) - 0.1 * (max(m$predicted.response, na.rm=TRUE) - min(m$predicted.response, na.rm=TRUE)))
  g <- g + geom_text(data = d.eqn, aes(x = x, y = y, label = label), parse=TRUE)
  
  g <- g + geom_smooth(method = "lm", se = FALSE)
  
  g <- g + ggtitle(paste0("alpha = ", target.alpha, " drug = ", target.drug))
  cat("Returning result\n")
  print(g)
  g
}


remove.outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}


bootstrap.elastic.net.fimm <- function(drc.df, expr.df, drugs, response.col, alphas, num.bootstraps = 100, seed = 1234, remove.outlier = FALSE) {

  lst <- prepare.fimm.drug.response.and.expr.matrices(drc.df, expr.df, drugs, response.col)
  drc.df <- lst[["drc.df"]]
  expr.df <- lst[["expr.df"]]

  if(remove.outlier) {
    for(i in 1:nrow(drc.df)) {
      drc.df[i,] <- remove.outliers(drc.df[i, ])
    }
  }

  ## FIMM has relatively few samples--49.  Let's demand at least 35 of them be non-NAs.
  cat(paste0("Original number of drugs: ", nrow(drc.df), "\n"))
  flag <- unlist(apply(drc.df, 1, function(row) length(which(!is.na(row))) >= 35))
  drc.df <- drc.df[flag, ]
  cat(paste0("Number of drugs with >= 35 samples: ", nrow(drc.df), "\n"))
  cat("Done preparing data\n")

  drugs <- rownames(drc.df)

  common.cols <- intersect(colnames(drc.df), colnames(expr.df))
  drc.df <- drc.df[, common.cols]
  expr.df <- expr.df[, common.cols]
  ## For each drug
  res <- llply(drugs, .parallel = TRUE,
               .fun = function(drug) {
                        cat(paste0("Fitting elastic net to ", drug, "\n"))
                        ## For each bootstrap
                        res.drug <- llply(1:num.bootstraps, .parallel = FALSE,
                                          .fun = function(boot.i) {
                                                   set.seed(boot.i)
                                                   
                                                   sample.indices <- which(!is.na(drc.df[drug,]))
                                                   n.samples <- length(sample.indices)
                                                   train.indices <- sample(x = sample.indices, size = n.samples, replace = TRUE)
                                                   test.indices <- sample.indices[!(sample.indices %in% train.indices)] 

                                                   if( (length(train.indices) == 0) && (length(test.indices) == 0) ) {
                                                     null.vec <- rep(NULL, length(alphas))
                                                     fit <- list(alphas = alphas, model = NULL, coeffs = null.vec, predictions = null.vec, test.response = null.vec)
                                                   } else {
                                                     boot.name <- paste0("V", 1:length(train.indices))
                                                     drc.train <- drc.df[drug, train.indices]
                                                     expr.train <- expr.df[, train.indices]

                                                     colnames(drc.train) <- boot.name
                                                     colnames(expr.train) <- boot.name

                                                     drc.test <- drc.df[drug, test.indices]
                                                     expr.test <- expr.df[, test.indices]

                                                     fit <- fit.elastic.net_(drc.train, expr.train, drc.test, expr.test, alphas, train.alphas.in.parallel = FALSE)
                                                   }
                                                   fit[["drug"]] <- drug
                                                   fit[["boot"]] <- boot.i
                                                   fit
                                          })
                        res.drug

               })
  res
}

bootstrap.elastic.net.ohsu <- function(drc.df, expr.df, rnaseq.summary.df, drugs, response.col, alphas, num.bootstraps = 100, seed = 1234) {

  lst <- prepare.ohsu.drug.response.and.expr.matrices(drc.df, expr.df, rnaseq.summary.df, drugs, response.col)
  drc.df <- lst[["drc.df"]]
  expr.df <- lst[["expr.df"]]

  common.cols <- intersect(colnames(drc.df), colnames(expr.df))
  drc.df <- drc.df[, common.cols]
  expr.df <- expr.df[, common.cols]

  drugs <- rownames(drc.df)

  ## For each drug
  res <- llply(drugs, .parallel = TRUE,
               .fun = function(drug) {
                        cat(paste0("Fitting elastic net to ", drug, "\n"))
                        ## For each bootstrap
                        res.drug <- llply(1:num.bootstraps, .parallel = FALSE,
                                          .fun = function(boot.i) {
                                                   set.seed(boot.i)
                                                   
                                                   sample.indices <- which(!is.na(drc.df[drug,]))
                                                   n.samples <- length(sample.indices)
                                                   train.indices <- sample(x = sample.indices, size = n.samples, replace = TRUE)
                                                   test.indices <- sample.indices[!(sample.indices %in% train.indices)] 

                                                   if( (length(train.indices) == 0) && (length(test.indices) == 0) ) {
                                                     null.vec <- rep(NULL, length(alphas))
                                                     fit <- list(alphas = alphas, model = NULL, coeffs = null.vec, predictions = null.vec, test.response = null.vec)
                                                   } else {

                                                     boot.name <- paste0("V", 1:length(train.indices))
                                                     drc.train <- drc.df[drug, train.indices]
                                                     expr.train <- expr.df[, train.indices]
                                                     colnames(drc.train) <- boot.name
                                                     colnames(expr.train) <- boot.name

                                                     drc.test <- drc.df[drug, test.indices]
                                                     expr.test <- expr.df[, test.indices]

                                                     fit <- fit.elastic.net_(drc.train, expr.train, drc.test, expr.test, alphas, train.alphas.in.parallel = FALSE)
                                                   }
                                                   fit[["drug"]] <- drug
                                                   fit[["boot"]] <- boot.i
                                                   fit
                                          })
                        res.drug

               })
  res
}

## alpha = 1: LASSO (shrink coefficients to zero)
## alpha = 0: ridge
fit.elastic.net_ <- function(y.train, X.train, y.test, X.test, alphas, train.alphas.in.parallel = FALSE) {

  ## Independently z-score the training and test data sets
  ## No--don't zscore.  If we do, it will be difficult to plot actual vs predicted responses,
  ## since actual will change as a function of bootstrap (as will predicted, of course)
##                        y.train <- t(scale(t(y.train), center = TRUE, scale = TRUE))
##                        y.test <- t(scale(t(y.test), center = TRUE, scale = TRUE))
  X.train <- scale(X.train, center = TRUE, scale = TRUE)
  X.test <- scale(X.test, center = TRUE, scale = TRUE)
 
  nfolds <- 5
  N <- length(y.train)
  foldid <- sample(rep(seq(nfolds), length = N))
  train.models <- llply(alphas, .parallel = train.alphas.in.parallel, 
                        .fun = function(alpha) { 
                                  fit.elastic.net.alpha(alpha, x = t(X.train), y = as.numeric(y.train), foldid = foldid, nfolds = nfolds)
                        })

  preds <- (lapply(train.models, function(model) predict(model, newx=t(X.test), s="lambda.1se", type="link")))
  ## Return the model, the predictions, and the test response.
  ## Don't return the model -- it takes up too much memory.
  coeffs <- lapply(train.models, function(model) {
                                    df <- data.frame(coef.name = dimnames(coef(model))[[1]], coef.value = matrix(coef(model)))
                                    if(any(df$coef.value > 0)) {
                                      df <- df[df$coef.value > 0,,drop=F]
                                    } else {
                                      df <- NULL
                                    }
                                    df
                                 })
                                    
  rm(train.models)
  train.models <- NULL
  gc()
  vec <- list(alphas = alphas, model = train.models, coeffs = coeffs, predictions = preds, test.response = y.test)
  return(vec)

}

plot.predicted.vs.actual.response.all <- function(df) {
  df$response <- as.numeric(df$response)
  predicted.flag <- df$response.type == "predicted"
  predicted.df <- df[predicted.flag, c("sample", "response", "alpha", "drug")]
  colnames(predicted.df) <- c("sample", "predicted.response", "alpha", "drug")
  actual.df <- unique(df[!predicted.flag, c("sample", "response", "alpha", "drug")])
  colnames(actual.df) <- c("sample", "actual.response", "alpha", "drug")
  m <- merge(predicted.df, actual.df)
  m$predicted.response <- as.numeric(m$predicted.response)
  m$actual.response <- as.numeric(m$actual.response)
  g <- ggplot(data = m, aes(x = actual.response, y = predicted.response))
  g <- g + geom_point()
  g <- g + ggtitle(paste0("Drug: ", df$drug[1], " Alpha: ", df$alpha[1]))
  g
}

post.process.fits <- function(lst, prefix) {
  ## Flatten the predicted responses
  predicted.responses <- ldply(1:length(lst), .parallel = FALSE, 
                               .fun = function(drug.i) {
                                 drug.lst <- lst[[drug.i]]
                                 ret.drug <- ldply(drug.lst, .parallel = FALSE,
                                                   .fun = function(fit) {
                                                     ret.alpha <- ldply(1:length(fit$alphas), .parallel = FALSE,
                                                                        .fun = function(i) {
                                                                          resps <- fit$predictions[[i]][,1]
                                                                          samples <- rownames(fit$predictions[[i]])
                                                                          if(length(resps) <= 1) { return(c(alpha = fit$alphas[i], predicted.response = NA, predicted.samples = NA)) }
                                                                          response <- data.frame(sample = samples, response = resps)
                                                                          resp.str <- paste0(resps, collapse=",")
                                                                          sample.str <- paste0(samples, collapse=",")
                                                                          vec <- c(fit$alphas[i], resp.str, sample.str)
                                                                          names(vec) <- c("alpha", "predicted.response", "predicted.samples")
                                                                          vec
                                                                        })
                                                     ret.alpha$drug <- fit$drug
                                                     ret.alpha$boot <- fit$boot
                                                     ret.alpha
                                                   })
                                 ret.drug
                               })
  
  ## Flatten the actual responses
  actual.responses <- ldply(1:length(lst), .parallel = FALSE, 
                            .fun = function(drug.i) {
                              drug.lst <- lst[[drug.i]]
                              ret.drug <- ldply(drug.lst, .parallel = FALSE,
                                                .fun = function(fit) {
                                                  resps <- as.vector(t(fit$test.response))
                                                  samples <- colnames(fit$test.response)
                                                  if(length(resps) <= 1) { return(c(drug = fit$drug, actual.response = NA, actual.samples = NA)) }
                                                  resp.str <- paste0(resps, collapse=",")
                                                  sample.str <- paste0(samples, collapse=",")
                                                  vec <- c(fit$drug, fit$boot, resp.str, sample.str)
                                                  names(vec) <- c("drug", "boot", "actual.response", "actual.samples")
                                                  vec
                                                })
                              ret.drug
                            })
  
  actual.and.predicted.responses <- merge(predicted.responses, actual.responses)
  actual.and.predicted.responses <- subset(actual.and.predicted.responses, !is.na(predicted.response))
  actual.and.predicted.responses <- subset(actual.and.predicted.responses, !is.na(actual.response))
  
  cat("Done with assembling responses\n")
  
  ## 5. For each drug
  ## 6.   Calculate Pearson and Spearman correlations for actual drug response (AUC) vs median (within each sample) predicted auc
  ## 7.   Plot each of these
  ## Create a short-form matrix having columns response.type (= "predicted", "actual"), sample, response, drug, bootstrap, and alpha
  short.df <- 
    ldply(1:nrow(actual.and.predicted.responses), .parallel = FALSE,
          .fun = function(i) {
            drug <- actual.and.predicted.responses$drug[i]
            bootstrap <- actual.and.predicted.responses$bootstrap[i]
            alpha <- actual.and.predicted.responses$alpha[i]
            vec <- unlist(strsplit(actual.and.predicted.responses$predicted.response[i], split=","))
            name.vec <- unlist(strsplit(actual.and.predicted.responses$predicted.samples[i], split=","))
            pred.ret <- ldply(1:length(vec), .parallel = FALSE,
                              .fun = function(vec.i) {
                                res <- c(sample = name.vec[vec.i], response = vec[vec.i], response.type = "predicted")
                                res
                              })
            vec <- unlist(strsplit(actual.and.predicted.responses$actual.response[i], split=","))
            name.vec <- unlist(strsplit(actual.and.predicted.responses$actual.samples[i], split=","))
            actual.ret <- ldply(1:length(vec), .parallel = FALSE,
                                .fun = function(vec.i) {
                                  res <- c(sample = name.vec[vec.i], response = vec[vec.i], response.type = "actual")
                                  res
                                })
            ret <- rbind(pred.ret, actual.ret)
            ret$drug <- drug
            ret$bootstrap <- bootstrap
            ret$alpha <- alpha
            ret
          })
  cat("Done creating short.df\n")
  
  do.plot <- FALSE
  actual.predicted.correlations.df <- 
    ddply(short.df, .variables = c("alpha", "drug"), .parallel = FALSE,
          .fun = function(df) {
            df$response <- as.numeric(df$response)
            predicted.flag <- df$response.type == "predicted"
            median.pred <- ddply(df[predicted.flag,], .variables = "sample",
                                 .fun = function(df.sample) {
                                   vec <- c(predicted.response = median(as.numeric(df.sample$response)), sample = df.sample$sample[1])
                                   vec
                                 })
            actual.df <- unique(df[!predicted.flag, c("sample", "response", "alpha", "drug")])
            colnames(actual.df) <- c("sample", "actual.response", "alpha", "drug")
            m <- merge(median.pred, actual.df)
            m$predicted.response <- as.numeric(m$predicted.response)
            m$actual.response <- as.numeric(m$actual.response)
            if(do.plot) {
              g <- ggplot(data = m, aes(x = actual.response, y = predicted.response))
              g <- g + geom_point()
              g <- g + ggtitle(paste0("Drug: ", df$drug[1], " Alpha: ", df$alpha[1]))
              print(g)
            }
            ct.spearman <- cor.test(m$actual.response, m$predicted.response, method="spearman")
            ct.pearson <- cor.test(m$actual.response, m$predicted.response, method="pearson")
            vec <- c(ct.spearman$p.value, unname(ct.spearman$estimate), ct.pearson$p.value, unname(ct.pearson$estimate))
            names(vec) <- c("spearman.p", paste0("spearman.", names(ct.spearman$estimate)), "pearson.p", paste0("pearson.", names(ct.pearson$estimate)))
            vec
          })
  cat("Done creating actual predicted correlations\n")
  
  actual.predicted.correlations.df$spearman.ap <- rep(NA, nrow(actual.predicted.correlations.df))
  actual.predicted.correlations.df$pearson.ap <- rep(NA, nrow(actual.predicted.correlations.df))
  for(alpha in unique(actual.predicted.correlations.df$alpha)) {
    alpha.flag <- actual.predicted.correlations.df$alpha == alpha
    actual.predicted.correlations.df$spearman.ap[alpha.flag] <- p.adjust(actual.predicted.correlations.df$spearman.p[alpha.flag], method = "bonferroni")
    actual.predicted.correlations.df$pearson.ap[alpha.flag] <- p.adjust(actual.predicted.correlations.df$pearson.p[alpha.flag], method = "bonferroni")
  }
  
  pearson.num.sig <- ddply(actual.predicted.correlations.df, .variables = "alpha",
                           .fun = function(df) {
                             alpha <- unique(df$alpha)
                             num.sig <- length(which(df$pearson.ap < 0.05))
                             tot <- nrow(df)
                             c(alpha = alpha, num.sig = num.sig, tot.tested = tot)
                           })
  pearson.num.sig$test <- "Pearson"
  
  spearman.num.sig <- ddply(actual.predicted.correlations.df, .variables = "alpha",
                            .fun = function(df) {
                              alpha <- unique(df$alpha)
                              num.sig <- length(which(df$spearman.ap < 0.05))
                              tot <- nrow(df)
                              c(alpha = alpha, num.sig = num.sig, tot.tested = tot)
                            })
  spearman.num.sig$test <- "Spearman"
  num.sig <- rbind(pearson.num.sig, spearman.num.sig)
  
  g <- ggplot(data = num.sig, aes(x = alpha, y = num.sig))
  g <- g + geom_point()
  g <- g + facet_wrap(~ test)
  g <- g + xlab("")
  g <- g + ylab("Number Significant (Bonferonni-corrected p < 0.05)")
  g <- g + ggtitle("Number of Drugs (n = 79) with\nSignificant Prediction vs Median Response")
  pdf(paste0(prefix, "num-significant-drugs.pdf"))
  print(g)
  d <- dev.off()
  
  plts <- 
    dlply(actual.predicted.correlations.df, .variables = "alpha", .parallel = FALSE,
          .fun = function(tbl) {
            alpha <- unique(tbl$alpha)
            pval.col <- "pearson.p"
            tbl <- tbl[order(tbl[,pval.col]),]
            x.col <- "drug"
            tbl[,x.col] <- factor(tbl[,x.col], levels=unique(tbl[,x.col]))
            file <- paste0("pearson-alpha-", alpha, ".pdf")
            pdf(file, onefile=FALSE)
            g1 <- plot.stacked.pval.correlation.figures(tbl, "pearson.p", "pearson.cor", "drug", y.cor.lab = "Pearson's Correlation", main = paste0("Alpha = ", alpha), text.size = 5)
            plot(g1)
            d <- dev.off()
            file <- paste0("spearman-alpha-", alpha, ".pdf")
            pdf(file, onefile=FALSE)
            g2 <- plot.stacked.pval.correlation.figures(tbl, "spearman.p", "spearman.rho", "drug", y.cor.lab = "Spearman's Rho", main = paste0("Alpha = ", alpha), text.size = 5)
            plot(g2)
            d <- dev.off()
            return(list(g1, g2))
          })
  
  ## Plot the distribution of correlation values as a function of alpha
  g1 <- ggplot(data = actual.predicted.correlations.df, aes(x = factor(alpha), y = spearman.rho))
  g1 <- g1 + geom_violin()
  g1 <- g1 + ggtitle("Spearman's rho")
  g1 <- g1 + geom_boxplot(width = 0.5)
  g1 <- g1 + xlab("alpha")
  g1 <- g1 + ylab("Spearman's rho")
  
  g2 <- ggplot(data = actual.predicted.correlations.df, aes(x = factor(alpha), y = pearson.cor))
  g2 <- g2 + geom_violin()
  g2 <- g2 + ggtitle("Pearson's correlation")
  g2 <- g2 + geom_boxplot(width = 0.5)
  g2 <- g2 + xlab("alpha")
  g2 <- g2 + ylab("Pearson's r")
  
  pdf(paste0(prefix, "correlation-vs-alpha.pdf"))
  grid.arrange(g1, g2)
  d <- dev.off()
  
  corr.alpha.g <- arrangeGrob(grobs = list(g1, g2))
  
  tmp <- subset(actual.predicted.correlations.df, alpha == 0)
  tmp <- tmp[order(tmp$pearson.p,decreasing=FALSE),]
  
  individual.drug.fits <- list()
  for(i in 1:min(4,nrow(tmp))) {
    target.alpha <- tmp$alpha[i]
    target.drug <- tmp$drug[i]
    pdf(paste0(prefix, target.drug, "-alpha-", target.alpha, ".pdf"))
    drug.g <- plot.predicted.vs.actual.response(short.df, target.alpha = target.alpha, target.drug = target.drug)
    individual.drug.fits[[target.drug]] <- drug.g
    print(drug.g)
    d <- dev.off()
  }
  
  mydoc <- pptx()
  
  mydoc <- addSlide( mydoc, slide.layout = 'Title Slide' )
  mydoc <- addTitle( mydoc, 'Elastic Net Analysis of OHSU (AUC vs expr)' )
  mydoc <- addSubtitle( mydoc , 'June 28, 2017')
  
  mydoc <- addSlide( mydoc, slide.layout = "Title and Content" )
  ## title <- "FIMM (t = 0)"
  ## mydoc <- addTitle( mydoc, title )
  mydoc <- addPlot( doc = mydoc, fun = plot, x = corr.alpha.g)
  
  for(i in 1:length(plts)) {
    mydoc <- addSlide( mydoc, slide.layout = "Two Content" )
    ## title <- paste0("OHSU (t = ", dss.t, "): ", mechanism)
    ## mydoc <- addTitle( mydoc, title )
    mydoc <- addPlot( doc = mydoc, fun = plot, x = plts[[i]][[1]] )
    mydoc <- addPlot( doc = mydoc, fun = plot, x = plts[[i]][[2]] )
  }
  
  for(drug in names(individual.drug.fits)) {
    mydoc <- addSlide( mydoc, slide.layout = "Title and Content" )
    title <- paste0("Drug: ", drug, " Alpha: ", 0)
    mydoc <- addTitle( mydoc, title )
    mydoc <- addPlot( doc = mydoc, fun = plot, x = individual.drug.fits[[drug]])
  }
  
  
  writeDoc(mydoc, paste0(prefix, "070417-elastic-net.pptx"))
  
  cat("Done with creating correlation table\n")
  lst <- list(short.df = short.df, plts = plts, actual.predicted.correlations.df = actual.predicted.correlations.df)
  lst
}

