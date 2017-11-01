## For each drug in common between OHSU and FIMM, fit elastic net to response (auc)
## Specifically:
## 1. Split into 70%/30% training/test
## 2. Standardize training/test sets separately for X (gene expr) and y (drug response)
## 3. Optimize elastic net parameters on standardized training set with cv.glmnet
## 4. Run prediction on test


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


bootstrap.elastic.net.fimm <- function(drc.df, expr.df, drugs, response.col, alphas, num.bootstraps = 100, seed = 1234) {

  lst <- prepare.fimm.drug.response.and.expr.matrices(drc.df, expr.df, drugs, response.col)
  drc.df <- lst[["drc.df"]]
  expr.df <- lst[["expr.df"]]

  cat("Done preparing data\n")

  common.cols <- intersect(colnames(drc.df), colnames(expr.df))
  drc.df <- drc.df[, common.cols]
  expr.df <- expr.df[, common.cols]

  ## For each drug
  res <- llply(drugs, .parallel = FALSE,
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

                                                   boot.name <- paste0("V", 1:length(train.indices))
                                                   drc.train <- drc.df[drug, train.indices]
                                                   expr.train <- expr.df[, train.indices]
                                                   colnames(drc.train) <- boot.name
                                                   colnames(expr.train) <- boot.name

                                                   drc.test <- drc.df[drug, test.indices]
                                                   expr.test <- expr.df[, test.indices]

                                                   fit <- fit.elastic.net_(drc.train, expr.train, drc.test, expr.test, alphas, train.alphas.in.parallel = FALSE)
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

                                                   boot.name <- paste0("V", 1:length(train.indices))
                                                   drc.train <- drc.df[drug, train.indices]
                                                   expr.train <- expr.df[, train.indices]
                                                   colnames(drc.train) <- boot.name
                                                   colnames(expr.train) <- boot.name

                                                   drc.test <- drc.df[drug, test.indices]
                                                   expr.test <- expr.df[, test.indices]

                                                   fit <- fit.elastic.net_(drc.train, expr.train, drc.test, expr.test, alphas, train.alphas.in.parallel = FALSE)
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
  rm(train.models)
  train.models <- NULL
  gc()
  vec <- list(alphas = alphas, model = train.models, predictions = preds, test.response = y.test)
  return(vec)

}

