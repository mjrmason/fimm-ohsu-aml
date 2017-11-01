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

bootstrap.elastic.net.ohsu <- function(drc.df, expr.df, rnaseq.summary.df, drugs, response.col, alphas, num.bootstraps = 100) {

  lst <- prepare.ohsu.drug.response.and.expr.matrices(drc.df, expr.df, rnaseq.summary.df, drugs, response.col)
  drc.df <- lst[["drc.df"]]
  expr.df <- lst[["expr.df"]]

  common.cols <- intersect(colnames(drc.df), colnames(expr.df))
  drc.df <- drc.df[, common.cols]
  expr.df <- expr.df[, common.cols]

  res <- llply(1:num.bootstraps, .parallel = TRUE,
               .fun = function(boot.i) {
                        cat(paste0("Bootstrap iteration ", boot.i, "\n"))
                        set.seed(boot.i)

                        ## Generate a bootstrap sample and a map from bootstrap index to sample
                        boot.indices <- sample.int(n = ncol(drc.df), replace=TRUE)
                        boot.name <- paste0("V", 1:length(boot.indices))
                        bootstrap.sample.map <- data.frame(bootstrap.sample = boot.name, orig.sample = colnames(drc.df)[boot.indices])

                        drc.boot <- drc.df[, boot.indices]
                        expr.boot <- expr.df[, boot.indices]
                        colnames(drc.boot) <- boot.name
                        colnames(expr.boot) <- boot.name

                        fits <- fit.elastic.net.ohsu_(drc.boot, expr.boot, rnaseq.summary.df, drugs, alphas, seed = boot.i, do.drugs.in.parallel = FALSE, train.alphas.in.parallel = FALSE)

                        return(list(bootstrap.sample.map = bootstrap.sample.map, fits = fits))
                      })
  res
}

fit.elastic.net.ohsu <- function(drc.df, expr.df, rnaseq.summary.df, drugs, response.col, regression.method = c("elastic", "lasso", "ridge"), seed = 1234, do.drugs.in.parallel = TRUE, train.alphas.in.parallel = FALSE) {

  lst <- prepare.ohsu.drug.response.and.expr.matrices(drc.df, expr.df, rnaseq.summary.df, drugs, response.col)
  drc.df <- lst[["drc.df"]]
  expr.df <- lst[["expr.df"]]

  regression.method <- match.arg(regression.method)

  res <- llply(drugs, .parallel = do.drugs.in.parallel,
               .fun = function(drug) {
                        y <- drc.df[drug, ]
  
                        ## Drop samples that have an NA response
                        y <- y[, !is.na(y), drop=FALSE]
                        expr.y <- expr.df[, colnames(y)]
  
                        ## Use caret to create a 70%/30% split into training and test sets
                        set.seed(seed)
                        train.indices <- createDataPartition(as.numeric(y), list = FALSE, p=0.7)
                        y.train <- y[, train.indices]
                        y.test <- y[, -train.indices]

                        X.train <- expr.y[, train.indices]
                        X.test <- expr.y[, -train.indices]

                        ## Independently z-score the training and test data sets
                        y.train <- t(scale(t(y.train), center = TRUE, scale = TRUE))
                        y.test <- t(scale(t(y.test), center = TRUE, scale = TRUE))
                        X.train <- scale(X.train, center = TRUE, scale = TRUE)
                        X.test <- scale(X.test, center = TRUE, scale = TRUE)

                        ## alpha = 1: LASSO (shrink coefficients to zero)
                        ## alpha = 0: ridge
                        alphas <- switch(regression.method,
                                         lasso = 1,
                                         ridge = 0,
                                         elastic = seq(from = 0, to = 1, by = 0.05))
                        nfolds <- 5
                        N <- length(y.train)
                        foldid <- sample(rep(seq(nfolds), length = N))
##                        num.processes <- detectCores() - 1
##                        num.processes <- 1
##                        train.models <- mclapply(alphas, fit.elastic.net.alpha, x = t(X.train), y = as.numeric(y.train), foldid = foldid, nfolds = nfolds, seed = seed, mc.cores = num.processes)
                        train.models <- llply(alphas, .parallel = train.alphas.in.parallel, 
                                              .fun = function(alpha) { 
                                                       fit.elastic.net.alpha(alpha, x = t(X.train), y = as.numeric(y.train), foldid = foldid, nfolds = nfolds, seed = seed)
                                                     })

                        mse <- unlist(lapply(train.models, function(x) x$cvm[which(x$lambda == x$lambda.1se)]))
                        ## plot(alphas, mse, main="lambda.1se")

                        ## Choose the alpha parameter that gives the lowest MSE
                        best.alpha <- alphas[which(mse == min(mse))]
                        best.model <- train.models[[which(alphas == best.alpha)]]
                        preds <- predict(best.model, newx=t(X.test), s="lambda.1se", type="link")
                        ## Return the model, the predictions, and the test response.
                        vec <- list(drug = drug, model = best.model, predictions = preds, test.response = y.test)
                        return(vec)
                      })
  res
}


## alpha = 1: LASSO (shrink coefficients to zero)
## alpha = 0: ridge
fit.elastic.net.ohsu_ <- function(drc.df, expr.df, rnaseq.summary.df, drugs, alphas, seed = 1234, do.drugs.in.parallel = TRUE, train.alphas.in.parallel = FALSE) {


  res <- llply(drugs, .parallel = do.drugs.in.parallel,
               .fun = function(drug) {
                        cat(paste0("    Fitting elastic net to: ", drug, "\n"))
                        y <- drc.df[drug, ]
  
                        ## Drop samples that have an NA response
                        y <- y[, !is.na(y), drop=FALSE]
                        expr.y <- expr.df[, colnames(y)]
  
                        ## Use caret to create a 70%/30% split into training and test sets
                        set.seed(seed)
                        train.indices <- createDataPartition(as.numeric(y), list = FALSE, p=0.7)
                        y.train <- y[, train.indices]
                        y.test <- y[, -train.indices]

                        X.train <- expr.y[, train.indices]
                        X.test <- expr.y[, -train.indices]

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
                                                       fit.elastic.net.alpha(alpha, x = t(X.train), y = as.numeric(y.train), foldid = foldid, nfolds = nfolds, seed = seed)
                                                     })

                        preds <- (lapply(train.models, function(model) predict(model, newx=t(X.test), s="lambda.1se", type="link")))
                        ## Return the model, the predictions, and the test response.
                        ## Don't return the model -- it takes up too much memory.
                        rm(train.models)
                        train.models <- NULL
                        gc()
                        vec <- list(drug = drug, alphas = alphas, model = train.models, predictions = preds, test.response = y.test)
                        return(vec)
                      })
  res
}

