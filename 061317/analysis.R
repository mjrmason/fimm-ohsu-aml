## For each drug in common between OHSU and FIMM, fit elastic net to response (auc)
## Specifically:
## 1. Split into 70%/30% training/test
## 2. Standardize training/test sets separately for X (gene expr) and y (drug response)
## 3. Optimize elastic net parameters on standardized training set with cv.glmnet
## 4. Run prediction on test

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
  cv <- cv.glmnet(x, y, family = "gaussian", type.measure = "mse", foldid = foldid, nfolds = nfolds, alpha = alpha, standardize = FALSE)
  cv
}

fit.elastic.net.ohsu <- function(drc.df, expr.df, rnaseq.summary.df, drugs, response.col, regression.method = c("elastic", "lasso", "ridge"), seed = 1234) {
  regression.method <- match.arg(regression.method)

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

  res <- llply(drugs, .parallel = TRUE,
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
                        num.processes <- detectCores() - 1
                        num.processes <- 1
                        train.models <- mclapply(alphas, fit.elastic.net.alpha, x = t(X.train), y = as.numeric(y.train), foldid = foldid, nfolds = nfolds, seed = seed, mc.cores = num.processes)

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

