suppressPackageStartupMessages(library("hdi"))

model.fimm.and.ohsu <- function(ohsu.dss.arg, ohsu.expr.arg, ohsu.genomic.arg, ohsu.clinical.arg, ohsu.common.drugs, ohsu.drugs, fimm.dss.arg, fimm.expr.arg, fimm.genomic.arg, fimm.clinical.arg, fimm.common.drugs, fimm.drugs, response.col, seed = 1234, train.on.ohsu = TRUE) {

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
                         alphas <- c(0, 0.25, 0.5, 0.7, 1)
                         alphas <- c(0, 1)
                         svalues <- c("lambda.min", "lambda.1se")
                         svalues <- c("lambda.min")
                         for(alpha in alphas) {
                           for(s in svalues) {
                             set.seed(seed)
                             fit <- tryCatch({cv.glmnet(x = t(x), y, family = "gaussian", type.measure = "mse", nfolds = nfolds, alpha = alpha, standardize = FALSE, parallel = TRUE)}, error = function(e) { cat("glmnet err\n"); return(NA) })
                             pred <- tryCatch({predict(fit, newx=t(newx), s = s)}, error = function(e) { cat("glmnet pred err\n"); return(NA) })
                             ret.list[[length(ret.list)+1]] <- list(train.drug = train.drug, test.drug = test.drug, predicted = pred, actual = newy, fit = fit, alpha = alpha, s = s, model = "glmnet", n.train = n.train, n.test = n.test)
                           }
                         }

                         ## Use RF
                         set.seed(seed)
                         fit <- tryCatch({randomForest(x = t(x), y)}, error = function(e) { cat("rf error\n"); return(NA) })
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

model.train.and.test <- function(train.dss.arg, train.expr.arg, train.genomic.arg, train.clinical.arg, train.common.drugs, train.drugs, train.drug.col, train.patient.col, train.response.col, test.dss.arg, test.expr.arg, test.genomic.arg, test.clinical.arg, test.common.drugs, test.drugs, test.drug.col, test.patient.col, test.response.col, seed = 1234, use.rf = FALSE, use.svm = FALSE) {

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

                         ## Use glmnet
                         alphas <- c(0, 0.25, 0.5, 0.7, 1)
                         alphas <- c(0, 1)
                         svalues <- c("lambda.min", "lambda.1se")
                         svalues <- c("lambda.min")
                         for(alpha in alphas) {
                           for(s in svalues) {
                             set.seed(seed)
                             fit <- tryCatch({cv.glmnet(x = t(x), y, family = "gaussian", type.measure = "mse", nfolds = nfolds, alpha = alpha, standardize = FALSE, parallel = TRUE)}, error = function(e) { cat("glmnet err\n"); return(NA) })
                             pred <- tryCatch({predict(fit, newx=t(newx), s = s)}, error = function(e) { cat("glmnet pred err\n"); return(NA) })
                             ret.list[[length(ret.list)+1]] <- list(train.drug = train.drug, test.drug = test.drug, predicted = pred, actual = newy, fit = fit, alpha = alpha, s = s, model = "glmnet", n.train = n.train, n.test = n.test)
                           }
                         }

                         ## Use RF
                         if(use.rf) {
                           set.seed(seed)
                           fit <- tryCatch({randomForest(x = t(x), y)}, error = function(e) { cat("rf error\n"); return(NA) })
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
                         len <- length(newy) 
                         pred <- rep(mean(y, na.rm = TRUE), len)
                         ret.list[[length(ret.list)+1]] <- list(train.drug = train.drug, test.drug = test.drug, predicted = pred, actual = newy, fit = NA, alpha = NA, s = NA, model = "mean", n.train = n.train, n.test = n.test)

                         ret.list
                })
   ret
} 

train.model <- function(train.dss.arg, train.expr.arg, train.genomic.arg, train.clinical.arg, train.common.drugs, train.drugs, train.drug.col, train.patient.col, train.response.col, seed = 1234, use.rf = FALSE, use.svm = FALSE) {

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
     train.genomic <- train.genomic[, train.samples]
   }
   if(!is.null(train.clinical)) {
     train.clinical <- train.clinical[, train.samples]
   }
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

                         ## Use glmnet
                         alphas <- c(0, 0.25, 0.5, 0.7, 1)
                         alphas <- c(0, 1)
                         svalues <- c("lambda.min", "lambda.1se")
                         svalues <- c("lambda.min")
                         for(alpha in alphas) {
                           pvals <- NA
                           if(alpha == 1) {
                             ## lasso
                             ## Calculate the p-values for the individual components
                             fit.lasso <- lasso.proj(t(x), y)
                             pvals <- fit.lasso$pval
                           }
                           for(s in svalues) {
                             set.seed(seed)
                             fit <- tryCatch({cv.glmnet(x = t(x), y, family = "gaussian", type.measure = "mse", nfolds = nfolds, alpha = alpha, standardize = FALSE, parallel = TRUE)}, error = function(e) { cat("glmnet err\n"); return(NA) })
                             ret.list[[length(ret.list)+1]] <- list(train.drug = train.drug, fit = fit, alpha = alpha, s = s, model = "glmnet", n.train = n.train, pvals = pvals)
                           }
                         }

                         ## Use RF
                         if(use.rf) {
                           set.seed(seed)
                           fit <- tryCatch({randomForest(x = t(x), y)}, error = function(e) { cat("rf error\n"); return(NA) })
                           pred <- tryCatch({predict(fit, t(newx))}, error = function(e) { cat("rf pred error\n"); return(NA) })
                           ret.list[[length(ret.list)+1]] <- list(train.drug = train.drug, fit = fit, alpha = NA, s = NA, model = "rf", n.train = n.train, pvals = NA)
                         }

                         ## Use SVM
                         if(use.svm) {
                           set.seed(seed)
                           fit <- tryCatch({svm(x = t(x), y, scale = FALSE)}, error = function(e) { cat("svm error\n"); return(NA) })
                           pred <- tryCatch({predict(fit, t(newx))}, error = function(e) { cat("svm pred error\n"); return(NA) })
                           ret.list[[length(ret.list)+1]] <- list(train.drug = train.drug, fit = fit, alpha = NA, s = NA, model = "svm", n.train = n.train, pvals = NA)
                         }

                         ## Return a baseline in which the prediction is just the mean of training data set
                         ret.list[[length(ret.list)+1]] <- list(train.drug = train.drug, fit = NA, alpha = NA, s = NA, model = "mean", n.train = n.train, pvals = NA)

                         ret.list
                })
   ret
} 

                                              

model.train.and.test.with.downsampling <- function(train.dss.arg, train.expr.arg, train.genomic.arg, train.clinical.arg, train.common.drugs, train.drugs, train.drug.col, train.patient.col, train.response.col, test.dss.arg, test.expr.arg, test.genomic.arg, test.clinical.arg, test.common.drugs, test.drugs, test.drug.col, test.patient.col, test.response.col, seed = 1234, use.rf = FALSE, use.svm = FALSE, num.samples.per.drug, num.iterations = 100, with.replacement = FALSE) {  

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
       train.drc <- t(scale(t(train.drc), center = TRUE, scale = TRUE))
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
                                        alphas <- c(0, 0.25, 0.5, 0.7, 1)
                                        alphas <- c(0, 1)
                                        svalues <- c("lambda.min", "lambda.1se")
                                        svalues <- c("lambda.min")
                                        for(alpha in alphas) {
                                          for(s in svalues) {
                                            set.seed(seed)
                                            fit <- tryCatch({cv.glmnet(x = t(x), y, family = "gaussian", type.measure = "mse", nfolds = nfolds, alpha = alpha, standardize = FALSE, parallel = TRUE)}, error = function(e) { cat("glmnet err\n"); return(NA) })
                                            pred <- tryCatch({predict(fit, newx=t(newx), s = s)}, error = function(e) { cat("glmnet pred err\n"); return(NA) })
                                            ret.list[[length(ret.list)+1]] <- list(train.drug = train.drug, test.drug = test.drug, predicted = pred, actual = newy, fit = fit, alpha = alpha, s = s, model = "glmnet", n.train = n.train, n.test = n.test, iter = iter)
                                          }
                                        }

                                        ## Use RF
                                        if(use.rf) {
                                          set.seed(seed)
                                          fit <- tryCatch({randomForest(x = t(x), y)}, error = function(e) { cat("rf error\n"); return(NA) })
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

extract.results <- function(ret) {
  tbl <-
     ldply(ret,
             .fun = function(foo) {
                ldply(foo,
             .fun = function(df) {
                ret <- NULL

                if((class(df$predicted) != "logical") && !is.na(df$predicted)) {
                  method <- "pearson"
                  ct <- cor.test(df$predicted, df$actual, method = method)
                  ret <- data.frame(val = ct$estimate, pval = ct$p.value, train.drug = df$train.drug, test.drug = df$test.drug, alpha = df$alpha, s = df$s, model = df$model, metric = method, n.train = df$n.train, n.test = df$n.test)
                  method <- "spearman"
                  ct <- cor.test(df$predicted, df$actual, method = method)
                  ret2 <- data.frame(val = ct$estimate, pval = ct$p.value, train.drug = df$train.drug, test.drug = df$test.drug, alpha = df$alpha, s = df$s, model = df$model, metric = method, n.train = df$n.train, n.test = df$n.test)
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
                  ct <- cor.test(df$predicted, df$actual, method = method)
                  ret <- data.frame(val = ct$estimate, pval = ct$p.value, train.drug = df$train.drug, test.drug = df$test.drug, alpha = df$alpha, s = df$s, model = df$model, metric = method, n.train = df$n.train, n.test = df$n.test)
                  method <- "spearman"
                  ct <- cor.test(df$predicted, df$actual, method = method)
                  ret2 <- data.frame(val = ct$estimate, pval = ct$p.value, train.drug = df$train.drug, test.drug = df$test.drug, alpha = df$alpha, s = df$s, model = df$model, metric = method, n.train = df$n.train, n.test = df$n.test)
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
