model.fimm.and.ohsu <- function(ohsu.dss.arg, ohsu.expr.arg, ohsu.common.drugs, ohsu.drugs, fimm.dss.arg, fimm.expr.arg, fimm.common.drugs, fimm.drugs, response.col) {

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

   common.genes <- intersect(rownames(ohsu.expr), rownames(fimm.expr))
   ohsu.expr <- ohsu.expr[common.genes,]
   fimm.expr <- fimm.expr[common.genes,]

   ## Z-score the genes/rows (individually)
   z.score <- TRUE
   if(z.score) {
     ohsu.expr <- t(scale(t(ohsu.expr), center = TRUE, scale = TRUE))
     fimm.expr <- t(scale(t(fimm.expr), center = TRUE, scale = TRUE))
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

   ## Fit ridge to OHSU data and apply to FIMM for each drug independently.
   ret <- llply(1:length(ohsu.drugs),
                .parallel = TRUE,
                .fun = function(drug.i) {
                         ohsu.drug <- ohsu.drugs[drug.i]
                         fimm.drug <- fimm.drugs[drug.i]
                         cat(paste0("Modeling ", ohsu.drug, " using ", response.col, "\n"))
                         x <- ohsu.expr
                         y <- as.numeric(ohsu.drc[ohsu.drug,])
                         flag <- is.na(y)
                         x <- x[, !flag]
                         y <- y[!flag]

                         newx <- fimm.expr
                         newy <- as.numeric(fimm.drc[fimm.drug,])
                         flag <- is.na(newy)
                         newx <- newx[, !flag]
                         newy <- newy[!flag]

                         ## Exclude any genes that have no variation in either data set (or were set to NA/NaN by
                         ## z-scoring above)
                         constant.genes <- rownames(x)[unlist(apply(x, 1, function(row) any(!is.finite(row)) | all(row == row[1])))]
                         constant.genes <- unique(c(constant.genes, rownames(newx)[unlist(apply(newx, 1, function(row) any(!is.finite(row)) | all(row == row[1])))]))
                         x <- x[!(rownames(x) %in% constant.genes),]
                         newx <- newx[!(rownames(newx) %in% constant.genes),]

                         nfolds <- 5 
                         ret.list <- list()
                         ## Apply fit to FIMM

                         ## Use glmnet
                         alphas <- c(0, 0.25, 0.5, 0.7, 1)
                         alphas <- c(0, 1)
                         svalues <- c("lambda.min", "lambda.1se")
                         svalues <- c("lambda.min")
                         for(alpha in alphas) {
                           for(s in svalues) {
                             fit <- cv.glmnet(x = t(x), y, family = "gaussian", type.measure = "mse", nfolds = nfolds, alpha = alpha, standardize = FALSE, parallel = TRUE)
                             pred <- predict(fit, newx=t(newx), s = s)
                      
                             ret.list[[length(ret.list)+1]] <- list(ohsu.drug = ohsu.drug, fimm.drug = fimm.drug, predicted = pred, actual = newy, fit = fit, alpha = alpha, s = s, model = "glmnet")
                           }
                         }

                         ## Use RF
                         rf.fit <- randomForest(x = t(x), y)
                         pred <- predict(rf.fit, t(newx))
                         ret.list[[length(ret.list)+1]] <- list(ohsu.drug = ohsu.drug, fimm.drug = fimm.drug, predicted = pred, actual = newy, fit = fit, alpha = NA, s = NA, model = "rf")

                         ## Use SVM
                         svm.fit <- svm(x = t(x), y, scale = FALSE)
                         pred <- predict(svm.fit, t(newx))
                         ret.list[[length(ret.list)+1]] <- list(ohsu.drug = ohsu.drug, fimm.drug = fimm.drug, predicted = pred, actual = newy, fit = fit, alpha = NA, s = NA, model = "svm")
                         ret.list
                })
   ret
} 

plot.correlation.vs.drug <- function(tbl) {
  rftbl <- subset(tbl, model == "rf")
  rftbl <- rftbl[order(rftbl$cor, decreasing=TRUE),]
  ordered.levels <- unique(rftbl$ohsu.drug)
  tbl$model <- factor(tbl$model)
  tbl$ohsu.drug <- factor(tbl$ohsu.drug, levels = ordered.levels)
  tbl$cor <- as.numeric(tbl$cor)
  g <- ggplot(data = tbl, aes(x = ohsu.drug, y = cor, colour = model, group = model))
  g <- g + geom_point() + geom_line()
  g <- g + xlab("Drug")
  g <- g + ylab("Correlation")
  g <- g + theme(axis.text.x=element_text(angle = 45, hjust = 1))
  g
} 

plot.correlation.vs.model <- function(tbl) {
  tbl$model <- factor(tbl$model)
  tbl$ohsu.drug <- factor(tbl$ohsu.drug)
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

