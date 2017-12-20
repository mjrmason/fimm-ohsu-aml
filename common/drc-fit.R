suppressPackageStartupMessages(library("drc"))

## Fit _log-logistic_ function (fct = "LL.4") of _logistic_ function (fct = "L.4") to OHSU data using DRC.
## Also, compute AUC.
## response.is.viability = TRUE implies that the drug response is percent viability
## otherwise, it is percent inhibition.
## Concentrations are assumed to be in nanomolar in conc.nM.col
## drug.screen.cols should be the columns that uniquely identify a drug x patient x replicate experiment
## e.g., in OHSU this is c("patient_id", "inhibitor", "lab_id", "replicant", "time_of_read"),
## outlier.col gives the column within data that indicates whether the individual response should be excluded
## as an outlier (TRUE) or not (FALSE).  If outlier.col is NULL, no outliers are excluded
fit.drc <- function(data, fct = "LL.4", drug.screen.cols = NULL, conc.nM.col = NULL, response.col = NULL, response.is.viability = FALSE, outlier.col = NULL) {
  ddply(data,
        drug.screen.cols,
        .parallel = TRUE,
        .fun = function(df) {
          fit <- NULL
          if(response.is.viability) {
            df[, response.col] <- 100 - df[, response.col]
          }
          ## response values are in first column and dose values in the second column
          df.no.outliers <- df
          if(!is.null(outlier.col)) {
            df.no.outliers <- df[!(df[, outlier.col]), , drop=FALSE]
          }
          if(nrow(df.no.outliers) > 0) {
            drm.df <- df.no.outliers[, c(response.col, conc.nM.col)]
            if(fct == "LL.4") {
              tryCatch({
                fit <- suppressWarnings(drm(drm.df, fct=LL.4(fixed = c(NA, NA, NA,NA), names = c("SLOPE","MIN","MAX","IC50"))))
              }, 
              error = function(e) {
                cat("Trapped error in LL.4\n")
                NULL
              })
            } else {
              tryCatch({
                fit <- suppressWarnings(drm(drm.df, fct=L.4(fixed = c(NA, NA, NA,NA), names = c("SLOPE","MIN","MAX","IC50"))))
              }, 
              error = function(e) {
                cat("Trapped error in L.4\n")
                NULL
              })
            }
          }
          min.conc <- min(df.no.outliers[, conc.nM.col])
          max.conc <- max(df.no.outliers[, conc.nM.col])
          o <- order(as.numeric(df.no.outliers[, conc.nM.col]))
          all.concs <- paste(as.numeric(df.no.outliers[o, conc.nM.col]), collapse=",")
          all.responses <- paste(as.numeric(df.no.outliers[o, response.col]), collapse=",")
          uniq.concs <- paste(unique(as.numeric(df.no.outliers[o, conc.nM.col])), collapse=",")

          min.conc.w.outliers <- min(df[, conc.nM.col])
          max.conc.w.outliers <- max(df[, conc.nM.col])
          o <- order(as.numeric(df[, conc.nM.col]))
          all.concs.w.outliers <- paste(as.numeric(df[o, conc.nM.col]), collapse=",")
          all.responses.w.outliers <- paste(as.numeric(df[o, response.col]), collapse=",")
          uniq.concs.w.outliers <- paste(unique(as.numeric(df[o, conc.nM.col])), collapse=",")

          gof <- 0
          auc <- NA
          ## Regardless of using L.4 and LL.4, the parameters have the following interpretation
          ## b.param: slope
          ## c.param: minimum asymptote
          ## d.param: maximum asymptote
          ## e.param: ic50 (not log ic50--even for LL.4)
          b.param <- 0
          c.param <- 0
          d.param <- 0
          e.param <- 0
          converged <- 0
          if(!is.null(fit) && (fit$fit$convergence || fit$convergence)) {
            converged <- 1
            ssres <- sum(residuals(fit)^2)
            resp <- na.omit(fit$dataList$resp)
            sstot <- sum((mean(resp)-resp)^2)
            gof <- 1 - (ssres/sstot)
            b.param <- coef(fit)[1]
            c.param <- coef(fit)[2]
            d.param <- coef(fit)[3]
            e.param <- coef(fit)[4]
            if(fct == "LL.4") {
              auc <- compute.ll.4.auc(b.param, c.param, d.param, e.param, min.conc, max.conc)
            } else {
              auc <- compute.l.4.auc(b.param, c.param, d.param, e.param, min.conc, max.conc, numerical = FALSE) 
            }
          }
          vec <- c(converged, b.param, c.param, d.param, e.param, gof, auc, min.conc, max.conc, all.concs, uniq.concs, all.responses, 
                   min.conc.w.outliers, max.conc.w.outliers, all.concs.w.outliers, uniq.concs.w.outliers, all.responses.w.outliers)
          names(vec) <- c("converged", "b", "c", "d", "e", "gof", "auc", "min.conc.nM", "max.conc.nM", "all.concs.nM", "uniq.concs.nM", "all.responses",
                          "min.conc.w.outliers.nM", "max.conc.w.outliers.nM", "all.concs.w.outliers.nM", "uniq.concs.w.outliers.nM", "all.responses.w.outliers")


          all.concs <- paste(as.numeric(df.no.outliers[o, conc.nM.col]), collapse=",")
          all.responses <- paste(as.numeric(df.no.outliers[o, response.col]), collapse=",")
          uniq.concs <- paste(unique(as.numeric(df.no.outliers[o, conc.nM.col])), collapse=",")

          min.conc.w.outliers <- min(df[, conc.nM.col])
          max.conc.w.outliers <- max(df[, conc.nM.col])
          o <- order(as.numeric(df[, conc.nM.col]))
          all.concs.w.outliers <- paste(as.numeric(df[o, conc.nM.col]), collapse=",")
          all.responses.w.outliers <- paste(as.numeric(df[o, response.col]), collapse=",")
          uniq.concs.w.outliers <- paste(unique(as.numeric(df[o, conc.nM.col])), collapse=",")


          vec
        })
  }
