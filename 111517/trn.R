trn.tbl <- ldply(train.set.names, .parallel = FALSE,
             .fun = function(train.set) {
                      ldply(names(gene.sets), .parallel = FALSE,
                        .fun = function(gene.set) {
                                 ldply(res.aucs[[1]][["all.fits"]][[train.set]][[gene.set]], .parallel = FALSE,
                                   .fun = function(lst) {
                                     ldply(lst, .parallel = FALSE,
                                       .fun = function(df) {
                                         ## Only pull out the regression and rf results
                                         if((df$model != "glmnet") && (df$model != "rf")) { return(NULL) }
                                         if(is.null(df$fit) || is.na(df$fit)) { return(NULL) }
print(df)
                                         coeffs <- NULL
                                         switch(df$model,
                                           "glmnet" = {
                                             coeffs <- coefficients(df$fit)
                                           },
                                           "rf" = {
                                             col <- colnames(importance(df$fit))[1]
                                             if(!grepl(col, pattern="IncMSE")) { stop(paste0("Was expected ", col, " to be %IncMSE\n")) }
                                             coeffs <- importance(df$fit)[,1,drop=FALSE]
                                           },
                                           { stop(paste0("Unknown model ", df$model, "\n")) })
                                         ns <- rownames(coeffs)
                                         coeffs <- as.vector(coeffs)
                                         names(coeffs) <- ns
                                         coeffs <- coeffs[!grepl(pattern="Intercept", names(coeffs))]
                                         n.features <- length(coeffs) 
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
