do.pca_ <- function(mat, nPcs, scale = "none", center = FALSE, pca.methods = c("nipals", "bpca", "ppca", "svdImpute")) {
  names(pca.methods) <- pca.methods
  pc.res <- llply(pca.methods, .parallel = TRUE, .fun = function(method) {
    res <- tryCatch({pca(mat, method = method, nPcs = nPcs, scale = scale, center = center)}, error = function(e) { return(NULL) })
    res    
  })
}

impute.pca_ <- function(mat, mice.imputation.methods = c("pmm", "midastouch", "sample", "cart", "rf", "mean", "norm", "norm.nob", "norm.boot", "norm.predict", "quadratic", "ri"), impute.row = FALSE) {
  names(mice.imputation.methods) <- mice.imputation.methods
  imps <- llply(mice.imputation.methods, .parallel = TRUE, .fun = function(imp) {
    mat.impute <- NULL
    if(impute.row) {
      mi <- tryCatch({mice(t(mat), method = imp, m = 1)}, error = function(e) { return(NULL) })
      if(!is.null(mi)) { mat.impute <- t(complete(mi)) }
    } else {
      mi <- tryCatch({mice(mat, method = imp, m = 1)}, error = function(e) { return(NULL) })
      if(!is.null(mi)) { mat.impute <- complete(mi) }
    }
    mat.impute
  })
  imps
}

