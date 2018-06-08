
## For each of FIMM and OHSU independently, correlate PCs of expression matrix with PCs of drug response
for(train.indx in 1:length(train.set.names)) {
  train.set <- train.set.names[[train.indx]]
  l <- prepare.drug.response.and.expr.matrices(train.dss.args[[train.indx]], train.expr.args[[train.indx]], drugs = train.common.drugs[[train.indx]], 
                                               drug.col = train.drug.cols[[train.indx]], patient.col = train.patient.cols[[train.indx]], 
                                               response.col = train.response.cols[[train.indx]])
  drc.df <- l[["drc.df"]]
  expr.df <- l[["expr.df"]]

  expr.scales <- c(TRUE, FALSE)
  centers <- c(TRUE, FALSE)
  for(indx in 1:length(expr.scales)) {
  expr.scale <- expr.scales[indx]
  center <- centers[indx]
  drug.scale <- "none"
  if(expr.scale == TRUE) { drug.scale <- "uv" }

  ## Calculate the PCs of the (complete) expression matrix
##  expr.df <- scale(expr.df, center = center, scale = expr.scale)
##  expr.pca <- prcomp(t(expr.df), center = FALSE, scale. = FALSE)
  expr.pca <- prcomp(t(expr.df), center = center, scale. = expr.scale)
##  varplot(expr.pca,cols=Cols, Main = "Uncorrected")

  ## Calculate the PCs of the (incomplete) drug response matrix
  pca.methods <- c("nipals", "bpca", "ppca", "svdImpute")
##  pca.methods <- "nipals"
  mice.imputation.methods <- c("pmm", "midastouch", "sample", "cart", "rf", "mean", "norm", "norm.nob", "norm.boot", "norm.predict", "quadratic", "ri")
  library(mice)
  do.impute <- TRUE
  drc.pcas <- do.pca_(t(drc.df), nPcs = 10, scale = drug.scale, center = center, pca.methods = pca.methods)
  drc.pcas <- drc.pcas[!unlist(lapply(drc.pcas, is.null))]

  for(method in names(drc.pcas)) {
    drug.pcs <- scores(drc.pcas[[method]])
##    drug.pcs <- loadings(drc.pcas[[method]])
    ## Could use rotation if did not use transpose for prcomp
    ## Don't know what you would use in place of scores -- probably loadings
##    drug.pcs <- drc.pcas[[method]]$rotation
    n.expr <- 3
    n.drug <- 3
##    p1 <- expr.pca$rotation[,1:n.expr]
    p1 <- expr.pca$x[,1:n.expr]
    p2 <- drug.pcs[,1:n.drug]
    common <- intersect(rownames(p1), rownames(p2))
    mat <- cbind(p1[common,], p2[common,])
    colnames(mat) <- c(paste0("expr.PC", 1:n.expr), paste0("drug.PC", 1:n.drug))
    g <- ggpairs(as.data.frame(mat))
    pdf(paste(train.set, "scale", expr.scale, "center", center, method, "pairs.pdf", sep="-"))
    print(g)
    d <- dev.off()
  }

  if(!do.impute) { next }
  imps <- impute.pca_(t(drc.df), mice.imputation.methods = mice.imputation.methods)
  imps <- imps[!unlist(lapply(imps, is.null))]
  drc.pcas <- llply(imps, .parallel = TRUE, 
                    .fun = function(m) {
                             prcomp(m, center = center, scale. = expr.scale)
                    })

  drc.pcas <- drc.pcas[!unlist(lapply(drc.pcas, is.null))]
  for(method in names(drc.pcas)) {
##    drug.pcs <- drc.pcas[[method]]$rotation
    drug.pcs <- drc.pcas[[method]]$x
    n.expr <- 3
    n.drug <- 3
##    p1 <- expr.pca$rotation[,1:n.expr]
    p1 <- expr.pca$x[,1:n.expr]
    p2 <- drug.pcs[,1:n.drug]
    common <- intersect(rownames(p1), rownames(p2))
    mat <- cbind(p1[common,], p2[common,])
    colnames(mat) <- c(paste0("expr.PC", 1:n.expr), paste0("drug.PC", 1:n.drug))
    g <- ggpairs(as.data.frame(mat))
    pdf(paste(train.set, "scale", expr.scale, "center", center, "imp", method, "pairs.pdf", sep="-"))
    print(g)
    d <- dev.off()
  }

  }
}
