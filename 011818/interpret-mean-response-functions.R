
do.pca <- function(mat, mean.response, main.prefix) {
  nPcs <- min(20, min(nrow(mat), ncol(mat)))
  pca.methods <- c("nipals", "bpca", "ppca", "svdImpute")
  names(pca.methods) <- pca.methods
  pc.res <- llply(pca.methods, .parallel = TRUE, .fun = function(method) {
    res <- tryCatch({pca(t(mat), method = method, nPcs = nPcs, scale = "none", center = FALSE)}, error = function(e) { return(NULL) })
    res    
  })
  for(method in pca.methods) {
    res <- pc.res[[method]]
    if(!is.null(res)) {
      pcs <- scores(res)
##      var <- sDev(res)^2
      var <- 100 * res@R2
      cm <- intersect(names(mean.response), rownames(pcs))
      pcs <- pcs[cm, ]
      resp <- mean.response[cm]
      resp.rank <- rank(abs(resp))
      col <- unlist(lapply(resp, function(r) ifelse(r < 0, "black", "steelblue2")))
      percent.nas <- 100 * length(which(is.na(mat))) / (nrow(mat) * ncol(mat))
      main <- paste0(main.prefix, " pca method = ", method, ": ", round(percent.nas, digits=1), "% NAs")
      str <- "row"
      pdf(paste0(main.prefix, "-pca-method-", method, ".pdf"))
      par(mfrow=c(1+3,2), oma=c(0,0,2,0))
      plot(var, xlab = "Index", ylab = "Percent Variance Explained")
      for(i in seq(from=1, to=5, by=2)) {
        j <- i + 1
        symbols(pcs[,i], pcs[,j], circles = resp.rank, inches=1/20, bg = col, xlab = paste0("PC", i, " (", round(as.numeric(var[i]), digits=2), "%)"),
                ylab = paste0("PC", j, " (", round(as.numeric(var[j]), digits=2), "%)"))
      }
      for(i in seq(from=1, to=4, by=1)) {
        plot.r2.base(resp, pcs[,i], xlab = "Mean Response", ylab = paste0("PC", i, " (", round(as.numeric(var[i]), digits=2), "%)"))
      }
      title(main, outer=TRUE)
      d <- dev.off()
    }
  }
}

impute.pca <- function(mat, mean.response, main.prefix, impute.row = TRUE) {
  mice.imputation.methods <- c("pmm", "midastouch", "sample", "cart", "rf", "mean", "norm", "norm.nob", "norm.boot", "norm.predict", "quadratic", "ri")
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
  print(names(imps))
  for(imp in mice.imputation.methods) {
    mat.impute <- imps[[imp]]
    if(!is.null(mat.impute)) {
      pr <- tryCatch({prcomp(t(mat.impute), scale = FALSE, center = FALSE)}, error = function(e) { return(NULL) })
      if(!is.null(pr)) {
        pcs <- pr$x
        var <- summary(pr)$importance[2, ]
        cm <- intersect(names(mean.response), rownames(pcs))
        pcs <- pcs[cm, ]
        resp <- mean.response[cm]
        resp.rank <- rank(abs(resp))
        col <- unlist(lapply(resp, function(r) ifelse(r < 0, "black", "steelblue2")))
        percent.nas <- 100 * length(which(is.na(mat))) / (nrow(mat) * ncol(mat))
        main <- paste0(main.prefix, " ", imp, " imputation: ", round(percent.nas, digits=1), "% NAs")
        str <- "row"
        if(!impute.row) { str <- "col" }
        pdf(paste0(file.prefix, "-", main.prefix, "-imp-", imp, "-pca-", str, ".pdf"))
        par(mfrow=c(3,2), oma=c(0,0,2,0))
        plot(var*100, xlab = "Index", ylab = "Percent Variance Explained")
        for(i in seq(from=1, to=9, by=2)) {
          j <- i + 1
          symbols(pcs[,i], pcs[,j], circles = resp.rank, inches=1/20, bg = col, xlab = paste0("PC", i, " (", 100 * round(as.numeric(var[i]), digits=2), "%)"),
                  ylab = paste0("PC", j, " (", 100 * round(as.numeric(var[j]), digits=2), "%)"))
        }
        title(main, outer=TRUE)
        d <- dev.off()
      }
    }
  }
}

overlap.vec.with.matrix <- function(vec, mat) {
  common.samples <- intersect(names(vec), colnames(mat))
  vec <- vec[common.samples]
  mat <- mat[, common.samples]
  indices <- 1:nrow(mat)
  names(indices) <- rownames(mat)
  corrs <- ldply(indices, .parallel = TRUE,
                        .fun = function(i) {
                                 ct <- cor.test(as.numeric(vec), as.numeric(mat[i,]))
                                 data.frame(pval = ct$p.value, cor = ct$estimate)             
                        })
  colnames(corrs) <- c("term", "pval", "cor")
  corrs
}


overlap.with.gene.set <- function(vec, expr.mat, gene.set) {
  suppressPackageStartupMessages(library("GSVA"))
  es <- gsva(as.matrix(expr.mat), gene.set, parallel.sz=num.processes, verbose=TRUE)$es.obs
  return(overlap.vec.with.matrix(vec, es))
  common.samples <- intersect(names(vec), colnames(es))
  vec <- vec[common.samples]
  es <- es[, common.samples]
  indices <- 1:nrow(es)
  names(indices) <- rownames(es)
  corrs <- ldply(indices, 
                        .fun = function(i) { 
                                 ct <- cor.test(as.numeric(vec), es[i,,drop=T])
                                 data.frame(pval = ct$p.value, cor = ct$estimate)             
                        })
  colnames(corrs) <- c("term", "pval", "cor")
  corrs
}

na.dist <- function(x,...) {
 t.dist <- dist(x,...)
 t.dist <- as.matrix(t.dist)
 t.limit <- 1.1*max(t.dist,na.rm=T)
 t.dist[is.na(t.dist)] <- t.limit
 t.dist <- as.dist(t.dist)
 return(t.dist)
}

interpret.mean.response <- function(drc.mat, expr.mat, prefix, gene.sets) {

  colnames(drc.mat) <- unlist(lapply(colnames(drc.mat), function(nm) gsub("X(.*)\\.(.*)", "\\1-\\2", nm)))
  colnames(expr.mat) <- unlist(lapply(colnames(expr.mat), function(nm) gsub("X(.*)\\.(.*)", "\\1-\\2", nm)))
  frac.na <- 0.25
  tmp <- drc.mat
  flag <- unlist(apply(tmp, 1, function(row) all(is.na(row))))
  tmp <- tmp[!flag,]
  flag <- unlist(apply(tmp, 1, function(row) length(which(is.na(row))) > frac.na * length(row)))
  tmp <- tmp[!flag,]
  flag <- unlist(apply(tmp, 2, function(col) all(is.na(col))))
  tmp <- tmp[, !flag]
  flag <- unlist(apply(tmp, 2, function(col) length(which(is.na(col))) > frac.na * length(col)))
  tmp <- tmp[, !flag]
  tmp.scaled <- t(scale(t(tmp)))

  library(gplots)
  pdf(paste0(prefix, "-heatmaps.pdf"))
  par(mfrow=c(1,2))
print(head(as.matrix(tmp)))
##  heatmap.2(tmp, scale = "none", dendrogram = "none", main = paste0(prefix, " unscaled"), trace="none", Rowv = TRUE, Colv = FALSE, dist = na.dist)
  heatmap.2(as.matrix(tmp), scale = "none", dendrogram = "none", main = paste0(prefix, " scaled"), trace="none", Rowv = TRUE, Colv = FALSE, dist = na.dist)
  heatmap.2(tmp.scaled, scale = "none", dendrogram = "none", main = paste0(prefix, " scaled"), trace="none", Rowv = TRUE, Colv = FALSE, dist = na.dist)
  d <- dev.off()
  return()

  mean.response <- colMeans(tmp.scaled, na.rm=TRUE)
  common.samples <- intersect(names(mean.response), colnames(expr.mat))
  genes <- c("ABCB1")
  ensg.genes <- c("ENSG00000085563")
  for(i in 1:length(genes)) {
    gene <- genes[i]
    ensg.gene <- ensg.genes[i]
    pdf(paste0(prefix, "-", gene, "-vs-mean-response.pdf"))
    plot.r2.base(as.numeric(expr.mat[ensg.gene, common.samples]), mean.response[common.samples], ylab = "Mean Response", xlab = paste0(gene, " Expression"), main = prefix)
    d <- dev.off()
  }
  return()

  corrs <- llply(gene.sets, .fun = function(gene.set) overlap.with.gene.set(mean.response, expr.mat, gene.set))
  corrs[["gene"]] <- overlap.vec.with.matrix(mean.response, expr.mat)
  main.prefix <- prefix
##  do.pca(tmp.scaled, mean.response, main.prefix)
##  impute.pca(tmp.scaled, mean.response, main.prefix, impute.row = TRUE)
##  impute.pca(tmp.scaled, mean.response, main.prefix, impute.row = FALSE)
  corrs

}