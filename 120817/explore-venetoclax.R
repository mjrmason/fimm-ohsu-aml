library(ggplot2)

## load(".Rdata.overlap.120817")

plot.coeffs <- function(res, ds1, drug1, ds2, drug2, save.to.pdf = TRUE, use.smooth = FALSE, alpha = 0, model = "glmnet", pts.to.label = NULL, labels = NULL, ...) {
  par(pty = "s")  
  foo <- extract.fit(res[["all.fits"]][[ds1]][["gene"]], train.drug = drug1, alpha = alpha, model = model)
  bar <- extract.fit(res[["all.fits"]][[ds2]][["gene"]], train.drug = drug2, alpha = alpha, model = model)
  c1 <- get.fit.coeff(foo, model)
  c2 <- get.fit.coeff(bar, model)
  c1 <- c1[!grepl(pattern="Intercept", names(c1))]
  c2 <- c2[!grepl(pattern="Intercept", names(c2))]
  c2 <- c2[names(c1)]

  if(save.to.pdf) { pdf(paste0(ds1, "-", drug11, "-vs-", ds2, "-", drug2, "-scatter.pdf")) }
  
  df <- data.frame(x = c1, y = c2, colour = "black", label = "", stringsAsFactors = FALSE)
  if(!is.null(pts.to.label)) {
    flag <- names(c1) %in% pts.to.label
    if(any(flag)) {
      df$colour[flag] <- "blue"
      names(labels) <- pts.to.label
      df$label[flag] <- labels[names(c1)[flag]]
    }
  }
##  df$colour <- factor(df$colour)

  g <- ggplot(df, aes(x = x, y = y, colour = colour, label = label))
  for(col in unique(c("black", unique(df$colour)))) {
print(col)
    g <- g + geom_point(data = df[df$colour == col,], aes(x = x, y = y, colour = colour))
  }
##  g <- g + geom_point()
  g <- g + geom_text(hjust=0, vjust=0)
  g <- g + xlab(paste0(ds1, " ", drug1))
  g <- g + ylab(paste0(ds2, " ", drug2))

  if(save.to.pdf) { d <- dev.off() }
  g
}

source("../common/utils.R")
venetoclax.feature.syms <- c("BCL2", "CD14", "FLT3", "MAFB", "LRP1")

venetoclax.feature.gene.tbl <- symbols.to.ensg.mapping(venetoclax.feature.syms)

ds1 <- "ohsu"
ds2 <- "fimm"
drug1 <- "Venetoclax"
drug2 <- drug1
pts.to.label <- venetoclax.feature.gene.tbl$ensg
labels <- venetoclax.feature.gene.tbl$symbol

## ridge
g <- plot.coeffs(res.aucs[["mean.feature"]], ds1, drug1, ds2, drug2, save.to.pdf = FALSE, pts.to.label = pts.to.label, labels = labels, alpha = 0, model = "glmnet")

## lasso
g <- plot.coeffs(res.aucs[["mean.feature"]], ds1, drug1, ds2, drug2, save.to.pdf = FALSE, pts.to.label = pts.to.label, labels = labels, alpha = 1, model = "glmnet")

## random forest
g <- plot.coeffs(res.aucs[["mean.feature"]], ds1, drug1, ds2, drug2, save.to.pdf = FALSE, pts.to.label = pts.to.label, labels = labels, alpha = NA, model = "rf")
