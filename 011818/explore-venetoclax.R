library(ggplot2)
library(plyr)
library(glmnet)
library(randomForest)

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
  
  flag <- !is.na(c1) & !is.na(c2) & ( (c1 != 0) | (c2 != 0))
  c1 <- c1[flag]
  c2 <- c2[flag]
  df <- data.frame(x = c1, y = c2, colour = "black", label = "", stringsAsFactors = FALSE)
  ## exclude zeros in both
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
  g <- g + scale_colour_manual(values = c("black" = "black", "blue" = "blue"))
  g <- g + guides(colour = FALSE)
  var <- "feature coefficient"
  if(model == "rf") { var <- "feature importance" }

  g <- g + xlab(paste0(drug1, " ", ds1, " ", model, " ", var))
  g <- g + ylab(paste0(drug2, " ", ds2, " ", model, " ", var))

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

pts.to.label <- c("mean.resp", as.character(pts.to.label))
labels <- c("mean.resp", labels)


## ridge
g <- plot.coeffs(res.aucs[["mean.feature"]], ds1, drug1, ds2, drug2, save.to.pdf = FALSE, pts.to.label = pts.to.label, labels = labels, alpha = 0, model = "glmnet")
pdf("venetoclax-ridge-ohsu-vs-fimm.pdf")
print(g)
d <- dev.off()

g <- plot.coeffs(res.aucs[["mean.feature"]], ds1, drug1, ds2, drug2, save.to.pdf = FALSE, pts.to.label = NULL, labels = NULL, alpha = 0, model = "glmnet")
pdf("venetoclax-ridge-ohsu-vs-fimm-no-labels.pdf")
print(g)
d <- dev.off()

## lasso
g <- plot.coeffs(res.aucs[["mean.feature"]], ds1, drug1, ds2, drug2, save.to.pdf = FALSE, pts.to.label = pts.to.label, labels = labels, alpha = 1, model = "glmnet")
pdf("venetoclax-lasso-ohsu-vs-fimm.pdf")
print(g)
d <- dev.off()

## random forest
g <- plot.coeffs(res.aucs[["mean.feature"]], ds1, drug1, ds2, drug2, save.to.pdf = FALSE, pts.to.label = pts.to.label, labels = labels, alpha = NA, model = "rf")
pdf("venetoclax-rf-ohsu-vs-fimm.pdf")
print(g)
d <- dev.off()

g <- plot.coeffs(res.aucs[["mean.feature"]], ds1, drug1, ds2, drug2, save.to.pdf = FALSE, pts.to.label = NULL, labels = NULLS, alpha = NA, model = "rf")
pdf("venetoclax-rf-ohsu-vs-fimm-no-labels.pdf")
print(g)
d <- dev.off()
