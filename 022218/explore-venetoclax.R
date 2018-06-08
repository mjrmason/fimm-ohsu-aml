library(ggplot2)
library(plyr)
library(glmnet)
library(randomForest)

rdata.file <- ".Rdata.overlap.120817"
load(rdata.file)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

plot.coeffs <- function(res, ds1, drug1, ds2, drug2, save.to.pdf = TRUE, use.smooth = FALSE, alpha = 0, model = "glmnet", pts.to.label = NULL, labels = NULL, exclude.mean.response = FALSE, ...) {
  par(pty = "s")  
  foo <- extract.fit(res[["all.fits"]][[ds1]][["gene"]], train.drug = drug1, alpha = alpha, model = model)
  bar <- extract.fit(res[["all.fits"]][[ds2]][["gene"]], train.drug = drug2, alpha = alpha, model = model)
  c1 <- get.fit.coeff(foo, model)
  c2 <- get.fit.coeff(bar, model)
  c1 <- c1[!grepl(pattern="Intercept", names(c1))]
  if(exclude.mean.response) {
    c1 <- c1[!grepl(pattern="mean", ignore.case=TRUE, names(c1))]
  }
  c2 <- c2[names(c1)]

  if(save.to.pdf) { pdf(paste0(ds1, "-", drug11, "-vs-", ds2, "-", drug2, "-scatter.pdf")) }
  
  flag <- !is.na(c1) & !is.na(c2) & ( (c1 != 0) | (c2 != 0))
  c1 <- c1[flag]
  c2 <- c2[flag]
  df <- data.frame(x = c1, y = c2, colour = "default", label = "", stringsAsFactors = FALSE)
  ## exclude zeros in both
  if(!is.null(pts.to.label)) {
    flag <- names(c1) %in% pts.to.label
    if(any(flag)) {
      df$colour[flag] <- "black"
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
  df.txt <- df[df$label != "",]
##  g <- g + geom_text(data = df[flag,], aes(x = x, y = y, colour = colour, label = label), hjust=0, vjust=0, ...)
  fl <- df.txt$label == "LOC200772"
  g <- g + geom_text(data = df.txt[!fl,,drop=F], aes(x = x, y = y, colour = colour, label = label), hjust=0, vjust=0, ...)
  g <- g + geom_text(data = df.txt[fl,,drop=F], aes(x = x, y = y, colour = colour, label = label), hjust=0, vjust=1, ...)
## g <- g + geom_text(data = labels, aes(x = val.1, y = val.2, label = train.drug), nudge_x = - 0.02, hjust = "right", size = 6)
  cols <- gg_color_hue(2)
  g <- g + scale_colour_manual(values = c("default" = cols[1], "black" = "black"))
  g <- g + guides(colour = FALSE)
  var <- "feature coefficient"
  if(model == "rf") { var <- "feature importance" }

  g <- g + xlab(paste0(drug1, " ", ds1, " ", model, " ", var))
  g <- g + ylab(paste0(drug2, " ", ds2, " ", model, " ", var))

  lb <- min(c(df$x, df$y))
  ub <- max(c(df$x, df$y))
##  lb <- min(lb, -ub)
##  ub <- max(-lb, ub)
  g <- g + xlim(c(lb, ub)) + ylim(c(lb, ub))

  if(save.to.pdf) { d <- dev.off() }
  g
}

source("../common/utils.R")
venetoclax.feature.syms <- c("BCL2", "CD14", "FLT3", "MAFB", "LRP1", "CREG1", "LOC200772")

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
g <- plot.coeffs(res.aucs[["unpenalized.mean.feature"]], ds1, drug1, ds2, drug2, save.to.pdf = FALSE, pts.to.label = pts.to.label, labels = labels, alpha = 0, model = "glmnet", exclude.mean.response = TRUE, size = 6)
g <- g + theme(text = element_text(size = 20))
g <- g + xlab(paste(toupper(ds1), drug1, "Ridge", "Coefficient", sep=" "))
g <- g + ylab(paste(toupper(ds2), drug2, "Ridge", "Coefficient", sep=" "))
pdf("venetoclax-ridge-ohsu-vs-fimm.pdf")
print(g)
d <- dev.off()
