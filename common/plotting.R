suppressPackageStartupMessages(library(gtable))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(gridExtra))

scatterplot <- function(ds1, drug1, ds2, drug2, save.to.pdf = TRUE, use.smooth = FALSE, alpha = 0, model = "glmnet", ...) {
  par(pty = "s")  
  foo <- extract.fit(res.aucs[[1]][["all.fits"]][[ds1]][["gene"]], train.drug = drug1, alpha = alpha, model = model)
  bar <- extract.fit(res.aucs[[1]][["all.fits"]][[ds2]][["gene"]], train.drug = drug2, alpha = alpha, model = model)
  c1 <- get.fit.coeff(foo, model)
  c2 <- get.fit.coeff(bar, model)
  c1 <- c1[!grepl(pattern="Intercept", names(c1))]
  c2 <- c2[!grepl(pattern="Intercept", names(c2))]
  c2 <- c2[names(c1)]

  print(drug1)
  r1 <- rank(c1)
  r2 <- rank(c2)
  top <- 20
  flag <- (r1 < top) & (r2 < top)
  if(any(flag)) { 
    print(c1[flag])
    cat(paste(names(c1)[flag], collapse=", "), "\n") 
  }
  r1 <- rank(-c1)
  r2 <- rank(-c2)
  flag <- (r1 < top) & (r2 < top)
  if(any(flag)) { 
    print(c1[flag])
    cat(paste(names(c1)[flag], collapse=", "), "\n") 
  }

  if(save.to.pdf) { pdf(paste0(ds1, "-", drug11, "-vs-", ds2, "-", drug2, "-scatter.pdf")) }
  x = c1
  y = c2

  all = c(x,y)
##  range = c(min(all, na.rm=TRUE), max(all, na.rm=TRUE))
##  plot(c1, c2, xlab = paste0(ds1, " ", drug1, " ridge coefficients"), ylab = paste0(ds2, " ", drug2, " ridge coefficients"), xlim=range, ylim=range)
  if(use.smooth) {
    smoothScatter(c1, c2, xlab = paste0(ds1, " ", drug1, " ridge coefficients"), ylab = paste0(ds2, " ", drug2, " ridge coefficients"), xlim=range, ylim=range)
  } else {
##    plot(c1, c2, xlab = paste0(ds1, " ", drug1, " ridge coefficients"), ylab = paste0(ds2, " ", drug2, " ridge coefficients"), xlim=range, ylim=range)
    plot(c1, c2, xlab = paste0(ds1, " ", drug1, " ridge coefficients"), ylab = paste0(ds2, " ", drug2, " ridge coefficients"), ...)
  }
  reg <- lm(c2 ~ c1)
  abline(reg, col = "blue")
  if(save.to.pdf) { d <- dev.off() }
}

drug.scatterplot <- function(res.auc, gene.set, ds1, drug1, ds2, drug2, save.to.pdf = TRUE, use.smooth = FALSE, alpha = 0, model = "glmnet", ...) {
  par(pty = "s")  
  foo <- extract.fit(res.auc[["all.fits"]][[ds1]][[gene.set]], train.drug = drug1, alpha = alpha, model = model)
  bar <- extract.fit(res.auc[["all.fits"]][[ds2]][[gene.set]], train.drug = drug2, alpha = alpha, model = model)
  c1 <- get.fit.coeff(foo, model)
  c2 <- get.fit.coeff(bar, model)
  c1 <- c1[!grepl(pattern="Intercept", names(c1))]
  c2 <- c2[!grepl(pattern="Intercept", names(c2))]
  c2 <- c2[names(c1)]

  print(drug1)
  r1 <- rank(c1)
  r2 <- rank(c2)
  top <- 20
  flag <- (r1 < top) & (r2 < top)
  if(any(flag)) { 
    print(c1[flag])
    cat(paste(names(c1)[flag], collapse=", "), "\n") 
  }
  r1 <- rank(-c1)
  r2 <- rank(-c2)
  flag <- (r1 < top) & (r2 < top)
  if(any(flag)) { 
    print(c1[flag])
    cat(paste(names(c1)[flag], collapse=", "), "\n") 
  }

  if(save.to.pdf) { pdf(paste0(ds1, "-", drug11, "-vs-", ds2, "-", drug2, "-scatter.pdf")) }
  x = c1
  y = c2

  all = c(x,y)
##  range = c(min(all, na.rm=TRUE), max(all, na.rm=TRUE))
##  plot(c1, c2, xlab = paste0(ds1, " ", drug1, " ridge coefficients"), ylab = paste0(ds2, " ", drug2, " ridge coefficients"), xlim=range, ylim=range)
  if(use.smooth) {
    smoothScatter(c1, c2, xlab = paste0(ds1, " ", drug1, " ridge coefficients"), ylab = paste0(ds2, " ", drug2, " ridge coefficients"), xlim=range, ylim=range)
  } else {
##    plot(c1, c2, xlab = paste0(ds1, " ", drug1, " ridge coefficients"), ylab = paste0(ds2, " ", drug2, " ridge coefficients"), xlim=range, ylim=range)
    plot(c1, c2, xlab = paste0(ds1, " ", drug1, " ridge coefficients"), ylab = paste0(ds2, " ", drug2, " ridge coefficients"), ...)
  }
  reg <- lm(c2 ~ c1)
  abline(reg, col = "blue")
  if(save.to.pdf) { d <- dev.off() }
}


my.qqnorm <- function(vec, ...) {
  qq = qqnorm(vec, ...)
  outliers <- boxplot.stats(vec)$out
  outlier.names <- c()
  if(length(outliers) > 0) {
    indices <- which(vec %in% outliers)
    text(qq$x[indices] - 0.2, qq$y[indices], names(qq$y)[indices])
    outlier.names <- names(qq$y)[indices]
  }
  outlier.names
}

varplot = function(pca, Main=NULL, cols = "slategrey")
{
  pcaSum = summary(pca)$importance; 
  ## make the axes square and have the same range
  par(pty="s")
  xlim <- c(min(pca$x[,1], pca$x[,2]) - 0.1 * abs(min(pca$x[,1], pca$x[,2])), max(pca$x[,1], pca$x[,2]) + 0.1 * abs(max(pca$x[,1], pca$x[,2])))
  ylim <- xlim
  library(scales)
##  plot(pca$x[,1], pca$x[,2], xlab = paste("PC1", signif(pcaSum[2,1],3)),ylab = paste("PC2", signif(pcaSum[2,2],3)), pch=16, col = alpha(cols, 0.5), cex = .65, main=Main,
##       xlim = xlim, ylim = ylim)
  df <- data.frame(x = pca$x[,1], y = pca$x[,2], col = cols)
  g <- ggplot(data = df, aes(x = x, y = y, fill = col, colour = col))
  g <- g + geom_point(alpha = 0.6)
##  g <- g + geom_point(alpha = 0.2, shape = 21)
  g <- g + xlab(paste("PC1", signif(pcaSum[2,1],3)))
  g <- g + ylab(paste("PC2", signif(pcaSum[2,2],3)))
  g <- g + xlim(xlim)
  g <- g + ylim(ylim)
  g
}

varplot2 = function(pca, Main=NULL, cols = "slategrey")
{
  pcaSum = summary(pca)$importance; 
  plot(pca$x[,1], pca$x[,3], xlab = paste("PC1", signif(pcaSum[2,1],3)),ylab = paste("PC3", signif(pcaSum[2,3],3)), pch=16, col =cols,cex = .65, main=Main)
}

capwords <- function(s, strict = FALSE) {
    cap <- function(s) paste(toupper(substring(s, 1, 1)),
                  {s <- substring(s, 2); if(strict) tolower(s) else s},
                             sep = "", collapse = " " )
    sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}

lm_eqn <- function(df){
    m <- lm(y ~ x, df);
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
         list(a = format(coef(m)[1], digits = 2), 
              b = format(coef(m)[2], digits = 2), 
             r2 = format(summary(m)$r.squared, digits = 3)))
    eq <- substitute(italic(r)^2~"="~r2, 
             list(r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(eq));                 
}

lm_eqn_bq <- function(df){
    m <- lm(y ~ x, df);
    r2 = format(summary(m)$r.squared, digits = 3)
    eq <- bquote(italic(r)^2~"="~.(r2)) 
    eq
}

## Plot a few examples of predicted vs actual 
plot.r2 <- function(x, y, labels = NULL) {
  df <- data.frame(x = x, y = y)
  if(!is.null(labels)) {
    df$labels <- labels
  }
  g <- NULL
  if(is.null(labels)) {
    g <- ggplot(df, aes(x = x, y = y))
  } else {
    g <- ggplot(df, aes(x = x, y = y, label = labels))
  }
  g <- g + geom_point()
  if(!is.null(labels)) {
    g <- g + geom_text(vjust = "inward", hjust = "inward")
  }
  g <- g + theme(legend.position="none")
  g <- g + geom_smooth(method='lm')
  x.min <- min(df$x, na.rm=TRUE)
  x.max <- max(df$x, na.rm=TRUE)
  y.min <- min(df$y, na.rm=TRUE)
  y.max <- max(df$y, na.rm=TRUE)
  g <- g + geom_text(x = x.min + 0.5 * (x.max - x.min), y = y.min + 1 * (y.max - y.min), label = lm_eqn(df), parse=TRUE)
  ## g <- g + ggtitle(lm_eqn_bq(df))
  g
}

plot.table <- function(tbl, main = "", suppress.rows = TRUE, ...) {
  t1 <- NULL
  if(suppress.rows == TRUE) {
    t1 <- tableGrob(tbl, rows=NULL, ...)
  } else {
    t1 <- tableGrob(tbl, ...)
  }
  title <- textGrob(main,gp=gpar(fontsize=15))
  padding <- unit(5,"mm")

  table <- gtable_add_rows(
       t1, 
       heights = grobHeight(title) + padding,
       pos = 0)
  table <- gtable_add_grob(
      table, 
      title, 
      1, 1, 1, ncol(table))

  grid.newpage()
  grid.draw(table)
}

panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...) 
{
  usr <- par("usr"); on.exit(par(usr)) 
  par(usr = c(0, 1, 0, 1)) 
  ok <- is.finite(x) & is.finite(y) 
  r <- abs(cor(x[ok], y[ok])) 
  txt <- format(c(r, 0.123456789), digits=digits)[1] 
  txt <- paste(prefix, txt, sep="") 
  if(missing(cex.cor)) cex <- 0.8/strwidth(txt) 
  
  test <- cor.test(x[ok],y[ok]) 
  # borrowed from printCoefmat
  Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                   symbols = c("***", "**", "*", ".", " ")) 
  
  text(0.5, 0.5, txt, cex = cex * r) 
  text(.8, .8, Signif, cex=cex, col=2) 
}

panel.rank.cor <- function(x, y, digits=2, prefix="", cex.cor, ...) 
{
  usr <- par("usr"); on.exit(par(usr)) 
  par(usr = c(0, 1, 0, 1)) 
  ok <- is.finite(x) & is.finite(y) 
  if(length(which(ok)) < 3) { return() }
  test <- cor.test(x[ok], y[ok], method="spearman")
  r <- test$estimate 
  txt <- format(c(r, 0.123456789), digits=digits)[1] 
  txt <- paste(prefix, txt, sep="") 
  if(missing(cex.cor)) cex <- 0.8/strwidth(txt) 
  
  ##test <- cor.test(x[ok],y[ok]) 
  # borrowed from printCoefmat
  Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                   symbols = c("***", "**", "*", ".", " ")) 
  
  text(0.5, 0.5, txt, cex = cex * r) 
  text(.8, .8, Signif, cex=cex, col=2) 
}

panel.regression <- function (x, y, col = par("col"), bg = NA, pch = par("pch"), 
                              cex = 1, col.regres = "red", ...) 
{ 
  points(x, y, pch = pch, col = col, bg = bg, cex = cex) 
  ok <- is.finite(x) & is.finite(y) 
  if (length(which(ok)) > 3) 
    abline(stats::lm(y[ok] ~ x[ok]), col = col.regres, ...) 
  ## Plot y = x
  abline(coef=c(0,1), lty=2)
} 

panel.rank.regression <- function (x, y, col = par("col"), bg = NA, pch = par("pch"), 
                              cex = 1, col.regres = "red", ...) 
{ 
  ok <- is.finite(x) & is.finite(y) 
  points(rank(x[ok]), rank(y[ok]), pch = pch, col = col, bg = bg, cex = cex) 
  if (any(ok)) 
    abline(stats::lm(rank(y[ok]) ~ rank(x[ok])), col = col.regres, ...) 
} 
