suppressPackageStartupMessages(library(gtable))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(gridExtra))

output.correlations <- function(ds, processed.exprs, processed.drcs, drug.to.feature.map, drug.feature.name.to.id.map, drug.ranges, prefix, title) {
  samples <- intersect(colnames(processed.exprs[[ds]]), colnames(processed.drcs[[ds]]))
  l_ply(1:nrow(drug.to.feature.map),
        .fun = function(i) {
                 drug.id <- drug.to.feature.map$drug[i]
                 drug.name <- drug.to.feature.map$drug.name[i]
                 feature.name <- drug.to.feature.map$features[i]
                 if(!(feature.name %in% names(drug.feature.name.to.id.map))) { return() }
                 features <- drug.feature.name.to.id.map[[feature.name]]
                 features <- na.omit(features)
                 if(nrow(features) == 0) { return() }
                 features <- features[features$id %in% rownames(processed.exprs[[ds]]),]
                 if(nrow(features) == 0) { return() }
  print(features)
                 lst <- llply(1:nrow(features),
                              .fun = function(indx) {
                                       gene.id <- features$id[indx]
                                       if(!(gene.id %in% rownames(processed.exprs[[ds]]))) { return(NULL) }
                                       if(!(drug.id %in% rownames(processed.drcs[[ds]]))) { return(NULL) }
                                       gene.name <- features$name[indx]
                                       flag <- !is.na(processed.exprs[[ds]][gene.id, samples]) & !is.na(processed.drcs[[ds]][drug.id, samples])
                                       if(length(which(flag)) == 0) { return(NULL) }
                                       gene.expr <- scale(as.numeric(processed.exprs[[ds]][gene.id,samples]))
                                       drug.resp <- scale(as.numeric(processed.drcs[[ds]][drug.id,samples]))
                                       g <- plot.correlation(gene.expr, drug.resp, display.pval = TRUE, size = 6, colour = "blue")
                                       g <- g + ylab(paste("Drug Response", collapse=" "))
                                       g <- g + xlab(paste(gene.name, "Expression", collapse=" "))
                                       g <- g + theme(text = element_text(size = 15))
                                       g
                              })
                 pg <- plot_grid(plotlist = lst, labels = LETTERS[1:length(lst)])
                 min.conc <- as.numeric(drug.ranges$min.conc.nM[drug.ranges$drug.name == drug.name])
                 max.conc <- as.numeric(drug.ranges$max.conc.nM[drug.ranges$drug.name == drug.name])
                 title <- ggdraw() + draw_label(paste(title, drug.name, "vs", feature.name, "biomarkers (", min.conc, "-", max.conc, "nM)", sep = " "))
                 pg <- plot_grid(title, pg, ncol=1, rel_heights=c(0.1, 1)) # rel_heights values control title margin
                 pg
                 file <- paste0(prefix, "-", make.names(drug.name), "-vs-", make.names(feature.name), "-biomarkers.pdf")
                 ggsave(file, width = 14)
         })
}      

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

lm_corr_eqn <- function(df, method = "pearson", display.r2 = FALSE, display.pval = FALSE){
    m <- lm(y ~ x, df);
    ct <- cor.test(df$x, df$y, method = method)
    estimate <- ct$estimate
    if(display.r2 == TRUE) { estimate <- estimate*estimate }
    pval <- ct$p.value
    eq <- NULL
    if((method == "pearson") && (display.r2 == TRUE)) { 
      if(display.pval) { 
        eq <- substitute(italic(r)^2~"="~est*","~~italic(p)~"="~pval, 
                         list(est = format(estimate, digits=3, scientific=0),
                              pval = format(pval, digits=3, scientific=0)))
      } else {
        eq <- substitute(italic(r)^2~"="~est, 
                         list(est = format(estimate, digits=3, scientific=0)))

      }
    } else if((method == "pearson") && (display.r2 == FALSE)) {
      if(display.pval) { 
        eq <- substitute(italic(r)~"="~est*","~~italic(p)~"="~pval, 
                         list(est = format(estimate, digits=3, scientific=0),
                              pval = format(pval, digits=3, scientific=0)))
      } else {
        eq <- substitute(italic(r)~"="~est, 
                         list(est = format(estimate, digits=3, scientific=0)))

      }
    } else if((method == "spearman") && (display.r2 == FALSE)) {
      if(display.pval) { 
        eq <- substitute(rho~"="~est*","~~italic(p)~"="~pval, 
                         list(est = format(estimate, digits=3, scientific=0),
                              pval = format(pval, digits=3, scientific=0)))
      } else {
        eq <- substitute(rho~"="~est, 
                         list(est = format(estimate, digits=3, scientific=0)))

      }
    } else {
      stop(paste("lm_corr_eqn does not know how to handle method = ", method,  " display.r2 = ", display.r2, "\n"))
    }
    as.character(as.expression(eq));                 
}

lm_eqn_bq <- function(df){
    m <- lm(y ~ x, df);
    r2 = format(summary(m)$r.squared, digits = 3)
    eq <- bquote(italic(r)^2~"="~.(r2)) 
    eq
}

plot.r2.base <- function(x, y, ...) {
  df <- data.frame(x = x, y = y)
  plot(x, y, ...)
  abline(lm(y~x)) ## regression line (y~x)
  xdim <- c(min(x, na.rm=TRUE), max(x, na.rm=TRUE))
  ydim <- c(min(y, na.rm=TRUE), max(y, na.rm=TRUE))

  ct <- cor.test(x, y)
  r.pearson <- format(ct$estimate, digits=3)
  ## data.frame(pval = ct$p.value, cor = ct$estimate)             

##  text(xdim[1] + 0.4 * (xdim[2] - xdim[1]), ydim[2] - 0.1 * (ydim[2] - ydim[1]), labels = lm_eqn_bq(df))
  text(xdim[1] + 0.4 * (xdim[2] - xdim[1]), ydim[2] - 0.1 * (ydim[2] - ydim[1]), labels = bquote(italic(r)~"="~.(r.pearson)))
}

plot.correlation <- function(x, y, labels = NULL, colors = NULL, display.r2 = FALSE, method = "pearson", display.pval = FALSE, ...) {
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
  if(!is.null(colors)) {
    g <- g + geom_point(aes(colour = colors))
  } else {
    g <- g + geom_point()
  }
  if(!is.null(labels)) {
    g <- g + geom_text(vjust = "inward", hjust = "inward")
  }
##  g <- g + theme(legend.position="none")
  g <- g + geom_smooth(data = df, aes(x = x, y = y), method='lm')
  x.min <- min(df$x, na.rm=TRUE)
  x.max <- max(df$x, na.rm=TRUE)
  y.min <- min(df$y, na.rm=TRUE)
  y.max <- max(df$y, na.rm=TRUE)

  ylimits <- NULL
  use.ggplot.2.2.1.limit.code <- TRUE
  if(use.ggplot.2.2.1.limit.code) {
    ylimits <- ggplot_build(g)$layout$panel_ranges[[1]]$y.range
  } else {
    ylimits <- ggplot_build(g)$layout$panel_params[[1]]$y.range
  }

##  g <- g + geom_text(x = x.min + 0.5 * (x.max - x.min), y = y.min + 1 * (y.max - y.min), label = lm_corr_eqn(df, method = method, display.r2 = display.r2, display.pval = display.pval), parse=TRUE, ...)
  g <- g + geom_text(x = x.min + 0.5 * (x.max - x.min), y = 0.8 * ylimits[2], label = lm_corr_eqn(df, method = method, display.r2 = display.r2, display.pval = display.pval), parse=TRUE, ...)
  g
}

plot.r2 <- function(x, y, labels = NULL, colors = NULL) {
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
  if(!is.null(colors)) {
    g <- g + geom_point(aes(colour = colors))
  } else {
    g <- g + geom_point()
  }
  if(!is.null(labels)) {
    g <- g + geom_text(vjust = "inward", hjust = "inward")
  }
##  g <- g + theme(legend.position="none")
  g <- g + geom_smooth(data = df, aes(x = x, y = y), method='lm')
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

bee.dodge.width <- 0.75
bee.sz <- 3

data2npc <- function(x, range) scales::rescale(c(range, x), c(0,1))[-c(1,2)]

pval.to.text <- function(pval, plot.pvals.as.stars = TRUE, stat.name = "p") {
    if(is.na(pval)) { return("n.s.") }
    if(!plot.pvals.as.stars) {
        eq <- substitute(italic(var)~"="~pval, list(var = stat.name, pval = ifelse(pval < 0.001, format(pval, digits = 1, scientific=TRUE), format(pval, digits = 1))))
        return(eq)
##        return(as.character(signif(pval, digits=2)))
    }
    if(pval >= 0.05) {
        return("n.s.")
    } else if(pval < 0.0001) {
        return("****")
    } else if(pval < 0.001) {
        return("***")
    } else if(pval < 0.01) {
        return("**")
    } else if(pval < 0.05) {
        return("*")
    }
    die("Should not be here\n")
}

## Draw an error bar between two facets (that may or may not be the same).
## panel.maxs[faceti-xj] gives the max value for x value xj in facet faceti
draw.err.bar <- function(pval, ranges, panel.maxs, g, facet1, x1, facet2, x2, facetIndx1, facetIndx2, xis, yoffset, plot.pvals.as.stars = TRUE, stat.name = "p", star.size = 35, ns.size = 22) {

    ## Get the maximum values across the two facets.
    m <- panel.maxs[paste0(facet1, "-", xis[1])]
    for(xi in xis) {
      m <- max(m, panel.maxs[paste0(facet1, "-", xi)])
      m <- max(m, panel.maxs[paste0(facet2, "-", xi)])
    }

    ## Get the length of the yaxis and define a step size, yoffset, with
    ## respect to it.  We will position the error bars in units of this
    ## step size.
##    ymax <- max(ranges[[facetIndx1]][["y.range"]]) - min(ranges[[facetIndx1]][["y.range"]])
##    yoffset <- 0.01 * ymax
    
    ## Use data2npc to translate from data space coordinates to
    ## coordinates used by ggplot within the facet.

    ## This is the vertical line from the top of the KRAS mutant or WT expression
    ## data, which are x offset 1 or 2 within the CMS/facet, to what will be the horizontal line of the error bar.
    x1.index <- which(x1 == xis)
    x2.index <- which(x2 == xis)
    start <- c(data2npc(x1.index,ranges[[facetIndx1]][["x.range"]]),
               data2npc(m + yoffset, ranges[[facetIndx1]][["y.range"]]))
    
    end <- c(data2npc(x1.index,ranges[[facetIndx1]][["x.range"]]),
             data2npc(m + 2 * yoffset, ranges[[facetIndx1]][["y.range"]]))

    ## The first grob/panel is 4.
    l1 <- 2 + 2 * facetIndx1
    l2 <- 2 + 2 * facetIndx2
    ## This used to work
    t <- 4
    ## This is required in ggplot 2.2.0
    t <- 8
    
    ## Give the grob a random/unique name.  I don't know whether this is
    ## required.
    name=stringi::stri_rand_strings(1,5)
    ## I don't know why this delta is necessary--ggplot2 seems to be
    ## incorrectly confusing/reordering the grobs if I do not make them unique.
    delta <- runif(n=1, min=10^-5, max=10^-3)

    ## Set the current position to the start of the vertical line.
    g <- gtable_add_grob(g, grid.move.to(start[1],start[2],draw=FALSE,name=name), z = Inf, t = t + delta, l = l1, b = 4, r = l1)

    ## Draw line from the current position to the end of the vertical line
    ## (and set that end point as the current position)
    name=stringi::stri_rand_strings(1,5)
    delta <- runif(n=1, min=10^-5, max=10^-3)    
    g <- gtable_add_grob(g, lineToGrob(end[1],end[2],name=name), z = Inf, t = t + delta, l = l1, r = l1, b = 4)

    ## Similarly, draw the horizontal line of the error bar--note that this spans from
    ## one facet to another 
    start <- c(data2npc(x1.index,ranges[[facetIndx1]][["x.range"]]),
               data2npc(m + 2 * yoffset, ranges[[facetIndx1]][["y.range"]]))
    
    end <- c(data2npc(x2.index,ranges[[facetIndx2]][["x.range"]]),
             data2npc(m + 2 * yoffset, ranges[[facetIndx2]][["y.range"]]))

    name=stringi::stri_rand_strings(1,5)
    delta <- runif(n=1, min=10^-5, max=10^-3)    
    g <- gtable_add_grob(g, lineToGrob(end[1],end[2],name=name), z = Inf, t = t + delta, l = l2, r = l2, b = 4)
    
    ## Finally, draw the vertical line on the "right side" of the error bar--
    ## this goes from the error bar to the maximum of the MT KRAS
    ## expression values in the second facet/CMS 2.
    start <- c(data2npc(x2.index,ranges[[facetIndx2]][["x.range"]]),
               data2npc(m + 2 * yoffset, ranges[[facetIndx2]][["y.range"]]))
    
    end <- c(data2npc(x2.index,ranges[[facetIndx2]][["x.range"]]),
             data2npc(m + 1 * yoffset, ranges[[facetIndx2]][["y.range"]]))
    
    name=stringi::stri_rand_strings(1,5)
    delta <- runif(n=1, min=10^-5, max=10^-3)    
    g <- gtable_add_grob(g, lineToGrob(end[1],end[2],name=name), z = Inf, t = t + delta, l = l2, r = l2, b = 4)

    ## Update the maximum values used within each facet to reflect the
    ## newly added error bars
    for(xi in xis) {
      panel.maxs[paste0(facet1,"-",xi)] <- m + 2 * 6 * yoffset
      panel.maxs[paste0(facet2,"-",xi)] <- m + 2 * 6 * yoffset
    }

    ## Add the asterisk designation of the pvalue.
    text <- pval.to.text(pval, plot.pvals.as.stars = plot.pvals.as.stars, stat.name = stat.name)

    ## Position the text in the middle of the error bar.
    xrange <- ranges[[facetIndx1]][["x.range"]]
    ## I believe this 3 comes from the fact that the high range is 2.6 and I was
    ## padding for the space between the facets
    xrange[2] <- xrange[2] + 3 * abs(facetIndx2 - facetIndx1)
    ## pos <- c(0.5 * ( data2npc(x1.index,ranges[[facetIndx1]][["x.range"]]) + data2npc(x2.index,ranges[[facetIndx2]][["x.range"]]) ), data2npc(m + 4 * yoffset, ranges[[facetIndx1]][["y.range"]]))
    num.dy <- 2
    sz <- star.size
    num.dy <- 0
    vjust <- 0.5
    if((text == "n.s.") || (plot.pvals.as.stars)) {
        num.dy <- 4
        sz <- ns.size
    }
    if(text == "n.s.") {
        vjust <- 0
    }
    pos <- c(data2npc(0.5 * ( x1.index + 3 * abs(facetIndx2 - facetIndx1) + x2.index), xrange),
             data2npc(m + num.dy * yoffset, ranges[[facetIndx1]][["y.range"]]))
    ## pos <- c(data2npc(0.5, xrange), data2npc(m + 4 * yoffset, ranges[[facetIndx1]][["y.range"]]))
    


    delta <- runif(n=1, min=10^-5, max=10^-3)
    g <- gtable_add_grob(g, textGrob(text, pos[1], pos[2], gp=gpar(fontsize = sz), vjust = vjust), t = t + delta, l = l1, b = 4, r = l2)
    return(list(g=g, panel.maxs=panel.maxs))
}

