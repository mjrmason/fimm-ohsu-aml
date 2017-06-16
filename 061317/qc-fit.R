library(synapseClient)
library(data.table)
library(ggplot2)
library(gridExtra)

synapseLogin()

plot.smooth.scatter <- function(data, x.col, y.col, x.label, y.label) {
  g <- ggplot(data, aes_string(x = x.col, y = y.col))
  g <- g + stat_density2d(aes(fill = ..density..^0.25), geom = "tile", contour = FALSE, n = 200)
  g <- g + scale_fill_continuous(low = "white", high = "dodgerblue4", guide=FALSE)
  g <- g + geom_point(alpha = 0.1, shape = 20)
  g <- g + xlab(x.label)
  g <- g + ylab(y.label)
  g
}

plot.scatter <- function(data, x.col, y.col, x.label, y.label) {
  g <- ggplot(data, aes_string(x = x.col, y = y.col))
  g <- g + geom_point()
  g <- g + xlab(x.label)
  g <- g + ylab(y.label)
  g
}

remove.outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

## Load in the FIMM and OHSU fits.

## path <- "fimm.dss.t0.tsv"
synId <- "syn10083488"
obj <- synGet(synId, downloadFile = TRUE)
file <- getFileLocation(obj)
fimm.dss.t0 <- as.data.frame(fread(file))

## The e parameter in both L.4 and LL.4 is IC50 (not log IC50).
## Remove outliers and then plot their comparison
fimm.dss.t0$ic50.ll4 <- remove.outliers(fimm.dss.t0$e.ll4)
fimm.dss.t0$ic50.l4 <- remove.outliers(fimm.dss.t0$e.l4)

## path <- "fimm.dss.t10.tsv"
synId <- "syn10083490"
obj <- synGet(synId, downloadFile = TRUE)
file <- getFileLocation(obj)
fimm.dss.t10 <- as.data.frame(fread(file))
fimm.dss.t10$ic50.ll4 <- remove.outliers(fimm.dss.t10$e.ll4)
fimm.dss.t10$ic50.l4 <- remove.outliers(fimm.dss.t10$e.l4)

## path <- "ohsu.dss.t0.tsv"
synId <- "syn10083494"
obj <- synGet(synId, downloadFile = TRUE)
file <- getFileLocation(obj)
ohsu.dss.t0 <- as.data.frame(fread(file))
ohsu.dss.t0$ic50.ll4 <- remove.outliers(ohsu.dss.t0$e.ll4)
ohsu.dss.t0$ic50.l4 <- remove.outliers(ohsu.dss.t0$e.l4)

## path <- "ohsu.dss.t10.tsv"
synId <- "syn10083499"
obj <- synGet(synId, downloadFile = TRUE)
file <- getFileLocation(obj)
ohsu.dss.t10 <- as.data.frame(fread(file))
ohsu.dss.t10$ic50.ll4 <- remove.outliers(ohsu.dss.t10$e.ll4)
ohsu.dss.t10$ic50.l4 <- remove.outliers(ohsu.dss.t10$e.l4)


## Plot the metrics computed in L.4 (logistic) vs LL.4 (log-logistic) fits.
plot.all.metrics <- function(data, title) {
  glist <- list()
  for(metric in c("ic50", "auc", "dss1", "dss2", "dss3")) {
    x.col <- paste0(metric, ".l4")
    y.col <- paste0(metric, ".ll4")
    x.lab <- paste0(toupper(metric), " (L4)")
    y.lab <- paste0(toupper(metric), " (LL4)")
    g <- plot.smooth.scatter(data, x.col, y.col, x.lab, y.lab)
    num <- nrow(na.omit(data[,c(x.col, y.col)]))
    g <- g + ggtitle(paste0("n = ", num))
    ## Add regression and identity lines
    g <- g + geom_smooth(method = "lm", se = FALSE)
    g <- g + geom_abline(intercept = 0, slope = 1, linetype = 2, color = "blue")
    glist[[length(glist)+1]] <- g
  }
  do.call("grid.arrange", c(glist, top = title))
}

pdf("fimm-dss-t0-metric-l4-vs-ll4-concordance.pdf")
plot.all.metrics(fimm.dss.t0, "FIMM (t = 0)")
d <- dev.off()

pdf("fimm-dss-t10-metric-l4-vs-ll4-concordance.pdf")
plot.all.metrics(fimm.dss.t10, "FIMM (t = 10)")
d <- dev.off()

pdf("ohsu-dss-t0-metric-l4-vs-ll4-concordance.pdf")
plot.all.metrics(ohsu.dss.t0, "OHSU (t = 0)")
d <- dev.off()

pdf("ohsu-dss-t10-metric-l4-vs-ll4-concordance.pdf")
plot.all.metrics(ohsu.dss.t10, "OHSU (t = 10)")
d <- dev.off()