plot.aligned <- function(tbl1, x1, y1, tbl2, x2, y2, y1lab = NULL, x2lab = NULL, show.x.labels = TRUE, g1hline = NULL, title = NULL) {
    g1 <- ggplot(data = tbl1)
    g1 <- g1 + geom_point(aes_string(x = x1, y = y1))
    g1 <- g1 + theme(axis.text.x = element_blank())
    if(!is.null(y1lab)) {
      g1 <- g1 + ylab(y1lab)
    }
    if(!is.null(g1hline)) {
      g1 <- g1 + geom_hline(yintercept = g1hline)
    }

    g2 <- ggplot(data = tbl2)
    g2 <- g2 + geom_violin(aes_string(x = x2, y = y2))
    if(!is.null(x2lab)) {
      g2 <- g2 + xlab(x2lab)
    }
    if(show.x.labels) { 
      g2 <- g2 + theme(axis.text.x = element_text(angle = 90))
    } else {
      g2 <- g2 + theme(axis.text.x = element_blank())
    }
    g1.table <- ggplotGrob(g1)
    g1.ylab <- gtable::gtable_filter(g1.table, "ylab-l")  ## extract the y axis label
    g1.yaxis <- gtable::gtable_filter(g1.table, "axis-l") ## extract the y axis
    g1.plt <- gtable::gtable_filter(g1.table, "panel")    ## extract the plot
    g1.spacer <- gtable::gtable_filter(g1.table, "spacer") ## extract the spacer between the plot and axis label

    g2.table <- ggplotGrob(g2)
    g2.ylab <- gtable::gtable_filter(g2.table, "ylab-l")  ## extract the y axis label
    g2.yaxis <- gtable::gtable_filter(g2.table, "axis-l") ## extract the y axis
    g2.xlab <- gtable::gtable_filter(g2.table, "xlab-b")  ## extract the x axis label
    g2.xaxis <- gtable::gtable_filter(g2.table, "axis-b") ## extract the x axis
    g2.plt <- gtable::gtable_filter(g2.table, "panel")    ## extract the plot
    g2.spacer <- gtable::gtable_filter(g2.table, "spacer") ## extract the spacer between the plot and axis label

    grobs <- list(g1.ylab, g1.yaxis, g1.plt, g2.ylab, g2.yaxis, g2.plt, g2.xaxis, g2.xlab, g1.spacer, g2.spacer)
    layout <- rbind(c(1,2,3),c(NA, NA, 9), c(4,5,6), c(NA, NA, 10), c(NA, NA, 7), c(NA,NA,8))
    heights <- c(10,1,10,1,1,1)
    widths <- c(1,1,10)
    grob.table <- grid.arrange(grobs=grobs, layout_matrix=layout, heights=heights, widths=widths, top=title)
}

## Plot densities of responses as a function of patient, drug, plate, etc. (on bottom panel)
## Plot number of responses as a function of patient, drug, plate, etc. (on top panel)
## Also return a list, with a slot for each experimental factor, that gives outliers (to exclude from data) based on that experimental factor.
## In particular, these samples to exclude have an outlying (and high) number of cases outside the min.response to max.response range (relative
## to other levels of that experimental factor)
correlate.drug.response.with.experimental.factors <- function(drug.data.long, experimental.factors, response.col, min_response, max_response, prefix) {
  exclude <- list()
  for(group in experimental.factors) {
      res <- ddply(.data = drug.data.long, .variables = group, .parallel = FALSE,
                   .fun = function(df) {
                            flag <- !is.na(df[, response.col])
                            df <- df[flag, ]
                            vals <- as.numeric(df[, response.col])
                            spread <- max(df[, response.col], na.rm = TRUE) - min(df[, response.col], na.rm = TRUE)
                            short.group <- substr(as.character(df[1, group]), 1, min(10, nchar(as.character(df[1, group]))))
                            tot <- nrow(df)
                            num.extremal <- length(which( ( vals <= min_response ) | ( vals >= max_response ) ))
                            frac.extremal <- num.extremal / tot
                            data.frame(group = group, short.group = short.group, spread = spread, num = tot, frac.extremal = frac.extremal)
                   })

      df <- drug.data.long
      df$short.group <- unlist(lapply(as.character(df[,group]), function(grp) substr(grp, 1, min(10, nchar(grp)))))

      res <- res[order(res$spread),]
      levels <- unique(res[, group])
      df[, group] <- factor(df[, group], levels = levels)
      res[, group] <- factor(res[, group], levels = levels)
      levels <- unique(res$short.group)
      df$short.group <- factor(df$short.group, levels = levels)
      res$short.group <- factor(res$short.group, levels = levels)

      file <- paste0(prefix, "-", response.col, "-vs-", group, "-cnt-and-range.pdf")
      pdf(file, onefile=FALSE)
      plot.aligned(tbl1 = res, x1 = group, y1 = "num", tbl2 = df, x2 = group, y2 = response.col, x2lab = group, show.x.labels = FALSE,
                   title = paste0(prefix, " ", response.col, " vs ", group))
      d <- dev.off()

      res <- res[order(res$frac.extremal),]
      levels <- unique(res[, group])
      df[, group] <- factor(df[, group], levels = levels)
      res[, group] <- factor(res[, group], levels = levels)
      levels <- unique(res$short.group)
      df$short.group <- factor(df$short.group, levels = levels)
      res$short.group <- factor(res$short.group, levels = levels)
      print(head(res))

      ## Detect cases with outlying fraction of responses outside of range.
      ## Don't bother to do this for time_of_read, which we will subset to 24.
      outliers <- boxplot.stats(res$frac.extremal)$out
      outliers <- outliers[outliers > mean(res$frac.extremal)]
      hline <- NULL
      if( (length(outliers) > 0) && (group != "time_of_read") ) {
        hline <- min(outliers)
        exclude[[group]] <- res[res$frac.extremal >= hline, group]
      }

      file <- paste0(prefix, "-", response.col, "-vs-", group, "-frac-extremal-and-range.pdf")
      pdf(file, onefile=FALSE)
      plot.aligned(tbl1 = res, x1 = group, y1 = "frac.extremal", y1lab = paste0("Frac outside [", round(min_response), ", ", round(max_response), "]"), 
                   tbl2 = df, x2 = group, y2 = response.col, x2lab = group, show.x.labels = FALSE, g1hline = hline,
                   title = paste0(prefix, " ", response.col, " vs ", group))
      d <- dev.off()
  }
  exclude
}

approximately.monotonic <- function(vec, tolerance = 0, direction = 'inc') {
  if(length(vec) <= 1) { return(TRUE) }
  if(tolerance == 0) { return(monotonic(vec, direction = direction)) }
  switch(direction,
    inc = { },
    dec = { vec <- rev(vec) },
    { stop(paste0("Unknown direction", direction, "\n")) })
  for(i in 1:(length(vec)-1)) {
    if( (vec[i] - vec[i+1]) > tolerance ) { return(FALSE) }
  }  
  return(TRUE)
}

