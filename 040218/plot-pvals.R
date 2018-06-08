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
