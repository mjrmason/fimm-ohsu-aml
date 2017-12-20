
## Make up some data -- for each patient we will have some random number of random values
group <- "patient"
response.col <- "val"
patients <- 1:10

set.seed(1234)
df <- ldply(patients,
            .fun = function(pt) {
                     num <- sample(500:1000, 1)
                     data.frame(patient = pt, val = rnorm(num))
            })

## Count the number of values per patient
df.cnt <- ddply(.data = df, .variables = group,
                .fun = function(tbl) {
                         data.frame(patient = tbl[1, group], num = nrow(tbl))
                })

## Make sure the x axes of both plots are aligned (by making them factors
## with the same levels)
df.cnt <- df.cnt[order(df.cnt$num), ]
levels <- unique(df.cnt[, group])
df.cnt[, group] <- factor(df.cnt[, group], levels = levels)
df[, group] <- factor(df[, group], levels = levels)

## In the top panel, plot the number of responses vs patients
g1 <- ggplot(data = df.cnt)
g1 <- g1 + geom_point(aes_string(x = group, y = "num"))
g1 <- g1 + theme(axis.text.x = element_blank())

## In the bottom panel, plot the distribution of responses vs patients
g2 <- ggplot(data = df)
g2 <- g2 + geom_violin(aes_string(x = group, y = response.col))
g2 <- g2 + theme(axis.text.x = element_text(angle = 90))

## Extract components of the top plot
g1.table <- ggplotGrob(g1)
g1.ylab <- gtable::gtable_filter(g1.table, "ylab-l")  ## extract the y axis label
g1.yaxis <- gtable::gtable_filter(g1.table, "axis-l") ## extract the y axis
g1.plt <- gtable::gtable_filter(g1.table, "panel")    ## extract the plot
g1.spacer <- gtable::gtable_filter(g1.table, "spacer") ## extract the spacer between the plot and axis label

## Extract components of the bottom plot
g2.table <- ggplotGrob(g2)
g2.ylab <- gtable::gtable_filter(g2.table, "ylab-l")  ## extract the y axis label
g2.yaxis <- gtable::gtable_filter(g2.table, "axis-l") ## extract the y axis
g2.xlab <- gtable::gtable_filter(g2.table, "xlab-b")  ## extract the x axis label
g2.xaxis <- gtable::gtable_filter(g2.table, "axis-b") ## extract the x axis
g2.plt <- gtable::gtable_filter(g2.table, "panel")    ## extract the plot
g2.spacer <- gtable::gtable_filter(g2.table, "spacer") ## extract the spacer between the plot and axis label

## Re-assemble components of both plots such that they are aligned
grobs <- list(g1.ylab, g1.yaxis, g1.plt, g2.ylab, g2.yaxis, g2.plt, g2.xaxis, g2.xlab, g1.spacer, g2.spacer)
layout <- rbind(c(1,2,3),c(NA, NA, 9), c(4,5,6), c(NA, NA, 10), c(NA, NA, 7), c(NA,NA,8))
heights <- c(10,1,10,1,1,1)
widths <- c(1,1,10)
grob.table <- grid.arrange(grobs=grobs, layout_matrix=layout, heights=heights, widths=widths)

## Compare above with
## grid.arrange(g1, g2)
## which does not align axes