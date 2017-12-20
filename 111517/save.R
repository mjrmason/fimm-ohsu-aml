suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(synapseClient))
suppressPackageStartupMessages(library("foreach"))
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("ggplot2"))
##suppressPackageStartupMessages(library("dplyr"))

source("../common/dss.R")
source("../common/drc-fit.R")
source("../common/plotting.R")

num.cores <- detectCores()
if(!is.na(num.cores) && (num.cores > 1)) {
  suppressPackageStartupMessages(library("doMC"))
  cat(paste("Registering ", num.cores-1, " cores.\n", sep=""))
  registerDoMC(cores=(num.cores-1))
}

synapseLogin()

cat("Use DRC to fit L.4 (4-parameter logistic) and LL.4 (4-paramter log logistic) functions to OHSU data.\n")

## BEGIN setup

data.set <- "OHSU"

## Read in the raw OHSU inhibitor data 
## inhibitor_data_points_2017_01_12.txt
## This latest file seems to have data from the prior releases
synId <- "syn8149180"
obj <- synGet(synId, downloadFile=TRUE)
file <- getFileLocation(obj)
drug.data.long <- fread(file)
drug.data.long <- as.data.frame(drug.data.long)

## NB: OHSU data are viability; "invert" (subtract from 100%) to define inhbition.
drug.data.long$PERCENT_INHIBITION <- 100 - drug.data.long$normalized_viability
## OHSU concentrations are in uM, convert to nM
drug.data.long$conc.nM <- as.numeric(drug.data.long$well_concentration) * 10^3

conc.col <- "conc.nM"
response.col <- "PERCENT_INHIBITION"
response.is.viability <- FALSE

drug.screen.cols <- c("DRUG_ID", "DRUG_NAME", "DRUG_SET", "SCREEN_ID", "SAMPLE_DATE", "PATIENT_ID")
drug.screen.cols <- c("inhibitor", "lab_id", "patient_id", "replicant", "time_of_read")
experimental.factors <- c("inhibitor", "lab_id", "patient_id", "time_of_read")
experimental.factors <- c("time_of_read")

aml.diagnoses <- c("ACUTE MYELOID LEUKAEMIA (AML) AND RELATED PRECURSOR NEOPLASMS")
aml.flag <- drug.data.long$diagnosis %in% aml.diagnoses

## END setup

## At min concentration, percent inhibition should be zero; percent viability should be 100.
## Define mean and std dev of this peak as mean_min_conc and sd_min_conc.
## Exclude response at min concentration that is outside mean_min_conc +/- 3 * sd_min_conc.
## Or, more generally, outside min_response_min_conc and max_response_min_conc

## At max concentration, max percent inhibition should be 100; min percent viability should be 0.
## But lesser inhibitions/greater viabilities are possible.  Hence, fit a gaussian to these peaks.
## Define mean and std dev of this peak as mean_max_conc and sd_max_conc.
## Exclude response at max concentration that has a percent inhibition greater than mean_max_conc + 3 * sd_max_conc or
##                                       that has a percent viability less than mean_max_conc - 3 * sd_max_conc.

## Exclude _any_ inhibition that is outside mean_min_conc - 3 * sd_min_conc < inhibition < mean_max_conc + 3 * sd_max_conc
## or, equivalently,
##         _any_ viability that is outside mean_max_conc - 3 * sd_max_conc < viability < mean_min_conc + 3 * sd_min_conc

## Fit gaussian to response at min drug concentration.  This is easy, just take mean and std of all response 
## (after excluding outliers).

## Return the response at the min concentration for each screen
responses.at.min.conc <- unlist(dlply(drug.data.long, .variables = drug.screen.cols,
                                      .fun = function(df) { (df[which(df[, conc.col] == min(df[, conc.col]))[1], response.col]) }))
## Define min_response_min_conc and max_response_min_conc as outliers in a qqplot
## Limit to a reasonable range of responses
resp.flag <- (responses.at.min.conc > -100) & (responses.at.min.conc < 100)
vec <- responses.at.min.conc[resp.flag]

mean_response_min_conc <- mean(vec)
sd_response_min_conc <- sd(vec)

## Define min/max responses at min concentration based on outliers.
## e.g., min_response_min_conc is the _maximum_ response _less than_ the mean that is an outlier
## NB: then exclude any response that is less than or equal to this min_response_min_conc
outliers <- boxplot.stats(vec)$out
min_response_min_conc <- max(outliers[outliers < mean_response_min_conc])
max_response_min_conc <- min(outliers[outliers > mean_response_min_conc])

names(vec) <- 1:length(vec)
qq <- qqnorm(vec, main = paste0(data.set, " ", response.col, " at min concentration"))
qqline(vec)
abline(h = min_response_min_conc)
abline(h = max_response_min_conc)

g <- ggplot(data.frame(responses.at.min.conc = responses.at.min.conc[resp.flag]))
g <- g + geom_density(aes(x = responses.at.min.conc))
g <- g + xlab(paste0(response.col, " at min concentration"))
g <- g + ggtitle(paste0(data.set, " ", response.col, " at min concentration"))
g <- g + geom_vline(xintercept = min_response_min_conc)
g <- g + geom_vline(xintercept = max_response_min_conc)

## The minimum response at the minimum concentration is also the minimum allowed response at _any_ concentration
min_response <- min_response_min_conc

## Return the response at the max concentration for each screen
## Unfortunately, OHSU data is pinned at min viability = 0; max inhibition = 100
## Hence, there is no cutoff at the high end -- make 100.1 the max response (.1 so that we don't exclude the response == 100 values)
responses.at.max.conc <- unlist(dlply(drug.data.long, .variables = drug.screen.cols,
                                .fun = function(df) { (df[which(df[, conc.col] == max(df[, conc.col]))[1], response.col]) }))
max_response <- 100.1
if(max(responses.at.max.conc) != 100) {
  stop(paste0("For OHSU was expecting max response to be 100, was instead: ", max(responses.at.max.conc), "\n"))
}

## TODO
## Mark outliers to be excluded
## Change fit to include all concs/response; non excluded concs/response
## Exclude fits that have < median pts - 1 points
## Exclude fits that have poor fit (based on what?)
## Exclude fits based on time of read == 24
## Exclude fits based on patient (those with a large fraction outside of range)
## Exclude fits based on drug (those with a large fraction outside of range)

plot.aligned <- function(tbl1, x1, y1, tbl2, x2, y2, y1lab = NULL, x2lab = NULL) {
    g1 <- ggplot(data = tbl1)
    g1 <- g1 + geom_point(aes_string(x = x1, y = y1))
    g1 <- g1 + theme(axis.text.x = element_blank())
    if(!is.null(y1lab)) {
      g1 <- g1 + ylab(y1lab)
    }

    g2 <- ggplot(data = df)
    g2 <- g2 + geom_violin(aes_string(x = "short.group", y = response.col))
    if(!is.null(x2lab)) {
      g2 <- g2 + xlab(x2lab)
    }
    g2 <- g2 + theme(axis.text.x = element_text(angle = 90))

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
    grob.table <- grid.arrange(grobs=grobs, layout_matrix=layout, heights=heights, widths=widths)
}

## Plot densities of responses as a function of patient, drug, plate, etc. (on bottom panel)
## Plot number of responses as a function of patient, drug, plate, etc. (on top panel)
for(group in experimental.factors) {
    res <- ddply(.data = drug.data.long[aml.flag, ], .variables = group, .parallel = FALSE,
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
    df <- drug.data.long[aml.flag,]
    df$short.group <- unlist(lapply(as.character(df[,group]), function(grp) substr(grp, 1, min(10, nchar(grp)))))

    res <- res[order(res$spread),]
    levels <- unique(res[, group])
    df[, group] <- factor(df[, group], levels = levels)
    res[, group] <- factor(res[, group], levels = levels)
    levels <- unique(res$short.group)
    df$short.group <- factor(df$short.group, levels = levels)
    res$short.group <- factor(res$short.group, levels = levels)

    g1 <- ggplot(data = res)
    g1 <- g1 + geom_point(aes_string(x = "short.group", y = "num"))
    g1 <- g1 + theme(axis.text.x = element_blank())

    g3 <- ggplot(data = res)
    g3 <- g3 + geom_point(aes_string(x = "short.group", y = "frac.extremal"))
    g3 <- g3 + theme(axis.text.x = element_blank())
    g3 <- g3 + ylab(paste0("Frac outside [", round(min_response), ", ", round(max_response), "]"))

    g2 <- ggplot(data = df)
    g2 <- g2 + geom_violin(aes_string(x = "short.group", y = response.col))
    g2 <- g2 + xlab(group)
    g2 <- g2 + theme(axis.text.x = element_text(angle = 90))

    g1.table <- ggplotGrob(g1)
    g1.ylab <- gtable::gtable_filter(g1.table, "ylab-l")  ## extract the y axis label
    g1.yaxis <- gtable::gtable_filter(g1.table, "axis-l") ## extract the y axis
    g1.plt <- gtable::gtable_filter(g1.table, "panel")    ## extract the plot
    g1.spacer <- gtable::gtable_filter(g1.table, "spacer") ## extract the spacer between the plot and axis label

    g3.table <- ggplotGrob(g3)
    g3.ylab <- gtable::gtable_filter(g3.table, "ylab-l")  ## extract the y axis label
    g3.yaxis <- gtable::gtable_filter(g3.table, "axis-l") ## extract the y axis
    g3.plt <- gtable::gtable_filter(g3.table, "panel")    ## extract the plot
    g3.spacer <- gtable::gtable_filter(g3.table, "spacer") ## extract the spacer between the plot and axis label

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
    grob.table <- grid.arrange(grobs=grobs, layout_matrix=layout, heights=heights, widths=widths)
##    grid.arrange(g1,g2)

    grobs <- list(g3.ylab, g3.yaxis, g3.plt, g2.ylab, g2.yaxis, g2.plt, g2.xaxis, g2.xlab, g3.spacer, g2.spacer)
    layout <- rbind(c(1,2,3),c(NA, NA, 9), c(4,5,6), c(NA, NA, 10), c(NA, NA, 7), c(NA,NA,8))
    heights <- c(10,1,10,1,1,1)
    widths <- c(1,1,10)
    grob.table <- grid.arrange(grobs=grobs, layout_matrix=layout, heights=heights, widths=widths)
}

stop("stop")



l <- dlply(drug.data.long, .variables = drug.screen.cols,
      .fun = function(df) { any((df[, response.col] < min.cutoff) | (df[, response.col] > max.cutoff)) })
table(unlist(l))
stop("stop")

stop("stop")

save.image(".Rdata")

## Scale concentrations from micro- to nano-molar
drug.data.long$conc.nM <- drug.data.long$conc.uM * 10^3 


all.cols <- c(drug.screen.cols, conc.col, response.col)

cat("Fitting drug response using LL4\n")
ll4.fits <- fit.drc(drug.data.long[,all.cols], fct = "LL.4", drug.screen.cols = drug.screen.cols, conc.nM.col = conc.col, response.col = response.col, 
                    response.is.viability = response.is.viability)
cat("Done fitting drug response using LL4\n")
save.image(".Rdata")

if(nrow(ll4.fits) != nrow(drug.data)) {
  stop(paste0("Number of LL4 fits (", nrow(ll4.fits), ") != Number of original entries (", nrow(drug.data), ")\n"))
}

cat("Fitting drug response using L4\n")
l4.fits <- fit.drc(drug.data.long[,all.cols], fct = "L.4", drug.screen.cols = drug.screen.cols, conc.nM.col = conc.col, response.col = response.col, response.is.viability = FALSE)
cat("Done fitting drug response using L4\n")
save.image(".Rdata")

if(nrow(l4.fits) != nrow(drug.data)) {
  stop(paste0("Number of L4 fits (", nrow(l4.fits), ") != Number of original entries (", nrow(drug.data), ")\n"))
}

common.cols <- c(drug.screen.cols, "min.conc.nM", "max.conc.nM", "all.concs.nM", "uniq.concs.nM")
fits <- merge(ll4.fits, l4.fits, by = common.cols, suffixes = c(".ll4", ".l4"), all = TRUE)

## Store in
## FIMM_BEAT AML Collaboration/Files/Sanger Data
parentId <- "syn11287866"

path <- paste0("sanger.l4.and.ll4.fits.tsv")
write.table(file=path, fits, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
f <- File(path, parentId = parentId, synapseStore = TRUE)
synStore(f, executed = NULL, forceVersion = FALSE)

cat("Fit drug response data using L.4 and LL.4\n")

q()
