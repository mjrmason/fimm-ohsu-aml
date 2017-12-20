suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(synapseClient))
suppressPackageStartupMessages(library("foreach"))
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("ggplot2"))
##suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("MonoInc"))
suppressPackageStartupMessages(library("data.table"))

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
prefix <- "ohsu"

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
## experimental.factors <- c("time_of_read")
## experimental.factors <- c("inhibitor")

all.cols <- c(drug.screen.cols, conc.col, response.col)

aml.diagnoses <- c("ACUTE MYELOID LEUKAEMIA (AML) AND RELATED PRECURSOR NEOPLASMS")
aml.flag <- drug.data.long$diagnosis %in% aml.diagnoses
mono.flag <- !grepl(pattern=" - ", x = as.character(drug.data.long$inhibitor))

## Limit to AML and to monotherapies
drug.data.long <- drug.data.long[aml.flag & mono.flag, ]

exclude <- list()

## Exclude fits based on time of read == 24
flag <- drug.data.long$time_of_read != 24
exclude[["time_of_read"]] <- unique(drug.data.long$time_of_read[flag])

## Store in
## FIMM_BEAT AML Collaboration/Files/BEAT AML Data/Processed Data
parentId <- "syn10083332"

output.file <- paste0("ohsu.l4.and.ll4.fits.112017.tsv")

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
responses.at.min.conc <- ddply(drug.data.long, .variables = drug.screen.cols,
                                      .fun = function(df) { (df[which(df[, conc.col] == min(df[, conc.col]))[1], ]) })
## Define min_response_min_conc and max_response_min_conc as outliers in a qqplot
## Limit to a reasonable range of responses
resp.flag <- (responses.at.min.conc[, response.col] > -100) & (responses.at.min.conc[, response.col] < 100)
vec <- responses.at.min.conc[resp.flag, response.col]

mean_response_min_conc <- mean(vec)
sd_response_min_conc <- sd(vec)

## Define min/max responses at min concentration based on outliers.
## e.g., min_response_min_conc is the _maximum_ response _less than_ the mean that is an outlier
## NB: then exclude any response that is less than or equal to this min_response_min_conc
outliers <- boxplot.stats(vec)$out
min_response_min_conc <- round(max(outliers[outliers < mean_response_min_conc]), digits = 2)
max_response_min_conc <- round(min(outliers[outliers > mean_response_min_conc]), digits = 2)

resp.flag <- (responses.at.min.conc[, response.col] <= min_response_min_conc) | (responses.at.min.conc[, response.col] >= max_response_min_conc)
outlier.responses.at.min.conc <- responses.at.min.conc[resp.flag, ]

file <- paste0(prefix, "-", response.col, "-qq.pdf")
pdf(file)
names(vec) <- 1:length(vec)
qq <- qqnorm(vec, main = paste0(data.set, " ", response.col, " at min concentration"))
qqline(vec)
abline(h = min_response_min_conc)
abline(h = max_response_min_conc)
d <- dev.off()

file <- paste0(prefix, "-", response.col, "-density.pdf")
pdf(file)
g <- ggplot(data.frame(responses.at.min.conc = vec))
g <- g + geom_density(aes(x = responses.at.min.conc))
g <- g + xlab(paste0(response.col, " at min concentration"))
g <- g + ggtitle(paste0(data.set, " ", response.col, " at min concentration"))
g <- g + geom_vline(xintercept = min_response_min_conc)
g <- g + geom_vline(xintercept = max_response_min_conc)
print(g)
d <- dev.off()

## The minimum response at the minimum concentration is also the minimum allowed response at _any_ concentration
min_response <- min_response_min_conc

## Return the response at the max concentration for each screen
## Unfortunately, OHSU data is pinned at min viability = 0; max inhibition = 100
## Hence, there is no cutoff at the high end -- make 100.1 the max response (.1 so that we don't exclude the response == 100 values)
responses.at.max.conc <- (ddply(drug.data.long, .variables = drug.screen.cols,
                                .fun = function(df) { (df[which(df[, conc.col] == max(df[, conc.col]))[1], ]) }))
max_response_max_conc <- 100.1
max_response <- max_response_max_conc

cat(paste0("min_response_min_conc: ", min_response_min_conc, "\n"))
cat(paste0("max_response_min_conc: ", max_response_min_conc, "\n"))
cat(paste0("min_response: ", min_response, "\n"))
cat(paste0("max_response_max_conc: ", max_response_max_conc, "\n"))
cat(paste0("max_response: ", max_response, "\n"))

if(max(responses.at.max.conc[, response.col]) != 100) {
  stop(paste0("For OHSU was expecting max response to be 100, was instead: ", max(responses.at.max.conc[, response.col]), "\n"))
}

resp.flag <- (responses.at.max.conc[, response.col] >= max_response_max_conc)
outlier.responses.at.max.conc <- responses.at.max.conc[resp.flag, ]

## Finally, exclude any response that is <= min_response or >= max_response
resp.flag <- (drug.data.long[, response.col] <= min_response) | (drug.data.long[, response.col] >= max_response)
outlier.responses.at.any.conc <- drug.data.long[resp.flag, ]

all.outliers <- unique(rbind(outlier.responses.at.min.conc, outlier.responses.at.max.conc, outlier.responses.at.any.conc))

## Mark outliers in the original raw data
all.outliers$outlier.response <- TRUE
drug.data.long <- merge(drug.data.long, all.outliers, all.x = TRUE)
flag <- is.na(drug.data.long$outlier.response)
drug.data.long$outlier.response[flag] <- FALSE

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

fits <- list()
fcts <- c("LL.4", "L.4")
outlier.col <- "outlier.response"
for(fct in fcts) {
  cat(paste0("Fitting drug response using ", fct, "\n"))
  fits[[fct]] <- fit.drc(drug.data.long[,c(all.cols, outlier.col)], fct = fct, drug.screen.cols = drug.screen.cols, conc.nM.col = conc.col, 
                         response.col = response.col, response.is.viability = response.is.viability, outlier.col = outlier.col)

  cat(paste0("Done fitting drug response using ", fct, "\n"))

  ## Allow at most one outlier to be dropped per curve
  ## For some reason, some of these curves have many more responses -- possibly multiple replicates have
  ## been lumped together.  Filter these out. 
  num.responses <- unlist(lapply(fits[[fct]]$uniq.concs.nM, function(str) length(unlist(strsplit(str, split=",[ ]*")))))
  tbl <- as.data.frame(table(num.responses))
  expected.num.responses <- as.numeric(as.character(tbl$num.responses[tbl$Freq == max(tbl$Freq)]))
  accepted.responses <- (expected.num.responses-1):expected.num.responses
  ## Require at least 5 pts
  accepted.responses <- accepted.responses[accepted.responses >= 5]
  cat(paste0("Limiting curves to those with the following number of responses: ", paste(accepted.responses, collapse=", "), "\n"))
  exclusion.flag <- !(as.numeric(num.responses) %in% accepted.responses)
  exclusion.col <- "num.pts.exclude"
  fits[[fct]][, exclusion.col] <- FALSE
  fits[[fct]][exclusion.flag, exclusion.col] <- TRUE
 
  ## Exclude based on cutoffs established above
  for(nm in names(exclude)) {
    exclusion.col <- paste0(nm, ".exclude")
    exclusion.flag <- !is.na(fits[[fct]][, nm]) & (fits[[fct]][, nm] %in% exclude[[nm]])
    fits[[fct]][, exclusion.col] <- FALSE
    fits[[fct]][exclusion.flag, exclusion.col] <- TRUE
  }

  if(FALSE) {
  any.neg.response <- unlist(lapply(fits[[fct]]$all.responses, function(str) any(as.numeric(unlist(strsplit(str, split=",[ ]*"))) < -30)))
  monotone <- unlist(lapply(fits[[fct]]$all.responses, function(str) {
                                                          vec <- as.numeric(unlist(strsplit(str, split=",[ ]*")))
                                                          monotonic(vec, direction='inc') || monotonic(vec, direction='dec')
                                                        }))
  }

  ## Exclude based on gof
  ## Establish a gof cutoff by comparing approximately monotonic data vs non-approximately monotonic
  ## Consider two consecutive components of a vector v[i] and v[i+1] to be approximately monotonic increasing
  ## if (v[i] - v[i+1]) < 10 -- i.e., allowing a tolerance of 10%
  any.exclude <- unlist(apply(fits[[fct]][, grepl(pattern="exclude", x=colnames(fits[[fct]]))], 1, function(row) any(row)))
  tmp <- fits[[fct]][!any.exclude,]
  approx.monotone <- unlist(lapply(tmp$all.responses, function(str) {
                                                          vec <- as.numeric(unlist(strsplit(str, split=",[ ]*")))
                                                          approximately.monotonic(vec, direction='inc', tolerance = 10) || 
                                                             approximately.monotonic(vec, direction='dec', tolerance = 10)
                                                        }))

  file <- paste0(prefix, "-", fct, "-gof-cutoff.pdf")
  pdf(file, onefile=FALSE)
  tbl <- data.frame(gof = as.numeric(tmp$gof), condition = approx.monotone)
  gof.cutoff <- round(median(tbl$gof[tbl$condition == FALSE]), digits=2)
  g1 <- ggplot(tbl)
  g1 <- g1 + geom_boxplot(aes(x = condition, y = gof))
  g1 <- g1 + xlab("Approximately Monotonic")

  g2 <- ggplot(tbl)
  g2 <- g2 + geom_density(aes(x = gof))
  g2 <- g2 + geom_vline(xintercept = gof.cutoff)
  grid.arrange(g1, g2, nrow=1, top = paste0(prefix, ": GOF cutoff = ", gof.cutoff))  
  d <- dev.off()

  exclusion.flag <- as.numeric(fits[[fct]]$converged) == 0 
  exclusion.flag <- exclusion.flag | ( !is.na(fits[[fct]]$gof) & (as.numeric(fits[[fct]]$gof) < gof.cutoff) )
  exclusion.col <- "gof.exclude"
  fits[[fct]][, exclusion.col] <- FALSE
  fits[[fct]][exclusion.flag, exclusion.col] <- TRUE

  any.exclude <- unlist(apply(fits[[fct]][, grepl(pattern="exclude", x=colnames(fits[[fct]]))], 1, function(row) any(row)))
  exclusion.flag <- any.exclude
  exclusion.col <- "any.exclude"
  fits[[fct]][, exclusion.col] <- FALSE
  fits[[fct]][exclusion.flag, exclusion.col] <- TRUE

  for(exclusion.col in colnames(fits[[fct]])[grepl(pattern="exclude", x=colnames(fits[[fct]]))]) {
      cat(paste0("Num excluded: ", exclusion.col, "\n"))
      print(table(fits[[fct]][, exclusion.col]))
  }
}

cat("Saving image\n")
rdata.file <- paste0(".Rdata.", prefix)
save.image(rdata.file)

common.cols <- c(drug.screen.cols, "min.conc.nM", "max.conc.nM", "all.concs.nM", "uniq.concs.nM", "all.responses",
                 "min.conc.w.outliers.nM", "max.conc.w.outliers.nM", "all.concs.w.outliers.nM", "uniq.concs.w.outliers.nM", "all.responses.w.outliers")
fits <- merge(fits[["LL.4"]], fits[["L.4"]], by = common.cols, suffixes = c(".ll4", ".l4"), all = TRUE)

cat("Saving image\n")
save.image(rdata.file)

write.table(file=output.file, fits, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
f <- File(output.file, parentId = parentId, synapseStore = TRUE)
synStore(f, executed = NULL, forceVersion = FALSE)

cat("Done fitting drug response data using L.4 and LL.4\n")
