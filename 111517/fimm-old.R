suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(synapseClient))
suppressPackageStartupMessages(library("foreach"))
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("ggplot2"))
##suppressPackageStartupMessages(library("dplyr"))

num.cores <- detectCores()
if(!is.na(num.cores) && (num.cores > 1)) {
  suppressPackageStartupMessages(library("doMC"))
  cat(paste("Registering ", num.cores-1, " cores.\n", sep=""))
  registerDoMC(cores=(num.cores-1))
}

synapseLogin()

cat("Use DRC to fit L.4 (4-parameter logistic) and LL.4 (4-paramter log logistic) functions to FIMM data.\n")

## BEGIN setup

## Read in FIMM raw drug response data
synId <- "syn8488824"
obj <- synGet(id=synId, downloadFile = TRUE)
file <- getFileLocation(obj)
drug.data.long <- read.xlsx(file)

## drug.data.long has one concentration/response per row.
## NB: concentrations are in nM.  Let's make that explicit.
drug.data.long$conc.nM <- drug.data.long$CONCENTRATION

## NB: PERCENT_INHBITION is ... inhibition!
## NB: PERCENT_INHIBITION really is a percent -- it runs from 0 to 100 (well, actually, slightly below and above).

conc.col <- "conc.nM"
response.col <- "PERCENT_INHIBITION"
response.is.viability <- FALSE

drug.screen.cols <- c("DRUG_ID", "DRUG_NAME", "DRUG_SET", "SCREEN_ID", "SAMPLE_DATE", "PATIENT_ID")
experimental.factors <- c("DRUG_NAME", "DRUG_SET", "SCREEN_ID", "PATIENT_ID")

## END setup

## Plot densities of responses as a function of patient, drug, plate, etc.
for(group in experimental.factors) {
    res <- ddply(.data = drug.data.long, .variables = group, .parallel = FALSE,
                 .fun = function(df) {
                          flag <- !is.na(df[, response.col])
                          df <- df[flag, ]
                          vals <- as.numeric(df[, response.col])
                          spread <- max(df[, response.col], na.rm = TRUE) - min(df[, response.col], na.rm = TRUE)
                          data.frame(group = df[1, group], spread = spread)
                 })
    df <- drug.data.long
    res <- res[order(res$spread),]
    levels <- unique(res[, group])
    df[, group] <- factor(df[, group], levels = levels)
    g <- ggplot(data = df)
    g <- g + geom_violin(aes_string(x = group, y = response.col))
    print(g)
}

min.cutoff <- -25
max.cutoff <- 125

min.cutoff <- -10
max.cutoff <- 110

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

## Drop at most one point per curve and do not let it be an extremal one.

## Return the response at the min concentration for each screen
responses.at.min.conc <- unlist(dlply(drug.data.long, .variables = drug.screen.cols,
                                      .fun = function(df) { (df[which(df[, conc.col] == min(df[, conc.col]))[1], response.col]) }))
## Define min_response_min_conc and max_response_min_conc as outliers in a qqplot
plot(density(responses.at.min.conc))
outliers <- boxplot.stats(responses.at.min.conc)$out
min(outliers[outliers > 0])
max(outliers[outliers < 0])
mean(responses.at.min.conc)
sd(responses.at.min.conc)

## Return the response at the max concentration for each screen
responses.at.max.conc <- unlist(dlply(drug.data.long, .variables = drug.screen.cols,
                                .fun = function(df) { (df[which(df[, conc.col] == max(df[, conc.col]))[1], response.col]) }))
plot(density(responses.at.max.conc))

## Find the peak of the gaussian near the expected max response
expected.max.response.at.max.conc <- 100
d <- density(responses.at.max.conc)
peak.indx <- which(abs(d$x - expected.max.response.at.max.conc) == min(abs(d$x - expected.max.response.at.max.conc)))
x.peak <- d$x[peak.indx]
## Reflect the distribution about the peak
vals <- responses.at.max.conc[responses.at.max.conc >= x.peak]
reflected.vals <- x.peak - abs(responses.at.max.conc[responses.at.max.conc > x.peak] - x.peak)
all.vals <- c(vals, reflected.vals)
plot(density(responses.at.max.conc))
lines(density(all.vals))


l <- dlply(drug.data.long, .variables = drug.screen.cols,
      .fun = function(df) { any((df[, response.col] < min.cutoff) | (df[, response.col] > max.cutoff)) })
table(unlist(l))
stop("stop")

stop("stop")

save.image(".Rdata")

## Scale concentrations from micro- to nano-molar
drug.data.long$conc.nM <- drug.data.long$conc.uM * 10^3 

source("../common/dss.R")
source("../common/drc-fit.R")

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
