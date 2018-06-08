library(glmnet)
suppressPackageStartupMessages(library("randomForest"))
suppressPackageStartupMessages(library("e1071"))
suppressPackageStartupMessages(library("foreach"))
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(tidyr))

suppressPackageStartupMessages(library("minpack.lm"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("gridExtra"))
suppressPackageStartupMessages(library("caTools"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("pcaMethods"))
## suppressPackageStartupMessages(library("GGally"))
suppressPackageStartupMessages(library("VennDiagram"))
suppressPackageStartupMessages(library("caret"))

library(ggbeeswarm)
library(plyr)
library(stringr)
library(synapseClient)

num.cores <- detectCores()
if(!is.na(num.cores) && (num.cores > 1)) {
  suppressPackageStartupMessages(library("doMC"))
  cat(paste("Registering ", num.cores-1, " cores.\n", sep=""))
  registerDoMC(cores=(num.cores-1))
}

## Get the FIMM raw data

## Read in FIMM raw drug response data
synId <- "syn8488824"
obj <- synGet(id=synId, downloadFile = TRUE)
file <- getFileLocation(obj)
fimm.raw.dr <- read.xlsx(file)

## Read in FIMM drug annotations
synId <- "syn8434109"
obj <- synGet(id=synId, downloadFile = TRUE)
file <- getFileLocation(obj)
fimm.drug.annotations <- read.xlsx(file)

## Read in the subset of drug annotations in both FIMM and OHSU
synId <- "syn9731315"
obj <- synGet(id=synId, downloadFile = TRUE)
ohsu.fimm.drugs <- openxlsx:::read.xlsx(getFileLocation(obj), sheet=1)

mek.ohsu.fimm.drugs <- ohsu.fimm.drugs[grepl(ohsu.fimm.drugs$`Mechanism/Targets`, pattern="MEK"),]

## Read in the raw OHSU inhibitor data 
## inhibitor_data_points_2017_01_12.txt
## This latest file seems to have data from the prior releases
synId <- "syn8149180"
inh.obj <- synGet(synId, downloadFile=TRUE)
ohsu.inh.tbl <- fread(getFileLocation(inh.obj))
ohsu.inh.tbl <- as.data.frame(ohsu.inh.tbl)

mek.ohsu.raw.tbl <- subset(ohsu.inh.tbl, inhibitor == mek.ohsu.fimm.drugs$ID_Drug.ohsu[1])
mek.fimm.raw.tbl <- subset(fimm.raw.dr, DRUG_ID == mek.ohsu.fimm.drugs$ID_Drug[1])

plot.LL.4.curve <- function(slope, min.asymptote, max.asymptote, ic50, x.min, x.max) {
  x.seq <- seq(from = x.min, to = x.max, by = 0.01)
  b <- slope
  c <- min.asymptote
  d <- max.asymptote
  e <- ic50
  y <- unlist(lapply(x.seq, function(x) {
    c + ( (d-c)/(1+exp(b*(log(x)-log(e)))) )
  }))
  plot(log10(x.seq), y, type="l")
}

lines.LL.4.curve <- function(slope, min.asymptote, max.asymptote, ic50, x.min, x.max) {
  x.seq <- seq(from = x.min, to = x.max, by = 0.01)
  b <- slope
  c <- min.asymptote
  d <- max.asymptote
  e <- ic50
  y <- unlist(lapply(x.seq, function(x) {
    c + ( (d-c)/(1+exp(b*(log(x)-log(e)))) )
  }))
  lines(log10(x.seq), y, type="l")
}

shade.LL.4.curve <- function(slope, min.asymptote, max.asymptote, ic50, x.min, x.max) {
  x.seq <- seq(from = x.min, to = x.max, by = 0.01)
  b <- slope
  c <- min.asymptote
  d <- max.asymptote
  e <- ic50
  y <- unlist(lapply(x.seq, function(x) {
    c + ( (d-c)/(1+exp(b*(log(x)-log(e)))) )
  }))
  polygon(c(min(log10(x.seq)), log10(x.seq)), c(min(y), y))
}

library(drc)
dat <- subset(mek.fimm.raw.tbl, SCREEN_ID == mek.fimm.raw.tbl$SCREEN_ID[1])
## drc.fit <- drm(normalized_viability ~ well_concentration, data = h, fct=LL.4(fixed = c(NA, NA, NA,NA), names = c("SLOPE","MIN","MAX","IC50")))

drc.fit <- drm(100 - PERCENT_INHIBITION ~ CONCENTRATION, data = dat, fct=LL.4(fixed = c(NA, NA, NA,NA), names = c("SLOPE","MIN","MAX","IC50")))

source("../common/dss.R")
pdf("mek-fimm-raw-data.pdf")
plot(log10(dat$CONCENTRATION), 100 - dat$PERCENT_INHIBITION, xlab = "log10 Concentration (nM)", ylab = "Viability")
d <- dev.off()

pdf("mek-fimm-raw-data-fit.pdf")
plot(log10(dat$CONCENTRATION), 100 - dat$PERCENT_INHIBITION, xlab = "log10 Concentration (nM)", ylab = "Viability")
lines.LL.4.curve(coef(drc.fit)[1], coef(drc.fit)[2], coef(drc.fit)[3], coef(drc.fit)[4], x.min=0.02, x.max=10000)
d <- dev.off()

pdf("mek-fimm-raw-data-fit-ics50.pdf")
plot(log10(dat$CONCENTRATION), 100 - dat$PERCENT_INHIBITION, xlab = "log10 Concentration (nM)", ylab = "Viability")
lines.LL.4.curve(coef(drc.fit)[1], coef(drc.fit)[2], coef(drc.fit)[3], coef(drc.fit)[4], x.min=0.02, x.max=10000)
ic50 <- as.numeric((coef(drc.fit)[4]))
ic50.via <- ll.4.func(ic50, coef(drc.fit)[1], coef(drc.fit)[2], coef(drc.fit)[3], coef(drc.fit)[4])
segments(x0 = log10(ic50), y0 = ic50.via, x1 = log10(0.02), y1 = ic50.via, lty = 2)
segments(x0 = log10(ic50), y0 = ic50.via, x1 = log10(ic50), y1 = 0, lty = 2)
d <- dev.off()

dat <- subset(mek.ohsu.raw.tbl, patient_id == mek.ohsu.raw.tbl$patient_id[1])
dat$well_concentration <- dat$well_concentration * 10^3
ohsu.drc.fit <- drm(normalized_viability ~ well_concentration, data = dat, fct=LL.4(fixed = c(NA, NA, NA,NA), names = c("SLOPE","MIN","MAX","IC50")))

abline(h = ic50.via, lty=2)
ll.4.func <- function(x, b, c, d, e) {
  c + ( (d-c)/(1+exp(b*(log(x)-log(e)))) )
}


stop("stop")

## Scale concentrations from micro- to nano-molar
ctrp.data$conc.nM <- ctrp.data$cpd_conc_umol * 10^3 

## Note that viabilities generally max out at 1. 
## i.e., they are fractional viabilities.
## Convert to percent viability.
ctrp.data$percent_viability <- 100 * (2^ctrp.data$bsub_value_log2)

source("../common/dss.R")
source("../common/drc-fit.R")

drug.screen.cols <- c("experiment_id", "master_cpd_id", "master_ccl_id")
conc.col <- "conc.nM"
response.col <- "percent_viability"
all.cols <- c(drug.screen.cols, conc.col, response.col)

cat("Fitting CTRPv2 using LL4\n")
ll4.fits <- fit.drc(ctrp.data[,all.cols], fct = "LL.4", drug.screen.cols = drug.screen.cols, conc.nM.col = conc.col, response.col = response.col, response.is.viability = TRUE)
cat("Done fitting CTRPv2 using LL4\n")
save.image(".Rdata")

