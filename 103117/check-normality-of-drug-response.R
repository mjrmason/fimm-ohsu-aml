## rm(list = ls())

suppressPackageStartupMessages(library(glmnet))
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
suppressPackageStartupMessages(library("VennDiagram"))
suppressPackageStartupMessages(library("caret"))

suppressPackageStartupMessages(library(ggbeeswarm))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(synapseClient))

suppressPackageStartupMessages(library("GSVA"))
# use GSA to read in a gmt file
suppressPackageStartupMessages(library("GSA"))

## Bring in common code
source("../common/data-preprocessing.R")
source("../common/models.R")
source("../common/plotting.R")

## Register parallelization
num.cores <- detectCores()
if(!is.na(num.cores) && (num.cores > 1)) {
  suppressPackageStartupMessages(library("doMC"))
  cat(paste("Registering ", num.cores-1, " cores.\n", sep=""))
  registerDoMC(cores=(num.cores-1))
}
num.processes <- num.cores - 1

## Log in to Synapse
synapseLogin()

## Begin setup (this will need to be changed as the data sets in the intersection change)

## Purpose:
## This file models 2 data sets (OHSU and FIMM)

## Data sets to be analyzed:
data.sets <- c("ohsu", "fimm")
data.sets.to.analyze <- data.sets

## These fit files are assumed to be for AUCs/DSS values calculated across shared concentration ranges across data sets.
## These files may be created by a script such as
## calculate-dss-fimm-ohsu.R
data.set.dss.fit.synIds <- c("syn11362885", "syn11362891")
names(data.set.dss.fit.synIds) <- data.sets
data.set.dss.fit.files <- c("ohsu-fimm-ohsu-outlier1-dss-t0.tsv", "fimm-fimm-ohsu-outlier1-dss-t0.tsv")

## The columns within each respective data set that identifies the drug (in the name space of that data set)
data.set.drug.id.cols <- list("ohsu" = "inhibitor", "fimm" = "DRUG_ID")
patient.id.cols <- list("ohsu" = "patient_id", "fimm" = "PATIENT_ID")

## A string that will be included in any output files.
## Here "foc-aml-s-aml" = "fimm-ohsu-ctrp-aml-s-aml"
file.prefix <- "fimm-ohsu-drug-response-normality"

## End setup

cat(paste0("Check normality of drug response in data sets: ", paste(data.sets, collapse=","), "\n"))

## Load in the fit files
fits <- llply(data.set.dss.fit.synIds, .parallel = FALSE,
              .fun = function(synId) {
                       obj <- synGet(synId, downloadFile = TRUE)
                       file <- getFileLocation(obj)
                       read.table(file, sep="\t", header=TRUE, as.is=TRUE, comment.char="", quote="\"")
              })

responses <- c("dss.auc.ll4", "dss.auc.l4")

p.val.cutoff <- 0.05

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

tbl <- c()
for(response in responses) {
  for(ds in data.sets) {

    transform <- "none"
    file <- paste0(paste(file.prefix, ds, response, transform, "qq", sep="-"), ".pdf")
    pdf(file)
    l <- dlply(fits[[ds]], .variables = data.set.drug.id.cols[[ds]], .parallel = FALSE,
                 .fun = function(df) {
                          vals <- df[, response]
                          names(vals) <- df[, patient.id.cols[[ds]]]
                          vals <- vals[!is.na(vals)]
                          outliers <- my.qqnorm(vals, main = shapiro.test(vals)$p.value)
                          qqline(vals)
                          list("patients" = names(vals), "outliers" = outliers)
                        })
    d <- dev.off()

    file <- paste0(paste(file.prefix, ds, response, transform, "pt-outlier", sep="-"), ".pdf")
    pdf(file)
    pt.tally <- table(unlist(lapply(l, function(f) f$patients)))
    tmp <- table(unlist(lapply(l, function(f) f$outliers)))
    outlier.tally <- rep(0, length(pt.tally))
    names(outlier.tally) <- names(pt.tally)
    outlier.tally[names(tmp)] <- tmp

    outlier.ratio <- as.numeric(outlier.tally / pt.tally[names(outlier.tally)])
    df <- data.frame(ratio = outlier.ratio, patient = names(outlier.tally))
    g <- ggplot(data = df, aes(x = patient, y = ratio))
    g <- g + geom_point()
    g <- g + ylab("fraction outliers")
    g <- g + theme(axis.text.x = element_text(angle = 90))
    print(g)
    d <- dev.off()

    res <- ddply(fits[[ds]], .variables = data.set.drug.id.cols[[ds]], 
                 .fun = function(df) {
                          flag <- !is.na(df[, response])
                          df <- df[flag, ]
                          vals <- df[, response]
                          df$rank <- rank(vals)/length(vals)
                          df$neg <- vals < 0
                          df
                 })
    plot(res$rank, res$gof.ll4)
    g <- ggplot(data = res)
    g <- g + geom_boxplot(aes(x = neg, y = gof.ll4))
    print(g)

    plot(density(res$gof.ll4))
    plot(density(res$gof.ll4[res$neg == FALSE]))
    plot(density(res$gof.ll4[res$neg == TRUE]))

    res <- dlply(fits[[ds]], .variables = data.set.drug.id.cols[[ds]], 
                 .fun = function(df) {
                          vals <- df[, response]
                          vals <- vals[!is.na(vals)]
                          vals <- vals[vals > 0]
                          shapiro.test(vals)$p.value
                 })
    num.non.sig <- length(which(unlist(res) > p.val.cutoff))
    num.tot <- length(res)
    frac <- num.non.sig / num.tot
    tbl <- rbind(tbl, c(ds, response, transform, num.non.sig, num.tot, frac))
    
    ## Shift and log transform
    res <- ddply(fits[[ds]], .variables = data.set.drug.id.cols[[ds]], .parallel = FALSE,
                 .fun = function(df) {
                          flag <- !is.na(df[, response]) & (df[, response] > 0)
                          vals <- df[flag, response]
                          patient <- df[flag, patient.id.cols[[ds]]]
                          shift <- 0
                          if(any(vals < 0)) {
                            shift <- 0.001 + -1 * min(vals, na.rm=TRUE)
##			    shift <- -1 * min(vals, na.rm=TRUE)
##                            shift <- shift + median(vals + shift, na.rm = TRUE)
                          }
                          data.frame(patient = patient, transformed = log(vals + shift))
                        })
    transform <- "shift-and-log"
    file <- paste0(paste(file.prefix, ds, response, transform, "qq", sep="-"), ".pdf")
    pdf(file)
    l <- dlply(res, .variables = data.set.drug.id.cols[[ds]], .parallel = FALSE,
                 .fun = function(df) {
                          vals <- df[, "transformed"]
                          names(vals) <- df[, "patient"]
                          vals <- vals[!is.na(vals)]
                          outliers <- my.qqnorm(vals, main = shapiro.test(vals)$p.value)
                          qqline(vals)
                          list("patients" = names(vals), "outliers" = outliers)
                        })
    d <- dev.off()

    file <- paste0(paste(file.prefix, ds, response, transform, "pt-outlier", sep="-"), ".pdf")
    pdf(file)
    pt.tally <- table(unlist(lapply(l, function(f) f$patients)))
    tmp <- table(unlist(lapply(l, function(f) f$outliers)))
    outlier.tally <- rep(0, length(pt.tally))
    names(outlier.tally) <- names(pt.tally)
    outlier.tally[names(tmp)] <- tmp

    outlier.ratio <- as.numeric(outlier.tally / pt.tally[names(outlier.tally)])
    df <- data.frame(ratio = outlier.ratio, patient = names(outlier.tally))
    g <- ggplot(data = df, aes(x = patient, y = ratio))
    g <- g + geom_point()
    g <- g + ylab("fraction outliers")
    g <- g + theme(axis.text.x = element_text(angle = 90))
    print(g)
    d <- dev.off()

    res <- dlply(res, .variables = data.set.drug.id.cols[[ds]], 
                 .fun = function(df) {
                          vals <- df[, "transformed"]
                          vals <- vals[!is.na(vals)]
                          shapiro.test(vals)$p.value
                 })
    num.non.sig <- length(which(unlist(res) > p.val.cutoff))
    num.tot <- length(res)
    frac <- num.non.sig / num.tot
    tbl <- rbind(tbl, c(ds, response, transform, num.non.sig, num.tot, frac))
    
    ## Shift and probit (qnorm) transform
    res <- ddply(fits[[ds]], .variables = data.set.drug.id.cols[[ds]], .parallel = FALSE,
                 .fun = function(df) {
                          flag <- !is.na(df[, response]) & (df[, response] > 0)
                          vals <- df[flag, response]
                          patient <- df[flag, patient.id.cols[[ds]]]
                          shift <- 0
                          if(any(vals < 0)) {
                            shift <- 0.001 + -1 * min(vals, na.rm=TRUE)
##			    shift <- -1 * min(vals, na.rm=TRUE)
##                            shift <- shift + median(vals + shift, na.rm = TRUE)

                          }
                          data.frame(patient = patient, transformed = qnorm((vals + shift)/(max(vals + shift + 0.001))))
                        })
    transform <- "shift-and-probit"
    file <- paste0(paste(file.prefix, ds, response, transform, "qq", sep="-"), ".pdf")
    pdf(file)
    l <- dlply(res, .variables = data.set.drug.id.cols[[ds]], .parallel = FALSE,
                 .fun = function(df) {
                          vals <- df[, "transformed"]
                          names(vals) <- df[, "patient"]
                          vals <- vals[!is.na(vals)]
                          outliers <- my.qqnorm(vals, main = shapiro.test(vals)$p.value)
                          qqline(vals)
                          list("patients" = names(vals), "outliers" = outliers)
                        })
    d <- dev.off()

    file <- paste0(paste(file.prefix, ds, response, transform, "pt-outlier", sep="-"), ".pdf")
    pdf(file)
    pt.tally <- table(unlist(lapply(l, function(f) f$patients)))
    tmp <- table(unlist(lapply(l, function(f) f$outliers)))
    outlier.tally <- rep(0, length(pt.tally))
    names(outlier.tally) <- names(pt.tally)
    outlier.tally[names(tmp)] <- tmp

    outlier.ratio <- as.numeric(outlier.tally / pt.tally[names(outlier.tally)])
    df <- data.frame(ratio = outlier.ratio, patient = names(outlier.tally))
    g <- ggplot(data = df, aes(x = patient, y = ratio))
    g <- g + geom_point()
    g <- g + ylab("fraction outliers")
    g <- g + theme(axis.text.x = element_text(angle = 90))
    print(g)
    d <- dev.off()

    res <- dlply(res, .variables = data.set.drug.id.cols[[ds]], 
                 .fun = function(df) {
                          vals <- df[, "transformed"]
                          vals <- vals[!is.na(vals)]
                          shapiro.test(vals)$p.value
                 })
    num.non.sig <- length(which(unlist(res) > p.val.cutoff))
    num.tot <- length(res)
    frac <- num.non.sig / num.tot
    tbl <- rbind(tbl, c(ds, response, transform, num.non.sig, num.tot, frac))
    
    ## Shift and logit transform
    library(car)
    res <- ddply(fits[[ds]], .variables = data.set.drug.id.cols[[ds]], .parallel = FALSE,
                 .fun = function(df) {
                          flag <- !is.na(df[, response]) & (df[, response] > 0)
                          vals <- df[flag, response]
                          patient <- df[flag, patient.id.cols[[ds]]]
                          shift <- 0
                          if(any(vals < 0)) {
                            shift <- 0.001 + -1 * min(vals, na.rm=TRUE)
##			    shift <- -1 * min(vals, na.rm=TRUE)
##                            shift <- shift + median(vals + shift, na.rm = TRUE)
                          }
                          data.frame(patient = patient, transformed = logit((vals + shift)/(max(vals + shift + 0.001))))
                        })
    transform <- "shift-and-logit"
    file <- paste0(paste(file.prefix, ds, response, transform, "qq", sep="-"), ".pdf")
    pdf(file)
    l <- dlply(res, .variables = data.set.drug.id.cols[[ds]], .parallel = FALSE,
                 .fun = function(df) {
                          vals <- df[, "transformed"]
                          names(vals) <- df[, "patient"]
                          vals <- vals[!is.na(vals)]
                          outliers <- my.qqnorm(vals, main = shapiro.test(vals)$p.value)
                          qqline(vals)
                          list("patients" = names(vals), "outliers" = outliers)
                        })
    d <- dev.off()

    file <- paste0(paste(file.prefix, ds, response, transform, "pt-outlier", sep="-"), ".pdf")
    pdf(file)
    pt.tally <- table(unlist(lapply(l, function(f) f$patients)))
    tmp <- table(unlist(lapply(l, function(f) f$outliers)))
    outlier.tally <- rep(0, length(pt.tally))
    names(outlier.tally) <- names(pt.tally)
    outlier.tally[names(tmp)] <- tmp

    outlier.ratio <- as.numeric(outlier.tally / pt.tally[names(outlier.tally)])
    df <- data.frame(ratio = outlier.ratio, patient = names(outlier.tally))
    g <- ggplot(data = df, aes(x = patient, y = ratio))
    g <- g + geom_point()
    g <- g + ylab("fraction outliers")
    g <- g + theme(axis.text.x = element_text(angle = 90))
    print(g)
    d <- dev.off()

    res <- dlply(res, .variables = data.set.drug.id.cols[[ds]], 
                 .fun = function(df) {
                          vals <- df[, "transformed"]
                          vals <- vals[!is.na(vals)]
                          shapiro.test(vals)$p.value
                 })
    num.non.sig <- length(which(unlist(res) > p.val.cutoff))
    num.tot <- length(res)
    frac <- num.non.sig / num.tot
    tbl <- rbind(tbl, c(ds, response, transform, num.non.sig, num.tot, frac))
  }
}
colnames(tbl) <- c("data.set", "response", "transform", "num.norm", "num.tot", "frac.norm")
tbl <- as.data.frame(tbl)

tbl$frac.norm <- as.numeric(format(as.numeric(as.character(tbl$frac.norm)), digits=2))
file <- paste0(paste(file.prefix, "table", sep="-"), ".pdf")
pdf(file)
plot.table(tbl, main = "Normal Drug Response (Shapiro Test)")
d <- dev.off()