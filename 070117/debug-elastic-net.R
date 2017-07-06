library(synapseClient)
library(data.table)
library(ggplot2)
library(gridExtra)

suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(caret))
suppressPackageStartupMessages(library(glmnet))
library( ReporteRs )
library(scales)

suppressPackageStartupMessages(library("openxlsx"))

synapseLogin()

remove.outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

synId <- "syn8488824"
obj <- synGet(id=synId, downloadFile = TRUE)
file <- getFileLocation(obj)
fimm.raw.dr <- read.xlsx(file)

## path <- "ohsu.dss.t0.tsv"
synId <- "syn10083494"
obj <- synGet(synId, downloadFile = TRUE)
file <- getFileLocation(obj)
ohsu.dss.t0 <- as.data.frame(fread(file))
## ohsu.dss.t0$ic50.ll4 <- remove.outliers(ohsu.dss.t0$e.ll4)
## ohsu.dss.t0$ic50.l4 <- remove.outliers(ohsu.dss.t0$e.l4)

synId <- "syn10083488"
obj <- synGet(synId, downloadFile = TRUE)
file <- getFileLocation(obj)
fimm.dss.t0 <- as.data.frame(fread(file))

bad.drug <- "FIMM100365"
good.drug <- "FIMM100364"
## head(fp[,c("alpha","drug","pearson.cor")])
## alpha       drug pearson.cor
## 50      0 FIMM100365  -0.6040647
## 22      0 FIMM003751  -0.5379798
## 104  0.25 FIMM023795  -0.5315151
## 42      0 FIMM003793  -0.5291159
## 154   0.5 FIMM003782  -0.5118911
## 95   0.25 FIMM003782  -0.5052962

##  tail(fp[,c("alpha","drug","pearson.cor")])
##alpha       drug pearson.cor
##173   0.5 FIMM100395   0.6429484
##68   0.25 FIMM003711   0.6449133
##114  0.25 FIMM100395   0.6602826
###9       0 FIMM003711   0.6622430
##48      0 FIMM100362   0.6666116
##49      0 FIMM100364   0.6767241

fp <- read.table(file="fimm.pp.actual.predicted.correlations.tsv", sep="\t", header=TRUE)
short.df <- read.table(file="fimm.pp.short.df", sep="\t", header=TRUE)

load(".Rdata.fimm")

fp <- unique(fimm.pp$actual.predicted.correlations.df)
short.df <- unique(fimm.pp$short.df)

library(tidyr)

my.dup <- function(df) { duplicated(df, fromLast=TRUE) | duplicated(df, fromLast=FALSE) }

short.df <- unique(short.df)

dups <- short.df[my.dup(short.df$sample) & short.df$drug == "FIMM000160" & short.df$alpha == 0,]
dups[order(dups$sample),]


bad.drug.tbl <- subset(fp, pearson.cor < 0 & pearson.p < 0.01 & alpha == 0)
good.drug.tbl <- subset(fp, pearson.cor > 0 & pearson.p < 0.01 & alpha == 0)
bad.drugs <- bad.drug.tbl$drug
good.drugs <- good.drug.tbl$drug

pdf("bad-drugs.pdf")
for(target.drug in unique(bad.drugs)) {
  print(target.drug)
  sb <- subset(short.df, (drug == target.drug) & (alpha == 0))
  fm <- subset(fimm.dss.t0, DRUG_ID == target.drug)
  print(table(fm$b.l4 < 0))
  sb.pred <- unique(sb[sb$response.type == "predicted",c("sample", "response")])
  colnames(sb.pred) <- c("sample", "predicted.response")
  sb.resp <- unique(sb[sb$response.type == "actual",c("sample", "response")])
  colnames(sb.resp) <- c("sample", "actual.response")
  m <- merge(sb.pred, sb.resp, by="sample")
  plot(m$actual.response, m$predicted.response, main=target.drug, xlab = "Actual Response", ylab = "Median Predicted Response")
##  o <- m
##  o$actual.response <- remove.outliers(o$actual.response)
##  o <- na.omit(o)
##  plot(o$actual.response, o$predicted.response)
  plot(density(fm$gof.l4))
}
d <- dev.off()

pdf("good-drugs.pdf")
good.drugs <- unique(good.drugs)
for(target.drug in unique(good.drugs)) {
  print(target.drug)
  sb <- subset(short.df, drug == target.drug & alpha == 0)
  fm <- subset(fimm.dss.t0, DRUG_ID == target.drug)
  print(table(fm$b.l4 < 0))
  sb.pred <- unique(sb[sb$response.type == "predicted",c("sample", "response")])
  colnames(sb.pred) <- c("sample", "predicted.response")
  sb.resp <- unique(sb[sb$response.type == "actual",c("sample", "response")])
  colnames(sb.resp) <- c("sample", "actual.response")
  m <- merge(sb.pred, sb.resp, by="sample")
  plot(m$actual.response, m$predicted.response, main=target.drug, xlab = "Actual Response", ylab = "Median Predicted Response")
  plot(density(fm$gof.l4))
}
d <- dev.off()


plot(density(fimm.dss.t0$gof.l4[fimm.dss.t0$b.l4 < 0]), col="blue")
lines(density(fimm.dss.t0$gof.l4[fimm.dss.t0$b.l4 > 0]))

plot(density(ohsu.dss.t0$gof.l4[ohsu.dss.t0$b.l4 < 0]), col="blue")
lines(density(ohsu.dss.t0$gof.l4[ohsu.dss.t0$b.l4 > 0]))


## Look at outliers and max negative response

library(plyr)
bad.drug.stats <- ddply(short.df[(short.df$response.type == "actual") & (short.df$drug %in% bad.drugs),], .variables = c("drug"),
                        .fun = function(df) max(df$response)[1])
good.drug.stats <- ddply(short.df[(short.df$response.type == "actual") & (short.df$drug %in% good.drugs),], .variables = c("drug"),
                        .fun = function(df) max(df$response)[1])
plot(density(bad.drug.stats$V1))
lines(density(good.drug.stats$V1), col="blue")

unique(subset(short.df, drug == "FIMM003782" & alpha == 0 & response.type == "actual"))

unq <- unique(subset(short.df, response < 0 & alpha == 0 & response.type == "actual"))

screen.id <- "FHRB_560_02052012_0800_BM"
sub <- subset(fimm.raw.dr, DRUG_ID == "FIMM003782" & SCREEN_ID == screen.id)
plot(sub$CONCENTRATION, sub$PERCENT_INHIBITION)
plot(log10(sub$CONCENTRATION), sub$PERCENT_INHIBITION)


unq <- unique(subset(short.df, response > 0 & alpha == 0 & response.type == "actual"))
indx <- 5
drug.id <- unq$drug[indx]
screen.id <- unq$sample[indx]
sub <- subset(fimm.raw.dr, DRUG_ID == "FIMM003782" & SCREEN_ID == screen.id)
plot(log10(sub$CONCENTRATION), sub$PERCENT_INHIBITION)


screen.id <- "FHRB_370_24032015_9999_BM"
screen.id <- "FHRB_784_16052011_9999_BM" 

bad.drug.stats <- ddply(short.df[(short.df$response.type == "actual") & (short.df$drug %in% bad.drugs),], .variables = c("drug"),
                        .fun = function(df) min(df$response)[1])
good.drug.stats <- ddply(short.df[(short.df$response.type == "actual") & (short.df$drug %in% good.drugs),], .variables = c("drug"),
                         .fun = function(df) min(df$response)[1])
plot(density(good.drug.stats$V1), col="blue")
lines(density(bad.drug.stats$V1))

bad.drug.stats <- ddply(short.df[(short.df$response.type == "actual") & (short.df$drug %in% bad.drugs),], .variables = c("drug"),
                        .fun = function(df) max(df$response)[1] - min(df$response)[1])
good.drug.stats <- ddply(short.df[(short.df$response.type == "actual") & (short.df$drug %in% good.drugs),], .variables = c("drug"),
                         .fun = function(df) max(df$response)[1] - min(df$response)[1])
plot(density(bad.drug.stats$V1))
lines(density(good.drug.stats$V1), col="blue")

bad.drug.stats <- ddply(short.df[(short.df$alpha == 0) & (short.df$response.type == "actual") & (short.df$drug %in% bad.drugs),], .variables = c("drug"),
                        .fun = function(df) dim(df))
good.drug.stats <- ddply(short.df[(short.df$alpha == 0) & (short.df$response.type == "actual") & (short.df$drug %in% good.drugs),], .variables = c("drug"),
                         .fun = function(df) dim(df))
plot(density(bad.drug.stats$V1))
lines(density(good.drug.stats$V1), col="blue")

bad.drugs <- c("FIMM100365", "FIMM003751", "FIMM023795", "FIMM003793")

bad.drugs <- bad.drug.tbl$drug
good.drugs <- good.drug.tbl$drug

from=-5*10^5
to=5*10^5
plot(density(subset(fimm.dss.t0, DRUG_ID %in% bad.drugs)$auc.l4, from=from, to=to))
lines(density(subset(fimm.dss.t0, DRUG_ID %in% good.drugs)$auc.l4, from=from,to=to))
plot(density(subset(fimm.dss.t0, DRUG_ID %in% good.drugs)$gof.l4))
lines(density(subset(fimm.dss.t0, DRUG_ID %in% bad.drugs)$gof.l4))
