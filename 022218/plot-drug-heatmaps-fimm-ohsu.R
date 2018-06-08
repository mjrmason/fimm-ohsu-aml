rm(list = ls())

suppressPackageStartupMessages(library("foreach"))
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library(tidyr))

suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("gridExtra"))
suppressPackageStartupMessages(library("caTools"))
suppressPackageStartupMessages(library("data.table"))

suppressPackageStartupMessages(library(ggbeeswarm))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(synapseClient))
suppressPackageStartupMessages(library(mygene))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library("cowplot"))

library(rJava)
library(rcdk)
library(fingerprint)
library(plyr)
library(tidyverse)

library(dynamicTreeCut)
library(MCL)
library(clusteval)

library(WeightedCluster)
library(igraph)

library("dendextend")
## library("dendextendRcpp")

library(corrplot)

## Log in to Synapse
synapseLogin()

## Bring in common code
source("../common/data-preprocessing.R")
source("../common/models.R")
source("../common/plotting.R")
source("../common/utils.R")
source("heatmap-utils.R")

## These are defnied in heatmap-utils
orig.ohsu.metadata <- ohsu.metadata
orig.fimm.metadata <- fimm.metadata

## Register parallelization
num.cores <- detectCores()
if(!is.na(num.cores) && (num.cores > 1)) {
  suppressPackageStartupMessages(library("doMC"))
  cat(paste("Registering ", num.cores-1, " cores.\n", sep=""))
  registerDoMC(cores=(num.cores-1))
}
num.processes <- num.cores - 1

## Begin setup (this will need to be changed as the data sets in the intersection change)

source("fimm-ohsu-setup-120817.R")

rdata.file <- ".Rdata.model.mean.response.120817"

## Over-write some of the variables
## These fit files are assumed to be for AUCs/DSS values calculated across shared concentration ranges across data sets.
## These files may be created by a script such as
## calculate-dss-fimm-ohsu.R
## data.set.dss.fit.synIds <- c("syn11362885", "syn11362891")
## names(data.set.dss.fit.synIds) <- data.sets
## data.set.dss.fit.files <- c("ohsu-fimm-ohsu-outlier1-dss-t0.tsv", "fimm-fimm-ohsu-outlier1-dss-t0.tsv")

## A string that will be included in any output files.
file.prefix <- "fimm-ohsu-validate-trained-models-112217-one-ohsu-no-filtering-all-mean-response"

## End setup

## Load the map of drugs between data sets

source("process-drc-and-expr.R")
source("plot-drug-heatmap-utils.R")

## make.drug.heatmaps(orig.scaled.drcs)
make.drug.correlation.plots(orig.scaled.drcs, file.prefix = "fimm-ohsu-drug-correlations", multiple.plots.in.one.figure = TRUE)
## make.drug.correlation.plots(orig.scaled.drcs, file.prefix = NULL, multiple.plots.in.one.figure = TRUE)