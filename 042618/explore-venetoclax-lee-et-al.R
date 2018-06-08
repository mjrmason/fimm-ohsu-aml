library(ggplot2)
library(plyr)
library(glmnet)
library(randomForest)
suppressPackageStartupMessages(library("cowplot"))


source("../common/utils.R")
## venetoclax.feature.syms <- c("BCL2", "CD14", "FLT3", "MAFB", "LRP1", "CREG1", "LOC200772")
venetoclax.feature.syms <- c("BCL2", "CD14", "FLT3", "MAFB", "LRP1")
venetoclax.feature.gene.tbl <- symbols.to.ensg.mapping(venetoclax.feature.syms)

