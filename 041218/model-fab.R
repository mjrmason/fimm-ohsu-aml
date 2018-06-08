source("model-fab-local-init.R")

my.dup <- function(x) duplicated(x, fromLast = TRUE) | duplicated(x, fromLast = FALSE)

ohsu.patient.fab <- unique(ohsu.metadata[, c("patient_id", "fab")])

## One patient has inconsistent FABs
ohsu.patient.fab[my.dup(ohsu.patient.fab$patient_id),]

ohsu.patient.fab <- subset(ohsu.patient.fab, !is.na(fab))
ohsu.patient.fab <- subset(ohsu.patient.fab, !is.na(patient_id))
ohsu.patient.fab <- subset(ohsu.patient.fab, !my.dup(patient_id))
rownames(ohsu.patient.fab) <- ohsu.patient.fab$patient_id
ohsu.patient.fab$fab <- factor(ohsu.patient.fab$fab)
fab <- ohsu.patient.fab$fab
names(fab) <- rownames(ohsu.patient.fab)

common.samples <- intersect(rownames(ohsu.patient.fab), colnames(orig.drcs[["ohsu"]]))

ohsu.scaled.drcs <- t(scale(t(orig.drcs[["ohsu"]])))

ohsu.glds <- colMeans(t(scale(t(orig.drcs[["ohsu"]][, common.samples]))), na.rm=TRUE)

## Look for correlation between mean response and FAB
fit <- lm(as.numeric(ohsu.glds[common.samples]) ~ fab[common.samples])
summary(fit)

## Look for correlation across drugs and FAB
library(dplyr)
library(plyr)
mat <- ohsu.scaled.drcs
indices <- 1:nrow(mat)
names(indices) <- rownames(mat)
fits <- llply(indices, .parallel = FALSE,
              .fun = function(i) {
                  print(rownames(mat)[i])
                  cols <- !is.na(as.numeric(mat[i, ])) & (colnames(mat) %in% names(fab))
                  ## Remove intercept
                  fit <- tryCatch({lm(as.numeric(mat[i, cols]) ~ 0 + fab[cols])}, error = function(e) NULL)
                  ## fit <- tryCatch({lm(as.numeric(mat[i, cols]) ~  fab[cols])}, error = function(e) NULL)
                  fit
              })

## Try this without an intercept--M0 is often different from others

## Multinomial logistic regression
library(nnet)

flag <- !my.dup(ohsu.metadata$lab_id)
fab.sample <- ohsu.metadata[flag, "fab"]
names(fab.sample) <- ohsu.metadata[flag, "SeqID"]
fab.sample <- fab.sample[!is.na(fab.sample)]
fab.sample <- as.factor(fab.sample)

expr.mat <- orig.exprs[["ohsu"]]
colnames(expr.mat) <- gsub("X(.*)", "\\1", colnames(expr.mat))
colnames(expr.mat) <- gsub("\\.", "-", colnames(expr.mat))
cs <- intersect(names(fab.sample), colnames(expr.mat))
fab.sample <- fab.sample[cs]
expr.mat <- expr.mat[, cs]

df <- as.data.frame(t(expr.mat))
df$fab <- fab.sample

stop("stop")

mod <- multinom(fab ~ ., data = df)

library(randomForest)
rf <- randomForest(fab ~ ., data = df)


llply(fits, .fun = function(fit) summary(fit))
## Trametinib
## SNS-032 (BMS-387032)
## Palbociclib
## Neratinib (HKI-272)
## MGCD-265
## LY-333531
## Lapatinib
## INK-128
## Doramapimod (BIRB 796)
## Crizotinib (PF-2341066)
## BEZ235
## Azacytidine
## AGI-5198
## Afatinib (BIBW-2992)

