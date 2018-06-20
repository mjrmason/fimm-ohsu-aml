suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(synapseClient))
suppressPackageStartupMessages(library("foreach"))
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("ggplot2"))
##suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("MonoInc"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library(tidyr))

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

cat("Use DRC to fit L.4 (4-parameter logistic) and LL.4 (4-paramter log logistic) functions to FIMM data CM.\n")

## BEGIN setup

## plan
## what is concentration range?
## how often do we drop the first point? last point?
## make sure we drop these outside of range from auc calculation

source("fimm-cm-setup-050318.R")

prefix <- "fimm-cm"
data.set <- "fimm-cm"

conc.col <- conc.cols[[data.set]]
response.col <- response.cols[[data.set]]
response.is.viability <- responses.are.viabilities[[data.set]]
drug.col <- data.set.drug.id.cols[[data.set]]
drug.screen.cols <- screen.id.cols[[data.set]]
experimental.factors <- experimental.factor.cols[[data.set]]
parentId <- processed.folder.synIds[[data.set]]

## Read in FIMM CM raw drug response data (FIMM_DSRT_Data_CM.txt)
synId <- raw.drug.synIds[[data.set]]
obj <- synGet(id=synId, downloadFile = TRUE)
file <- getFileLocation(obj)
## These data have the same columns and concentration ranges (hence in nM) as the
## FIMM MCM data ("syn8488824").  Just drop the extraneous "X" column.
## And "DSRT_SET" in these CM data correspond to "DRUG_SET" in MCM data.
drug.data.long <- read.table(file, sep="\t", header=TRUE, as.is=TRUE)
colnames(drug.data.long)[colnames(drug.data.long) == "DSRT_SET"] <- "DRUG_SET"
drug.data.long <- drug.data.long[, !(colnames(drug.data.long) == "X")]

## Read in the expression data
## synId <- data.set.orig.expr.synIds[[data.set]]
synId <- "syn12498090"
obj <- synGet(id=synId, downloadFile = TRUE)
file <- getFileLocation(obj)
## expr <- fread(file)
expr <- read.table(file, sep="\t", header=TRUE, as.is=TRUE)

## Load FIMM metadata (including relapse/refractory/diagnosis)
metadata <- NULL
if(!is.null(metadata.synIds[[data.set]])) { ## Don't yet have metadata for CM
  obj <- synGet(id=metadata.synIds[[data.set]], downloadFile = TRUE)
  file <- getFileLocation(obj)
  metadata <- openxlsx:::read.xlsx(file, sheet = 1)
}

## Subset to samples that are listed in metadata as CM, but are _not_ in the MCM data.
## This ensures that the valiation is completely orthogonal -- not just different
## media (CM vs MCM), but different patients.
if(is.null(metadata)) { die("Was not expecting NULL metadata\n") }

## Define the CM patients as those in our CM data
tmp <- unique(drug.data.long[, c("SCREEN_ID", "PATIENT_ID")])
## Ensure that they all associated SCREEN_IDs with _CM suffixes
if(!(all(grepl(tmp$SCREEN_ID, pattern="_CM$")))) {
  stop("Some CM patients have a SCREEN_ID without a _CM postfix\n")
}

cm.patients <- unique(tmp$PATIENT_ID)

## Now, exclude patients that we would have used in our MCM analysis.
source("../common/utils.R")
source("../common/data-preprocessing.R")
source("../common/process-drc-and-expr.R")

ds <- "fimm"
l <- prepare.drug.response.and.expr.matrices(orig.fits[[ds]], orig.exprs[[ds]], drugs = orig.drug.name.tbl[, drug.name.col],
                                             drug.col = drug.name.col, patient.col = expr.patient.id.cols[[ds]],
                                             response.col = "auc.ll4")
mcm.drug.screens.and.pts <- unique(orig.fits[["fimm"]][, c("SCREEN_ID", "PATIENT_ID")])
mcm.expr.screens.and.pts <- mcm.drug.screens.and.pts[mcm.drug.screens.and.pts$SCREEN_ID %in% colnames(exprs[["fimm"]]),]
mcm.drug.pts <- unique(mcm.drug.screens.and.pts$PATIENT_ID)
mcm.expr.pts <- unique(mcm.expr.screens.and.pts$PATIENT_ID)
mcm.pts <- intersect(mcm.drug.pts, mcm.expr.pts)
mcm.screens.and.pts <- subset(mcm.drug.screens.and.pts, PATIENT_ID %in% mcm.pts)

uniquely.cm.patients <- cm.patients[!(cm.patients %in% mcm.pts)]

## drug.data.long has one concentration/response per row.
## NB: concentrations are in nM.  Let's make that explicit.
drug.data.long$conc.nM <- drug.data.long$CONCENTRATION

## NB: PERCENT_INHBITION is ... inhibition!
## NB: PERCENT_INHIBITION really is a percent -- it runs from 0 to 100 (well, actually, slightly below and above).

all.drugs <- unique(drug.data.long[, drug.col])
drugs.to.exclude <- all.drugs[!grepl(all.drugs, pattern="FIMM")]
cat(paste0("Excluding the following (control?) drugs from outlier detection: ", paste(drugs.to.exclude, collapse=", "), "\n"))

all.cols <- c(drug.screen.cols, conc.col, response.col)

## Subset the drug data to those patients were only assayed in CM media -- i.e., were not
## used in training our models with MCM media

orig.drug.data.long <- drug.data.long
drug.data.long <- subset(drug.data.long, PATIENT_ID %in% uniquely.cm.patients)

## Limit fo AML
## We assume all of these samples are from AML

output.file <- paste0("fimm.cm.l4.and.ll4.fits.050318.tsv")

dir.create(output.dir)

## Restrict to samples for which we have both expression and drug response
common.samples <- intersect(drug.data.long[, expr.patient.id.cols[[data.set]]], colnames(expr))
drug.data.long <- drug.data.long[drug.data.long[, expr.patient.id.cols[[data.set]]] %in% common.samples, ]

## END setup

## Plot distribution of these control drugs
for(drug.id in drugs.to.exclude) {
  sub <- subset(drug.data.long, DRUG_ID == drug.id)
  if(nrow(sub) > 0) {
    pdf(paste0(output.dir, "/", "control-", drug.id, "-distribution.pdf"))
    cat(paste0(drug.id, " has ", nrow(sub), " entries\n"))
    plot(density(sub$PERCENT_INHIBITION), main = drug.id)
    d <- dev.off()
  }
}

## Confirm that the LOWn (n = 1 to 10) and HIGHn are specific to one patient.  So let's drop those.
low.high.drugs <- all.drugs[grepl(pattern="LOW", all.drugs) | grepl(pattern="HIGH", all.drugs)]
low.high.drugs <- low.high.drugs[!(low.high.drugs %in% c("LOW", "HIGH"))]

pts <- unique(subset(drug.data.long, DRUG_ID %in% low.high.drugs)$PATIENT_ID)
cat(paste0("Patients with LOWn/HIGHn controls: ", paste(pts, collapse=", "), "\n"))
if(length(pts) > 1) {
  stop("Was expecting at most one patient with LOWn/HIGHn controls\n")
}

ctrl.drugs <- drugs.to.exclude[!(drugs.to.exclude %in% low.high.drugs)]
ctrl.data <- subset(drug.data.long, DRUG_ID %in% drugs.to.exclude)

## Is there a patient, screen, or plate-specific affect on the control drugs?

## From Swapnil's email on May 11, 2018, we expect that there may be a batch effect in the BLANK or EMPTY
## wells because students sometimes did or did not add CTG:
## The ‘Sometimes’ tag comes because in the initial days students would sometimes skip adding CTG in these wells because they are blank anyway. 
## The FO4 layouts however should be more consistent.

## So, let's not test BLANK or EMPTY.

exclude <- list()
ctrl.drugs <- ctrl.drugs[!(ctrl.drugs %in% c("BLANK", "EMPTY"))]
for(ctrl in ctrl.drugs) {
  exclude[[ctrl]] <- list()
  df <- subset(drug.data.long, DRUG_ID == ctrl)
  covariates <- drug.screen.cols[!(drug.screen.cols %in% c("DRUG_ID", "DRUG_NAME"))]
  ## For BLANK and EMPTY, don't include SCREEN_ID or PATIENT_ID as covariates.
  ## There are no cells in the wells for BLANK or EMPTY, so these covariates do not make sense.
  if(ctrl %in% c("BLANK", "EMPTY")) {
    covariates <- covariates[!(covariates %in% c("SCREEN_ID", "PATIENT_ID"))]
  }
  keep <- unlist(lapply(covariates, function(covariate) length(unique(df[, covariate])) > 1))
  if(!any(keep)) { next }
  covariates <- covariates[keep]
  form <- paste(response.col, "~ 1. + ", paste(covariates, collapse = " + "))
  lmo <- lm(formula = form, data = df)
  sm <- summary(lmo)
  pval <- 1 - pf(sm$fstatistic[1], sm$fstatistic[2], sm$fstatistic[3]) 
  cat(paste0("Testing ", form, " for DRUG_ID = ", ctrl, "; pval = ", pval, "\n"))
  if(pval < 0.05) {
    coefs <- coefficients(sm)
    sig.flag <- !is.na(coefs[,4]) & (coefs[,4] < 0.05)
    if(any(sig.flag)) {
##      write.table(coefs[sig.flag, , drop = FALSE], sep="\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
      print(coefs[sig.flag, , drop = FALSE])
      sig.coef.names <- rownames(coefs)[sig.flag]
      for(var in covariates) {
        var.flag <- grepl(pattern = var, sig.coef.names)
        if(!(any(var.flag))) { next }
        sig.levels <- unlist(lapply(sig.coef.names[var.flag], function(str) gsub(str, pattern=paste0("^", var), replacement="")))
        exclude[[ctrl]][[var]] <- sig.levels
      }
    }
  }
  cat("\n")
}

for(ctrl in names(exclude)) {
  for(var in names(exclude[[ctrl]])) {
    flag <- drug.data.long[, var] %in% exclude[[ctrl]][[var]]
    cat(paste0("Excluding the following ", var, " because of ", ctrl, " (n = ", length(which(flag)), 
               " rows): ", paste(exclude[[ctrl]][[var]], collapse=", "), "\n"))
    exclude.col <- paste0(ctrl, ".", var, ".exclude")
    drug.data.long[, exclude.col] <- FALSE
    drug.data.long[flag, exclude.col] <- TRUE
  }
}

## Define the acceptable response range at the minimum concentration (min.response.at.min.conc, max.response.at.min.conc)

source("../common/annotate-drug-response-outliers.R")

## Define the distribution of the response at the minimum concentration
flag <- !(drug.data.long[, drug.col] %in% drugs.to.exclude)
drug.data.long <- drug.data.long[flag, ]
res <- define.response.at.min.conc.distribution(drug.data.long, drug.screen.cols, conc.col, response.col)
responses.at.min.conc <- res[["responses.at.min.conc"]]$PERCENT_INHIBITION
## min.conc.mean <- res[["mean"]]
## min.conc.sd <- res[["sd"]]
## min.response.at.min.conc <- min.conc.mean - 3 * min.conc.sd
## max.response.at.min.conc <- min.conc.mean + 3 * min.conc.sd

## Define the distribution of the response of the "HIGH" controls
high.responses <- subset(ctrl.data, DRUG_ID == "HIGH")$PERCENT_INHIBITION
lst <- list("Response at min conc" = responses.at.min.conc, "HIGH" = high.responses)

## Define the distribution of the response of the "DMSO" controls
if(any(ctrl.data$DRUG_ID == "DMSO")) {
  dmso.responses <- subset(ctrl.data, DRUG_ID == "DMSO")$PERCENT_INHIBITION
##  dmso.responses <- subset(orig.drug.data.long, DRUG_ID == "DMSO")$PERCENT_INHIBITION
  lst[["DMSO"]] <- dmso.responses
}

indices <- 1:length(lst)
names(indices) <- names(lst)
file <- "min-conc-distributions-fimm-cm.pdf"
pdf(paste0(output.dir, "/", file))
df <- Reduce("rbind", llply(indices, .fun = function(i) data.frame(color = names(lst)[i], x = lst[[i]])))
g <- ggplot(data = df, aes(x = x, color = color))
g <- g + geom_density()
g <- g + xlab("PERCENT_INHIBITION")
g <- g + ggtitle("FIMM CM") + theme(legend.title=element_blank())
g
d <- dev.off()
l_ply(indices, 
      .fun = function(i) {
         qnt <- quantile(lst[[i]], probs = c(0.05, 0.95))
         cat(paste0(names(lst)[i], " mean = ", mean(lst[[i]]), " sd = ", sd(lst[[i]]), " ",
             paste(unlist(lapply(1:length(qnt), function(i) paste0(names(qnt)[i], ": ", unname(qnt[i])))), collapse = " "), "\n"))

      })

## Use the HIGH controls to define the distribution at the minimum concentration
min.conc.mean <- mean(lst[["HIGH"]])
min.conc.sd <- sd(lst[["HIGH"]])
min.response.at.min.conc <- min.conc.mean - 3 * min.conc.sd
max.response.at.min.conc <- min.conc.mean + 3 * min.conc.sd

## Define the acceptable response range at the minimum concentration (max.response)

## Define the distribution of the response at the minimum concentration, shifted to 100
shifted.responses.at.min.conc <- responses.at.min.conc - mean(responses.at.min.conc) + 100

## Define the distribution of the response of the "LOW" controls
low.responses <- subset(ctrl.data, DRUG_ID == "LOW")$PERCENT_INHIBITION

lst <- list("Response at min conc (shifted)" = shifted.responses.at.min.conc, "LOW" = low.responses)

## Define the distribution of the response of the "BZCL" controls
if(any(ctrl.data$DRUG_ID == "BZCL")) {
  bzcl.responses <- subset(ctrl.data, DRUG_ID == "BZCL")$PERCENT_INHIBITION
##   bzcl.responses <- subset(orig.drug.data.long, DRUG_ID == "BZCL")$PERCENT_INHIBITION
  lst[["BzCl"]] <- bzcl.responses
}

indices <- 1:length(lst)
names(indices) <- names(lst)
file <- "max-conc-distributions-fimm-cm.pdf"
pdf(paste0(output.dir, "/", file))
df <- Reduce("rbind", llply(indices, .fun = function(i) data.frame(color = names(lst)[i], x = lst[[i]])))
g <- ggplot(data = df, aes(x = x, color = color))
g <- g + geom_density()
g <- g + xlab("PERCENT_INHIBITION")
g <- g + ggtitle("FIMM CM") + theme(legend.title=element_blank())
g
d <- dev.off()
l_ply(indices, 
      .fun = function(i) {
         qnt <- quantile(lst[[i]], probs = c(0.05, 0.95))
         cat(paste0(names(lst)[i], " mean = ", mean(lst[[i]]), " sd = ", sd(lst[[i]]), " ",
             paste(unlist(lapply(1:length(qnt), function(i) paste0(names(qnt)[i], ": ", unname(qnt[i])))), collapse = " "), "\n"))

      })

## Use the BzCl controls to define the distribution at the maximum concentration
cntrl.name <- "BzCl"
cntrl.name <- "LOW"
## Use the LOW controls to define the distribution at the maximum concentration
max.conc.mean <- mean(lst[[cntrl.name]])
max.conc.sd <- sd(lst[[cntrl.name]])
max.response.at.max.conc <- max.conc.mean + 3 * max.conc.sd
max.response <- max.response.at.max.conc

## Add a flag outlier.response (== TRUE) if the concentration/inhibition pair appears to be an outlier

if(TRUE) {
## Annotate outliers based on min/max responses
drug.data.long <- annotate.drug.response.outliers_(drug.data.long, drug.screen.cols, conc.col, response.col, 
                                                   min_response_min_conc = min.response.at.min.conc, 
                                                   max_response_min_conc = max.response.at.min.conc, 
                                                   max_response = max.response, prefix = paste0(output.dir, "/", prefix))
}

if(FALSE) {
## Annotate outliers based on pvals defined from distributions
drug.data.long <- annotate.drug.response.outlier.pvals_(drug.data.long, response.col, min.conc.mean, min.conc.sd, max.conc.mean, max.conc.sd)
}

## Plot responses vs each of the experimental factors.  Detect any outliers that we should exclude based on these experimental factors.
## In particular, for each level of each factor, calculate the fraction of samples that have response outside the range min.response to max.response.
## If that fraction is an outlying (and high) fraction across all such fractions for that factor, exclude the level. 
source("../common/drug-response-quality-control.R")
min_response <- min.response.at.min.conc
max_response <- max(max.response.at.max.conc, max.response.at.min.conc)
exclude <- correlate.drug.response.with.experimental.factors(drug.data.long, experimental.factors, response.col, min_response, max_response, prefix = paste0(output.dir, "/", prefix))

## Perform the fits, after excluding outliers
fits <- list()
fcts <- c("LL.4", "L.4")
outlier.col <- "outlier.response"
for(fct in fcts) {
  cat(paste0("Fitting drug response using ", fct, "\n"))
  flag <- !(drug.data.long[, drug.col] %in% drugs.to.exclude)
  fits[[fct]] <- fit.drc(drug.data.long[flag, c(all.cols, outlier.col)], fct = fct, drug.screen.cols = drug.screen.cols, conc.nM.col = conc.col, 
                         response.col = response.col, response.is.viability = response.is.viability, outlier.col = outlier.col, compute.auc = FALSE)

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
  exclusion.cols <- c(exclusion.col)
  factor.cols <- c(NA)
  fits[[fct]][, exclusion.col] <- FALSE
  fits[[fct]][exclusion.flag, exclusion.col] <- TRUE
 
  ## Exclude based on cutoffs established above
  for(nm in names(exclude)) {
    exclusion.col <- paste0(nm, ".exclude")
    exclusion.flag <- !is.na(fits[[fct]][, nm]) & (fits[[fct]][, nm] %in% exclude[[nm]])
    fits[[fct]][, exclusion.col] <- FALSE
    fits[[fct]][exclusion.flag, exclusion.col] <- TRUE
    exclusion.cols <- c(exclusion.cols, exclusion.col)
    factor.cols <- c(factor.cols, nm)
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
  ## Consider two consecutive components of a vector v[i] and v[i+1] to be approximately _non_-monotonic increasing
  ## if (v[i] - v[i+1]) > 30 -- i.e., allowing a tolerance of 30%
  any.exclude <- unlist(apply(fits[[fct]][, grepl(pattern="exclude", x=colnames(fits[[fct]]))], 1, function(row) any(row)))
  tmp <- fits[[fct]][!any.exclude,]
  approx.monotone <- unlist(lapply(tmp$all.responses, function(str) {
                                                          vec <- as.numeric(unlist(strsplit(str, split=",[ ]*")))
                                                          approximately.monotonic(vec, direction='inc', tolerance = 10) || 
                                                             approximately.monotonic(vec, direction='dec', tolerance = 10)
                                                        }))
  approx.non.monotone <- unlist(lapply(tmp$all.responses, function(str) {
                                                          vec <- as.numeric(unlist(strsplit(str, split=",[ ]*")))
                                                          !approximately.monotonic(vec, direction='inc', tolerance = 30) &&
                                                             !approximately.monotonic(vec, direction='dec', tolerance = 30)
                                                        }))
  monotone <- rep(NA, nrow(tmp))
  flag <- approx.monotone == TRUE
  monotone[flag] <- TRUE
  flag <- approx.non.monotone == TRUE
  monotone[flag] <- FALSE
  flag <- !is.na(monotone)

  file <- paste0(prefix, "-", fct, "-gof-cutoff.pdf")
  pdf(paste0(output.dir, "/", file), onefile=FALSE)
  tbl <- data.frame(gof = as.numeric(tmp$gof[flag]), condition = monotone[flag])
##  gof.cutoff <- round(median(tbl$gof[tbl$condition == FALSE]), digits=2)
  suppressPackageStartupMessages(library(maxstat))
  formula <- as.formula("condition ~ gof")
  mt <- maxstat.test(formula, data=tbl, smethod="Wilcoxon", pmethod="none")
  gof.cutoff <- mt$estimate
  cat(paste0("Using GOF cutoff: ", gof.cutoff, "\n"))
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
  exclusion.cols <- c(exclusion.cols, exclusion.col)
  factor.cols <- c(factor.cols, NA)
  fits[[fct]][, exclusion.col] <- FALSE
  fits[[fct]][exclusion.flag, exclusion.col] <- TRUE

  any.exclude <- unlist(apply(fits[[fct]][, grepl(pattern="exclude", x=colnames(fits[[fct]]))], 1, function(row) any(row)))
  exclusion.flag <- any.exclude
  exclusion.col <- "any.exclude"
  exclusion.cols <- c(exclusion.cols, exclusion.col)
  factor.cols <- c(factor.cols, NA)
  fits[[fct]][, exclusion.col] <- FALSE
  fits[[fct]][exclusion.flag, exclusion.col] <- TRUE

##  for(exclusion.col in colnames(fits[[fct]])[grepl(pattern="exclude", x=colnames(fits[[fct]]))]) {
  for(indx in 1:length(exclusion.cols)) {
      exclusion.col <- exclusion.cols[indx]
      factor.col <- factor.cols[indx]
      cat(paste0("\n\nNum excluded: ", exclusion.col, "\n"))
      print(table(fits[[fct]][, exclusion.col]))
      exclude.flag <- fits[[fct]][, exclusion.col] == TRUE
      if(any(exclude.flag)) {
        if(!is.na(factor.col)) {
          excluded <- unique(fits[[fct]][exclude.flag, factor.col])
          cat(paste0("Num ", factor.col, " excluded: ", length(excluded), "\n"))
          cat(paste0(paste(excluded, collapse=", "), "\n"))
        }
      }
  }
}

cat("Saving image\n")
rdata.file <- paste0(".Rdata.", prefix)
save.image(rdata.file)

common.cols <- c(drug.screen.cols, "min.conc.nM", "max.conc.nM", "all.concs.nM", "uniq.concs.nM", "all.responses",
                 "min.conc.w.outliers.nM", "max.conc.w.outliers.nM", "all.concs.w.outliers.nM", "uniq.concs.w.outliers.nM", "all.responses.w.outliers")
fits <- merge(fits[["LL.4"]], fits[["L.4"]], by = common.cols, suffixes = c(".ll4", ".l4"), all = TRUE)

my.dup <- function(x) duplicated(x, fromLast = TRUE) | duplicated(x, fromLast = FALSE)

## Make sure each drug is tested over a common range across patients

cols <- c(drug.col, "min.conc.w.outliers.nM", "max.conc.w.outliers.nM")
df <- fits[!fits$any.exclude.ll4, cols]
uniq.conc.orig <- ddply(df, .variables = cols, .fun = function(df) df$cnt <- nrow(df))
colnames(uniq.conc.orig) <- c(cols, "cnt")

cat("Note that the original ranges (before dropping any points) are the same for all of the drugs (modulo rounding errors)\n")
print(uniq.conc.orig[my.dup(uniq.conc.orig[, drug.col]),])

cols <- c(drug.col, "min.conc.nM", "max.conc.nM")
df <- fits[!fits$any.exclude.ll4, cols]
uniq.conc <- ddply(df, .variables = cols, .fun = function(df) df$cnt <- nrow(df))
colnames(uniq.conc) <- c(cols, "cnt")
uniq.conc <- uniq.conc[order(uniq.conc[, drug.col]),]
print(uniq.conc[my.dup(uniq.conc[, drug.col]),])


source("../common/restrict-concentration-to-common-range.R")

dummy.drug.map <- data.frame("inhibitor" = unique(fits[, drug.col]))
colnames(dummy.drug.map)[colnames(dummy.drug.map) == "inhibitor"] <- drug.col
dummy.drug.map.drug.id.cols <- list("data.set" = drug.col)
names(dummy.drug.map.drug.id.cols) <- c(data.set)
fit.list <- list("data.set" = fits)
names(fit.list) <- data.set
ret <- determine.drug.ranges(fit.list, screen.id.cols, data.set.drug.id.cols, 
                             dummy.drug.map, dummy.drug.map.drug.id.cols, shared.conc.postfix = ".auc.conc") 

print(uniq.conc[my.dup(uniq.conc[, drug.col]),])
flag <- ret$shared.drug.ranges[["LL.4"]][, drug.col]  %in% uniq.conc[my.dup(uniq.conc[, drug.col]),drug.col]
m <- merge(uniq.conc[my.dup(uniq.conc[, drug.col]), ], ret$shared.drug.ranges[["LL.4"]], by = drug.col)

## Compute the auc
tmp <- ret$fits.with.concs[[data.set]]
for(nm in names(tmp)) {
  ## Drop the columns with .ll4 prefix in the L.4 fits and vice versa
  pattern = "\\.l4"
  if(nm == "L.4") { pattern = "\\.ll4" }
  tmp[[nm]] <- tmp[[nm]][, !grepl(x = colnames(tmp[[nm]]), pattern = pattern)]
}

fits <- merge(tmp[["LL.4"]], tmp[["L.4"]], by = common.cols, all = TRUE)

source("../common/dss.R")

fcts <- list("LL.4" = "LL.4", "L.4" = "L.4")
postfixes <- list("LL.4" = ".ll4", "L.4" = ".l4")

for(nm in names(fcts)) {  
  min.conc.col <- paste0("min.auc.conc", postfixes[[nm]])
  max.conc.col <- paste0("max.auc.conc", postfixes[[nm]])
  b.param.col <- paste0("b", postfixes[[nm]])
  c.param.col <- paste0("c", postfixes[[nm]])
  d.param.col <- paste0("d", postfixes[[nm]])
  e.param.col <- paste0("e", postfixes[[nm]])
  converged.col <- paste0("converged", postfixes[[nm]])
  auc.output.col <- paste0("auc", postfixes[[nm]])
  shifted.auc.output.col <- paste0("shifted.auc", postfixes[[nm]])
  fits <- compute.all.auc(fits, fct = fcts[[nm]], min.conc.col = min.conc.col, max.conc.col = max.conc.col, 
                          b.param.col = b.param.col, c.param.col = c.param.col, d.param.col = d.param.col, e.param.col = e.param.col, 
                          converged.col = converged.col, auc.output.col = auc.output.col, shifted.auc.output.col = shifted.auc.output.col)
}


cat("Saving image\n")
save.image(rdata.file)

write.table(file=output.file, fits, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
f <- File(output.file, parentId = parentId, synapseStore = TRUE)
synStore(f, executed = NULL, forceVersion = FALSE)

cat("Done fitting drug response data using L.4 and LL.4\n")

q()

