suppressPackageStartupMessages(library("nplr"))
suppressPackageStartupMessages(library("foreach"))
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(tidyr))

suppressPackageStartupMessages(library("drc"))
suppressPackageStartupMessages(library("minpack.lm"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("gridExtra"))
suppressPackageStartupMessages(library("caTools"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("maftools"))

library(ggbeeswarm)
library(plyr)
library(stringr)
library(synapseClient)

## Add ensembl identifiers to the drug annotation table
## Translate gene symbols to ensg identifiers
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(Biobase))

suppressPackageStartupMessages(library(gplots))

load("../050517/.RData.fits")
cat("Loaded fits\n")

synapseLogin()

cat("This script in 050817/run-models.R performs 10-fold CV of z-scored IC50 vs z-scored expression\n")

## Inputs
zscore.expr.data <- TRUE
zscore.drug.data <- TRUE
gof.cutoff <- 0.4

output.path <- "output"
if (!file.exists(output.path)) {
  dir.create(output.path)
}

my.dup <- function(x) {
  duplicated(x, fromLast=TRUE) | duplicated(x, fromLast=FALSE)
}

suppressPackageStartupMessages(library("scales")) ## for scientific

Sys.setenv(R_ZIPCMD= "/usr/bin/zip")  

num.cores <- detectCores()
if(!is.na(num.cores) && (num.cores > 1)) {
  suppressPackageStartupMessages(library("doMC"))
  cat(paste("Registering ", num.cores-1, " cores.\n", sep=""))
  registerDoMC(cores=(num.cores-1))
}

## mutation.types.to.total: list of mutation types that will be considered a mutation (i.e., a "1")
##                          in the final 0/1 matrix.
## e.g., c("Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Splice_Site", "Translation_Start_Site")
read.mafs <- function(maf.names, maf.tbl, mutation.types.to.total = c("Missense_Mutation", "Nonsense_Mutation")) {
  ## Read in the MAFs
  maf.list  <- list()
  for(i in 1:nrow(maf.tbl)){ maf.list[[maf.names[i]]] <-  read.maf(getFileLocation(synGet(maf.tbl$file.id[i])), verbose = F); print(i)}
  
  maf.tbl <- ldply(maf.list, getGeneSummary)
  rm(maf.list); gc()
  
  ## Tally mutation.types.to.total (e.g., Missense_Mutation and Nonsense_Mutation)
  maf.tbl$total <- unlist(apply(maf.tbl[,mutation.types.to.total], 1, function(row) sum(row, na.rm=TRUE)))
  print(colnames(maf.tbl))
  maf.tbl <- maf.tbl[,c(".id","Hugo_Symbol","total")]
  
  maf.tbl                  <- spread(maf.tbl, Hugo_Symbol,total); 
  rownames(maf.tbl)        <- maf.tbl$.id; 
  maf.tbl                  <- maf.tbl[,-1]; 
  maf.tbl[is.na(maf.tbl)] <- 0
  maf.tbl[maf.tbl > 0] <- 1
  maf.tbl
}

## Filter fits
## min.gof:  exclude fits having a gof < min.gof
## remove.non.finite.ic50:  remove ic50 values that are NA or NaN.
## require.ic50.be.in.range:  require that ic50 be between min and max concentration tested
## max.replicate.cv:  require that the coefficient of variation (sd/mean) of the IC50 across replicates be < max.replicate.cv
## max.specimen.cv:  require that the coefficient of variation (sd/mean) of the IC50 across specimens (i.e., labids) be < max.replicate.cv
## specimen.col: the column holding the replicate identifier
## patient.col: the column holding the patient identifier
## inhibitor.col: the column holding the inhibitor identifier
filter.drug.response.fits <- function(fits.df, min.gof = 0, remove.non.finite.ic50 = TRUE, require.ic50.be.in.range = TRUE, max.replicate.cv = NULL, max.specimen.cv = NULL, specimen.col = "lab_id", patient.col = "patient_id", inhibitor.col = "inhibitor") {
  if(remove.non.finite.ic50) {
    old.nrow <- nrow(fits.df)
    flag <- is.finite(fits.df$IC50)
    fits.df <- fits.df[flag,]
    new.nrow <- nrow(fits.df)
    cat(paste0("Filtering by finite IC50 reduces ", old.nrow, " fits to ", new.nrow, "\n"))
  }

  if(require.ic50.be.in.range) {
    old.nrow <- nrow(fits.df)
    flag <- is.finite(fits.df$IC50) & is.finite(fits.df$Max.Conc.tested) & is.finite(fits.df$Min.Conc.tested) & ((10^fits.df$IC50) >= fits.df$Min.Conc.tested) & ((10^fits.df$IC50) <= fits.df$Max.Conc.tested)
    fits.df <- fits.df[flag,]
    new.nrow <- nrow(fits.df)
    cat(paste0("Filtering by requiring IC50 be in range reduces ", old.nrow, " fits to ", new.nrow, "\n"))
  }
  
  if(min.gof > 0) {
    old.nrow <- nrow(fits.df)
    flag <- is.finite(fits.df$gof) & (fits.df$gof >= min.gof)
    fits.df <- fits.df[flag,]
    new.nrow <- nrow(fits.df)
    cat(paste0("Filtering by GOF < ", min.gof, " reduces ", old.nrow, " fits to ", new.nrow, "\n"))
  }
  
  if(!is.null(max.replicate.cv)) {
    ## Calculate CV across replicates within a specimen
    old.nrow <- nrow(fits.df)
    cv <- ddply(fits.df[,c("IC50", specimen.col, patient.col, inhibitor.col)], c(specimen.col, patient.col, inhibitor.col),
                .parallel = TRUE,
                .fun = function(df) {
                  cv <- max.replicate.cv + 1
                  ic50s <- df$IC50
                  num.reps <- length(which(is.finite(ic50s)))
                  m <- NA
                  if(num.reps >= 1) {
                    ic50s <- 10^(ic50s[is.finite(ic50s)])
                    cv <- 0
                    m <- mean(ic50s)
                    if((num.reps > 1) && (m != 0)) {
                      cv <- sd(ic50s) / m
                    }
                  }
                  c(m, cv, num.reps)
                })
    colnames(cv) <- c(specimen.col, patient.col, inhibitor.col, "ic50.repl.mean", "ic50.repl.cv", "num.repl")
    flag <- (cv$ic50.repl.cv < max.replicate.cv) | (cv$num.repl == 1)
    cv <- cv[flag,]
    fits.df <- merge(fits.df, cv)
    new.nrow <- nrow(fits.df)
    cat(paste0("Filtering by max replicate-level CV = ", max.replicate.cv, " reduces ", old.nrow, " fits to ", new.nrow, "\n"))
  }
  
  if(!is.null(max.specimen.cv)) {
    ## Calculate CV across specimens within a patient -- i.e., mean of replicate means and sd of replicate means
    old.nrow <- nrow(fits.df)
    cv <- ddply(fits.df[,c("ic50.repl.mean", patient.col, inhibitor.col)], c(patient.col, inhibitor.col),
                .parallel = TRUE,
                .fun = function(df) {
                  cv <- max.specimen.cv + 1
                  ic50s <- df$ic50.repl.mean
                  num.reps <- length(which(is.finite(ic50s)))
                  m <- NA
                  if(num.reps >= 1) {
                    ## NB: these mean ic50s have already been exponentiated above
                    ic50s <- (ic50s[is.finite(ic50s)])
                    cv <- 0
                    m <- mean(ic50s)
                    if((num.reps > 1) && (m != 0)) {
                      cv <- sd(ic50s) / m
                    }
                  }
                  c(m, cv, num.reps)
                })
    colnames(cv) <- c(patient.col, inhibitor.col, "ic50.specimen.mean", "specimen.cv", "num.specimens")
    flag <- (cv$specimen.cv < max.specimen.cv) | (cv$num.specimens == 1)
    cv <- cv[flag,]
    fits.df <- merge(fits.df, cv)
    new.nrow <- nrow(fits.df)
    cat(paste0("Filtering by max specimen-level CV = ", max.specimen.cv, " reduces ", old.nrow, " fits to ", new.nrow, "\n"))
  }

  fits.df
}



## WORKING

## We will combine columns, i.e., replicates should have the same length as the number of columns.
## Columns i and j will be combined if replicates[i] == replicates[j]
combine_columns_across_replicates <- function(matrix, replicates, agg.fun = "mean", debug = FALSE) {
  suppressPackageStartupMessages(library("dplyr"))
  
  supported.agg.funcs <- c("mean", "max", "min")
  if(!(agg.fun %in% supported.agg.funcs)) {
    stop(paste0("Supported aggregate functions are: ", paste0(supported.agg.funcs, collapse=","), "\n"))
  }
  
  if(ncol(matrix) != length(replicates)) {
    stop(paste0("Number of columns (", ncol(matrix), ") and length of replicates (", length(replicates), ") should be the same\n"))
  }
  
  df <- as.data.frame(t(matrix))
  rownames(df) <- 1:nrow(df)
  colnames(df) <- 1:ncol(df)
  names(replicates) <- 1:length(replicates)
  df$group <- replicates
  agg <- t((df %>% group_by(group) %>% summarise_each(funs_(agg.fun))))
  if(debug == TRUE) {
    ## Output a column in the original matrix that is duplicated in the output matrix
    dups <- replicates[duplicated(replicates, fromLast=TRUE) | duplicated(replciates, fromLast=FALSE)]
    if(length(dups) > 1) {
      ## Take the first duplicate case
      dups <- dups[dups == dups[1]]
      cat("These columns:\n")
      print(head((t(df))[, names(dups)]))
      cat(paste0("Were combined using ", agg.fun, " into:\n"))
      print(head(agg[, (agg[1,] %in% dups)]))
    }
  }
  colnames(agg) <- agg[1,]
  agg <- (agg[-1,])
  class(agg) <- "numeric"
  df <- as.data.frame(agg)
  df
}

## specimen ids are columns of expr.mat
## specimen.id.col is col of drug data frame holding the specimen id
## patient.id.col is col of drug data frame holding the patient id
## Perform 10-fold CV.
correlate.expression.with.ic50 <- function(drug.df, drug.id.col, drugs.to.test, zscore.response = FALSE, zscore.expression = FALSE, 
                                           genes.to.test = NULL, expr.mat, specimen.id.col, patient.id.col) {
  df <- drug.df[drug.df[,drug.id.col] %in% drugs.to.test,]
  if(is.null(genes.to.test)) {
    df <- df[!is.na(df$Ensg.Targets),]
  }
  if(is.null(genes.to.test)) {
    flag <- unlist(lapply(df$Ensg.Targets, function(x) any(unlist(strsplit(x, split=",[ ]*")) %in% rownames(expr.mat))))
    df <- df[flag,]
  }
  if(!is.null(genes.to.test)) {
    if(!(any(genes.to.test %in% rownames(expr.mat)))) {
      warning("UNEXPECTED no genes to test are in data\n")
    }
  }
  ret <- ddply(df, drug.id.col, .parallel = TRUE,
               .fun = function(df.drug) {
                 cat(paste0("Correlating gene expr with IC50 for drug ", unique(df.drug[,drug.id.col]), "\n"))
                 drug.targets <- df.drug$Ensg.Targets[1]
                 drug.targets <- unlist(strsplit(drug.targets, split=",[ ]*"))
                 if(!is.null(genes.to.test)) {
                   drug.targets <- genes.to.test
                 }
                 drug.targets <- drug.targets[drug.targets %in% rownames(expr.mat)]
                 resp <- data.frame(IC50 = df.drug$IC50, PATIENT_ID = df.drug[,patient.id.col], id = df.drug[,specimen.id.col])
                 resp <- resp[!is.infinite(resp$IC50) & !is.nan(resp$IC50) & !is.na(resp$IC50),]
                 common.samples <- intersect(resp$id, colnames(expr.mat))
                 num.resp.samples <- length(resp$id)
                 num.unique.resp.samples <- length(unique(resp$id))
                 flag <- resp$id %in% common.samples
                 resp <- resp[flag,]
                 tmp <- as.data.frame(table(resp$PATIENT_ID))
                 colnames(tmp) <- c("PATIENT_ID", "weight")
                 tmp$weight <- 1/sqrt(tmp$weight)
                 resp <- merge(resp, tmp)
                 common.samples <- as.character(resp$id)
                 ldply(drug.targets, .parallel = FALSE,
                       .fun = function(gene) {
                         if(!(gene %in% rownames(expr.mat))) {
                           warning(paste0("UNEXPECTED: ", gene, " not in expr.mat\n"))
                         }
                         if(!(all(common.samples %in% colnames(expr.mat)))) {
                           warning("UNEXPECTED columns not in expr.mat\n")
                         }
                         gene.expr <- as.vector(unlist(expr.mat[gene, common.samples]))
                         if(zscore.response) {
                           ic50.sd <- sd(resp$IC50[is.finite(resp$IC50)], na.rm=TRUE)
                           if((ic50.sd == 0) || !is.finite(ic50.sd)) { ic50.sd <- 1 }
                           resp$IC50 <- ( resp$IC50 - mean(resp$IC50[is.finite(resp$IC50)], na.rm=TRUE) ) / ic50.sd
                         }
                         if(zscore.expression) {
                           gene.sd <- sd(gene.expr[is.finite(gene.expr)], na.rm=TRUE)
                           if((gene.sd == 0) || !is.finite(gene.sd)) { gene.sd <- 1 }
                           gene.expr <- ( gene.expr - mean(gene.expr[is.finite(gene.expr)], na.rm=TRUE) ) / gene.sd
                         }
                         ## Perform 10-fold CV
                         n_folds <- 10
                         folds_i <- sample(rep(1:n_folds, length.out = length(gene.expr)))
                         res.names <- c("gene", "pval", "r2", "k", "predicted", "actual")
                         res <- matrix(NA, nrow = n_folds, ncol = length(res.names))
                         colnames(res) <- res.names
                         for (k in 1:n_folds) {
                           test_i <- which(folds_i == k)
                           df <- data.frame(x = gene.expr[-test_i], y = resp$IC50[-test_i])
                           lm.fit <- lm(data = df, y ~ x)
                           sum <- summary(lm.fit)
                           f <- sum$fstatistic
                           r2 <- sum$r.squared
                           p <- pf(f[1],f[2],f[3],lower.tail=F)
                           num.unique.samples <- length(unique(common.samples))
                           num.samples <- length(common.samples)
                           actual <- resp$IC50[test_i]
                           pred <- unname(predict(lm.fit, data.frame(x=gene.expr[test_i])))
                           res[k, ] <- c(gene, p, r2, k, paste(pred, collapse=","), paste(actual, collapse=","))
                         }
                         res
                       })
               })
  ret
}

common.genes <- intersect(rownames(ohsu.expr.data), rownames(fimm.expr.data))

ohsu.expr.data <- ohsu.expr.data[common.genes, ]
fimm.expr.data <- fimm.expr.data[common.genes, ]
gc()

cat("Filtering OHSU drug response\n")
ohsu.fits.gof.cutoff.df <- filter.drug.response.fits(ohsu.fits.df, min.gof = gof.cutoff, remove.non.finite.ic50 = TRUE, require.ic50.be.in.range = TRUE, max.replicate.cv = 0.5, max.specimen.cv = 0.5)
cat("Filtering FIMM drug response\n")
fimm.fits.gof.cutoff.df <- filter.drug.response.fits(fimm.fits.df, min.gof = gof.cutoff, remove.non.finite.ic50 = TRUE, require.ic50.be.in.range = TRUE, max.replicate.cv = NULL, max.specimen.cv = NULL)

rm(ohsu.fits.df)
rm(fimm.fits.df)
gc()
ohsu.fits.df <- ohsu.fits.gof.cutoff.df
fimm.fits.df <- fimm.fits.gof.cutoff.df
rm(ohsu.fits.gof.cutoff.df)
rm(fimm.fits.gof.cutoff.df)
gc()


if(zscore.drug.data) {
  cat("Zscoring OHSU drug data\n")
  drug.id.col <- "ID_Drug"
  ohsu.fits.df <- ddply(ohsu.fits.df, drug.id.col, 
                        .fun = function(df) {
                          ic50.sd <- sd(df$IC50, na.rm = TRUE)
                          if( (ic50.sd == 0) || !is.finite(ic50.sd) ) { ic50.sd <- 1 }
                          ret <- ( df$IC50 - mean(df$IC50, na.rm = TRUE) ) / ic50.sd
                          df$IC50.z <- ret
                          df
                        })
  ohsu.fits.df$IC50 <- ohsu.fits.df$IC50.z

  cat("Zscoring FIMM drug data\n")
  drug.id.col <- "DRUG_ID"
  fimm.fits.df <- ddply(fimm.fits.df, drug.id.col, 
                        .fun = function(df) {
                          ic50.sd <- sd(df$IC50, na.rm = TRUE)
                          if( (ic50.sd == 0) || !is.finite(ic50.sd) ) { ic50.sd <- 1 }
                          ret <- ( df$IC50 - mean(df$IC50, na.rm = TRUE) ) / ic50.sd
                          df$IC50.z <- ret
                          df
                        })
  fimm.fits.df$IC50 <- fimm.fits.df$IC50.z
}

min.avg.log.cpm <- 1
cat(paste0("Ensuring avg log CPM > ", min.avg.log.cpm, "\n"))
flag <- unlist(apply(ohsu.expr.data, 1, function(row) mean(row, na.rm=TRUE) > min.avg.log.cpm))
ohsu.expr.data <- ohsu.expr.data[flag, ]

flag <- unlist(apply(fimm.expr.data, 1, function(row) mean(row, na.rm=TRUE) > min.avg.log.cpm))
fimm.expr.data <- fimm.expr.data[flag, ]

if(zscore.expr.data) {
  cat("Zscoring OHSU expression data\n")
  ohsu.expr.data <- t(scale(t(ohsu.expr.data), center = TRUE, scale = TRUE))
  
  cat("Zscoring FIMM expression data\n")
  fimm.expr.data <- t(scale(t(fimm.expr.data), center = TRUE, scale = TRUE))
}

ohsu.rna.replicates <- rnaseq.sample.summary.tbl$PatientID
names(ohsu.rna.replicates) <- rnaseq.sample.summary.tbl$SeqID
ohsu.rna.replicates <- ohsu.rna.replicates[!duplicated(names(ohsu.rna.replicates))]
ohsu.rna.replicates <- ohsu.rna.replicates[names(ohsu.rna.replicates) %in% colnames(ohsu.expr.data)]

cat(paste0(length(which(colnames(ohsu.expr.data) %in% ohsu.fits.df$SeqID)), " of ", ncol(ohsu.expr.data), " OHSU specimens with expression have drug response data.\n"))
flag <- ohsu.fits.df$SeqID %in% colnames(ohsu.expr.data)
cat(paste0(length(unique(ohsu.fits.df$patient_id[flag])), " unique OSHU patients have at least one specimen with both RNA-seq and drug response data.\n"))
cat(paste0(length(which(unique(ohsu.rna.replicates) %in% ohsu.fits.df$patient_id)), " unique OSHU patients have at least one specimen with RNA-seq and at least one specimen with drug response data.\n"))

ohsu.expr.data <- ohsu.expr.data[, colnames(ohsu.expr.data) %in% names(ohsu.rna.replicates)]
ohsu.rna.replicates <- ohsu.rna.replicates[colnames(ohsu.expr.data)]

## ohsu.expr.data <- combine_columns_across_replicates(ohsu.expr.data, ohsu.rna.replicates, agg.fun = "mean", debug = TRUE)
ohsu.expr.data <- ohsu.expr.data[, colnames(ohsu.expr.data) %in% ohsu.fits.df$SeqID]
ohsu.rna.replicates <- ohsu.rna.replicates[colnames(ohsu.expr.data)]
ohsu.expr.data <- ohsu.expr.data[, names(ohsu.rna.replicates)[!duplicated(ohsu.rna.replicates)]]

fimm.rna.replicates <- gsub(x=fimm.metadata$person, pattern="\\.", replacement="_")
print(fimm.rna.replicates)
names(fimm.rna.replicates) <- rownames(fimm.metadata)
fimm.rna.replicates <- fimm.rna.replicates[!duplicated(names(fimm.rna.replicates))]
fimm.rna.replicates <- fimm.rna.replicates[names(fimm.rna.replicates) %in% colnames(fimm.expr.data)]

cat(paste0(length(which(colnames(fimm.expr.data) %in% fimm.fits.df$SCREEN_ID)), " of ", ncol(fimm.expr.data), " FIMM specimens with expression have drug response data.\n"))
flag <- fimm.fits.df$SCREEN_ID %in% colnames(fimm.expr.data)
u1 <- unique(fimm.fits.df$PATIENT_ID[flag])
cat(paste0(length(unique(fimm.fits.df$PATIENT_ID[flag])), " unique FIMM patients have at least one specimen with both RNA-seq and drug response data.\n"))
print(unique(fimm.fits.df$PATIENT_ID[flag]))
u <- unique(fimm.rna.replicates)
print(u1[!(u1 %in% u)])
## "FHRB_370" "FHRB_250"
cat(paste0(length(which(u %in% fimm.fits.df$PATIENT_ID)), " unique FIMM patients have at least one specimen with RNA-seq and at least one specimen with drug response data.\n"))
print(u[u %in% fimm.fits.df$PATIENT_ID])

fimm.expr.data <- fimm.expr.data[, colnames(fimm.expr.data) %in% names(fimm.rna.replicates)]
fimm.rna.replicates <- fimm.rna.replicates[colnames(fimm.expr.data)]

## fimm.expr.data <- combine_columns_across_replicates(fimm.expr.data, fimm.rna.replicates, agg.fun = "mean")
fimm.expr.data <- fimm.expr.data[, colnames(fimm.expr.data) %in% fimm.fits.df$SCREEN_ID]
fimm.rna.replicates <- fimm.rna.replicates[colnames(fimm.expr.data)]
fimm.expr.data <- fimm.expr.data[, names(fimm.rna.replicates)[!duplicated(fimm.rna.replicates)]]



common.drugs.tested <- intersect(unique(ohsu.fits.df$ID_Drug), unique(fimm.fits.df$DRUG_ID))
common.drugs.tested <- intersect(common.drugs.tested, ohsu.fimm.drugs$ID_Drug[!is.na(ohsu.fimm.drugs$Gene.Targets)])
common.drugs.tested <- intersect(common.drugs.tested, ohsu.fimm.drugs$ID_Drug[!is.na(ohsu.fimm.drugs$Ensg.Targets)])

save.image(".RData")

ohsu.expr.corr <- correlate.expression.with.ic50(drug.df = ohsu.fits.df, drug.id.col = "ID_Drug", genes.to.test = common.genes, drugs.to.test = common.drugs.tested, expr.mat = ohsu.expr.data, specimen.id.col = "SeqID", patient.id.col = "patient_id", zscore.response = FALSE, zscore.expression = FALSE)
cat("Done correlating expression with IC50 for OHSU\n")
save.image(".RData")

## common.drugs.tested <- intersect(common.drugs.tested, c("FIMM136387", "FIMM133832", "FIMM133867", "FIMM133902"))

fimm.expr.corr <- correlate.expression.with.ic50(drug.df = fimm.fits.df, drug.id.col = "DRUG_ID", genes.to.test = common.genes, drugs.to.test = common.drugs.tested, expr.mat = fimm.expr.data, specimen.id.col = "SCREEN_ID", patient.id.col = "PATIENT_ID", zscore.response = FALSE, zscore.expression = FALSE)
cat("Done correlating expression with IC50 for FIMM\n")
cat("Done with expression\n")
save.image(".RData")
cat("Saved image\n")
stop("stop")

ohsu.fimm.expr.corr <- merge(ohsu.expr.corr, fimm.expr.corr, by.x = c("ID_Drug", "gene"), by.y = c("DRUG_ID", "gene"), suffixes = c(".ohsu", ".fimm"))

subset(ohsu.fimm.expr.corr, (pval.ohsu < 0.05) & (pval.fimm < 0.05))
## ID_Drug            gene          pval.ohsu             r2.ohsu          pval.fimm           r2.fimm
## 100 FIMM003783 ENSG00000094631 0.0419410434254526  0.0494985601145971 0.0330721048250537 0.228393094777627
## 127 FIMM003794 ENSG00000122025 0.0140628306932017 0.00833458654791967 0.0449806673025468  0.24180451968676
ohsu.fimm.expr.corr$pval.fimm <- as.numeric(ohsu.fimm.expr.corr$pval.fimm)
ohsu.fimm.expr.corr$pval.ohsu <- as.numeric(ohsu.fimm.expr.corr$pval.ohsu)

## Plot pvalue distributions of FIMM and OHSU
pdf("ohsu-fimm-pval-dist.pdf")
g1 <- ggplot(ohsu.fimm.expr.corr, aes(x=pval.fimm)) + geom_histogram() + xlab("FIMM p-value")
g2 <- ggplot(ohsu.fimm.expr.corr, aes(x=pval.ohsu)) + geom_histogram() + xlab("OHSU p-value")
grid.arrange(g1, g2)
d <- dev.off()

ohsu.fimm.expr.corr$drug.gene <- unlist(apply(ohsu.fimm.expr.corr[,c("ID_Drug", "gene")], 1, function(row) paste0(row[1], "-", row[2])))
## Make Venn diagram of FIMM vs OHSU significant
vennList = list("FIMM"=ohsu.fimm.expr.corr$drug.gene[ohsu.fimm.expr.corr$pval.fimm < 0.05], "OHSU"=ohsu.fimm.expr.corr$drug.gene[ohsu.fimm.expr.corr$pval.ohsu < 0.05])
pdf("ohsu-fimm-sig-venn.pdf")
venn(vennList)
d <- dev.off()

old.nrow <- nrow(ohsu.fimm.expr.corr)
ohsu.fimm.expr.corr <- merge(ohsu.fimm.expr.corr, unique(fimm.ohsu.drug.annotations[,c("ID_Drug", "DRUG_NAME")]))
if(old.nrow != nrow(ohsu.fimm.expr.corr)) {
  warning("UNEXPECTED drug correlation results grew\n")
}

old.nrow <- nrow(ohsu.fimm.expr.corr)
ohsu.fimm.expr.corr <- merge(ohsu.fimm.expr.corr, unique(sym.to.ensg), by.x = "gene", by.y = "ENSG")
if(old.nrow != nrow(ohsu.fimm.expr.corr)) {
  warning("UNEXPECTED drug correlation results grew\n")
}

stop("stop")

common.genes <- intersect(rownames(beat.mafs), rownames(fimm.genomic.bin))
ohsu.mut.corr <- correlate.genomic.with.ic50(drug.df = ohsu.fits.gof.0.7.df, drug.id.col = "ID_Drug", drugs.to.test = common.drugs.tested, genes.to.test = common.genes, mutation.mat = beat.mafs, specimen.id.col = "SeqID", patient.id.col = "patient_id")
fimm.mut.corr <- correlate.genomic.with.ic50(drug.df = fimm.fits.gof.0.7.df, drug.id.col = "DRUG_ID", drugs.to.test = common.drugs.tested, genes.to.test = common.genes, mutation.mat = fimm.genomic.bin, specimen.id.col = "SCREEN_ID", patient.id.col = "PATIENT_ID")

cat("Done with genomic\n")
save.image(".RData")
cat("Saved image\n")
stop("stop")

drug = "FIMM003774"
gene = "FLT3"
pdf(paste0(drug, "-", gene, "-ohsu-fimm-mut.pdf"))
g1 <- plot.genomic.vs.ic50(drug.df = ohsu.fits.gof.0.7.df, drug = drug, gene = gene, drug.id.col = "ID_Drug", mutation.mat = beat.mafs, specimen.id.col = "SeqID", patient.id.col = "patient_id")
g1 <- g1 + ggtitle("OHSU")
g2 <- plot.genomic.vs.ic50(drug.df = fimm.fits.gof.0.7.df, drug = drug, gene = gene, drug.id.col = "DRUG_ID", mutation.mat = fimm.genomic.bin, specimen.id.col = "SCREEN_ID", patient.id.col = "PATIENT_ID")
g2 <- g2 + ggtitle("FIMM")
grid.arrange(g1, g2)
d <- dev.off()

##for(drug in unique(fimm.mut.corr$DRUG_ID)) {
##  print(drug)
##  ohsu.mut.corr <- correlate.genomic.with.ic50(drug.df = ohsu.fits.gof.0.7.df, drug.id.col = "ID_Drug", drugs.to.test = drug, genes.to.test = common.genes, mutation.mat = beat.mafs, specimen.id.col = "SeqID", patient.id.col = "patient_id")
##}

ohsu.fimm.mut.corr <- merge(ohsu.mut.corr, fimm.mut.corr, by.x = c("ID_Drug", "gene"), by.y = c("DRUG_ID", "gene"), suffixes = c(".ohsu", ".fimm"))

ohsu.fimm.mut.corr$drug.gene <- unlist(apply(ohsu.fimm.mut.corr[,c("ID_Drug", "gene")], 1, function(row) paste0(row[1], "-", row[2])))
vennList = list("FIMM"=ohsu.fimm.mut.corr$drug.gene[!is.na(ohsu.fimm.mut.corr$pval.fimm) & (ohsu.fimm.mut.corr$pval.fimm < 0.05)], "OHSU"=ohsu.fimm.mut.corr$drug.gene[!is.na(ohsu.fimm.mut.corr$pval.ohsu) & (ohsu.fimm.mut.corr$pval.ohsu < 0.05)])
pdf("ohsu-fimm-sig-mut-venn.pdf")
venn(vennList)
d <- dev.off()

ohsu.fimm.mut.corr$pval.fimm <- as.numeric(ohsu.fimm.mut.corr$pval.fimm)
ohsu.fimm.mut.corr$pval.ohsu <- as.numeric(ohsu.fimm.mut.corr$pval.ohsu)

## Plot pvalue distributions of FIMM and OHSU
pdf("ohsu-fimm-mut-pval-dist.pdf")
g1 <- ggplot(ohsu.fimm.mut.corr, aes(x=pval.fimm)) + geom_histogram() + xlab("FIMM p-value")
g2 <- ggplot(ohsu.fimm.mut.corr, aes(x=pval.ohsu)) + geom_histogram() + xlab("OHSU p-value")
grid.arrange(g1, g2)
d <- dev.off()


## specimen ids are columns of mutation.mat
## gene symbols are rows of mutation.mat
correlate.genomic.with.ic50 <- function(drug.df, drug.id.col, drugs.to.test, genes.to.test = NULL, mutation.mat, specimen.id.col, patient.id.col) {
  df <- drug.df[drug.df[,drug.id.col] %in% drugs.to.test,]
  if(is.null(genes.to.test)) {
    df <- df[!is.na(df$Gene.Targets),]
  }
  if(is.null(genes.to.test)) {
    flag <- unlist(lapply(df$Gene.Targets, function(x) any(unlist(strsplit(x, split=",[ ]*")) %in% rownames(mutation.mat))))
    df <- df[flag,]
  }
  if(!is.null(genes.to.test)) {
    if(!(any(genes.to.test %in% rownames(mutation.mat)))) {
      warning("UNEXPECTED no genes to test are in data\n")
    }
  }
  ret <- ddply(df, drug.id.col, .parallel = TRUE,
               .fun = function(df.drug) {
                 drug.targets <- df.drug$Gene.Targets[1]
                 drug.targets <- unlist(strsplit(drug.targets, split=",[ ]*"))
                 if(!is.null(genes.to.test)) {
                   drug.targets <- genes.to.test
                 }
                 drug.targets <- drug.targets[drug.targets %in% rownames(mutation.mat)]
                 resp <- data.frame(IC50 = df.drug$IC50, PATIENT_ID = df.drug[,patient.id.col], id = df.drug[,specimen.id.col])
                 resp <- resp[!is.infinite(resp$IC50) & !is.nan(resp$IC50) & !is.na(resp$IC50),]
                 common.samples <- intersect(resp$id, colnames(mutation.mat))
                 num.resp.samples <- length(resp$id)
                 num.unique.resp.samples <- length(unique(resp$id))
                 if(length(common.samples) == 0) {
                   return(ldply(drug.targets, .parallel = FALSE,
                                .fun = function(gene) {
                                  p <- NA
                                  num.unique.samples <- length(unique(common.samples))
                                  num.samples <- length(common.samples)
                                  res <- c(gene, p, num.unique.samples, num.samples, num.unique.resp.samples, num.resp.samples)
                                  names(res) <- c("gene", "pval", "num.uniq.samples", "num.samples", "num.uniq.resp.samples", "num.resp.samples")
                                  res
                                }))
                 }
                 flag <- resp$id %in% common.samples
                 resp <- resp[flag,]
                 tmp <- as.data.frame(table(resp$PATIENT_ID))
                 colnames(tmp) <- c("PATIENT_ID", "weight")
                 tmp$weight <- 1/sqrt(tmp$weight)
                 resp <- merge(resp, tmp)
                 common.samples <- as.character(resp$id)
                 ldply(drug.targets, .parallel = FALSE,
                       .fun = function(gene) {
                         if(!(gene %in% rownames(mutation.mat))) {
                           warning(paste0("UNEXPECTED: ", gene, " not in mutation.mat\n"))
                         }
                         if(!(all(common.samples %in% colnames(mutation.mat)))) {
                           warning("UNEXPECTED columns not in mutation.mat\n")
                         }
                         gene.mutation <- as.vector(unlist(mutation.mat[gene, common.samples]))
                         if(!(all(gene.mutation %in% c(0,1)))) {
                           warning(paste0("WARNING: non-binary gene mutations: ", paste(gene.mutation, collapse=", "), "\n"))
                         }
                         p <- NA
                         if(all((c(0,1) %in% gene.mutation))) {
                           test.df <- data.frame(x = factor(gene.mutation), y = resp$IC50)
                           wt <- wilcox.test(y ~ x, data = test.df)
                           p <- wt$p.value
                           print(wt)
                         }
                         num.unique.samples <- length(unique(common.samples))
                         num.samples <- length(common.samples)
                         res <- c(gene, p, num.unique.samples, num.samples, num.unique.resp.samples, num.resp.samples)
                         names(res) <- c("gene", "pval", "num.uniq.samples", "num.samples", "num.uniq.resp.samples", "num.resp.samples")
                         res
                       })
               })
  ret
}

plot.genomic.vs.ic50 <- function(drug.df, drug, gene, drug.id.col, mutation.mat, specimen.id.col, patient.id.col) {
  df.drug <- drug.df[drug.df[,drug.id.col] == drug,]
  resp <- data.frame(IC50 = df.drug$IC50, PATIENT_ID = df.drug[,patient.id.col], id = df.drug[,specimen.id.col])
  resp <- resp[!is.infinite(resp$IC50) & !is.nan(resp$IC50) & !is.na(resp$IC50),]
  common.samples <- intersect(resp$id, colnames(mutation.mat))
  flag <- resp$id %in% common.samples
  resp <- resp[flag,]
  common.samples <- as.character(resp$id)
  gene.mutation <- factor(as.vector(unlist(mutation.mat[gene, common.samples])))
  df <- data.frame(mutation = gene.mutation, resp = resp$IC50)  
  ggplot(df, aes(x = mutation, y = resp)) + 
    geom_beeswarm() +
    xlab(paste0(gene, "\nmutation status")) + ylab(paste0("Log10 ", drug, " IC50"))
}



## BEGIN load OHSU dna data


## beat.mafs <- read.mafs(maf.names = ohsu.beat.maf.names, maf.tbl = ohsu.beat.mafs, mutation.types.to.total = c("Missense_Mutation", "Nonsense_Mutation"))
beat.mafs <- read.mafs(maf.names = ohsu.beat.maf.names, maf.tbl = ohsu.beat.mafs, mutation.types.to.total = c("Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Splice_Site", "Translation_Start_Site"))
## Make the samples the columns
beat.mafs <- t(beat.mafs)

## END load OHSU dna data 


## HERE


cat("Fitting FIMM curves\n")
## This works with plyr 1.7.1 and 1.8.1, but not 1.8.4
fimm.fits <- dlply(fimm.raw.dr.subset[,c("DRUG_ID", "DRUG_NAME", "CONCENTRATION", "PERCENT_INHIBITION", "SCREEN_ID", "SAMPLE_DATE", "PATIENT_ID")],
                   c("SCREEN_ID", "DRUG_ID"),
                   .parallel = TRUE,
                   .fun = 
                     function(df) {
                       print(nrow(df))
                       id <- unique(df[, !(colnames(df) %in% c("CONCENTRATION", "PERCENT_INHIBITION")),drop=F])
                       if(nrow(df) < 2) {
                         return(list(id=id, tbl=NULL))
                       }
                       print(df)
                       dose <- df$CONCENTRATION
                       inhibition <- df$PERCENT_INHIBITION
                       patient.id <- df$PATIENT_ID[1]
                       inhibitor <- df$DRUG_NAME[1]
                       file <- paste0(patient.id, "-", inhibitor)
                       tbl <- NULL
                       fit <- NULL
                       res <- NULL
                       res <- suppressWarnings(drc.calc.ic50.robust(dose, inhibition, file, patient.id, inhibitor, return.fit = FALSE, return.residuals = FALSE, apply.constraints = FALSE, plot.fit = FALSE, robust = "mean", force.min.to.zero = FALSE))
                       if(!is.null(res) & !is.null(res$tbl)) {
                         tbl <- res$tbl
                       }
                       if(!is.null(res) & !is.null(res$res)) {
                         tbl <- res$tbl
                       }
                       list(id=id, tbl=tbl)
                    })


## WORKING

cat("Fitting OHSU curves\n")
## Fit curves to OHSU data (AML only, those drugs that overlap with FIMM drugs)
ohsu.inh.tbl.subset <- subset(ohsu.inh.tbl, (ohsu.inh.tbl$inhibitor %in% ohsu.fimm.drugs$ID_Drug.ohsu))
rm(ohsu.inh.tbl); gc()

ohsu.inh.tbl.subset <- ohsu.inh.tbl.subset[order(ohsu.inh.tbl.subset$inhibitor, ohsu.inh.tbl.subset$patient_id, ohsu.inh.tbl.subset$lab_id),]
ohsu.fits <- dlply(ohsu.inh.tbl.subset[,c("inhibitor", "well_concentration", "normalized_viability", "time_of_read", "lab_id", "replicant", "patient_id")],
                   c("patient_id", "inhibitor", "lab_id", "replicant", "time_of_read"),
                   .parallel = TRUE,
                   .fun = 
                     function(df) {
                       print(nrow(df))
                       id <- unique(df[, !(colnames(df) %in% c("well_concentration", "normalized_viability")),drop=F])
                       if(nrow(df) < 2) {
                         return(list(id=id, tbl=NULL))
                       }
                       print(df)
                       dose <- df$well_concentration
                       inhibition <- 100 - df$normalized_viability
                       patient.id <- df$patient_id[1]
                       inhibitor <- df$inhibitor[1]
                       file <- paste0(patient.id, "-", inhibitor)
                       tbl <- NULL
                       fit <- NULL
                       res <- NULL
                       res <- suppressWarnings(drc.calc.ic50.robust(dose, inhibition, file, patient.id, inhibitor, return.fit = FALSE, return.residuals = FALSE, apply.constraints = FALSE, plot.fit = FALSE, robust = "mean", force.min.to.zero = FALSE))
                       if(!is.null(res) & !is.null(res$tbl)) {
                         tbl <- res$tbl
                       }
                       if(!is.null(res) & !is.null(res$res)) {
                         tbl <- res$tbl
                       }
                       list(id=id, tbl=tbl)
                     })

cat("Done\n")
save.image(".RData")
stop("stop")

## Collect the results into a single table
tmp <- lapply(fimm.fits, function(x) x$tbl)
fimm.fits.tbl <- do.call("rbind", tmp)

tmp <- lapply(fimm.fits, function(x) x$id)
fimm.fits.id <- do.call("rbind", tmp)

fimm.fits.df <- cbind(fimm.fits.tbl, fimm.fits.id)
if(!(all(rownames(fimm.fits.tbl) == rownames(fimm.fits.id)))) {
  warning("UNEXPECTED fimm.fits.tbl rownames != fimm.fits.id rownames\n")
}
if(!(all(rownames(fimm.fits.tbl) == rownames(fimm.fits.df)))) {
  warning("UNEXPECTED fimm.fits.tbl rownames != fimm.fits.df rownames\n")
}
rm(fimm.fits.tbl)
rm(fimm.fits.id)
gc()

tmp <- lapply(ohsu.fits, function(x) x$tbl)
ohsu.fits.tbl <- do.call("rbind", tmp)

tmp <- lapply(ohsu.fits, function(x) x$id)
ohsu.fits.id <- do.call("rbind", tmp)

ohsu.fits.df <- cbind(ohsu.fits.tbl, ohsu.fits.id)
if(!(all(rownames(ohsu.fits.tbl) == rownames(ohsu.fits.id)))) {
  warning("UNEXPECTED ohsu.fits.tbl rownames != ohsu.fits.id rownames\n")
}
if(!(all(rownames(ohsu.fits.tbl) == rownames(ohsu.fits.df)))) {
  warning("UNEXPECTED ohsu.fits.tbl rownames != ohsu.fits.df rownames\n")
}
rm(ohsu.fits.tbl)
rm(ohsu.fits.id)
rm(tmp)
gc()

## Translate gene symbols to ensg identifiers
symbols.to.ensg.mapping <- function(symbols) {
  # Get a mapping from gene symbol to ensembl id
  ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
  bm <- getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id'), 
              filters = 'hgnc_symbol', 
              values = symbols, 
              mart = ensembl)
  names(bm) <- c("SYMBOL", "ENSG")
  bm <- bm[!(bm$ENSG %in% c("")),]
  bm
}

all.targets <- unique(unlist(lapply(ohsu.fimm.drugs$Gene.Targets, function(x) unlist(strsplit(x, split=",[ ]*")))))
all.targets <- all.targets[!is.na(all.targets)]

sym.to.ensg <- symbols.to.ensg.mapping(all.targets)

ohsu.fimm.drugs$Ensg.Targets <- unlist(lapply(ohsu.fimm.drugs$Gene.Targets, 
                                           function(target.str) {
                                             targets <- unlist(strsplit(target.str, split=",[ ]*"))
                                             if(!(any(targets %in% sym.to.ensg$SYMBOL))) { 
                                               return(NA)
                                             }
                                             targets <- targets[targets %in% sym.to.ensg$SYMBOL]
                                             targets <- unique(sym.to.ensg$ENSG[sym.to.ensg$SYMBOL %in% targets])
                                             paste(targets, collapse=", ")
                                           }))

## Look at correlation of drug IC50 with drug target gene expression for both
## FIMM and OHSU.  Compare p-values between the two.

old.nrow <- nrow(ohsu.fits.df)
ohsu.fits.df <- merge(ohsu.fits.df, unique(ohsu.fimm.drugs[,c("ID_Drug", "ID_Drug.ohsu", "Gene.Targets", "Ensg.Targets")]), by.x = "inhibitor", by.y = "ID_Drug.ohsu")
if(old.nrow != nrow(ohsu.fits.df)) {
  warning("UNEXPECTED size of ohsu.fits.df changed\n")
}

old.nrow <- nrow(fimm.fits.df)
fimm.fits.df <- merge(fimm.fits.df, unique(ohsu.fimm.drugs[,c("ID_Drug", "ID_Drug.ohsu", "Gene.Targets", "Ensg.Targets")]), by.x = "DRUG_ID", by.y = "ID_Drug")
if(old.nrow != nrow(fimm.fits.df)) {
  warning("UNEXPECTED size of fimm.fits.df changed\n")
}

g1 <- ggplot(data = ohsu.fits.df, aes(x = gof))
g1 <- g1 + geom_density() + xlab("OHSU GOF")
g2 <- ggplot(data = fimm.fits.df, aes(x = gof))
g2 <- g2 + geom_density() + xlab("FIMM GOF")
pdf("gof-density.pdf")
grid.arrange(g1, g2)
d <- dev.off()

all.ensg.targets <- unique(unlist(lapply(ohsu.fimm.drugs$Ensg.Targets, function(x) unlist(strsplit(x, split=",[ ]*")))))
all.ensg.targets <- all.ensg.targets[!is.na(all.ensg.targets)]

common.genes <- intersect(rownames(fimm.expr), rownames(ohsu.expr))
common.genes <- intersect(common.genes, all.ensg.targets)

fimm.expr.common <- fimm.expr[common.genes,]
ohsu.expr.common <- ohsu.expr[common.genes,]

## Correlate gene target expression with target IC50 for FIMM
df <- fimm.fits.df[fimm.fits.df$DRUG_ID %in% common.drugs.tested,]
flag <- unlist(lapply(df$Ensg.Targets, function(x) any(unlist(strsplit(x, split=",[ ]*")) %in% common.genes)))
df <- df[flag,]
fimm.expr.df <- merge(df, ohsu.fimm.drugs[,c("ID_Drug", "Gene.Targets")], by.x = "DRUG_ID", by.y = "ID_Drug")
fimm.expr.corr <- ddply(fimm.expr.df, c("DRUG_ID"), .parallel = TRUE,
                        .fun = function(df.drug) {
                          drug.targets <- df.drug$Ensg.Targets[1]
                          drug.targets <- unlist(strsplit(drug.targets, split=",[ ]*"))
                          drug.targets <- drug.targets[drug.targets %in% rownames(fimm.expr.common)]
                          resp <- data.frame(IC50 = df.drug$IC50, PATIENT_ID = df.drug$PATIENT_ID, id = df.drug$SCREEN_ID)
                          rownames(resp) <- resp$id
                          resp <- resp[!is.infinite(resp$IC50) & !is.nan(resp$IC50) & !is.na(resp$IC50),]
                          common.samples <- intersect(rownames(resp), colnames(fimm.expr.common))
                          resp <- resp[common.samples,]
                          tmp <- as.data.frame(table(resp$PATIENT_ID))
                          colnames(tmp) <- c("PATIENT_ID", "weight")
                          tmp$weight <- 1/sqrt(tmp$weight)
                          resp <- merge(resp, tmp)
                          rownames(resp) <- resp$id
                          if(!(all(common.samples %in% rownames(resp)))) {
                            warning("UNEXPECTED row names")
                          }
                          resp <- resp[common.samples,]
                          ldply(drug.targets, .parallel = FALSE,
                                .fun = function(gene) {
                                  if(!(gene %in% rownames(fimm.expr.common))) {
                                    warning(paste0("UNEXPECTED: ", gene, " not in fimm.expr\n"))
                                  }
                                  if(!(all(common.samples %in% colnames(fimm.expr.common)))) {
                                    warning("UNEXPECTED columns not in fimm expr\n")
                                  }
                                  gene.expr <- fimm.expr.common[gene, common.samples]
                                  lm.fit <- lm(resp$IC50 ~ gene.expr, weights = resp$weight)
                                  sum <- summary(lm.fit)
                                  f <- sum$fstatistic
                                  r2 <- sum$r.squared
                                  p <- pf(f[1],f[2],f[3],lower.tail=F)
                                  res <- c(gene, p, r2)
                                  num.unique.samples <- length(unique(common.samples))
                                  num.samples <- length(common.samples)
                                  res <- c(gene, p, r2, num.unique.samples, num.samples)
                                  names(res) <- c("gene", "pval", "r2", "num.uniq.samples", "num.samples")
                                  res
                                })
                        })

plot.fimm.expr.vs.ic50 <- function(df, drug, ensg.gene, expr.mat, xlab, ylab) {
  df.drug <- subset(df, DRUG_ID == drug)
  resp <- data.frame(IC50 = df.drug$IC50, PATIENT_ID = df.drug$PATIENT_ID, id = df.drug$SCREEN_ID)
  rownames(resp) <- resp$id
  resp <- resp[!is.infinite(resp$IC50) & !is.nan(resp$IC50) & !is.na(resp$IC50),]
  common.samples <- intersect(rownames(resp), colnames(expr.mat))
  resp <- resp[common.samples,]
  tmp <- as.data.frame(table(resp$PATIENT_ID))
  colnames(tmp) <- c("PATIENT_ID", "weight")
  tmp$weight <- 1/sqrt(tmp$weight)
  resp <- merge(resp, tmp)
  rownames(resp) <- resp$id
  if(!(all(common.samples %in% rownames(resp)))) {
    warning("UNEXPECTED row names")
  }
  if(!(ensg.gene %in% rownames(expr.mat))) {
    warning(paste0("UNEXPECTED: ", ensg.gene, " not in expr.mat\n"))
  }
  if(!(all(common.samples %in% colnames(expr.mat)))) {
    warning("UNEXPECTED columns not in expr.mat\n")
  }
  gene.expr <- fimm.expr.common[ensg.gene, common.samples]
  df <- data.frame(expr = gene.expr, resp = resp$IC50, weights=resp$weight)
  ggplot(df, aes(x = expr, y = resp, weight=weights)) + 
    geom_point() + 
    xlab(xlab) + ylab(ylab) +
    stat_smooth(method = "lm", col = "red")  
}


## Correlate gene target expression with target IC50 for OHSU
df <- ohsu.fits.df[ohsu.fits.df$ID_Drug %in% common.drugs.tested,]
flag <- unlist(lapply(df$Ensg.Targets, function(x) any(unlist(strsplit(x, split=",[ ]*")) %in% common.genes)))
df <- df[flag,]
df <- merge(df, ohsu.fimm.drugs[,c("ID_Drug", "Gene.Targets")], by.x = "ID_Drug", by.y = "ID_Drug")
## What is this lab_id
## Need to translate lab_id to seq_id
old.nrow <- nrow(df)
new.df <- merge(df, unique(aml.rnaseq.sample.summary.tbl[,c("SeqID", "Original_LabID")]), by.x = "lab_id", by.y = "Original_LabID")
ohsu.expr.df <- new.df
ohsu.expr.corr <- ddply(ohsu.expr.df, c("ID_Drug"), .parallel = TRUE,
                        .fun = function(df.drug) {
                          drug.targets <- df.drug$Ensg.Targets[1]
                          drug.targets <- unlist(strsplit(drug.targets, split=",[ ]*"))
                          drug.targets <- drug.targets[drug.targets %in% rownames(ohsu.expr.common)]
                          ids <- df.drug$SeqID
                          ids <- gsub(pattern="-", x=ids, replacement=".")
                          ids <- unlist(lapply(ids, function(x) paste0("X", x)))
                          resp <- data.frame(IC50 = df.drug$IC50, PATIENT_ID = df.drug$patient_id, id = ids)
                          resp <- resp[!is.infinite(resp$IC50) & !is.nan(resp$IC50) & !is.na(resp$IC50),]
                          common.samples <- intersect(resp$id, colnames(ohsu.expr.common))
                          flag <- resp$id %in% common.samples
                          resp <- resp[flag,]
                          tmp <- as.data.frame(table(resp$PATIENT_ID))
                          colnames(tmp) <- c("PATIENT_ID", "weight")
                          tmp$weight <- 1/sqrt(tmp$weight)
                          resp <- merge(resp, tmp)
                          common.samples <- as.character(resp$id)
                          ldply(drug.targets, .parallel = FALSE,
                                .fun = function(gene) {
                                  if(!(gene %in% rownames(ohsu.expr.common))) {
                                    warning(paste0("UNEXPECTED: ", gene, " not in ohsu.expr\n"))
                                  }
                                  if(!(all(common.samples %in% colnames(ohsu.expr.common)))) {
                                    warning("UNEXPECTED columns not in ohsu expr\n")
                                  }
                                  gene.expr <- as.vector(unlist(ohsu.expr.common[gene, common.samples]))
                                  lm.fit <- lm(resp$IC50 ~ gene.expr, weights = resp$weight)
                                  sum <- summary(lm.fit)
                                  f <- sum$fstatistic
                                  r2 <- sum$r.squared
                                  p <- pf(f[1],f[2],f[3],lower.tail=F)
                                  num.unique.samples <- length(unique(common.samples))
                                  num.samples <- length(common.samples)
                                  res <- c(gene, p, r2, num.unique.samples, num.samples)
                                  names(res) <- c("gene", "pval", "r2", "num.uniq.samples", "num.samples")
                                  res
                                })
                        })

plot.ohsu.expr.vs.ic50 <- function(df, drug, ensg.gene, expr.mat, xlab, ylab) {
  df.drug <- subset(df, ID_Drug == drug)
  ids <- df.drug$SeqID
  ids <- gsub(pattern="-", x=ids, replacement=".")
  ids <- unlist(lapply(ids, function(x) paste0("X", x)))
  resp <- data.frame(IC50 = df.drug$IC50, PATIENT_ID = df.drug$patient_id, id = ids)
  resp <- resp[!is.infinite(resp$IC50) & !is.nan(resp$IC50) & !is.na(resp$IC50),]
  common.samples <- intersect(resp$id, colnames(expr.mat))
  flag <- resp$id %in% common.samples
  resp <- resp[flag,]
  tmp <- as.data.frame(table(resp$PATIENT_ID))
  colnames(tmp) <- c("PATIENT_ID", "weight")
  tmp$weight <- 1/sqrt(tmp$weight)
  resp <- merge(resp, tmp)
  common.samples <- as.character(resp$id)
  if(!(ensg.gene %in% rownames(expr.mat))) {
    warning(paste0("UNEXPECTED: ", ensg.gene, " not in expr.mat\n"))
  }
  if(!(all(common.samples %in% colnames(expr.mat)))) {
    warning("UNEXPECTED columns not in expr.mat\n")
  }
  gene.expr <- as.vector(unlist(expr.mat[ensg.gene, common.samples]))
  df <- data.frame(expr = gene.expr, resp = resp$IC50, weights=resp$weight)
  ggplot(df, aes(x = expr, y = resp, weight=weights)) + 
    geom_point() + 
    xlab(xlab) + ylab(ylab) +
    stat_smooth(method = "lm", col = "red")  
}


## Plot few cases that have significant correlation in both
both.sig <- subset(ohsu.fimm.expr.corr, (pval.ohsu < 0.05) & (pval.fimm < 0.05))
for(i in 1:nrow(both.sig)) {
  ensg.gene <- both.sig$gene[i]
  drug.id <- both.sig$ID_Drug[i]
  pval.ohsu <- both.sig$pval.ohsu[i]
  pval.fimm <- both.sig$pval.fimm[i]
  
  gene.sym <- sym.to.ensg$SYMBOL[sym.to.ensg$ENSG == ensg.gene]
  drug.name <- fimm.ohsu.drug.annotations$DRUG_NAME[fimm.ohsu.drug.annotations$ID_Drug == drug.id]
  g1 <- plot.fimm.expr.vs.ic50(fimm.expr.df, drug.id, ensg.gene, fimm.expr.common, paste0("Log2 CPM (", gene.sym, ")"), paste0("Log10 IC50 (", drug.name, ")"))
  g1 <- g1 + ggtitle(paste0("FIMM (pval = ", pval.fimm, ")"))
  g2 <- plot.ohsu.expr.vs.ic50(ohsu.expr.df, drug.id, ensg.gene, ohsu.expr.common, paste0("Log2 CPM (", gene.sym, ")"), paste0("Log10 IC50 (", drug.name, ")"))
  g2 <- g2 + ggtitle(paste0("OHSU (pval = ", pval.ohsu, ")"))
  pdf(paste0("both-sig-", ensg.gene, "-", drug.id, ".pdf"))
  grid.arrange(g1, g2)
  d <- dev.off()
}

## Plot the top hits in each
sig <- subset(ohsu.fimm.expr.corr, pval.ohsu < 0.01)
for(i in 1:nrow(sig)) {
  ensg.gene <- sig$gene[i]
  drug.id <- sig$ID_Drug[i]
  pval.ohsu <- sig$pval.ohsu[i]
  pval.fimm <- sig$pval.fimm[i]
  
  gene.sym <- sym.to.ensg$SYMBOL[sym.to.ensg$ENSG == ensg.gene]
  drug.name <- fimm.ohsu.drug.annotations$DRUG_NAME[fimm.ohsu.drug.annotations$ID_Drug == drug.id]
  g1 <- plot.fimm.expr.vs.ic50(fimm.expr.df, drug.id, ensg.gene, fimm.expr.common, paste0("Log2 CPM (", gene.sym, ")"), paste0("Log10 IC50 (", drug.name, ")"))
  g1 <- g1 + ggtitle(paste0("FIMM (pval = ", pval.fimm, ")"))
  g2 <- plot.ohsu.expr.vs.ic50(ohsu.expr.df, drug.id, ensg.gene, ohsu.expr.common, paste0("Log2 CPM (", gene.sym, ")"), paste0("Log10 IC50 (", drug.name, ")"))
  g2 <- g2 + ggtitle(paste0("OHSU (pval = ", pval.ohsu, ")"))
  pdf(paste0("ohsu-sig-", ensg.gene, "-", drug.id, ".pdf"))
  grid.arrange(g1, g2)
  d <- dev.off()
}

sig <- subset(ohsu.fimm.expr.corr, pval.fimm < 0.01)
for(i in 1:nrow(sig)) {
  ensg.gene <- sig$gene[i]
  drug.id <- sig$ID_Drug[i]
  pval.ohsu <- sig$pval.ohsu[i]
  pval.fimm <- sig$pval.fimm[i]
  
  gene.sym <- sym.to.ensg$SYMBOL[sym.to.ensg$ENSG == ensg.gene]
  drug.name <- fimm.ohsu.drug.annotations$DRUG_NAME[fimm.ohsu.drug.annotations$ID_Drug == drug.id]
  g1 <- plot.fimm.expr.vs.ic50(fimm.expr.df, drug.id, ensg.gene, fimm.expr.common, paste0("Log2 CPM (", gene.sym, ")"), paste0("Log10 IC50 (", drug.name, ")"))
  g1 <- g1 + ggtitle(paste0("FIMM (pval = ", pval.fimm, ")"))
  g2 <- plot.ohsu.expr.vs.ic50(ohsu.expr.df, drug.id, ensg.gene, ohsu.expr.common, paste0("Log2 CPM (", gene.sym, ")"), paste0("Log10 IC50 (", drug.name, ")"))
  g2 <- g2 + ggtitle(paste0("OHSU (pval = ", pval.ohsu, ")"))
  pdf(paste0("fimm-sig-", ensg.gene, "-", drug.id, ".pdf"))
  grid.arrange(g1, g2)
  d <- dev.off()
}


## Plot the best cases in both

## - for intersected drugs, use symbols/ensg for genes


## - look at correlation of expr in gene targets with IC50
## - look at correlation of mutant in gene targets with IC50 (only those that have a lot of mutations)

## - 10-fold cross validation of ohsu using expr
## - 10-fold cross validation of ohsu using genomic
