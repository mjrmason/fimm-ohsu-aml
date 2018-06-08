## Load the map of drugs between data sets
obj <- synGet(drug.mapping.synId, downloadFile = TRUE)
drug.map <- read.table(getFileLocation(obj), sep="\t", header=TRUE, quote="\"", comment.char="")

## Drop any drugs that show up in multiple rows of the drug map--which would indicate an ambiguous mapping.
my.dup <- function(x) duplicated(x, fromLast = TRUE) | duplicated(x, fromLast = FALSE)
for(col in drug.map.drug.id.cols) {
  flag <- !my.dup(drug.map[, col])
  drug.map <- drug.map[flag, ]
}

process.fits <- function(fits, drug.map) {
  ## Restrict the fits to those drugs in common across the data sets -- do so by merging.
  cat("Restricting fits to those involving common drugs\n")
  for(ds in names(fits)) {
    fits[[ds]] <- merge(fits[[ds]], drug.map, by.x = data.set.drug.id.cols[[ds]], by.y = drug.map.drug.id.cols[[ds]], all = FALSE)
  }

  common.drugs <- Reduce("intersect", lapply(fits, function(fit) unique(fit[, drug.name.col])))
  for(i in 1:length(fits)) {
    ds <- names(fits)[i]
    flag <- fits[[ds]][, drug.name.col] %in% common.drugs
    fits[[ds]] <- fits[[ds]][flag, ]
  }


  ## Drop any fits that are flagged as outliers
  if(TRUE) {
    fct.postfixes <- list("LL.4" = ".ll4", "L.4" = ".l4")
    fct <- "LL.4"
    for(ds in names(fits)) {
        exclude.cols <- colnames(fits[[ds]])[grepl(pattern=paste0("exclude", fct.postfixes[[fct]]), x=colnames(fits[[ds]]))]
        any.excluded <- unlist(apply(fits[[ds]][, exclude.cols], 1, function(row) any(row)))
        orig.nrow <- nrow(fits[[ds]])
        fits[[ds]] <- fits[[ds]][!any.excluded, ]
        new.nrow <- nrow(fits[[ds]])
        cat(paste0(ds, ": filtered ", orig.nrow - new.nrow, " of the original ", orig.nrow, " fits, leaving ", new.nrow, " fits.\n"))
    }
  }

  ## Merge metadata to fit files
  cat("Merging metadata to fits\n")
  for(ds in names(fits)) {
    switch(ds,
      "ohsu" = {
        ## Add the SeqID to the OHSU drug data.
        ## Read in the rna-seq summary table, which provides a map from the lab_id's used
        ## in the drug response data to the SeqID's used in the expression data.
        synId <- "syn10083817"
        obj <- synGet(synId, downloadFile = TRUE)
        file <- getFileLocation(obj)
        ohsu.rnaseq.sample.summary <- read.table(file, header=TRUE, sep="\t")
        if(!("SeqID" %in% colnames(fits[[ds]]))) {
          fits[[ds]] <- merge(fits[[ds]], ohsu.rnaseq.sample.summary, by.x = "lab_id", by.y = "Original_LabID")
        }
      },
      "fimm" = {
      },
      "ctrp" = {
        fits[[ds]] <- merge(fits[[ds]], ctrp.cell.line.metadata[, c("master_ccl_id", "ccl_name")])
      },
      "ntap" = {
      },
      "sanger" = {
      },
      {
        stop("Unrecognized data set\n")
      })
  }
  fits
}

process.exprs <- function(exprs) {
  ## Sanitize some of the column headers
  for(nm in names(exprs)) {
    switch(nm,
      "ohsu" = {
        ## Convert the column names from X20.00347 to 20-00347
        colnames(exprs[[nm]]) <- gsub("X(.*)\\.(.*)", "\\1-\\2", colnames(exprs[[nm]]))
  
        ## Drop two outliers
        outliers <- c("14-00800", "20-00062")
##        if(!(all(outliers %in% colnames(exprs[[nm]])))) {
##          stop(paste0("Was expecting outliers: ", paste(outliers, collapse=","), " in expr matrix\n"))
##        }
        exprs[[nm]] <- exprs[[nm]][, !(colnames(exprs[[nm]]) %in% outliers)]
      },
      "fimm" = {
      },
      "ctrp" = {
      },
      "ntap" = {
      },
      "sanger" = {
        colnames(exprs[[nm]]) <- gsub("X(.*)", "\\1", colnames(exprs[[nm]]))
      },
      {
        stop("Unrecognized data set\n")
      })
  }

  common.genes <- Reduce("intersect", lapply(exprs, rownames))

  cat(paste0("Num expression features: ", length(unique(common.genes)), "\n"))

  ## Restrict the matrices to those genes in common across data sets
  for(i in 1:length(exprs)) {
    ds <- names(exprs)[i]
    exprs[[ds]] <- exprs[[ds]][common.genes, ]
  }
  exprs
}

orig.fits <- llply(data.set.drug.fit.synIds, .parallel = FALSE,
              .fun = function(synId) {
                       obj <- synGet(synId, downloadFile = TRUE)
                       file <- getFileLocation(obj)
                       read.table(file, sep="\t", header=TRUE, as.is=TRUE, comment.char="", quote="\"")
              })

orig.fits <- process.fits(orig.fits, drug.map)

orig.common.drugs <- Reduce("intersect", lapply(orig.fits, function(fit) unique(fit[, drug.name.col])))

## Limit the drug map to common drugs (which, by design of the drug map file, should already be the case)
orig.drug.name.tbl <- drug.map[drug.map[, drug.name.col] %in% orig.common.drugs,]

## Load in the expr files
cat("Reading expression data sets\n")
exprs <- llply(data.set.expr.synIds, .parallel = FALSE,
              .fun = function(synId) {
                       obj <- synGet(synId, downloadFile = TRUE)
                       file <- getFileLocation(obj)
                       read.table(file, sep="\t", header=TRUE, as.is=TRUE, comment.char="", quote="\"")
              })

cat("Reading original expression data sets\n")
orig.exprs <- llply(data.set.orig.expr.synIds, .parallel = FALSE,
              .fun = function(synId) {
                       obj <- synGet(synId, downloadFile = TRUE)
                       file <- getFileLocation(obj)
                       read.table(file, sep="\t", header=TRUE, as.is=TRUE, comment.char="", quote="\"")
              })

exprs <- process.exprs(exprs)
orig.exprs <- process.exprs(orig.exprs)

## We may have filtered samples out of the combat-normalized data (i.e., in exprs) but not in orig.exprs.
## Apply that filtering here.

for(ds in names(exprs)) {
  samples <- colnames(exprs[[ds]])
  orig.exprs[[ds]] <- orig.exprs[[ds]][, colnames(orig.exprs[[ds]]) %in% samples]
}

orig.drcs <- llply(data.sets,
                   .fun = function(ds) {
                            drug.col <- "OHSU_DRUG_NAME"
                            if(ds == "ohsu") { drug.col <- "inhibitor" }
                            prepare.drug.response.matrix(orig.fits[[ds]], drug.col = drug.col,
                                                         patient.col = patient.id.cols[[ds]], response.col = "auc.ll4")                            
                          })

## Restrict to drugs that are non-NA in a fractional number of rows and columns
orig.scaled.drcs <- llply(data.sets,
                   .fun = function(ds) {
                            subset.and.zscore.matrix(orig.drcs[[ds]], col.frac.cutoff = 0.25, row.frac.cutoff = 0.25)
                          })

## Get the drugs that are included in all matrices
non.na.drugs <- Reduce(intersect, lapply(orig.scaled.drcs, function(x) rownames(x)))

## Restrict to those drugs
orig.scaled.drcs <- llply(data.sets,
                   .fun = function(ds) {
                            orig.scaled.drcs[[ds]][non.na.drugs,]
                          })

