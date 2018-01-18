prepare.ohsu.drug.response.and.expr.matrices <- function(drc.df, expr.df, rnaseq.summary.df = NULL, drugs, response.col) {

  drc.df <- subset(drc.df, inhibitor %in% drugs)

  ## Convert lab_id as used in drug response to SeqIDs used in expr matrix.
  if(!is.null(rnaseq.summary.df)) {
    drc.df <- merge(drc.df, rnaseq.summary.df, by.x = "lab_id", by.y = "Original_LabID")
  }

  ## Exclude any NAs in response
  flag <- !is.na(drc.df[, response.col])
  drc.df <- drc.df[flag, ]

  ## Take the median of replicates in the X and y data
  drc.df <- ddply(drc.df, c("inhibitor", "SeqID"), .fun = function(df) median(df[, response.col]))
  colnames(drc.df) <- c("inhibitor", "SeqID", response.col)

  ## Convert drug response table into a matrix in which columns are samples and
  ## rows are drug response
  drc.df <- spread_(drc.df, "SeqID", response.col)  
  rownames(drc.df) <- drc.df$inhibitor
  drc.df <- drc.df[, !(colnames(drc.df) == "inhibitor")]

  ## Restrict X and response to be over the same samples
  if(!is.null(expr.df)) {
    common.samples <- intersect(colnames(drc.df), colnames(expr.df))
    drc.df <- drc.df[, common.samples]
    expr.df <- expr.df[, common.samples]
  }

  return(list(drc.df = drc.df, expr.df = expr.df))
}

prepare.fimm.drug.response.and.expr.matrices <- function(drc.df, expr.df, drugs, response.col) {

  drug.col <- "DRUG_ID"
  patient.col <- "SCREEN_ID"

  flag <- drc.df[, drug.col] %in% drugs
  drc.df <- drc.df[flag, ]

  ## Exclude any NAs in response
  flag <- !is.na(drc.df[, response.col])
  drc.df <- drc.df[flag, ]

  ## Take the median of replicates in the X and y data
  drc.df <- drc.df[, c(drug.col, patient.col, response.col)]

  ## Convert drug response table into a matrix in which columns are samples and
  ## rows are drug response
  drc.df <- spread_(drc.df, patient.col, response.col)  
  rownames(drc.df) <- drc.df[, drug.col]
  drc.df <- drc.df[, !(colnames(drc.df) == drug.col)]

  ## Restrict X and response to be over the same samples
  if(!is.null(expr.df)) {
    common.samples <- intersect(colnames(drc.df), colnames(expr.df))
    drc.df <- drc.df[, common.samples]
    expr.df <- expr.df[, common.samples]
  }

  return(list(drc.df = drc.df, expr.df = expr.df))
}

## For FIMM:
##   drug.col <- "DRUG_ID"
##   patient.col <- "SCREEN_ID"
## For OHSU:
##   drug.col <- "inhibitor"
##   patient.col <- "SeqID"
## For CTRP:
##   drug.col <- "master_cpd_id"
##   patient.col <- "ccl_name"
## For NTAP:
##   drug.col <- "SMI"
##   patient.col <- "sampleIdentifier"
prepare.drug.response.matrix <- function(drc.df, drug.col, patient.col, response.col, drugs = NULL) {

  if(!is.null(drugs)) { 
    flag <- drc.df[, drug.col] %in% drugs
    drc.df <- drc.df[flag, ]
  }

  ## Exclude any NAs in response
  flag <- !is.na(drc.df[, response.col])
  drc.df <- drc.df[flag, ]

  ## Take the median of replicates in the X and y data
  drc.df <- drc.df[, c(drug.col, patient.col, response.col)]
  drc.df <- ddply(drc.df, c(drug.col, patient.col), .fun = function(df) median(df[, response.col]))
  colnames(drc.df) <- c(drug.col, patient.col, response.col)

  ## Convert drug response table into a matrix in which columns are samples and
  ## rows are drug response
  drc.df <- spread_(drc.df, patient.col, response.col)  
  rownames(drc.df) <- drc.df[, drug.col]
  drc.df <- drc.df[, !(colnames(drc.df) == drug.col)]

  drc.df
}

prepare.drug.response.and.expr.matrices <- function(drc.df, expr.df, drugs, drug.col, patient.col, response.col) {

  drc.df <- prepare.drug.response.matrix(drc.df, drug.col, patient.col, response.col, drugs = drugs)
  ## Restrict X and response to be over the same samples
  if(!is.null(expr.df) && !is.na(expr.df)) {
    colnames(drc.df) <- make.names(colnames(drc.df))
    colnames(expr.df) <- make.names(colnames(expr.df))
    common.samples <- intersect(colnames(drc.df), colnames(expr.df))
    cat(paste0("Samples with both expression and drug assays: ", length(common.samples), "\n"))
    drc.df <- drc.df[, as.character(common.samples)]
    expr.df <- expr.df[, as.character(common.samples)]
  }

  return(list(drc.df = drc.df, expr.df = expr.df))

  if(FALSE) {
  flag <- drc.df[, drug.col] %in% drugs
  drc.df <- drc.df[flag, ]

  ## Exclude any NAs in response
  flag <- !is.na(drc.df[, response.col])
  drc.df <- drc.df[flag, ]

  ## Take the median of replicates in the X and y data
  drc.df <- drc.df[, c(drug.col, patient.col, response.col)]
  drc.df <- ddply(drc.df, c(drug.col, patient.col), .fun = function(df) median(df[, response.col]))
  colnames(drc.df) <- c(drug.col, patient.col, response.col)

  ## Convert drug response table into a matrix in which columns are samples and
  ## rows are drug response
  drc.df <- spread_(drc.df, patient.col, response.col)  
  rownames(drc.df) <- drc.df[, drug.col]
  drc.df <- drc.df[, !(colnames(drc.df) == drug.col)]

  ## Restrict X and response to be over the same samples
  if(!is.null(expr.df)) {
    common.samples <- intersect(colnames(drc.df), colnames(expr.df))
    drc.df <- drc.df[, common.samples]
    expr.df <- expr.df[, common.samples]
  }

  return(list(drc.df = drc.df, expr.df = expr.df))
  }
}

