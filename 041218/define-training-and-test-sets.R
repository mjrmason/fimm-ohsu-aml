common.drugs.by.data.set <- list()
for(ds in names(fits)) {
  common.drugs.by.data.set[[ds]] <- drug.name.tbl[, drug.name.col]
}

train.set.names <- list()
train.dss.args <- list()
train.expr.args <- list()
train.genomic.args <- list()
train.clinical.args <- list()
train.common.drugs <- list()
train.drugs <- list()
train.drug.cols <- list()
train.patient.cols <- list()
train.response.cols <- list()

test.set.names <- list()
test.dss.args <- list()
test.expr.args <- list()
test.genomic.args <- list()
test.clinical.args <- list()
test.common.drugs <- list()
test.drug.cols <- list()
test.patient.cols <- list()
test.response.cols <- list()

ds <- "ohsu"
if(ds %in% data.sets.to.analyze) {

  nxt.indx <- length(test.set.names) + 1
  test.set.names[[nxt.indx]] <- ds
  test.dss.args[[nxt.indx]] <- fits[[ds]]
  test.expr.args[[nxt.indx]] <- expr.filt.matrices[[ds]][common.genes, ]
  test.genomic.args[nxt.indx] <- list(NULL)
  test.clinical.args[nxt.indx] <- list(NULL)
  test.common.drugs[[nxt.indx]] <- common.drugs.by.data.set[[ds]]
  test.drug.cols[[nxt.indx]] <- drug.name.col
  test.patient.cols[[nxt.indx]] <- expr.patient.id.cols[[ds]]
  test.response.cols[[nxt.indx]] <- "dss.auc.ll4"

  nxt.indx <- length(train.set.names) + 1
  train.set.names[[nxt.indx]] <- ds
  train.dss.args[[nxt.indx]] <- fits[[ds]]
  train.expr.args[[nxt.indx]] <- expr.filt.matrices[[ds]][common.genes, ] 
  train.genomic.args[nxt.indx] <- list(NULL)
  train.clinical.args[nxt.indx] <- list(NULL)
  train.common.drugs[[nxt.indx]] <- common.drugs.by.data.set[[ds]]
  train.drugs[[nxt.indx]] <- common.drugs.by.data.set[[ds]]
  train.drug.cols[[nxt.indx]] <- drug.name.col
  train.patient.cols[[nxt.indx]] <- expr.patient.id.cols[[ds]]
  train.response.cols[[nxt.indx]] <- "dss.auc.ll4"

  set.seed(1234)
  half.index <- createDataPartition(fits[[ds]][, data.set.drug.id.cols[[ds]]], p = .5, list = FALSE)
  fits.set1 <- fits[[ds]][half.index, ]
  fits.set2 <- fits[[ds]][-half.index, ]
if(FALSE) {
  nxt.indx <- length(train.set.names) + 1
  train.set.names[[nxt.indx]] <- "ohsu.set1"
  train.dss.args[[nxt.indx]] <- fits.set1
  train.expr.args[[nxt.indx]] <- expr.filt.matrices[[ds]][common.genes, ] 
  train.genomic.args[nxt.indx] <- list(NULL)
  train.clinical.args[nxt.indx] <- list(NULL)
  train.common.drugs[[nxt.indx]] <- common.drugs.by.data.set[[ds]]
  train.drugs[[nxt.indx]] <- common.drugs.by.data.set[[ds]]
  train.drug.cols[[nxt.indx]] <- drug.name.col
  train.patient.cols[[nxt.indx]] <- expr.patient.id.cols[[ds]]
  train.response.cols[[nxt.indx]] <- "dss.auc.ll4"

  nxt.indx <- length(train.set.names) + 1
  train.set.names[[nxt.indx]] <- "ohsu.set2"
  train.dss.args[[nxt.indx]] <- fits.set2
  train.expr.args[[nxt.indx]] <- expr.filt.matrices[[ds]][common.genes, ] 
  train.genomic.args[nxt.indx] <- list(NULL)
  train.clinical.args[nxt.indx] <- list(NULL)
  train.common.drugs[[nxt.indx]] <- common.drugs.by.data.set[[ds]]
  train.drugs[[nxt.indx]] <- common.drugs.by.data.set[[ds]]
  train.drug.cols[[nxt.indx]] <- drug.name.col
  train.patient.cols[[nxt.indx]] <- expr.patient.id.cols[[ds]]
  train.response.cols[[nxt.indx]] <- "dss.auc.ll4"

  nxt.indx <- length(test.set.names) + 1
  test.set.names[[nxt.indx]] <- "ohsu.set2"
  test.dss.args[[nxt.indx]] <- fits.set2
  test.expr.args[[nxt.indx]] <- expr.filt.matrices[[ds]][common.genes, ] 
  test.genomic.args[nxt.indx] <- list(NULL)
  test.clinical.args[nxt.indx] <- list(NULL)
  test.common.drugs[[nxt.indx]] <- common.drugs.by.data.set[[ds]]
  test.drug.cols[[nxt.indx]] <- drug.name.col
  test.patient.cols[[nxt.indx]] <- expr.patient.id.cols[[ds]]
  test.response.cols[[nxt.indx]] <- "dss.auc.ll4"
}

} ## if("ohsu" %in% data.sets.to.analyze) 

ds <- "fimm"
if(ds %in% data.sets.to.analyze) {

  nxt.indx <- length(test.set.names) + 1
  test.set.names[[nxt.indx]] <- ds
  test.dss.args[[nxt.indx]] <- fits[[ds]]
  test.expr.args[[nxt.indx]] <- expr.filt.matrices[[ds]][common.genes, ]
  test.genomic.args[nxt.indx] <- list(NULL)
  test.clinical.args[nxt.indx] <- list(NULL)
  test.common.drugs[[nxt.indx]] <- common.drugs.by.data.set[[ds]]
  test.drug.cols[[nxt.indx]] <- drug.name.col
  test.patient.cols[[nxt.indx]] <- expr.patient.id.cols[[ds]]
  test.response.cols[[nxt.indx]] <- "dss.auc.ll4"

  nxt.indx <- length(train.set.names) + 1
  train.set.names[[nxt.indx]] <- ds
  train.dss.args[[nxt.indx]] <- fits[[ds]]
  train.expr.args[[nxt.indx]] <- expr.filt.matrices[[ds]][common.genes, ]
  train.genomic.args[nxt.indx] <- list(NULL)
  train.clinical.args[nxt.indx] <- list(NULL)
  train.common.drugs[[nxt.indx]] <- common.drugs.by.data.set[[ds]]
  train.drugs[[nxt.indx]] <- common.drugs.by.data.set[[ds]]
  train.drug.cols[[nxt.indx]] <- drug.name.col
  train.patient.cols[[nxt.indx]] <- expr.patient.id.cols[[ds]]
  train.response.cols[[nxt.indx]] <- "dss.auc.ll4"

  if(FALSE) {
  set.seed(1234)
  half.index <- createDataPartition(fits[[ds]][, data.set.drug.id.cols[[ds]]], p = .5, list = FALSE)
  fits.set1 <- fits[[ds]][half.index, ]
  fits.set2 <- fits[[ds]][-half.index, ]

  nxt.indx <- length(train.set.names) + 1
  train.set.names[[nxt.indx]] <- "fimm.set1"
  train.dss.args[[nxt.indx]] <- fits.set1
  train.expr.args[[nxt.indx]] <- expr.filt.matrices[[ds]][common.genes, ] 
  train.genomic.args[nxt.indx] <- list(NULL)
  train.clinical.args[nxt.indx] <- list(NULL)
  train.common.drugs[[nxt.indx]] <- common.drugs.by.data.set[[ds]]
  train.drugs[[nxt.indx]] <- common.drugs.by.data.set[[ds]]
  train.drug.cols[[nxt.indx]] <- drug.name.col
  train.patient.cols[[nxt.indx]] <- expr.patient.id.cols[[ds]]
  train.response.cols[[nxt.indx]] <- "dss.auc.ll4"

  nxt.indx <- length(train.set.names) + 1
  train.set.names[[nxt.indx]] <- "fimm.set2"
  train.dss.args[[nxt.indx]] <- fits.set2
  train.expr.args[[nxt.indx]] <- expr.filt.matrices[[ds]][common.genes, ] 
  train.genomic.args[nxt.indx] <- list(NULL)
  train.clinical.args[nxt.indx] <- list(NULL)
  train.common.drugs[[nxt.indx]] <- common.drugs.by.data.set[[ds]]
  train.drugs[[nxt.indx]] <- common.drugs.by.data.set[[ds]]
  train.drug.cols[[nxt.indx]] <- drug.name.col
  train.patient.cols[[nxt.indx]] <- expr.patient.id.cols[[ds]]
  train.response.cols[[nxt.indx]] <- "dss.auc.ll4"
  }

} ## if("fimm" %in% data.sets.to.analyze) 

ds <- "ctrp"
if(ds %in% data.sets.to.analyze) {

  nxt.indx <- length(train.set.names) + 1
  train.set.names[[nxt.indx]] <- ds
  train.dss.args[[nxt.indx]] <- fits[[ds]]
  train.expr.args[[nxt.indx]] <- expr.filt.matrices[[ds]][common.genes, ]
  train.genomic.args[nxt.indx] <- list(NULL)
  train.clinical.args[nxt.indx] <- list(NULL)
  train.common.drugs[[nxt.indx]] <- common.drugs.by.data.set[[ds]]
  train.drugs[[nxt.indx]] <- common.drugs.by.data.set[[ds]]
  train.drug.cols[[nxt.indx]] <- drug.name.col
  train.patient.cols[[nxt.indx]] <- expr.patient.id.cols[[ds]]
  train.response.cols[[nxt.indx]] <- "dss.auc.ll4"

} ## if("ctrp" %in% data.sets.to.analyze) 

ds <- "ntap"
if(ds %in% data.sets.to.analyze) {

  nxt.indx <- length(train.set.names) + 1
  train.set.names[[nxt.indx]] <- ds
  train.dss.args[[nxt.indx]] <- fits[[ds]]
  train.expr.args[[nxt.indx]] <- expr.filt.matrices[[ds]][common.genes, ]
  train.genomic.args[nxt.indx] <- list(NULL)
  train.clinical.args[nxt.indx] <- list(NULL)
  train.common.drugs[[nxt.indx]] <- common.drugs.by.data.set[[ds]]
  train.drugs[[nxt.indx]] <- common.drugs.by.data.set[[ds]]
  train.drug.cols[[nxt.indx]] <- drug.name.col
  train.patient.cols[[nxt.indx]] <- expr.patient.id.cols[[ds]]
  train.response.cols[[nxt.indx]] <- "dss.auc.ll4"
} ## if("ntap" %in% data.sets.to.analyze) 

ds <- "sanger"
if(ds %in% data.sets.to.analyze) {

  nxt.indx <- length(train.set.names) + 1
  train.set.names[[nxt.indx]] <- ds
  train.dss.args[[nxt.indx]] <- fits[[ds]]
  train.expr.args[[nxt.indx]] <- expr.filt.matrices[[ds]][common.genes, ]
  train.genomic.args[nxt.indx] <- list(NULL)
  train.clinical.args[nxt.indx] <- list(NULL)
  train.common.drugs[[nxt.indx]] <- common.drugs.by.data.set[[ds]]
  train.drugs[[nxt.indx]] <- common.drugs.by.data.set[[ds]]
  train.drug.cols[[nxt.indx]] <- drug.name.col
  train.patient.cols[[nxt.indx]] <- expr.patient.id.cols[[ds]]
  train.response.cols[[nxt.indx]] <- "dss.auc.ll4"
} ## if("sanger" %in% data.sets.to.analyze) 

