orig.fits <- list()
fits <- list()
exprs <- list()
orig.exprs <- list()
orig.drcs <- list()
prep <- list()
prep[["train.x"]] <- list()
prep[["train.drcs"]] <- list()

for(ds in c("fimm", "ohsu")) {
  orig.fits[[ds]] <- read.table(file = paste0(ds, "-orig-fits.tsv"), sep="\t", header = TRUE, as.is=TRUE)
  fits[[ds]] <- read.table(file = paste0(ds, "-fits.tsv"), sep="\t", header = TRUE, as.is=TRUE)
  orig.exprs[[ds]] <- read.table(file = paste0(ds, "-orig-exprs.tsv"), sep="\t", header = TRUE, as.is=TRUE)
  orig.drcs[[ds]] <- read.table(file = paste0(ds, "-orig-drcs.tsv"), sep="\t", header = TRUE, as.is=TRUE)
  exprs[[ds]] <- read.table(file = paste0(ds, "-exprs.tsv"), sep="\t", header = TRUE, as.is=TRUE)

  prep[["train.x"]][[ds]] <- list()
  prep[["train.x"]][[ds]][["gene"]] <- read.table(file = paste0(ds, "-prep-expr.tsv"), sep="\t", header = TRUE, as.is=TRUE)

  prep[["train.drcs"]][[ds]] <- read.table(file = paste0(ds, "-prep-drcs.tsv"), sep="\t", header = TRUE, as.is=TRUE)
}

ohsu.metadata <- read.table(file = "ohsu-metadata.tsv", sep="\t", header = TRUE, as.is=TRUE)
patient.age.tbl <- read.table(file = "ohsu-patient-age.tsv", sep="\t", header = TRUE, as.is=TRUE)
