n.genes <- 5000
source("corr-of-corrs.R")
corrs <- list()
for(i.sample in 1:1) {
  set.seed(i.sample)

  fimm.indices <- sample.int(n = fimm.n.samples, size = n.samples, replace = FALSE)
  ohsu.indices <- sample.int(n = ohsu.n.samples, size = n.samples, replace = FALSE)
  ctrp.aml.indices <- sample.int(n = ctrp.aml.n.samples, size = n.samples, replace = FALSE)
  ctrp.aml2.indices <- sample.int(n = ctrp.aml2.n.samples, size = n.samples, replace = FALSE)
  ctrp.non.aml.indices <- sample.int(n = ctrp.non.aml.n.samples, size = n.samples, replace = FALSE)
  ctrp.non.aml.heme.indices <- sample.int(n = ctrp.non.aml.heme.n.samples, size = n.samples, replace = FALSE)

  do.downsample <- FALSE
  if(!do.downsample) {
    fimm.indices <- 1:fimm.n.samples
    ohsu.indices <- 1:ohsu.n.samples
    ctrp.aml.indices <- 1:ctrp.aml.n.samples
    ctrp.aml2.indices <- 1:ctrp.aml2.n.samples
    ctrp.non.aml.indices <- 1:ctrp.non.aml.n.samples
    ctrp.non.aml.heme.indices <- 1:ctrp.non.aml.heme.n.samples
  }

  all.expr <- cbind(fimm.expr[common.genes, fimm.indices],
                    ohsu.expr[common.genes, ohsu.indices],
                    ctrp.aml.log.cpm.expr[common.genes, ctrp.aml.indices],
                    ctrp.aml2.log.cpm.expr[common.genes, ctrp.aml2.indices],
                    ctrp.non.aml.log.cpm.expr[common.genes, ctrp.non.aml.indices],
                    ctrp.non.aml.heme.log.cpm.expr[common.genes, ctrp.non.aml.heme.indices])

  data.sets <- c(rep("fimm", n.samples), rep("ohsu", n.samples), rep("ctrp.aml", n.samples), rep("ctrp.aml2", n.samples), rep("ctrp.non.aml", n.samples), rep("ctrp.non.aml.heme", n.samples))

  if(FALSE) {
  all.expr <- cbind(fimm.expr[common.genes, fimm.indices],
                    ohsu.expr[common.genes, ohsu.indices])

  data.sets <- c(rep("fimm", length(fimm.indices)), rep("ohsu", length(ohsu.indices)))
  }
  study.expressed.genes <- llply(unique(data.sets), .fun = function(data.set) {
                                   flag <- data.sets == data.set
                                   avg.log.cpm <- unlist(apply(all.expr[, flag], 1, function(row) mean(row, na.rm=TRUE)))
                                   expressed.genes <- rownames(all.expr)[avg.log.cpm > 0]
                                   expressed.genes
                                 })
  expressed.genes <- Reduce(intersect, study.expressed.genes)

  study.variable.genes <- llply(unique(data.sets), .fun = function(data.set) { 
                                                     flag <- data.sets == data.set
                                                     sd.log.cpm <- unlist(apply(all.expr[, flag], 1, function(row) sd(row, na.rm=TRUE)))
                                                     var.genes <- rownames(all.expr); var.genes <- var.genes[order(sd.log.cpm, decreasing=TRUE)]
                                                     var.genes[1:min(length(var.genes),n.genes)]
                                                   })
  variable.genes <- Reduce(intersect, study.variable.genes)

  fimm.flag <- data.sets == "fimm"
  ohsu.flag <- data.sets == "ohsu"
  ctrp.aml.flag <- data.sets == "ctrp.aml"
  ctrp.aml2.flag <- data.sets == "ctrp.aml2"
  ctrp.non.aml.flag <- data.sets == "ctrp.non.aml"
  ctrp.non.aml.heme.flag <- data.sets == "ctrp.non.aml.heme"

  cat(paste0("Doing downsample iteration ", i.sample, "\n"))

  ## corrs[[length(corrs)+1]] <- do.corrs.of.corrs(all.expr, data.sets, common.genes, ohsu.flag, fimm.flag, ctrp.flag, postfix = paste0("ds-", i.sample))
  ## corrs[[length(corrs)+1]] <- do.corrs.of.corrs(all.expr, data.sets, common.genes, c("ohsu", "fimm", "ctrp"), c(ohsu.flag, fimm.flag, ctrp.flag), postfix = paste0("ds-", i.sample))
  genes <- intersect(common.genes, expressed.genes)
  genes <- intersect(genes, variable.genes)
  genes <- common.genes
  print(length(genes))
##  genes <- intersect(genes, cancer.genes)
  if(FALSE) {
  flags <- list(ohsu.flag, fimm.flag, ctrp.aml.flag, ctrp.aml2.flag, ctrp.non.aml.flag, ctrp.non.aml.heme.flag)
  names <- list("ohsu", "fimm", "ctrp.aml", "ctrp.aml2", "ctrp.non.aml", "ctrp.non.aml.heme")

  all.expr <- cbind(fimm.expr[common.genes, fimm.indices],
                    ohsu.expr[common.genes, ohsu.indices])

  data.sets <- c(rep("fimm", length(fimm.indices)), rep("ohsu", length(ohsu.indices)))

  flags <- list(ohsu.flag, fimm.flag)
  names <- list("ohsu", "fimm")
  }
  corrs[[length(corrs)+1]] <- do.corrs.of.corrs(all.expr, data.sets, genes, names,
                                                flags, postfix = paste0("ds-", i.sample))
}
