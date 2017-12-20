## see 070617/fit-ohsu-and-fimm-drug-response-data.R 

## Filter FIMM results to include only those drugs with at least 5 uniq concentration points
for(ds in names(fits)) {
  switch(ds,
    "fimm" = {
      min.conc.pts <- 5
      keep <- sapply(fits[[ds]][, c("uniq.concs.nM")], function(str) length(unique(unlist(strsplit(str, split=",")))) >= min.conc.pts)
      cat(paste0("Keeping FIMM screens with at least ", min.conc.pts, " unique concentrations\n"))
      print(table(keep))
      fits[[ds]] <- fits[[ds]][keep,]
    }, 
    "ohsu" = {
      min.conc.pts <- 5
      max.conc.pts <- 7
      keep <- unlist(apply(fits[[ds]][ ,c("all.concs.nM", "uniq.concs.nM")], 1, function(row) {
                 ( length(unique(unlist(strsplit(row[1], split=",")))) <= max.conc.pts ) && ( length(unique(unlist(strsplit(row[2], split=",
")))) >= min.conc.pts )
              }))
      cat(paste0("Keeping OHSU screens with at least ", min.conc.pts, " unique concentrations and at most ", max.conc.pts, " total concentra
tions\n"))
      print(table(keep))
      fits[[ds]] <- fits[[ds]][keep,]
    },
    {
    })
}

## Plot densities of responses as a function of patient, drug, plate, etc.
groups <- c("DRUG_NAME", "DRUG_SET", "SCREEN_ID", "PATIENT_ID")
for(group in groups) {
    res <- ddply(.data = drug.data.long, .variables = group, .parallel = FALSE,
                 .fun = function(df) {
                          flag <- !is.na(df[, response.col])
                          df <- df[flag, ]
                          vals <- as.numeric(df[, response.col])
                          tot <- nrow(df)
                          num.below <- length(which(vals < min.cutoff))
                          num.above <- length(which(vals > max.cutoff))
                          frac.below <- num.below / tot
                          frac.above <- num.above / tot
                          spread <- max(df[, response.col], na.rm = TRUE) - min(df[, response.col], na.rm = TRUE)
                          data.frame(group = df[1, group], frac.below = frac.below, frac.above = frac.above, spread = spread)
                 })
    df <- drug.data.long
##    df <- df[order(df[, response.col]), ]
##    levels <- unique(df[, group])
    res <- res[order(res$spread),]
    levels <- unique(res[, group])
    print(head(df))
    df[, group] <- factor(df[, group], levels = levels)
    g <- ggplot(data = df)
    g <- g + geom_violin(aes_string(x = group, y = response.col))
    print(g)
    g <- ggplot(data = res)
    g <- g + geom_point(aes_string(x = group, y = "frac.below"))
##    g <- g + geom_bar(aes(x = below.lb))
##    g <- g + facet_wrap(as.formula(paste("~", group)))
    print(g)
}

