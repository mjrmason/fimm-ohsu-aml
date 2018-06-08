
decimalplaces <- function(x) {
    if ((x %% 1) != 0) {
##        nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed=TRUE)[[1]][[2]])
        nchar(strsplit(as.character(x), ".", fixed=TRUE)[[1]][[2]])
    } else {
        return(0)
    }
}

## Calculate drug ranges for each data set and each fit type (L.4 or LL.4)
## Drug ranges may differ between fits since we do different filtering (e.g., based on GOF)
determine.drug.ranges <- function(fits, screen.id.cols, data.set.drug.id.cols, drug.map, drug.map.drug.id.cols, shared.conc.postfix = ".shared.conc", min.conc.col = "min.conc.nM", max.conc.col = "max.conc.nM", min.conc.output.col = "min.conc", max.conc.output.col = "max.conc") {

  data.sets <- names(fits)
  names(data.sets) <- data.sets
  drug.ranges <- list()
  fcts <- c("LL.4", "L.4")
  exclude.patterns <- list("LL.4" = "\\.l4", "L.4" = "\\.ll4")
  include.patterns <- list("LL.4" = ".ll4", "L.4" = ".l4")
  fct.postfixes <- list("LL.4" = ".ll4", "L.4" = ".l4")

  for(ds in names(fits)) {
    drug.ranges[[ds]] <- list()
    for(fct in fcts) {
      ## Calculate the concentration range for each screen
      ## Do not consider filtered fits
      exclude.cols <- colnames(fits[[ds]])[grepl(pattern=paste0("exclude", fct.postfixes[[fct]]), x=colnames(fits[[ds]]))]
      any.excluded <- unlist(apply(fits[[ds]][, exclude.cols], 1, function(row) any(row)))
      screen.ranges <- ddply(fits[[ds]][!any.excluded, ], .variables = screen.id.cols[[ds]],
                                  .fun = function(df) {
                                           v <- c(min(df[, min.conc.col]), max(df[, max.conc.col]), length(unique(unlist(strsplit(df$uniq.concs.nM[1], split=",")))))  
                                           names(v) <- c(min.conc.col, max.conc.col, "num.concs")
                                           v 
                                  })

      options(scipen=999) ## disable scientific notation so we can compute number of decimals
      ## Round off insignificant differences between concentrations for a given drug.
      ## For any _min_ concentrations that are the same (to within rounding error), set to the max of the rounding-equivalent concentrations
      ## For any _max_ concentrations that are the same (to within rounding error), set to the min of the rounding-equivalent concentrations
      screen.ranges <- ddply(screen.ranges, .variables = c(data.set.drug.id.cols[[ds]]),
                                  .fun = function(df) {
                                           df[, min.conc.col] <- as.numeric(df[, min.conc.col])
                                           df[, max.conc.col] <- as.numeric(df[, max.conc.col])
                                           for(j in 1:nrow(df)) {
                                             num.decimal <- decimalplaces(df[j, min.conc.col])
                                             flag <- round(df[j, min.conc.col], digits = num.decimal) == round(df[, min.conc.col], digits = num.decimal)
                                             df[flag, min.conc.col] <- max(df[flag, min.conc.col])
                                             num.decimal <- decimalplaces(df[j, max.conc.col])
                                             flag <- round(df[j, max.conc.col], digits = num.decimal) == round(df[, max.conc.col], digits = num.decimal)
                                             df[flag, max.conc.col] <- min(df[flag, max.conc.col])
                                           } 
                                           df        
                                  })
      ## Calculate the concentration range for each drug
      drug.ranges[[ds]][[fct]] <- ddply(screen.ranges, .variables = c(data.set.drug.id.cols[[ds]], min.conc.col, max.conc.col),
                                .fun = function(df) {
                                         vec <- c(df[1, data.set.drug.id.cols[[ds]]], df$num.concs[1], df[1, min.conc.col], df[1, max.conc.col], freq = nrow(df))  
                                         names(vec) <- c(data.set.drug.id.cols[[ds]], "num.concs", min.conc.col, max.conc.col, "freq")
                                         vec
                                       })
    }
  }

  ## Define a single range for each drug (in each data set) and each fit type (LL.4 or L.4)
  ## If any range occurs across more than 80% of the screens, use that range.
  ## Otherwise, drop the 20% most narrow ranges and use the most restrictive remaining range.                                 
  ## Define range by the number of concentration points (not the absolute range)
  final.screen.ranges <- list()
  min.frac.for.range <- 0.8
  frac.to.discard <- 0.2
  for(ds in names(fits)) {
    final.screen.ranges[[ds]] <- list()
    for(fct in fcts) {
      final.screen.ranges[[ds]][[fct]] <- ddply(drug.ranges[[ds]][[fct]], .variables = c(data.set.drug.id.cols[[ds]]),
                                  .fun = function(df) {
                                           df[, min.conc.col] <- as.numeric(df[, min.conc.col])
                                           df[, max.conc.col] <- as.numeric(df[, max.conc.col])
                                           df$freq <- as.numeric(df$freq)
                                           df$rel.freq <- df$freq / sum(df$freq)
                                           flag <- ( df$rel.freq == max(df$rel.freq) ) & ( df$rel.freq >= min.frac.for.range )
                                           ret <- NULL
                                           if(any(flag)) {
                                             indx <- which(flag)[1]
                                             ret <- df[indx, c(data.set.drug.id.cols[[ds]], min.conc.col, max.conc.col)]
                                           } else {
                                             ## Discard the smallest frac.to.discard ranges and then define the restrictive range over what's left
                                             ##  df$range <- abs(df[, max.conc.col] - df[, min.conc.col])
                                             df$range <- as.numeric(df$num.concs)
                                             df <- df[order(df$range, df$freq, decreasing=FALSE),]
                                             cum.freq <- 0
                                             first.indx.to.keep <- 1
                                             for(j in 1:(nrow(df)-1)) {
                                               if( ( cum.freq + df$rel.freq[j] ) > frac.to.discard ) { break }
                                               cum.freq <- cum.freq + df$rel.freq[j]
                                               first.indx.to.keep <- j + 1
                                             } 
                                             df <- df[first.indx.to.keep:nrow(df),,drop=F]
                                             ret <- df[1, data.set.drug.id.cols[[ds]], drop=F]
                                             ret[, min.conc.col] <- max(df[, min.conc.col])
                                             ret[, max.conc.col] <- min(df[, max.conc.col])
                                           }
                                         return(ret)
                                  })
    }
  }

  ## Merge the ranges for all data sets together
  shared.drug.ranges <- list()
  for(fct in fcts) {
    shared.drug.ranges[[fct]] <- drug.map
    for(ds in names(fits)) {
      tmp <- final.screen.ranges[[ds]][[fct]][, c(data.set.drug.id.cols[[ds]], min.conc.col, max.conc.col)]
      colnames(tmp) <- c(data.set.drug.id.cols[[ds]], paste0(data.sets[[ds]], ".", min.conc.col), paste0(data.sets[[ds]], ".", max.conc.col))
      cat(paste0("Merging by ", data.set.drug.id.cols[[ds]], " and ", drug.map.drug.id.cols[[ds]], "\n"))
      shared.drug.ranges[[fct]] <- merge(shared.drug.ranges[[fct]], tmp, by.y = data.set.drug.id.cols[[ds]], by.x = drug.map.drug.id.cols[[ds]], all = FALSE)
    }
  }

  min.conc.cols <- unlist(lapply(data.sets, function(ds) paste0(ds, ".", min.conc.col)))
  max.conc.cols <- unlist(lapply(data.sets, function(ds) paste0(ds, ".", max.conc.col)))
  ## DSS code will use min.conc and max.conc columns
  for(fct in fcts) {
    shared.drug.ranges[[fct]][, min.conc.output.col] <- unlist(apply(shared.drug.ranges[[fct]][, min.conc.cols, drop = FALSE], 1, function(row) max(as.numeric(row))))
    shared.drug.ranges[[fct]][, max.conc.output.col] <- unlist(apply(shared.drug.ranges[[fct]][, max.conc.cols, drop = FALSE], 1, function(row) min(as.numeric(row))))
  }

  ## Merge the processed drug response tables with the range tables and exclude any screens that do 
  ## not completely cover the min/max range
  fits.with.concs <- list()
  for(ds in names(fits)) {
    fits.with.concs[[ds]] <- list()
    for(fct in fcts) {
      tmp <- shared.drug.ranges[[fct]][, c(drug.map.drug.id.cols[[ds]], min.conc.output.col, max.conc.output.col)]
      colnames(tmp)[colnames(tmp) == min.conc.output.col] <- paste0("min", shared.conc.postfix, fct.postfixes[[fct]])
      colnames(tmp)[colnames(tmp) == max.conc.output.col] <- paste0("max", shared.conc.postfix, fct.postfixes[[fct]])
      ## ## Exclude min.conc.col and max.conc.col from the fits, to avoid confusion
      ## fits.with.concs[[ds]][[fct]] <- fits[[ds]][, !(colnames(fits[[ds]]) %in% c(min.conc.col, max.conc.col))]
      fits.with.concs[[ds]][[fct]] <- fits[[ds]]
      cat(paste0("Merging by ", drug.map.drug.id.cols[[ds]], " and ", data.set.drug.id.cols[[ds]], "\n"))
      fits.with.concs[[ds]][[fct]] <- merge(fits.with.concs[[ds]][[fct]], tmp, by.y = drug.map.drug.id.cols[[ds]], by.x = data.set.drug.id.cols[[ds]], all = FALSE)
      keep <- unlist(apply(fits.with.concs[[ds]][[fct]][, c("uniq.concs.nM", paste0("min", shared.conc.postfix, fct.postfixes[[fct]]), paste0("max", shared.conc.postfix, fct.postfixes[[fct]]))], 1,   
                           function(row) {
                             concs <- as.numeric(unlist(strsplit(row[1], split=",")))
                             min.conc.nM <- as.numeric(row[2])
                             max.conc.nM <- as.numeric(row[3])
                             (min.conc.nM < max.conc.nM) & (min(concs) < max(concs)) & (min(concs) <= min.conc.nM) & (max(concs) >= max.conc.nM)
                           }))
      fits.with.concs[[ds]][[fct]] <- fits.with.concs[[ds]][[fct]][keep, ]
    }
  }
  return(list("fits.with.concs" = fits.with.concs, "shared.drug.ranges" = shared.drug.ranges))
}

