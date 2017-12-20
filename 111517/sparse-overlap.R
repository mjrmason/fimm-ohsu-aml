## transform = c("none", "abs", "negate")
## top = # of features to consider in overall
## top = 0 -> use all
find.overlap.of.sparse.predictors <- function(tbl.input, universe, transform, top, gene.sets = NULL, verbose = FALSE) {
  ret.tbl <- c()
  for(gs in unique(tbl.input$gene.set)) {
    tmp <- subset(tbl.input, gene.set == gs)
    tr.sets <- as.character(unique(tmp$train.set))
    for(i in 1:(length(tr.sets)-1)) {
      tr1 <- tr.sets[i]
      tmp1 <- subset(tmp, train.set == tr1)
      merge.col <- drug.name.col
      tmp1 <- merge(drug.name.tbl, tmp1, by.x = merge.col, by.y = "train.drug")
      for(j in (i+1):length(tr.sets)) {
        tr2 <- tr.sets[j]
        if(((tr1 == "ohsu") || (tr2 == "ohsu")) && (grepl(pattern="ohsu.set1", paste0(tr1,tr2)) || grepl(pattern="ohsu.set2", paste0(tr1,tr2)))) { next }
        if(((tr1 == "fimm") || (tr2 == "fimm")) && (grepl(pattern="fimm.set1", paste0(tr1,tr2)) || grepl(pattern="fimm.set2", paste0(tr1,tr2)))) { next }
        tmp2 <- subset(tmp, train.set == tr2)
        merge.col <- drug.name.col
        tmp2 <- merge(drug.name.tbl, tmp2, by.x = merge.col, by.y = "train.drug")
        tmp1$id.1 <- 1:nrow(tmp1)
        m <- merge(tmp1, tmp2, by = merge.col, suffixes = c(".1", ".2"))
        m <- m[, c(merge.col, "id.1", "n.features.1", "n.features.2", "coeffs.1", "coeffs.2", "coeff.vals.1", "coeff.vals.2", "train.set.1", "train.set.2", "gene.set.1")]
        m$id <- 1:nrow(m)
        hyp <- ddply(m, .variables = "id.1", .parallel = FALSE,
                     .fun = function(r1) {
                        n.f1 <- as.numeric(r1$n.features.1)
                        coeff.vals.1 <- as.numeric(unlist(strsplit(as.character(r1$coeff.vals.1), split=",")))
                        coeffs.1 <- (unlist(strsplit(as.character(r1$coeffs.1), split=",")))
                        names(coeff.vals.1) <- coeffs.1
                        coeffs.1 <- coeff.vals.1
                        switch(transform,
                          "drop.zero" = { 
                             coeffs.1 <- coeffs.1[coeffs.1 > 0]
                          },
                          "none" = { },
                          "abs" = {
                             coeffs.1 <- abs(coeffs.1)
                          },
                          "abs.drop.zero" = {
                             coeffs.1 <- abs(coeffs.1)
                             coeffs.1 <- coeffs.1[coeffs.1 != 0]
                          },
                          "negate" = {
                             coeffs.1 <- - coeffs.1
                          })

                        ddply(r1, .variables = "id", .parallel = FALSE,
                              .fun = function(r) {
                                 coeff.vals.2 <- as.numeric(unlist(strsplit(as.character(r$coeff.vals.2), split=",")))
                                 coeffs.2 <- (unlist(strsplit(as.character(r$coeffs.2), split=",")))
                                 names(coeff.vals.2) <- coeffs.2
                                 coeffs.2 <- coeff.vals.2
                                 switch(transform,
                                   "drop.zero" = { 
                                      coeffs.2 <- coeffs.2[coeffs.2 > 0]
                                   },
                                   "none" = { },
                                   "abs" = {
                                      coeffs.2 <- abs(coeffs.2)
                                   },
                                   "abs.drop.zero" = {
                                      coeffs.2 <- abs(coeffs.2)
                                      coeffs.2 <- coeffs.2[coeffs.2 != 0]
                                   },
                                   "negate" = {
                                      coeffs.2 <- - coeffs.2
                                   })

                                 ldply(top, .parallel = FALSE,
                                       .fun = function(sz) {
                                                c1 <- coeffs.1
                                                c2 <- coeffs.2
                                                if((sz != 0) && (length(c1) > 0)) {
                                                  c1 <- sort(c1, decreasing=TRUE)
                                                  c1 <- c1[1:min(sz, length(c1))]
                                                }
                                                if(length(c1) > 0) {
                                                  c1 <- c1[!duplicated(names(c1))]
                                                }
                                                
                                                if((sz != 0) && (length(c2) > 0)) {
                                                  c2 <- sort(c2, decreasing=TRUE)
                                                  c2 <- c2[1:min(sz, length(c2))]
                                                }
                                                if(length(c2) > 0) {
                                                  c2 <- c2[!duplicated(names(c2))]
                                                }
                                                ## Find overlap between genes and gene sets--the gene sets will be our features.
                                                if(!is.null(gene.sets)) {
                        
                                                }
                                                n.1 <- length(c1)
                                                n.2 <- length(c2)
                                                both <- intersect(names(c1), names(c2))
                                                n.both <- length(both)
                                                both <- paste(both, collapse=",")
                                                ## pval <- phyper(q = n.both - 1, m = n.1, n = n.f1 - n.1, k = n.2, lower.tail = FALSE) 
n.1 <- length(setdiff(names(c1),names(c2)))
n.2 <- length(setdiff(names(c2),names(c1)))
neither <- length(setdiff(universe, union(names(c1),names(c2))))
mat <- cbind(c(n.both, n.1), c(n.2, neither))
##                                                mat <- cbind(c(n.both, n.1 - n.both), c(n.2 - n.both, n.f1 - n.1))
                                                fet <- fisher.test(mat, alternative = "greater")
                                                pval <- fet$p.value
if(verbose) {
  print(mat)
  print(fet)
  print(pval)
}
                                                vec <- c(r$id[1], n.1, n.2, n.both, both, pval, sz)
                                                names(vec) <- c("id", "n.1", "n.2", "n.both", "both", "pval", "top")
                                                vec
                                  })
                               })
                             })
        hyp <- hyp[, !grepl(pattern="id.1", x=colnames(hyp))]
        m <- m[, !grepl(pattern="id.1", x=colnames(m))]
        hyp <- merge(hyp, m, by = "id")
        ret.tbl <- rbind(ret.tbl, hyp)
      }
    }
  }
  ret.tbl
}
