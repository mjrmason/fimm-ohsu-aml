find.overlap.of.sparse.predictors <- function(tbl.input, transform, top) {
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

        m <- merge(tmp1, tmp2, by = merge.col, suffixes = c(".1", ".2"))
        m <- m[, c(merge.col, "n.features.1", "coeffs.1", "coeffs.2", "coeff.vals.1", "coeff.vals.2", "train.set.1", "train.set.2", "gene.set.1")]
        m$id <- 1:nrow(m)
        hyp <- ddply(m, .variables = "id", 
               .fun = function(r) {
                        n.f <- as.numeric(r$n.features.1)
                        coeff.vals.1 <- as.numeric(unlist(strsplit(as.character(r$coeff.vals.1), split=",")))
                        coeffs.1 <- (unlist(strsplit(as.character(r$coeffs.1), split=",")))
                        names(coeff.vals.1) <- coeffs.1
                        coeffs.1 <- coeff.vals.1
                        coeff.vals.2 <- as.numeric(unlist(strsplit(as.character(r$coeff.vals.2), split=",")))
                        coeffs.2 <- (unlist(strsplit(as.character(r$coeffs.2), split=",")))
                        names(coeff.vals.2) <- coeffs.2
                        coeffs.2 <- coeff.vals.2
                        switch(transform,
                          "drop.zero" = { 
                             coeffs.1 <- coeffs.1[coeffs.1 > 0]
                             coeffs.2 <- coeffs.2[coeffs.2 > 0]
                          },
                          "none" = { },
                          "abs" = {
                             coeffs.1 <- abs(coeffs.1)
                             coeffs.2 <- abs(coeffs.2)
                          },
                          "negate" = {
                             coeffs.1 <- - coeffs.1
                             coeffs.2 <- - coeffs.2
                          })
                        if(top != 0) {
                          coeffs.1 <- sort(coeffs.1)
                          coeffs.2 <- sort(coeffs.2)
                          coeffs.1 <- coeffs.1[1:min(top, length(coeffs.1))]
                          coeffs.2 <- coeffs.2[1:min(top, length(coeffs.2))]
                        }
                        n.1 <- length(coeffs.1)
                        n.2 <- length(coeffs.2)
                        both <- intersect(names(coeffs.1), names(coeffs.2))
                        n.both <- length(both)
                        both <- paste(both, collapse=",")
                        ## pval <- phyper(q = n.both - 1, m = n.1, n = n.f - n.1, k = n.2, lower.tail = FALSE) 
                        fet <- fisher.test(cbind(c(n.both, n.1), c(n.2, n.f - n.1 - n.2 + n.both)), alternative = "greater")
                        pval <- fet$p.value
                        vec <- c(r$id[1], n.1, n.2, n.both, both, pval)
                        names(vec) <- c("id", "n.1", "n.2", "n.both", "both", "pval")
                        vec
                      })
        hyp <- merge(hyp, m, by = "id")
        ret.tbl <- rbind(ret.tbl, hyp)
      }
    }
  }
  ret.tbl
}