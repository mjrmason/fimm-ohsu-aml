
## for(mdl in unique(tbl$model)) {
for(mdl in c("LASSO", "Ridge")) {
  for(resp in unique(renamed.tbl$response)) {
    for(gset in c("gene")) {
      for(met in c("pearson", "spearman")) {
        train.ds <- "ohsu"; test.ds <- "fimm"
        drugs <- rownames(all.resps[[train.ds]])
        names(drugs) <- drugs
        train.drug.vs.mean.resp <- ldply(drugs, .fun = function(drug) cor(all.resps[[train.ds]][drug,,drop=TRUE], mean.resps[[train.ds]], method=met, use="pairwise.complete.obs"))
        colnames(train.drug.vs.mean.resp) <- c("drug", "train.mean.cor")
        test.drug.vs.mean.resp <- ldply(drugs, .fun = function(drug) cor(all.resps[[test.ds]][drug,,drop=TRUE], mean.resps[[test.ds]], method=met, use="pairwise.complete.obs"))
        colnames(test.drug.vs.mean.resp) <- c("drug", "test.mean.cor")
        for(indx in 1:length(pt.responses1)) {
          pt.response1 <- pt.responses1[indx]; pt.response2 <- pt.responses2[indx]
          sub1 <- subset(renamed.tbl, metric == met & response == resp & gene.set == gset & train.set == "ohsu" & test.set == "fimm" & model == mdl & pt.response == pt.response1)
          sub2 <- subset(renamed.tbl, metric == met & response == resp & gene.set == gset & train.set == "ohsu" & test.set == "fimm" & model == mdl & pt.response == pt.response2)
          sub <- merge(sub1, sub2, by = "train.drug", suffixes = c(".1", ".2"))
          sub <- merge(sub, train.drug.vs.mean.resp, by.x = "train.drug", by.y = "drug")
          sub <- merge(sub, test.drug.vs.mean.resp, by.x = "train.drug", by.y = "drug")
          sub$diff <- sub$val.1 - sub$val.2
          ## Plot difference in pred/actual correlation (test data) vs drug/mean correlation (train or test data)
          for(ds in c("train", "test")) {  
            flag <- is.na(sub$val.1) | (sub$val.1 == 0)
            flag <- flag | is.na(sub$val.2) | (sub$val.2 == 0)
            g <- ggplot(data = sub[!flag,])
            g <- g + geom_point(aes_string(x = "diff", y = paste0(ds, ".mean.cor"))) 
            g <- g + xlab(paste(pt.response1, "-", pt.response2, "pred/actual", capwords(met), "correlation", sep=" "))
            g <- g + ylab(paste("Drug/mean response", capwords(met), "correlation (", capwords(ds), "data", ")", sep=" "))
            g <- g + ggtitle(paste0(resp, " ", capwords(met), " ", mdl, " ", gset))
            pdf(paste(file.prefix, mdl, resp, met, gset, make.names(pt.response1), "vs", make.names(pt.response2), "volcano.pdf", sep="-"))
            print(g)
            d <- dev.off()
          }
##          flag <- !is.na(sub$val.1) & !is.na(sub$val.2) & ( abs(sub$val.1 - sub$val.2) > 0.1 )
          flag <- !is.na(sub$val.1) & !is.na(sub$val.2) & ( (sub$val.1 - sub$val.2) < -0.1 )
          flag <- flag | ( sub$train.drug == "Venetoclax" )
          sub$color <- "black"
          sub$color[flag] <- "blue"
          g <- ggplot(data = sub)
          g <- g + geom_abline(slope = 1, intercept = 0, linetype = "dashed")
          if(any(flag)) {
            labels <- sub[flag,]
            extremal.drugs <- merge(labels, drug.map[, c("OHSU_DRUG_NAME", "FIMM_DRUG_NAME")], by.x = "train.drug", by.y = "FIMM_DRUG_NAME")
            if(nrow(labels) != nrow(extremal.drugs)) {
              stop("Unexpected change in extremal drug table")
            }
            write.table(file = paste(file.prefix, mdl, resp, met, gset, make.names(pt.response1), "vs", make.names(pt.response2), "extremal-points.tsv", sep="-"), extremal.drugs, sep="\t", row.names=FALSE, col.names=TRUE)
            nudge <- 0.025
set.seed(1234)
##            g <- g + geom_text(data = labels, aes(x = val.1, y = val.2, label = train.drug), position = position_jitter(width=nudge, height=nudge), hjust = "right", size = 8)
##            g <- g + geom_text(data = labels, aes(x = val.1, y = val.2, label = train.drug), nudge_x = - 0.02, hjust = "right", size = 6)
              left.labels <- labels
              right.labels <- NULL
              if((mdl == "Ridge") && (met == "pearson")) {
                right.labels <- subset(labels, train.drug %in% c("XAV-939"))
              }
              if(!is.null(right.labels)) {
                left.labels <- subset(left.labels, !(train.drug %in% right.labels$train.drug))
              }
              if(!is.null(left.labels) && (nrow(left.labels) >= 1)) { g <- g + geom_text(data = left.labels, aes(x = val.1, y = val.2, label = train.drug), nudge_x = - 0.02, hjust = "right", size = 6) }
              if(!is.null(right.labels) && (nrow(right.labels) >= 1)) { g <- g + geom_text(data = right.labels, aes(x = val.1, y = val.2, label = train.drug), nudge_x = 0.02, hjust = "left", size = 6) }
          }
          g <- g + geom_point(data = sub, aes(x = val.1, y = val.2, col = color), show.legend = FALSE)
          lb <- min(min(c(sub$val.1, sub$val.2)) - 0.1, -0.2)
          g <- g + xlim(c(lb,1)) + ylim(c(lb,1))
          g <- g + xlab(capwords(paste(pt.response1, met, "correlation", sep=" ")))
          g <- g + ylab(capwords(paste(pt.response2, met, "correlation", sep=" ")))
          g <- g + theme(text = element_text(size = 20))
##          g <- g + ggtitle(paste0(resp, " ", capwords(met), " ", mdl, " ", gset))
          pdf(paste(file.prefix, mdl, resp, met, gset, make.names(pt.response1), "vs", make.names(pt.response2), "comparison.pdf", sep="-"))
          print(g)
          d <- dev.off()
          ## Plot the correlation of each drug vs mean response and group according to relation between the two pt.responses:
          ## (1) cov(pred_mean.feature.only, actual = 0); (2) cor(pred_1, actual) - cor(pred_2, actual) > 0.1; (3) cor(pred_1, actual) - cor(pred_2, actual) < -0.1; other
          sub$relation <- unlist(apply(sub[, c("val.1", "val.2")], 1, 
                                 function(row) {
                                   x <- row[1]; y <- row[2]; relation <- "error"
                                   if( x == 0 ) { relation <- paste0(pt.response1, " = 0")
                                   } else if (y == 0) { relation <- paste0(pt.response2, " = 0")
                                   } else if ( ( x - y ) > 0.1 ) { relation <- paste0(pt.response1, " >>")
                                   } else if ( ( x - y ) < -0.1 ) { relation <- paste0(pt.response2, " >>")
                                   } else { relation <- "similar" }
                                   relation
                                 }))
          sub$coarse.relation <- unlist(lapply(sub$relation, function(str) ifelse(grepl(str, pattern="= 0"), str, "non-degenerate model")))
          sub$relation <- factor(sub$relation, levels = c(paste0(pt.response1, " = 0"), paste0(pt.response2, " = 0"), "similar", paste0(pt.response1, " >>"), paste0(pt.response2, " >>")))
          sub$coarse.relation <- factor(sub$coarse.relation, levels = c(paste0(pt.response1, " = 0"), paste0(pt.response2, " = 0"), "non-degenerate model"))
          g <- ggplot(data = sub[!grepl(sub$relation, pattern="= 0"),], aes(x = relation, y = test.mean.cor))
          g <- g + geom_violin()
          g <- g + geom_beeswarm()
          g <- g + xlab(paste0("Relation between ", pt.response1, " and ", pt.response2, " correlations (Test: ", capwords(test.ds), ")"))
          g <- g + ylab(paste0("Drug correlation with mean response (Test: ", capwords(test.ds), ")")) 
          g <- g + ggtitle(paste0(resp, " ", capwords(met), " ", mdl, " ", gset))
          pdf(paste(file.prefix, mdl, resp, met, gset, make.names(pt.response1), "vs", make.names(pt.response2), "test-mean-response-correlation.pdf", sep="-"))
          print(g)
          d <- dev.off()
          g <- ggplot(data = sub, aes(x = coarse.relation, y = train.mean.cor))
          g <- g + geom_violin()
          g <- g + geom_beeswarm()
          g <- g + xlab(paste0("Relation between ", pt.response1, " and ", pt.response2, " correlations (Test: ", capwords(test.ds), ")"))
          g <- g + ylab(paste0("Drug correlation with mean response (Training: ", capwords(train.ds), ")")) 
          g <- g + ggtitle(paste0(resp, " ", capwords(met), " ", mdl, " ", gset))
          pdf(paste(file.prefix, mdl, resp, met, gset, make.names(pt.response1), "vs", make.names(pt.response2), "train-mean-response-correlation.pdf", sep="-"))
          print(g)
          d <- dev.off()
        }
      }
    }
  } }

