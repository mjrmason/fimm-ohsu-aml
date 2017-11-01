source("models.R")

fit.ComBat <-
function (dat, batch, mod = NULL, par.prior = TRUE, prior.plots = FALSE, 
    mean.only = FALSE) 
{
    if (mean.only == TRUE) {
        cat("Using the 'mean only' version of ComBat\n")
    }
    if (length(dim(batch)) > 1) {
        stop("This version of ComBat only allows one batch variable")
    }
    batch <- as.factor(batch)
    batchmod <- model.matrix(~-1 + batch)
    cat("Found", nlevels(batch), "batches\n")
    n.batch <- nlevels(batch)
    batches <- list()
    for (i in 1:n.batch) {
        batches[[i]] <- which(batch == levels(batch)[i])
    }
    n.batches <- sapply(batches, length)
    if (any(n.batches == 1)) {
        mean.only = TRUE
        cat("Note: one batch has only one sample, setting mean.only=TRUE\n")
    }
    n.array <- sum(n.batches)
    design <- cbind(batchmod, mod)
    check <- apply(design, 2, function(x) all(x == 1))
    design <- as.matrix(design[, !check])
    cat("Adjusting for", ncol(design) - ncol(batchmod), "covariate(s) or covariate level(s)\n")
    if (qr(design)$rank < ncol(design)) {
        if (ncol(design) == (n.batch + 1)) {
            stop("The covariate is confounded with batch! Remove the covariate and rerun ComBat")
        }
        if (ncol(design) > (n.batch + 1)) {
            if ((qr(design[, -c(1:n.batch)])$rank < ncol(design[, 
                -c(1:n.batch)]))) {
                stop("The covariates are confounded! Please remove one or more of the covariates so the design is not confounded")
            }
            else {
                stop("At least one covariate is confounded with batch! Please remove confounded covariates and rerun ComBat")
            }
        }
    }
    NAs = any(is.na(dat))
    if (NAs) {
        cat(c("Found", sum(is.na(dat)), "Missing Data Values\n"), 
            sep = " ")
    }
    cat("Standardizing Data across genes\n")
    if (!NAs) {
        B.hat <- solve(t(design) %*% design) %*% t(design) %*% 
            t(as.matrix(dat))
    }
    else {
        B.hat = apply(dat, 1, Beta.NA, design)
    }
    grand.mean <- t(n.batches/n.array) %*% B.hat[1:n.batch, ]
    if (!NAs) {
        var.pooled <- ((dat - t(design %*% B.hat))^2) %*% rep(1/n.array, 
            n.array)
    }
    else {
        var.pooled <- apply(dat - t(design %*% B.hat), 1, var, 
            na.rm = T)
    }
    stand.mean <- t(grand.mean) %*% t(rep(1, n.array))
    if (!is.null(design)) {
        tmp <- design
        tmp[, c(1:n.batch)] <- 0
        stand.mean <- stand.mean + t(tmp %*% B.hat)
    }
    s.data <- (dat - stand.mean)/(sqrt(var.pooled) %*% t(rep(1, 
        n.array)))
    cat("Fitting L/S model and finding priors\n")
    batch.design <- design[, 1:n.batch]
    if (!NAs) {
        gamma.hat <- solve(t(batch.design) %*% batch.design) %*% 
            t(batch.design) %*% t(as.matrix(s.data))
    }
    else {
        gamma.hat = apply(s.data, 1, Beta.NA, batch.design)
    }
    delta.hat <- NULL
    for (i in batches) {
        if (mean.only == TRUE) {
            delta.hat <- rbind(delta.hat, rep(1, nrow(s.data)))
        }
        else {
            delta.hat <- rbind(delta.hat, apply(s.data[, i], 
                1, var, na.rm = T))
        }
    }
    gamma.bar <- apply(gamma.hat, 1, mean)
    t2 <- apply(gamma.hat, 1, var)
    a.prior <- apply(delta.hat, 1, sva:::aprior)
    b.prior <- apply(delta.hat, 1, sva:::bprior)
    if (prior.plots & par.prior) {
        par(mfrow = c(2, 2))
        tmp <- density(gamma.hat[1, ])
        plot(tmp, type = "l", main = "Density Plot")
        xx <- seq(min(tmp$x), max(tmp$x), length = 100)
        lines(xx, dnorm(xx, gamma.bar[1], sqrt(t2[1])), col = 2)
        qqnorm(gamma.hat[1, ])
        qqline(gamma.hat[1, ], col = 2)
        tmp <- density(delta.hat[1, ])
        invgam <- 1/rgamma(ncol(delta.hat), a.prior[1], b.prior[1])
        tmp1 <- density(invgam)
        plot(tmp, typ = "l", main = "Density Plot", ylim = c(0, 
            max(tmp$y, tmp1$y)))
        lines(tmp1, col = 2)
        qqplot(delta.hat[1, ], invgam, xlab = "Sample Quantiles", 
            ylab = "Theoretical Quantiles")
        lines(c(0, max(invgam)), c(0, max(invgam)), col = 2)
        title("Q-Q Plot")
    }
    gamma.star <- delta.star <- NULL
    if (par.prior) {
        cat("Finding parametric adjustments\n")
        for (i in 1:n.batch) {
            if (mean.only) {
                gamma.star <- rbind(gamma.star, postmean(gamma.hat[i, 
                  ], gamma.bar[i], 1, 1, t2[i]))
                delta.star <- rbind(delta.star, rep(1, nrow(s.data)))
            }
            else {
                temp <- sva:::it.sol(s.data[, batches[[i]]], gamma.hat[i, 
                  ], delta.hat[i, ], gamma.bar[i], t2[i], a.prior[i], 
                  b.prior[i])
                gamma.star <- rbind(gamma.star, temp[1, ])
                delta.star <- rbind(delta.star, temp[2, ])
            }
        }
    }
    else {
        cat("Finding nonparametric adjustments\n")
        for (i in 1:n.batch) {
            if (mean.only) {
                delta.hat[i, ] = 1
            }
            temp <- int.eprior(as.matrix(s.data[, batches[[i]]]), 
                gamma.hat[i, ], delta.hat[i, ])
            gamma.star <- rbind(gamma.star, temp[1, ])
            delta.star <- rbind(delta.star, temp[2, ])
        }
    }
    cat("Adjusting the Data\n")
    bayesdata <- s.data
    j <- 1
    for (i in batches) {
        bayesdata[, i] <- (bayesdata[, i] - t(batch.design[i, 
            ] %*% gamma.star))/(sqrt(delta.star[j, ]) %*% t(rep(1, 
            n.batches[j])))
        j <- j + 1
    }
    bayesdata <- (bayesdata * (sqrt(var.pooled) %*% t(rep(1, 
        n.array)))) + stand.mean
##    return(bayesdata)
    lst <- list("bayesdata" = bayesdata, "var.pooled" = var.pooled, "stand.mean" = stand.mean, "gamma.star" = gamma.star, "delta.star" = delta.star, "B.hat" = B.hat)
    return(lst)
}

apply.ComBat3 <-
function (dat, batch, mod = NULL, B.hat, stand.mean, var.pooled, gamma.star, delta.star, par.prior = TRUE, prior.plots = FALSE, 
    mean.only = FALSE) 
{
    if (mean.only == TRUE) {
        cat("Using the 'mean only' version of ComBat\n")
    }
    if (length(dim(batch)) > 1) {
        stop("This version of ComBat only allows one batch variable")
    }
    batch <- as.factor(batch)
    batchmod <- model.matrix(~-1 + batch)
    cat("Found", nlevels(batch), "batches\n")
    n.batch <- nlevels(batch)
    batches <- list()
    for (i in 1:n.batch) {
        batches[[i]] <- which(batch == levels(batch)[i])
    }
    n.batches <- sapply(batches, length)
    if (any(n.batches == 1)) {
        mean.only = TRUE
        cat("Note: one batch has only one sample, setting mean.only=TRUE\n")
    }
    n.array <- sum(n.batches)
    design <- cbind(batchmod, mod)
    check <- apply(design, 2, function(x) all(x == 1))
    design <- as.matrix(design[, !check])
    s.data <- (dat - stand.mean)/(sqrt(var.pooled) %*% t(rep(1, 
        n.array)))
    batch.design <- design[, 1:n.batch]
    cat("Adjusting the Data\n")
    bayesdata <- s.data
    j <- 1
    for (i in batches) {
        bayesdata[, i] <- (bayesdata[, i] - t(batch.design[i, 
            ] %*% gamma.star))/(sqrt(delta.star[j, ]) %*% t(rep(1, 
            n.batches[j])))
        j <- j + 1
    }
    bayesdata <- (bayesdata * (sqrt(var.pooled) %*% t(rep(1, 
        n.array)))) + stand.mean
    return(bayesdata)
}

apply.ComBat2 <-
function (dat, batch, mod = NULL, gamma.star, delta.star, par.prior = TRUE, prior.plots = FALSE, 
    mean.only = FALSE) 
{
    if (mean.only == TRUE) {
        cat("Using the 'mean only' version of ComBat\n")
    }
    if (length(dim(batch)) > 1) {
        stop("This version of ComBat only allows one batch variable")
    }
    batch <- as.factor(batch)
    batchmod <- model.matrix(~-1 + batch)
    cat("Found", nlevels(batch), "batches\n")
    n.batch <- nlevels(batch)
    batches <- list()
    for (i in 1:n.batch) {
        batches[[i]] <- which(batch == levels(batch)[i])
    }
    n.batches <- sapply(batches, length)
    if (any(n.batches == 1)) {
        mean.only = TRUE
        cat("Note: one batch has only one sample, setting mean.only=TRUE\n")
    }
    n.array <- sum(n.batches)
    design <- cbind(batchmod, mod)
    check <- apply(design, 2, function(x) all(x == 1))
    design <- as.matrix(design[, !check])
    cat("Adjusting for", ncol(design) - ncol(batchmod), "covariate(s) or covariate level(s)\n")
    if (qr(design)$rank < ncol(design)) {
        if (ncol(design) == (n.batch + 1)) {
            stop("The covariate is confounded with batch! Remove the covariate and rerun ComBat")
        }
        if (ncol(design) > (n.batch + 1)) {
            if ((qr(design[, -c(1:n.batch)])$rank < ncol(design[, 
                -c(1:n.batch)]))) {
                stop("The covariates are confounded! Please remove one or more of the covariates so the design is not confounded")
            }
            else {
                stop("At least one covariate is confounded with batch! Please remove confounded covariates and rerun ComBat")
            }
        }
    }
    NAs = any(is.na(dat))
    if (NAs) {
        cat(c("Found", sum(is.na(dat)), "Missing Data Values\n"), 
            sep = " ")
    }
    cat("Standardizing Data across genes\n")
    if (!NAs) {
        B.hat <- solve(t(design) %*% design) %*% t(design) %*% 
            t(as.matrix(dat))
    }
    else {
        B.hat = apply(dat, 1, Beta.NA, design)
    }
    grand.mean <- t(n.batches/n.array) %*% B.hat[1:n.batch, ]
    if (!NAs) {
        var.pooled <- ((dat - t(design %*% B.hat))^2) %*% rep(1/n.array, 
            n.array)
    }
    else {
        var.pooled <- apply(dat - t(design %*% B.hat), 1, var, 
            na.rm = T)
    }
    stand.mean <- t(grand.mean) %*% t(rep(1, n.array))
    if (!is.null(design)) {
        tmp <- design
        tmp[, c(1:n.batch)] <- 0
        stand.mean <- stand.mean + t(tmp %*% B.hat)
    }
    s.data <- (dat - stand.mean)/(sqrt(var.pooled) %*% t(rep(1, 
        n.array)))
    batch.design <- design[, 1:n.batch]
    cat("Adjusting the Data\n")
    bayesdata <- s.data
    j <- 1
    for (i in batches) {
        bayesdata[, i] <- (bayesdata[, i] - t(batch.design[i, 
            ] %*% gamma.star))/(sqrt(delta.star[j, ]) %*% t(rep(1, 
            n.batches[j])))
        j <- j + 1
    }
    bayesdata <- (bayesdata * (sqrt(var.pooled) %*% t(rep(1, 
        n.array)))) + stand.mean
    return(bayesdata)
}

apply.ComBat <-
function (dat, batch, mod = NULL, B.hat, gamma.star, delta.star, par.prior = TRUE, prior.plots = FALSE, 
    mean.only = FALSE) 
{
    if (mean.only == TRUE) {
        cat("Using the 'mean only' version of ComBat\n")
    }
    if (length(dim(batch)) > 1) {
        stop("This version of ComBat only allows one batch variable")
    }
    batch <- as.factor(batch)
    batchmod <- model.matrix(~-1 + batch)
    cat("Found", nlevels(batch), "batches\n")
    n.batch <- nlevels(batch)
    batches <- list()
    for (i in 1:n.batch) {
        batches[[i]] <- which(batch == levels(batch)[i])
    }
    n.batches <- sapply(batches, length)
    if (any(n.batches == 1)) {
        mean.only = TRUE
        cat("Note: one batch has only one sample, setting mean.only=TRUE\n")
    }
    n.array <- sum(n.batches)
    design <- cbind(batchmod, mod)
    check <- apply(design, 2, function(x) all(x == 1))
    design <- as.matrix(design[, !check])
    cat("Adjusting for", ncol(design) - ncol(batchmod), "covariate(s) or covariate level(s)\n")
    if (qr(design)$rank < ncol(design)) {
        if (ncol(design) == (n.batch + 1)) {
            stop("The covariate is confounded with batch! Remove the covariate and rerun ComBat")
        }
        if (ncol(design) > (n.batch + 1)) {
            if ((qr(design[, -c(1:n.batch)])$rank < ncol(design[, 
                -c(1:n.batch)]))) {
                stop("The covariates are confounded! Please remove one or more of the covariates so the design is not confounded")
            }
            else {
                stop("At least one covariate is confounded with batch! Please remove confounded covariates and rerun ComBat")
            }
        }
    }
    NAs = any(is.na(dat))
    if (NAs) {
        cat(c("Found", sum(is.na(dat)), "Missing Data Values\n"), 
            sep = " ")
    }
    cat("Standardizing Data across genes\n")
    grand.mean <- t(n.batches/n.array) %*% B.hat[1:n.batch, ]
    if (!NAs) {
        var.pooled <- ((dat - t(design %*% B.hat))^2) %*% rep(1/n.array, 
            n.array)
    }
    else {
        var.pooled <- apply(dat - t(design %*% B.hat), 1, var, 
            na.rm = T)
    }
    stand.mean <- t(grand.mean) %*% t(rep(1, n.array))
    if (!is.null(design)) {
        tmp <- design
        tmp[, c(1:n.batch)] <- 0
        stand.mean <- stand.mean + t(tmp %*% B.hat)
    }
    s.data <- (dat - stand.mean)/(sqrt(var.pooled) %*% t(rep(1, 
        n.array)))
    cat("Adjusting the Data\n")
    bayesdata <- s.data
    batch.design <- design[, 1:n.batch]
    j <- 1
    for (i in batches) {
        bayesdata[, i] <- (bayesdata[, i] - t(batch.design[i, 
            ] %*% gamma.star))/(sqrt(delta.star[j, ]) %*% t(rep(1, 
            n.batches[j])))
        j <- j + 1
    }
    bayesdata <- (bayesdata * (sqrt(var.pooled) %*% t(rep(1, 
        n.array)))) + stand.mean
    return(bayesdata)
}

do.combat <- function(all.expr, batch, postfix) {

##  batch           <- as.factor(data.sets)
  modcombat       <- model.matrix(~1,data=batch)

  all.expr.combat <- ComBat(all.expr, batch = batch, mod = modcombat)

  all.pca <- prcomp(t(all.expr))
  all.combat.pca <- prcomp(t(all.expr.combat))

  list("all.expr.combat" = all.expr.combat, "pca" = all.pca, "combat.pca" = all.combat.pca)
}


varplot = function(pca, Main=NULL, cols = "slategrey")
{
  pcaSum = summary(pca)$importance; 
  plot(pca$x[,1], pca$x[,2], xlab = paste("PC1", signif(pcaSum[2,1],3)),ylab = paste("PC2", signif(pcaSum[2,2],3)), pch=16, col =cols,cex = .65, main=Main)
}

varplot2 = function(pca, Main=NULL, cols = "slategrey")
{
  pcaSum = summary(pca)$importance; 
  plot(pca$x[,1], pca$x[,3], xlab = paste("PC1", signif(pcaSum[2,1],3)),ylab = paste("PC3", signif(pcaSum[2,3],3)), pch=16, col =cols,cex = .65, main=Main)
}

plot.correction <- function() {

  tmp <- all.expr.combat
  offset <- min(tmp)
  if(offset < 0) {
    tmp <- tmp - offset + 0.001
  }
  mu <- unlist(apply(tmp, 1, mean, na.rm=TRUE))
  std.dev <- unlist(apply(tmp, 1, sd, na.rm=TRUE))
  png(paste0("mu-vs-sd-combat-", postfix, ".png"))
  smoothScatter(mu, std.dev)
  dev.off()

  cov <- std.dev / mu
  names(cov) <- common.genes
  cov <- cov[order(cov, decreasing=TRUE)]

  Cols = as.numeric(as.factor(data.sets))

  png(paste0("pca-", postfix, ".png"))
  par(mfrow=c(2,1))
  varplot(all.pca,cols=Cols, Main = "Uncorrected")
  varplot(all.combat.pca,cols=Cols, Main = "Study Corrected")
  dev.off()


}


do.corrs.of.corrs.old <- function(all.expr, data.sets, common.genes, names, flags, postfix) {

  print(table(data.sets))

  max.y <- 0
  max.x <- 0
  for(i in 1:length(names)) {
    dn <- density(unlist(all.expr[common.genes, flags[[i]]]))
    max.y <- max(max.y, max(dn$y))
    max.x <- max(max.x, max(dn$x))
  }

  png(paste0("exprDensity-", postfix, ".png"))
  plot(density(unlist(all.expr[common.genes, flags[[1]]])), main = "Expression Density", xlab="Expression", ylim=c(0, 1.1*max.y))
  labels <- names
  n.labels <- length(labels)
  legend(0.6*max.x,0.9*max.y,legend=labels, fill=1:n.labels,border=1:n.labels, box.col="white")
  for(i in 2:length(names)){lines(density(unlist(all.expr[common.genes, flags[[i]]])), col = i);}
  dev.off()

  gene.subset <- common.genes

  use.combat <- FALSE
  use.combat <- FALSE
  all.expr.combat <- NULL
  if(use.combat) { 

    batch           <- as.factor(data.sets)
    modcombat       <- model.matrix(~1,data=batch)

    all.expr.combat  <- ComBat(all.expr, batch = batch, mod=modcombat)
    offset <- min(all.expr.combat)
    if(offset < 0) {
      all.expr.combat <- all.expr.combat - offset + 0.001
    }
    mu <- unlist(apply(all.expr.combat[common.genes, ], 1, mean, na.rm=TRUE))
    std.dev <- unlist(apply(all.expr.combat[common.genes, ], 1, sd, na.rm=TRUE))
    png(paste0("mu-vs-sd-combat-", postfix, ".png"))
    smoothScatter(mu, std.dev)
    dev.off()

    cov <- std.dev / mu
    names(cov) <- common.genes
    cov <- cov[order(cov, decreasing=TRUE)]
##    gene.subset <- names(cov)[1:min(1300,length(cov))]

    all.pca <- prcomp(t(all.expr))
    all.combat.pca <- prcomp(t(all.expr.combat))

    Cols = as.numeric(as.factor(data.sets))

    png(paste0("pca-", postfix, ".png"))
    par(mfrow=c(2,1))
    varplot(all.pca,cols=Cols, Main = "Uncorrected")
    varplot(all.combat.pca,cols=Cols, Main = "Study Corrected")
    dev.off()

  }

  cat(paste0(length(gene.subset), " genes in subset.\n"))

  expr.mats <- list()
##  if(use.combat) { 
##    for(i in 1:length(names)) { expr.mats[[i]] <- all.expr.combat[gene.subset, flags[[i]]] }
##  } else {
    for(i in 1:length(names)) { expr.mats[[i]] <- all.expr[gene.subset, flags[[i]]] }
##  }
  names(expr.mats) <- names
  ## expr.mats <- list("ohsu" = all.expr.combat[gene.subset, ohsu.flag], "fimm" = all.expr.combat[gene.subset, fimm.flag], "ctrp" = all.expr.combat[gene.subset, ctrp.flag])

  sampleCors <- list()
  indices <- 1:length(expr.mats)
  names(indices) <- names(expr.mats)
  sampleCors = llply(.data=expr.mats, .fun = function(x){cor(apply(x,2,rank))}, .parallel = T)

  names(sampleCors) <- names(expr.mats)

  if(use.combat) { 
    png(paste0("exprDensity-combat-", postfix, ".png"))
    max.y <- 0
    max.x <- 0
    for(i in 1:length(expr.mats)) {
      dn <- density(unlist(expr.mats[[i]]))
      max.y <- max(max.y, max(dn$y))
      max.x <- max(max.x, max(dn$x))
    }
    plot(density(unlist(expr.mats[[1]])), main = "Expression Density", xlab="Expression", ylim=c(0, 1.1*max.y))
    labels <- names(expr.mats)
    n.labels <- length(labels)
    legend(0.6*max.x,0.9*max.y,legend=labels, fill=1:n.labels,border=1:n.labels, box.col="white")
    for(i in 2:length(expr.mats)){lines(density(unlist(expr.mats[[i]])), col = i);}
    dev.off()
  }

  png(paste0("studyCors-", postfix, ".png"))
  use.fisher <- TRUE
  use.dist <- TRUE
  max.y <- 0
  max.x <- 0
  for(i in 1:length(studyCors)) {
    dn <- NULL
    if(use.fisher) { 
      if(use.dist) {
        dn <- density(fisherz(unlist(as.dist(studyCors[[i]]))), na.rm=TRUE)
      } else {
        dn <- density(fisherz(unlist((studyCors[[i]]))), na.rm=TRUE)
      }
    } else { 
      if(use.dist) { 
        dn <- density((unlist(as.dist(studyCors[[i]]))), na.rm=TRUE)
      } else {
        dn <- density((unlist((studyCors[[i]]))), na.rm=TRUE)
      }
    }
    max.y <- max(max.y, max(dn$y))
    max.x <- max(max.x, max(dn$x))
  }

  if(use.fisher) { 
    if(use.dist) { 
      plot(density(fisherz(unlist(as.dist(studyCors[[1]]))),na.rm=T), main = "Intra Study Correlations", xlab="Fisher z transform of\nPearson Correlations", ylim=c(0, 1.1*max.y))
    } else {
      plot(density(fisherz(unlist((studyCors[[1]]))),na.rm=T), main = "Intra Study Correlations", xlab="Fisher z transform of\nPearson Correlations", ylim=c(0, 1.1*max.y))
    }
  } else { 
    if(use.dist) { 
      plot(density((unlist(as.dist(studyCors[[1]]))),na.rm=T), main = "Intra Study Correlations", xlab="Pearson Correlations", ylim=c(0, 1.1*max.y))
    } else {
      plot(density((unlist((studyCors[[1]]))),na.rm=T), main = "Intra Study Correlations", xlab="Pearson Correlations", ylim=c(0, 1.1*max.y))
    }
  }
  labels <- names(studyCors)
  n.labels <- length(labels)
  legend(0.6*max.x,0.9*max.y,legend=labels, fill=1:n.labels,border=1:n.labels, box.col="white")
  if(use.fisher) { 
    if(use.dist) { 
      for(i in 2:length(studyCors)){lines(density(fisherz(unlist(as.dist(studyCors[[i]]))),na.rm=T), col = i);}
    } else {
      for(i in 2:length(studyCors)){lines(density(fisherz(unlist((studyCors[[i]]))),na.rm=T), col = i);}
    }
  } else { 
    if(use.dist) { 
      for(i in 2:length(studyCors)){lines(density((unlist(as.dist(studyCors[[i]]))),na.rm=T), col = i);}
    } else {
      for(i in 2:length(studyCors)){lines(density((unlist((studyCors[[i]]))),na.rm=T), col = i);}
    }
  }
  abline(v = 0)
  dev.off()

  corMat     <- matrix(NA, choose(length(gene.subset),2), length(expr.mats))

  use.fast <- TRUE
  if(!use.fast) { 
    indices <- 1:length(expr.mats)
    names(indices) <- names(expr.mats)
    geneCors <- llply(indices, .parallel = TRUE,
                       .fun = function(i) {
                         cor(t(expr.mats[[i]]), t(expr.mats[[i]]), use = "pairwise")
                       })
    for(i in 1:length(expr.mats)) {
      corMat[,i] <- as.vector(as.dist(geneCors[[i]]))
    }
  } else { 
    for(i in 1:length(expr.mats))
    {
##      geneCors[[i]]  <- fastCor(expr.mats[[i]], Split = 1500,Method = "spearman")
##      geneCors[[i]] <- fastCor(expr.mats[[i]], Split = 1500,Method = "pearson")
      corMat[,i] <- as.vector(as.dist(fastCor(expr.mats[[i]], Split = 1500,Method = "spearman")))
    }
  }

  # Use Fisher's z transform to remove hard (-1,1) bounds on correlation distributions
  corMat <- apply(corMat,2, fisherz)

  corMat[corMat == Inf]  <- NA
  cors = cor(corMat, use = "pairwise")
  rm(corMat);gc()
  colnames(cors) <- labels
  rownames(cors) <- labels

  print(postfix)
  print(cors)
  ## write.table(data.frame(cors), sep="\t", col.names=T, row.names=T, quote=F)
  cat("\n")

  write.table(cors, file=paste0("corOfCors-", postfix, ".txt"),sep="\t", col.names=T, row.names=T,quote=F)
  list("cors" = cors, "all.expr.combat" = all.expr.combat)
}


do.corrs.of.corrs <- function(all.expr, data.sets, common.genes, names, flags, postfix) {

  print(table(data.sets))

  max.y <- 0
  max.x <- 0
  for(i in 1:length(names)) {
    dn <- density(unlist(all.expr[common.genes, flags[[i]]]))
    max.y <- max(max.y, max(dn$y))
    max.x <- max(max.x, max(dn$x))
  }

  png(paste0("exprDensity-", postfix, ".png"))
  plot(density(unlist(all.expr[common.genes, flags[[1]]])), main = "Expression Density", xlab="Expression", ylim=c(0, 1.1*max.y))
  labels <- names
  n.labels <- length(labels)
  legend(0.6*max.x,0.9*max.y,legend=labels, fill=1:n.labels,border=1:n.labels, box.col="white")
  for(i in 2:length(names)){lines(density(unlist(all.expr[common.genes, flags[[i]]])), col = i);}
  dev.off()

  gene.subset <- common.genes

  cat(paste0(length(gene.subset), " genes in subset.\n"))

  expr.mats <- list()
  for(i in 1:length(names)) { expr.mats[[i]] <- all.expr[gene.subset, flags[[i]]] }
  names(expr.mats) <- names

  sampleCors <- list()
  indices <- 1:length(expr.mats)
  names(indices) <- names(expr.mats)
  sampleCors = llply(.data=expr.mats, .fun = function(x){cor(apply(x,2,rank))}, .parallel = T)

  names(sampleCors) <- names(expr.mats)

  png(paste0("sampleCors-", postfix, ".png"))
  use.fisher <- TRUE
  use.dist <- TRUE
  max.y <- 0
  max.x <- 0
  for(i in 1:length(sampleCors)) {
    dn <- NULL
    if(use.fisher) { 
      if(use.dist) {
        dn <- density(fisherz(unlist(as.dist(sampleCors[[i]]))), na.rm=TRUE)
      } else {
        dn <- density(fisherz(unlist((sampleCors[[i]]))), na.rm=TRUE)
      }
    } else { 
      if(use.dist) { 
        dn <- density((unlist(as.dist(sampleCors[[i]]))), na.rm=TRUE)
      } else {
        dn <- density((unlist((sampleCors[[i]]))), na.rm=TRUE)
      }
    }
    max.y <- max(max.y, max(dn$y))
    max.x <- max(max.x, max(dn$x))
  }

  if(use.fisher) { 
    if(use.dist) { 
      plot(density(fisherz(unlist(as.dist(sampleCors[[1]]))),na.rm=T), main = "Intra Study Correlations", xlab="Fisher z transform of\nPearson Correlations", ylim=c(0, 1.1*max.y))
    } else {
      plot(density(fisherz(unlist((sampleCors[[1]]))),na.rm=T), main = "Intra Study Correlations", xlab="Fisher z transform of\nPearson Correlations", ylim=c(0, 1.1*max.y))
    }
  } else { 
    if(use.dist) { 
      plot(density((unlist(as.dist(sampleCors[[1]]))),na.rm=T), main = "Intra Study Correlations", xlab="Pearson Correlations", ylim=c(0, 1.1*max.y))
    } else {
      plot(density((unlist((sampleCors[[1]]))),na.rm=T), main = "Intra Study Correlations", xlab="Pearson Correlations", ylim=c(0, 1.1*max.y))
    }
  }
  labels <- names(sampleCors)
  n.labels <- length(labels)
  legend(0.6*max.x,0.9*max.y,legend=labels, fill=1:n.labels,border=1:n.labels, box.col="white")
  if(use.fisher) { 
    if(use.dist) { 
      for(i in 2:length(sampleCors)){lines(density(fisherz(unlist(as.dist(sampleCors[[i]]))),na.rm=T), col = i);}
    } else {
      for(i in 2:length(sampleCors)){lines(density(fisherz(unlist((sampleCors[[i]]))),na.rm=T), col = i);}
    }
  } else { 
    if(use.dist) { 
      for(i in 2:length(sampleCors)){lines(density((unlist(as.dist(sampleCors[[i]]))),na.rm=T), col = i);}
    } else {
      for(i in 2:length(sampleCors)){lines(density((unlist((sampleCors[[i]]))),na.rm=T), col = i);}
    }
  }
  abline(v = 0)
  dev.off()

  corMat     <- matrix(NA, choose(length(gene.subset),2), length(expr.mats))

  use.fast <- TRUE
  if(!use.fast) { 
    indices <- 1:length(expr.mats)
    names(indices) <- names(expr.mats)
    geneCors <- llply(indices, .parallel = TRUE,
                       .fun = function(i) {
                         cor(t(expr.mats[[i]]), t(expr.mats[[i]]), use = "pairwise")
                       })
    for(i in 1:length(expr.mats)) {
      corMat[,i] <- as.vector(as.dist(geneCors[[i]]))
    }
  } else { 
    for(i in 1:length(expr.mats))
    {
##      geneCors[[i]]  <- fastCor(expr.mats[[i]], Split = 1500,Method = "spearman")
##      geneCors[[i]] <- fastCor(expr.mats[[i]], Split = 1500,Method = "pearson")
      corMat[,i] <- as.vector(as.dist(fastCor(expr.mats[[i]], Split = 1500,Method = "spearman")))
    }
  }

  # Use Fisher's z transform to remove hard (-1,1) bounds on correlation distributions
  corMat <- apply(corMat,2, fisherz)

  corMat[corMat == Inf]  <- NA
  cors = cor(corMat, use = "pairwise")
  rm(corMat);gc()
  colnames(cors) <- labels
  rownames(cors) <- labels

  print(postfix)
  print(cors)
  ## write.table(data.frame(cors), sep="\t", col.names=T, row.names=T, quote=F)
  cat("\n")

  write.table(cors, file=paste0("corOfCors-", postfix, ".txt"),sep="\t", col.names=T, row.names=T,quote=F)
  file=paste0("corOfCors-", postfix, ".png")
  png(file)
  plot.table(round(cors, digits=3), main="Correlation of\nCorrelations", suppress.rows = FALSE)
  cors
}

