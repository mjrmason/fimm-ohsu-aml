library(ggplot2)
library(ggbeeswarm)

## This from https://stats.stackexchange.com/questions/15011/generate-a-random-variable-with-a-defined-correlation-to-an-existing-variable/15040#15040
## Create a vector x2 that has the given correlation from the vector x1
calculate.correlated.vector <- function(x1, rho) {
  n <- length(x1)
  theta <- acos(rho)             # corresponding angle
  x2    <- rnorm(n, 2, 0.5)      # new random data
  X     <- cbind(x1, x2)         # matrix
  Xctr  <- scale(X, center=TRUE, scale=FALSE)   # centered columns (mean 0)

  Id   <- diag(n)                               # identity matrix
  Q    <- qr.Q(qr(Xctr[ , 1, drop=FALSE]))      # QR-decomposition, just matrix Q
  P    <- tcrossprod(Q)          # = Q Q'       # projection onto space defined by x1
  x2o  <- (Id-P) %*% Xctr[ , 2]                 # x2ctr made orthogonal to x1ctr
  Xc2  <- cbind(Xctr[ , 1], x2o)                # bind to matrix
  Y    <- Xc2 %*% diag(1/sqrt(colSums(Xc2^2)))  # scale columns to length 1

  x <- Y[ , 2] + (1 / tan(theta)) * Y[ , 1]     # final new vector
  ## cor(x1, x)                 
  x
}

## Generate 100 500 vectors and 500 targeted correlations (around a mean correlation)
vec.len <- 500
num.samples <- 100
mean.corr <- 0.3
corrs <- rnorm(num.samples, mean = mean.corr, sd = 0.1 * mean.corr)
x1s <- lapply(1:num.samples, function(i) rnorm(n = vec.len))
x2s <- lapply(1:num.samples, function(i) calculate.correlated.vector(x1s[[i]], corrs[[i]]))

vec.sizes <- c(50, 100, 250, 500)

corrs.vs.sample.size <- ldply(vec.sizes,
                              .fun = function(vec.size) {
                                       indices <- lapply(1:num.samples, function(i) sample.int(vec.len, size = vec.size, replace = FALSE))
                                       x1.ds <- lapply(1:num.samples, function(i) x1s[[i]][indices[[i]]])
                                       x2.ds <- lapply(1:num.samples, function(i) x2s[[i]][indices[[i]]])
                                       cs <- unlist(lapply(1:num.samples, function(i) cor(x1.ds[[i]], x2.ds[[i]])))
                                       df <- data.frame(cor = cs)
                                       df$vec.size <- vec.size
                                       df
                              })

corrs.vs.sample.size$vec.size <- factor(corrs.vs.sample.size$vec.size)
g <- ggplot(data = corrs.vs.sample.size, aes(x = vec.size, y = cor))
g <- g + geom_violin()
## g <- g + geom_beeswarm()
g <- g + geom_boxplot(width = 0.1)
g <- g + xlab("Vector Length")
g <- g + ylab("Spearman Correlation")
g