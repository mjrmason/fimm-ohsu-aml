## At min concentration, percent inhibition should be zero; percent viability should be 100.
## Define mean and std dev of this peak as mean_min_conc and sd_min_conc.
## Exclude response at min concentration that is outside mean_min_conc +/- 3 * sd_min_conc.
## Or, more generally, outside min_response_min_conc and max_response_min_conc

## At max concentration, max percent inhibition should be 100; min percent viability should be 0.
## But lesser inhibitions/greater viabilities are possible.  Hence, fit a gaussian to these peaks.
## Define mean and std dev of this peak as mean_max_conc and sd_max_conc.
## Exclude response at max concentration that has a percent inhibition greater than mean_max_conc + 3 * sd_max_conc or
##                                       that has a percent viability less than mean_max_conc - 3 * sd_max_conc.

## Exclude _any_ inhibition that is outside mean_min_conc - 3 * sd_min_conc < inhibition < mean_max_conc + 3 * sd_max_conc
## or, equivalently,
##         _any_ viability that is outside mean_max_conc - 3 * sd_max_conc < viability < mean_min_conc + 3 * sd_min_conc

## Fit gaussian to response at min drug concentration.  This is easy, just take mean and std of all response 
## (after excluding outliers).

## NB: this function assumes that the response is inhibition (not survival).  And that it is a percent (e.g., 100 = 100%).
## e.g., below we define a "reasonable range of responses" as those between -100 and 100% at lowest concentration.  i.e., they straddle 0% (inhibition).
## If these responses were survival, we would instead expect a reasonable range to straddle 100%. (survival).
get.responses.at.min.conc <- function(drug.data.long, drug.screen.cols, conc.col, response.col) {
  ## Return the response at the min concentration for each screen
  responses.at.min.conc <- ddply(drug.data.long, .variables = drug.screen.cols,
                                        .fun = function(df) { (df[which(df[, conc.col] == min(df[, conc.col]))[1], ]) })
  responses.at.min.conc
}

get.responses.at.max.conc <- function(drug.data.long, drug.screen.cols, conc.col, response.col) {
  ## Return the response at the max concentration for each screen
  responses.at.max.conc <- ddply(drug.data.long, .variables = drug.screen.cols,
                                        .fun = function(df) { (df[which(df[, conc.col] == max(df[, conc.col]))[1], ]) })
  responses.at.max.conc
}

define.response.at.min.conc.distribution <- function(drug.data.long, drug.screen.cols, conc.col, response.col) {

  responses.at.min.conc <- get.responses.at.min.conc(drug.data.long, drug.screen.cols, conc.col, response.col)

  ## Limit to a reasonable range of responses
  resp.flag <- (responses.at.min.conc[, response.col] > -100) & (responses.at.min.conc[, response.col] < 100)
  vec <- responses.at.min.conc[resp.flag, response.col]

  mean_response_min_conc <- mean(vec)
  sd_response_min_conc <- sd(vec)

  return(list("responses.at.min.conc" = responses.at.min.conc, "mean" = mean_response_min_conc, "sd" = sd_response_min_conc))
}

annotate.drug.response.outliers_ <- function(drug.data.long, drug.screen.cols, conc.col, response.col, min_response_min_conc, max_response_min_conc, max_response, prefix = prefix) {

  responses.at.min.conc <- get.responses.at.min.conc(drug.data.long, drug.screen.cols, conc.col, response.col)

  ## Limit to a reasonable range of responses
  resp.flag <- (responses.at.min.conc[, response.col] > -100) & (responses.at.min.conc[, response.col] < 100)
  vec <- responses.at.min.conc[resp.flag, response.col]

  mean_response_min_conc <- mean(vec)
  sd_response_min_conc <- sd(vec)

  resp.flag <- (responses.at.min.conc[, response.col] <= min_response_min_conc) | (responses.at.min.conc[, response.col] >= max_response_min_conc)
  outlier.responses.at.min.conc <- responses.at.min.conc[resp.flag, ]

  file <- paste0(prefix, "-", response.col, "-qq.pdf")
  pdf(file)
  names(vec) <- 1:length(vec)
  qq <- qqnorm(vec, main = paste0(data.set, " ", response.col, " at min concentration"))
  qqline(vec)
  abline(h = min_response_min_conc)
  abline(h = max_response_min_conc)
  d <- dev.off()

  file <- paste0(prefix, "-", response.col, "-density.pdf")
  pdf(file)
  g <- ggplot(data.frame(responses.at.min.conc = vec))
##  g <- g + geom_density(aes(x = responses.at.min.conc))
  g <- g + geom_histogram(aes(x = responses.at.min.conc), binwidth=5)
  g <- g + xlab(paste0(response.col, " at min concentration"))
  g <- g + ggtitle(paste0(data.set, " ", response.col, " at min concentration"))
  g <- g + geom_vline(xintercept = min_response_min_conc)
  g <- g + geom_vline(xintercept = max_response_min_conc)
  print(g)
  d <- dev.off()

  ## The minimum response at the minimum concentration is also the minimum allowed response at _any_ concentration
  ## This only makes sense for inhibition, not cell survival
  min_response <- min_response_min_conc

  ## Return the response at the max concentration for each screen
  ## This only makes sense for inhibition, not cell survival
  responses.at.max.conc <- get.responses.at.max.conc(drug.data.long, drug.screen.cols, conc.col, response.col)

  if(FALSE) {
    ## Find the peak of the gaussian near the expected max response
    expected.max.response.at.max.conc <- 100
    flag <- (responses.at.max.conc[, response.col] > 90) & (responses.at.max.conc[, response.col] < 110)
    d <- density(responses.at.max.conc[flag, response.col])
    flag <- d$x > 50
    dx <- d$x[flag]
    dy <- d$y[flag]
    x.peak <- dx[dy == max(dy)]

    num.sigmas <- 1
    ## Probability density between mean + num.sigmas * sigma and mean is:
    prob.density <- pnorm(num.sigmas) - pnorm(0)
    ## So go out from x.peak until we have covered prob.density -- that will tell us sigma
    flag <- (responses.at.max.conc[, response.col] >= x.peak) & (responses.at.max.conc[, response.col] < 150)
    ox <- sort(responses.at.max.conc[flag, response.col])
    sigma <- (ox[floor(prob.density*length(ox))] - x.peak)/num.sigmas

    emp <- rnorm(n = 1000, mean = x.peak, sd = sigma)

    ## peak.indx <- which(abs(d$x - expected.max.response.at.max.conc) == min(abs(d$x - expected.max.response.at.max.conc)))
    ## x.peak <- d$x[peak.indx]
    ## Reflect the distribution about the peak
    vals <- responses.at.max.conc[(responses.at.max.conc[, response.col] >= x.peak) & (responses.at.max.conc[, response.col] <= 250), response.col]
    reflected.vals <- x.peak - abs(responses.at.max.conc[(responses.at.max.conc[, response.col] > x.peak) & (responses.at.max.conc[, response.col] <= 250), response.col] - x.peak)
    all.vals <- c(vals, reflected.vals)
    plot(d)
    dnew <- density(all.vals)
    dnew <- density(emp)
    dnew$y <- dnew$y * max(d$y) / max(dnew$y)
    lines(dnew)
  } ## END FALSE

  cat(paste0("min_response_min_conc: ", min_response_min_conc, "\n"))
  cat(paste0("max_response_min_conc: ", max_response_min_conc, "\n"))
  cat(paste0("min_response: ", min_response, "\n"))
  cat(paste0("max_response: ", max_response, "\n"))

  vec <- responses.at.max.conc[, response.col]
  file <- paste0(prefix, "-max-", response.col, "-density.pdf")
  pdf(file)
  g <- ggplot(data.frame(responses.at.max.conc = vec))
##  g <- g + geom_density(aes(x = responses.at.max.conc))
  g <- g + geom_histogram(aes(x = responses.at.max.conc), bindwidth = 5)
  g <- g + geom_vline(xintercept = max_response)
  g <- g + xlab(paste0(response.col, " at max concentration"))
  g <- g + ggtitle(paste0(data.set, " ", response.col, " at max concentration"))
  print(g)
  d <- dev.off()

  resp.flag <- (responses.at.max.conc[, response.col] >= max_response)
  outlier.responses.at.max.conc <- responses.at.max.conc[resp.flag, ]

  ## Finally, exclude any response that is <= min_response or >= max_response
  resp.flag <- (drug.data.long[, response.col] <= min_response) | (drug.data.long[, response.col] >= max_response)
  outlier.responses.at.any.conc <- drug.data.long[resp.flag, ]

  all.outliers <- unique(rbind(outlier.responses.at.min.conc, outlier.responses.at.max.conc, outlier.responses.at.any.conc))

  ## Mark outliers in the original raw data
  all.outliers$outlier.response <- TRUE
  drug.data.long <- merge(drug.data.long, all.outliers, all.x = TRUE)
  flag <- is.na(drug.data.long$outlier.response)
  drug.data.long$outlier.response[flag] <- FALSE
  return(drug.data.long)
}

annotate.drug.response.outliers <- function(drug.data.long, drug.screen.cols, conc.col, response.col, prefix = prefix) {

  responses.at.min.conc <- get.responses.at.min.conc(drug.data.long, drug.screen.cols, conc.col, response.col)

  ## Limit to a reasonable range of responses
  resp.flag <- (responses.at.min.conc[, response.col] > -100) & (responses.at.min.conc[, response.col] < 100)
  vec <- responses.at.min.conc[resp.flag, response.col]

  mean_response_min_conc <- mean(vec)
  sd_response_min_conc <- sd(vec)

  ## Define min/max responses at min concentration based on outliers.
  ## e.g., min_response_min_conc is the _maximum_ response _less than_ the mean that is an outlier
  ## NB: then exclude any response that is less than or equal to this min_response_min_conc
  outliers <- boxplot.stats(vec)$out
  min_response_min_conc <- round(max(outliers[outliers < mean_response_min_conc]), digits = 2)
  max_response_min_conc <- round(min(outliers[outliers > mean_response_min_conc]), digits = 2)

  ## The response is sharply peaked near 100 -- so just set max response to 105.
  ## No. Use the spread near the 0 peak to determine this spread.
  max_response_max_conc <- 105
  max_response_max_conc <- round(100 + 0.5 * ( max_response_min_conc - min_response_min_conc ), digits = 2)
  max_response <- max_response_max_conc

  drug.data.long <- annotate.drug.response.outliers_(drug.data.long, drug.screen.cols, conc.col, response.col, min_response_min_conc, max_response_min_conc, max_response, prefix = prefix)

  return(list("drug.data.long" = drug.data.long, "min_response" = min_response, "max_response" = max_response, 
              "responses.at.min.conc" = responses.at.min.conc, "responses.at.max.conc" = responses.at.max.conc))
}

annotate.drug.response.outlier.pvals_ <- function(drug.data.long, response.col, min.conc.mean, min.conc.sd, max.conc.mean, max.conc.sd) {
  drug.data.long$lower.outlier.pval <- unlist(lapply(drug.data.long[, response.col], function(x) pnorm(x, mean = min.conc.mean, sd = min.conc.sd, lower.tail = TRUE)))
  drug.data.long$upper.outlier.pval <- unlist(lapply(drug.data.long[, response.col], function(x) pnorm(x, mean = max.conc.mean, sd = max.conc.sd, lower.tail = FALSE)))
  drug.data.long$outlier.response <- FALSE
  flag <- ( drug.data.long$lower.outlier.pval < 0.05 ) | ( drug.data.long$upper.outlier.pval < 0.05 )
  drug.data.long$outlier.response[flag] <- TRUE
  return(drug.data.long)

}

find.extremal.responses.qqplots_ <- function(drug.data.long, response.col) {
  vec <- drug.data.long[, response.col]
  d <- density(vec)
  max.indx <- which(d$y == max(d$y))
  reflection.pt <- d$x[max.indx]
  pts.to.reflect <- vec[vec < reflection.pt]
  reflected.pts <- reflection.pt + abs(pts.to.reflect - reflection.pt)
  all.pts <- c(vec[vec <= reflection.pt], reflected.pts)
  outliers <- boxplot.stats(all.pts)$out
  min_response_min_conc <- round(max(outliers[outliers < 0]), digits = 2)
  qqnorm(all.pts)
  qqline(all.pts)
  abline(h = min_response_min_conc)
  min.response <- min_response_min_conc
  
  max.response <- NA
  if(max(vec) < 110) {
    max.response <- max(vec) + 10
  } else {
    stop("Have not implemented find.extremal.responses.qqplots_ for max end\n")
  }
  return(list("min.response" = min.response, "max.response" = max.response))
}

annotate.drug.response.outliers.min.max_ <- function(drug.data.long, response.col, min.response, max.response) {
  drug.data.long$outlier.response <- FALSE
  flag <- is.na(drug.data.long[, response.col]) | ( drug.data.long[, response.col] < min.response ) | ( drug.data.long[, response.col] > max.response )
  drug.data.long$outlier.response[flag] <- TRUE
  return(drug.data.long)
}

