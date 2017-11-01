suppressPackageStartupMessages(library("plyr"))

## Concentration parameters (ic50, x.min, and x.max are in real, not log, space)
## LL.4 uses llogistic to fit this function:
## f(x) = c + \frac{d-c}{1+\exp(b(\log(x)-\log(e)))}
## i.e., the parameter e is the IC50 (not log IC50)
## and x is in real (not log) space
## b = slope
## c = min.asymptote
## d = max.asymptote
## e = IC50
ll.4.func <- function(x, b, c, d, e) {
  c + ( (d-c)/(1+exp(b*(log(x)-log(e)))) )
}

## L.4 uses logistic to fit this function:
## f(x) = c + \frac{d-c}{1+\exp(b(x-e))}
## i.e., the parameter e is the IC50 
## and x is in real (not log) space
## b = slope
## c = min.asymptote
## d = max.asymptote
## e = IC50
l.4.func <- function(x, b, c, d, e) {
  c + ( (d-c)/(1+exp(b*(x-e))) )
}

## CTRP uses a (log-)logistic-like function.  
## This 4-parameter model is given in the supplement of
## Harnessing connectivity in a large-scale small-molecule sensitivity dataset
## Seashore-Ludlow et al (2015) Cancer Discovery
## Note some differences between l4.func above:
## 1. denominator is h, not the difference of the min and max asymptotes
## 2. slope beta is in the exponent's denominator
## The first commented out function is the one given in the manuscript.
## However, empirically, x is on a log2 scale.  Hence, I have pushed the log2 into the function.
## Hence, the function is actually related to a log-logistic function.
## As such, there is a Jacobian when we integrate and we can not use the
## analytical integral from the logistic function.  This should instead
## be integrated numerically.
## Note that integrating the first function (with x) from log2(min.conc) to log2(max.conc)
## gives the AUC computed by Seashore-Ludlow.  But, this is incorrect, we should integrate
## the function with log2(x) from min.conc to max.conc.  And this gives different results.
## ctrp.curve <- function(x, b, h, alpha, beta) { b + ( h / ( 1 + exp(- ( ( x - alpha ) / beta ) ) ) ) }
ctrp.curve <- function(x, b, h, alpha, beta) { b + ( h / ( 1 + exp(- ( ( log2(x) - alpha ) / beta ) ) ) ) }

## Integrate the CTRP log-logistic-like function numerically.  See above
## comments to understand why we don't evaluate the (apparently) logistic function
## used by Seashore-Ludlow analytically.
## NB: min.conc/max.conc are in uMol (real, not log, space)
compute.ctrp.auc <- function(b, h, alpha, beta, min.conc, max.conc) {
  auc <- NA
  fnc <- function(x) { ctrp.curve(x, b = b, h = h, alpha = alpha, beta = beta) }
  val <- integrate(fnc, min.conc, max.conc)
  if(val$message == "OK") { auc <- val$value }
  auc
}

## Compute the AUC for a _log-logistic_ fit using LL.4/llogistic with drc
## All concentrations/IC50s are in real space
## NB: FIMM/Yadav fits a logistic function
## NB: min.conc/max.conc are in nMol (real, not log, space)
compute.ll.4.auc <- function(b.param, c.param, d.param, e.param, min.conc, max.conc) {
  ## FIMM defines logistic (not log-logistic) function as:
  ## f(x) = d + ( a - d ) * [ 1 + 10^(b * c - b * x) ]^-1   (Yadav Eq. 1)
  ## with (indefinite) integral
  ## Y(x) = a * x + { ( 1 / b ) * (a - d) * log10 [ 1 + 10^(b * c - b * x) ] }
  ## 
  ## LL.4 in DRC is:
  ## f(x) = c' + ( d' - c' ) / ( 1 + exp(b' * log(x) - b' * log(e')) )  (LL.4 Eq. 1)
  ## Input parameters are defined in terms of LL.4.  i.e., b.param = b'

  ## a = d' = max.asymptote
  ## d = c' = min.asymptote
  ## c = e' = ic50
  ## The interpretation of the slope/b/b' is different between the two
  ## models, but not needed to calculate dss.

  ## LL.4 in DRC is:
  ## f(x) = c + \frac{d-c}{1+\exp(b(\log(x)-\log(e)))}
  ## f(x) = c + ( d - c ) / ( 1 + exp(b * log(x) - b * log(f)) )
  ## FIMM/Yadav defines x1 as (Eq. 4), but this is defined for the logistic function (Yadav Eq. 1)
  ## x1 = c - (1/b) * [ log10(a - t) - log10(t - d) ]   (Yadav Eq. 4)
  ## This comes from inverting (Yadav Eq. 1) above.  
  ## Similarly, inverting LL.4 Eq. 1 above gives:
  ## x1 = e' * exp{ (1/b') * [ log(d' - t) - log(t - c') ] }
  ## Begin the integration at x1(y=t)

  fnc <- function(x) { ll.4.func(x, b.param, c.param, d.param, e.param) }

  auc <- NA

  val <- tryCatch({integrate(fnc, min.conc, max.conc)}, error = function(e) { NA })
  if(!(class(val) == "integrate")) { return(NA) }
  if(val$message == "OK") { auc <- val$value }
  if(!is.finite(auc)) {
    auc <- NA
  }
  return(auc)
}

## Compute the AUC for a _logistic_ fit using L.4/logistic with drc
## NB: min.conc/max.conc are in nMol (real, not log, space)
compute.l.4.auc <- function(b.param, c.param, d.param, e.param, min.conc, max.conc, numerical = TRUE) {
  ## FIMM defines logistic (not log-logistic) function as:
  ## f(x) = d + ( a - d ) * [ 1 + 10^(b * c - b * x) ]^-1   (Yadav Eq. 1)
  ##
  ## DRC L.4 is
  ## f(x) = c' + ( d' - c' ) * [ 1 + exp(b' * x - b' * e' ) ]^-1
  ## i.e., 
  ## a = d' = max.asymptote
  ## b = - ln(10) b' = -ln(10) slope
  ## c = e' = ic50
  ## d = c' = min.asymptote
  ## Here, the primed values are those from L.4 and correspond to the parameters of this function.
  ## e.g., b' = b.param
  auc <- NA
  if(numerical) {  
    fnc <- function(x) { l.4.func(x, b.param, c.param, d.param, e.param) }
    val <- integrate(fnc, min.conc, max.conc)
    if(val$message == "OK") { auc <- val$value }
  } else {
    ## Here a, b, c, and d are defined in terms of 
    a <- d.param
    b <- - log(10) * b.param
    c <- e.param
    d <- c.param

    ## Implement Yadav Eq. 3:
    ## Y(x) = a * x + { ( 1 / b ) * (a - d) * log10 [ 1 + 10^(b * c - b * x) ] }  (Yadav Eq. 3)
    y.int <- function(x, a, b, c, d) {
      ## This exponential sometimes blows up
      res <- 0
      exponent <- (b * c - b * x)
      pwr <- 10^exponent
      ## This exponential sometimes blows up
      ## res <- a * x + ( ( 1 / b ) * (a - d) * log10(1 + 10^exponent ) )
      ## If exponent >> 1, then approximate log10(1 + 10^exponent) ~ log10(10^exponent) = exponent
      if((pwr > 0) && is.infinite(pwr)) { 
        res <- a * x + ( ( 1 / b ) * (a - d) * exponent )
      } else {
        res <- a * x + ( ( 1 / b ) * (a - d) * log10(1 + pwr ) )
      }
      res
    }
    auc <- y.int(max.conc, a, b, c, d) - y.int(min.conc, a, b, c, d)
    if(!is.finite(auc)) { auc <- NA }
  }
  auc
}

## Compute DSS1, DSS2, and DSS3 as described by Yadav et al. (2014) Scientific Reports
## t: min activity level (between 0 and 100)
## c.min, c.max: min and max concentrations (real, not log, scales)
## x1, x2: min and max of selected (real, not log) concentration range over which to integrate.
##         Yadav states that x2 = c.max is the default.
##         If x1 is NULL, then start the integration at the point at which 
##         t is crossed.  x1 is defined by Yadav Eqn 4.
## top.asymptote is the max response asymptote from the curve fit.
## Yadav sets DSS = 0 when the ic50 is at or beyond the max dose level tested c.max.
## We won't do that here.
## I'm not sure why Yadav parameterizes this by both x1 and c.min (and x2 and c.max),
## I believe x1 = c.min and x2 = c.max.
compute.dss <- function(auc, t = 10, c.min, c.max, x1, x2, top.asymptote) {
  if(x1 != c.min) {
    stop(paste0("x1 (", x1, ") != c.min (", c.min, ")\n"))
  }
  if(x2 != c.max) {
    stop(paste0("x2 (", x2, ") != c.max (", c.max, ")\n"))
  }
  if(c.max < c.min) { 
    warning("c.max should be > c.min\n") 
    ret <- c(NA, NA, NA)
    names(ret) <- c("dss1", "dss2", "dss3")
    return(ret)
  }
  if(x2 < x1) { 
    warning("x2 should be > x1\n") 
    ret <- c(NA, NA, NA)
    names(ret) <- c("dss1", "dss2", "dss3")
    return(ret)
  }
  ## Distinguish between cases in which it makes sense for DSS to be zero and when
  ## it is an error condition and should be NA.
  if(!is.finite(auc)) {
    ret <- c(NA, NA, NA)
    names(ret) <- c("dss1", "dss2", "dss3")
    return(ret)
  }
  if(auc == 0) {
    ret <- c(0, 0, 0)
    names(ret) <- c("dss1", "dss2", "dss3")
    return(ret)
  }
  if(any(is.na(c(auc, t, c.min, c.max, x1, x2, top.asymptote)))) {
    ret <- c(NA, NA, NA)
    names(ret) <- c("dss1", "dss2", "dss3")
    return(ret)
  }
 
  ## Define DSS1
  dss1 <- auc - t * (x2 - x1)
  dss1 <- dss1 / ( (100 - t) * (c.max - c.min) )
  ## This catches dividing by zero (above)
  if(c.max == c.min) { dss1 <- NA }

  ## Define DSS2
  dss2 <- dss1 / log10(top.asymptote)
  ## This catches taking the log of zero/negative (above).
  if(top.asymptote <= 0) { dss2 <- NA }

  ## Define DSS3
  dss3 <- dss2 * (x2 - x1) / (c.max - c.min)
  ## This catches dividing by zero (above)
  if(c.max == c.min) {
    ## Numerator and denominator cancel here
    if( (x2 - x1) == (c.max - c.min) ) {
      dss3 <- dss2
    } else {
      dss3 <- NA 
    }
  }
  ret <- c(dss1, dss2, dss3)
  names(ret) <- c("dss1", "dss2", "dss3")
  ret
}

## Compute DSS1, DSS2, and DSS3 for the _log-logistic_ fit for CTRP
## Return a vector holding:
## the AUC value used to calculate DSS ("dss.auc"),
## the "shifted" AUC value where the lower asymptote is considered to have response 0 ("shifted.auc"),
## and DSS1, DSS2, and DSS3 as "dss1", "dss2", and "dss3", respectively.
## NB: dss.auc is the DSS calculated between min.conc and max.conc
## NB: min.conc/max.conc are in uMol (real, not log, space)
## Is t in [0, 1] or [0, 100] (generally)
compute.ctrp.dss <- function(b, h, alpha, beta, t = 0, min.conc, max.conc) {
  ## CTRP fits the function (see ctrp.curve above)
  ## f(x) = b + ( h / ( 1 + exp(- ( ( log2(x) - alpha ) / beta ) ) ) ) 
  ## Invert this function for f(x) = t gives:
  ## x1 = 2^alpha * 2^{ beta * [ log(t - b - h) - log(t - b) ] }

  ## Determine the lower limit of integration.
  x1 <- min.conc
  if(t != 0) {
    x1 <- tryCatch({(2^alpha) * 2^( b * (log(t - b - h) - log(t - b)) )},
                   error = function(e) { NA })
  }
  if(!is.finite(x1)) {
    vec <- c(NA, NA, NA, NA, NA)
    names(vec) <- c("dss.auc", "shifted.auc", "dss1", "dss2", "dss3")
    return(vec)
  }
  x1 <- max(x1, min.conc)

  ## Compute AUC to be used for DSS.  It should be calculated between x1 and max.conc.
  dss.auc <- tryCatch({compute.ctrp.auc(b, h, alpha, beta, x1, max.conc)},
                      error = function(e) { NA })
  if(!is.finite(dss.auc)) {
    vec <- c(NA, NA, NA, NA, NA)
    names(vec) <- c("dss.auc", "shifted.auc", "dss1", "dss2", "dss3")
    return(vec)
  }

  ## Subtract off the area/rectangle "below" the lower asymptote (between y = 0 and
  ## the y value of the lower asymptote) and between the x limits of integration
  ## x1 to max.conc
  ## NB: we are taking the lesser of the "lower" and "upper" asymptotes from the
  ## original fit (not at the truncated range of integration)
  y.low <- b
  y.high <- b + h
  y.asym <- min(y.low, y.high)
  y0 <- 0
  shifted.auc <- dss.auc - ( ( y.asym - y0 ) * ( max.conc - x1 ) )

  ## Compute DSS1, 2, and 3.
  dss <- compute.dss(dss.auc, t = t, c.min = x1, c.max = max.conc, x1 = x1, x2 = max.conc, top.asymptote = b + h)
  vec <- c(dss.auc, shifted.auc, dss)
  names(vec) <- c("dss.auc", "shifted.auc", names(dss))
  return(vec)
}

## Compute DSS1, DSS2, and DSS3 for a _log-logistic_ fit using LL.4/llogistic with drc
## Return a vector holding:
## the AUC value used to calculate DSS ("dss.auc")
## the "shifted" AUC value where the lower asymptote is considered to have response 0 ("shifted.auc"),
## and DSS1, DSS2, and DSS3 as "dss1", "dss2", and "dss3", respectively.
## NB: dss.auc is the DSS calculated between min.conc and max.conc
## NB: min.conc/max.conc are in nMol (real, not log, space)
compute.ll.4.dss <- function(b.param, c.param, d.param, e.param, t = 0, min.conc, max.conc) {
  ## FIMM defines logistic (not log-logistic) function as:
  ## f(x) = d + ( a - d ) * [ 1 + 10^(b * c - b * x) ]^-1   (Yadav Eq. 1)
  ## 
  ## LL.4 in DRC is:
  ## f(x) = c' + ( d' - c' ) / ( 1 + exp(b' * log(x) - b' * log(e')) )  (LL.4 Eq. 1)
  ## Input parameters are defined in terms of LL.4.  i.e., b.param = b'

  ## a = d' = max.asymptote
  ## d = c' = min.asymptote
  ## c = e' = ic50
  ## The interpretation of the slop/b/b' is different between the two
  ## models, but not needed to calculate dss.

  ## LL.4 in DRC is:
  ## f(x) = c + \frac{d-c}{1+\exp(b(\log(x)-\log(e)))}
  ## f(x) = c + ( d - c ) / ( 1 + exp(b * log(x) - b * log(f)) )
  ## FIMM/Yadav defines x1 as (Eq. 4), but this is defined for the logistic function:
  ## x1 = c - (1/b) * [ log10(a - t) - log10(t - d) ] 
  ## This comes from inverting (Yadav Eq. 1) above.  
  ## Similarly, inverting LL.4 Eq. 1 above gives:
  ## x1 = e' * exp{ (1/b') * [ log(d' - t) - log(t - c') ] }
  ## Begin the integration at x1(y=t)

  ## Determine the lower limit of integration.
  x1 <- min.conc
  if(t != 0) {
    x1 <- tryCatch({e.param * exp( (1/b.param) * (log(d.param - t) - log(t - c.param)) )},
                   error = function(e) { NA })
  }
  if(!is.finite(x1)) {
    vec <- c(NA, NA, NA, NA, NA)
    names(vec) <- c("dss.auc", "shifted.auc", "dss1", "dss2", "dss3")
    return(vec)
  }
  x1 <- max(x1, min.conc)

  ## Compute AUC to be used for DSS.  It should be calculated between x1 and max.conc.
  dss.auc <- tryCatch({compute.ll.4.auc(b.param, c.param, d.param, e.param, x1, max.conc)},
                      error = function(e) { NA })
  if(!is.finite(dss.auc)) {
    vec <- c(NA, NA, NA, NA, NA)
    names(vec) <- c("dss.auc", "shifted.auc", "dss1", "dss2", "dss3")
    return(vec)
  }

  ## Subtract off the area/rectangle "below" the lower asymptote (between y = 0 and
  ## the y value of the lower asymptote) and between the x limits of integration
  ## x1 to max.conc
  ## NB: we are taking the lesser of the "lower" and "upper" asymptotes from the
  ## original fit (not at the truncated range of integration)
  y.low <- c.param
  y.high <- d.param
  y.asym <- min(y.low, y.high)
  y0 <- 0
  shifted.auc <- dss.auc - ( ( y.asym - y0 ) * ( max.conc - x1 ) )

  ## Compute DSS1, 2, and 3.
  dss <- compute.dss(dss.auc, t = t, c.min = x1, c.max = max.conc, x1 = x1, x2 = max.conc, top.asymptote = d.param)
  vec <- c(dss.auc, shifted.auc, dss)
  names(vec) <- c("dss.auc", "shifted.auc", names(dss))
  return(vec)
}

## Compute DSS1, DSS2, and DSS3 for a _logistic_ fit using L.4/logistic with drc
## Return a vector holding:
## the AUC value used to calculate DSS ("dss.auc")
## the "shifted" AUC value where the lower asymptote is considered to have response 0 ("shifted.auc"),
## and DSS1, DSS2, and DSS3 as "dss1", "dss2", and "dss3", respectively.
## NB: dss.auc is the DSS calculated between min.conc and max.conc
## NB: min.conc/max.conc are in nMol (real, not log, space)
compute.l.4.dss <- function(b.param, c.param, d.param, e.param, t = 0, min.conc, max.conc, numerical = TRUE) {
  ## FIMM defines logistic (not log-logistic) function as:
  ## f(x) = d + ( a - d ) * [ 1 + 10^(b * c - b * x) ]^-1   (Yadav Eq. 1)
  ## 
  ## L.4 in DRC is:
  ## f(x) = c' + ( d' - c' ) / ( 1 + exp(b' * x - b' * e') )  (L.4 Eq. 1)
  ## Input parameters are defined in terms of L.4.  i.e., b.param = b'
  ## a = d' = max.asymptote
  ## d = c' = min.asymptote
  ## c = e' = ic50
  ## b' = - b log(10)

  ## FIMM/Yadav defines x1 as (Eq. 4), which is defined for the logistic function (Yadav Eq. 1)
  ## x1 = c - (1/b) * [ log10(a - t) - log10(t - d) ]   (Yadav Eq. 4)
  ## This comes from inverting (Yadav Eq. 1) above.  
  ## 
  ## Similarly, inverting L.4 Eq. 1 (and/or plugging the above substitutions into Yadav Eq. 4 gives:
  ## x1 = e' + (1/b') * [ log(d' - t) - log(t - c') ]
  ## (i.e., x1 = e' + [ log(10) / b' ] * [ log10(d' - t) - log10(t - c') ], recalling that log10(x) = log(x) / log(10) )
  ## Begin the integration at x1(y=t)

  ## Determine the lower limit of integration.
  x1 <- min.conc
  if(t != 0) {
    x1 <- tryCatch({e.param + (1/b.param) * (log(d.param - t) - log(t - c.param))},
                   error = function(e) { NA })
  }
  if(!is.finite(x1)) {
    vec <- c(NA, NA, NA, NA, NA)
    names(vec) <- c("dss.auc", "shifted.auc", "dss1", "dss2", "dss3")
    return(vec)
  }
  x1 <- max(x1, min.conc)

  ## Compute AUC to be used for DSS.  It should be calculated between x1 and max.conc.
  dss.auc <- tryCatch({compute.l.4.auc(b.param, c.param, d.param, e.param, x1, max.conc, numerical = numerical)},
                      error = function(e) { NA })
  if(!is.finite(dss.auc)) {
    vec <- c(NA, NA, NA, NA, NA)
    names(vec) <- c("dss.auc", "shifted.auc", "dss1", "dss2", "dss3")
    return(vec)
  }

  ## Subtract off the area/rectangle "below" the lower asymptote (between y = 0 and
  ## the y value of the lower asymptote) and between the x limits of integration
  ## x1 to max.conc.
  ## NB: we are taking the lesser of the "lower" and "upper" asymptotes from the
  ## original fit (not at the truncated range of integration)
  y.low <- c.param
  y.high <- d.param
  y.asym <- min(y.low, y.high)
  y0 <- 0
  shifted.auc <- dss.auc - ( ( y.asym - y0 ) * ( max.conc - x1 ) )

  ## Compute DSS1, 2, and 3.
  dss <- compute.dss(dss.auc, t = t, c.min = x1, c.max = max.conc, x1 = x1, x2 = max.conc, top.asymptote = d.param)
  vec <- c(dss.auc, shifted.auc, dss)
  names(vec) <- c("dss.auc", "shifted.auc", names(dss))
  return(vec)
}

compute.all.dss <- function(data, t = 0, fct = "LL.4") {
  res <- ddply(data,
               colnames(data), 
               .parallel = TRUE,
               .fun = function(df) {
                 if(nrow(df) != 1) { stop("Expected to process one row\n") }
                 converged <- as.numeric(df$converged[1])
                 vec <- c(NA, NA, NA, NA, NA)
                 names(vec) <- c("dss.auc", "shifted.auc", "dss1", "dss2", "dss3")
                 if(converged == 1) {
                   b.param <- as.numeric(df$b[1])
                   c.param <- as.numeric(df$c[1])
                   d.param <- as.numeric(df$d[1])
                   e.param <- as.numeric(df$e[1])
                   min.conc <- as.numeric(df$min.conc[1])
                   max.conc <- as.numeric(df$max.conc[1])
                   vec <- tryCatch({
                              if(fct == "LL.4") {
                                compute.ll.4.dss(b.param, c.param, d.param, e.param, t = t, min.conc, max.conc)
                              } else {
                                compute.l.4.dss(b.param, c.param, d.param, e.param, t = t, min.conc, max.conc, numerical = FALSE) 
                              }
                            }, 
                            error = function(e) {
                              vec <- c(NA, NA, NA, NA, NA)
                              names(vec) <- c("dss.auc", "shifted.auc", "dss1", "dss2", "dss3")
                              vec    
                            })
                   }
                   vec
                 })
  res
}

compute.all.ctrp.dss <- function(data, t = 0) {
  res <- ddply(data,
               colnames(data), 
               .parallel = TRUE,
               .fun = function(df) {
                 if(nrow(df) != 1) { stop("Expected to process one row\n") }
                 vec <- c(NA, NA, NA, NA)
                 names(vec) <- c("dss.auc", "dss1", "dss2", "dss3")
                 b <- as.numeric(df$p4_baseline[1])
                 h <- as.numeric(df$p3_total_decline[1])
                 alpha <- as.numeric(df$p1_center[1])
                 beta <- as.numeric(df$p2_slope[1])
                 ## min/max concentration are in nMol, but compute.ctrp.dss expects values in uMol
                 min.conc <- as.numeric(df$min.conc[1]) / 10^3
                 max.conc <- as.numeric(df$max.conc[1]) / 10^3
                 vec <- tryCatch({
                            compute.ctrp.dss(b = b, h = h, alpha = alpha, beta = beta, t = t, min.conc, max.conc)
                          }, 
                          error = function(e) {
                            vec <- c(NA, NA, NA, NA)
                            names(vec) <- c("dss.auc", "dss1", "dss2", "dss3")
                            vec    
                          })
                   vec
                 })
  res
}

