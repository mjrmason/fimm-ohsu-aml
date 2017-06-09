bioc.packages <- c("nplr", "openxlsx", "drc", "minpack.lm", "maftools", "vcd", "binom", "DMwR", "WGCNA", "VennDiagram")

install.bioc.packages.auto <- function(x) { 
  x <- as.character(x)
  print(x)
##  if(isTRUE(x %in% .packages(all.available=TRUE))) { 
##    eval(parse(text = sprintf("require(\"%s\")", x)))
##  } else { 
##    #update.packages(ask= FALSE) #update installed packages.
##    eval(parse(text = sprintf("install.packages(\"%s\", dependencies = TRUE)", x)))
##  }

  if(isTRUE(x %in% .packages(all.available=TRUE))) { 
    eval(parse(text = sprintf("require(\"%s\")", x)))
  } else {
    source("http://bioconductor.org/biocLite.R")
    #biocLite(character(), ask=FALSE) #update installed packages.
    eval(parse(text = sprintf("biocLite(\"%s\")", x)))
    eval(parse(text = sprintf("require(\"%s\")", x)))
  }
}

d <- unlist(lapply(bioc.packages, install.bioc.packages.auto))

suppressPackageStartupMessages(library("nplr"))
suppressPackageStartupMessages(library("foreach"))
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(tidyr))

suppressPackageStartupMessages(library("drc"))
suppressPackageStartupMessages(library("minpack.lm"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("gridExtra"))
suppressPackageStartupMessages(library("caTools"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("maftools"))
suppressPackageStartupMessages(library("pcaMethods"))
## suppressPackageStartupMessages(library("GGally"))
suppressPackageStartupMessages(library("VennDiagram"))

library(ggbeeswarm)
library(plyr)
library(stringr)
library(synapseClient)
synapseLogin()

## The following comes from Deducer (https://github.com/ifellows/Deducer/blob/master/Deducer/R/util.R)
## which I could not install.
## Renamed this function from d
ddf<-function(..., row.names = NULL, check.rows = FALSE,
            check.names = FALSE,
            stringsAsFactors = FALSE){
  data.frame(...,row.names=row.names,check.rows=check.rows,check.names=check.names,stringsAsFactors=stringsAsFactors)
}

cor.matrix<-function(variables,with.variables,data=NULL,test=cor.test,...){
  arguments <- as.list(match.call()[-1])
  variables<-eval(substitute(variables),data,parent.frame())
  if(length(dim(variables))<1.5){
    variables<-ddf(variables)
    fn <- arguments$variables
    names(variables)<-if(is.call(fn)) format(fn) else as.character(fn)
  }
  if(missing(with.variables))
    with.variables <-variables
  else{
    with.variables<-eval(substitute(with.variables),data,parent.frame())
    if(length(dim(with.variables))<1.5){
      with.variables<-ddf(with.variables)
      fn <- arguments$with.variables
      names(with.variables)<-if(is.call(fn)) format(fn) else as.character(fn)
    }		
  }
  cors<-list()
  for(var1 in colnames(variables)){
    cors[[var1]]<-list()
    for(var2  in colnames(with.variables)){
      tmp<-na.omit(data.frame(as.numeric(variables[[var1]]),as.numeric(with.variables[[var2]])))
      names(tmp)<-c(var1,var2)
      cors[[var1]][[var2]]<-test(tmp[[1]],tmp[[2]],...)
      attr(cors[[var1]][[var2]],"N")<-nrow(tmp)
    }
  }
  class(cors)<-"cor.matrix"
  cors
}

ggcorplot <- function(cor.mat,data=NULL,lines=TRUE,line.method=c("lm","loess"),type="points",
                      alpha=.25,main="auto",var_text_size=5,
                      cor_text_limits=c(5,25),level=.05){
  x_var <- y_var <- trans <- rsq <- p <- x_label <- NULL
  #define a helper function (borrowed from the "ez" package)
  ezLev<-function(x,new_order){
    for(i in rev(new_order)){
      x<-relevel(x,ref=i)
    }
    return(x)
  }							
  
  if(all(line.method==c("lm","loess")))
    line.method<-"lm"	
  
  nm <- names(cor.mat)
  for(i in 1:length(nm))
    dat <- if(i==1) ddf(eval(parse(text=nm[i]),data,parent.frame())) else ddf(dat, eval(parse(text=nm[i]),data,parent.frame()))
  data <- dat
  names(data) <- nm
  # normalize data
  for(i in 1:length(data)){
    data[,i]<-as.numeric(data[,i])
    data[,i]<-(data[,i]-mean(data[,i],na.rm=TRUE))/sd(data[,i],na.rm=TRUE)
  }
  # obtain new data frame
  z<-data.frame()
  i <- 1
  j <- i
  while(i<=length(data)){
    if(j>length(data)){
      i<-i+1
      j<-i
    }else{
      x <- data[,i]
      y <- data[,j]
      temp<-as.data.frame((cbind(x,y)))
      temp<-cbind(temp,names(data)[i],names(data)[j])
      z<-rbind(z,temp)
      
      j<-j+1
    }
  }
  z<-cbind(z,alpha)
  names(z)=c('x_var','y_var','x_label','y_label','trans')
  z$x_label <- ezLev(factor(z$x_label),names(data))
  z$y_label <- ezLev(factor(z$y_label),names(data))
  z=z[z$x_label!=z$y_label,]
  #obtain correlation values
  z_cor <- data.frame()
  i <- 1
  j <- i
  while(i<=length(data)){
    if(j>length(data)){
      i<-i+1
      j<-i
    }else{
      x <- na.omit(data[,i])
      y <- na.omit(data[,j])
      x_mid <- min(x)+diff(range(x))/2
      y_mid <- min(y)+diff(range(y))/2
      this_cor <- cor.mat[[i]][[j]]$estimate
      this_cor.test <- cor.mat[[i]][[j]]
      this_col <- ifelse(this_cor.test$p.value<level,"red"
                         ,"blue")
      this_size <- (this_cor)^2
      cor_text <- ifelse(
        this_cor>0
        ,substr(format(c(this_cor,.123456789),digits=2)[1],2,4)
        ,paste('-',substr(format(c(this_cor,.123456789),digits=2)[1],3,5),sep='')
      )
      b<-as.data.frame(cor_text)
      b<-cbind(b,x_mid,y_mid,this_col,this_size,names(data)[j],names(data)[i])
      z_cor<-rbind(z_cor,b)
      j<-j+1
    }
  }
  names(z_cor)<-c('cor','x_mid','y_mid','p','rsq','x_label','y_label')
  z_cor$x_label <- ezLev(factor(z_cor$x_label),names(data))
  z_cor$y_label <- ezLev(factor(z_cor$y_label),names(data))
  diag <- z_cor[z_cor$x_label==z_cor$y_label,]
  z_cor<-z_cor[z_cor$x_label!=z_cor$y_label,]
  
  #start creating layers	
  points_layer <- geom_point(aes(x = x_var, y = y_var, alpha=trans),data = z)
  bin_layer<-geom_hex(data = z, mapping = aes(x = x_var, y = y_var,alpha=trans),bins=10)
  lm_line_layer <- stat_smooth(aes(x=x_var, y=y_var), method=line.method)
  cor_text <- geom_text(aes(x=y_mid, y=x_mid, label=cor, size = rsq, colour = p), data=z_cor)
  var_text <- geom_text(aes(x=y_mid, y=x_mid, label=x_label), data=diag, size=var_text_size)
  
  f <- facet_grid(y_label~x_label,scales='free')
  o <- theme(
    panel.grid.minor = element_blank()
    ,panel.grid.major = element_blank()
    ,axis.ticks = element_blank()
    ,axis.text.y = element_blank()
    ,axis.text.x = element_blank()
    ,axis.title.y = element_blank()
    ,axis.title.x = element_blank()
    ,legend.position='none'
  )
  
  size_scale <- scale_size(limits = c(0,1),range=cor_text_limits)
  the.plot<-ggplot(data=z)
  if(type=="bins")
    the.plot<-the.plot+bin_layer
  else if(type=="points")
    the.plot<-the.plot+points_layer + scale_alpha_identity()
  the.plot<-the.plot+var_text+
    cor_text+
    f+
    o+
    size_scale
  if(type=="bins")
    the.plot<-the.plot+scale_fill_gradient(low="grey", high="black")
  if(lines)
    the.plot<-the.plot+lm_line_layer
  if(main=="auto")
    main<-cor.mat[[1]][[1]]$method
  the.plot<-the.plot+ggtitle(main)
  return(the.plot)
}

output.path <- "output"
if (!file.exists(output.path)) {
  dir.create(output.path)
}

my.dup <- function(x) {
  duplicated(x, fromLast=TRUE) | duplicated(x, fromLast=FALSE)
}

suppressPackageStartupMessages(library("scales")) ## for scientific

Sys.setenv(R_ZIPCMD= "/usr/bin/zip")  

num.cores <- detectCores()
if(!is.na(num.cores) && (num.cores > 1)) {
  suppressPackageStartupMessages(library("doMC"))
  cat(paste("Registering ", num.cores-1, " cores.\n", sep=""))
  registerDoMC(cores=(num.cores-1))
}

cat("Fit _all_ OHSU drug response data\n")

## BSW add min response (min) as a parameter
dss <- function(ic50,slope,max,min.conc.tested,max.conc.tested,y.arg=10,DSS.type=2,concn_scale=1e-9, min.response=0){
  #rdata should be in in format containing IC50, SLOPE, MAX,MIN.Concentration,MAX.Concentration
  
  ## BSW added:
  ## If the min response is 0, then make the activity threshold = y.arg (e.g,. 10% by default)
  ## If the min response is not 0, then make the activity threshold 10% (by default) of the range (max - min.response) _above_ the min.response
  y <- as.numeric(unname(y.arg))
  if(max < as.numeric(unname(min.response))) { 
    dss <- 0
    return(dss)
  }
  if(as.numeric(unname(min.response)) != 0) {
    range <- max - as.numeric(unname(min.response))
    y <- min(max, as.numeric(unname(min.response)) + (((as.numeric(unname(y.arg)))/100) * range))
  }
  
  a=as.numeric(unname(max))
  
  b=as.numeric(unname(slope))
  d=as.numeric(unname(min.response)) # min response
  ic50 = as.numeric(unname(ic50))
  min.conc.tested = as.numeric(unname(min.conc.tested))
  max.conc.tested = as.numeric(unname(max.conc.tested))
  Min.Conc<- log10(min.conc.tested*concn_scale) #
  Max.Conc<- max.conc.tested
  x2<-log10(Max.Conc*concn_scale)  
  
  
  if(is.na(ic50)||is.na(b)||is.na(a)||is.na(Min.Conc)||is.na(Max.Conc)){
    dss<-NA
  }
  else if(isTRUE(ic50>=Max.Conc)){
    dss<-0
  }
  ### BSW added
  else if(isTRUE(ic50<=min.conc.tested)){
    dss<-0
  }
  else if(isTRUE(b==0)){
    dss<-0
  }
  else{
    if(a>100){ a<-100  }
    if(isTRUE(b<0)){ b<--b  }
    c<-log10(ic50*concn_scale)
    if(a>y){
      if(y!=0){
        x1<-(c - ((log(a-y)-log(y-d))/(b*log(10))))
        if(isTRUE(x1 < Min.Conc)){x1<-Min.Conc}
        else if(isTRUE(x1 > x2)){x1<-x2}
      }
      else {x1<-Min.Conc}
      
      # This is a logistic function used in Dotmatics.com
      # y = d+(a-d)/(1+10^(b*(c-x)))
      #inverse function
      # x = c - ((log(a-y)-log(d-y))/(b*log(10)))
      
      int_y=(((((a-d)*log(1+10^(b*(c-x2))))/(b*log(10)))+a*x2)-((((a-d)*log(1+10^(b*(c-x1))))/(b*log(10)))+a*x1)) - (y*(x2-x1))
      
      total_area<-(x2-Min.Conc)*(100-y)
      
      if(DSS.type==1){
        norm_area<-((int_y/total_area)*100)#DSS1
      }
      if(DSS.type==2){
        #       if(a>100){a<-100}
        norm_area<-((int_y/total_area)*100)/log10(a)#DSS2 #AUC1
        if(isTRUE(norm_area > 50)){ norm_area <- 0}
      }
      if(DSS.type==3){
        #       if(a>100){a<-100}
        norm_area<-((int_y/total_area)*100)*(log10(100)/log10(a))*((x2-x1)/(x2-Min.Conc)) #DSS3 #AUC5
      }
      if(isTRUE(norm_area < 0|norm_area > 100)){
        dss<-0
      }else{
        dss<-round(norm_area,digits=4)}
    } else {dss<-0} 
  } 
  return (dss)
}

drc.calc.ic50.robust <- function(dose, inhibition, file, patient.id, inhibitor, return.fit = TRUE, return.residuals = FALSE, apply.constraints = FALSE, plot.fit = FALSE, robust = c("mean", "median", "lts", "lms", "trimmed", "winsor", "tukey"), force.min.to.zero = FALSE) {
      
      ## Jitter any duplicated data to avoid singularities during fitting
      dupes <- duplicated(inhibition)
      if(any(dupes)) {
        inhibition[dupes] <- inhibition[dupes] + rnorm(n = length(which(dupes)), mean = 0, sd = 0.0001)
      }
      
      dupes <- duplicated(dose)
      if(any(dupes)) {
        dose[dupes] <- dose[dupes] + rnorm(n = length(which(dupes)), mean = 0, sd = 0.0001)
      }
      
      robust <- match.arg(robust)
      
      ## Combine the data and sort by dose
      mat_tbl <- data.frame(logconc = log10(dose), dose = dose, inhibition = inhibition)
      mat_tbl <- mat_tbl[order(mat_tbl$logconc), ]
      
      ## Calculate the IC50.
      
      ## Use drc to get the initial parameters.  
      ## FIMM function:
      ##   y = d + ( a - d ) / ( 1 + ( 10^( b * ( c - x ) ) ) )
      ## Asymptotics for positive b: 
      ##              x -> infty --> y -> a
      ##              x -> -infty --> y -> d
      ## In an inhibition setting, FIMM sets b = 1, d (min asymptote) = 0, and determines a and c via the linear regression self starter
      ## drc function:
      ##   y = c' + ( d' - c' ) / ( 1 + ( exp( b' * ( log(x) - log(e') ) ) ) )
      ## here c' is the lower asymptote, d' is the upper asymptote, and b' is the slope (see Ritz et al). 
      ## Note that the drc slope b' has the opposite sign as the FIMM slope b.
      ## Hence, to ensure that y is a decreasing function (e.g., suitable for viability response where viability
      ## decreases with concentration) and c' is the lower asymptote, we need b' to be positive.
      ## Asymptotics:
      ##    x -> infy, b' < 0 --> y -> d'
      ##    x -> infy, b' > 0 --> y -> c'
      ##    x ->    0, b' < 0 --> y -> c'
      ##    x ->    0, b' > 0 --> y -> d'
      ## Hence, for an increasing function (x increases -> y increases), we will have b' < 0
      ## so that d' is the upper asymptote.  For a decreasing function (x increases -> y decreases),
      ## we will have b' > 0, so that c' is the lower asymptote.
      ## And this is stated by Ritz et al. here:
      ## 'the slope paramter b [what I'm called b'] should be negative for an increasing dose-response relationship
      ## to ensure the interpretation that c is the lower asymptote and d the upper asymptote.'
      
      ## Limits for SLOPE, MIN, MAX, and IC50  
      lowerl <- c(-Inf, -Inf, -Inf, -Inf)
      upperl <- c(Inf, Inf, Inf, Inf)
      
      x <- mat_tbl$logconc
      yProp <- mat_tbl$inhibition
      
      if(apply.constraints) {
        ## Constrain IC50 to be between the min and max concentrations
        lowerl[4] <- min(x)
        upperl[4] <- max(x)
        
        ## Constrain bottom to be between [min response - 10% of range, max response]
        if(max(yProp) != min(yProp)) {
          lowerl[2] <- min(yProp) - 0.1 * abs(max(yProp)-min(yProp))
          upperl[2] <- max(yProp)
        }
        
        ## Constrain top to be bewteen [min response, max response + 10% of range]
        if(max(yProp) != min(yProp)) {
          lowerl[3] <- min(yProp)
          upperl[3] <- max(yProp) + 0.1 * abs(max(yProp)-min(yProp))
        }
      }
      if(force.min.to.zero == TRUE) {
        lowerl[2] <- 0
        upperl[2] <- 0.0001
      }
      ## For inhibition data (which is what we fit with drc), the function y will be increasing.  
      ## Hence, b' should be negative.
      ## In any case, c' and d' are the min and max, respectively.
      ## Parameters are in the order names = c("b'", "c'", "d'", "e'")
      ##  print(mat_tbl)
      ##  print(lowerl)
      ##  print(upperl)
      ##  print(robust)
      drc.fit <- NULL
      tryCatch({
        drc.fit <- suppressWarnings(drm(inhibition ~ logconc, data = mat_tbl, 
                                        fct = LL.4(fixed = c(NA, NA, NA,NA), names = c("SLOPE","MIN","MAX","IC50")),
                                        lowerl = lowerl,
                                        upperl = upperl,
                                        robust = robust,
                                        logDose = 10, 
                                        control = drmc(errorm = FALSE)))
      }, 
      error = function(e) {
        cat("error function\n")
        drc.fit <- tryCatch({drm(inhibition ~ logconc, data = mat_tbl, 
                                 fct = L.4(fixed = c(NA, NA, NA,NA), names = c("SLOPE","MIN","MAX","IC50")),
                                 lowerl = lowerl,
                                 upperl = upperl,
                                 robust = robust,
                                 logDose=10)
        }, error = function(e) { NULL })
      }
      ##  ,warning = function(w) {
      ##    cat("warn function\n")
      ##    drm(inhibition ~ logconc, data = mat_tbl, 
      ##        fct = L.4(fixed = c(NA, NA, NA,NA), names = c("SLOPE","MIN","MAX","IC50")),
      ##        lowerl = lowerl,
      ##        upperl = upperl,
      ##        robust = robust,
      ##        logDose=10)
      ##  }
      )
      converged <- TRUE
      if(is.null(drc.fit)) {
        converged <- FALSE
      } else {
        if(!is.null(drc.fit$convergence)) {
          converged <- drc.fit$convergence
        }
        if(!is.null(drc.fit$fit$convergence)) {
          converged <- drc.fit$fit$convergence
        }
      }
      
      if(converged == FALSE) {
        drc.fit <- tryCatch({
          drm(inhibition ~ logconc, data = mat_tbl, 
              fct = L.4(fixed = c(NA, NA, NA,NA), names = c("SLOPE","MIN","MAX","IC50")),
              lowerl = lowerl,
              upperl = upperl,
              robust = robust,
              logDose=10)
        },
        error = function(e) {
          NULL
        })
      }
      converged <- TRUE
      if(is.null(drc.fit)) {
        converged <- FALSE
      } else {
        if(!is.null(drc.fit$convergence)) {
          converged <- drc.fit$convergence
        }
        if(!is.null(drc.fit$fit$convergence)) {
          converged <- drc.fit$fit$convergence
        }
      }
      if(converged == FALSE) {
        l <- list(tbl = NULL)
        if(return.fit) {
          cat("Reurnin fit\n")
          l <- list(tbl = NULL, res = NULL)
        }
        return(l)
      }
      
      ## Extract and name the coefficients
      coef_estim <- coef(drc.fit)
      ## b is slope; c is min; d is max, and e is IC50
      names(coef_estim) <- c("SLOPE","MIN","MAX","IC50")
      
      ## DSS inverts the sign of the slope [see above and compare Yadav et al. supp Eq(1) to Ritz et al. (2)] 
      coef_estim["SLOPE"] <- coef_estim["SLOPE"] * -1 
      
      do.corner.cases <- FALSE
      if(do.corner.cases == FALSE) {
        # NB: switching IC50 to log10 scale here
        coef_estim["IC50"] <- log10(coef_estim["IC50"])
      } else {
        ## Handle a bunch of corner cases.  This is taken verbatim from Swapnil's curve_Fitting.R.  I don't necessarily follow the logic.
        ## BEGIN CORNER CASES.
        # if curve decreases or IC50 is higher than max (i.e. IC50 is "outlier"), set IC50 to max conc.
        coef_estim["IC50"] <- ifelse(coef_estim["MAX"]<=coef_estim["MIN"] | coef_estim["IC50"]>max(mat_tbl$dose,na.rm=T), max(mat_tbl$dose,na.rm=T),coef_estim["IC50"])
        # if IC50 is less than 0 set it to min. conc. and if even min. conc. < 0, then set IC50 to mean of all conc.
        coef_estim["IC50"] <- ifelse(coef_estim["IC50"]<0,min(mat_tbl$dose,na.rm=T),coef_estim["IC50"])
        coef_estim["IC50"] <- ifelse(coef_estim["IC50"]<0,mean(mat_tbl$dose,na.rm=T),coef_estim["IC50"])
      
        # NB: switching IC50 to log10 scale here
        coef_estim["IC50"] <- log10(coef_estim["IC50"])
      
        ## Brian: check this--I believe that dss will not calculate dss if IC50 is set to max log concentration, which is why
        ## this is done here.  But, this is not pass directly to DSS, but to nls.
        # similar to previous step but now compare log10(IC50) with log(min. conc.).
        coef_estim["IC50"] <- ifelse(coef_estim["IC50"]<min(mat_tbl$logconc),max(mat_tbl$logconc),coef_estim["IC50"])
        # if all inhib. < 0 set IC50 to max. log. conc !!!!! not obvious why!
        ## Brian: again, may like above?
        coef_estim["IC50"] <- ifelse(all(mat_tbl$inhibition<0),max(mat_tbl$logconc,na.rm=T),coef_estim["IC50"])
        #(Trying to fix curves that need outlier kickout)
        coef_estim["MIN"] <- 0; coef_estim["MAX"] <- max(mat_tbl$inhibition,na.rm=T)
        #(Fix off minimums) Find lowest inhibition value. If it is not in (0:100), fix it whether to 0 or 99.  
        min_lower <- ifelse(min(mat_tbl$inhibition,na.rm=T) > 0,min(mat_tbl$inhibition,na.rm=T),0)
        min_lower <- ifelse(min_lower >= 100,99,min_lower)
        #similar to previous step but for MAX
        coef_estim["MAX"] <- ifelse(coef_estim["MAX"]>100,100,coef_estim["MAX"])
        coef_estim["MAX"] <- ifelse(coef_estim["MAX"]<0,100,coef_estim["MAX"])
        #max_lower and max_upper - lower and upper bounds for 'nl2sol' algorithm in nonlinear least-squares
        max_lower <- ifelse(max(mat_tbl$inhibition,na.rm=T)>100,coef_estim["MAX"],max(mat_tbl$inhibition,na.rm=T))
        max_lower <- ifelse(max_lower < 0,coef_estim["MAX"],max(mat_tbl$inhibition,na.rm=T))
        max_lower <- ifelse(max_lower < 0,0,max_lower)
        max_lower <- ifelse(max_lower > 100,100,max_lower)
        #(Fix upper maximum for negative slopes)
        run_avg <- runmean(mat_tbl$inhibition, 10)
        max_upper <- ifelse(any(run_avg[-nrow(mat_tbl)]>run_avg[nrow(mat_tbl)]),max(mat_tbl$inhibition[run_avg>run_avg[nrow(mat_tbl)]]),coef_estim["MAX"])
        max_upper <- ifelse(any(mat_tbl$inhibition > max_upper),mean(mat_tbl$inhibition[mat_tbl$inhibition > max_upper])+5,max_upper)
        max_upper <- ifelse(max_upper < 0,coef_estim["MAX"],max_upper)
        max_upper <- ifelse(max_upper > 100,100,max_upper) #coef_estim["MAX"]
        max_upper <- ifelse(max_lower > max_upper,coef_estim["MAX"],max_upper)
        # left it as it was, just rewritten a bit (ALEKS). not clear how values 25, 60 and 5 are chosen. 
        mean_inh_last = mean(tail(mat_tbl$inhibition,2),na.rm=T)
        if(mean_inh_last < 60) {
          if(mean_inh_last > 25) coef_estim["IC50"] <- mean(mat_tbl$logconc,na.rm=T)
          else if(mean_inh_last < 25) coef_estim["IC50"] <- max(mat_tbl$logconc,na.rm=T)}
        if(mean(mat_tbl$inhibition[1:3],na.rm=T)<5) coef_estim["IC50"] <- max(mat_tbl$logconc,na.rm=T)
        #add a bit of positive noise to MAX if it is the same as MIN. 
        if(unname(coef_estim["MIN"]) == unname(coef_estim["MAX"])) coef_estim["MAX"] <- coef_estim["MAX"] + rnorm(n=1, mean = 0, sd = 0.001)
      
      ## END CORNER CASES.
      }    
      ## Regardless of the theoretical lower asymptote (with non-log scale x -> -Inf),
      ## the min value that can be obtained occurs as (non-log scale) x -> 0.
      coef_estim["MIN"] <- drc.fit$fct$fct(0, t(as.matrix(drc.fit$fit$par)))
      coef_estim["MAX"] <- drc.fit$fct$fct(Inf, t(as.matrix(drc.fit$fit$par)))
      
      #Calculate the standard error scores
      sumIC50 = summary(drc.fit); 
      ic50std_Error <- sumIC50$coefficients["IC50:(Intercept)","Std. Error"]
      ## residual standard error
      ic50std_resid <- round(sqrt(sumIC50$resVar),1)  
      
      # plot IC50
      if(plot.fit) {
        png(paste0(output.path, "/", file, ".png"))
        x <- seq(min(mat_tbl$logconc),max(mat_tbl$logconc), length=100)
        y <- predict(drc.fit, data.frame(logconc=x))
        icpl <- ggplot(mat_tbl, aes(logconc, inhibition)) + scale_x_continuous(breaks=mat_tbl$logconc,labels=mat_tbl$dose) +
          geom_point(color = "blue", size = 2.8) + geom_line(data = data.frame(x = x, y = y), aes(x, y), color="blue", size = 0.8) +
          geom_vline(xintercept = log10(coef_ic50["IC50"]), colour="grey", size = 0.8) 
        icpl <- icpl + ggtitle(file)
        ## icpl <- icpl + ggtitle(paste0(drug_name," (dss:",round(IC50_dataframe$DSS,1),")\n"))
        icpl <- icpl + theme(legend.title = element_text(size = 10)) + theme_bw() + labs(y = "% inhibition", x = "conc(nM)")  +  ylim(min(-25, 1.1*min(inhibition), 1.1*min(y)),max(125,1.1*max(inhibition),1.1*max(y))) +
          geom_text(mapping=aes(x2,y2,label = text2), data=data.frame(x2=log10(coef_ic50["IC50"])*0.95, y2=115, text2="IC50"), color="grey", parse=T) +
          theme(plot.background = element_rect(fill = "transparent",colour = NA),
                panel.background =element_rect(fill = "transparent",colour = NA))
        
        print(icpl)
        d <- dev.off()
      }
      
      ## mean squared error
      max_signal <- max(mat_tbl$dose,na.rm=T); min_signal <- min(mat_tbl$dose,na.rm=T)
      mse <- sum(residuals(drc.fit)^2)/length(residuals(drc.fit))
      ssres <- sum(residuals(drc.fit)^2)
      resp <- na.omit(drc.fit$dataList$resp)
      sstot <- sum((mean(resp)-resp)^2)
      gof <- 1 - (ssres/sstot)
      coef_ic50 <- coef_estim
      dss_score <- dss(coef_ic50["IC50"],coef_ic50["SLOPE"],coef_ic50["MAX"],min_signal,max_signal, concn_scale=1e-6, min.response=coef_ic50["MIN"]) 
      ##  IC50_dataframe <- data.frame(ID=file,PATIENT_ID=patient.id,DRUG_NAME=inhibitor,t(as.matrix(coef_estim)), MIN.lb=lowerl[2], MIN.ub=upperl[2], MAX.lb=lowerl[3], MAX.ub=upperl[3], IC50.lb=lowerl[4], IC50.ub=upperl[4], SLOPE.lb=lowerl[1], SLOPE.ub=upperl[1], mse=mse)
      ####coef_ic50 <- c(tmp,Max.Conc.tested=max_signal,Min.Conc.tested=min_signal,DSS=dss_score,IC50_std_error=ic50std_Error,S_est = ic50std_resid)
      ###working
      IC50_dataframe <- data.frame(t(as.matrix(coef_estim)), Max.Conc.tested=max_signal,Min.Conc.tested=min_signal, DSS=dss_score, MIN.lb=lowerl[2], MIN.ub=upperl[2], MAX.lb=lowerl[3], MAX.ub=upperl[3], IC50.lb=lowerl[4], IC50.ub=upperl[4], SLOPE.lb=lowerl[1], SLOPE.ub=upperl[1], gof=gof)
      
      l <- list(tbl = IC50_dataframe)
      if(return.fit) {
        l <- list(tbl = IC50_dataframe, res = NULL)
      }
      if(return.residuals) {
        l$residuals <- data.frame(residuals = residuals(drc.fit), logconc = drc.fit$data$logconc, inhibition = drc.fit$data$inhibition, fitted = fitted(drc.fit))    
      }
      return(l)
      
}
    
## end drc.calc.ic50.robust

## BEGIN load FIMM data

## Utility function to read in a file that has rownames (i.e., without a corresponding column)
fread.with.rownames <- function(file, ...) {
  data.set <- fread(file, header=TRUE, fill=TRUE, ...)
  df <- as.data.frame(data.set)
  rownames(df) <- df[,1]
  cols <- colnames(df)[1:(ncol(df)-1)]
  df <- df[,-1]
  colnames(df) <- cols
  data.set <- df
  rm(df)
  data.set
}



## Load (batch corrected) FIMM expression data
obj <- synGet(id="syn8270602", downloadFile = TRUE, downloadLocation = ".")
file <- getFileLocation(obj)
##fimm.expr <- read.table(file, header=TRUE, sep=",")
fimm.expr <- fread.with.rownames(file, sep=",")
fimm.expr <- t(fimm.expr)

## Load FIMM metadata (including relapse/refractory/diagnosis)
synId <- "syn8270594"
obj <- synGet(id=synId, downloadFile = TRUE, downloadLocation = ".")
file <- getFileLocation(obj)
fimm.metadata <- read.table(file, header=TRUE, sep=",")

## Load FIMM genomic data
synId <- "syn8270591"
obj <- synGet(id=synId, downloadFile = TRUE, downloadLocation = ".")
file <- getFileLocation(obj)
fimm.genomic.detail <- read.table(file, header=TRUE, sep=",", comment.char="", quote="\"")

## The "bin" file simply is 0/1 matrix indicating whether the sample (row) has a mutation
## in the gene (column).  These use gene symbols.
synId <- "syn8270590"
obj <- synGet(id=synId, downloadFile = TRUE, downloadLocation = ".")
file <- getFileLocation(obj)
fimm.genomic.bin <- read.table(file, header=TRUE, sep=",", comment.char="", quote="\"")
## Make the samples the columns
fimm.genomic.bin <- t(fimm.genomic.bin)

## Read in FIMM raw drug response data
synId <- "syn8488824"
obj <- synGet(id=synId, downloadFile = TRUE, downloadLocation = ".")
file <- getFileLocation(obj)
fimm.raw.dr <- read.xlsx(file)

## Read in FIMM drug annotations
synId <- "syn8434109"
obj <- synGet(id=synId, downloadFile = TRUE, downloadLocation = ".")
file <- getFileLocation(obj)
fimm.drug.annotations <- read.xlsx(file)

## Read in the subset of drug annotations in both FIMM and OHSU
synId <- "syn9731315"
obj <- synGet(id=synId, downloadFile = TRUE, downloadLocation = ".")
ohsu.fimm.drugs <- openxlsx:::read.xlsx(getFileLocation(obj), sheet=1)

## END load FIMM data

## BEGIN load OHSU data

## Load the OHSU diagnostics
obj <- synGet(id="syn7488457", downloadFile = TRUE, downloadLocation = ".")
file <- getFileLocation(obj)
ohsu.response.tbl <- read.xlsx(file, sheet = 5)
flag <- !is.na(ohsu.response.tbl$response) & grepl(pattern="Complete Response", ohsu.response.tbl$response)
ohsu.response.tbl$response[flag] <- "Complete Response"

ohsu.simplified.response.tbl <- ddply(na.omit(ohsu.response.tbl[,c("patient_id", "response")]), "patient_id", .parallel=TRUE,
                                      .fun = function(df) {
                                        if(any(df$response == "Complete Response")) { return("Complete Response")}
                                        if(any(df$response == "Hematologic CR")) { return("Hematologic CR")}
                                        if(any(df$response == "Refractory")) { return("Refractory")}
                                        if(any(df$response == "Supportive/Palliative Care")) { return("Supportive/Palliative Care")}
                                        if(any(df$response == "Unknown")) { return("Unknown")}
                                        warning("What!?")
                                      })
colnames(ohsu.simplified.response.tbl) <- c("patient_id", "response")
ohsu.responses <- ohsu.simplified.response.tbl$response
names(ohsu.responses) <- ohsu.simplified.response.tbl$patient_id

obj <- synGet(id="syn8441027", downloadFile = TRUE, downloadLocation = ".")
inhibitor.plate.number.file <- getFileLocation(obj)
inhibitor.plate.tbl <- unique(read.xlsx(inhibitor.plate.number.file))
inhibitor.plate.tbl <- inhibitor.plate.tbl[order(inhibitor.plate.tbl$Drug, inhibitor.plate.tbl$Replicant),]

## Read in the raw OHSU inhibitor data inhibitor_data_points_2017_01_12.txt
inh.obj <- synGet("syn8149180", downloadFile=TRUE)
ohsu.inh.tbl <- fread(getFileLocation(inh.obj))
ohsu.inh.tbl <- as.data.frame(ohsu.inh.tbl)

## Load OHSU expression data
## Download the AML CPM expression data
## NB: linking to the original Beat AML project (syn2942337), since this has more updated files/releases
## (as well as the raw fastqs)
## Read in BeatAML_DR1_RNASeq_log2_cpm_2016_08_02.csv
rna.dr1.obj <- synGet("syn7124177", downloadFile=TRUE)
# Load the data
## rna.dr1.tbl <- read.table(getFileLocation(rna.dr1.obj), sep=",", header=TRUE, as.is=TRUE)
rna.dr1.tbl <- as.data.frame(fread(getFileLocation(rna.dr1.obj), sep=",", header=TRUE))
rna.dr1.lab.ids <- gsub("\\.", "-", gsub("^X", "", names(rna.dr1.tbl[-c(1:9)])))

cat("Downloading and reading DR2 Q1 data\n")
## Read in BeatAML_DR2q1_RNASeq_log2_cpm_2016_12_16.csv
rna.dr2q1.obj <- synGet("syn8149220", downloadFile=TRUE)
## rna.dr2q1.tbl <- read.table(getFileLocation(rna.dr2q1.obj), sep=",", header=TRUE, as.is=TRUE)
rna.dr2q1.tbl <- as.data.frame(fread(getFileLocation(rna.dr2q1.obj), sep=",", header=TRUE))
rna.dr2q1.lab.ids <- gsub("\\.", "-", gsub("^X", "", names(rna.dr2q1.tbl[-c(1:9)])))

cat("Downloading and reading DR2 Q2 data\n")
## Read in BeatAML_DR2q2_RNASeq_log2_cpm_2017_01_11.csv
rna.dr2q2.obj <- synGet("syn8149259", downloadFile=TRUE)
## rna.dr2q2.tbl <- read.table(getFileLocation(rna.dr2q2.obj), sep=",", header=TRUE, as.is=TRUE)
rna.dr2q2.tbl <- as.data.frame(fread(getFileLocation(rna.dr2q2.obj), sep=",", header=TRUE))
rna.dr2q2.lab.ids <- gsub("\\.", "-", gsub("^X", "", names(rna.dr2q2.tbl[-c(1:9)])))

## Make sure we have a common set of genes across all 3 expression data sets
common.genes <- intersect(rna.dr1.tbl$Gene, intersect(rna.dr2q1.tbl$Gene, rna.dr2q2.tbl$Gene))
rownames(rna.dr1.tbl) <- rna.dr1.tbl$Gene
rownames(rna.dr2q1.tbl) <- rna.dr2q1.tbl$Gene
rownames(rna.dr2q2.tbl) <- rna.dr2q2.tbl$Gene
rna.dr1.tbl <- rna.dr1.tbl[common.genes,]
rna.dr2q1.tbl <- rna.dr2q1.tbl[common.genes,]
rna.dr2q2.tbl <- rna.dr2q2.tbl[common.genes,]

## Drop the annotation columns from the expression matrices
rna.dr1.tbl <- rna.dr1.tbl[,-c(1:9)]
rna.dr2q1.tbl <- rna.dr2q1.tbl[,-c(1:9)]
rna.dr2q2.tbl <- rna.dr2q2.tbl[,-c(1:9)]
gc()

all.samples <- c(colnames(rna.dr1.tbl), colnames(rna.dr2q1.tbl), colnames(rna.dr2q2.tbl))
if(any(duplicated(all.samples))) {
  warning("Some RNA-seq samples are duplicated across the 3 releases\n")
} else {
  cat("As expected, no samples are duplicated across the 3 releases\n")
}
## Since no samples are duplicated, we can just merge the columns
ohsu.expr <- cbind(rna.dr1.tbl, rna.dr2q1.tbl, rna.dr2q2.tbl)
rm(rna.dr1.tbl)
rm(rna.dr2q1.tbl)
rm(rna.dr2q2.tbl)
gc()

## Read in the latest rna-seq dashboard
## BeatAML_rnaseq_2017_01_19_public_dashboard.xlsx (syn8149218)
## NB: this includes data release 1, release 2 q1, and release 2 q2.
rna.obj <- synGet("syn8149218", downloadFile=TRUE)

## The samples summary is the 4th sheet
## This sample.summary table provides a map from patient id to sequence id.
## Except that since the sample sequence ids where column headers in our expression
## data, they were corrupted from AA-BBBBB (where A and B are digits) to XAA.BBBBB.
## So change that.
cat(paste0("Reading sample summary file ", getFileLocation(rna.obj), "\n"))
rnaseq.sample.summary.tbl <- openxlsx:::read.xlsx(getFileLocation(rna.obj), sheet=4)

## Read in the OHSU patient diagnosis table (corresponding to the drug response data).
obj <- synGet("syn8149174", downloadFile=TRUE)
ohsu.diagnosis.tbl <- openxlsx:::read.xlsx(getFileLocation(obj), sheet=1)

## Load OHSU genomic data (just mutect for now)
origBeatMafs     <- synQuery("SELECT name FROM file WHERE parentId=='syn5522834'")
beatMafs     <- origBeatMafs[grep("mutect",origBeatMafs$file.name),,drop=FALSE]
beatMafNames <- gsub("^.*_(.*?)_AML_.*maf$","\\1",beatMafs$file.name)

## Read in the OHSU seqcap dashboard file (so we can subset to those patients that have AML)
## For DNA sheet 4 is all samples (tumor and normal)
## Sheet 5 is tumor only.
## Tumor-only means that there was no normal comparator.
## I believe that sheet 5 is a strict subset of sheet 4, but to be be sure
## merge the two.
## Read in the latest dna-seq dashboard
## BeatAML_seqcap_2017_01_19_public_dashboard.xlsx (syn8149211)
dna.obj <- synGet("syn8149211", downloadFile=TRUE)
cat(paste0("Reading sample summary file ", getFileLocation(dna.obj), "\n"))
all.dnaseq.sample.summary.tbl <- openxlsx:::read.xlsx(getFileLocation(dna.obj), sheet=4)
tumor.only.dnaseq.sample.summary.tbl <- openxlsx:::read.xlsx(getFileLocation(dna.obj), sheet=5)
dnaseq.sample.summary.tbl <- unique(rbind(all.dnaseq.sample.summary.tbl, tumor.only.dnaseq.sample.summary.tbl))

## Make sure all of the beatMaf labids are in the diagnosis table
if(!(all(beatMafNames %in% dnaseq.sample.summary.tbl$SeqID))) {
  warning("UNEXPECTED: some maf lab ids were not in DNA dashboard\n")
} else {
  cat("As expected, all maf lab ids were in DNA dashboard\n")
}

## END load OHSU data

## BEGIN sample filtering

## Filter based on sample diagnosis (e.g., AML only)
## restrict.ohsu.to.aml: whether to restrict analysis to only those specimens with an AML diagnosis
##                       ACUTE MYELOID LEUKAEMIA (AML) AND RELATED PRECURSOR NEOPLASMS
## restrict.fimm.to.aml: whether to restrict analysis to only those specimens with an AML diagnosis: "AML"
## restrict.to.paired.dna.seq: whether to restrict DNA-seq analysis to only those specimens for which
##                       we have a paired normal
## restrict.fimm.disease.stage:  array of strings to which FIMM disease stage should be restricted
##                       (subset of "Diagnosis", "Refractory", and "Relapse") or NULL (indicating
##                       no restriction)
## match.at.specimen.level: in requiring concordance between RNA-seq and drug data (or DNA-seq and
##                          drug data) require specimens be assayed by both (if TRUE), otherwise
##                          only require that patients be assayed by both (even if the specimens
##                          assayed differ).
filter.data.sets <- function(ohsu.dna.summary.tbl, oshu.tumor.only.dna.summary.tbl, ohsu.rna.summary.tbl, ohsu.beat.maf.names, ohsu.beat.mafs, ohsu.expr.data, fimm.expr.data, ohsu.diagnosis.tbl, ohsu.raw.drug.data, fimm.raw.drug.data, fimm.patient.metadata, restrict.ohsu.to.aml = TRUE, restrict.fimm.to.aml = TRUE, restrict.fimm.disease.stage = NULL, restrict.to.paired.dna.seq = TRUE, match.at.specimen.level = FALSE) {
  
  aml.diagnosis <- "ACUTE MYELOID LEUKAEMIA (AML) AND RELATED PRECURSOR NEOPLASMS"
  
  ## Only keep those mafs and drug response data that have an AML diagnosis (as indicated in the dashboard)
  if(restrict.ohsu.to.aml == TRUE) {
    ## Restrict drug response data to AML diagnosis (as annotated in the raw inhibitor file)
    cat(paste0("Only analyzing OHSU drug response specimens with AML diagnosis: ", aml.diagnosis, "\n"))
    ohsu.raw.drug.data <- ohsu.raw.drug.data[ohsu.raw.drug.data$diagnosis == aml.diagnosis,]
    
    ## Further, restrict it to AML in the diagnostic file (as either the most recent or the time of acquisition specimen)
    ohsu.aml.diagnosis.tbl <- subset(ohsu.diagnosis.tbl, most_recent_diagnosis %in% aml.diagnosis)
    ohsu.aml.diagnosis.tbl <- subset(ohsu.aml.diagnosis.tbl, diagnosis_at_time_of_specimen_acquisition %in% aml.diagnosis)
    ohsu.raw.drug.data <- ohsu.raw.drug.data[ohsu.raw.drug.data$lab_id %in% ohsu.aml.diagnosis.tbl$lab_id,]
  }

  if(restrict.to.paired.dna.seq == TRUE) {
    cat(paste0("Only analyzing OHSU DNA-seq specimens having a paired normal\n"))
    flag <- !(ohsu.dna.summary.tbl$SeqID %in% oshu.tumor.only.dna.summary.tbl$SeqID)
    ohsu.dna.summary.tbl <- ohsu.dna.summary.tbl[flag, ]
  }
  
  ## Resrict DNA-seq data to those for which we have drug response (after potential filter based on AML diagnosis)
  cat(paste0("Only analyzing OHSU DNA-seq specimens having (diagnosis-filtered) drug response\n"))
  if(match.at.specimen.level == TRUE) {
    ## NB: we are requiring the same _specimen_ be assayed by DNA-seq and for drug response,
    ## not solely the same patient.
    ## NB: the lab_id in this diagnosis table corresponds to Original_LabID in the dna dashboard file
    cat("Requiring concordance between DNA-seq and drug response based on specimen/lab ids\n")
    ohsu.dna.summary.tbl <- ohsu.dna.summary.tbl[(ohsu.dna.summary.tbl$Original_LabID %in% ohsu.raw.drug.data$lab_id),]
  } else {
    cat("Requiring concordance between DNA-seq and drug response based on patient ids\n")
    ohsu.dna.summary.tbl <- ohsu.dna.summary.tbl[(ohsu.dna.summary.tbl$PatientID %in% ohsu.raw.drug.data$patient_id),]
  }  
  flag <- ohsu.beat.maf.names %in% ohsu.dna.summary.tbl$SeqID
  ohsu.beat.maf.names <- ohsu.beat.maf.names[flag]
  ohsu.beat.mafs <- ohsu.beat.mafs[flag,]
  
  cat(paste0("Only analyzing OHSU RNA-seq specimens having (diagnosis-filtered) drug response\n"))
  ## Resrict RNA-seq data to those for which we have drug response (after potential filter based on AML diagnosis)
  if(match.at.specimen.level == TRUE) {
    cat("Requiring concordance between RNA-seq and drug response based on specimen/lab ids\n")
    ohsu.rna.summary.tbl <- ohsu.rna.summary.tbl[(ohsu.rna.summary.tbl$Original_LabID %in% ohsu.raw.drug.data$lab_id),]
  } else {
    ohsu.rna.summary.tbl <- ohsu.rna.summary.tbl[(ohsu.rna.summary.tbl$PatientID %in% ohsu.raw.drug.data$patient_id),]
  }
  if(!is.null(restrict.fimm.disease.stage)) {
    restricted.diseases <- restrict.fimm.disease.stage[restrict.fimm.disease.stage %in% fimm.metadata$diseases.stage]
    cat(paste0("Restricting FIMM disease stage to: ", paste(restricted.diseases, collapse=", "), "\n"))
    fimm.patient.metadata <- subset(fimm.patient.metadata, diseases.stage %in% restricted.diseases)
  }
  
  if(restrict.fimm.to.aml) {
    cat(paste("Restricting FIMM diagnosis to having the pattern \"AML\"\n"))
    fimm.patient.metadata <- subset(fimm.patient.metadata, grepl(diagnosis, pattern="AML"))
  }
  
  fimm.raw.drug.data <- subset(fimm.raw.drug.data, (SCREEN_ID %in% rownames(fimm.patient.metadata)))

  ## Filter the expression data based on the RNA-seq filtering above
  ohsu.expr.lab.ids <- gsub("\\.", "-", gsub("^X", "", names(ohsu.expr.data)))
  flag <- ohsu.expr.lab.ids %in% ohsu.rna.summary.tbl$SeqID
  ohsu.expr.data <- ohsu.expr.data[, flag]
  
  if(match.at.specimen.level == TRUE) {
    cat("Requiring concordance between FIMM RNA-seq and drug response based on specimen/lab ids\n")
    flag <- colnames(fimm.expr.data) %in% fimm.raw.drug.data$SCREEN_ID
    fimm.expr.data <- fimm.expr.data[, flag]
  } else {
    cat("Requiring concordance between FIMM RNA-seq and drug response based on patient ids\n")
    flag <- rownames(fimm.patient.metadata) %in% fimm.raw.drug.data$SCREEN_ID
    flag <- fimm.patient.metadata$person %in% fimm.patient.metadata$person[flag]
    flag <- colnames(fimm.expr.data) %in% rownames(fimm.patient.metadata)[flag]
    fimm.expr.data <- fimm.expr.data[, flag]
  }
  ret <- list(ohsu.expr.data = ohsu.expr.data, fimm.expr.data = fimm.expr.data, ohsu.beat.maf.names = ohsu.beat.maf.names, ohsu.beat.mafs = ohsu.beat.mafs, ohsu.raw.drug.data = ohsu.raw.drug.data, fimm.raw.drug.data = fimm.raw.drug.data)
  ret
}

## END sample filtering

## BEGIN drug response outlier sample/individual concentration filtering

## remove.ohsu.zero.responses: whether to remove any OHSU normalized_viabilities
##                             that are exactly zero (this is an error condition)
## max.ohsu.response: exclude anything in the OHSU data above this response
## require.intersection:  require drug be in both FIMM and OHSU

## OHSU filters sent by Dave Edwards
## 1.  Outliers (norm viability > 150)
## 2.  Low quality replicates (>= 3 with norm viability > 150) -- i.e., min.ohsu.responses.in.screen = 5 (since there should be 7 pts)
## Rest will be post-fit filtering
## 1.  Oversized AUC (AUC > 1100)
## 2.  Low replicant R^2 (R^2 < 0.4)
## 3a. High percent std deviation - replicates (CV = sd/mean > 50%): 
## 3b. High percent std deviation - specimen (CV = sd/mean > 50%)
filter.drug.response.data.sets <- function(ohsu.raw.drug.data, fimm.raw.drug.data, remove.ohsu.zero.responses = TRUE, max.ohsu.response = 150, min.ohsu.responses.in.screen = 5, max.ohsu.responses.in.screen = 7, require.intersection = TRUE) {
  
  if(remove.ohsu.zero.responses) {
    flag <- ohsu.raw.drug.data$normalized_viability == 0
    cat(paste0("Removing OHSU normalized viability == 0: ", length(which(flag)), " of ", length(flag), " cases\n"))
    ohsu.raw.drug.data <- ohsu.raw.drug.data[!flag, ]
  }
     
  flag <- ohsu.raw.drug.data$normalized_viability > max.ohsu.response
  cat(paste0("Removing OHSU normalized viability > ", max.ohsu.response, ": ", length(which(flag)), " of ", length(flag), " cases\n"))
  ohsu.raw.drug.data <- ohsu.raw.drug.data[!flag, ]

  if(require.intersection) {
    cat(paste("Requiring drugs be in both FIMM and OHSU\n"))
    ## First do this quick subset before the merge (and defining ids below) for efficiency
    ohsu.raw.drug.data <- subset(ohsu.raw.drug.data, inhibitor %in% ohsu.fimm.drugs$ID_Drug.ohsu)
  }
  
  ## Filter specimen screens to ensure we have at least min.ohsu.responses.in.screen responses/concentrations 
  ## per specimen and at most max.ohsu.responses.in.screen.
  ## The upper bound is in place because we sometimes see a screen mapped to multiple plates; this
  ## effectively excludes those cases.
  old.nrow <- nrow(ohsu.raw.drug.data)
  ohsu.raw.drug.data$ids <- unlist(apply(ohsu.raw.drug.data[,c("inhibitor", "replicant", "lab_id", "patient_id", "time_of_read")], 1, function(row) paste(row, collapse="-")))
  tbl <- as.data.frame(table(ohsu.raw.drug.data$ids))
  colnames(tbl) <- c("ids", "freq")
  tbl <- subset(tbl, (freq >= min.ohsu.responses.in.screen) & (freq <= max.ohsu.responses.in.screen))
  ohsu.raw.drug.data <- subset(ohsu.raw.drug.data, ids %in% tbl$ids)
  cat(paste0("After ensuring num responses per screen >= ", min.ohsu.responses.in.screen, " and <= ", max.ohsu.responses.in.screen, " left with ", nrow(ohsu.raw.drug.data), " data pts out of ", old.nrow, "\n"))

  if(require.intersection) {
    cat(paste("Requiring drugs be in both FIMM and OHSU\n"))
    ## First do this quick subset before the merge for efficiency
    ohsu.raw.drug.data <- subset(ohsu.raw.drug.data, inhibitor %in% ohsu.fimm.drugs$ID_Drug.ohsu)
    ## ID_Drug is the FIMM drug identifier
    ohsu.raw.drug.data <- merge(ohsu.raw.drug.data, unique(ohsu.fimm.drugs[,c("ID_Drug", "ID_Drug.ohsu")]), by.x = "inhibitor", by.y = "ID_Drug.ohsu")
    ohsu.raw.drug.data <- subset(ohsu.raw.drug.data, ID_Drug %in% fimm.raw.drug.data$DRUG_ID)
    fimm.raw.drug.data <- subset(fimm.raw.drug.data, DRUG_ID %in% ohsu.raw.drug.data$ID_Drug)
  }
  
  ret <- list(ohsu.raw.drug.data = ohsu.raw.drug.data, fimm.raw.drug.data = fimm.raw.drug.data)  
  ret
}
  
## END drug response outlier sample/individual concentration filtering

## mutation.types.to.total: list of mutation types that will be considered a mutation (i.e., a "1")
##                          in the final 0/1 matrix.
## e.g., c("Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Splice_Site", "Translation_Start_Site")
read.mafs <- function(maf.names, maf.tbl, mutation.types.to.total = c("Missense_Mutation", "Nonsense_Mutation")) {
  ## Read in the MAFs
  maf.list  <- list()
  for(i in 1:nrow(maf.tbl)){ maf.list[[maf.names[i]]] <-  read.maf(getFileLocation(synGet(maf.tbl$file.id[i])), verbose = F); print(i)}
  
  maf.tbl <- ldply(maf.list, getGeneSummary)
  rm(maf.list); gc()
  
  ## Tally mutation.types.to.total (e.g., Missense_Mutation and Nonsense_Mutation)
  maf.tbl$total <- unlist(apply(maf.tbl[,mutation.types.to.total], 1, function(row) sum(row, na.rm=TRUE)))
  print(colnames(maf.tbl))
  maf.tbl <- maf.tbl[,c(".id","Hugo_Symbol","total")]
  
  maf.tbl                  <- spread(maf.tbl, Hugo_Symbol,total); 
  rownames(maf.tbl)        <- maf.tbl$.id; 
  maf.tbl                  <- maf.tbl[,-1]; 
  maf.tbl[is.na(maf.tbl)] <- 0
  maf.tbl[maf.tbl > 0] <- 1
  maf.tbl
}

cat(paste0("Before filtering raw drug data\n"))
cat(paste0("   FIMM has: ", nrow(fimm.raw.dr), "\n"))
cat(paste0("   OHSU has: ", nrow(ohsu.inh.tbl), "\n"))
filtered.res <- filter.drug.response.data.sets(ohsu.raw.drug.data = ohsu.inh.tbl, fimm.raw.drug.data = fimm.raw.dr, remove.ohsu.zero.responses = TRUE, max.ohsu.response = 150, min.ohsu.responses.in.screen = 5, max.ohsu.responses.in.screen = 7, require.intersection = FALSE)
fimm.raw.drug.data <- filtered.res[["fimm.raw.drug.data"]]
ohsu.raw.drug.data <- filtered.res[["ohsu.raw.drug.data"]]
rm(ohsu.inh.tbl)
rm(fimm.raw.dr)
rm(filtered.res); gc()


cat(paste0("After filtering raw drug data\n"))
cat(paste0("   FIMM has: ", nrow(fimm.raw.drug.data), "\n"))
cat(paste0("   OHSU has: ", nrow(ohsu.raw.drug.data), "\n"))

filtered.res <- filter.data.sets(ohsu.dna.summary.tbl = dnaseq.sample.summary.tbl, oshu.tumor.only.dna.summary.tbl = tumor.only.dnaseq.sample.summary.tbl, ohsu.rna.summary.tbl = rnaseq.sample.summary.tbl, ohsu.beat.maf.names = beatMafNames, ohsu.beat.mafs = beatMafs, ohsu.expr.data = ohsu.expr, fimm.expr.data = fimm.expr, ohsu.diagnosis.tbl = ohsu.diagnosis.tbl, ohsu.raw.drug.data = ohsu.raw.drug.data, fimm.raw.drug.data = fimm.raw.drug.data, fimm.patient.metadata = fimm.metadata, restrict.ohsu.to.aml = TRUE, restrict.fimm.to.aml = TRUE, restrict.fimm.disease.stage = NULL, restrict.to.paired.dna.seq = FALSE, match.at.specimen.level = FALSE)
ohsu.expr.data <- filtered.res[["ohsu.expr.data"]]
fimm.expr.data <- filtered.res[["fimm.expr.data"]]
ohsu.beat.maf.names <- filtered.res[["ohsu.beat.maf.names"]]
ohsu.beat.mafs <- filtered.res[["ohsu.beat.mafs"]]
ohsu.raw.drug.data <- filtered.res[["ohsu.raw.drug.data"]]
fimm.raw.drug.data <- filtered.res[["fimm.raw.drug.data"]]
rm(ohsu.expr)
rm(fimm.expr)
rm(filtered.res); gc()

## This works with plyr 1.7.1 and 1.8.1, but not 1.8.4
fit.fimm.drug.response.data <- function(fimm.drug.response.data) {
fimm.fits <- dlply(fimm.drug.response.data[,c("DRUG_ID", "DRUG_NAME", "CONCENTRATION", "PERCENT_INHIBITION", "SCREEN_ID", "SAMPLE_DATE", "PATIENT_ID")],
                   c("SCREEN_ID", "DRUG_ID"),
                   .parallel = TRUE,
                   .fun = 
                     function(df) {
                       id <- unique(df[, !(colnames(df) %in% c("CONCENTRATION", "PERCENT_INHIBITION")),drop=F])
                       if(nrow(df) < 2) {
                         return(list(id=id, tbl=NULL))
                       }
                       dose <- df$CONCENTRATION
                       inhibition <- df$PERCENT_INHIBITION
                       patient.id <- df$PATIENT_ID[1]
                       inhibitor <- df$DRUG_NAME[1]
                       file <- paste0(patient.id, "-", inhibitor)
                       print(paste0(patient.id, "-", df$DRUG_ID[1]))
                       tbl <- NULL
                       fit <- NULL
                       res <- suppressWarnings(drc.calc.ic50.robust(dose, inhibition, file, patient.id, inhibitor, return.fit = FALSE, return.residuals = FALSE, apply.constraints = FALSE, plot.fit = FALSE, robust = "mean", force.min.to.zero = FALSE))
                       if(!is.null(res) & !is.null(res$tbl)) {
                         tbl <- res$tbl
                       }
                       if(!is.null(res) & !is.null(res$res)) {
                         fit <- res$tbl
                       }
                       list(id=id, tbl=tbl, fit=fit)
                     })
fimm.fits
}

fit.ohsu.drug.response.data <- function(ohsu.drug.response.data) {
ohsu.fits <- dlply(ohsu.drug.response.data[,c("inhibitor", "well_concentration", "normalized_viability", "time_of_read", "lab_id", "replicant", "patient_id")],
                   c("patient_id", "inhibitor", "lab_id", "replicant", "time_of_read"),
                   .parallel = TRUE,
                   .fun = 
                     function(df) {
                       id <- unique(df[, !(colnames(df) %in% c("well_concentration", "normalized_viability")),drop=F])
                       if(nrow(df) < 2) {
                         return(list(id=id, tbl=NULL))
                       }
                       dose <- df$well_concentration
                       inhibition <- 100 - df$normalized_viability
                       patient.id <- df$patient_id[1]
                       inhibitor <- df$inhibitor[1]
                       file <- paste0(patient.id, "-", inhibitor)
                       tbl <- NULL
                       fit <- NULL
                       res <- suppressWarnings(drc.calc.ic50.robust(dose, inhibition, file, patient.id, inhibitor, return.fit = FALSE, return.residuals = FALSE, apply.constraints = FALSE, plot.fit = FALSE, robust = "mean", force.min.to.zero = FALSE))
                       if(!is.null(res) & !is.null(res$tbl)) {
                         tbl <- res$tbl
                       }
                       if(!is.null(res) & !is.null(res$res)) {
                         fit <- res$tbl
                       }
                       list(id=id, tbl=tbl, fit=fit)
                     })
ohsu.fits
}

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

l.4.func <- function(x, b, c, d, e) {
  c + ( (d-c)/(1+exp(b*(x-e))) )
}


## Do correlation of all genes vs eigendrugs
panel.rank.cor <- function(x, y, digits=2, prefix="", cex.cor, ...) 
{
  usr <- par("usr"); on.exit(par(usr)) 
  par(usr = c(0, 1, 0, 1)) 
  ok <- is.finite(x) & is.finite(y) 
  if(length(which(ok)) < 3) { return() }
  test <- cor.test(x[ok], y[ok], method="spearman")
  r <- test$estimate 
  txt <- format(c(r, 0.123456789), digits=digits)[1] 
  txt <- paste(prefix, txt, sep="") 
  if(missing(cex.cor)) cex <- 0.8/strwidth(txt) 
  
  ##test <- cor.test(x[ok],y[ok]) 
  # borrowed from printCoefmat
  Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                   symbols = c("***", "**", "*", ".", " ")) 
  
  text(0.5, 0.5, txt, cex = cex * r) 
  text(.8, .8, Signif, cex=cex, col=2) 
}

panel.regression <- function (x, y, col = par("col"), bg = NA, pch = par("pch"), 
                              cex = 1, col.regres = "red", ...) 
{ 
  points(x, y, pch = pch, col = col, bg = bg, cex = cex) 
  ok <- is.finite(x) & is.finite(y) 
  if (length(which(ok)) > 3) 
    abline(stats::lm(y[ok] ~ x[ok]), col = col.regres, ...) 
} 

panel.rank.regression <- function (x, y, col = par("col"), bg = NA, pch = par("pch"), 
                                   cex = 1, col.regres = "red", ...) 
{ 
  ok <- is.finite(x) & is.finite(y) 
  points(rank(x[ok]), rank(y[ok]), pch = pch, col = col, bg = bg, cex = cex) 
  if (any(ok)) 
    abline(stats::lm(rank(y[ok]) ~ rank(x[ok])), col = col.regres, ...) 
} 



plot.pairs <- function(data, metrics, drug.col = "drug", prefix = NULL, log.transform.IC = TRUE, ...) {
  for(metric in metrics) {
    print(metric)
    ## Take the median of the metric across replicants
    tmp <- ddply(data[,c("patient_id", drug.col, metric)], c(drug.col, "patient_id"), 
                 .fun = function(df) { median(df[,metric]) })
    colnames(tmp) <- c(drug.col, "patient_id", metric)
    tmp <- spread_(tmp, key = "patient_id", value = metric)
    rownames(tmp) <- tmp[,drug.col]
    tmp <- tmp[,-1]
    ## Remove any outliers
    tmp <- as.matrix(tmp)
    if(grepl(metric, pattern="IC", ignore.case=TRUE)) { 
      ## Remove any IC50s that were set to the max
      for(row in 1:nrow(tmp)) {
        flag <- tmp[row,] == max(tmp[row,], na.rm=TRUE)
        tmp[row,flag] <- NA
      }
    }
    for(row in 1:nrow(tmp)) {
      flag <- !is.na(tmp[row,]) & unname(abs( (tmp[row,] - median(tmp[row,], na.rm=TRUE)) / mad(tmp[row,], na.rm=TRUE) ) > 3)
      tmp[row,flag] <- NA
    }
    if(grepl(metric, pattern="IC", ignore.case=TRUE) && log.transform.IC) { 
      tmp <- log(tmp) 
      for(row in 1:nrow(tmp)) {
        flag <- !is.finite(tmp[row,])
        tmp[row,flag] <- NA
        tmp[row,] <- scale(tmp[row,])
      }
      for(row in 1:nrow(tmp)) {
        flag <- !is.na(tmp[row,]) & ( abs(tmp[row,]) > 3 )
        tmp[row,flag] <- NA
      }
      
    }
    if(!is.null(prefix)) {
      pdf(paste0("ohsu-mek-", metric, ".pdf"))
    }
    pairs(t(tmp), ...)
    if(!is.null(prefix)) {
      d <- dev.off()
    }
  }
}

plot.LL.4.curve <- function(slope, min.asymptote, max.asymptote, ic50, x.min, x.max) {
  x.seq <- seq(from = x.min, to = x.max, by = 0.01)
  b <- slope
  c <- min.asymptote
  d <- max.asymptote
  e <- ic50
  y <- unlist(lapply(x.seq, function(x) {
    c + ( (d-c)/(1+exp(b*(log(x)-log(e)))) )
  }))
  plot(log10(x.seq), y, type="l")
}

## L.4 fits the function
## f(x) = c + \frac{d-c}{(1+\exp(b(x - e)))}
## i.e., x is on log, not real, scale
plot.L.4.curve <- function(slope, min.asymptote, max.asymptote, ic50, x.min, x.max, add = FALSE) {
  x.seq <- seq(from = x.min, to = x.max, by = 0.01)
  b <- slope
  c <- min.asymptote
  d <- max.asymptote
  e <- ic50
  y <- unlist(lapply(x.seq, function(x) {
    c + ( (d-c)/(1+exp(b*(x-e))) )
  }))
  if(add==FALSE) {
    plot(log10(x.seq), y, type="l")
  } else {
    lines(log10(x.seq), y, type="l")
  }
}

## All concentrations/IC50s are in real space
compute.ll.4.auc <- function(slope, min.asymptote, max.asymptote, ic50, min.conc, max.conc, force.ic50.in.range = TRUE) {
  ## FIMM defines logistic (not log-logistic) function as:
  ## f(x) = d + ( a - d ) * [ 1 + 10^(b * c - b * x) ]^-1
  ## with (indefinite) integral
  ## Y(x) = a * x + { ( 1 / b ) * (a - d) * log10 [ 1 + 10^(b * c - b * x) ] }
  ## 
  ## LL.4 in DRC is:
  ## f(x) = c + \frac{d-c}{1+\exp(b(\log(x)-\log(e)))}
  ## f(x) = c + ( d - c ) / ( 1 + exp(b * log(x) - b * log(f)) )
  ## The analog of Yadav Eqn 4 for this eqn is:
  ## x1 = e + (1/b) * [ log(d - y) - log(y - c) ]
  ## Begin the integration at x1(y=t)
  b <- slope
  c <- min.asymptote
  d <- max.asymptote
  e <- ic50
  fnc <- function(x) { ll.4.func(x, b, c, d, e) }

  x1 <- min.conc
  t <- 10
  x1 <- e + (1/b) * (log(d - t) - log(t - c))
  if(!is.infinite(x1)) {
    vec <- c(NA, NA, NA, NA)
    names(vec) <- c("auc", "dss1", "dss2", "dss3")
    return(vec)
  }
  val <- integrate(fnc, min.conc, max.conc)
  auc <- NA
  if(val$message == "OK") { auc <- val$value }
  if(!is.infinite(auc)) {
    vec <- c(NA, NA, NA, NA)
    names(vec) <- c("auc", "dss1", "dss2", "dss3")
    return(vec)
  }
  dss <- compute.dss(auc, ic50, t = t, c.min = min.conc, c.max = max.conc, x1, x2 = max.conc, top.asymptote = max.asymptote, force.ic50.in.range = force.ic50.in.range)
  vec <- c(auc, dss)
  names(vec) <- c("auc", names(dss))
  return(vec)
}

compute.l.4.auc <- function(slope, min.asymptote, max.asymptote, ic50, min.conc, max.conc, numerical = TRUE) {
  ## FIMM defines logistic (not log-logistic) function as:
  ## f(x) = d + ( a - d ) * [ 1 + 10^(b * c - b * x) ]^-1
  ## with (indefinite) integral
  ## Y(x) = a * x + { ( 1 / b ) * (a - d) * log10 [ 1 + 10^(b * c - b * x) ] }
  ## DRC L.4 is
  ## f(x) = c' + ( d' - c' ) * [ 1 + exp(b' * x - b' * e' ) ]^-1
  ## i.e., 
  ## a = d' = max.asymptote
  ## b = - ln(10) b' = -ln(10) slope
  ## c = e' = ic50
  ## d = c' = min.asymptote
  stop("this needs to be fixed for AUC/DSS/x1")
  auc <- 0
  b <- slope
  c <- min.asymptote
  d <- max.asymptote
  e <- ic50
  if(numerical) {  
    fnc <- function(x) { l.4.func(x, b, c, d, e) }
    val <- integrate(fnc, min.conc, max.conc)
    if(val$message == "OK") { auc <- val$value }
  } else {
    a <-  max.asymptote
    b <- - log(10) * slope
    c <- ic50
    d <- min.asymptote
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
      ## res <- a * x + ( ( 1 / b ) * (a - d) * log10(1 + 10^(b * c - b * x) ) )
      res
    }
    auc <- y.int(max.conc, a, b, c, d) - y.int(min.conc, a, b, c, d)
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
## Yadav sets DSS = 0when the ic50 is at or beyond the max dose level tested c.max.
compute.dss <- function(auc, ic50, t = 10, c.min, c.max, x1, x2, top.asymptote, force.ic50.in.range = TRUE) {
  if(c.max < c.min) { warning("c.max should be > c.min\n") }
  if(x2 < x1) { warning("x2 should be > x1\n") }
  ## Distinguish between cases in which it makes sense for DSS to be zero and when
  ## it is an error condition and should be NA.
  if(auc == 0) {
    ret <- c(0, 0, 0)
    names(ret) <- c("dss1", "dss2", "dss3")
    return(ret)
  }
  if(any(is.na(c(auc, ic50, t, c.min, c.max, x1, x2, top.asymptote)))) {
    ret <- c(NA, NA, NA)
    names(ret) <- c("dss1", "dss2", "dss3")
    return(ret)
  }
  dss1 <- auc - t * (x2 - x1)
  dss1 <- dss1 / ( (100 - t) * (c.max - c.min) )
  if(force.ic50.in.range && (ic50 >= c.max)) {
    ## Set DSS1 to 0 if the ic50 is at or exceeds the max does level tested.
    ## Since DSS2 and DSS3 are multiplies of DSS1, this will effectively set 
    ## DSS2 and DSS3 to zero as well (below).
    dss1 <- 0
  }
  ## This catches dividing by zero (above)
  if(c.max == c.min) { dss1 <- NA }
  dss2 <- dss1 / log10(top.asymptote)
  ## This catches taking the log of zero/negative.
  if(top.asymptote <= 0) { dss2 <- NA }
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

fit.fct.to.ohsu <- function(data, fct = "LL.4") {
  ddply(data[,c("inhibitor", "well_concentration", "normalized_viability", "time_of_read", "lab_id", "replicant", "patient_id")],
        c("patient_id", "inhibitor", "lab_id", "replicant", "time_of_read"),
        .parallel = TRUE,
        .fun = function(df) {
          fit <- NULL
          if(fct == "LL.4") {
            tryCatch({
              fit <- suppressWarnings(drm(normalized_viability ~ well_concentration, data = df, fct=LL.4(fixed = c(NA, NA, NA,NA), names = c("SLOPE","MIN","MAX","IC50"))))
            }, 
            error = function(e) {
              cat("Trapped error in LL.4\n")
              NULL
            })
          } else {
            tryCatch({
              fit <- suppressWarnings(drm(normalized_viability ~ well_concentration, data = df, fct=L.4(fixed = c(NA, NA, NA,NA), names = c("SLOPE","MIN","MAX","IC50"))))
            }, 
            error = function(e) {
              cat("Trapped error in L.4\n")
              NULL
            })
          }
##          fit
          gof <- 0
          ic50 <- 0
          auc <- 0
          slope <- 0
          min.asymptote <- 0
          max.asymptote <- 0
          min.conc <- min(df$well_concentration)
          max.conc <- max(df$well_concentration)
          if(!is.null(fit) && (fit$fit$convergence || fit$convergence)) {
            ssres <- sum(residuals(fit)^2)
            resp <- na.omit(fit$dataList$resp)
            sstot <- sum((mean(resp)-resp)^2)
            gof <- 1 - (ssres/sstot)
            slope <- coef(fit)[1]
            min.asymptote <- coef(fit)[2]
            max.asymptote <- coef(fit)[3]
            ic50 <- coef(fit)[4]
            if(fct == "LL.4") {
              auc <- compute.ll.4.auc(slope, min.asymptote, max.asymptote, ic50, min.conc, max.conc)
            } else {
              auc <- compute.l.4.auc(slope, min.asymptote, max.asymptote, ic50, min.conc, max.conc, numerical = FALSE) 
            }
          }
          vec <- c(slope, min.asymptote, max.asymptote, ic50, gof, auc, min.conc, max.conc)
          names(vec) <- c("slope", "min.asymptote", "max.asymptote", "ic50", "gof", "auc", "min.conc", "max.conc")
          vec
        })
  }

fit.fct.to.fimm <- function(data, fct = "LL.4") {
  ddply(data[,c("DRUG_ID", "DRUG_NAME", "CONCENTRATION", "PERCENT_INHIBITION", "SCREEN_ID", "SAMPLE_DATE", "PATIENT_ID")],
        c("SCREEN_ID", "DRUG_ID"),
        .parallel = TRUE,
        .fun = function(df) {
          fit <- NULL
          if(fct == "LL.4") {
            tryCatch({
              fit <- suppressWarnings(drm(PERCENT_INHIBITION ~ CONCENTRATION, data = df, fct=LL.4(fixed = c(NA, NA, NA,NA), names = c("SLOPE","MIN","MAX","IC50"))))
            }, 
            error = function(e) {
              cat("Trapped error in LL.4\n")
              NULL
            })
          } else {
            tryCatch({
              fit <- suppressWarnings(drm(PERCENT_INHIBITION ~ CONCENTRATION, data = df, fct=L.4(fixed = c(NA, NA, NA,NA), names = c("SLOPE","MIN","MAX","IC50"))))
            }, 
            error = function(e) {
              cat("Trapped error in L.4\n")
              NULL
            })
          }
          ##          fit
          gof <- 0
          ic50 <- 0
          auc <- 0
          slope <- 0
          min.asymptote <- 0
          max.asymptote <- 0
          min.conc <- min(df$CONCENTRATION)
          max.conc <- max(df$CONCENTRATION)
          if(!is.null(fit) && (fit$fit$convergence || fit$convergence)) {
            ssres <- sum(residuals(fit)^2)
            resp <- na.omit(fit$dataList$resp)
            sstot <- sum((mean(resp)-resp)^2)
            gof <- 1 - (ssres/sstot)
            slope <- coef(fit)[1]
            min.asymptote <- coef(fit)[2]
            max.asymptote <- coef(fit)[3]
            ic50 <- coef(fit)[4]
            if(fct == "LL.4") {
              auc <- compute.ll.4.auc(slope, min.asymptote, max.asymptote, ic50, min.conc, max.conc)
            } else {
              auc <- compute.l.4.auc(slope, min.asymptote, max.asymptote, ic50, min.conc, max.conc, numerical = FALSE) 
            }
          }
          vec <- c(slope, min.asymptote, max.asymptote, ic50, gof, auc, min.conc, max.conc)
          names(vec) <- c("slope", "min.asymptote", "max.asymptote", "ic50", "gof", "auc", "min.conc", "max.conc")
          vec
        })
}


cat("Fitting LL4\n")
ll4.fits <- fit.fct.to.ohsu(ohsu.raw.drug.data, fct = "LL.4")
cat("Done fitting LL4\n")
save.image(".RData")
cat("Done saving .RData\n")

ll4.dss.fits.dss.fixed <- compute.ohsu.dss(ll4.fits, force.ic50.in.range=TRUE)
ll4.dss.fits.no.force.dss.fixed <- compute.ohsu.dss(ll4.fits, force.ic50.in.range=FALSE)

remove.outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

plot.ohsu.panel.of.pairs <- function(data, main = NULL) {
  gs <- list()
  for(metric in c("ic50", "auc", "dss1", "dss2")) {
    tmp <- data[,c("lab_id", "inhibitor", metric, "replicant", "time_of_read")]
    if(metric == "ic50") {
      tmp[,metric] <- log(tmp[,metric])
    }
    flag <- is.finite(tmp[,metric])
    tmp <- tmp[flag,]
    dat <- spread_(data=tmp, key="inhibitor", value=metric)
    cols <- unique(tmp$inhibitor)
    df <- dat[,cols]
    for(col in 1:ncol(df)) {
      df[,col] <- remove.outliers(df[,col])
    }
    
    colnames(df) <- make.names(colnames(df))
    corr.mat1<-tryCatch({cor.matrix(variables=do.call("ddf", as.list(df)),,
                          dat=df,
                          test=cor.test,
                          method='spearman',
                          alternative="two.sided",exact=FALSE)},
                        error = function(e) { 
                          cat("Caught error in cor.matrix")
                          return(NULL) 
                        })
    if(is.null(corr.mat1)) { next }
    
    p <- ggcorplot(corr.mat1,data = df, var_text_size = 3)
    if(is.null(main)) {
      p <- p + ggtitle(metric)
    } else {
      p <- p + ggtitle(paste0(main, ":\n", metric))
    }
    gs[[metric]] <- p
  }
  ## do.call("grid.arrange", gs, ncol=2)
  return(do.call("arrangeGrob", list(grobs = gs, ncol = 2)))
}

plot.fimm.panel.of.pairs <- function(data, main = NULL) {
  gs <- list()
  for(metric in c("ic50", "auc", "dss1", "dss2")) {
    tmp <- data[,c("SCREEN_ID", "DRUG_ID", metric)]
    if(metric == "ic50") {
      tmp[,metric] <- log(tmp[,metric])
    }
    flag <- is.finite(tmp[,metric])
    tmp <- tmp[flag,]
    dat <- spread_(data=tmp, key="DRUG_ID", value=metric)
    cols <- unique(tmp$DRUG_ID)
    df <- dat[,cols]
    for(col in 1:ncol(df)) {
      df[,col] <- remove.outliers(df[,col])
    }
##    print(head(df))
##    print(length(which(complete.cases(df))))
##    if(length(which(complete.cases(df))) < 5) { return(NULL) }
    colnames(df) <- make.names(colnames(df))
    corr.mat1<-tryCatch({cor.matrix(variables=do.call("ddf", as.list(df)),,
                          dat=df,
                          test=cor.test,
                          method='spearman',
                          alternative="two.sided",exact=FALSE)},
                        error = function(e) {
                          cat("Caught error in cor.matrix\n")
                          return(NULL)
                        })
    if(is.null(corr.mat1)) { next }
    p <- ggcorplot(corr.mat1,data = df, var_text_size = 3)
    if(is.null(main)) {
      p <- p + ggtitle(metric)
    } else {
      p <- p + ggtitle(paste0(main, ":\n", metric))
    }
    gs[[metric]] <- p
  }
  ## do.call("grid.arrange", gs)
  return(do.call("arrangeGrob", list(grobs = gs, ncol = 2)))
}

mechanisms <- as.data.frame(table(ohsu.fimm.drugs$`Mechanism/Targets`))
colnames(mechanisms)[1] <- "mechanism"
mechanisms <- mechanisms$mechanism[mechanisms$Freq > 1]
for(mechanism in mechanisms) {
  print(mechanism)
  mechanism.annotations <- ohsu.fimm.drugs[ohsu.fimm.drugs$`Mechanism/Targets` == mechanism,]

  ohsu.mechanism.data <- subset(ll4.dss.fits, inhibitor %in% mechanism.annotations$ID_Drug.ohsu)
  fname <- paste0(make.names(mechanism), "-correlation-ohsu-restricted-ic50.pdf")
  pdf(fname, onefile = FALSE)
  plt <- plot.ohsu.panel.of.pairs(ohsu.mechanism.data, main = paste0("OHSU: ", mechanism))
  ## ggsave(fname, plt)
  plot(plt)
  d <- dev.off()

  ohsu.mechanism.data <- subset(ll4.dss.fits.no.force, inhibitor %in% mechanism.annotations$ID_Drug.ohsu)
  fname <- paste0(make.names(mechanism), "-correlation-ohsu-unrestricted-ic50.pdf")
  pdf(fname, onefile = FALSE)
  plt <- plot.ohsu.panel.of.pairs(ohsu.mechanism.data, main = paste0("OHSU ", mechanism, " (unrestricted IC50)"))
  ## ggsave(fname, plt)
  plot(plt)
  d <- dev.off()
}

cat("Fitting FIMM LL4\n")
fimm.ll4.fits <- fit.fct.to.fimm(fimm.raw.drug.data, fct = "LL.4")
cat("Done fitting FIMM LL4\n")
save.image(".RData")

fimm.ll4.dss.fits <- compute.fimm.dss(fimm.ll4.fits, force.ic50.in.range=TRUE)
fimm.ll4.dss.fits.no.force <- compute.fimm.dss(fimm.ll4.fits, force.ic50.in.range=FALSE)

for(mechanism in mechanisms) {
  print(mechanism)
  mechanism.annotations <- ohsu.fimm.drugs[ohsu.fimm.drugs$`Mechanism/Targets` == mechanism,]
  
  fimm.mechanism.data <- subset(fimm.ll4.dss.fits, DRUG_ID %in% mechanism.annotations$ID_Drug)
  fname <- paste0(make.names(mechanism), "-correlation-fimm-restricted-ic50.pdf")
  pdf(fname, onefile=FALSE)
  plt <- plot.fimm.panel.of.pairs(fimm.mechanism.data, main = paste0("FIMM ", mechanism))
  plot(plt)
  d <- dev.off()
  
  fimm.mechanism.data <- subset(fimm.ll4.dss.fits.no.force, DRUG_ID %in% mechanism.annotations$ID_Drug)
  fname <- paste0(make.names(mechanism), "-correlation-fimm-unrestricted-ic50.pdf")
  pdf(fname, onefile=FALSE)
  plt <- plot.fimm.panel.of.pairs(fimm.mechanism.data, main = paste0("FIMM ", mechanism, " (unrestricted IC50)"))
  plot(plt)
  d <- dev.off()
}

## Univariate analysis of MEK1/2 inhibitors

## Add the translation from drug response labId to sequencing ID to the ohsu data
ll4.dss.fits.no.force.seq <- merge(ll4.dss.fits.no.force, unique(rnaseq.sample.summary.tbl[,c("SeqID", "Original_LabID")]), by.x = "lab_id", by.y = "Original_LabID")
ll4.dss.fits.no.force.seq <- subset(ll4.dss.fits.no.force.seq, gof > 0)

ll4.dss.fits.no.force.seq.dss.fixed <- merge(ll4.dss.fits.no.force.dss.fixed, unique(rnaseq.sample.summary.tbl[,c("SeqID", "Original_LabID")]), by.x = "lab_id", by.y = "Original_LabID")
ll4.dss.fits.no.force.seq.dss.fixed <- subset(ll4.dss.fits.no.force.seq.dss.fixed, gof > 0)

## Do each inhibitor separately.

common.genes <- intersect(rownames(beat.mafs), rownames(fimm.genomic.bin))
common.genes <- intersect(rownames(fimm.expr.data), rownames(ohsu.expr.data))
ohsu.mek.expr.auc.corr <- correlate.expression.with.drug.response(drug.df = ll4.dss.fits.no.force.seq, drug.id.col = "inhibitor", genes.to.test = common.genes, drugs.to.test = mek.annotations$ID_Drug.ohsu, expr.mat = ohsu.expr.data, specimen.id.col = "SeqID", patient.id.col = "patient_id", metric.col = "auc")
cat("Done correlating expression with AUC for OHSU\n")
save.image(".RData")

screen.pt.tbl <- unique(fimm.raw.drug.data[,c("PATIENT_ID", "SCREEN_ID")])
fimm.ll4.dss.fits.no.force.pt <- merge(fimm.ll4.dss.fits.no.force, screen.pt.tbl)
fimm.ll4.dss.fits.no.force.pt <- subset(fimm.ll4.dss.fits.no.force.pt, gof > 0)
  
fimm.mek.expr.auc.corr <- correlate.expression.with.drug.response(drug.df = fimm.ll4.dss.fits.no.force.pt, drug.id.col = "DRUG_ID", genes.to.test = common.genes, drugs.to.test = mek.annotations$ID_Drug, expr.mat = fimm.expr.data, specimen.id.col = "SCREEN_ID", patient.id.col = "PATIENT_ID", metric.col = "auc")
cat("Done correlating expression with AUC for FIMM\n")
plot(hist(as.numeric(na.omit(fimm.mek.expr.auc.corr$pval))))

ohsu.mek.expr.all.auc.corr <- correlate.expression.with.all.drug.response(drug.df = ll4.dss.fits.no.force.seq, drug.id.col = "inhibitor", genes.to.test = common.genes, drugs.to.test = mek.annotations$ID_Drug.ohsu, expr.mat = ohsu.expr.data, specimen.id.col = "SeqID", patient.id.col = "patient_id", metric.col = "auc")
fimm.mek.expr.all.auc.corr <- correlate.expression.with.all.drug.response(drug.df = fimm.ll4.dss.fits.no.force.pt, drug.id.col = "DRUG_ID", genes.to.test = common.genes, drugs.to.test = mek.annotations$ID_Drug, expr.mat = fimm.expr.data, specimen.id.col = "SCREEN_ID", patient.id.col = "PATIENT_ID", metric.col = "auc")
pdf("fimm-ohsu-auc-expr-gene-drug-pval.pdf")
par(mfrow=c(1,2))
hist(as.numeric(na.omit(fimm.mek.expr.all.auc.corr$pval)), main="FIMM MEK:\nAUC ~ gene + drug_id", xlab = "p-value")
hist(as.numeric(na.omit(ohsu.mek.expr.all.auc.corr$pval)), main="OHSU MEK:\nAUC ~ gene + drug_id", xlab = "p-value")
d <- dev.off()

ohsu.mek.expr.all.dss1.corr <- correlate.expression.with.all.drug.response(drug.df = ll4.dss.fits.no.force.seq.dss.fixed, drug.id.col = "inhibitor", genes.to.test = common.genes, drugs.to.test = mek.annotations$ID_Drug.ohsu, expr.mat = ohsu.expr.data, specimen.id.col = "SeqID", patient.id.col = "patient_id", metric.col = "dss1")
fimm.mek.expr.all.dss1.corr <- correlate.expression.with.all.drug.response(drug.df = fimm.ll4.dss.fits.no.force.pt, drug.id.col = "DRUG_ID", genes.to.test = common.genes, drugs.to.test = mek.annotations$ID_Drug, expr.mat = fimm.expr.data, specimen.id.col = "SCREEN_ID", patient.id.col = "PATIENT_ID", metric.col = "dss1")
pdf("fimm-ohsu-dss1-expr-gene-drug-pval.pdf")
par(mfrow=c(1,2))
hist(as.numeric(na.omit(fimm.mek.expr.all.dss1.corr$pval)), main="FIMM MEK:\nDSS1 ~ gene + drug_id", xlab = "p-value")
hist(as.numeric(na.omit(ohsu.mek.expr.all.dss1.corr$pval)), main="OHSU MEK:\nDSS1 ~ gene + drug_id", xlab = "p-value")
d <- dev.off()

ohsu.mek.expr.all.ic50.corr <- correlate.expression.with.all.drug.response(drug.df = ll4.dss.fits.no.force.seq, drug.id.col = "inhibitor", genes.to.test = common.genes, drugs.to.test = mek.annotations$ID_Drug.ohsu, expr.mat = ohsu.expr.data, specimen.id.col = "SeqID", patient.id.col = "patient_id", metric.col = "ic50", log.transform = TRUE)
fimm.mek.expr.all.ic50.corr <- correlate.expression.with.all.drug.response(drug.df = fimm.ll4.dss.fits.no.force.pt, drug.id.col = "DRUG_ID", genes.to.test = common.genes, drugs.to.test = mek.annotations$ID_Drug, expr.mat = fimm.expr.data, specimen.id.col = "SCREEN_ID", patient.id.col = "PATIENT_ID", metric.col = "ic50", log.transform = TRUE)
pdf("fimm-ohsu-ic50-expr-gene-drug-pval.pdf")
par(mfrow=c(1,2))
hist(as.numeric(na.omit(fimm.mek.expr.all.ic50.corr$pval)), main="FIMM MEK:\nlog IC50 ~ gene + drug_id", xlab = "p-value")
hist(as.numeric(na.omit(ohsu.mek.expr.all.ic50.corr$pval)), main="OHSU MEK:\nlog IC50 ~ gene + drug_id", xlab = "p-value")
d <- dev.off()

my.plot.venn <- function(df, pval.threshold, set1.name, set2.name, pval.col1, pval.col2, gene.col1, gene.col2, title=NULL) {
  sig.genes1 <- df[!is.na(df[,pval.col1]) & (as.numeric(df[,pval.col1]) < pval.threshold),gene.col1]
  sig.genes2 <- df[!is.na(df[,pval.col2]) & (as.numeric(df[,pval.col2]) < pval.threshold),gene.col2]
  vennList <- list("set1"=sig.genes1, "set2"=sig.genes2)
  names(vennList) <- c(set1.name, set2.name)

  num.sig1 <- length(sig.genes1)
  num.sig2 <- length(sig.genes2)
  num.sig.intersect <- length(intersect(sig.genes1, sig.genes2))
  do.sim <- FALSE
  pval <- 0
  if(do.sim == FALSE) {
    x1 <- !is.na(df[,pval.col1]) & (as.numeric(df[,pval.col1]) < pval.threshold)
    x2 <- !is.na(df[,pval.col2]) & (as.numeric(df[,pval.col2]) < pval.threshold)
    ft <- fisher.test(table(x1, x2))
    pval <- ft$p.val
  } else {
    num.iters <- 1000
  
    res.pval <- sum(unlist(llply(1:num.iters, .parallel = FALSE, .fun = function(i) {
      indices1 <- sample.int(size = num.sig1, n = nrow(df), replace=FALSE)
      indices2 <- sample.int(size = num.sig2, n = nrow(df), replace=FALSE)
      df1 <- df[indices1,]
      df2 <- df[indices2,]
      tmp.sig.genes1 <- df1[!is.na(df1[,pval.col1]) & (as.numeric(df1[,pval.col1]) < pval.threshold),gene.col1]
      tmp.sig.genes2 <- df2[!is.na(df2[,pval.col2]) & (as.numeric(df2[,pval.col2]) < pval.threshold),gene.col2]
      tmp.num.sig.intersect <- length(intersect(tmp.sig.genes1, tmp.sig.genes2))
      if(tmp.num.sig.intersect >= num.sig.intersect) {
        return(1)
      } else {
        return(0)
      }
    })))
    pval <- res.pval/num.iters
    if(res.pval == 0) { pval <- 1/num.iters }
  }
  pval <- signif(pval, digits=2)
  cat(paste0("# sig (p-val < ", pval.threshold, ") in ", set1.name, ": ", num.sig1, "\n"))
  cat(paste0("# sig (p-val < ", pval.threshold, ") in ", set2.name, ": ", num.sig2, "\n"))
  cat(paste0("# sig (p-val < ", pval.threshold, ") in both: ", num.sig.intersect, "\n"))
  cat(paste0("p-val of this intersection:", ifelse(res.pval == 0, " <", " "), pval, "\n"))
  pval.str <- paste0("intersection pval:", ifelse(res.pval == 0, " < ", " "), pval)
  main <- pval.str
  if(!is.null(title)) {
    main <- paste0(title, "\n", pval.str)
  }
  
  ## plot(venn(vennList), main=main)
  g <- venn.diagram(vennList, filename=NULL, main=main)
  g
}

m <- merge(fimm.mek.expr.all.auc.corr, ohsu.mek.expr.all.auc.corr, by = "gene", suffixes = c(".fimm", ".ohsu"))
pdf("fimm-ohsu-auc-expr-gene-drug-venn.pdf")
g <- my.plot.venn(m, 0.05, "FIMM", "OHSU", "pval.fimm", "pval.ohsu", "gene", "gene", title = "significant genes (p < 0.05)\nauc ~ expr + drug")
grid.draw(g)
d <- dev.off()

## NB: in.ras.sig == 0 means it is in RAS signature.
## This allows us to use the venn diagram code above.
m$in.ras.sig <- 1
flg <- m$gene %in% ras.sig.genes$ENSG
m$in.ras.sig[flg] <- 0
g1 <- my.plot.venn(m, 0.05, "FIMM", "RAS", "pval.fimm", "in.ras.sig", "gene", "gene", title = "significant genes (p < 0.05)\nMEK auc ~ expr + drug")
g2 <- my.plot.venn(m, 0.05, "OHSU", "RAS", "pval.ohsu", "in.ras.sig", "gene", "gene", title = "significant genes (p < 0.05)\nMEK auc ~ expr + drug")
pdf("fimm-ohsu-auc-expr-gene-drug-ras-venn.pdf")
grid.arrange(gTree(children=g1),gTree(children=g2), ncol=2)
d <- dev.off()

m <- merge(fimm.mek.expr.all.dss1.corr, ohsu.mek.expr.all.dss1.corr, by = "gene", suffixes = c(".fimm", ".ohsu"))
pdf("fimm-ohsu-dss1-expr-gene-drug-venn.pdf")
g <- my.plot.venn(m, 0.05, "FIMM", "OHSU", "pval.fimm", "pval.ohsu", "gene", "gene", title = "significant genes (p < 0.05)\nDSS1 ~ expr + drug")
grid.draw(g)
d <- dev.off()

## NB: in.ras.sig == 0 means it is in RAS signature.
## This allows us to use the venn diagram code above.
m$in.ras.sig <- 1
flg <- m$gene %in% ras.sig.genes$ENSG
m$in.ras.sig[flg] <- 0
g1 <- my.plot.venn(m, 0.05, "FIMM", "RAS", "pval.fimm", "in.ras.sig", "gene", "gene", title = "significant genes (p < 0.05)\nMEK auc ~ expr + drug")
g2 <- my.plot.venn(m, 0.05, "OHSU", "RAS", "pval.ohsu", "in.ras.sig", "gene", "gene", title = "significant genes (p < 0.05)\nMEK auc ~ expr + drug")
pdf("fimm-ohsu-dss1-expr-gene-drug-ras-venn.pdf")
grid.arrange(gTree(children=g1),gTree(children=g2), ncol=2)
d <- dev.off()

## Plot a few of the genes that are in RAS signature
gene <- "ENSG00000104972"
gene <- "ENSG00000145416"
gene <- ras.sig.genes$ENSG[5]
gene <- "ENSG00000141505"
gene <- "ENSG00000139318"
gene <- "ENSG00000197956"

ras.sig <- subset(m, (in.ras.sig == 0) & (as.numeric(pval.fimm) < 0.05) & (as.numeric(pval.ohsu) < 0.05))
## ras.sig <- subset(m, (in.ras.sig == 0) & (as.numeric(pval.ohsu) < 0.05))
ras.sig <- ras.sig[order(as.numeric(ras.sig$pval.ohsu)),]
head(ras.sig[order(as.numeric(ras.sig$pval.ohsu)),c("gene","pval.ohsu","r2.ohsu")])
for(i in 1:5) {
  gene <- ras.sig$gene[i]
  pdf(paste0("fimm-ohsu-ras-", gene, "-vs-dss1.pdf"))
  ##par(mfrow=c(1,2))
  t.gene <- data.frame(gene = unname(t(ohsu.expr.data[gene,])), seqid = colnames(ohsu.expr.data))
  d.tmp <- subset(ll4.dss.fits.no.force.seq.dss.fixed, inhibitor %in% mek.annotations$ID_Drug.ohsu)
  m.tmp <- merge(d.tmp, t.gene, by.x = "SeqID", by.y = "seqid")
  m.tmp$auc2 <- ( 100 * (m.tmp$max.conc - m.tmp$min.conc ) ) - m.tmp$auc
  fit <- lm(m.tmp$dss1 ~ m.tmp$gene + m.tmp$inhibitor)
  rmse <- round(sqrt(mean(resid(fit)^2)), 2)
  r2 <- round(summary(fit)$r.squared, 2)
  eqn <- bquote(r^2 == .(r2) * "," ~~ RMSE == .(rmse))
  plot(m.tmp$gene, m.tmp$dss1, xlab = paste0(gene, " Expr"), ylab = "AUC", main = "OHSU")
  abline(fit)
  text(ceiling(min(m.tmp$gene)), max(m.tmp$auc) - 100, eqn, pos = 4)
  
  t.gene <- data.frame(gene = unname(t(fimm.expr.data[gene,,drop=F])), SCREEN_ID = colnames(fimm.expr.data))
  d.tmp <- subset(fimm.ll4.dss.fits.no.force.pt, DRUG_ID %in% mek.annotations$ID_Drug)
  m.tmp <- merge(d.tmp, t.gene, by = "SCREEN_ID")
  m.tmp$DRUG_ID <- as.factor(m.tmp$DRUG_ID)
  fit <- lm(formula = dss1 ~ gene + DRUG_ID, data = m.tmp)
  plot(m.tmp$dss1 ~ m.tmp$gene)
  for(drug.id in levels(m.tmp$DRUG_ID)) {
    new.x <- unique(m.tmp$gene)
    pred <- predict(fit, newdata = data.frame(gene = new.x, DRUG_ID = drug.id))
    pred <- data.frame(gene = new.x, dss1 = unname(pred))
    pred <- pred[order(pred$gene),]
    lines(pred$gene, pred$dss1)
  }
  rmse <- round(sqrt(mean(resid(fit)^2)), 2)
  r2 <- round(summary(fit)$r.squared, 2)
  eqn <- bquote(r^2 == .(r2) * "," ~~ RMSE == .(rmse))
  plot(m.tmp$gene, m.tmp$dss1, xlab = paste0(gene, " Expr"), ylab = "AUC", main = "FIMM")
  abline(fit)
  text(ceiling(min(m.tmp$gene)), max(m.tmp$auc) - 100, eqn, pos = 4)
  d <- dev.off()
}

ras.sig <- subset(m, (in.ras.sig == 0) & (as.numeric(pval.fimm) < 0.05) & (as.numeric(pval.ohsu) < 0.05))
## ras.sig <- subset(m, (in.ras.sig == 0) & (as.numeric(pval.ohsu) < 0.05))
ras.sig <- ras.sig[order(as.numeric(ras.sig$pval.ohsu)),]
for(i in 1:5) {
  gene <- ras.sig$gene[i]
  sym <- ensg.to.symbols.mapping(gene)$SYMBOL[1]
  pdf(paste0("fimm-ohsu-ras-", sym, "-vs-dss1.pdf"))
  ##par(mfrow=c(1,2))
  glist1 <- list()
  glist2 <- list()
  
  t.gene <- data.frame(gene = unname(t(ohsu.expr.data[gene,])), seqid = colnames(ohsu.expr.data))
  d.tmp <- subset(ll4.dss.fits.no.force.seq.dss.fixed, inhibitor %in% mek.annotations$ID_Drug.ohsu)
  m.tmp <- merge(d.tmp, t.gene, by.x = "SeqID", by.y = "seqid")
  m.tmp$auc2 <- ( 100 * (m.tmp$max.conc - m.tmp$min.conc ) ) - m.tmp$auc
  drug.col <- "inhibitor"
  df <- m.tmp[,c("gene", "dss1", drug.col)]
  colnames(df) <- c("gene", "metric", "drug")
  df$drug <- factor(df$drug)
  fit <- lm(data = df, metric ~ gene + drug)
  for(drug.id in unique(df$drug)) {
    flag <- df$drug == drug.id
    g <- ggplot(data = df[flag,], aes(x = gene, y = metric))
    g <- g + geom_point()
    g <- g + ggtitle(paste0("OHSU: ", drug.id))
    new.x <- unique(df$gene)
    pred <- predict(fit, newdata = data.frame(gene = new.x, drug = drug.id))
    pred <- data.frame(gene = new.x, metric = pred)
    g <- g + geom_line(data = pred, aes(x = gene, y = metric))
    g <- g + xlab(sym)
    g <- g + ylab(metric)
    rmse <- round(sqrt(mean(resid(fit)^2)), 2)
    r2 <- round(summary(fit)$r.squared, 2)
    eqn <- bquote(r^2 == .(r2) * "," ~~ RMSE == .(rmse))
    eqn <- r2
    eqn <- substitute(italic(r)^2~"="~r2, 
                     list(r2 = format(summary(fit)$r.squared, digits = 3)))
    eqn <- as.character(as.expression(eqn));
    d.eqn <- data.frame(x = min(df$gene,na.rm=TRUE) + 0.1 * (max(df$gene, na.rm=TRUE) - min(df$gene, na.rm=TRUE)), y = 0.8, label = eqn)
    g <- g + geom_text(data = d.eqn, aes(x = x, y = y, label = label), parse=TRUE)
    glist1[[length(glist1)+1]] <- g
  }
  
  t.gene <- data.frame(gene = unname(t(fimm.expr.data[gene,,drop=F])), SCREEN_ID = colnames(fimm.expr.data))
  d.tmp <- subset(fimm.ll4.dss.fits.no.force.pt, DRUG_ID %in% mek.annotations$ID_Drug)
  m.tmp <- merge(d.tmp, t.gene, by = "SCREEN_ID")
  drug.col <- "DRUG_ID"
  df <- m.tmp[,c("gene", "dss1", drug.col)]
  colnames(df) <- c("gene", "metric", "drug")
  df$drug <- factor(df$drug)
  fit <- lm(data = df, metric ~ gene + drug)
  for(drug.id in unique(df$drug)) {
    flag <- df$drug == drug.id
    g <- ggplot(data = df[flag,], aes(x = gene, y = metric))
    g <- g + geom_point()
    g <- g + ggtitle(paste0("FIMM: ", drug.id))
    new.x <- unique(m.tmp$gene)
    pred <- predict(fit, newdata = data.frame(gene = new.x, drug = drug.id))
    pred <- data.frame(gene = new.x, metric = pred)
    g <- g + geom_line(data = pred, aes(x = gene, y = metric))
    g <- g + xlab(sym)
    g <- g + ylab(metric)
    rmse <- round(sqrt(mean(resid(fit)^2)), 2)
    r2 <- round(summary(fit)$r.squared, 2)
    eqn <- bquote(r^2 == .(r2) * "," ~~ RMSE == .(rmse))
    eqn <- r2
    eqn <- substitute(italic(r)^2~"="~r2, 
                      list(r2 = format(summary(fit)$r.squared, digits = 3)))
    eqn <- as.character(as.expression(eqn));
    d.eqn <- data.frame(x = min(df$gene,na.rm=TRUE) + 0.1 * (max(df$gene, na.rm=TRUE) - min(df$gene, na.rm=TRUE)), y = 0.8, label = eqn)
    g <- g + geom_text(data = d.eqn, aes(x = x, y = y, label = label), parse=TRUE)
    glist2[[length(glist2)+1]] <- g
  }

  grid.arrange(grobs = c(glist1, glist2), ncol=2, as.table=FALSE)
  ##do.call("grid.arrange", glist)
  
  d <- dev.off()
}

## Repeat above, but for those not in RAS-signature
non.ras.sig <- subset(m, (in.ras.sig == 1) & (as.numeric(pval.fimm) < 0.05) & (as.numeric(pval.ohsu) < 0.05))
non.ras.sig <- non.ras.sig[order(as.numeric(non.ras.sig$pval.ohsu)),]
for(i in 1:5) {
  gene <- non.ras.sig$gene[i]
  sym <- ensg.to.symbols.mapping(gene)$SYMBOL[1]
  pdf(paste0("fimm-ohsu-non-ras-", sym, "-vs-dss1.pdf"))
  ##par(mfrow=c(1,2))
  glist1 <- list()
  glist2 <- list()
  
  t.gene <- data.frame(gene = unname(t(ohsu.expr.data[gene,])), seqid = colnames(ohsu.expr.data))
  d.tmp <- subset(ll4.dss.fits.no.force.seq.dss.fixed, inhibitor %in% mek.annotations$ID_Drug.ohsu)
  m.tmp <- merge(d.tmp, t.gene, by.x = "SeqID", by.y = "seqid")
  m.tmp$auc2 <- ( 100 * (m.tmp$max.conc - m.tmp$min.conc ) ) - m.tmp$auc
  drug.col <- "inhibitor"
  df <- m.tmp[,c("gene", "dss1", drug.col)]
  colnames(df) <- c("gene", "metric", "drug")
  df$drug <- factor(df$drug)
  fit <- lm(data = df, metric ~ gene + drug)
  for(drug.id in unique(df$drug)) {
    flag <- df$drug == drug.id
    g <- ggplot(data = df[flag,], aes(x = gene, y = metric))
    g <- g + geom_point()
    g <- g + ggtitle(paste0("OHSU: ", drug.id))
    new.x <- unique(df$gene)
    pred <- predict(fit, newdata = data.frame(gene = new.x, drug = drug.id))
    pred <- data.frame(gene = new.x, metric = pred)
    g <- g + geom_line(data = pred, aes(x = gene, y = metric))
    g <- g + xlab(sym)
    g <- g + ylab(metric)
    rmse <- round(sqrt(mean(resid(fit)^2)), 2)
    r2 <- round(summary(fit)$r.squared, 2)
    eqn <- bquote(r^2 == .(r2) * "," ~~ RMSE == .(rmse))
    eqn <- r2
    eqn <- substitute(italic(r)^2~"="~r2, 
                      list(r2 = format(summary(fit)$r.squared, digits = 3)))
    eqn <- as.character(as.expression(eqn));
    d.eqn <- data.frame(x = min(df$gene,na.rm=TRUE) + 0.1 * (max(df$gene, na.rm=TRUE) - min(df$gene, na.rm=TRUE)), y = 0.8, label = eqn)
    g <- g + geom_text(data = d.eqn, aes(x = x, y = y, label = label), parse=TRUE)
    glist1[[length(glist1)+1]] <- g
  }
  
  t.gene <- data.frame(gene = unname(t(fimm.expr.data[gene,,drop=F])), SCREEN_ID = colnames(fimm.expr.data))
  d.tmp <- subset(fimm.ll4.dss.fits.no.force.pt, DRUG_ID %in% mek.annotations$ID_Drug)
  m.tmp <- merge(d.tmp, t.gene, by = "SCREEN_ID")
  drug.col <- "DRUG_ID"
  df <- m.tmp[,c("gene", "dss1", drug.col)]
  colnames(df) <- c("gene", "metric", "drug")
  df$drug <- factor(df$drug)
  fit <- lm(data = df, metric ~ gene + drug)
  for(drug.id in unique(df$drug)) {
    flag <- df$drug == drug.id
    g <- ggplot(data = df[flag,], aes(x = gene, y = metric))
    g <- g + geom_point()
    g <- g + ggtitle(paste0("FIMM: ", drug.id))
    new.x <- unique(m.tmp$gene)
    pred <- predict(fit, newdata = data.frame(gene = new.x, drug = drug.id))
    pred <- data.frame(gene = new.x, metric = pred)
    g <- g + geom_line(data = pred, aes(x = gene, y = metric))
    g <- g + xlab(sym)
    g <- g + ylab(metric)
    rmse <- round(sqrt(mean(resid(fit)^2)), 2)
    r2 <- round(summary(fit)$r.squared, 2)
    eqn <- bquote(r^2 == .(r2) * "," ~~ RMSE == .(rmse))
    eqn <- r2
    eqn <- substitute(italic(r)^2~"="~r2, 
                      list(r2 = format(summary(fit)$r.squared, digits = 3)))
    eqn <- as.character(as.expression(eqn));
    d.eqn <- data.frame(x = min(df$gene,na.rm=TRUE) + 0.1 * (max(df$gene, na.rm=TRUE) - min(df$gene, na.rm=TRUE)), y = 0.8, label = eqn)
    g <- g + geom_text(data = d.eqn, aes(x = x, y = y, label = label), parse=TRUE)
    glist2[[length(glist2)+1]] <- g
  }
  
  grid.arrange(grobs = c(glist1, glist2), ncol=2, as.table=FALSE)
  ##do.call("grid.arrange", glist)
  
  d <- dev.off()
}



m <- merge(fimm.mek.expr.all.ic50.corr, ohsu.mek.expr.all.ic50.corr, by = "gene", suffixes = c(".fimm", ".ohsu"))
pdf("fimm-ohsu-ic50-expr-gene-drug-venn.pdf")
g <- my.plot.venn(m, 0.05, "FIMM", "OHSU", "pval.fimm", "pval.ohsu", "gene", "gene", title = "significant genes (p < 0.05)\nlog(ic50) ~ expr + drug")
grid.draw(g)
d <- dev.off()

m$in.ras.sig <- 1
flg <- m$gene %in% ras.sig.genes$ENSG
m$in.ras.sig[flg] <- 0
g1 <- my.plot.venn(m, 0.05, "FIMM", "RAS", "pval.fimm", "in.ras.sig", "gene", "gene", title = "significant genes (p < 0.05)\nMEK log(ic50) ~ expr + drug")
g2 <- my.plot.venn(m, 0.05, "OHSU", "RAS", "pval.ohsu", "in.ras.sig", "gene", "gene", title = "significant genes (p < 0.05)\nMEK log(ic50) ~ expr + drug")
pdf("fimm-ohsu-ic50-expr-gene-drug-ras-venn.pdf")
grid.arrange(gTree(children=g1),gTree(children=g2), ncol=2)
d <- dev.off()

## 1. Plot a few of the genes that are in RAS signature
## 5. Do elastic net -- train on OHSU
## 6. Elastic net AUC venn
## 7. Elastic net IC50 venn
## 6. Elastic net AUC pval dist
## 7. Elastic net IC50 pval dist

## Intersect the results with Justin's RAS signature genes, which he showed were
## correlated with MEK inhibition.

d.tmp <- subset(ll4.dss.fits.no.force.seq, inhibitor %in% mek.annotations$ID_Drug.ohsu)

tmp.res <- ohsu.mek.expr.all.auc.corr
tmp.res <- tmp.res[tmp.res$pval != -1,]
head(tmp.res[order(as.numeric(tmp.res$pval)),c("pval","r2","gene")])

ras.res <- subset(ohsu.mek.expr.all.auc.corr, gene %in% ras.sig.genes$ENSG)
num.sig.ras.res <- length(which(!is.na(ras.res$pval) & (as.numeric(ras.res$pval < 0.05))))
frac.sig <- num.sig.ras.res / nrow(ras.res)
cat(paste0(num.sig.ras.res, " of ", nrow(ras.res), " (", round(100 * frac.sig, digits=2), "%) of genes in RAS signature are sig associated with MEK\n"))

## Let's do a simulation to see if this number is significant
num.iters <- 10000
res.pval <- sum(unlist(llply(1:num.iters, .parallel = TRUE, .fun = function(i) {
  indices <- sample.int(size = nrow(ras.res), n = nrow(ohsu.mek.expr.all.auc.corr), replace=FALSE)
  tmp.res <- ohsu.mek.expr.all.auc.corr[indices,]
  num.sig <- length(which(!is.na(tmp.res$pval) & (as.numeric(tmp.res$pval < 0.05))))
  if(num.sig >= num.sig.ras.res) {
    return(1)
  } else {
    return(0)
  }
})))
## This pvalue is not significant (around 0.17)
res.pval/num.iters

gene <- "ENSG00000104972"
gene <- "ENSG00000145416"
gene <- ras.sig.genes$ENSG[5]
gene <- "ENSG00000141505"
gene <- "ENSG00000139318"
gene <- "ENSG00000197956"
t.gene <- data.frame(gene = unname(t(ohsu.expr.data[gene,])), seqid = colnames(ohsu.expr.data))
m.tmp <- merge(d.tmp, t.gene, by.x = "SeqID", by.y = "seqid")
lmf <- lm(m.tmp$auc ~ m.tmp$gene + m.tmp$inhibitor)
plot(m.tmp$gene, m.tmp$auc)

## Model all MEK drugs simultaneously

## common.drugs.tested <- intersect(common.drugs.tested, c("FIMM136387", "FIMM133832", "FIMM133867", "FIMM133902"))

fimm.expr.corr <- correlate.expression.with.drug.response(drug.df = fimm.fits.gof.0.7.df, drug.id.col = "DRUG_ID", genes.to.test = common.genes, drugs.to.test = common.drugs.tested, expr.mat = fimm.expr.data, specimen.id.col = "SCREEN_ID", patient.id.col = "PATIENT_ID")



save.image(".RData")

cat("Fitting L4\n")
l4.fits <- fit.fct.to.ohsu(ohsu.raw.drug.data, fct = "L.4")
cat("Done fitting L4\n")

## Might need to exponentiate ic50?  yeah--log (or log10) is model
## l4.dss.fits <- compute.ohsu.dss(l4.fits)

save.image(".RData")
cat("Done saving .RData\n")
cat("Exiting successfully\n")
q()


## NB: ohsu data are viability, whereas compute.dss is written assuming inhibition data.
## So, need compute.dss max.asymptote_dss = 100 - min.asymptote_ohsu_fit
compute.all.dss <- function(data, invert.response = TRUE, ...) {
  ddply(data,
        colnames(data), 
        .parallel = TRUE,
        .fun = function(df) {
          auc <- df$auc[1]
          ic50 <- df$ic50[1]
          c.min <- df$min.conc[1]
          c.max <- df$max.conc[1]
          x1 <- c.min
          x2 <- c.max
          t <- 20
          if(df$min.asymptote[1] > df$max.asymptote[1]) {
            warning("min asymptote > max asymptote\n")
          }
          top.asymptote <- df$max.asymptote[1]
          if(invert.response) {
            top.asymptote <- 100 - df$min.asymptote[1]
            auc <- ( 100 * ( c.max - c.min ) ) - auc
          }
          dss <- compute.dss(auc, ic50, t, c.min, c.max, x1, x2, top.asymptote, ...)
          dss
        })
}

compute.ohsu.dss <- function(data, ...) {
  compute.all.dss(data, invert.response=TRUE, ...)
}

compute.fimm.dss <- function(data, ...) {
  compute.all.dss(data, invert.response=FALSE, ...)
}


## Make scatterplot of drugs sharing the same mechanism
mechanisms <- as.data.frame(table(ohsu.fimm.drugs$`Mechanism/Targets`))
colnames(mechanisms)[1] <- "mechanism"
mechanisms <- mechanisms$mechanism[mechanisms$Freq > 1]
for(mechanism in mechanisms) {
  mechanism.annotations <- ohsu.fimm.drugs[ohsu.fimm.drugs$`Mechanism/Targets` == mechanism,]
  ohsu.mechanism.data <- subset(ll4.fits, inhibitor %in% mechanism.annotations$ID_Drug.ohsu)
  for(metric in c("ic50", "auc")) {
    pdf(paste0("mechanism-", make.names(mechanism), "-", metric, "-pairs.pdf"))
    plot.pairs(ohsu.mechanism.data, metric, drug.col = "inhibitor", prefix = NULL, log.transform.IC = TRUE, lower.panel=panel.regression, upper.panel=panel.rank.cor, main = paste0(mechanism, ": IC50"))
    d <- dev.off()
  }
}

corDist <- function(x) {
  dst <- suppressWarnings(as.dist((1-cor(t(x), use="pairwise", method="spearman"))))
  dst.mat <- as.matrix(dst)
  dst.mat[!is.finite(dst.mat)] <- 1
  as.dist(dst.mat)
}

tmp <- ddply(ll4.fits, c("inhibitor", "patient_id"), 
             .parallel = TRUE,
             .fun = function(df) {
               median(df$auc)
             })


colnames(tmp) <- c("inhibitor", "patient_id", "auc")
ohsu.aucs <- spread(tmp, key = patient_id, value = auc)
rownames(ohsu.aucs) <- ohsu.aucs$inhibitor

## Restrict to those for which we have RNA-seq
ohsu.aucs <- ohsu.aucs[, colnames(ohsu.aucs) %in% rnaseq.sample.summary.tbl$PatientID]

mat <- t(ohsu.aucs)
flag <- unlist(apply(mat, 2, function(col) length(which(!is.na(col))) >= 10))
na.mat <- mat[,flag]

compound <- grepl(colnames(na.mat), pattern=" - ")
na.mat <- na.mat[,!compound]


## Cluster drugs based on AUC
distance <- corDist(t(na.mat))
hc <- hclust(distance)
attr(distance, "Labels") <- unlist(lapply(attr(distance, "Labels"), function(str) substr(str, 1, max(5:length(str)))))
cutoff <- 0.7
clusters <- cutree(hc, h=cutoff)
num.clusters <- length(unique(clusters))
pdf("drug-auc-spearman-clustering-monotherapies.pdf")
main <- paste0("Drug AUC Clustering\n(", num.clusters, " clusters at cutoff = ", cutoff, ")")
plot(hclust(distance), ylab = "1 - Spearman Correlation", main = main)
abline(h=cutoff)
d <- dev.off()

nipals <- pca(good.mat, method="nipals", center=TRUE, nPcs=5)


## Cluster patients based on AUC
distance <- corDist((na.mat))
hc <- hclust(distance)
attr(distance, "Labels") <- unlist(lapply(attr(distance, "Labels"), function(str) substr(str, 1, max(5:length(str)))))
cutoff <- 1
clusters <- cutree(hc, h=cutoff)
num.clusters <- length(unique(clusters))
pdf("patient-auc-spearman-clustering-monotherapies.pdf")
main <- paste0("Patient AUC Clustering\n(", num.clusters, " clusters at cutoff = ", cutoff, ")")
plot(hclust(distance), ylab = "1 - Spearman Correlation", main = main)
abline(h=cutoff)
d <- dev.off()

## Plot the patients, annotated with response
distance <- corDist((na.mat))
dend <- as.dendrogram(hclust(distance))

pdf("ohsu-patient-auc-annotated-response.pdf")
plot(dend, main = "Patient AUC clustering (Annotated Response)")
cols <- as.numeric(as.factor(ohsu.responses[labels(dend)]))
colored_bars(colors = cols, dend = dend, sort_by_labels_order = FALSE)
legend("topright", legend=unique(ohsu.responses[labels(dend)]), fill=as.numeric(as.factor(unique(ohsu.responses[labels(dend)]))))
d <- dev.off()

plot.missing <- function(mat) {
  frac.missing <- length(which(is.na(as.vector(mat))))/length(as.vector(mat))
  mat[!is.na(mat)] <- 1
  mat[is.na(mat)] <- 0
  per.missing = round(100 * frac.missing, digits=0)
  heatmap.2(mat, scale="none", trace="none", main = paste0(per.missing, "% missing data (red)"))
}

synId <- "syn7440451"
obj <- synGet(synId, downloadFile=TRUE)
ohsu.interpreted.data <- read.table(getFileLocation(obj), sep="\t", header=TRUE)

extract.monotherapy.matrix <- function(tbl, metric.col = "auc", drug.col = "inhibitor") {
  tmp <- ddply(tbl, c(drug.col, "patient_id"), 
               .parallel = TRUE,
               .fun = function(df) {
                 median(df[,metric.col])
              })

  colnames(tmp) <- c(drug.col, "patient_id", metric.col)
  mat <- spread_(tmp, key = "patient_id", metric.col)
  rownames(mat) <- mat[, drug.col]

  ## Restrict to those for which we have RNA-seq
  mat <- mat[, colnames(mat) %in% rnaseq.sample.summary.tbl$PatientID]

  mat <- t(mat)
  flag <- unlist(apply(mat, 2, function(col) length(which(!is.na(col))) >= 10))
  na.mat <- mat[,flag]

  compound <- grepl(colnames(na.mat), pattern=" - ")
  na.mat <- na.mat[,!compound]
  na.mat
}

interp.mat  <- extract.monotherapy.matrix(ohsu.interpreted.data[,c("drug", "patient_id", "Area_under_the_curve")], metric.col = "Area_under_the_curve", drug.col = "drug")


compound <- grepl(colnames(mat), pattern=" - ")
tmp <- mat[,!compound]

common.drugs <- intersect(colnames(tmp), colnames(interp.mat))
common.samples <- intersect(rownames(tmp), rownames(interp.mat))

pdf("monotherapy-missing-pre-fit-filtering.pdf")
plot.missing(tmp[as.character(common.samples), common.drugs])
d <- dev.off()

pdf("monotherapy-missing-interp.pdf")
plot.missing(interp.mat[as.character(common.samples), common.drugs])
d <- dev.off()


monotherapy <- !compound
ms <- unlist(apply(mat[,monotherapy], 2, function(col) length(which(is.na(col))) / length(col)))
plot(ms)
abline(h=0.25)
good.drugs <- names(ms)[ms < 0.25]

ms <- unlist(apply(t(mat[,monotherapy]), 2, function(col) length(which(is.na(col))) / length(col)))
plot(ms)
abline(h=0.4)
good.patients <- names(ms)[ms < 0.4]

tmp <- mat[,!compound]
good.mat <- tmp[as.character(good.patients), good.drugs]
plot.missing(good.mat)

plot.pcas <- function(mat, pca.methods = c("nipals", "bpca", "ppca", "svdImpute", "imputeSVD"), prefix = "prefix", k = 4, nPcs = 5) {
  
  ## Here, just using the various pca methods on the complete data
  smat <- scale(mat)
  ## nPcs <- min(5, ncol(mat))
  l <- llply(c(pca.methods), .parallel = TRUE,
             .fun = function(method) {
               print(method)
               if(method == "nlpca") {
                 return(pca(smat, method=method, center=FALSE, nPcs=nPcs, maxSteps=10000))
               }else if(method == "imputeSVD") {
                 imat <- llsImpute(mat, center=FALSE, k=k)
                 srmat.lo <- scale(completeObs(imat))
                 return(pca(srmat.lo, method="svd", center=FALSE, nPcs=nPcs))
               } else {
                 return(pca(smat, method=method, center=FALSE, nPcs=nPcs))
               }
             })
  names(l) <- c(pca.methods)

  ## Plot SVD vs all other PCA approaches
  gs <- list()
  frac <- 0
  for(method in pca.methods) {
    component <- 1
    df <- data.frame(x = l[[method]]@scores[,1], y = l[[method]]@scores[,2])
    print(length(df$x))
    g <- ggplot(df)
    g <- g + geom_point(aes(x = x, y = y))
    g <- g + xlab("PC1")
    g <- g + ylab("PC2")
    title <- method
    if(method == "imputeSVD") { title <- "impute then SVD" }
    g <- g + ggtitle(title)
    gs[[method]] <- g
  }
  pdf(paste0(prefix, "-missing.pdf"))
  do.call("grid.arrange", gs)
  d <- dev.off()
  
  ## Plot the eigenvalues
  gs <- list()
  frac <- 0
  for(method in pca.methods) {
    component <- 1
    eigs <- sDev(l[[method]])
    x <- 1:length(eigs)
    eigs <- 100 * eigs^2 
    if(nPcs == ncol(mat)) {
      eigs <- 100 * eigs^2 / sum(eigs^2)
    }
    
    df <- data.frame(x = x, y = eigs)
    g <- ggplot(df)
    g <- g + geom_point(aes(x = x, y = y))
    g <- g + xlab("Eigenvalue Number")
##    g <- g + ylab("Eigenvalue")
    if(nPcs == ncol(mat)) {
      g <- g + ylab("Variance Explained")
    } else {
      g <- g + ylab("Proportional to Variance Explained")
    }
    title <- method
    if(method == "imputeSVD") { title <- "impute then SVD" }
    g <- g + ggtitle(title)
    gs[[method]] <- g
  }
  pdf(paste0(prefix, "-eigenvalues-missing.pdf"))
  do.call("grid.arrange", gs)
  d <- dev.off()

}

plot.pcas(mat[,monotherapy], prefix="monotherapy-pcas", nPcs = ncol(mat[,monotherapy]), pca.methods = c("nipals", "bpca", "ppca", "svdImpute", "imputeSVD"))
plot.pcas(mat[,monotherapy], prefix="monotherapy-pcas", nPcs = ncol(mat[,monotherapy]), pca.methods = c("nipals", "svdImpute", "imputeSVD"))
plot.pcas(good.mat, prefix="good-monotherapy-pcas", nPcs = ncol(good.mat), pca.methods = c("nipals", "svdImpute", "imputeSVD"))


mek.annotations <- ohsu.fimm.drugs[ohsu.fimm.drugs$`Mechanism/Targets` == "MEK1/2 inhibitor",]
good.mek.mat <- good.mat[, colnames(good.mat) %in% mek.annotations$ID_Drug.ohsu]
flag <- unlist(apply(good.mek.mat, 1, function(row) any(is.na(row))))
plot.pcas(good.mek.mat[!flag,], nPcs = ncol(good.mek.mat[!flag,]), prefix="good-mek-pcas", k=2, pca.methods = c("svd", "nipals", "bpca", "ppca", "svdImpute"))

save.image(".RData")

nipals <- pca(good.mat, method="nipals", center=TRUE, nPcs=5)
plot(nipals@scores[,1], nipals@scores[,2])

## Make a drug-drug correlation plot
cor.mat <- cor((na.mat), use="pairwise.complete.obs")
cor.mat[lower.tri(cor.mat, diag=TRUE)] <- NA

## Now label the entries according to whether the two correlated drugs are in the same or different class
cor.melt.orig <- reshape2:::melt.matrix(cor.mat)
colnames(cor.melt.orig) <- c("inhibitor1", "inhibitor2", "cor")
cor.melt.orig <- cor.melt.orig[!is.na(cor.melt.orig$cor),]

drug.mechanisms <- ohsu.fimm.drugs[,c("ID_Drug.ohsu", "Mechanism/Targets")]

cor.melt <- merge(cor.melt.orig, ohsu.fimm.drugs[,c("ID_Drug.ohsu", "Mechanism/Targets")], by.x = "inhibitor1", by.y = "ID_Drug.ohsu", all=FALSE)
colnames(cor.melt)[4] <- "mechanism1"

cor.melt <- merge(cor.melt, ohsu.fimm.drugs[,c("ID_Drug.ohsu", "Mechanism/Targets")], by.x = "inhibitor2", by.y = "ID_Drug.ohsu", all=FALSE)
colnames(cor.melt)[5] <- "mechanism2"

flag <- cor.melt$mechanism1 == cor.melt$mechanism2
cor.melt$mechanism <- "Different"
cor.melt$mechanism[flag] <- "Same"

mean.diff <- abs(mean(cor.melt$cor[flag]) - mean(cor.melt$cor[!flag]))

## Calculate p value by simulation
num.samples <- 10000
p.val <- llply(1:num.samples, .parallel = TRUE,
              .fun = function(i) {
  drug.mechanisms <- ohsu.fimm.drugs[,c("ID_Drug.ohsu", "Mechanism/Targets")]
  drug.mechanisms.shuffle <- drug.mechanisms
  drug.mechanisms.shuffle$`Mechanism/Targets` <- drug.mechanisms.shuffle$`Mechanism/Targets`[sample.int(n = nrow(drug.mechanisms.shuffle), size = nrow(drug.mechanisms.shuffle), replace = FALSE)] 
    
  cor.melt.shuf <- merge(cor.melt.orig, drug.mechanisms.shuffle, by.x = "inhibitor1", by.y = "ID_Drug.ohsu", all=FALSE)
  colnames(cor.melt.shuf)[4] <- "mechanism1"
  
  cor.melt.shuf <- merge(cor.melt.shuf, drug.mechanisms.shuffle, by.x = "inhibitor2", by.y = "ID_Drug.ohsu", all=FALSE)
  colnames(cor.melt.shuf)[5] <- "mechanism2"
  
  f <- cor.melt.shuf$mechanism1 == cor.melt.shuf$mechanism2
  dff <- abs(mean(cor.melt.shuf$cor[f]) - mean(cor.melt.shuf$cor[!f]))
  p.val <- 0
  if(dff > mean.diff) {
    p.val <- 1
  }
  p.val
              })

p.val <- sum(unlist(p.val))/length(p.val)

g <- ggplot(cor.melt)
g <- g + geom_violin(aes(x = mechanism, y = cor))
g <- g + geom_boxplot(aes(x = mechanism, y = cor), width = 0.5)
g <- g + ylab("Drug-drug AUC correlation")
g <- g + xlab("Drug Target/Mechanisms")
g <- g + ggtitle(paste0("p-value: ", p.val))
pdf("monotherapy-auc-correlations-vs-mechanism.pdf")
print(g)
d <- dev.off()


## TODO
- do pca of patients based on auc
- do pca of drugs based on auc
- missing heatmap--but labeling according to missing, filtered, or present
- for each patients, drugs
for each clinical variable
cluster and color based on clinical
- repeat above but using pca

- email suleiman

mek.annotations <- ohsu.fimm.drugs[grepl(ohsu.fimm.drugs$`Mechanism/Targets`, pattern="MEK"),]
ohsu.mek.data <- subset(ohsu.interpreted.data, drug %in% mek.annotations$ID_Drug.ohsu)


## Filter fits
## min.gof:  exclude fits having a gof < min.gof
## remove.non.finite.ic50:  remove ic50 values that are NA or NaN.
## require.ic50.be.in.range:  require that ic50 be between min and max concentration tested
filter.drug.response.fits <- function(fits.df, min.gof = 0, remove.non.finite.ic50 = TRUE, require.ic50.be.in.range = TRUE) {
  if(remove.non.finite.ic50) {
    old.nrow <- nrow(fits.df)
    flag <- is.finite(fits.df$IC50)
    fits.df <- fits.df[flag,]
    new.nrow <- nrow(fits.df)
    cat(paste0("Filtering by finite IC50 reduces ", old.nrow, " fits to ", new.nrow, "\n"))
  }

  if(require.ic50.be.in.range) {
    old.nrow <- nrow(fits.df)
    flag <- is.finite(fits.df$IC50) & is.finite(fits.df$Max.Conc.tested) & is.finite(fits.df$Min.Conc.tested) & ((10^fits.df$IC50) >= fits.df$Min.Conc.tested) & ((10^fits.df$IC50) <= fits.df$Max.Conc.tested)
    fits.df <- fits.df[flag,]
    new.nrow <- nrow(fits.df)
    cat(paste0("Filtering by requiring IC50 be in range reduces ", old.nrow, " fits to ", new.nrow, "\n"))
  }
  
  if(min.gof > 0) {
    old.nrow <- nrow(fits.df)
    flag <- is.finite(fits.df$gof) & (fits.df$gof >= min.gof)
    fits.df <- fits.df[flag,]
    new.nrow <- nrow(fits.df)
    cat(paste0("Filtering by GOF < ", min.gof, " reduces ", old.nrow, " fits to ", new.nrow, "\n"))
  }
  
  fits.df
}

## common.drugs.tested <- intersect(unique(ohsu.raw.drug.data$ID_Drug), unique(fimm.raw.drug.data$DRUG_ID))

## cat("Common drugs\n")
## print(common.drugs.tested)

## sort( sapply(ls(),function(x){object.size(get(x))})) 

## Fit OHSU drug response curves
cat("Fitting OHSU drug response data\n")

print(colnames(ohsu.raw.drug.data))
## ohsu.raw.drug.data <- subset(ohsu.raw.drug.data, ID_Drug %in% common.drugs.tested)
ohsu.fits <- fit.ohsu.drug.response.data(ohsu.raw.drug.data)

cat("Done with OHSU drug response fits\n")

tmp <- lapply(ohsu.fits, function(x) x$tbl)
ohsu.fits.tbl <- do.call("rbind", tmp)

tmp <- lapply(ohsu.fits, function(x) x$id)
ohsu.fits.id <- do.call("rbind", tmp)

ohsu.fits.df <- cbind(ohsu.fits.tbl, ohsu.fits.id)
if(!(all(rownames(ohsu.fits.tbl) == rownames(ohsu.fits.id)))) {
  warning("UNEXPECTED ohsu.fits.tbl rownames != ohsu.fits.id rownames\n")
}
if(!(all(rownames(ohsu.fits.tbl) == rownames(ohsu.fits.df)))) {
  warning("UNEXPECTED ohsu.fits.tbl rownames != ohsu.fits.df rownames\n")
}
rm(ohsu.fits.tbl)
rm(ohsu.fits.id)
rm(tmp)
gc()

save.image(".RData.ohsu.fits")
cat("Done saving OHSU fits\n")
stop("stop")

## Fit FIMM drug response curves
## cat("Fitting FIMM drug response data\n")

## fimm.raw.drug.data <- subset(fimm.raw.drug.data, DRUG_ID %in% common.drugs.tested)
fimm.fits <- fit.fimm.drug.response.data(fimm.raw.drug.data)



## Collate the results
tmp <- lapply(fimm.fits, function(x) x$tbl)
fimm.fits.tbl <- do.call("rbind", tmp)

tmp <- lapply(fimm.fits, function(x) x$id)
fimm.fits.id <- do.call("rbind", tmp)

fimm.fits.df <- cbind(fimm.fits.tbl, fimm.fits.id)
if(!(all(rownames(fimm.fits.tbl) == rownames(fimm.fits.id)))) {
  warning("UNEXPECTED fimm.fits.tbl rownames != fimm.fits.id rownames\n")
}
if(!(all(rownames(fimm.fits.tbl) == rownames(fimm.fits.df)))) {
  warning("UNEXPECTED fimm.fits.tbl rownames != fimm.fits.df rownames\n")
}
rm(fimm.fits.tbl)
rm(fimm.fits.id)
gc()

## Add ensembl identifiers to the drug annotation table
## Translate gene symbols to ensg identifiers
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(Biobase))

symbols.to.ensg.mapping <- function(symbols) {
  # Get a mapping from gene symbol to ensembl id
  ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
  bm <- getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id'), 
              filters = 'hgnc_symbol', 
              values = symbols, 
              mart = ensembl)
  names(bm) <- c("SYMBOL", "ENSG")
  bm <- bm[!(bm$ENSG %in% c("")),]
  bm
}

all.targets <- unique(unlist(lapply(ohsu.fimm.drugs$Gene.Targets, function(x) unlist(strsplit(x, split=",[ ]*")))))
all.targets <- all.targets[!is.na(all.targets)]

sym.to.ensg <- symbols.to.ensg.mapping(all.targets)

ohsu.fimm.drugs$Ensg.Targets <- unlist(lapply(ohsu.fimm.drugs$Gene.Targets, 
                                              function(target.str) {
                                                targets <- unlist(strsplit(target.str, split=",[ ]*"))
                                                if(!(any(targets %in% sym.to.ensg$SYMBOL))) { 
                                                  return(NA)
                                                }
                                                targets <- targets[targets %in% sym.to.ensg$SYMBOL]
                                                targets <- unique(sym.to.ensg$ENSG[sym.to.ensg$SYMBOL %in% targets])
                                                paste(targets, collapse=", ")
                                              }))


## Merge the drug targets (ensembl identifiers and symbols) to the drug fits
old.nrow <- nrow(ohsu.fits.df)
ohsu.fits.df <- merge(ohsu.fits.df, unique(ohsu.fimm.drugs[,c("ID_Drug", "ID_Drug.ohsu", "Gene.Targets", "Ensg.Targets")]), by.x = "inhibitor", by.y = "ID_Drug.ohsu")
if(old.nrow != nrow(ohsu.fits.df)) {
  warning("UNEXPECTED size of ohsu.fits.df changed\n")
}

old.nrow <- nrow(fimm.fits.df)
fimm.fits.df <- merge(fimm.fits.df, unique(ohsu.fimm.drugs[,c("ID_Drug", "ID_Drug.ohsu", "Gene.Targets", "Ensg.Targets")]), by.x = "DRUG_ID", by.y = "ID_Drug")
if(old.nrow != nrow(fimm.fits.df)) {
  warning("UNEXPECTED size of fimm.fits.df changed\n")
}

## Add the translation from drug response labId to sequencing ID to the ohsu data
ohsu.fits.df <- merge(ohsu.fits.df, unique(rnaseq.sample.summary.tbl[,c("SeqID", "Original_LabID")]), by.x = "lab_id", by.y = "Original_LabID")

## When the sequence ids are used as column names, they will be converted from 12-00023 -> X12.00023.  Add this field so we can use to
## index columns so named.
ohsu.fits.df$ColSeqID <- gsub(pattern="-", x=ohsu.fits.df$SeqID, replacement=".")
ohsu.fits.df$ColSeqID <- unlist(lapply(ohsu.fits.df$ColSeqID, function(x) paste0("X", x)))

save.image(".RData.fits")
q()

ohsu.fits.gof.0.7.df <- filter.drug.response.fits(ohsu.fits.df, min.gof = 0.7, remove.non.finite.ic50 = TRUE, require.ic50.be.in.range = TRUE)
fimm.fits.gof.0.7.df <- filter.drug.response.fits(fimm.fits.df, min.gof = 0.7, remove.non.finite.ic50 = TRUE, require.ic50.be.in.range = TRUE)

common.drugs.tested <- intersect(unique(ohsu.fits.df$ID_Drug), unique(fimm.fits.df$DRUG_ID))
common.drugs.tested <- intersect(common.drugs.tested, ohsu.fimm.drugs$ID_Drug[!is.na(ohsu.fimm.drugs$Gene.Targets)])
common.drugs.tested <- intersect(common.drugs.tested, ohsu.fimm.drugs$ID_Drug[!is.na(ohsu.fimm.drugs$Ensg.Targets)])

## specimen ids are columns of expr.mat
## specimen.id.col is col of drug data frame holding the specimen id
## patient.id.col is col of drug data frame holding the patient id
## metric.col is IC50, AUC, DSS1, etc.
correlate.expression.with.drug.response <- function(drug.df, drug.id.col, drugs.to.test, genes.to.test = NULL, expr.mat, specimen.id.col, patient.id.col, metric.col = "IC50", log.transform = FALSE) {
  df <- drug.df[drug.df[,drug.id.col] %in% drugs.to.test,]
  if(is.null(genes.to.test)) {
    df <- df[!is.na(df$Ensg.Targets),]
  }
  if(is.null(genes.to.test)) {
    flag <- unlist(lapply(df$Ensg.Targets, function(x) any(unlist(strsplit(x, split=",[ ]*")) %in% rownames(expr.mat))))
    df <- df[flag,]
  }
  if(!is.null(genes.to.test)) {
    if(!(any(genes.to.test %in% rownames(expr.mat)))) {
      warning("UNEXPECTED no genes to test are in data\n")
    }
  }
  ret <- ddply(df, drug.id.col, .parallel = TRUE,
               .fun = function(df.drug) {
                 cat(paste0("Correlating gene expr with ", metric.col, " for drug ", unique(df.drug[,drug.id.col]), "\n"))
                 drug.targets <- NULL
                 if(!is.null(genes.to.test)) {
                   drug.targets <- genes.to.test
                 } else {
                   drug.targets <- df.drug$Ensg.Targets[1]
                   drug.targets <- unlist(strsplit(drug.targets, split=",[ ]*"))
                 }
                 drug.targets <- drug.targets[drug.targets %in% rownames(expr.mat)]
                 resp <- data.frame(metric = df.drug[,metric.col], PATIENT_ID = df.drug[,patient.id.col], id = df.drug[,specimen.id.col])
                 resp <- resp[!is.infinite(resp$metric) & !is.nan(resp$metric) & !is.na(resp$metric),]
                 if(log.transform) { resp$metric <- log(resp$metric) }
                 common.samples <- intersect(resp$id, colnames(expr.mat))
                 num.resp.samples <- length(resp$id)
                 num.unique.resp.samples <- length(unique(resp$id))
                 flag <- resp$id %in% common.samples
                 resp <- resp[flag,]
                 tmp <- as.data.frame(table(resp$PATIENT_ID))
                 colnames(tmp) <- c("PATIENT_ID", "weight")
                 tmp$weight <- 1/sqrt(tmp$weight)
                 resp <- merge(resp, tmp)
                 common.samples <- as.character(resp$id)
                 ldply(drug.targets, .parallel = FALSE,
                       .fun = function(gene) {
                         if(!(gene %in% rownames(expr.mat))) {
                           warning(paste0("UNEXPECTED: ", gene, " not in expr.mat\n"))
                         }
                         if(!(all(common.samples %in% colnames(expr.mat)))) {
                           warning("UNEXPECTED columns not in expr.mat\n")
                         }
                         gene.expr <- as.vector(unlist(expr.mat[gene, common.samples]))
                         lm.fit <- lm(resp$metric ~ gene.expr, weights = resp$weight)
                         sum <- summary(lm.fit)
                         f <- sum$fstatistic
                         r2 <- 0
                         p <- -1
                         if(is.numeric(f)) { 
                           r2 <- sum$r.squared
                           p <- pf(f[1],f[2],f[3],lower.tail=F)
                         }
                         num.unique.samples <- length(unique(common.samples))
                         num.samples <- length(common.samples)
                         res <- c(gene, p, r2, num.unique.samples, num.samples, num.unique.resp.samples, num.resp.samples)
                         names(res) <- c("gene", "pval", "r2", "num.uniq.samples", "num.samples", "num.uniq.resp.samples", "num.resp.samples")
                         res
                       })
                 })
  ret
}

correlate.expression.with.all.drug.response <- function(drug.df, drug.id.col, drugs.to.test, genes.to.test = NULL, expr.mat, specimen.id.col, patient.id.col, metric.col = "IC50", log.transform = FALSE) {
  df <- drug.df[drug.df[,drug.id.col] %in% drugs.to.test,]
  if(is.null(genes.to.test)) {
    df <- df[!is.na(df$Ensg.Targets),]
  }
  if(is.null(genes.to.test)) {
    flag <- unlist(lapply(df$Ensg.Targets, function(x) any(unlist(strsplit(x, split=",[ ]*")) %in% rownames(expr.mat))))
    df <- df[flag,]
  }
  if(!is.null(genes.to.test)) {
    if(!(any(genes.to.test %in% rownames(expr.mat)))) {
      warning("UNEXPECTED no genes to test are in data\n")
    }
  }
  cat(paste0("Correlating gene expr with ", metric.col, " for all drugs simultaneously: ", paste(unique(df[,drug.id.col]), collapse=","), "\n"))
  drug.targets <- NULL
  if(!is.null(genes.to.test)) {
    drug.targets <- genes.to.test
  } else {
    drug.targets <- df$Ensg.Targets
    drug.targets <- unlist(lapply(drug.targets, function(str) strsplit(str, split=",[ ]*")))
  }
  drug.targets <- drug.targets[drug.targets %in% rownames(expr.mat)]
  resp <- data.frame(metric = df[,metric.col], drug = df[,drug.id.col], PATIENT_ID = df[,patient.id.col], id = df[,specimen.id.col])
  resp <- resp[!is.infinite(resp$metric) & !is.nan(resp$metric) & !is.na(resp$metric),]
  if(log.transform) { resp$metric <- log(resp$metric) }
  resp$id2 <- unlist(apply(resp[,c("PATIENT_ID", "drug")], 1, function(row) paste(row, collapse="-"))) 
  common.samples <- intersect(resp$id, colnames(expr.mat))
  num.resp.samples <- length(resp$id2)
  num.unique.resp.samples <- length(unique(resp$id2))
  flag <- resp$id %in% common.samples
  resp <- resp[flag,]
##  tmp <- as.data.frame(table(resp$PATIENT_ID))
  tmp <- as.data.frame(table(resp$id2))
  colnames(tmp) <- c("id2", "weight")
  tmp$weight <- 1/sqrt(tmp$weight)
  resp <- merge(resp, tmp)
  common.samples <- as.character(resp$id)
  ret <- ldply(drug.targets, .parallel = TRUE,
               .fun = function(gene) {
                 if(!(gene %in% rownames(expr.mat))) {
                   warning(paste0("UNEXPECTED: ", gene, " not in expr.mat\n"))
                 }
                 if(!(all(common.samples %in% colnames(expr.mat)))) {
                   warning("UNEXPECTED columns not in expr.mat\n")
                 }
                 gene.expr <- as.vector(unlist(expr.mat[gene, common.samples]))
                 lm.fit <- lm(resp$metric ~ gene.expr + resp$drug, weights = resp$weight)
                 sum <- summary(lm.fit)
                 f <- sum$fstatistic
                 r2 <- 0
                 p <- -1
                 if(is.numeric(f)) { 
                   r2 <- sum$r.squared
                   ## Take the pvalue of the gene, not of the model
                   flg <- grepl(rownames(coefficients(sum)), pattern="gene")
                   if(length(which(flg)) == 1) {
                     p <- as.numeric(coefficients(sum)[flg,4])
                   }
                   ## p <- pf(f[1],f[2],f[3],lower.tail=F)
                 }
                 num.unique.samples <- length(unique(common.samples))
                 num.samples <- length(common.samples)
                 res <- c(gene, p, r2, num.unique.samples, num.samples, num.unique.resp.samples, num.resp.samples)
                 names(res) <- c("gene", "pval", "r2", "num.uniq.samples", "num.samples", "num.uniq.resp.samples", "num.resp.samples")
                 res
               })
        
  ret
}

## RAS signature genes from Justin's paper
ras.sig.genes <- data.frame(symbols = c("NPTX2", "LMAN2L", "IL33", "DUSP6", "MYL4", "PHLDA1", "EYA2", "PHIP", "KCNN4", "HSPB3", "SKP1", "LAMC2", "HOXB6", "SKAP1", "RNF39", "SLC17A1", "BANP", "TBXAS1", "KCNG1", "TRPV1", "ZBTB40", "EDEM2", "SNX3", "NINL", "HDLBP", "MYBPC1", "S100A6", "SYNJ1", "GABBR1", "APOBEC1", "FAM60A", "ZBTB48", "HOXB7", "ASGR1", "SYCP1", "EREG", "CDKAL1", "BNC1", "MFSD11", "C10orf68", "KRT81", "LRRC6", "IFNGR2", "KRT24", "STAMBP", "PRUNE2", "GABRQ", "TNNI1", "CAPN9"))
ras.sig.genes <- na.omit(symbols.to.ensg.mapping(ras.sig.genes$symbols))

common.genes <- intersect(rownames(ohsu.expr.data), rownames(fimm.expr.data))

ohsu.expr.corr <- correlate.expression.with.drug.response(drug.df = ohsu.fits.gof.0.7.df, drug.id.col = "ID_Drug", genes.to.test = common.genes, drugs.to.test = common.drugs.tested, expr.mat = ohsu.expr.data, specimen.id.col = "SeqID", patient.id.col = "patient_id")
cat("Done correlating expression with IC50 for OHSU\n")
save.image(".RData")

## common.drugs.tested <- intersect(common.drugs.tested, c("FIMM136387", "FIMM133832", "FIMM133867", "FIMM133902"))

fimm.expr.corr <- correlate.expression.with.drug.response(drug.df = fimm.fits.gof.0.7.df, drug.id.col = "DRUG_ID", genes.to.test = common.genes, drugs.to.test = common.drugs.tested, expr.mat = fimm.expr.data, specimen.id.col = "SCREEN_ID", patient.id.col = "PATIENT_ID")
cat("Done correlating expression with IC50 for FIMM\n")
cat("Done with expression\n")
save.image(".RData")
cat("Saved image\n")
stop("stop")

common.genes <- intersect(rownames(beat.mafs), rownames(fimm.genomic.bin))
ohsu.mut.corr <- correlate.genomic.with.ic50(drug.df = ohsu.fits.gof.0.7.df, drug.id.col = "ID_Drug", drugs.to.test = common.drugs.tested, genes.to.test = common.genes, mutation.mat = beat.mafs, specimen.id.col = "SeqID", patient.id.col = "patient_id")
fimm.mut.corr <- correlate.genomic.with.ic50(drug.df = fimm.fits.gof.0.7.df, drug.id.col = "DRUG_ID", drugs.to.test = common.drugs.tested, genes.to.test = common.genes, mutation.mat = fimm.genomic.bin, specimen.id.col = "SCREEN_ID", patient.id.col = "PATIENT_ID")

cat("Done with genomic\n")
save.image(".RData")
cat("Saved image\n")
stop("stop")

drug = "FIMM003774"
gene = "FLT3"
pdf(paste0(drug, "-", gene, "-ohsu-fimm-mut.pdf"))
g1 <- plot.genomic.vs.ic50(drug.df = ohsu.fits.gof.0.7.df, drug = drug, gene = gene, drug.id.col = "ID_Drug", mutation.mat = beat.mafs, specimen.id.col = "SeqID", patient.id.col = "patient_id")
g1 <- g1 + ggtitle("OHSU")
g2 <- plot.genomic.vs.ic50(drug.df = fimm.fits.gof.0.7.df, drug = drug, gene = gene, drug.id.col = "DRUG_ID", mutation.mat = fimm.genomic.bin, specimen.id.col = "SCREEN_ID", patient.id.col = "PATIENT_ID")
g2 <- g2 + ggtitle("FIMM")
grid.arrange(g1, g2)
d <- dev.off()

##for(drug in unique(fimm.mut.corr$DRUG_ID)) {
##  print(drug)
##  ohsu.mut.corr <- correlate.genomic.with.ic50(drug.df = ohsu.fits.gof.0.7.df, drug.id.col = "ID_Drug", drugs.to.test = drug, genes.to.test = common.genes, mutation.mat = beat.mafs, specimen.id.col = "SeqID", patient.id.col = "patient_id")
##}

ohsu.fimm.mut.corr <- merge(ohsu.mut.corr, fimm.mut.corr, by.x = c("ID_Drug", "gene"), by.y = c("DRUG_ID", "gene"), suffixes = c(".ohsu", ".fimm"))

library(gplots)
ohsu.fimm.mut.corr$drug.gene <- unlist(apply(ohsu.fimm.mut.corr[,c("ID_Drug", "gene")], 1, function(row) paste0(row[1], "-", row[2])))
vennList = list("FIMM"=ohsu.fimm.mut.corr$drug.gene[!is.na(ohsu.fimm.mut.corr$pval.fimm) & (ohsu.fimm.mut.corr$pval.fimm < 0.05)], "OHSU"=ohsu.fimm.mut.corr$drug.gene[!is.na(ohsu.fimm.mut.corr$pval.ohsu) & (ohsu.fimm.mut.corr$pval.ohsu < 0.05)])
pdf("ohsu-fimm-sig-mut-venn.pdf")
venn(vennList)
d <- dev.off()

ohsu.fimm.mut.corr$pval.fimm <- as.numeric(ohsu.fimm.mut.corr$pval.fimm)
ohsu.fimm.mut.corr$pval.ohsu <- as.numeric(ohsu.fimm.mut.corr$pval.ohsu)

## Plot pvalue distributions of FIMM and OHSU
pdf("ohsu-fimm-mut-pval-dist.pdf")
g1 <- ggplot(ohsu.fimm.mut.corr, aes(x=pval.fimm)) + geom_histogram() + xlab("FIMM p-value")
g2 <- ggplot(ohsu.fimm.mut.corr, aes(x=pval.ohsu)) + geom_histogram() + xlab("OHSU p-value")
grid.arrange(g1, g2)
d <- dev.off()


## specimen ids are columns of mutation.mat
## gene symbols are rows of mutation.mat
correlate.genomic.with.ic50 <- function(drug.df, drug.id.col, drugs.to.test, genes.to.test = NULL, mutation.mat, specimen.id.col, patient.id.col) {
  df <- drug.df[drug.df[,drug.id.col] %in% drugs.to.test,]
  if(is.null(genes.to.test)) {
    df <- df[!is.na(df$Gene.Targets),]
  }
  if(is.null(genes.to.test)) {
    flag <- unlist(lapply(df$Gene.Targets, function(x) any(unlist(strsplit(x, split=",[ ]*")) %in% rownames(mutation.mat))))
    df <- df[flag,]
  }
  if(!is.null(genes.to.test)) {
    if(!(any(genes.to.test %in% rownames(mutation.mat)))) {
      warning("UNEXPECTED no genes to test are in data\n")
    }
  }
  ret <- ddply(df, drug.id.col, .parallel = TRUE,
               .fun = function(df.drug) {
                 drug.targets <- df.drug$Gene.Targets[1]
                 drug.targets <- unlist(strsplit(drug.targets, split=",[ ]*"))
                 if(!is.null(genes.to.test)) {
                   drug.targets <- genes.to.test
                 }
                 drug.targets <- drug.targets[drug.targets %in% rownames(mutation.mat)]
                 resp <- data.frame(IC50 = df.drug$IC50, PATIENT_ID = df.drug[,patient.id.col], id = df.drug[,specimen.id.col])
                 resp <- resp[!is.infinite(resp$IC50) & !is.nan(resp$IC50) & !is.na(resp$IC50),]
                 common.samples <- intersect(resp$id, colnames(mutation.mat))
                 num.resp.samples <- length(resp$id)
                 num.unique.resp.samples <- length(unique(resp$id))
                 if(length(common.samples) == 0) {
                   return(ldply(drug.targets, .parallel = FALSE,
                                .fun = function(gene) {
                                  p <- NA
                                  num.unique.samples <- length(unique(common.samples))
                                  num.samples <- length(common.samples)
                                  res <- c(gene, p, num.unique.samples, num.samples, num.unique.resp.samples, num.resp.samples)
                                  names(res) <- c("gene", "pval", "num.uniq.samples", "num.samples", "num.uniq.resp.samples", "num.resp.samples")
                                  res
                                }))
                 }
                 flag <- resp$id %in% common.samples
                 resp <- resp[flag,]
                 tmp <- as.data.frame(table(resp$PATIENT_ID))
                 colnames(tmp) <- c("PATIENT_ID", "weight")
                 tmp$weight <- 1/sqrt(tmp$weight)
                 resp <- merge(resp, tmp)
                 common.samples <- as.character(resp$id)
                 ldply(drug.targets, .parallel = FALSE,
                       .fun = function(gene) {
                         if(!(gene %in% rownames(mutation.mat))) {
                           warning(paste0("UNEXPECTED: ", gene, " not in mutation.mat\n"))
                         }
                         if(!(all(common.samples %in% colnames(mutation.mat)))) {
                           warning("UNEXPECTED columns not in mutation.mat\n")
                         }
                         gene.mutation <- as.vector(unlist(mutation.mat[gene, common.samples]))
                         if(!(all(gene.mutation %in% c(0,1)))) {
                           warning(paste0("WARNING: non-binary gene mutations: ", paste(gene.mutation, collapse=", "), "\n"))
                         }
                         p <- NA
                         if(all((c(0,1) %in% gene.mutation))) {
                           test.df <- data.frame(x = factor(gene.mutation), y = resp$IC50)
                           wt <- wilcox.test(y ~ x, data = test.df)
                           p <- wt$p.value
                           print(wt)
                         }
                         num.unique.samples <- length(unique(common.samples))
                         num.samples <- length(common.samples)
                         res <- c(gene, p, num.unique.samples, num.samples, num.unique.resp.samples, num.resp.samples)
                         names(res) <- c("gene", "pval", "num.uniq.samples", "num.samples", "num.uniq.resp.samples", "num.resp.samples")
                         res
                       })
               })
  ret
}

plot.genomic.vs.ic50 <- function(drug.df, drug, gene, drug.id.col, mutation.mat, specimen.id.col, patient.id.col) {
  df.drug <- drug.df[drug.df[,drug.id.col] == drug,]
  resp <- data.frame(IC50 = df.drug$IC50, PATIENT_ID = df.drug[,patient.id.col], id = df.drug[,specimen.id.col])
  resp <- resp[!is.infinite(resp$IC50) & !is.nan(resp$IC50) & !is.na(resp$IC50),]
  common.samples <- intersect(resp$id, colnames(mutation.mat))
  flag <- resp$id %in% common.samples
  resp <- resp[flag,]
  common.samples <- as.character(resp$id)
  gene.mutation <- factor(as.vector(unlist(mutation.mat[gene, common.samples])))
  df <- data.frame(mutation = gene.mutation, resp = resp$IC50)  
  ggplot(df, aes(x = mutation, y = resp)) + 
    geom_beeswarm() +
    xlab(paste0(gene, "\nmutation status")) + ylab(paste0("Log10 ", drug, " IC50"))
}



## BEGIN load OHSU dna data


## beat.mafs <- read.mafs(maf.names = ohsu.beat.maf.names, maf.tbl = ohsu.beat.mafs, mutation.types.to.total = c("Missense_Mutation", "Nonsense_Mutation"))
beat.mafs <- read.mafs(maf.names = ohsu.beat.maf.names, maf.tbl = ohsu.beat.mafs, mutation.types.to.total = c("Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Splice_Site", "Translation_Start_Site"))
## Make the samples the columns
beat.mafs <- t(beat.mafs)

## END load OHSU dna data 


## HERE


cat("Fitting FIMM curves\n")
## This works with plyr 1.7.1 and 1.8.1, but not 1.8.4
fimm.fits <- dlply(fimm.raw.dr.subset[,c("DRUG_ID", "DRUG_NAME", "CONCENTRATION", "PERCENT_INHIBITION", "SCREEN_ID", "SAMPLE_DATE", "PATIENT_ID")],
                   c("SCREEN_ID", "DRUG_ID"),
                   .parallel = TRUE,
                   .fun = 
                     function(df) {
                       print(nrow(df))
                       id <- unique(df[, !(colnames(df) %in% c("CONCENTRATION", "PERCENT_INHIBITION")),drop=F])
                       if(nrow(df) < 2) {
                         return(list(id=id, tbl=NULL))
                       }
                       print(df)
                       dose <- df$CONCENTRATION
                       inhibition <- df$PERCENT_INHIBITION
                       patient.id <- df$PATIENT_ID[1]
                       inhibitor <- df$DRUG_NAME[1]
                       file <- paste0(patient.id, "-", inhibitor)
                       tbl <- NULL
                       fit <- NULL
                       res <- NULL
                       res <- suppressWarnings(drc.calc.ic50.robust(dose, inhibition, file, patient.id, inhibitor, return.fit = FALSE, return.residuals = FALSE, apply.constraints = FALSE, plot.fit = FALSE, robust = "mean", force.min.to.zero = FALSE))
                       if(!is.null(res) & !is.null(res$tbl)) {
                         tbl <- res$tbl
                       }
                       if(!is.null(res) & !is.null(res$res)) {
                         tbl <- res$tbl
                       }
                       list(id=id, tbl=tbl)
                    })


## WORKING

cat("Fitting OHSU curves\n")
## Fit curves to OHSU data (AML only, those drugs that overlap with FIMM drugs)
ohsu.inh.tbl.subset <- subset(ohsu.inh.tbl, (ohsu.inh.tbl$inhibitor %in% ohsu.fimm.drugs$ID_Drug.ohsu))
rm(ohsu.inh.tbl); gc()

ohsu.inh.tbl.subset <- ohsu.inh.tbl.subset[order(ohsu.inh.tbl.subset$inhibitor, ohsu.inh.tbl.subset$patient_id, ohsu.inh.tbl.subset$lab_id),]
ohsu.fits <- dlply(ohsu.inh.tbl.subset[,c("inhibitor", "well_concentration", "normalized_viability", "time_of_read", "lab_id", "replicant", "patient_id")],
                   c("patient_id", "inhibitor", "lab_id", "replicant", "time_of_read"),
                   .parallel = TRUE,
                   .fun = 
                     function(df) {
                       print(nrow(df))
                       id <- unique(df[, !(colnames(df) %in% c("well_concentration", "normalized_viability")),drop=F])
                       if(nrow(df) < 2) {
                         return(list(id=id, tbl=NULL))
                       }
                       print(df)
                       dose <- df$well_concentration
                       inhibition <- 100 - df$normalized_viability
                       patient.id <- df$patient_id[1]
                       inhibitor <- df$inhibitor[1]
                       file <- paste0(patient.id, "-", inhibitor)
                       tbl <- NULL
                       fit <- NULL
                       res <- NULL
                       res <- suppressWarnings(drc.calc.ic50.robust(dose, inhibition, file, patient.id, inhibitor, return.fit = FALSE, return.residuals = FALSE, apply.constraints = FALSE, plot.fit = FALSE, robust = "mean", force.min.to.zero = FALSE))
                       if(!is.null(res) & !is.null(res$tbl)) {
                         tbl <- res$tbl
                       }
                       if(!is.null(res) & !is.null(res$res)) {
                         tbl <- res$tbl
                       }
                       list(id=id, tbl=tbl)
                     })

cat("Done\n")
save.image(".RData")
stop("stop")

## Collect the results into a single table
tmp <- lapply(fimm.fits, function(x) x$tbl)
fimm.fits.tbl <- do.call("rbind", tmp)

tmp <- lapply(fimm.fits, function(x) x$id)
fimm.fits.id <- do.call("rbind", tmp)

fimm.fits.df <- cbind(fimm.fits.tbl, fimm.fits.id)
if(!(all(rownames(fimm.fits.tbl) == rownames(fimm.fits.id)))) {
  warning("UNEXPECTED fimm.fits.tbl rownames != fimm.fits.id rownames\n")
}
if(!(all(rownames(fimm.fits.tbl) == rownames(fimm.fits.df)))) {
  warning("UNEXPECTED fimm.fits.tbl rownames != fimm.fits.df rownames\n")
}
rm(fimm.fits.tbl)
rm(fimm.fits.id)
gc()

tmp <- lapply(ohsu.fits, function(x) x$tbl)
ohsu.fits.tbl <- do.call("rbind", tmp)

tmp <- lapply(ohsu.fits, function(x) x$id)
ohsu.fits.id <- do.call("rbind", tmp)

ohsu.fits.df <- cbind(ohsu.fits.tbl, ohsu.fits.id)
if(!(all(rownames(ohsu.fits.tbl) == rownames(ohsu.fits.id)))) {
  warning("UNEXPECTED ohsu.fits.tbl rownames != ohsu.fits.id rownames\n")
}
if(!(all(rownames(ohsu.fits.tbl) == rownames(ohsu.fits.df)))) {
  warning("UNEXPECTED ohsu.fits.tbl rownames != ohsu.fits.df rownames\n")
}
rm(ohsu.fits.tbl)
rm(ohsu.fits.id)
rm(tmp)
gc()

## Translate gene symbols to ensg identifiers
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(Biobase))

symbols.to.ensg.mapping <- function(symbols) {
  # Get a mapping from gene symbol to ensembl id
  ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
  bm <- getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id'), 
              filters = 'hgnc_symbol', 
              values = symbols, 
              mart = ensembl)
  names(bm) <- c("SYMBOL", "ENSG")
  bm <- bm[!(bm$ENSG %in% c("")),]
  bm
}

ensg.to.symbols.mapping <- function(ensg) {
  # Get a mapping from gene symbol to ensembl id
  ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
  bm <- getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id'), 
              filters = 'ensembl_gene_id', 
              values = ensg, 
              mart = ensembl)
  names(bm) <- c("SYMBOL", "ENSG")
  bm <- bm[!(bm$ENSG %in% c("")),]
  bm
}

all.targets <- unique(unlist(lapply(ohsu.fimm.drugs$Gene.Targets, function(x) unlist(strsplit(x, split=",[ ]*")))))
all.targets <- all.targets[!is.na(all.targets)]

sym.to.ensg <- symbols.to.ensg.mapping(all.targets)

ohsu.fimm.drugs$Ensg.Targets <- unlist(lapply(ohsu.fimm.drugs$Gene.Targets, 
                                           function(target.str) {
                                             targets <- unlist(strsplit(target.str, split=",[ ]*"))
                                             if(!(any(targets %in% sym.to.ensg$SYMBOL))) { 
                                               return(NA)
                                             }
                                             targets <- targets[targets %in% sym.to.ensg$SYMBOL]
                                             targets <- unique(sym.to.ensg$ENSG[sym.to.ensg$SYMBOL %in% targets])
                                             paste(targets, collapse=", ")
                                           }))

## Look at correlation of drug IC50 with drug target gene expression for both
## FIMM and OHSU.  Compare p-values between the two.

old.nrow <- nrow(ohsu.fits.df)
ohsu.fits.df <- merge(ohsu.fits.df, unique(ohsu.fimm.drugs[,c("ID_Drug", "ID_Drug.ohsu", "Gene.Targets", "Ensg.Targets")]), by.x = "inhibitor", by.y = "ID_Drug.ohsu")
if(old.nrow != nrow(ohsu.fits.df)) {
  warning("UNEXPECTED size of ohsu.fits.df changed\n")
}

old.nrow <- nrow(fimm.fits.df)
fimm.fits.df <- merge(fimm.fits.df, unique(ohsu.fimm.drugs[,c("ID_Drug", "ID_Drug.ohsu", "Gene.Targets", "Ensg.Targets")]), by.x = "DRUG_ID", by.y = "ID_Drug")
if(old.nrow != nrow(fimm.fits.df)) {
  warning("UNEXPECTED size of fimm.fits.df changed\n")
}

g1 <- ggplot(data = ohsu.fits.df, aes(x = gof))
g1 <- g1 + geom_density() + xlab("OHSU GOF")
g2 <- ggplot(data = fimm.fits.df, aes(x = gof))
g2 <- g2 + geom_density() + xlab("FIMM GOF")
pdf("gof-density.pdf")
grid.arrange(g1, g2)
d <- dev.off()

all.ensg.targets <- unique(unlist(lapply(ohsu.fimm.drugs$Ensg.Targets, function(x) unlist(strsplit(x, split=",[ ]*")))))
all.ensg.targets <- all.ensg.targets[!is.na(all.ensg.targets)]

common.genes <- intersect(rownames(fimm.expr), rownames(ohsu.expr))
common.genes <- intersect(common.genes, all.ensg.targets)

fimm.expr.common <- fimm.expr[common.genes,]
ohsu.expr.common <- ohsu.expr[common.genes,]

## Correlate gene target expression with target IC50 for FIMM
df <- fimm.fits.df[fimm.fits.df$DRUG_ID %in% common.drugs.tested,]
flag <- unlist(lapply(df$Ensg.Targets, function(x) any(unlist(strsplit(x, split=",[ ]*")) %in% common.genes)))
df <- df[flag,]
fimm.expr.df <- merge(df, ohsu.fimm.drugs[,c("ID_Drug", "Gene.Targets")], by.x = "DRUG_ID", by.y = "ID_Drug")
fimm.expr.corr <- ddply(fimm.expr.df, c("DRUG_ID"), .parallel = TRUE,
                        .fun = function(df.drug) {
                          drug.targets <- df.drug$Ensg.Targets[1]
                          drug.targets <- unlist(strsplit(drug.targets, split=",[ ]*"))
                          drug.targets <- drug.targets[drug.targets %in% rownames(fimm.expr.common)]
                          resp <- data.frame(IC50 = df.drug$IC50, PATIENT_ID = df.drug$PATIENT_ID, id = df.drug$SCREEN_ID)
                          rownames(resp) <- resp$id
                          resp <- resp[!is.infinite(resp$IC50) & !is.nan(resp$IC50) & !is.na(resp$IC50),]
                          common.samples <- intersect(rownames(resp), colnames(fimm.expr.common))
                          resp <- resp[common.samples,]
                          tmp <- as.data.frame(table(resp$PATIENT_ID))
                          colnames(tmp) <- c("PATIENT_ID", "weight")
                          tmp$weight <- 1/sqrt(tmp$weight)
                          resp <- merge(resp, tmp)
                          rownames(resp) <- resp$id
                          if(!(all(common.samples %in% rownames(resp)))) {
                            warning("UNEXPECTED row names")
                          }
                          resp <- resp[common.samples,]
                          ldply(drug.targets, .parallel = FALSE,
                                .fun = function(gene) {
                                  if(!(gene %in% rownames(fimm.expr.common))) {
                                    warning(paste0("UNEXPECTED: ", gene, " not in fimm.expr\n"))
                                  }
                                  if(!(all(common.samples %in% colnames(fimm.expr.common)))) {
                                    warning("UNEXPECTED columns not in fimm expr\n")
                                  }
                                  gene.expr <- fimm.expr.common[gene, common.samples]
                                  lm.fit <- lm(resp$IC50 ~ gene.expr, weights = resp$weight)
                                  sum <- summary(lm.fit)
                                  f <- sum$fstatistic
                                  r2 <- sum$r.squared
                                  p <- pf(f[1],f[2],f[3],lower.tail=F)
                                  res <- c(gene, p, r2)
                                  num.unique.samples <- length(unique(common.samples))
                                  num.samples <- length(common.samples)
                                  res <- c(gene, p, r2, num.unique.samples, num.samples)
                                  names(res) <- c("gene", "pval", "r2", "num.uniq.samples", "num.samples")
                                  res
                                })
                        })

plot.fimm.expr.vs.ic50 <- function(df, drug, ensg.gene, expr.mat, xlab, ylab) {
  df.drug <- subset(df, DRUG_ID == drug)
  resp <- data.frame(IC50 = df.drug$IC50, PATIENT_ID = df.drug$PATIENT_ID, id = df.drug$SCREEN_ID)
  rownames(resp) <- resp$id
  resp <- resp[!is.infinite(resp$IC50) & !is.nan(resp$IC50) & !is.na(resp$IC50),]
  common.samples <- intersect(rownames(resp), colnames(expr.mat))
  resp <- resp[common.samples,]
  tmp <- as.data.frame(table(resp$PATIENT_ID))
  colnames(tmp) <- c("PATIENT_ID", "weight")
  tmp$weight <- 1/sqrt(tmp$weight)
  resp <- merge(resp, tmp)
  rownames(resp) <- resp$id
  if(!(all(common.samples %in% rownames(resp)))) {
    warning("UNEXPECTED row names")
  }
  if(!(ensg.gene %in% rownames(expr.mat))) {
    warning(paste0("UNEXPECTED: ", ensg.gene, " not in expr.mat\n"))
  }
  if(!(all(common.samples %in% colnames(expr.mat)))) {
    warning("UNEXPECTED columns not in expr.mat\n")
  }
  gene.expr <- fimm.expr.common[ensg.gene, common.samples]
  df <- data.frame(expr = gene.expr, resp = resp$IC50, weights=resp$weight)
  ggplot(df, aes(x = expr, y = resp, weight=weights)) + 
    geom_point() + 
    xlab(xlab) + ylab(ylab) +
    stat_smooth(method = "lm", col = "red")  
}


## Correlate gene target expression with target IC50 for OHSU
df <- ohsu.fits.df[ohsu.fits.df$ID_Drug %in% common.drugs.tested,]
flag <- unlist(lapply(df$Ensg.Targets, function(x) any(unlist(strsplit(x, split=",[ ]*")) %in% common.genes)))
df <- df[flag,]
df <- merge(df, ohsu.fimm.drugs[,c("ID_Drug", "Gene.Targets")], by.x = "ID_Drug", by.y = "ID_Drug")
## What is this lab_id
## Need to translate lab_id to seq_id
old.nrow <- nrow(df)
new.df <- merge(df, unique(aml.rnaseq.sample.summary.tbl[,c("SeqID", "Original_LabID")]), by.x = "lab_id", by.y = "Original_LabID")
ohsu.expr.df <- new.df
ohsu.expr.corr <- ddply(ohsu.expr.df, c("ID_Drug"), .parallel = TRUE,
                        .fun = function(df.drug) {
                          drug.targets <- df.drug$Ensg.Targets[1]
                          drug.targets <- unlist(strsplit(drug.targets, split=",[ ]*"))
                          drug.targets <- drug.targets[drug.targets %in% rownames(ohsu.expr.common)]
                          ids <- df.drug$SeqID
                          ids <- gsub(pattern="-", x=ids, replacement=".")
                          ids <- unlist(lapply(ids, function(x) paste0("X", x)))
                          resp <- data.frame(IC50 = df.drug$IC50, PATIENT_ID = df.drug$patient_id, id = ids)
                          resp <- resp[!is.infinite(resp$IC50) & !is.nan(resp$IC50) & !is.na(resp$IC50),]
                          common.samples <- intersect(resp$id, colnames(ohsu.expr.common))
                          flag <- resp$id %in% common.samples
                          resp <- resp[flag,]
                          tmp <- as.data.frame(table(resp$PATIENT_ID))
                          colnames(tmp) <- c("PATIENT_ID", "weight")
                          tmp$weight <- 1/sqrt(tmp$weight)
                          resp <- merge(resp, tmp)
                          common.samples <- as.character(resp$id)
                          ldply(drug.targets, .parallel = FALSE,
                                .fun = function(gene) {
                                  if(!(gene %in% rownames(ohsu.expr.common))) {
                                    warning(paste0("UNEXPECTED: ", gene, " not in ohsu.expr\n"))
                                  }
                                  if(!(all(common.samples %in% colnames(ohsu.expr.common)))) {
                                    warning("UNEXPECTED columns not in ohsu expr\n")
                                  }
                                  gene.expr <- as.vector(unlist(ohsu.expr.common[gene, common.samples]))
                                  lm.fit <- lm(resp$IC50 ~ gene.expr, weights = resp$weight)
                                  sum <- summary(lm.fit)
                                  f <- sum$fstatistic
                                  r2 <- sum$r.squared
                                  p <- pf(f[1],f[2],f[3],lower.tail=F)
                                  num.unique.samples <- length(unique(common.samples))
                                  num.samples <- length(common.samples)
                                  res <- c(gene, p, r2, num.unique.samples, num.samples)
                                  names(res) <- c("gene", "pval", "r2", "num.uniq.samples", "num.samples")
                                  res
                                })
                        })

plot.ohsu.expr.vs.ic50 <- function(df, drug, ensg.gene, expr.mat, xlab, ylab) {
  df.drug <- subset(df, ID_Drug == drug)
  ids <- df.drug$SeqID
  ids <- gsub(pattern="-", x=ids, replacement=".")
  ids <- unlist(lapply(ids, function(x) paste0("X", x)))
  resp <- data.frame(IC50 = df.drug$IC50, PATIENT_ID = df.drug$patient_id, id = ids)
  resp <- resp[!is.infinite(resp$IC50) & !is.nan(resp$IC50) & !is.na(resp$IC50),]
  common.samples <- intersect(resp$id, colnames(expr.mat))
  flag <- resp$id %in% common.samples
  resp <- resp[flag,]
  tmp <- as.data.frame(table(resp$PATIENT_ID))
  colnames(tmp) <- c("PATIENT_ID", "weight")
  tmp$weight <- 1/sqrt(tmp$weight)
  resp <- merge(resp, tmp)
  common.samples <- as.character(resp$id)
  if(!(ensg.gene %in% rownames(expr.mat))) {
    warning(paste0("UNEXPECTED: ", ensg.gene, " not in expr.mat\n"))
  }
  if(!(all(common.samples %in% colnames(expr.mat)))) {
    warning("UNEXPECTED columns not in expr.mat\n")
  }
  gene.expr <- as.vector(unlist(expr.mat[ensg.gene, common.samples]))
  df <- data.frame(expr = gene.expr, resp = resp$IC50, weights=resp$weight)
  ggplot(df, aes(x = expr, y = resp, weight=weights)) + 
    geom_point() + 
    xlab(xlab) + ylab(ylab) +
    stat_smooth(method = "lm", col = "red")  
}

ohsu.fimm.expr.corr <- merge(ohsu.expr.corr, fimm.expr.corr, by.x = c("ID_Drug", "gene"), by.y = c("DRUG_ID", "gene"), suffixes = c(".ohsu", ".fimm"))

subset(ohsu.fimm.expr.corr, (pval.ohsu < 0.05) & (pval.fimm < 0.05))
## ID_Drug            gene          pval.ohsu             r2.ohsu          pval.fimm           r2.fimm
## 100 FIMM003783 ENSG00000094631 0.0419410434254526  0.0494985601145971 0.0330721048250537 0.228393094777627
## 127 FIMM003794 ENSG00000122025 0.0140628306932017 0.00833458654791967 0.0449806673025468  0.24180451968676
ohsu.fimm.expr.corr$pval.fimm <- as.numeric(ohsu.fimm.expr.corr$pval.fimm)
ohsu.fimm.expr.corr$pval.ohsu <- as.numeric(ohsu.fimm.expr.corr$pval.ohsu)

## Plot pvalue distributions of FIMM and OHSU
pdf("ohsu-fimm-pval-dist.pdf")
g1 <- ggplot(ohsu.fimm.expr.corr, aes(x=pval.fimm)) + geom_histogram() + xlab("FIMM p-value")
g2 <- ggplot(ohsu.fimm.expr.corr, aes(x=pval.ohsu)) + geom_histogram() + xlab("OHSU p-value")
grid.arrange(g1, g2)
d <- dev.off()

ohsu.fimm.expr.corr$drug.gene <- unlist(apply(ohsu.fimm.expr.corr[,c("ID_Drug", "gene")], 1, function(row) paste0(row[1], "-", row[2])))
## Make Venn diagram of FIMM vs OHSU significant
library(gplots)
vennList = list("FIMM"=ohsu.fimm.expr.corr$drug.gene[ohsu.fimm.expr.corr$pval.fimm < 0.05], "OHSU"=ohsu.fimm.expr.corr$drug.gene[ohsu.fimm.expr.corr$pval.ohsu < 0.05])
pdf("ohsu-fimm-sig-venn.pdf")
venn(vennList)
d <- dev.off()

old.nrow <- nrow(ohsu.fimm.expr.corr)
ohsu.fimm.expr.corr <- merge(ohsu.fimm.expr.corr, unique(fimm.ohsu.drug.annotations[,c("ID_Drug", "DRUG_NAME")]))
if(old.nrow != nrow(ohsu.fimm.expr.corr)) {
  warning("UNEXPECTED drug correlation results grew\n")
}

old.nrow <- nrow(ohsu.fimm.expr.corr)
ohsu.fimm.expr.corr <- merge(ohsu.fimm.expr.corr, unique(sym.to.ensg), by.x = "gene", by.y = "ENSG")
if(old.nrow != nrow(ohsu.fimm.expr.corr)) {
  warning("UNEXPECTED drug correlation results grew\n")
}

## Plot few cases that have significant correlation in both
both.sig <- subset(ohsu.fimm.expr.corr, (pval.ohsu < 0.05) & (pval.fimm < 0.05))
for(i in 1:nrow(both.sig)) {
  ensg.gene <- both.sig$gene[i]
  drug.id <- both.sig$ID_Drug[i]
  pval.ohsu <- both.sig$pval.ohsu[i]
  pval.fimm <- both.sig$pval.fimm[i]
  
  gene.sym <- sym.to.ensg$SYMBOL[sym.to.ensg$ENSG == ensg.gene]
  drug.name <- fimm.ohsu.drug.annotations$DRUG_NAME[fimm.ohsu.drug.annotations$ID_Drug == drug.id]
  g1 <- plot.fimm.expr.vs.ic50(fimm.expr.df, drug.id, ensg.gene, fimm.expr.common, paste0("Log2 CPM (", gene.sym, ")"), paste0("Log10 IC50 (", drug.name, ")"))
  g1 <- g1 + ggtitle(paste0("FIMM (pval = ", pval.fimm, ")"))
  g2 <- plot.ohsu.expr.vs.ic50(ohsu.expr.df, drug.id, ensg.gene, ohsu.expr.common, paste0("Log2 CPM (", gene.sym, ")"), paste0("Log10 IC50 (", drug.name, ")"))
  g2 <- g2 + ggtitle(paste0("OHSU (pval = ", pval.ohsu, ")"))
  pdf(paste0("both-sig-", ensg.gene, "-", drug.id, ".pdf"))
  grid.arrange(g1, g2)
  d <- dev.off()
}

## Plot the top hits in each
sig <- subset(ohsu.fimm.expr.corr, pval.ohsu < 0.01)
for(i in 1:nrow(sig)) {
  ensg.gene <- sig$gene[i]
  drug.id <- sig$ID_Drug[i]
  pval.ohsu <- sig$pval.ohsu[i]
  pval.fimm <- sig$pval.fimm[i]
  
  gene.sym <- sym.to.ensg$SYMBOL[sym.to.ensg$ENSG == ensg.gene]
  drug.name <- fimm.ohsu.drug.annotations$DRUG_NAME[fimm.ohsu.drug.annotations$ID_Drug == drug.id]
  g1 <- plot.fimm.expr.vs.ic50(fimm.expr.df, drug.id, ensg.gene, fimm.expr.common, paste0("Log2 CPM (", gene.sym, ")"), paste0("Log10 IC50 (", drug.name, ")"))
  g1 <- g1 + ggtitle(paste0("FIMM (pval = ", pval.fimm, ")"))
  g2 <- plot.ohsu.expr.vs.ic50(ohsu.expr.df, drug.id, ensg.gene, ohsu.expr.common, paste0("Log2 CPM (", gene.sym, ")"), paste0("Log10 IC50 (", drug.name, ")"))
  g2 <- g2 + ggtitle(paste0("OHSU (pval = ", pval.ohsu, ")"))
  pdf(paste0("ohsu-sig-", ensg.gene, "-", drug.id, ".pdf"))
  grid.arrange(g1, g2)
  d <- dev.off()
}

sig <- subset(ohsu.fimm.expr.corr, pval.fimm < 0.01)
for(i in 1:nrow(sig)) {
  ensg.gene <- sig$gene[i]
  drug.id <- sig$ID_Drug[i]
  pval.ohsu <- sig$pval.ohsu[i]
  pval.fimm <- sig$pval.fimm[i]
  
  gene.sym <- sym.to.ensg$SYMBOL[sym.to.ensg$ENSG == ensg.gene]
  drug.name <- fimm.ohsu.drug.annotations$DRUG_NAME[fimm.ohsu.drug.annotations$ID_Drug == drug.id]
  g1 <- plot.fimm.expr.vs.ic50(fimm.expr.df, drug.id, ensg.gene, fimm.expr.common, paste0("Log2 CPM (", gene.sym, ")"), paste0("Log10 IC50 (", drug.name, ")"))
  g1 <- g1 + ggtitle(paste0("FIMM (pval = ", pval.fimm, ")"))
  g2 <- plot.ohsu.expr.vs.ic50(ohsu.expr.df, drug.id, ensg.gene, ohsu.expr.common, paste0("Log2 CPM (", gene.sym, ")"), paste0("Log10 IC50 (", drug.name, ")"))
  g2 <- g2 + ggtitle(paste0("OHSU (pval = ", pval.ohsu, ")"))
  pdf(paste0("fimm-sig-", ensg.gene, "-", drug.id, ".pdf"))
  grid.arrange(g1, g2)
  d <- dev.off()
}


## Plot the best cases in both

## - for intersected drugs, use symbols/ensg for genes


## - look at correlation of expr in gene targets with IC50
## - look at correlation of mutant in gene targets with IC50 (only those that have a lot of mutations)

## - 10-fold cross validation of ohsu using expr
## - 10-fold cross validation of ohsu using genomic

## Get some patients for Cristina
uniq.ohsu <- unique(ohsu.raw.drug.data[,c("inhibitor", "patient_id", "diagnosis", "specific_diagnosis", "lab_id", "specimen_type", "run_type", "replicant", "time_of_read")])

flag <- (uniq.ohsu$lab_id %in% colnames(ohsu.expr.data)) & (uniq.ohsu$lab_id %in% ohsu.beat.maf.names)

write.table(uniq.ohsu[flag, ], sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE, file="ohsu-drug-screens-with-rna-and-dna-seq.tsv")