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

library(ggbeeswarm)
library(plyr)
library(stringr)
library(synapseClient)
synapseLogin()

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
beatMafs     <- origBeatMafs[grep("mutect",origBeatMafs$file.name),]
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
filtered.res <- filter.drug.response.data.sets(ohsu.raw.drug.data = ohsu.inh.tbl, fimm.raw.drug.data = fimm.raw.dr, remove.ohsu.zero.responses = TRUE, max.ohsu.response = 150, min.ohsu.responses.in.screen = 5, max.ohsu.responses.in.screen = 7, require.intersection = TRUE)
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

## Fit FIMM drug response curves
cat("Fitting FIMM drug response data\n")

common.drugs.tested <- intersect(unique(ohsu.raw.drug.data$ID_Drug), unique(fimm.raw.drug.data$DRUG_ID))

cat("Common drugs\n")
print(common.drugs.tested)

fimm.raw.drug.data <- subset(fimm.raw.drug.data, DRUG_ID %in% common.drugs.tested)
fimm.fits <- fit.fimm.drug.response.data(fimm.raw.drug.data)

## sort( sapply(ls(),function(x){object.size(get(x))})) 

## Fit OHSU drug response curves
cat("Fitting OHSU drug response data\n")

print(colnames(ohsu.raw.drug.data))
ohsu.raw.drug.data <- subset(ohsu.raw.drug.data, ID_Drug %in% common.drugs.tested)
ohsu.fits <- fit.ohsu.drug.response.data(ohsu.raw.drug.data)

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
correlate.expression.with.ic50 <- function(drug.df, drug.id.col, drugs.to.test, genes.to.test = NULL, expr.mat, specimen.id.col, patient.id.col) {
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
                 cat(paste0("Correlating gene expr with IC50 for drug ", unique(df.drug[,drug.id.col]), "\n"))
                 drug.targets <- df.drug$Ensg.Targets[1]
                 drug.targets <- unlist(strsplit(drug.targets, split=",[ ]*"))
                 if(!is.null(genes.to.test)) {
                   drug.targets <- genes.to.test
                 }
                 drug.targets <- drug.targets[drug.targets %in% rownames(expr.mat)]
                 resp <- data.frame(IC50 = df.drug$IC50, PATIENT_ID = df.drug[,patient.id.col], id = df.drug[,specimen.id.col])
                 resp <- resp[!is.infinite(resp$IC50) & !is.nan(resp$IC50) & !is.na(resp$IC50),]
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
                         lm.fit <- lm(resp$IC50 ~ gene.expr, weights = resp$weight)
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

common.genes <- intersect(rownames(ohsu.expr.data), rownames(fimm.expr.data))

ohsu.expr.corr <- correlate.expression.with.ic50(drug.df = ohsu.fits.gof.0.7.df, drug.id.col = "ID_Drug", genes.to.test = common.genes, drugs.to.test = common.drugs.tested, expr.mat = ohsu.expr.data, specimen.id.col = "SeqID", patient.id.col = "patient_id")
cat("Done correlating expression with IC50 for OHSU\n")
save.image(".RData")

## common.drugs.tested <- intersect(common.drugs.tested, c("FIMM136387", "FIMM133832", "FIMM133867", "FIMM133902"))

fimm.expr.corr <- correlate.expression.with.ic50(drug.df = fimm.fits.gof.0.7.df, drug.id.col = "DRUG_ID", genes.to.test = common.genes, drugs.to.test = common.drugs.tested, expr.mat = fimm.expr.data, specimen.id.col = "SCREEN_ID", patient.id.col = "PATIENT_ID")
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