library(synapseClient)
library(data.table)
library(ggplot2)
library(gridExtra)
library(tidyr)
library( ReporteRs )

## Create a ppt presentation holding these results
mydoc <- pptx()

mydoc <- addSlide( mydoc, slide.layout = 'Title Slide' )
mydoc <- addTitle( mydoc, 'L.4 vs LL.4 QC' )
mydoc <- addSubtitle( mydoc , 'June 20, 2017')

synapseLogin()

plot.smooth.scatter <- function(data, x.col, y.col, x.label, y.label) {
  g <- ggplot(data, aes_string(x = x.col, y = y.col))
  g <- g + stat_density2d(aes(fill = ..density..^0.25), geom = "tile", contour = FALSE, n = 200)
  g <- g + scale_fill_continuous(low = "white", high = "dodgerblue4", guide=FALSE)
##  g <- g + geom_point(alpha = 0.1, shape = 20)
  g <- g + xlab(x.label)
  g <- g + ylab(y.label)
  g
}

plot.scatter <- function(data, x.col, y.col, x.label, y.label) {
  g <- ggplot(data, aes_string(x = x.col, y = y.col))
  g <- g + geom_point()
  g <- g + xlab(x.label)
  g <- g + ylab(y.label)
  g
}

remove.outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

## Load in the FIMM and OHSU fits.

## path <- "fimm.dss.t0.tsv"
synId <- "syn10083488"
obj <- synGet(synId, downloadFile = TRUE)
file <- getFileLocation(obj)
fimm.dss.t0 <- as.data.frame(fread(file))

## The e parameter in both L.4 and LL.4 is IC50 (not log IC50).
## Remove outliers and then plot their comparison
fimm.dss.t0$ic50.ll4 <- remove.outliers(fimm.dss.t0$e.ll4)
fimm.dss.t0$ic50.l4 <- remove.outliers(fimm.dss.t0$e.l4)

## path <- "fimm.dss.t10.tsv"
synId <- "syn10083490"
obj <- synGet(synId, downloadFile = TRUE)
file <- getFileLocation(obj)
fimm.dss.t10 <- as.data.frame(fread(file))
fimm.dss.t10$ic50.ll4 <- remove.outliers(fimm.dss.t10$e.ll4)
fimm.dss.t10$ic50.l4 <- remove.outliers(fimm.dss.t10$e.l4)

## path <- "ohsu.dss.t0.tsv"
synId <- "syn10083494"
obj <- synGet(synId, downloadFile = TRUE)
file <- getFileLocation(obj)
ohsu.dss.t0 <- as.data.frame(fread(file))
ohsu.dss.t0$ic50.ll4 <- remove.outliers(ohsu.dss.t0$e.ll4)
ohsu.dss.t0$ic50.l4 <- remove.outliers(ohsu.dss.t0$e.l4)

## path <- "ohsu.dss.t10.tsv"
synId <- "syn10083499"
obj <- synGet(synId, downloadFile = TRUE)
file <- getFileLocation(obj)
ohsu.dss.t10 <- as.data.frame(fread(file))
ohsu.dss.t10$ic50.ll4 <- remove.outliers(ohsu.dss.t10$e.ll4)
ohsu.dss.t10$ic50.l4 <- remove.outliers(ohsu.dss.t10$e.l4)

## path <- "ohsu.fimm.drugs.tsv"
synId <- "syn10083888"
obj <- synGet(synId, downloadFile = TRUE)
file <- getFileLocation(obj)
ohsu.fimm.drugs <- as.data.frame(fread(file))

## Plot the metrics computed in L.4 (logistic) vs LL.4 (log-logistic) fits.
plot.all.metrics <- function(dat.arg, title) {
  glist <- list()
  data <- dat.arg
  for(metric in c("ic50", "auc", "dss1", "dss2")) {
    x.col <- paste0(metric, ".l4")
    y.col <- paste0(metric, ".ll4")
    x.lab <- paste0(toupper(metric), " (L4)")
    y.lab <- paste0(toupper(metric), " (LL4)")
    data[, x.col] <- remove.outliers(data[, x.col])
    data[, y.col] <- remove.outliers(data[, y.col])
    g <- plot.scatter(data, x.col, y.col, x.lab, y.lab)
    num <- nrow(na.omit(data[,c(x.col, y.col)]))
    total.num <- nrow(data)
    g <- g + ggtitle(paste0("n = ", num, " (total = ", total.num, ")"))
    ## Add regression and identity lines
    g <- g + geom_smooth(method = "lm", se = FALSE)
    g <- g + geom_abline(intercept = 0, slope = 1, linetype = 2, color = "blue")
    glist[[length(glist)+1]] <- g
  }
  ## do.call("grid.arrange", c(glist, top = title))
  do.call("arrangeGrob", c(glist, top = title))
}

## Create the plots and store in synapse
## Store the plots in synapse in the following "Drug Response Fit QC" folder
parentId <- "syn10089446"

## Attach the results to this script in github via the 'executed' provenance
executed.url <- "https://github.com/bswhite/fimm-ohsu-aml/blob/master/061317/qc-fit.R"

do.syn.store <- TRUE

fname <- "fimm-dss-t0-metric-l4-vs-ll4-correlation.pdf"
pdf(fname, onefile = FALSE)
plt <- plot.all.metrics(fimm.dss.t0, "FIMM (t = 0)")
plot(plt)
d <- dev.off()
if(do.syn.store) {
  f <- File(fname, parentId = parentId, synapseStore = TRUE)
  synStore(f, executed = list(list(url=executed.url)), forceVersion = FALSE)
}
## mydoc <- addSlide( mydoc, slide.layout = "Title and Content" )
## title <- "FIMM (t = 0)"
## mydoc <- addTitle( mydoc, title )
## mydoc <- addPlot( doc = mydoc, fun = plot, x = plt)
gc()

fname <- "fimm-dss-t10-metric-l4-vs-ll4-correlation.pdf"
pdf(fname, onefile = FALSE)
plt <- plot.all.metrics(fimm.dss.t10, "FIMM (t = 10)")
plot(plt)
d <- dev.off()
if(do.syn.store) {
  f <- File(fname, parentId = parentId, synapseStore = TRUE)
  synStore(f, executed = list(list(url=executed.url)), forceVersion = FALSE)
}
## mydoc <- addSlide( mydoc, slide.layout = "Title and Content" )
## title <- "FIMM (t = 10)"
## mydoc <- addTitle( mydoc, title )
## mydoc <- addPlot( doc = mydoc, fun = plot, x = plt)
gc()

fname <- "ohsu-dss-t0-metric-l4-vs-ll4-correlation.pdf"
pdf(fname, onefile = FALSE)
plt <- plot.all.metrics(ohsu.dss.t0, "OHSU (t = 0)")
plot(plt)
d <- dev.off()
if(do.syn.store) {
  f <- File(fname, parentId = parentId, synapseStore = TRUE)
  synStore(f, executed = list(list(url=executed.url)), forceVersion = FALSE)
}
## mydoc <- addSlide( mydoc, slide.layout = "Title and Content" )
## title <- "OHSU (t = 0)"
## mydoc <- addTitle( mydoc, title )
## mydoc <- addPlot( doc = mydoc, fun = plot, x = plt)
gc()

fname <- "ohsu-dss-t10-metric-l4-vs-ll4-correlation.pdf"
pdf(fname, onefile = FALSE)
plt <- plot.all.metrics(ohsu.dss.t10, "OHSU (t = 10)")
plot(plt)
d <- dev.off()
if(do.syn.store) {
  f <- File(fname, parentId = parentId, synapseStore = TRUE)
  synStore(f, executed = list(list(url=executed.url)), forceVersion = FALSE)
}
## mydoc <- addSlide( mydoc, slide.layout = "Title and Content" )
## title <- "OHSU (t = 10)"
## mydoc <- addTitle( mydoc, title )
## mydoc <- addPlot( doc = mydoc, fun = plot, x = plt)
gc()

q(status=0)

## Plot correlation of drugs within a drug class for each metric

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

plot.ohsu.panel.of.pairs <- function(data, main = NULL) {
  gs <- list()
  for(metric in c("ic50", "auc", "dss1", "dss2")) {
    tmp <- data[,c("lab_id", "inhibitor", metric, "replicant", "time_of_read")]
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
      p <- p + ggtitle(toupper(metric))
    } else {
      p <- p + ggtitle(paste0(main, ":\n", toupper(metric)))
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
      p <- p + ggtitle(toupper(metric))
    } else {
      p <- p + ggtitle(paste0(main, ":\n", toupper(metric)))
    }
    gs[[metric]] <- p
  }
  ## do.call("grid.arrange", gs)
  return(do.call("arrangeGrob", list(grobs = gs, ncol = 2)))
}

## Plot drug-drug correlations as a function of metric and data set within each drug class

mechanisms <- as.data.frame(table(ohsu.fimm.drugs$`Mechanism/Targets`))
colnames(mechanisms)[1] <- "mechanism"
mechanisms <- mechanisms$mechanism[mechanisms$Freq > 1]
for(mechanism in mechanisms) {
  print(mechanism)
  mechanism.annotations <- ohsu.fimm.drugs[ohsu.fimm.drugs$`Mechanism/Targets` == mechanism,]

  ## Store the plots in synapse in the following "Drug Response Fit QC" folder
  parentId <- "syn10089446"

  for(dss.t in c(0, 10)) {
    ohsu.mechanism.data <- NULL
    if(dss.t == 0) {
      ohsu.mechanism.data <- subset(ohsu.dss.t0, inhibitor %in% mechanism.annotations$ID_Drug.ohsu)
    } else if(dss.t == 10) {
      ohsu.mechanism.data <- subset(ohsu.dss.t10, inhibitor %in% mechanism.annotations$ID_Drug.ohsu)
    }  
    ## mydoc <- addSlide( mydoc, slide.layout = "Title and Content" )
    mydoc <- addSlide( mydoc, slide.layout = "Two Content" )
    ## title <- paste0("OHSU (", ifelse(ll.type == "ll4", "LL.4", "L.4"), "; t = ", dss.t, "): ", mechanism)
    title <- paste0("OHSU (t = ", dss.t, "): ", mechanism)
    mydoc <- addTitle( mydoc, title )
    for(ll.type in c("ll4", "l4")) {
      fname <- paste0(make.names(mechanism), "-correlation-", ll.type, "-ohsu-t", dss.t, ".pdf")
      pdf(fname, onefile = FALSE)
      data <- ohsu.mechanism.data[,c("lab_id", "inhibitor", "replicant", "time_of_read", paste0("ic50.", ll.type), paste0("auc.", ll.type), paste0("dss1.", ll.type), paste0("dss2.", ll.type))]
      colnames(data) <- c("lab_id", "inhibitor", "replicant", "time_of_read", "ic50", "auc", "dss1", "dss2")
      plt <- plot.ohsu.panel.of.pairs(data, main = paste0("OHSU (", ifelse(ll.type == "ll4", "LL.4", "L.4"), "; t = ", dss.t, "): ", mechanism))
      plot(plt)
      d <- dev.off()

      if(do.syn.store) {
        f <- File(fname, parentId = parentId, synapseStore = TRUE)
        synStore(f, executed = list(list(url=executed.url)), forceVersion = FALSE)
      }

      ## mydoc <- addSlide( mydoc, slide.layout = "Title and Content" )
      ## title <- paste0("OHSU (", ifelse(ll.type == "ll4", "LL.4", "L.4"), "; t = ", dss.t, "): ", mechanism)
      ## mydoc <- addTitle( mydoc, title )
      if(ll.type == "ll4") {
##        mydoc <- addPlot( doc = mydoc, fun = plot, x = plt, offx = 0.5, offy = 2, width = 6, height = 5)
      } else {
##        mydoc <- addPlot( doc = mydoc, fun = plot, x = plt, offx = 7.5, offy = 2, width = 6, height = 5)
      }
      mydoc <- addPlot( doc = mydoc, fun = plot, x = plt )
      gc()
      
    }
  }

  for(dss.t in c(0, 10)) {
    fimm.mechanism.data <- NULL
    if(dss.t == 0) {
      fimm.mechanism.data <- subset(fimm.dss.t0, DRUG_ID %in% mechanism.annotations$ID_Drug)
    } else if(dss.t == 10) {
      fimm.mechanism.data <- subset(fimm.dss.t10, DRUG_ID %in% mechanism.annotations$ID_Drug)
    }  
    
    ## mydoc <- addSlide( mydoc, slide.layout = "Title and Content" )
    mydoc <- addSlide( mydoc, slide.layout = "Two Content" )
    ## title <- paste0("FIMM (", ifelse(ll.type == "ll4", "LL.4", "L.4"), "; t = ", dss.t, "): ", mechanism)
    title <- paste0("FIMM (t = ", dss.t, "): ", mechanism)
    mydoc <- addTitle( mydoc, title )
    for(ll.type in c("ll4", "l4")) {
      fname <- paste0(make.names(mechanism), "-correlation-", ll.type, "-fimm-t", dss.t, ".pdf")
      pdf(fname, onefile = FALSE)
      data <- fimm.mechanism.data[, c("SCREEN_ID", "DRUG_ID", paste0("ic50.", ll.type), paste0("auc.", ll.type), paste0("dss1.", ll.type), paste0("dss2.", ll.type))]
      colnames(data) <- c("SCREEN_ID", "DRUG_ID", "ic50", "auc", "dss1", "dss2")
      plt <- plot.fimm.panel.of.pairs(data, main = paste0("FIMM (", ifelse(ll.type == "ll4", "LL.4", "L.4"), "; t = ", dss.t, "): ", mechanism))
      plot(plt)
      d <- dev.off()

      if(do.syn.store) {
        f <- File(fname, parentId = parentId, synapseStore = TRUE)
        synStore(f, executed = list(list(url=executed.url)), forceVersion = FALSE)
      }

      ## mydoc <- addSlide( mydoc, slide.layout = "Title and Content" )
      ## title <- paste0("FIMM (", ifelse(ll.type == "ll4", "LL.4", "L.4"), "; t = ", dss.t, "): ", mechanism)
      ## mydoc <- addTitle( mydoc, title )
      ## mydoc <- addPlot( doc = mydoc, fun = plot, x = plt)
      if(ll.type == "ll4") {
##        mydoc <- addPlot( doc = mydoc, fun = plot, x = plt, offx = 0.5, offy = 2, width = 6, height = 5)
      } else {
##        mydoc <- addPlot( doc = mydoc, fun = plot, x = plt, offx = 7.5, offy = 2, width = 6, height = 5)
      }
      mydoc <- addPlot( doc = mydoc, fun = plot, x = plt )
      gc()
      
    }
  }


}

writeDoc(mydoc, "062017-ll4-vs-l4-qc.pptx")
