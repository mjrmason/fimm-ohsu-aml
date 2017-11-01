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

## Set the download and output directory
download.path <- "../input/"
if (!file.exists(download.path)) {
  dir.create(download.path)
}

output.path <- "output/"
if (!file.exists(output.path)) {
  dir.create(output.path)
}

beat.aml.cache <- "../input/beat-aml-files/"

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

my.dup <- function(x) {
  duplicated(x, fromLast=TRUE) | duplicated(x, fromLast=FALSE)
}

suppressPackageStartupMessages(library("scales")) ## for scientific

Sys.setenv(R_ZIPCMD= "/usr/bin/zip")  

suppressPackageStartupMessages(library("parallel"))
num.cores <- detectCores()
if(!is.na(num.cores) && (num.cores > 1)) {
  suppressPackageStartupMessages(library("doMC"))
  cat(paste("Registering ", num.cores-1, " cores.\n", sep=""))
  registerDoMC(cores=(num.cores-1))
}

cat("Recalculate OHSU and FIMM DSS values using L.4 and LL.4\n")
cat("Note: this enforces the overlap between the OHSU and FIMM concentration ranges\n")

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

## BEGIN load FIMM data

## Load (batch corrected) FIMM expression data
obj <- synGet(id="syn8270602", downloadFile = TRUE, downloadLocation = download.path)
file <- getFileLocation(obj)
##fimm.expr <- read.table(file, header=TRUE, sep=",")
fimm.expr <- fread.with.rownames(file, sep=",")
fimm.expr <- t(fimm.expr)

## Load FIMM metadata (including relapse/refractory/diagnosis)
synId <- "syn8270594"
obj <- synGet(id=synId, downloadFile = TRUE, downloadLocation = download.path)
file <- getFileLocation(obj)
fimm.metadata <- read.table(file, header=TRUE, sep=",")

## Load FIMM genomic data
synId <- "syn8270591"
obj <- synGet(id=synId, downloadFile = TRUE, downloadLocation = download.path)
file <- getFileLocation(obj)
fimm.genomic.detail <- read.table(file, header=TRUE, sep=",", comment.char="", quote="\"")

## The "bin" file simply is 0/1 matrix indicating whether the sample (row) has a mutation
## in the gene (column).  These use gene symbols.
synId <- "syn8270590"
obj <- synGet(id=synId, downloadFile = TRUE, downloadLocation = download.path)
file <- getFileLocation(obj)
fimm.genomic.bin <- read.table(file, header=TRUE, sep=",", comment.char="", quote="\"")
## Make the samples the columns
fimm.genomic.bin <- t(fimm.genomic.bin)

## Read in FIMM raw drug response data
synId <- "syn8488824"
obj <- synGet(id=synId, downloadFile = TRUE, downloadLocation = download.path)
file <- getFileLocation(obj)
fimm.raw.dr <- read.xlsx(file)

## Read in FIMM drug annotations
synId <- "syn8434109"
obj <- synGet(id=synId, downloadFile = TRUE, downloadLocation = download.path)
file <- getFileLocation(obj)
fimm.drug.annotations <- read.xlsx(file)

## Read in the subset of drug annotations in both FIMM and OHSU
synId <- "syn9731315"
obj <- synGet(id=synId, downloadFile = TRUE, downloadLocation = download.path)
ohsu.fimm.drugs <- openxlsx:::read.xlsx(getFileLocation(obj), sheet=1)

## END load FIMM data

## BEGIN load OHSU data

## Load the OHSU diagnostics -- the latest appears to have all patient info for all data releases
## Diagnosis_Labs_Treatments_Outcomes_2017_01_12.xlsx
synId <- "syn8149174"
obj <- synGet(id=synId, downloadFile = TRUE, downloadLocation = download.path)
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

## Read in the raw OHSU inhibitor data 
## inhibitor_data_points_2017_01_12.txt
## This latest file seems to have data from the prior releases
synId <- "syn8149180"
inh.obj <- synGet(synId, downloadFile=TRUE, downloadLocation = download.path)
ohsu.inh.tbl <- fread(getFileLocation(inh.obj))
ohsu.inh.tbl <- as.data.frame(ohsu.inh.tbl)

## Load OHSU expression data
load.ohsu.expr.data <- function() {
  ## Download the AML CPM expression data
  ## NB: linking to the original Beat AML project (syn2942337), since this has more updated files/releases
  ## (as well as the raw fastqs)
  ## Read in BeatAML_DR1_RNASeq_log2_cpm_2016_08_02.csv
  cat("Downloading and reading DR1 data\n")
  rna.dr1.obj <- synGet("syn7124177", downloadFile=TRUE, downloadLocation = download.path)
  # Load the data
  rna.dr1.tbl <- as.data.frame(fread(getFileLocation(rna.dr1.obj), sep=",", header=TRUE))
  rna.dr1.lab.ids <- gsub("\\.", "-", gsub("^X", "", names(rna.dr1.tbl[-c(1:9)])))

  cat("Downloading and reading DR2 Q1 data\n")
  ## Read in BeatAML_DR2q1_RNASeq_log2_cpm_2016_12_16.csv
  rna.dr2q1.obj <- synGet("syn8149220", downloadFile=TRUE, downloadLocation = download.path)
  rna.dr2q1.tbl <- as.data.frame(fread(getFileLocation(rna.dr2q1.obj), sep=",", header=TRUE))

  rna.dr2q1.lab.ids <- gsub("\\.", "-", gsub("^X", "", names(rna.dr2q1.tbl[-c(1:9)])))

  cat("Downloading and reading DR2 Q2 data\n")
  ## Read in BeatAML_DR2q2_RNASeq_log2_cpm_2017_01_11.csv
  rna.dr2q2.obj <- synGet("syn8149259", downloadFile=TRUE, downloadLocation = download.path)
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

  ## Ensure that DR2 Q1 and Q2 are non-overlapping.
  dr2q1.and.2.samples <- c(colnames(rna.dr2q1.tbl), colnames(rna.dr2q2.tbl))
  if(any(duplicated(dr2q1.and.2.samples))) {
    stop("Duplicated samples between DR2 Q1 and Q2!\n")
  } else {
    cat("As expected, no duplicated samples between DR2 Q1 and Q2\n")
  }

  all.samples <- c(colnames(rna.dr1.tbl), colnames(rna.dr2q1.tbl), colnames(rna.dr2q2.tbl))
  if(any(duplicated(all.samples))) {
    warning("Some RNA-seq samples are duplicated across the releases\n")
  } else {
    cat("As expected, no samples are duplicated across the releases\n")
  }
  ## Since no samples are duplicated, we can just merge the columns
  ohsu.expr <- cbind(rna.dr1.tbl, rna.dr2q1.tbl, rna.dr2q2.tbl)
  rm(rna.dr1.tbl)
  rm(rna.dr2q1.tbl)
  rm(rna.dr2q2.tbl)
  gc()
  ohsu.expr
}

ohsu.expr <- load.ohsu.expr.data()

## Read in the latest rna-seq dashboard
## BeatAML_rnaseq_2017_01_19_public_dashboard.xlsx (syn8149218)
## As stated in the README (datawave-2-q3_rnaseq_README.pdf),
## this is a cumulative dashboard that covers all current releases (i.e., DR1 and DR2)
synId <- "syn8149218"
rna.obj <- synGet(synId, downloadFile=TRUE, downloadLocation = download.path)
file <- getFileLocation(rna.obj)
## The samples summary is the 4th sheet
## This sample.summary table provides a map from patient id to sequence id.
## Except that since the sample sequence ids where column headers in our expression
## data, they were corrupted from AA-BBBBB (where A and B are digits) to XAA.BBBBB.
## So change that.
cat(paste0("Reading sample summary file ", file, "\n"))
rnaseq.sample.summary.tbl <- openxlsx:::read.xlsx(file, sheet=4)

## Read in the OHSU patient diagnosis table (corresponding to the drug response data).
## Diagnosis_Labs_Treatments_Outcomes_2017_01_12.xlsx (syn8149174)
synId <- "syn8149174"
obj <- synGet(synId, downloadFile=TRUE, downloadLocation = download.path)
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
## BeatAML_seqcap_2017_01_19_public_dashboard.xlsx
## As stated in the README (datawave-2_seqcap_README.pdf), this dashboard is
## cumulative for all of the current releases (i.e., DR1 and DR2)
synId <- "syn8149211"
dna.obj <- synGet(synId, downloadFile=TRUE, downloadLocation = download.path)
file <- getFileLocation(dna.obj)
cat(paste0("Reading sample summary file ", file, "\n"))
all.dnaseq.sample.summary.tbl <- openxlsx:::read.xlsx(file, sheet=4)
tumor.only.dnaseq.sample.summary.tbl <- openxlsx:::read.xlsx(file, sheet=5)
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
  for(i in 1:nrow(maf.tbl)){ maf.list[[maf.names[i]]] <-  read.maf(getFileLocation(synGet(maf.tbl$file.id[i], downloadLocation = download.path)), verbose = F); print(i)}
  
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
rm(fimm.expr)
rm(filtered.res); gc()

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

## Compute the AUC for a _log-logistic_ fit using LL.4/llogistic with drc
## All concentrations/IC50s are in real space
## NB: FIMM/Yadav fits a logistic function
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
  ## The interpretation of the slop/b/b' is different between the two
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

  val <- integrate(fnc, min.conc, max.conc)
  if(val$message == "OK") { auc <- val$value }
  if(!is.finite(auc)) {
    auc <- NA
  }
  return(auc)
}

## Compute the AUC for a _logistic_ fit using L.4/logistic with drc
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

## Compute DSS1, DSS2, and DSS3 for a _log-logistic_ fit using LL.4/llogistic with drc
## Return a vector holding:
## the AUC value used to calculate DSS ("dss.auc")
## and DSS1, DSS2, and DSS3 as "dss1", "dss2", and "dss3", respectively.
## NB: dss.auc is the DSS calculated between 
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
    vec <- c(NA, NA, NA, NA)
    names(vec) <- c("dss.auc", "dss1", "dss2", "dss3")
    return(vec)
  }

  ## Compute AUC to be used for DSS.  It should be calculated between x1 and max.conc.
  dss.auc <- tryCatch({compute.ll.4.auc(b.param, c.param, d.param, e.param, x1, max.conc)},
                      error = function(e) { NA })
  if(!is.finite(dss.auc)) {
    vec <- c(NA, NA, NA, NA)
    names(vec) <- c("dss.auc", "dss1", "dss2", "dss3")
    return(vec)
  }

  ## Compute DSS1, 2, and 3.
  dss <- compute.dss(dss.auc, t = t, c.min = x1, c.max = max.conc, x1 = x1, x2 = max.conc, top.asymptote = d.param)
  vec <- c(dss.auc, dss)
  names(vec) <- c("dss.auc", names(dss))
  return(vec)
}

## Compute DSS1, DSS2, and DSS3 for a _logistic_ fit using L.4/logistic with drc
## Return a vector holding:
## the AUC value used to calculate DSS ("dss.auc")
## and DSS1, DSS2, and DSS3 as "dss1", "dss2", and "dss3", respectively.
## NB: dss.auc is the DSS calculated between 
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
    vec <- c(NA, NA, NA, NA)
    names(vec) <- c("dss.auc", "dss1", "dss2", "dss3")
    return(vec)
  }

  ## Compute AUC to be used for DSS.  It should be calculated between x1 and max.conc.
  dss.auc <- tryCatch({compute.l.4.auc(b.param, c.param, d.param, e.param, x1, max.conc, numerical = numerical)},
                      error = function(e) { NA })
  if(!is.finite(dss.auc)) {
    vec <- c(NA, NA, NA, NA)
    names(vec) <- c("dss.auc", "dss1", "dss2", "dss3")
    return(vec)
  }

  ## Compute DSS1, 2, and 3.
  dss <- compute.dss(dss.auc, t = t, c.min = x1, c.max = max.conc, x1 = x1, x2 = max.conc, top.asymptote = d.param)
  vec <- c(dss.auc, dss)
  names(vec) <- c("dss.auc", names(dss))
  return(vec)
}

## Fit _log-logistic_ function (fct = "LL.4") of _logistic_ function (fct = "L.4") to OHSU data using DRC.
## Also, compute AUC.
## NB: OHSU data are viability; "invert" (subtract from 100%) to define inhbition.
fit.fct.to.ohsu <- function(data, fct = "LL.4") {
  ddply(data[,c("inhibitor", "well_concentration", "normalized_viability", "time_of_read", "lab_id", "replicant", "patient_id")],
        c("patient_id", "inhibitor", "lab_id", "replicant", "time_of_read"),
        .parallel = TRUE,
        .fun = function(df) {
          fit <- NULL
          ## NB: OHSU data are viability; "invert" (subtract from 100%) to define inhbition.
          df$inhibition <- 100 - df$normalized_viability
          if(fct == "LL.4") {
            tryCatch({
              fit <- suppressWarnings(drm(inhibition ~ well_concentration, data = df, fct=LL.4(fixed = c(NA, NA, NA,NA), names = c("SLOPE","MIN","MAX","IC50"))))
            }, 
            error = function(e) {
              cat("Trapped error in LL.4\n")
              NULL
            })
          } else {
            tryCatch({
              fit <- suppressWarnings(drm(inhibition ~ well_concentration, data = df, fct=L.4(fixed = c(NA, NA, NA,NA), names = c("SLOPE","MIN","MAX","IC50"))))
            }, 
            error = function(e) {
              cat("Trapped error in L.4\n")
              NULL
            })
          }
          min.conc <- min(df$well_concentration)
          max.conc <- max(df$well_concentration)
          gof <- 0
          auc <- NA
          ## Regardless of using L.4 and LL.4, the parameters have the following interpretation
          ## b.param: slope
          ## c.param: minimum asymptote
          ## d.param: maximum asymptote
          ## e.param: ic50 (not log ic50--even for LL.4)
          b.param <- 0
          c.param <- 0
          d.param <- 0
          e.param <- 0
          converged <- 0
          if(!is.null(fit) && (fit$fit$convergence || fit$convergence)) {
            converged <- 1
            ssres <- sum(residuals(fit)^2)
            resp <- na.omit(fit$dataList$resp)
            sstot <- sum((mean(resp)-resp)^2)
            gof <- 1 - (ssres/sstot)
            b.param <- coef(fit)[1]
            c.param <- coef(fit)[2]
            d.param <- coef(fit)[3]
            e.param <- coef(fit)[4]
            if(fct == "LL.4") {
              auc <- compute.ll.4.auc(b.param, c.param, d.param, e.param, min.conc, max.conc)
            } else {
              auc <- compute.l.4.auc(b.param, c.param, d.param, e.param, min.conc, max.conc, numerical = FALSE) 
            }
          }
          vec <- c(converged, b.param, c.param, d.param, e.param, gof, auc, min.conc, max.conc)
          names(vec) <- c("converged", "b", "c", "d", "e", "gof", "auc", "min.conc", "max.conc")
          vec
        })
  }

## Fit _log-logistic_ function (fct = "LL.4") of _logistic_ function (fct = "L.4") to FIMM data using DRC.
## Also, compute AUC.
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
          min.conc <- min(df$CONCENTRATION)
          max.conc <- max(df$CONCENTRATION)
          gof <- 0
          auc <- NA
          ## Regardless of using L.4 and LL.4, the parameters have the following interpretation
          ## b.param: slope
          ## c.param: minimum asymptote
          ## d.param: maximum asymptote
          ## e.param: ic50 (not log ic50--even for LL.4)
          b.param <- 0
          c.param <- 0
          d.param <- 0
          e.param <- 0
          converged <- 0
          if(!is.null(fit) && (fit$fit$convergence || fit$convergence)) {
            converged <- 1
            ssres <- sum(residuals(fit)^2)
            resp <- na.omit(fit$dataList$resp)
            sstot <- sum((mean(resp)-resp)^2)
            gof <- 1 - (ssres/sstot)
            b.param <- coef(fit)[1]
            c.param <- coef(fit)[2]
            d.param <- coef(fit)[3]
            e.param <- coef(fit)[4]
            if(fct == "LL.4") {
              auc <- compute.ll.4.auc(b.param, c.param, d.param, e.param, min.conc, max.conc)
            } else {
              auc <- compute.l.4.auc(b.param, c.param, d.param, e.param, min.conc, max.conc, numerical = FALSE) 
            }
          }
          vec <- c(converged, b.param, c.param, d.param, e.param, gof, auc, min.conc, max.conc)
          names(vec) <- c("converged", "b", "c", "d", "e", "gof", "auc", "min.conc", "max.conc")
          vec
        })
  }

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

compute.all.auc <- function(data, fct = "LL.4") {
  res <- ddply(data,
               colnames(data), 
               .parallel = TRUE,
               .fun = function(df) {
                 if(nrow(df) != 1) { stop("Expected to process one row\n") }
                 converged <- df$converged[1]
                 auc <- NA
                 if(converged == 1) {
                   b.param <- df$b[1]
                   c.param <- df$c[1]
                   d.param <- df$d[1]
                   e.param <- df$e[1]
                   min.conc <- df$min.conc[1]
                   max.conc <- df$max.conc[1]
                   auc <- tryCatch({
                              if(fct == "LL.4") {
                                compute.ll.4.auc(b.param, c.param, d.param, e.param, min.conc, max.conc)
                              } else {
                                compute.l.4.auc(b.param, c.param, d.param, e.param, min.conc, max.conc, numerical = FALSE) 
                              }
                            }, 
                            error = function(e) {
                              NA
                            })
                   }
                   auc
                 })
  colnames(res)[length(colnames(res))] <- "auc"
  res
}

## TODO
## -- roll back to old data
## -- auc -> dss drop?
## -- auc vs gene + drug (over-represented drugs?  or patients?)


compute.all.dss <- function(data, t = 0, fct = "LL.4") {
  res <- ddply(data,
               colnames(data), 
               .parallel = FALSE,
               .fun = function(df) {
                 if(nrow(df) != 1) { stop("Expected to process one row\n") }
                 converged <- df$converged[1]
                 vec <- c(NA, NA, NA, NA)
                 names(vec) <- c("dss.auc", "dss1", "dss2", "dss3")
                 if(converged == 1) {
                   b.param <- df$b[1]
                   c.param <- df$c[1]
                   d.param <- df$d[1]
                   e.param <- df$e[1]
                   min.conc <- df$min.conc[1]
                   max.conc <- df$max.conc[1]
                   vec <- tryCatch({
                              if(fct == "LL.4") {
                                compute.ll.4.dss(b.param, c.param, d.param, e.param, t = t, min.conc, max.conc)
                              } else {
                                compute.l.4.dss(b.param, c.param, d.param, e.param, t = t, min.conc, max.conc, numerical = FALSE) 
                              }
                            }, 
                            error = function(e) {
                              vec <- c(NA, NA, NA, NA)
                              names(vec) <- c("dss.auc", "dss1", "dss2", "dss3")
                              vec    
                            })
                   }
                   vec
                 })
  res
}

## Calculate fits

cat("Fitting LL.4 to OHSU data\n")
ll4.ohsu.fits <- fit.fct.to.ohsu(ohsu.raw.drug.data, fct = "LL.4") 
cat("Done fitting LL.4 to OHSU data\n\n")

cat("Fitting L.4 to OHSU data\n")
l4.ohsu.fits <- fit.fct.to.ohsu(ohsu.raw.drug.data, fct = "L.4") 
cat("Done fitting L.4 to OHSU data\n\n")

save.image(".RData")

cat("Fitting LL.4 to FIMM data\n")
ll4.fimm.fits <- fit.fct.to.fimm(fimm.raw.drug.data, fct = "LL.4") 
cat("Done fitting LL.4 to FIMM data\n\n")

cat("Fitting L.4 to FIMM data\n")
l4.fimm.fits <- fit.fct.to.fimm(fimm.raw.drug.data, fct = "L.4") 
cat("Done fitting L.4 to FIMM data\n\n")

save.image(".RData")

## Calculate AUC
## This is already calculated above in the fit.fct.to.* functions
## ll4.ohsu.auc <- compute.all.auc(ll4.ohsu.fits, fct = "LL.4") 
## l4.ohsu.auc <- compute.all.auc(l4.ohsu.fits, fct = "L.4") 

## ll4.fimm.auc <- compute.all.auc(ll4.fimm.fits, fct = "LL.4") 
## l4.fimm.auc <- compute.all.auc(l4.fimm.fits, fct = "L.4") 

## Calculate DSS
cat("Computing LL.4 OHSU DSS\n")
ll4.ohsu.dss <- compute.all.dss(ll4.ohsu.fits, t = 0, fct = "LL.4") 

cat("Computing L.4 OHSU DSS\n")
l4.ohsu.dss <- compute.all.dss(l4.ohsu.fits, t = 0, fct = "L.4") 

cat("Computing LL.4 FIMM DSS\n")
ll4.fimm.dss <- compute.all.dss(ll4.fimm.fits, t = 0, fct = "LL.4") 

cat("Computing L.4 FIMM DSS\n")
l4.fimm.dss <- compute.all.dss(l4.fimm.fits, t = 0, fct = "L.4") 

save.image(".RData")

cat("Computing LL.4 OHSU DSS t = 10\n")
ll4.ohsu.dss.t10 <- compute.all.dss(ll4.ohsu.fits, t = 10, fct = "LL.4") 

cat("Computing L.4 OHSU DSS t = 10\n")
l4.ohsu.dss.t10 <- compute.all.dss(l4.ohsu.fits, t = 10, fct = "L.4") 

cat("Computing LL.4 FIMM DSS t = 10\n")
ll4.fimm.dss.t10 <- compute.all.dss(ll4.fimm.fits, t = 10, fct = "LL.4") 

cat("Computing L.4 FIMM DSS t = 10\n")
l4.fimm.dss.t10 <- compute.all.dss(l4.fimm.fits, t = 10, fct = "L.4") 

save.image(".RData")

## Merge the FIMM and OHSU results
common.ohsu.cols <- c("patient_id", "inhibitor", "lab_id", "replicant", "time_of_read", "min.conc", "max.conc")
ohsu.dss.t0 <- merge(ll4.ohsu.dss, l4.ohsu.dss, by = common.ohsu.cols, suffixes = c(".ll4", ".l4"), all = TRUE)
ohsu.dss.t10 <- merge(ll4.ohsu.dss.t10, l4.ohsu.dss.t10, by = common.ohsu.cols, suffixes = c(".ll4", ".l4"), all = TRUE)

common.fimm.cols <- c("SCREEN_ID", "DRUG_ID", "min.conc", "max.conc")
fimm.dss.t0 <- merge(ll4.fimm.dss, l4.fimm.dss, by = common.fimm.cols, suffixes = c(".ll4", ".l4"), all = TRUE)
fimm.dss.t10 <- merge(ll4.fimm.dss.t10, l4.fimm.dss.t10, by = common.fimm.cols, suffixes = c(".ll4", ".l4"), all = TRUE)

## Store the fits, etc in synapse and link them to this script in github

## Store the FIMM fits in FIMM Data/Processed Data (syn8270577)
parentId <- "syn8270577"
executed.url <- "https://github.com/bswhite/fimm-ohsu-aml/blob/master/061317/fit-ohsu-and-fimm-drug-response-data.R"

path <- "fimm.dss.t0.tsv"
write.table(file=path, fimm.dss.t0, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
f <- File(path, parentId = parentId, synapseStore = TRUE)
synStore(f, executed = list(list(url=executed.url)), forceVersion = FALSE)

path <- "fimm.dss.t10.tsv"
write.table(file=path, fimm.dss.t10, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
f <- File(path, parentId = parentId, synapseStore = TRUE)
synStore(f, executed = list(list(url=executed.url)), forceVersion = FALSE)


## Store the OHSU fits and expression in BEAT AML Data/Processed Data (syn10083332)
parentId <- "syn10083332"

path <- "ohsu.dss.t0.tsv"
write.table(file=path, ohsu.dss.t0, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
f <- File(path, parentId = parentId, synapseStore = TRUE)
synStore(f, executed = list(list(url=executed.url)), forceVersion = FALSE)

path <- "ohsu.dss.t10.tsv"
write.table(file=path, ohsu.dss.t10, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
f <- File(path, parentId = parentId, synapseStore = TRUE)
synStore(f, executed = list(list(url=executed.url)), forceVersion = FALSE)
 
path <- "ohsu.expr.tsv"
write.table(file=path, ohsu.expr, row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)
f <- File(path, parentId = parentId, synapseStore = TRUE)
synStore(f, executed = list(list(url=executed.url)), forceVersion = FALSE)

## Store the RNA- and DNA-seq sample tables to synapse (note that these are just XLS sheets)
path <- "ohsu.dnaseq.sample.summary.tsv"
write.table(file=path, dnaseq.sample.summary.tbl, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
f <- File(path, parentId = parentId, synapseStore = TRUE)
synStore(f, executed = list(list(url=executed.url)), forceVersion = FALSE)

path <- "ohsu.rnaseq.sample.summary.tsv"
write.table(file=path, rnaseq.sample.summary.tbl, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
f <- File(path, parentId = parentId, synapseStore = TRUE)
synStore(f, executed = list(list(url=executed.url)), forceVersion = FALSE)

path <- "ohsu.fimm.drugs.tsv"
write.table(file=path, ohsu.fimm.drugs, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
f <- File(path, parentId = parentId, synapseStore = TRUE)
synStore(f, executed = list(list(url=executed.url)), forceVersion = FALSE)

stop("stop")
