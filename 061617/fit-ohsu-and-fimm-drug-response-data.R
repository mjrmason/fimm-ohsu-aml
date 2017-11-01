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

num.cores <- detectCores()
if(!is.na(num.cores) && (num.cores > 1)) {
  suppressPackageStartupMessages(library("doMC"))
  cat(paste("Registering ", num.cores-1, " cores.\n", sep=""))
  registerDoMC(cores=(num.cores-1))
}

cat("Fit OHSU and FIMM drug response data using L.4 and LL.4\n")
cat("Calculate IC50, AUC, DSS1, DSS2, and DSS3 and upload to Synapse\n")
cat("Create scatter plot of L.4 vs LL.4 for each metric (IC50, AUC, DSS1, DSS2, and DSS3)\n")

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
fimm.expr <- read.table(file, header=TRUE, sep=",")
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
## Diagnosis_Labs_Treatments_Outcomes_2017_05_17.xlsx
obj <- synGet(id="syn10008777", downloadFile = TRUE, downloadLocation = download.path)
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

## Read in the raw OHSU inhibitor data inhibitor_data_points_2017_05_17.txt
## This latest file seems to have data from the prior releases
inh.obj <- synGet("syn10008778", downloadFile=TRUE, downloadLocation = download.path)
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

  cat("Downloading and reading DR2 rnaseq data\n")
  rna.dr2.obj <- synGet("syn10008812", downloadFile=TRUE, downloadLocation = download.path)
  rna.dr2.tbl <- as.data.frame(fread(getFileLocation(rna.dr2.obj), sep=",", header=TRUE))
  rna.dr2.lab.ids <- gsub("\\.", "-", gsub("^X", "", names(rna.dr2.tbl[-c(1:9)])))

  ## Make sure we have a common set of genes across all 3 expression data sets
  common.genes <- intersect(rna.dr1.tbl$Gene, intersect(rna.dr2q1.tbl$Gene, intersect(rna.dr2q2.tbl$Gene, rna.dr2.tbl$Gene)))
  rownames(rna.dr1.tbl) <- rna.dr1.tbl$Gene
  rownames(rna.dr2q1.tbl) <- rna.dr2q1.tbl$Gene
  rownames(rna.dr2q2.tbl) <- rna.dr2q2.tbl$Gene
  rownames(rna.dr2.tbl) <- rna.dr2.tbl$Gene
  rna.dr1.tbl <- rna.dr1.tbl[common.genes,]
  rna.dr2q1.tbl <- rna.dr2q1.tbl[common.genes,]
  rna.dr2q2.tbl <- rna.dr2q2.tbl[common.genes,]
  rna.dr2.tbl <- rna.dr2.tbl[common.genes,]

  ## Drop the annotation columns from the expression matrices
  rna.dr1.tbl <- rna.dr1.tbl[,-c(1:9)]
  rna.dr2q1.tbl <- rna.dr2q1.tbl[,-c(1:9)]
  rna.dr2q2.tbl <- rna.dr2q2.tbl[,-c(1:9)]
  rna.dr2.tbl <- rna.dr2.tbl[,-c(1:9)]
  gc()

  ## Ensure that DR2 Q1 and Q2 are non-overlapping.
  dr2q1.and.2.samples <- c(colnames(rna.dr2q1.tbl), colnames(rna.dr2q2.tbl))
  if(any(duplicated(dr2q1.and.2.samples))) {
    stop("Duplicated samples between DR2 Q1 and Q2!\n")
  } else {
    cat("As expected, no duplicated samples between DR2 Q1 and Q2\n")
  }

  ## DR2 is not a strict superset of DR2 Q1 and Q2:
  cat(paste0(length(which(!(dr2q1.and.2.samples %in% colnames(rna.dr2.tbl)))), " of ", length(dr2q1.and.2.samples), " samples from DR2 Q1 and DR2 Q2 are not in DR2\n"))
  ## Presumably those few missing from DR2 were excluded for a reason.  
  ## As described in the README (datawave-2-q3_rnaseq_README.pdf):
  ## "Note that the exact log2(cpm) values will differ if samples are added or removed 
  ## due to the scaling of the prior.count value."
  ## Hence, the values in DR2 Q1 and Q2 for sample X is slightly different than for that
  ## same sample in DR2.  Let's just use DR2.

  ## Ensure that there is no overlap between DR1 and DR2.
  rm(rna.dr2q1.tbl)
  rm(rna.dr2q2.tbl)
  gc()

  all.samples <- c(colnames(rna.dr1.tbl), colnames(rna.dr2.tbl))
  if(any(duplicated(all.samples))) {
    warning("Some RNA-seq samples are duplicated across the 2 releases\n")
  } else {
    cat("As expected, no samples are duplicated across the 2 releases\n")
  }
  ## Since no samples are duplicated, we can just merge the columns
  ohsu.expr <- cbind(rna.dr1.tbl, rna.dr2.tbl)
  rm(rna.dr1.tbl)
  rm(rna.dr2.tbl)
  gc()
  ohsu.expr
}

ohsu.expr <- load.ohsu.expr.data()

## working
## this MD5 is broken
## Read in the latest rna-seq dashboard
## BeatAML_rnaseq_2017_06_09_public_dashboard.xlsx (syn10008793)
## As stated in the README (datawave-2-q3_rnaseq_README.pdf),
## this is a cumulative dashboard that covers all current releases (i.e., DR1 and DR2)
## rna.obj <- synGet("syn10008793", downloadFile=TRUE, downloadLocation = download.path)
## file <- getFileLocation(rna.obj)
file <- paste0(beat.aml.cache, "BeatAML_rnaseq_2017_06_09_public_dashboard.xlsx")
## The samples summary is the 4th sheet
## This sample.summary table provides a map from patient id to sequence id.
## Except that since the sample sequence ids where column headers in our expression
## data, they were corrupted from AA-BBBBB (where A and B are digits) to XAA.BBBBB.
## So change that.
cat(paste0("Reading sample summary file ", file, "\n"))
rnaseq.sample.summary.tbl <- openxlsx:::read.xlsx(file, sheet=4)

## Read in the OHSU patient diagnosis table (corresponding to the drug response data).
## Diagnosis_Labs_Treatments_Outcomes_2017_05_17.xlsx (syn10008777)
obj <- synGet("syn10008777", downloadFile=TRUE, downloadLocation = download.path)
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
## BeatAML_seqcap_2017_06_09_public_dashboard.xlsx (syn10008791)
## As stated in the README (datawave-2_seqcap_README.pdf), this dashboard is
## cumulative for all of the current releases (i.e., DR1 and DR2)
## working
## This MD5 is broken
## dna.obj <- synGet("syn10008791", downloadFile=TRUE, downloadLocation = download.path)
## file <- getFileLocation(dna.obj)
file <- paste0(beat.aml.cache, "BeatAML_seqcap_2017_06_09_public_dashboard.xlsx")
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
rm(ohsu.expr)
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

stop("stop")


mechanisms <- as.data.frame(table(ohsu.fimm.drugs$`Mechanism/Targets`))
colnames(mechanisms)[1] <- "mechanism"
mechanisms <- mechanisms$mechanism[mechanisms$Freq > 1]
for(mechanism in mechanisms) {
  print(mechanism)
  mechanism.annotations <- ohsu.fimm.drugs[ohsu.fimm.drugs$`Mechanism/Targets` == mechanism,]

  ohsu.mechanism.data <- subset(ll4.fits.dss.fixed, inhibitor %in% mechanism.annotations$ID_Drug.ohsu)
  fname <- paste0(make.names(mechanism), "-correlation-ohsu-restricted-ic50.pdf")
  pdf(fname, onefile = FALSE)
  plt <- plot.ohsu.panel.of.pairs(ohsu.mechanism.data, main = paste0("OHSU: ", mechanism))
  ## ggsave(fname, plt)
  plot(plt)
  d <- dev.off()

  ohsu.mechanism.data <- subset(ll4.fits.dss.no.force.fixed, inhibitor %in% mechanism.annotations$ID_Drug.ohsu)
  fname <- paste0(make.names(mechanism), "-correlation-ohsu-unrestricted-ic50.pdf")
  pdf(fname, onefile = FALSE)
  plt <- plot.ohsu.panel.of.pairs(ohsu.mechanism.data, main = paste0("OHSU ", mechanism, " (unrestricted IC50)"))
  ## ggsave(fname, plt)
  plot(plt)
  d <- dev.off()
}


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
ll4.dss.fits.no.force.seq <- merge(ll4.fits.dss.no.force.fixed, unique(rnaseq.sample.summary.tbl[,c("SeqID", "Original_LabID")]), by.x = "lab_id", by.y = "Original_LabID")
ll4.dss.fits.no.force.seq <- subset(ll4.dss.fits.no.force.seq, gof > 0)

ll4.dss.fits.no.force.seq <- merge(ll4.fits.dss.no.force.fixed, unique(rnaseq.sample.summary.tbl[,c("SeqID", "Original_LabID")]), by.x = "lab_id", by.y = "Original_LabID")
ll4.dss.fits.no.force.seq <- subset(ll4.dss.fits.no.force.seq, gof > 0)

## Do each inhibitor separately.

## common.genes <- intersect(rownames(beat.mafs), rownames(fimm.genomic.bin))
common.genes <- intersect(rownames(fimm.expr.data), rownames(ohsu.expr.data))
mek.annotations <- ohsu.fimm.drugs[ohsu.fimm.drugs$`Mechanism/Targets` == "MEK1/2 inhibitor",]
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

## WORKING

m <- merge(fimm.mek.expr.all.auc.corr, ohsu.mek.expr.all.auc.corr, by = "gene", suffixes = c(".fimm", ".ohsu"))
ras.sig <- subset(m, (as.numeric(pval.fimm) < 0.05) & (as.numeric(pval.ohsu) < 0.05))
## ras.sig <- subset(m, (in.ras.sig == 0) & (as.numeric(pval.ohsu) < 0.05))
ras.sig <- ras.sig[order(as.numeric(ras.sig$pval.ohsu)),]
head(ras.sig[order(as.numeric(ras.sig$pval.ohsu)),c("gene","pval.ohsu","r2.ohsu")])
## ohsu.mek.expr.all.auc.corr <- correlate.expression.with.all.drug.response(drug.df = ll4.dss.fits.no.force.seq, drug.id.col = "inhibitor", genes.to.test = common.genes, drugs.to.test = mek.annotations$ID_Drug.ohsu, expr.mat = ohsu.expr.data, specimen.id.col = "SeqID", patient.id.col = "patient_id", metric.col = "auc")
## fimm.mek.expr.all.auc.corr <- correlate.expression.with.all.drug.response(drug.df = fimm.ll4.dss.fits.no.force.pt, drug.id.col = "DRUG_ID", genes.to.test = common.genes, drugs.to.test = mek.annotations$ID_Drug, expr.mat = fimm.expr.data, specimen.id.col = "SCREEN_ID", patient.id.col = "PATIENT_ID", metric.col = "auc")

# working
for(i in 1:5) {
  gene <- ras.sig$gene[i]
  pdf(paste0("fimm-ohsu-ras-", gene, "-vs-auc.pdf"))
  ##par(mfrow=c(1,2))
  t.gene <- data.frame(gene = unname(t(ohsu.expr.data[gene,])), seqid = colnames(ohsu.expr.data))
  d.tmp <- subset(ll4.dss.fits.no.force.seq, inhibitor %in% mek.annotations$ID_Drug.ohsu)
  m.tmp <- merge(d.tmp, t.gene, by.x = "SeqID", by.y = "seqid")
  m.tmp$auc2 <- ( 100 * (m.tmp$max.conc - m.tmp$min.conc ) ) - m.tmp$auc
  fit <- lm(m.tmp$auc ~ m.tmp$gene + m.tmp$inhibitor)
  rmse <- round(sqrt(mean(resid(fit)^2)), 2)
  r2 <- round(summary(fit)$r.squared, 2)
  eqn <- bquote(r^2 == .(r2) * "," ~~ RMSE == .(rmse))
  plot(m.tmp$gene, m.tmp$auc, xlab = paste0(gene, " Expr"), ylab = "AUC", main = "OHSU")
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

## ras.sig <- subset(m, (in.ras.sig == 0) & (as.numeric(pval.fimm) < 0.05) & (as.numeric(pval.ohsu) < 0.05))
## ras.sig <- subset(m, (in.ras.sig == 0) & (as.numeric(pval.ohsu) < 0.05))
ras.sig <- subset(m, (as.numeric(pval.fimm) < 0.05) & (as.numeric(pval.ohsu) < 0.05))
ras.sig <- ras.sig[order(as.numeric(ras.sig$pval.ohsu)),]
## ohsu.mek.expr.all.auc.corr <- correlate.expression.with.all.drug.response(drug.df = ll4.dss.fits.no.force.seq, drug.id.col = "inhibitor", genes.to.test = common.genes, drugs.to.test = mek.annotations$ID_Drug.ohsu, expr.mat = ohsu.expr.data, specimen.id.col = "SeqID", patient.id.col = "patient_id", metric.col = "auc")
## fimm.mek.expr.all.auc.corr <- correlate.expression.with.all.drug.response(drug.df = fimm.ll4.dss.fits.no.force.pt, drug.id.col = "DRUG_ID", genes.to.test = common.genes, drugs.to.test = mek.annotations$ID_Drug, expr.mat = fimm.expr.data, specimen.id.col = "SCREEN_ID", patient.id.col = "PATIENT_ID", metric.col = "auc")
for(i in 1:5) {
  gene <- ras.sig$gene[i]
  sym <- ensg.to.symbols.mapping(gene)$SYMBOL[1]
  pdf(paste0("fimm-ohsu-ras-", sym, "-vs-auc.pdf"))
  ##par(mfrow=c(1,2))
  glist1 <- list()
  glist2 <- list()
  
  t.gene <- data.frame(gene = unname(t(ohsu.expr.data[gene,])), seqid = colnames(ohsu.expr.data))
  d.tmp <- subset(ll4.dss.fits.no.force.seq, inhibitor %in% mek.annotations$ID_Drug.ohsu)
  m.tmp <- merge(d.tmp, t.gene, by.x = "SeqID", by.y = "seqid")
  m.tmp$auc2 <- ( 100 * (m.tmp$max.conc - m.tmp$min.conc ) ) - m.tmp$auc
  drug.col <- "inhibitor"
  df <- m.tmp[,c("gene", "auc", drug.col)]
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
  df <- m.tmp[,c("gene", "auc", drug.col)]
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

## working
synId <- "syn7440451"
obj <- synGet(synId, downloadFile=TRUE, downloadLocation = download.path)
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
## - do pca of patients based on auc
## - do pca of drugs based on auc
## - missing heatmap--but labeling according to missing, filtered, or present
## - for each patients, drugs
## for each clinical variable
## cluster and color based on clinical
## - repeat above but using pca

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

