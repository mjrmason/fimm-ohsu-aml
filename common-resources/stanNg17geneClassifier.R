rm(list=ls())

library(plyr)
library(tidyr)
library(data.table)
library(oligo)
library(doMC)
library(githubr) 
library(survival)
library(ggplot2)
library(GGally)
library(survcomp)
library(rSML)
library(powerMediation)
library(powerSurvEpi)
library(illuminaHumanv3.db)
library(illuminaHumanv4.db)
library(mice)
library(org.Hs.eg.db)
detach("package:synapseClient" , unload = T)


library(synapser)
synLogin()



clin <- read.csv(synGet('syn12046306')$path)
oldClin <- synTableQuery('SELECT * FROM syn5481721')$asDataFrame()

if(F)
{
  fLoc <- synGet("syn2355474")$path
  
  expr   <- read.csv(fLoc, sep =" ", row.names=1,header = T); # read in orig data
  genes  <- as.data.frame(illuminaHumanv4ENTREZID); 
  expr   <- expr[rownames(expr) %in% genes$probe_id,]
  probes <- rownames(expr)
  samps  <- colnames(expr)
  expr   <- complete(mice(expr))                              # impute missing data from 7 probes
  expr   <- preprocessCore::normalize.quantiles.robust(x=as.matrix(expr))
  
  rownames(expr) <- probes
  colnames(expr) <- samps
  
  write.csv(expr,file="~/SageOutput/AML/CARDINALimputedQuantilenormalizedProbelevelExpression.csv")
  tFile <- File(path="~/SageOutput/AML/CARDINALimputedQuantilenormalizedProbelevelExpression.csv", parent = "syn2354332")
  synStore(tFile)
  
  genes <- genes[match(rownames(expr), genes$probe_id),]
  tempEx <- data.frame(genes$gene_id,expr)
  
  compressed <- ddply(tempEx, .variables = "genes.gene_id", function(x){ apply(x[,-1],2,function(y){exp(mean(log(y)))})})
  rownames(compressed) <- compressed[,1]; compressed <- compressed[,-1]
  
  write.csv(compressed,file="~/SageOutput/AML/CARDINALimputedQuantilenormalizedEntrezIdExpression.csv")
  tFile <- File(path="~/SageOutput/AML/CARDINALimputedQuantilenormalizedEntrezIdExpression.csv", parent = "syn2354332")
  synStore(tFile)
}else{
  compressed <- read.csv(synGet('syn12046020')$path,row.names=1); colnames(compressed) <- gsub("^X","",colnames(compressed))
}

plot(density(compressed[,1], na.rm=T))
for(i in seq(2,ncol(compressed),5) ){lines(density(compressed[,i],na.rm=T))}

compressed <- compressed[,names(compressed)%in% clin$acquisition_id]
clin       <- clin[match(names(compressed), clin$acquisition_id),]
oldClin    <- oldClin[ match(names(compressed),oldClin$PatientID),]

genes <- as.data.frame(org.Hs.egSYMBOL)
LSC17coef <- data.frame(gene=c("DNMT3B","ZBTB46","NYNRIN","ARHGAP22","LAPTM4B","MMRN1","DPYSL3","FAM30A","CDK6","CPXM1","SOCS2","SMIM24","EMP1","BEX3","CD34","AKR1C3","ADGRG1"),beta=c(0.0874,-0.0347,0.00865,-0.0138,0.00582,0.0258,0.0284,0.0196,-0.0704,-0.0258,0.0271,-0.0226,0.0146,0.0465,0.0338,-0.0402,0.0501))
LSC17coef$entrez <- genes$gene_id[match(LSC17coef$gene,genes$symbol)]


tempExp <- t(scale(t(compressed[match(LSC17coef$entrez,rownames(compressed)),])))

lsc17_1 <- apply(compressed[match(LSC17coef$entrez,rownames(compressed)),],2,function(x){sum(x*LSC17coef$beta)})
lsc17_2 <- apply(tempExp,2,function(x){sum(x*LSC17coef$beta)})
lsc17_3 <- apply(compressed[match(LSC17coef$entrez,rownames(compressed)),],2,function(x){sum(log2(x)*LSC17coef$beta)})


pdf("~/SageOutput/AML/Ng17sigboxplot.pdf")
boxplot(lsc17_2[oldClin$response_class != "NA"]~oldClin$response_class[oldClin$response_class != "NA"])
beeswarm::beeswarm(lsc17_2[oldClin$response_class != "NA"]~oldClin$response_class[oldClin$response_class != "NA"], add=T,col="slategrey", pch =19)
pVal <- t.test(lsc17_2[oldClin$response_class != "NA"]~oldClin$response_class[oldClin$response_class != "NA"])$p.value
title(paste("CARDINAL Ng17 sig score\np.value =", signif(pVal,3)))


x = oldClin$response_class[oldClin$response_class != "NA"]
x = factor(x, levels = c("responder","non-responder"))
pred1 <- as.vector(lsc17_1[oldClin$response_class != "NA"])
pred2 <- as.vector(lsc17_2[oldClin$response_class != "NA"])
pred3 <- as.vector(lsc17_3[oldClin$response_class != "NA"])

roc1 <- roc(x~pred1)
roc2 <- roc(x~pred2)
roc3 <- roc(x~pred3)

plot(roc1)
plot(roc2, col = 2, add=T)
plot(roc3, col = 4, add=T)
legend("bottomright", legend=c("expression","z-scored norm", "log2 exp"), fill = c("black","red","blue"), border = NA, bty = "n",title= "CARDINAL\nROC Curves")
dev.off()

survfit(Surv(clin$survivial_time, clin$))








