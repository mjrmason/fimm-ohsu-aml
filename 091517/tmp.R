intraCors = llply(.data=pubEntrezExprDat, .fun = function(x){cor(apply(x,2,rank))}, .parallel = T)
plot(density(as.dist(as.matrix(intraCors[[1]]))))
for(i in 2:length(intraCors)){lines(density(as.dist(as.matrix(intraCors[[i]]))), col = i)}

png("~/SageOutput/CelgeneMM/intraStudyCors.png")
plot(density(fisherz(unlist(as.dist(intraCors[[1]]))),na.rm=T), main = "Intra Study Correlations", xlab="Fisher z transform of\nPearson Correlations", ylim= c(0,2.5),xlim=c(-1,5))
legend(2.5,2.5,legend=labels, fill=1:7,border=1:7, box.col="white")
for(i in 2:ncol(corMat)){lines(density(fisherz(unlist(as.dist(intraCors[[i]]))),na.rm=T), col = i); print(i)}
dev.off()


chipTypes              <- list(EMTAB4032= "hugene10", GSE19784HOVON65 = "hgu133plus2",  GSE24080UAMS = "hgu133plus2", GSE15695 = "hgu133plus2",EMTAB372 = "hgu133plus2", MMRF = "RNAseq", DFCI="RNAseq", M2Gen="RNAseq") 
colnames(corMat)       <- names(chipTypes)
corMat[corMat == Inf]  <- NA
cors = cor(corMat, use = "pairwise")
#write_feather(as.data.frame(corMat), "~/SageOutput/CelgeneMM/corMat.feather")
rm(corMat);gc()

cors

write.table(cors, "~/SageOutput/CelgeneMM/corOfCors.txt",sep="\t", col.names=T, row.names=T,quote=F)