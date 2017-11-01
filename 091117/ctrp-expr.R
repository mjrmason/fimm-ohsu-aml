## Read in the CTRP/CCLE expression (count!) data
synId <- "syn5616077"
obj <- synGet(synId, downloadFile = TRUE)
file <- getFileLocation(obj)
ctrp.cnt.expr <- read.table(file, header=TRUE, sep=",")
rownames(ctrp.cnt.expr) <- ctrp.cnt.expr[,1]
ctrp.cnt.expr <- ctrp.cnt.expr[,-1]
