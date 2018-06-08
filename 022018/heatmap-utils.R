## Get the OHSU covariates
synId <- "syn8149174"
obj <- synGet(id=synId, downloadFile = TRUE)
file <- getFileLocation(obj)
tmp <- read.xlsx(file, sheet = 1)
ohsu.metadata <- tmp[, c("patient_id", "gender", "lab_id", "specimen_type")]

tmp <- read.xlsx(file, sheet = 3)
tmp <- tmp[grepl(pattern="FAB", ignore.case=TRUE, tmp$lab_type),]
tmp <- unique(tmp[, c("patient_id", "lab_result")])
colnames(tmp) <- c("patient_id", "fab")
ohsu.metadata <- merge(ohsu.metadata, tmp, by = "patient_id", all = TRUE)

tmp <- read.xlsx(file, sheet = 5)
tmp <- tmp[as.numeric(tmp$age) > 0,]
tmp <- tmp[, c("patient_id", "lab_id", "age")]
tmp <- ddply(unique(tmp[, c("patient_id", "lab_id", "age")]), c("patient_id", "lab_id"), .fun = function(df) { df$age <- mean(df$age); return(unique(df)) })
ohsu.metadata <- merge(ohsu.metadata, tmp[, c("patient_id", "lab_id", "age")], by = c("patient_id", "lab_id"), all = TRUE)

tmp <- read.xlsx(file, sheet = 5)
tmp <- tmp[as.numeric(tmp$age) > 0,]
tmp <- tmp[, c("patient_id", "lab_id", "response")]
ohsu.simplified.response.tbl <- ddply(tmp[!is.na(tmp$response),], c("patient_id", "lab_id"), 
                                      .fun = function(df) {
                                        if(any(grepl(x=df$response, pattern="Complete Response"))) { return("Complete Response")}
                                        if(any(df$response == "Hematologic CR")) { return("Hematologic CR")}
                                        if(any(df$response == "Refractory")) { return("Refractory")}
                                        if(any(df$response == "Supportive/Palliative Care")) { return("Supportive/Palliative Care")}
                                        if(any(df$response == "Unknown")) { return("Unknown")}
                                        warning(paste0("What is ", df$response, "\n"))
                                      })
colnames(ohsu.simplified.response.tbl) <- c("patient_id", "lab_id", "response")
ohsu.metadata <- merge(ohsu.metadata, ohsu.simplified.response.tbl, by = c("patient_id", "lab_id"), all = TRUE)


## Remove NAs
ohsu.metadata <- ohsu.metadata[!is.na(ohsu.metadata$patient_id),]

## Map lab_id to seq_ids
synId <- "syn10083817"
obj <- synGet(synId, downloadFile = TRUE)
file <- getFileLocation(obj)
ohsu.rnaseq.sample.summary <- read.table(file, header=TRUE, sep="\t")
ohsu.metadata <- merge(ohsu.metadata, ohsu.rnaseq.sample.summary[, c("Original_LabID", "SeqID")], by.x = "lab_id", by.y = "Original_LabID", all = FALSE)

## One patient has 2 different FAB subtypes
ohsu.metadata <- ohsu.metadata[!(ohsu.metadata$SeqID == "20-00347"),]
rownames(ohsu.metadata) <- ohsu.metadata$SeqID
ohsu.metadata <- ohsu.metadata[, !(colnames(ohsu.metadata) %in% c("lab_id", "patient_id", "SeqID"))]

## Load FIMM metadata
synId <- "syn8270594"
obj <- synGet(id=synId, downloadFile = TRUE, downloadLocation = ".")
file <- getFileLocation(obj)
fimm.metadata <- read.table(file, header=TRUE, sep=",")
fimm.metadata <- fimm.metadata[, c("diseases.stage"), drop=F]

library(ComplexHeatmap)
library(circlize)

ohsu.metadata <- ohsu.metadata[, c("gender", "specimen_type", "fab", "response", "age")]

na.dist <- function(x,...) {
 t.dist <- dist(x,...)
 t.dist <- as.matrix(t.dist)
 t.limit <- 1.1*max(t.dist,na.rm=T)
 t.dist[is.na(t.dist)] <- t.limit
 t.dist <- as.dist(t.dist)
 return(t.dist)
}

plot.heatmap <- function(mat, annotations, row.annotations = NULL) {

   colnames(mat) <- gsub("X(.*)\\.(.*)", "\\1-\\2", colnames(mat))
   common <- intersect(rownames(annotations), colnames(mat))
   annotations <- annotations[common, , drop=F]
   mat <- mat[, common]

   col.list <- list()
   col.is.numeric <- list()
   for(col in colnames(annotations)) {
     flag <- is.na(annotations[, col]) | (annotations[, col] == "NA")
     annotations[flag, col] <- NA
   ##  annotations[flag,col] <- "NA"
     cols <- NULL
     vals <- annotations[, col]
     vals.no.na <- vals[!is.na(vals)]
   ##  vals.no.na <- vals[vals != "NA"]
     col.is.numeric[[col]] <- FALSE
     if(!is.numeric(vals.no.na)) {
       cols <- rainbow(length(unique(vals.no.na)))
       names(cols) <- sort(unique(vals.no.na))
       col.list[[col]] <- cols
     } else {
       col.is.numeric[[col]] <- TRUE
       vals <- as.numeric(vals)
       cols <- colorRamp2(c(min(0, floor(min(vals, na.rm=TRUE))), ceiling(max(vals, na.rm=TRUE))), c("white", "red"))
     }
   ##  col.list[[col]] <- cols
   }

   ## ha <- HeatmapAnnotation(df = annotations, col = col.list, na_col = "white", show_annotation_name = TRUE)

   colnames(mat) <- NULL
   ## ha <- HeatmapAnnotation(df = annotations[, !(colnames(annotations) %in% c("age"))], resp = anno_boxplot(mat), age = anno_points(annotations$age), col = col.list, na_col = "black", show_annotation_name = TRUE)
   ## ha <- HeatmapAnnotation(df = annotations[, !(colnames(annotations) %in% c("age"))], age = anno_points(annotations$age), col = col.list, na_col = "black", show_annotation_name = TRUE)
   ha <- NULL
   if("age" %in% colnames(annotations)) { 
     ha <- HeatmapAnnotation(df = annotations[, names(col.is.numeric)[col.is.numeric==FALSE], drop=F], age = anno_points(annotations$age), col = col.list, na_col = "black", show_annotation_name = TRUE)
   } else {
     ha <- HeatmapAnnotation(df = annotations[, names(col.is.numeric)[col.is.numeric==FALSE], drop=F], col = col.list, na_col = "black", show_annotation_name = TRUE)
   }
   ha_boxplot <- HeatmapAnnotation(response = anno_boxplot(mat, axis = TRUE, outline = FALSE), sensitive = anno_points(annotations$sensitive, axis=TRUE), insensitive = anno_points(annotations$insensitive, axis=TRUE), show_annotation_name = TRUE)

   if(!is.null(row.annotations)) {
     common <- intersect(rownames(row.annotations), rownames(mat))
     mat <- mat[common, ]
     row.annotations <- row.annotations[common, , drop=F]
   }
   rownames(mat) <- NULL
   hm <- Heatmap(mat, top_annotation = ha, clustering_distance_rows = function(x) na.dist(x), clustering_distance_columns = function(x) na.dist(x), na_col = "black", bottom_annotation = ha_boxplot, bottom_annotation_height = unit(3, "cm"))
   if(!is.null(row.annotations)) {
     row.col.list <- list()
     for(col in colnames(row.annotations)) {
       flag <- is.na(row.annotations[, col]) | (row.annotations[, col] == "NA")
       row.annotations[flag, col] <- NA
       cols <- NULL
       vals <- row.annotations[, col]
       vals.no.na <- vals[!is.na(vals)]
       cols <- rainbow(length(unique(vals.no.na)))
       names(cols) <- unique(vals.no.na)
       row.col.list[[col]] <- cols
     }

     ra <- rowAnnotation(df = row.annotations, col = row.col.list, show_annotation_name = TRUE)
     hm <- hm + ra
   }
   hm
}

plot.ordered.based.on.entropy <- function(mat) {
  tmp.scaled <- t(scale(t(mat)))
  ## heatmap.2(tmp.scaled, scale="none", trace="none", Rowv = TRUE, Colv = TRUE, dist = na.dist)
  ranked <- ldply(1:nrow(tmp.scaled), .fun = function(i) rank(tmp.scaled[i,])/max(rank(tmp.scaled[i,])))
  entropy <- unlist(apply(ranked, 2, function(col) -sum((col)*log(col))))
##  entropy <- unlist(apply(ranked, 2, function(col) max(col) - min(col)))
  nms <- names(sort(entropy))
  heatmap.2(tmp.scaled[, nms], scale="none", trace="none", Rowv = TRUE, Colv = FALSE, dist = na.dist)
}

