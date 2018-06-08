
subset.and.zscore.matrix <- function(df, row.frac.cutoff = 0.25, col.frac.cutoff = 0.25) {
##  frac.na <- 0.25
##  cutoff <- frac.na * nrow(df)
##  cutoff <- 20
  tmp <- df
  flag <- unlist(apply(tmp, 1, function(row) all(is.na(row))))
  tmp <- tmp[!flag,]
  flag <- unlist(apply(tmp, 1, function(row) length(which(!is.na(row))) < ( row.frac.cutoff * length(row) )))
  tmp <- tmp[!flag,]
  flag <- unlist(apply(tmp, 2, function(col) all(is.na(col))))
  tmp <- tmp[, !flag]
  flag <- unlist(apply(tmp, 2, function(col) length(which(!is.na(col))) < ( col.frac.cutoff * length(col) )))
  tmp <- tmp[, !flag]
  tmp.scaled <- t(scale(t(tmp)))
  tmp.scaled
}

safe.t.test <- function(...) {
  res <- tryCatch({t.test(...)}, error = function(e) { data.frame(p.value = 1) })
  res
}

safe.cor.test <- function(...) {
  ct <- tryCatch({cor.test(...)}, error = function(e) { data.frame(estimate = NA, p.value = NA) })
  ct
}

plot.drug.heatmap <- function(mat, drug.tbl, show.col.names = TRUE, show.row.names = TRUE) {
  if(!show.col.names) { colnames(mat) <- NULL }
  if(!show.row.names) { rownames(mat) <- NULL }
  hm <- Heatmap(mat, clustering_distance_rows = function(x) na.dist(x), clustering_distance_columns = function(x) na.dist(x), na_col = "black", cluster_rows = FALSE)
  row.annotation <- data.frame(cluster = as.character(drug.tbl$cluster))
  rownames(row.annotation) <- drug.tbl$OHSU_DRUG_NAME
  cols <- rainbow(length(unique(row.annotation$cluster)))
  names(cols) <- unique(row.annotation$cluster)
  row.col.list <- list("cluster" = cols)
  ra <- rowAnnotation(df = row.annotation, col = row.col.list, show_annotation_name = TRUE)
  hm <- hm + ra
  hm
}

plot.drug.heatmap.no.anno <- function(mat, show.col.names = TRUE, show.row.names = TRUE, ...) {
  if(!show.col.names) { colnames(mat) <- NULL }
  if(!show.row.names) { rownames(mat) <- NULL }
  hm <- Heatmap(mat, clustering_distance_rows = function(x) na.dist(x), clustering_distance_columns = function(x) na.dist(x), na_col = "black", cluster_rows = FALSE, ...)
  hm
}

plot.drug.heatmap.no.anno.no.cluster <- function(mat, show.col.names = TRUE, show.row.names = TRUE, ...) {
  if(!show.col.names) { colnames(mat) <- NULL }
  if(!show.row.names) { rownames(mat) <- NULL }
  hm <- Heatmap(mat, clustering_distance_rows = function(x) na.dist(x), clustering_distance_columns = function(x) na.dist(x), na_col = "black", cluster_cols = FALSE, cluster_rows = FALSE, ...)
  hm
}

## Plot heatmaps of the drugs
make.drug.heatmaps <- function(scaled.drcs) {
  heatmaps <-
    llply(data.sets,
          .fun = function(ds) {
                   title <- toupper(ds)
                   if(ds == "sanger") { title <- capwords(ds) }
                   hm <- plot.drug.heatmap.no.anno(scaled.drcs[[ds]], show.col.names = FALSE, column_title = title)
                 })
}

## pg <- plot_grid(l[[1]], l[[2]], labels = c("A", "B"), nrow=2, rel_widths = c(10,1.5), align="none")
## pg <- plot_grid(heatmaps[[1]], heatmaps[[2]], labels = c("A", "B"), nrow=2, align="none")

## assignInNamespace(x="xcmsRaw", value="my.xcmsRaw", ns=asNamespace("xcms"))

## Replace corrplot:::reorder_using_hclust

as.na.dist <- function(x,...) {
 t.dist <- as.dist(x,...)
 t.dist <- as.matrix(t.dist)
 t.limit <- 1.1*max(t.dist,na.rm=T)
 t.dist[is.na(t.dist)] <- t.limit
 t.dist <- as.dist(t.dist)
 return(t.dist)
}

my_reorder_using_hclust <- function (corr, hclust.method) 
{
    hc <- hclust(as.na.dist(1 - corr), method = hclust.method)
    order.dendrogram(as.dendrogram(hc))
}

assignInNamespace(x="reorder_using_hclust", value="my_reorder_using_hclust", ns=asNamespace("corrplot"))

reorder_using_hclust <- my_reorder_using_hclust


make.drug.correlation.plots_ <- function(orig.scaled.drcs, file.prefix = NULL, multiple.plots.in.one.figure = TRUE, common.drugs, drugs.to.label) {

  if(!is.null(file.prefix)) {
    png(paste0(file.prefix, ".png"))
  }
  if(multiple.plots.in.one.figure) {
    nrows <- ceiling(length(orig.scaled.drcs) / 2)
    par(mfrow=c(nrows,2))
  }
  indices <- 1:length(data.sets)
  names(indices) <- names(data.sets)
  corrplots <-
    llply(indices,
          .fun = function(indx) {
                   ds <- data.sets[[indx]]
                   title <- toupper(ds)
                   M <- cor(t(orig.scaled.drcs[[ds]][common.drugs, ]), use = "pairwise.complete.obs")
                   nz <- unlist(apply(M, 1, function(row) length(which(row < 0.1))))
                   flag <- !(rownames(M) %in% drugs.to.label)
                   rownames(M)[flag] <- ""
                   flag <- !(colnames(M) %in% drugs.to.label)
                   colnames(M)[flag] <- ""
                   ## order = "hclust" breaks when there are NAs, but we fixed with above function replacement
                   ## corrplot(M, type = "lower", tl.col = "black", tl.srt = 45, title = title)
                   ## mtext("A", line=2.5, at = 0, cex = 1.5)
                   mar <- c(5, 4, 4, 2) + 0.1
                   mar <- c(0, 1, 1, 1) + 0.1
                   mar <- c(0, 2, 1, 1) + 0.1
                   if(!is.null(file.prefix) && !multiple.plots.in.one.figure) {
                     png(paste0(file.prefix, "-", title, ".png"))
                   }
                   cp <- corrplot(M, type = "lower", order = "original", tl.col = "black", tl.srt = 45, na.label = "square", na.label.col = "yellow", mar = mar, tl.pos = "ld")
                   title(xlab = "Drug", line = -3, cex.lab = 1.5)
                   title(ylab = "Drug", line = 2, cex.lab = 1.5)
                   ## mtext("A", line=2.5, at = 0, cex = 1.5)
                   mtext(LETTERS[indx], line=-1, cex = 1.5, side = 3, at = -1)
                   mtext(title, line=-1, cex = 1.5, side = 3, at = ncol(M)/2)
                   v <- par("usr")
                   if(any(is.na(M))) { legend(v[1] + 0.6 * (v[2] - v[1]), v[3] + 0.75 * (v[4] - v[3]), legend="NA", fill="yellow", bty = "n") }
                   if(!is.null(file.prefix) && !multiple.plots.in.one.figure) {
                     d <- dev.off()
                   }
                 })
  if(!is.null(file.prefix) && multiple.plots.in.one.figure) {
    d <- dev.off()
  }
}

make.drug.correlation.plots <- function(orig.scaled.drcs, file.prefix = NULL, multiple.plots.in.one.figure = TRUE) {
  common.drugs <- Reduce("intersect", lapply(orig.scaled.drcs, function(mat) rownames(mat)))

  ## Determine the correlation order, which we will impose below
  ds <- names(data.sets)[1]
  if("fimm" %in% data.sets) {
    ds <- "fimm"
  }
  M <- cor(t(orig.scaled.drcs[[ds]][common.drugs, ]), use = "pairwise.complete.obs")
  nz <- unlist(apply(M, 1, function(row) length(which(row < 0.1))))
  ## order = "hclust" breaks when there are NAs, but we fixed with above function replacement
  ## corrplot(M, type = "lower", tl.col = "black", tl.srt = 45, title = title)
  ## mtext("A", line=2.5, at = 0, cex = 1.5)
  ord <- corrplot(M, type = "lower", order = "hclust")
  common.drugs <- rownames(ord)

  drugs.to.label <- "Venetoclax"
  drugs.to.label <- rownames(orig.scaled.drcs[[1]])
  drugs.to.label <- c()

  make.drug.correlation.plots_(orig.scaled.drcs, file.prefix = file.prefix, multiple.plots.in.one.figure = multiple.plots.in.one.figure, 
                               common.drugs = common.drugs, drugs.to.label = drugs.to.label)
}

## Plot correlation of each drug with GLDS
plot.drug.correlation.with.glds <- function(mat) {
  glds <- colMeans(mat, na.rm=TRUE)
  indices <- 1:nrow(mat)
  names(indices) <- row.names(mat)
  drug.glds.corr <- ldply(indices,
                          .fun = function(i) {
                                   drug.resp <- mat[i, ]
                                   ct <- safe.cor.test(drug.resp, glds)
                                   data.frame(corr = ct$estimate, pval = ct$p.value)
                          })
  colnames(drug.glds.corr) <- c("drug", "corr", "pval")
  drug.glds.corr <- na.omit(drug.glds.corr)
  drug.glds.corr <- drug.glds.corr[order(drug.glds.corr$corr, decreasing=FALSE), ]

  g <- ggplot()
  g <- g + geom_point(data = drug.glds.corr, aes(x = corr, y = -log10(pval)))
  df <- subset(drug.glds.corr, pval > 0.1)
  df <- subset(drug.glds.corr, corr < 0)
  df <- drug.glds.corr[1:min(5, nrow(drug.glds.corr)),]
  df$drug <- unlist(lapply(df$drug, function(str) gsub(str, pattern="([^(]+)[ ]*\\(.*", replacement="\\1")))
  g <- g + xlab("Pearson Correlation")
  g <- g + xlim(c(min(min(drug.glds.corr$corr, na.rm=TRUE) - 0.1, -0.2), max(drug.glds.corr$corr, na.rm=TRUE) + 0.1))
  g <- g + geom_text(data = df, aes(x = corr, y = -log10(pval), label = drug), position = position_jitter(width=0, height=1), hjust = 1)
  list("g" = g, "drug.glds.corr" = drug.glds.corr)
}

if(FALSE) {

## compares every smiles on list one across list 2
cat("Calculating all pairwise distances between drugs\n")
drug.struct.dist <- ldply(fp.common, .fun = function(i) {
  sim <- sapply(fp.common, function(j) {
    1 - distance(i, j, method = "tanimoto")
  })
  sim
})

rownames(drug.struct.dist) <- drug.struct.dist$.id
drug.struct.dist <- drug.struct.dist[, !(colnames(drug.struct.dist) == ".id")]


orig.ohsu.drc <- prepare.drug.response.matrix(orig.fits[["ohsu"]], drug.col = "inhibitor", patient.col = "patient_id", response.col = "auc.ll4")
orig.ohsu.drc.scaled <- subset.and.zscore.matrix(orig.ohsu.drc, col.frac.cutoff = 0.25, row.frac.cutoff = 0.25)
orig.ohsu.drug.response.dist <- as.matrix(na.dist(orig.ohsu.drc.scaled))

## tmp <- merge(orig.fits[["fimm"]], drug.map[, c("FIMM_DRUG_NAME", "OHSU_DRUG_NAME")], by.x = "DRUG_NAME", by.y = "FIMM_DRUG_NAME")
orig.fimm.drc <- prepare.drug.response.matrix(orig.fits[["fimm"]], drug.col = "OHSU_DRUG_NAME", patient.col = "PATIENT_ID", response.col = "auc.ll4")
orig.fimm.drc.scaled <- subset.and.zscore.matrix(orig.fimm.drc, col.frac.cutoff = 0.25, row.frac.cutoff = 0.25)
orig.fimm.drug.response.dist <- as.matrix(na.dist(orig.fimm.drc.scaled))

ohsu.drc <- prepare.drug.response.matrix(fits[["ohsu"]], drug.col = "inhibitor", patient.col = "patient_id", response.col = "auc.ll4")
ohsu.drc.scaled <- subset.and.zscore.matrix(ohsu.drc, col.frac.cutoff = 0.25, row.frac.cutoff = 0.25)
ohsu.drug.response.dist <- as.matrix(na.dist(ohsu.drc.scaled))

## tmp <- merge(fits[["fimm"]], drug.map[, c("FIMM_DRUG_NAME", "OHSU_DRUG_NAME")], by.x = "DRUG_NAME", by.y = "FIMM_DRUG_NAME")
fimm.drc <- prepare.drug.response.matrix(fits[["fimm"]], drug.col = "OHSU_DRUG_NAME", patient.col = "PATIENT_ID", response.col = "auc.ll4")
fimm.drc.scaled <- subset.and.zscore.matrix(fimm.drc, col.frac.cutoff = 0.25, row.frac.cutoff = 0.25)
fimm.drug.response.dist <- as.matrix(na.dist(fimm.drc.scaled))

orig.common.drugs <- Reduce("intersect", lapply(list(orig.ohsu.drc.scaled, orig.fimm.drc.scaled), function(mat) rownames(mat)))
common.drugs <- Reduce("intersect", lapply(list(ohsu.drc.scaled, fimm.drc.scaled), function(mat) rownames(mat)))

pdf("ohsu-all-drug-glds-corr.pdf")
ret <- plot.drug.correlation.with.glds(orig.ohsu.drc.scaled)
g <- ret$g
g <- g + ggtitle("OHSU: (All) Drug correlation with GLDS")
print(g)
d <- dev.off()

pdf("ohsu-common-drug-glds-corr.pdf")
ret <- plot.drug.correlation.with.glds(orig.ohsu.drc.scaled[orig.common.drugs, ])
g <- ret$g
common.ohsu.glds.corr <- ret$drug.glds.corr
g <- g + ggtitle("OHSU: (Common) Drug correlation with GLDS")
print(g)
d <- dev.off()

pdf("fimm-all-drug-glds-corr.pdf")
ret <- plot.drug.correlation.with.glds(orig.fimm.drc.scaled)
g <- ret$g
g <- g + ggtitle("FIMM: (All) Drug correlation with GLDS")
print(g)
d <- dev.off()

pdf("fimm-common-drug-glds-corr.pdf")
ret <- plot.drug.correlation.with.glds(orig.fimm.drc.scaled[orig.common.drugs, ])
g <- ret$g
common.fimm.glds.corr <- ret$drug.glds.corr
g <- g + ggtitle("FIMM: (Common) Drug correlation with GLDS")
print(g)
d <- dev.off()

common.glds.corr <- merge(common.ohsu.glds.corr, common.fimm.glds.corr, by = "drug", suffixes = c(".ohsu", ".fimm"))
g <- ggplot(data = common.glds.corr)
g <- g + geom_point(aes(x = corr.ohsu, y = corr.fimm))
g <- g + xlab("OHSU Pearson") + ylab("FIMM Pearson")
df <- subset(common.glds.corr, (corr.ohsu < 0) | (corr.fimm < 0))
df$drug <- unlist(lapply(df$drug, function(str) gsub(str, pattern="([^(]+)[ ]*\\(.*", replacement="\\1")))
df2 <- data.frame(x = common.glds.corr$corr.ohsu, y = common.glds.corr$corr.fimm)
g <- g + geom_text(data = df, aes(x = corr.ohsu, y = corr.fimm, label = drug), position = position_jitter(width=0, height=0), hjust = 0)
g <- g + geom_smooth(data = df2, aes(x = x, y = y), method='lm')
g <- g + geom_text(x = 0, y = 0.75, label = lm_eqn(df2), parse=TRUE)

pdf("fimm-ohsu-correlations-with-glds.pdf")
print(g)
d <- dev.off()

common <- intersect(rownames(drug.struct.dist), rownames(orig.ohsu.drug.response.dist))
tmp1 <- orig.ohsu.drug.response.dist[common, common]
tmp2 <- drug.struct.dist[common, common]

x <- tmp1[upper.tri(tmp1, diag=FALSE)]
y <- tmp2[upper.tri(tmp2, diag=FALSE)]
ct <- cor.test(x, y)

common <- intersect(rownames(orig.fimm.drug.response.dist), rownames(orig.ohsu.drug.response.dist))
tmp1 <- orig.ohsu.drug.response.dist[common, common]
tmp2 <- orig.fimm.drug.response.dist[common, common]

x <- tmp1[upper.tri(tmp1, diag=FALSE)]
y <- tmp2[upper.tri(tmp2, diag=FALSE)]
ct <- cor.test(x, y)

## pdf(paste0(file.prefix, "-ohsu-all-heatmap.pdf"))
## plot.heatmap(tmp.scaled, ohsu.metadata, drug.metadata)
## d <- dev.off()

mat <- orig.ohsu.drc.scaled
hm <- Heatmap(mat, clustering_distance_rows = function(x) na.dist(x), clustering_distance_columns = function(x) na.dist(x), na_col = "black")
ohsu.row.order <- rownames(mat)[row_order(hm)[[1]]]
ohsu.col.order <- colnames(mat)[column_order(hm)]

mat <- orig.fimm.drc.scaled
hm <- Heatmap(mat, clustering_distance_rows = function(x) na.dist(x), clustering_distance_columns = function(x) na.dist(x), na_col = "black")
fimm.row.order <- rownames(mat)[row_order(hm)[[1]]]
fimm.col.order <- colnames(mat)[column_order(hm)]

mats <- list("orig-ohsu-common-drugs" = orig.ohsu.drc.scaled[orig.common.drugs,], 
             "ohsu-common-drugs" = ohsu.drc.scaled[common.drugs,],
             "orig-ohsu-all-drugs" = orig.ohsu.drc.scaled,
             "orig-fimm-common-drugs" = orig.fimm.drc.scaled[orig.common.drugs,], 
             "fimm-common-drugs" = fimm.drc.scaled[common.drugs,],
             "orig-fimm-all-drugs" = orig.fimm.drc.scaled,
             "orig-fimm-ohsu-order-common-drugs" = orig.fimm.drc.scaled[orig.common.drugs,], 
             "fimm-ohsu-order-common-drugs" = fimm.drc.scaled[common.drugs,],
             "orig-fimm-ohsu-order-all-drugs" = orig.fimm.drc.scaled)

mat.rows <- list("orig-ohsu-common-drugs" = ohsu.row.order,
             "ohsu-common-drugs" = ohsu.row.order,
             "orig-ohsu-all-drugs" = ohsu.row.order,
             "orig-fimm-common-drugs" = fimm.row.order,
             "fimm-common-drugs" = fimm.row.order,
             "orig-fimm-all-drugs" = fimm.row.order,
             "orig-fimm-ohsu-order-common-drugs" = ohsu.row.order, 
             "fimm-ohsu-order-common-drugs" = ohsu.row.order,
             "orig-fimm-ohsu-order-all-drugs" = ohsu.row.order)


mat.cols <- list("orig-ohsu-common-drugs" = ohsu.col.order,
             "ohsu-common-drugs" = ohsu.col.order,
             "orig-ohsu-all-drugs" = ohsu.col.order,
             "orig-fimm-common-drugs" = fimm.col.order,
             "fimm-common-drugs" = fimm.col.order,
             "orig-fimm-all-drugs" = fimm.col.order,
             "orig-fimm-ohsu-order-common-drugs" = fimm.col.order, 
             "fimm-ohsu-order-common-drugs" = fimm.col.order,
             "orig-fimm-ohsu-order-all-drugs" = fimm.col.order)


for(nm in names(mats)) {
  pdf(paste0(nm, "-heatmap.pdf"))
  mat <- mats[[nm]]
##  hm <- Heatmap(mat, clustering_distance_rows = function(x) na.dist(x), clustering_distance_columns = function(x) na.dist(x), na_col = "black")
  rows <- mat.rows[[nm]][mat.rows[[nm]] %in% rownames(mat)]
  cols <- mat.cols[[nm]][mat.cols[[nm]] %in% colnames(mat)]
  mat <- mat[rows, cols]
  colnames(mat) <- NULL
  hm <- Heatmap(mat, cluster_rows = FALSE, cluster_columns = FALSE, na_col = "black")
  print(hm)
  d <- dev.off()
}

for(nm in c("ohsu-common-drugs", "fimm-common-drugs")) {
  pdf(paste0(nm, "-hc-heatmap.pdf"))
  mat <- mats[[nm]]
  colnames(mat) <- NULL
  hm <- Heatmap(mat, clustering_distance_rows = function(x) na.dist(x), clustering_distance_columns = function(x) na.dist(x), na_col = "black")
  print(hm)
  d <- dev.off()
}

plot.drug.clusters <- function(mat, k.char, drug.map, hc.file = NULL) {
  ## ct <- cutreeDynamic(hc, distM = as.matrix(na.dist(mat)), minClusterSize = 1)
  ## df <- data.frame(drug = labels(hc), cluster = ct)
  ## df <- data.frame(drug = labels(hc), cluster = ct)
  ## Plot the tree to see manually where we should make the cut
  ## plot(hc)
  if(!is.null(hc.file)) { pdf(hc.file) }
  dend <- mat %>% na.dist %>% hclust %>% as.dendrogram
  # plot + color the dend's branches before, based on 3 clusters:
  ## dend %>% color_branches(k=3) %>% plot(horiz=FALSE)
  dend %>% plot(horiz=FALSE)
  ## add horiz rect
  #dend %>% rect.dendrogram(k=3,horiz=TRUE)
  ## add horiz (well, vertical) line:
  #abline(v = heights_per_k.dendrogram(dend)["3"] + .6, lwd = 2, lty = 2, col = "blue")
  h <- unname(heights_per_k.dendrogram(dend)[k.char])
  abline(h = h, lwd = 2, lty = 2, col = "blue")
  if(!is.null(hc.file)) { d <- dev.off() }
  hc <- hclust(na.dist(mat))
  ct <- cutree(hc, k=as.numeric(k.char))
  df <- data.frame(drug = names(ct), cluster = unname(ct))
  drug.clusters <- dlply(df, .variables = "cluster", .fun = function(tbl) as.character(tbl$drug))
  drug.tbl <- merge(drug.map, df, by.x = "OHSU_DRUG_NAME", by.y = "drug")
  rownames(mat) <- unlist(lapply(rownames(mat), function(str) gsub(str, pattern="([^(]+)[ ]+\\(.*", replacement="\\1")))
  drug.tbl$OHSU_DRUG_NAME <- unlist(lapply(drug.tbl$OHSU_DRUG_NAME, function(str) gsub(str, pattern="([^(]+)[ ]+\\(.*", replacement="\\1")))
  rownames(drug.tbl) <- as.character(drug.tbl$OHSU_DRUG_NAME)
  common <- intersect(rownames(mat), rownames(drug.tbl))
  drug.tbl <- drug.tbl[common,]
  drug.tbl <- drug.tbl[order(drug.tbl$cluster),]                                       
  common <- rownames(drug.tbl)
  mat <- mat[common, ]
  hm <- plot.drug.heatmap(mat, drug.tbl, show.col.names = FALSE)
  list("heatmap" = hm, "clusters" = drug.clusters)
}


mat <- orig.ohsu.drc.scaled[orig.common.drugs,]
hc <- hclust(na.dist(mat))
ct <- cutreeDynamic(hc, distM = as.matrix(na.dist(mat)), minClusterSize = 1)
df <- data.frame(drug = labels(hc), cluster = ct)
orig.ohsu.drug.clusters <- dlply(df, .variables = "cluster", .fun = function(tbl) as.character(tbl$drug))

mat <- orig.fimm.drc.scaled[orig.common.drugs,]
hc <- hclust(na.dist(mat))
## ct <- cutreeDynamic(hc, distM = as.matrix(na.dist(mat)), minClusterSize = 1)
## ct <- cutreeHybrid(hc, distM = as.matrix(na.dist(mat)), minClusterSize = 1)
## df <- data.frame(drug = labels(hc), cluster = ct)
## Plot the tree to see manually where we should make the cut
## plot(hc)

mat <- ohsu.drc.scaled[common.drugs,]
ret <- plot.drug.clusters(mat, "12", drug.map, hc.file = "ohsu-common-drugs-hc.pdf")
ohsu.drug.clusters <- ret$clusters
pdf("ohsu-common-drugs-clustered-heatmap.pdf")
print(ret$heatmap)
d <- dev.off()

mat <- fimm.drc.scaled[common.drugs,]
ret <- plot.drug.clusters(mat, "12", drug.map, hc.file = "fimm-common-drugs.hc.pdf")
fimm.drug.clusters <- ret$clusters
pdf("fimm-common-drugs-clustered-heatmap.pdf")
print(ret$heatmap)
d <- dev.off()

jaccard.index <- function(set1, set2) {
  if(length(intersect(set1, set2)) == 0) { return(0) }
  length(intersect(set1, set2))/length(union(set1, set2))
}

## Perform meta-clustering using the Jaccard similarity as a distance metric/edge weight
meta.mcl_ <- function(cls, ...) {
  edges <- ldply(cls,
                .fun = function(cluster1) {
                         df.c1 <- ldply(cls,
                                        .fun = function(cluster2) { 
                                                 sim <- jaccard.index(cluster1, cluster2)
                                                 data.frame(jaccard = sim)
                                        })
                         colnames(df.c1) <- c("node2", "jaccard")   
                         df.c1 
                })
  colnames(edges) <- c("node1", "node2", "jaccard")
  ## Remove self-loops
##  flag <- (edges$method1 == edges$method2) & (edges$cluster1 == edges$cluster2)
##  edges <- edges[!flag,]

  adj.matrix <- as.matrix(get.adjacency(graph.data.frame(edges[, c("node1", "node2", "jaccard")]), attr="jaccard"))
  mcl.out <- mcl(adj.matrix, ...)
  df <- data.frame(meta.cluster = 0, name = colnames(adj.matrix))
  if("Cluster" %in% names(mcl.out)) {
    df$meta.cluster <- mcl.out$Cluster
  }
  df
}

meta.mcl <- function(clusters, ...) {
  cls <- unlist(clusters, recursive = FALSE)
  meta.mcl_(cls, ...)
}

## df <- meta.mcl(clusters, inflation = 1, addLoops = FALSE)

source("calc-mcl.R")

## Change format from item to cluster assigments to a list of all items in a cluster
get.cluster.items <- function(res, cluster.col, item.col) {
  cluster.members <- dlply(res[res[, cluster.col] != 0,,drop=FALSE], .variables = cluster.col, 
                           .fun = function(df) paste(df[, item.col], collapse = ","))
  cluster.members
}

## Expand clusters to items within that cluster
expand.cluster.to.items <- function(meta.cluster, cluster.to.item.map) {
  clusters <- unlist(strsplit(meta.cluster, split=",[ ]*"))
  items <- unlist(llply(clusters, .fun = function(cls) cluster.to.item.map[[cls]]))
  paste(sort(items), collapse=",")
}

## Finally, find drugs that occur twice in a cluster--i.e., are included from both data sets
find.consistently.coclustered.drugs <- function(cluster.of.drugs) {
  clusters.of.consistent.drugs <- llply(clusters.of.drugs, .parallel = FALSE,
                                         .fun = function(cluster) {
                                                  items <- unlist(strsplit(cluster, split=",[ ]*"))
                                                  if(!(any(duplicated(items)))) { return(NA) }
                                                  items <- items[duplicated(items)]
                                                  items
                                         })
  clusters.of.consistent.drugs <- clusters.of.consistent.drugs[!is.na(clusters.of.consistent.drugs)]

  ## Limit to cases in which there are at least 2 consistent drugs in a cluster
  clusters.of.consistent.drugs <- llply(clusters.of.consistent.drugs,
                                         .fun = function(cluster) {
                                                  items <- cluster
                                                  if(length(items) == 1) { return(NA) }
                                                  items
                                         })                                   
  clusters.of.consistent.drugs <- clusters.of.consistent.drugs[!is.na(clusters.of.consistent.drugs)]
  clusters.of.consistent.drugs

}    


inflations <- seq(from=1, to=10, by=0.1)
names(inflations) <- inflations

## clusters <- list("ohsu" = orig.ohsu.drug.clusters, "fimm" = orig.fimm.drug.clusters)
clusters <- list("ohsu" = ohsu.drug.clusters, "fimm" = fimm.drug.clusters)

sils <- ldply(inflations, 
              .fun = function(inflation) {
                       data.frame(sil = calc.mcl.weighted.silhouette(clusters, n.iters = 1000, frac.to.downsample = 0.8, addLoops = FALSE, inflation = inflation, allow1 = TRUE))
                     })
colnames(sils) <- c("inflation", "sil")

save.image(".Rdata.inflate")
stop("stop")

inflations2 <- seq(from=1.4, to=2.0, by=0.1)
inflations2 <- seq(from=1.4, to=1.8, by=0.4)
names(inflations2) <- inflations2

sils2 <- ldply(inflations2, 
              .fun = function(inflation) {
                       data.frame(sil = calc.mcl.weighted.silhouette(clusters, n.iters = 1000, frac.to.downsample = 0.8, addLoops = FALSE, inflation = inflation, allow1 = TRUE))
                     })
colnames(sils2) <- c("inflation", "sil")


cls <- unlist(clusters, recursive = FALSE)
cluster.to.item.map <- cls
mcl.all.out <- meta.mcl_(cls, addLoops = FALSE, allow1 = TRUE, inflation = 1.5)

clusters.of.clusters <- get.cluster.items(mcl.all.out, "meta.cluster", "name")
clusters.of.drugs <- unlist(llply(unlist(clusters.of.clusters), .fun = function(meta.cluster) expand.cluster.to.items(meta.cluster, cluster.to.item.map)))
clusters.of.consistent.drugs <- find.consistently.coclustered.drugs(cluster.of.drugs) 



## Calculate structural distances for intra- and inter-cluster drugs
df <- ldply(clusters.of.consistent.drugs, .parallel = FALSE,
            .fun = function(cluster1) {
                     cls1.res <- ldply(clusters.of.consistent.drugs, .parallel = FALSE,
                                       .fun = function(cluster2) {
                                                cls2.res <- ldply(cluster1, .parallel = FALSE,
                                                                  .fun = function(drug1) {
                                                                           drg1.res <- ldply(cluster2, .parallel = FALSE,
                                                                                             .fun = function(drug2) {
                                                                                                      data.frame(drug2 = drug2, distance =  drug.struct.dist[drug1, drug2])
                                                                                                    })
                                                                           colnames(drg1.res) <- c("drug2", "distance")
                                                                           drg1.res <- cbind(drug1 = drug1, drg1.res)
                                                                           drg1.res
                                                                         })
                                                cls2.res
                                              })
                     colnames(cls1.res) <- c("cluster2", "drug1", "drug2", "distance")
                     cls1.res
                   })
colnames(df)[1] <- "cluster1"
df <- df[df$drug1 != df$drug2,]
df$relation <- unlist(apply(df[, c("cluster1", "cluster2")], 1, function(row) ifelse(row[1] == row[2], "Intra-Cluster", "Inter-Cluster")))

g <- ggplot(data = df, aes(x = relation, y = distance))
g <- g + geom_violin()
## g <- g + geom_beeswarm()
g <- g + geom_boxplot(width = 0.5)
g <- g + ylab("Tanimoto Distance")
pdf("cluster-relation-vs-tanimoto.pdf")
print(g)
d <- dev.off()

save(clusters.of.consistent.drugs, file="clusters.of.consistent.drugs.Rd")

## Expand back out to format in which a drug is assigned to a cluster
consistent.drug.tbl <- ldply(clusters.of.consistent.drugs,
                                       .fun = function(cluster) {
                                                items <- cluster
                                                data.frame(drug = items)
                                       })
colnames(consistent.drug.tbl) <- c("cluster", "drug")

consistent.drug.tbl <- merge(drug.map, consistent.drug.tbl, by.x = "OHSU_DRUG_NAME", by.y = "drug")
consistent.drug.tbl <- consistent.drug.tbl[order(consistent.drug.tbl$cluster),]                                       
             
consistent.drug.tbl[, c("OHSU_DRUG_NAME", "Mechanism.Targets.1", "cluster")]

mat <- fimm.drc.scaled[as.character(consistent.drug.tbl$OHSU_DRUG_NAME),]   
hm <- plot.drug.heatmap(mat, consistent.drug.tbl, show.col.names = FALSE)
pdf("fimm-consistently-co-clustered-drugs-heatmap.pdf")
print(hm)
d <- dev.off()

mat <- ohsu.drc.scaled[as.character(consistent.drug.tbl$OHSU_DRUG_NAME),]   
hm <- plot.drug.heatmap(mat, consistent.drug.tbl, show.col.names = FALSE)
pdf("ohsu-consistently-co-clustered-drugs-heatmap.pdf")
print(hm)
d <- dev.off()



## Look at intra- vs inter-correlations for these drugs (in both data sets)
## Modify code to do drug classes
## Can we model each class in OHSU and apply to FIMM -- plot resp vs predicted, coloring by drug
## Plot correlation vs mean response


## TODO
## - cluster by drug response
## - look at other annotations
## - predict all within cluster, with covariate for each drug
## - look at overlap with drug covariate model
## - look at prediction accuracy with drug covariate model (limit post hoc to drugs in common)

## hm <- Heatmap(mat, clustering_distance_rows = function(x) na.dist(x), clustering_distance_columns = function(x) na.dist(x), na_col = "black")

pdf(paste0(file.prefix, "-ohsu-all-heatmap.pdf"))
plot.heatmap(tmp.scaled, ohsu.metadata, drug.metadata)
d <- dev.off()

tmp <- fimm.drc
flag <- unlist(apply(tmp, 1, function(row) all(is.na(row))))
tmp <- tmp[!flag,]
flag <- unlist(apply(tmp, 1, function(row) length(which(is.na(row))) > frac.na * length(row)))
tmp <- tmp[!flag,]
flag <- unlist(apply(tmp, 2, function(col) all(is.na(col))))
tmp <- tmp[, !flag]
flag <- unlist(apply(tmp, 2, function(col) length(which(is.na(col))) > frac.na * length(col)))
tmp <- tmp[, !flag]
tmp.scaled <- t(scale(t(tmp)))

##cutoff <- 10^-5
##sensitive <- unlist(apply(tmp.scaled, 2, function(col) safe.t.test(col, alternative="greater")$p.value < cutoff))
##insensitive <- unlist(apply(tmp.scaled, 2, function(col) safe.t.test(col, alternative="less")$p.value < cutoff))
sensitive <- unlist(apply(tmp.scaled, 2, function(col) log(safe.t.test(col, alternative="greater")$p.value)))
insensitive <- unlist(apply(tmp.scaled, 2, function(col) log(safe.t.test(col, alternative="less")$p.value)))
df <- data.frame(sensitive = sensitive, insensitive = insensitive)
rownames(df) <- colnames(tmp.scaled)
fimm.metadata <- merge(orig.fimm.metadata, df, by = "row.names", all = TRUE)
rownames(fimm.metadata) <- fimm.metadata$Row.names
fimm.metadata <- fimm.metadata[, !(colnames(fimm.metadata) %in% "Row.names")]

hc <- hclust(na.dist(t(tmp.scaled)))
k <- 3
cat(paste0("\nRidge quant cut k = ", k, "\n"))
## print(cutree(hc, k = k))
clusters <- t(t(cutree(hc, k=k)))
colnames(clusters) <- "cluster"
rownames(clusters) <- unlist(lapply(rownames(clusters), function(nm) gsub("X(.*)\\.(.*)", "\\1-\\2", nm)))
fimm.metadata <- merge(clusters, fimm.metadata, by = "row.names", all = TRUE)
rownames(fimm.metadata) <- fimm.metadata$Row.names
fimm.metadata <- fimm.metadata[, !(colnames(fimm.metadata) %in% "Row.names")]

tmp <- ddply(fimm.metadata, "cluster", .fun = function(df) mean(df$sensitive))
sensitive.cluster <- tmp$cluster[!is.na(tmp$V1) & (tmp$V1 == min(tmp$V1, na.rm=TRUE))]
tmp <- ddply(fimm.metadata, "cluster", .fun = function(df) mean(df$insensitive))
insensitive.cluster <- tmp$cluster[!is.na(tmp$V1) & (tmp$V1 == min(tmp$V1, na.rm=TRUE))]
fimm.metadata$status <- NA
fimm.metadata$status[!is.na(fimm.metadata$cluster) & (fimm.metadata$cluster == sensitive.cluster)] <- "sensitive"
fimm.metadata$status[!is.na(fimm.metadata$cluster) & (fimm.metadata$cluster == insensitive.cluster)] <- "insensitive"

fimm.sensitive.samples <- data.frame(sample = rownames(fimm.metadata)[!is.na(fimm.metadata$status) & (fimm.metadata$status == "sensitive")])
fimm.not.sensitive.samples <- data.frame(sample = rownames(fimm.metadata)[is.na(fimm.metadata$status) | (fimm.metadata$status != "sensitive")])
fimm.intermediate.sensitive.samples <- data.frame(sample = rownames(fimm.metadata)[is.na(fimm.metadata$status)])
fimm.insensitive.samples <- data.frame(sample = rownames(fimm.metadata)[!is.na(fimm.metadata$status) & (fimm.metadata$status == "insensitive")])
write.table(file="fimm.sensitive.samples.tsv", fimm.sensitive.samples, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
write.table(file="fimm.not.sensitive.samples.tsv", fimm.not.sensitive.samples, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
write.table(file="fimm.intermediate.sensitive.samples.tsv", fimm.intermediate.sensitive.samples, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
write.table(file="fimm.insensitive.samples.tsv", fimm.insensitive.samples, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)

pdf(paste0(file.prefix, "-fimm-all-heatmap.pdf"))
plot.heatmap(tmp.scaled, fimm.metadata, drug.metadata)
d <- dev.off()

save.image(rdata.file)

cat("Exiting\n")
q(status = 0)

} # if(FALSE)
