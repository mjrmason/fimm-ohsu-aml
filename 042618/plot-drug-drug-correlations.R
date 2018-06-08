drugs.of.interest <- c("XAV-939", "Trametinib", "Venetoclax", "Selumetinib", "Selinexor", "SNS-032", "Panobinostat", "Palbociclib", "PD184352", "Lovastatin", "BI 2536", "AGI-5198",
                       "Navitoclax")


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


