calc.mcl.weighted.silhouette <- function(clusters, n.iters = 1000, frac.to.downsample = 0.8, ...) {
  cls <- unlist(clusters, recursive = FALSE)
  mcl.all.out <- meta.mcl_(cls, ...)

  ## Create a dissimilarity matrix across all clusters (using Jaccard coefficient as similarity)
  all.clusters <- names(cls)
  diss <- as.data.frame(matrix(data = 0, nrow = length(all.clusters), ncol = length(all.clusters)))
  rownames(diss) <- all.clusters
  colnames(diss) <- all.clusters
  for(i in 1:length(cls)) {
    for(j in 1:length(cls)) {
      sim <- jaccard.index(cls[[i]], cls[[j]])
      diss[names(cls)[i], names(cls)[j]] <- 1 - sim
      diss[names(cls)[j], names(cls)[i]] <- 1 - sim
    }
  }

  if(length(unique(mcl.all.out$meta.cluster)) == 1) { return(-Inf) }
##  if(all(mcl.all.out$meta.cluster == 0)) { return(-Inf) }

  all.entries <- llply(clusters, .fun = function(lst) unname(unlist(lst)))
  common <- Reduce(intersect, all.entries)
  for(i in 1:length(all.entries)) {
    if(any(!(union(all.entries[[i]], common) %in% intersect(all.entries[[i]], common)))) {
      cat(paste0("clusters[[", i, "]] != common: ", union(setdiff(all.entries[[i]], common), setdiff(common, all.entries[[i]])), "\n"))
    }
  }

  iters <- 1:n.iters
  mcl.downsampled.clusters <- 
    llply(iters, .parallel = TRUE,
          .fun = function(i) {
                   set.seed(i)
                   downsampled.items <- common[sample.int(n = length(common), size = floor(frac.to.downsample * length(common)))]
                   downsampled.clusters <- llply(clusters,
                                                 .fun = function(method) {
                                                          llply(method,
                                                                .fun = function(cluster) {
                                                                         ds <- cluster[cluster %in% downsampled.items]
                                                                         ds
                                                                 })
                                                 })
                   ds.cls <- unlist(downsampled.clusters, recursive = FALSE)
                   meta.mcl_(ds.cls, ...)
          })

  ## Create a list, where each entry is a list of items in one of the (meta-)clusters
  downsampled.cluster.members <- unlist(llply(mcl.downsampled.clusters,
                                  .fun = function(mcl.res) {
                                           mcl.tmp <- mcl.res[mcl.res$meta.cluster != 0,,drop=FALSE]
                                           mcl.tmp <- mcl.res
                                           dlply(mcl.tmp, .variables = "meta.cluster", 
                                                 .fun = function(df) paste(df$name, collapse = ","))
                                         }))

  ## Create a consensus matrix, where entry i,j is the fraction of the iterations in which
  ## clusters i and j are clustered together
  all.clusters <- names(cls)
  cons.matrix <- as.data.frame(matrix(data = 0, nrow = length(all.clusters), ncol = length(all.clusters)))
  rownames(cons.matrix) <- all.clusters
  colnames(cons.matrix) <- all.clusters
  for(i in 1:length(downsampled.cluster.members)) {
    items <- unlist(strsplit(downsampled.cluster.members[[i]], split=",[ ]*"))
    if(length(items) < 2) { next }
    for(j in 1:(length(items)-1)) {
      itemj <- items[j]
      for(k in (j+1):length(items)) {
        itemk <- items[k]
        cons.matrix[itemj, itemk] <- cons.matrix[itemj, itemk] + 1
        cons.matrix[itemk, itemj] <- cons.matrix[itemk, itemj] + 1
      }
    }
  }
  cons.matrix <- cons.matrix / n.iters

  ## Similarly, create a co-(meta-)clustering matrix if two clusters were clustered together by
  ## mcl on _all_ of the data
  mcl.tmp <- mcl.all.out[mcl.all.out$meta.cluster != 0,,drop=FALSE]
  mcl.tmp <- mcl.all.out
  cluster.members <- dlply(mcl.tmp, .variables = "meta.cluster", 
                                                 .fun = function(df) paste(df$name, collapse = ","))
  co.cluster.matrix <- as.data.frame(matrix(data = 0, nrow = length(all.clusters), ncol = length(all.clusters)))
  rownames(co.cluster.matrix) <- all.clusters
  colnames(co.cluster.matrix) <- all.clusters
  for(i in 1:length(cluster.members)) {
    items <- unlist(strsplit(cluster.members[[i]], split=",[ ]*"))
    if(length(items) < 2) { next }
    for(j in 1:(length(items)-1)) {
      itemj <- items[j]
      for(k in (j+1):length(items)) {
        itemk <- items[k]
        co.cluster.matrix[itemj, itemk] <- co.cluster.matrix[itemj, itemk] + 1
        co.cluster.matrix[itemk, itemj] <- co.cluster.matrix[itemk, itemj] + 1
      }
    }
  }

  ## Calculate the stability scores, 
  ## From Guinney et al (2015) Nat Medicine:
  ## On the basis of the consensus matrix, we assessed the robustness of each subtype with a stability score, 
  ## which is the average frequency that its within-cluster association with other subtypes is the same as predicted 
  ## by MCL on the network generated with all samples. 
  ## Element-wise matrix multiplication
  tmp <- cons.matrix * co.cluster.matrix
  stability.scores <- unlist(apply(tmp, 1, function(row) {
                                             if(all(row == 0)) { return(0) }
                                             row <- row[row != 0]
                                             mean(row)
                                           }))
  ## Calculate the weighted silhouette distance
  ## From Guinney et al (2015) Nat Medicine: 
  ## For evaluation of clustering performance, we employed weighted Silhouette width (R package ‘WeightedCluster’), 
  ## which extends Silhouette width by giving more weights to subtypes that are more representative of their assigned 
  ## clusters. Here, we used stability scores as weights to calculate weighted Silhouette width and took the median 
  ## over all subtypes as a measure of clustering performance, which was used to evaluate the optimal number of clusters.
  ## sil <- wcSilhouetteObs(diss, clust5, weights=aggMvad$aggWeights, measure="ASWw")
  sil <- wcSilhouetteObs(as.matrix(diss), factor(mcl.all.out$meta.cluster), weights=stability.scores, measure="ASWw")
  sil[!is.finite(sil)] <- 0
  median(sil)
}