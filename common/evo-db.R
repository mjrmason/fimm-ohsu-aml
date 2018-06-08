

getGenesfromDrugs <- function(drugnames, input_struct, tanimoto_threshold, parallelized = T){
  if(parallelized == T){
    foo <- mclapply(input_struct, function(x){
      structure <- as.character(x)
      sims <- similarityFunction(structure)
      getSimMols(sims, tanimoto_threshold)
      }, mc.cores = detectCores())
    names(foo) <- drugnames
    bar <- ldply(foo)
    colnames(bar) <- (c("input_drug", "common_name", "Tanimoto Similarity"))
  }
  
  if(parallelized == F){
    foo <- lapply(input_struct, function(x){
      structure <- as.character(x)
      sims <- similarityFunction(structure)
      getSimMols(sims, tanimoto_threshold)
    })
    names(foo) <- drugnames
    bar <- ldply(foo)
    colnames(bar) <- (c("input_drug", "common_name", "Tanimoto Similarity"))
  }
  
  bar <- inner_join(db, bar, by = "common_name")
  
  
}

