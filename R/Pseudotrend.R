Pseudotrend <- function(pseudoMetaObject, method = "Spearman", min.cell.expression = 5) {
  pseudotime <- pseudoMetaObject@meta.data[["pseudotime"]]
  non_zero_indices <- which(pseudotime != 0)
  sorted_indices <- order(pseudotime[non_zero_indices])
  TrendMatrix <- pseudoMetaObject@assays[["var.matrix"]][, non_zero_indices][, sorted_indices]
  
  geneNames <- rownames(TrendMatrix)
  pb <- txtProgressBar(min = 0, max = length(geneNames), style = 3)
  
  method_function_map <- list(
    
    "Kendall" = function() kendall_test(TrendMatrix, geneNames, min.cell.expression, pb),
    "Spearman" = function() spearman_test(TrendMatrix, geneNames, min.cell.expression, pb),
    "Mann_Kendall" = function() mann_kendall_test(TrendMatrix, geneNames, min.cell.expression, pb),
    "Page" = function() page_test(TrendMatrix, geneNames, min.cell.expression, pb)
  )
  
  if (!method %in% names(method_function_map)) {
    stop("Invalid method. Choose from 'Spearman', 'Kendall', 'Mann_Kendall', 'Page'.")
  }
  
  results <- method_function_map[[method]]()  # 调用相应的方法
  
  close(pb)
  pseudoMetaObject@analysis[[method]] <- list(results = results)
  
  return(pseudoMetaObject)
}
