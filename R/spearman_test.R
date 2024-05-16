spearman_test <- function(TrendMatrix, geneNames, min.cell.expression, pb) {
  Spearman_results <- data.frame(n = numeric(length(geneNames)), p.value = numeric(length(geneNames)), 
                                 monotony = character(length(geneNames)), row.names = geneNames)
  
  for (i in seq_along(geneNames)) {
    setTxtProgressBar(pb, i)
    gene_data <- TrendMatrix[geneNames[i], ]
    non_zero_data_indices <- which(gene_data != 0)
    
    if (length(non_zero_data_indices) >= min.cell.expression) {
      gene_data <- gene_data[non_zero_data_indices]
      n <- length(gene_data)
      
      if (n > 1) {
        ranks <- rank(gene_data)
        d_i <- ranks - seq_along(non_zero_data_indices)
        sum_d_i_squared <- sum(d_i^2)
        
        rho <- 1 - (6 * sum_d_i_squared / (n * (n^2 - 1)))
        Spearman_results[geneNames[i], "p.value"] <- abs(rho)
        Spearman_results[geneNames[i], "monotony"] <- ifelse(rho > 0, "positive", "negative")
        Spearman_results[geneNames[i], "n"] <- n / ncol(TrendMatrix)
      } else {
        Spearman_results[geneNames[i], "p.value"] <- NA
        Spearman_results[geneNames[i], "monotony"] <- "NA"
      }
    } else {
      Spearman_results[geneNames[i], "p.value"] <- NA
      Spearman_results[geneNames[i], "monotony"] <- "NA"
    }
  }
  return(Spearman_results)
}