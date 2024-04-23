kendall_test <- function(TrendMatrix, geneNames, min.cell.expression, pb) {
  Kendall_results <- data.frame(n = numeric(length(geneNames)), p.value = numeric(length(geneNames)), 
                                monotony = character(length(geneNames)), row.names = geneNames)

    if (!requireNamespace("Kendall", quietly = TRUE)) {
      stop(paste("Kendall_pkg_not installed or loaded. Please install it to use kendall_test."))
    } else {
      library("Kendall", character.only = TRUE)
    }
    
  for (i in seq_along(geneNames)) {
    setTxtProgressBar(pb, i)
    gene_data <- TrendMatrix[geneNames[i], ]
    non_zero_data_indices <- which(gene_data != 0)
    gene_data <- gene_data[non_zero_data_indices]
    
    if (length(non_zero_data_indices) >= min.cell.expression) {
      time_data <- 1:length(gene_data)
      kendall_result <- Kendall(time_data, gene_data)
      
      Kendall_results[geneNames[i], "p.value"] <- abs(kendall_result$tau)
      Kendall_results[geneNames[i], "monotony"] <- ifelse(kendall_result$tau > 0, "positive", "negative")
      Kendall_results[geneNames[i], "n"] <- length(gene_data)
    } else {
      Kendall_results[geneNames[i], "p.value"] <- NA
      Kendall_results[geneNames[i], "monotony"] <- "NA"
    }
  }
  return(Kendall_results)
}