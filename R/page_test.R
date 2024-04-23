page_test <- function(TrendMatrix, geneNames, min.cell.expression, pb) {
  Page_results <- data.frame(L = numeric(length(geneNames)), n = numeric(length(geneNames)),
                             p.value = numeric(length(geneNames)), monotony = character(length(geneNames)),
                             Z = numeric(length(geneNames)), row.names = geneNames)
  
  for (i in seq_along(geneNames)) {
    setTxtProgressBar(pb, i)
    gene_data <- TrendMatrix[geneNames[i], ]
    non_zero_data_indices <- which(gene_data != 0)
    
    if (length(non_zero_data_indices) >= min.cell.expression) {
      n <- length(gene_data)
      ranks <- rank(gene_data)
      
      L <- sum((1:n) * ranks)
      mu_L <- n * (n + 1)^2 / 4
      sigma_L <- sqrt(n * (n + 1)^2 * (2 * n + 1) / 72)
      Z <- (L - mu_L) / sigma_L
      
      p_value <- 2 * pnorm(-abs(Z))
      
      Page_results[geneNames[i], "L"] <- L
      Page_results[geneNames[i], "n"] <- n
      Page_results[geneNames[i], "p.value"] <- p_value
      Page_results[geneNames[i], "monotony"] <- ifelse(Z > 0, "positive", "negative")
      Page_results[geneNames[i], "Z"] <- Z
    } else {
      Page_results[geneNames[i], ] <- NA
    }
  }
  return(Page_results)
}