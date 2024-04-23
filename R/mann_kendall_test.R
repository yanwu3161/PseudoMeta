mann_kendall_test <- function(TrendMatrix, geneNames, min.cell.expression, pb) {
  MannKendall_results <- data.frame(n = numeric(length(geneNames)), p.value = numeric(length(geneNames)), 
                                    monotony = character(length(geneNames)), S = numeric(length(geneNames)),
                                    varS = numeric(length(geneNames)), row.names = geneNames)
  if (!requireNamespace("Kendall", quietly = TRUE)) {
    stop(paste("Kendall_pkg_not installed or loaded. Please install it to use mann_kendall_test."))
  } else {
    library("Kendall", character.only = TRUE)
  }
  for (i in seq_along(geneNames)) {
    setTxtProgressBar(pb, i)
    gene_data <- TrendMatrix[geneNames[i], ]
    non_zero_data_indices <- which(gene_data != 0)
    
    if (length(non_zero_data_indices) >= min.cell.expression) {
      mk_result <- Kendall::MannKendall(gene_data)
      
      MannKendall_results[geneNames[i], "p.value"] <- abs(mk_result$tau)
      MannKendall_results[geneNames[i], "monotony"] <- ifelse(mk_result$tau > 0, "positive", "negative")
      MannKendall_results[geneNames[i], "n"] <- length(gene_data)
      MannKendall_results[geneNames[i], "S"] <- mk_result$S
      MannKendall_results[geneNames[i], "varS"] <- mk_result$varS
    } else {
      MannKendall_results[geneNames[i], "p.value"] <- NA
      MannKendall_results[geneNames[i], "monotony"] <- NA
      MannKendall_results[geneNames[i], "n"] <- NA
      MannKendall_results[geneNames[i], "S"] <- NA
      MannKendall_results[geneNames[i], "varS"] <- NA
    }
  }
  return(MannKendall_results)
}