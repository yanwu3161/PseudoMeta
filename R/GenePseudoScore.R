GenePseudoScore <- function(Pseudometaobj, method = "Spearman", Gene_p_TH = 0.25, Gene_n_TH = 0.05, Minmaxnorm = TRUE) {
  Varexpr_matrix <- Pseudometaobj@assays[["var.matrix"]]
  trendtest_results <- Pseudometaobj@analysis[[method]][["results"]][c("n", "p.value", "monotony")]
  
  trendtest_results <- trendtest_results[trendtest_results$'p.value' >= Gene_p_TH, ]
  trendtest_results <- trendtest_results[trendtest_results$'n' >= Gene_n_TH, ]
  Varexpr_matrix <- Varexpr_matrix[rownames(Varexpr_matrix) %in% rownames(trendtest_results), ]
  
  cat("Calculating the colrank of var.matrix.\n")
  Varexpr_matrix[Varexpr_matrix == 0] <- NA
  Varranks_matrix <- apply(Varexpr_matrix, 2, rank, na.last = "keep")
  cat("Calculating The timescore of each cell.\n")
  
  adjusted_CellScore_matrix <- data.frame(matrix(nrow = ncol(Varranks_matrix), ncol = 1))
  rownames(adjusted_CellScore_matrix) <- colnames(Varranks_matrix)
  colnames(adjusted_CellScore_matrix) <- "CellScore"
  pb <- txtProgressBar(min = 0, max = ncol(Varranks_matrix), style = 3)
  
  for (Cell_to_cal_Score in colnames(Varranks_matrix)) {
    OneCell_EachGene_Score <- data.frame(matrix(nrow = nrow(Varranks_matrix), ncol = 1))
    rownames(OneCell_EachGene_Score) <- rownames(Varranks_matrix)
    for (Var_gene_to_score in rownames(Varranks_matrix)) {
      if (Var_gene_to_score %in% rownames(trendtest_results)) {
        p_value <- trendtest_results[Var_gene_to_score, "p.value"]
        mono <- trendtest_results[Var_gene_to_score, "monotony"]
        scores <- Varranks_matrix[Var_gene_to_score, Cell_to_cal_Score] * p_value
        
        if (mono == "negative") {
          scores <- -scores
        }
        OneCell_EachGene_Score[Var_gene_to_score,] <- scores
      }
    }
    adjusted_CellScore_matrix[Cell_to_cal_Score,] <- sum(OneCell_EachGene_Score[, 1], na.rm = TRUE)
    setTxtProgressBar(pb, which(colnames(Varranks_matrix) == Cell_to_cal_Score))
  }
  close(pb)
  
  columnName <- paste(method, "TimeScore", sep = "_")
  if (Minmaxnorm) {
    min_val <- min(adjusted_CellScore_matrix$CellScore)
    max_val <- max(adjusted_CellScore_matrix$CellScore)
    normalized_scores <- (adjusted_CellScore_matrix$CellScore - min_val) / (max_val - min_val)
    Pseudometaobj@meta.data[[columnName]] <- normalized_scores
  } else {
    Pseudometaobj@meta.data[[columnName]] <- adjusted_CellScore_matrix$CellScore
  }
  
  return(Pseudometaobj)
}
