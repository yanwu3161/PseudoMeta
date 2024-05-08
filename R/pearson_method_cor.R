pearson_method_cor <- function(pseudoMetaObject, method_to_cor = "ALL") {
  methods <- c("Spearman", "Kendall", "Mann_Kendall", "Page")
  if (method_to_cor != "ALL") {
    methods <- intersect(methods, method_to_cor)
  }

  cor_data <- lapply(methods, function(method) {
    pseudoMetaObject@analysis[[method]][["results"]][["p.value"]]
  })
  names(cor_data) <- methods
  method_corr_results <- matrix(NA, nrow = length(methods), ncol = length(methods), 
                                dimnames = list(methods, methods))
  
  for (i in 1:length(methods)) {
    for (j in i:length(methods)) {
      if (i == j) {
        method_corr_results[i, j] <- 1
      } else {
        method_corr_results[i, j] <- cor(cor_data[[methods[i]]], cor_data[[methods[j]]], use = "complete.obs")
        method_corr_results[j, i] <- method_corr_results[i, j]
      }
    }
  }
  
  method_corr_results <- as.data.frame(method_corr_results)

  pseudoMetaObject@analysis[["method_corr_results"]] <- method_corr_results

  return(pseudoMetaObject)
}
