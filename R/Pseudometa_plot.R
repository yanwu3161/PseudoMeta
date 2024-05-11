Pseudometa_plot <- function(PseudometaObj, plot = "method_cor", method = "Spearman") {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is not installed. Please install it using install.packages('ggplot2').")
  }
  if (!requireNamespace("ggcorrplot", quietly = TRUE)) {
    stop("Package 'ggcorrplot' is not installed. Please install it using install.packages('ggcorrplot').")
  }
  
  if (plot == "method_cor") {
    if (!is.null(PseudometaObj@analysis[["method_corr_results"]])) {
      ggcorrplot::ggcorrplot(PseudometaObj@analysis[["method_corr_results"]])
    } else {
      stop("Pearson Corr results are not available in the PseudometaObj, please run pearson_method_cor().")
    }
  } else if (plot == "Gene_Volcano") {
    if (!method %in% c("Spearman", "Kendall", "Mann_Kendall", "Page")) {
      stop("Invalid method specified. Choose from 'Spearman', 'Kendall', 'Mann_Kendall', 'Page'.")
    }
    
    # 获取指定方法的结果数据框
    results_df <- PseudometaObj@analysis[[method]][["results"]]
    if (is.null(results_df)) {
      stop(paste("Results of", method, "are not available in the PseudometaObj."))
    }
    
    # 修改x值以根据monotony正负反映n值方向
    results_df$x_position <- ifelse(results_df$monotony == "positive", results_df$n, -results_df$n)
    
    # 使用 ggplot2 绘制火山图
    volcano_plot <- ggplot2::ggplot(results_df, ggplot2::aes(x = x_position, y = p.value, color = monotony)) +
      ggplot2::geom_point() +
      ggplot2::scale_color_manual(values = c("positive" = "red", "negative" = "blue")) +
      ggplot2::labs(title = paste("Volcano Plot for", method, "Method"),
                    x = "Count of Non-zero Expressions (n), direction by monotony",
                    y = "P-value") +
      ggplot2::theme_minimal() +
      ggplot2::xlim(min(results_df$x_position) - 10, max(results_df$x_position) + 10) +
      ggplot2::ylim(0, 1.05) +
      ggplot2::scale_x_continuous(labels = function(x) format(abs(x), big.mark = ",", scientific = FALSE), 
                                  breaks = scales::pretty_breaks()) +
      ggplot2::scale_y_continuous(labels = scales::label_number(), breaks = scales::pretty_breaks()) +
      ggplot2::theme(axis.line = ggplot2::element_line(color = "black"),
                     axis.ticks = ggplot2::element_line(color = "black"),
                     axis.text = ggplot2::element_text(color = "black"))
    
    print(volcano_plot)
  } else {
    stop("Unsupported plot type. Available types: 'method_cor', 'Gene_Volcano'.")
  }
}
