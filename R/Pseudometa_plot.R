Pseudometa_plot <- function(PseudometaObj, plot = "Gene_Volcano", method = "Spearman", 
                            title_size = 14, x_axis_text_size = 12, y_axis_text_size = 12, 
                            x_axis_title_size = 14, y_axis_title_size = 14, 
                            legend_title_size = 12, legend_text_size = 10, legend_point_size = 3) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is not installed. Please install it using install.packages('ggplot2').")
  }
  
  if (plot == "method_cor") {
    if (!is.null(PseudometaObj@analysis[["method_corr_results"]])) {
      print(ggcorrplot::ggcorrplot(PseudometaObj@analysis[["method_corr_results"]]))
    } else {
      stop("Pearson Corr results are not available in the PseudometaObj, please run pearson_method_cor().")
    }
  } else if (plot == "Gene_Volcano") {
    if (!method %in% c("Spearman", "Kendall", "Mann_Kendall", "Page")) {
      stop("Invalid method specified. Choose from 'Spearman', 'Kendall', 'Mann_Kendall', 'Page'.")
    }
    
    results_df <- PseudometaObj@analysis[[method]][["results"]]
    if (is.null(results_df)) {
      stop(paste("Results of", method, "are not available in the PseudometaObj."))
    }
    
    results_df$x_position <- ifelse(results_df$monotony == "positive", results_df$n, -results_df$n)
    
    volcano_plot <- ggplot2::ggplot(results_df, ggplot2::aes(x = x_position, y = p.value, color = monotony)) +
      ggplot2::geom_point(size = legend_point_size) +
      ggplot2::scale_color_manual(values = c("positive" = "red", "negative" = "blue")) +
      ggplot2::labs(title = paste("Gene Plot of", method, "Method"),
                    x = "Proportion of non-zero expression",
                    y = "P-value") +
      ggplot2::theme_minimal() +
      ggplot2::xlim(min(results_df$x_position) - 10, max(results_df$x_position) + 10) +
      ggplot2::ylim(0, 1.05) +
      ggplot2::scale_x_continuous(labels = function(x) format(abs(x), big.mark = ",", scientific = FALSE),
                                  breaks = scales::pretty_breaks()) +
      ggplot2::scale_y_continuous(labels = scales::label_number(), breaks = scales::pretty_breaks()) +
      ggplot2::theme(
        axis.line = ggplot2::element_line(color = "black"),
        axis.ticks = ggplot2::element_line(color = "black"),
        axis.text.x = ggplot2::element_text(color = "black", size = x_axis_text_size),
        axis.text.y = ggplot2::element_text(color = "black", size = y_axis_text_size),
        axis.title.x = ggplot2::element_text(size = x_axis_title_size),
        axis.title.y = ggplot2::element_text(size = y_axis_title_size),
        plot.title = ggplot2::element_text(size = title_size),
        legend.title = ggplot2::element_text(size = legend_title_size),
        legend.text = ggplot2::element_text(size = legend_text_size)
      )
    
    print(volcano_plot)
  } else {
    stop("Unsupported plot type. Available types: 'method_cor', 'Gene_Volcano'.")
  }
}

