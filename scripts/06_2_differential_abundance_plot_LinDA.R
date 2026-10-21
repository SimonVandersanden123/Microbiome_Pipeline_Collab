# ==============================================================================
# FIXED 06_DIFFERENTIAL_ABUNDANCE_LINDA.R (INTERACTION COMPATIBLE RETRIEVAL)
# ==============================================================================
library(dplyr)
library(ggplot2)

# Global variables protection fallback checks (fixes the 'object not found' errors)
if (!exists("var_to_plot")) {
  var_to_plot <- "time_continuous:sample_depth_information0.60-0.80"
}
if (!exists("log_fold_cutoff_plotting")) {
  log_fold_cutoff_plotting <- 1 # Standard 5-fold change cutoff boundary
}

# 1. Safely check if the target contrast actually exists in the LinDA output list
if (!var_to_plot %in% names(linda_res$output)) {
  stop(paste("Error: The requested variable", var_to_plot, "does not exist in the LinDA output structure."))
}

# 2. Extract raw statistics directly bypassing the restrictive 'linda.plot()' engine
raw_metrics <- linda_res$output[[var_to_plot]]
# 2.1 Define the colours for the positive and negative values for the logratio
color_positive ="#2A6F97"
color_negative ="#FF2400"

# 3. Format matrix rows into an accessible data frame matching standard plot metrics
# Convert raw LinDA matrix into clean dataframe
plot_data_base <- data.frame(
  Taxa           = rownames(raw_metrics),
  Log2FoldChange = raw_metrics$log2FoldChange,
  lfcSE          = raw_metrics$lfcSE, 
  padj           = raw_metrics$padj,
  reject         = raw_metrics$reject
) %>%
  filter(Taxa != "27F-1492R") %>% 
  filter(reject == TRUE)


# ------------------------------------------------------------------------------
# 2. DYNAMIC FILTERING ENGINE
# ------------------------------------------------------------------------------
if (filter_mode == "top_n") {
  
  # MODE A: Directional Top N Slicing (Top N Positive + Top N Negative)
  plot_data_clean <- plot_data_base %>%
    group_by(direction = ifelse(Log2FoldChange > 0, "Positive", "Negative")) %>%
    slice_max(order_by = abs(Log2FoldChange), n = top_n_taxa, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(Taxa = reorder(Taxa, Log2FoldChange))
  
  subtitle_text <- paste0("Top ", top_n_taxa, " positively & negatively enriched taxa (LinDA, FDR < 0.05)")
  
} else if (filter_mode == "log_fold") {
  
  # MODE B: Absolute Log2 Fold-Change Thresholding
  plot_data_clean <- plot_data_base %>%
    filter(abs(Log2FoldChange) >= log_fold_cutoff_plotting) %>%
    mutate(Taxa = reorder(Taxa, Log2FoldChange))
  
  subtitle_text <- paste0("Significant taxa with FDR < 0.05 and |Log2FC| ≥ ", log_fold_cutoff_plotting)
  
} else {
  stop("Invalid 'filter_mode' specified. Please choose either 'top_n' or 'log_fold'.")
}

# 5. Check if any taxa passed your filters before attempting to plot
if (nrow(plot_data_clean) == 0) {
  warning("No microbial families met the significant p-value and fold-change threshold filters for this interaction.")
} else {
  
  # 6. Generate the custom Plot
  good_plot <- ggplot(plot_data_clean, aes(x = Taxa, y = Log2FoldChange)) +
    geom_pointrange(aes(
      ymin = Log2FoldChange - lfcSE,
      ymax = Log2FoldChange + lfcSE,
      color = Log2FoldChange < 0
    ), size = 0.6) +
    scale_color_manual(values = c("FALSE" = color_positive, "TRUE" = color_negative)) + # Safe distinct ecological colors
    coord_flip() +
    theme_minimal() +
    labs(
      title = var_to_plot,
      subtitle = paste(" "),
      x = "Taxon", 
      y = "Log2 Fold change (Debiased)"
    ) +
    theme(
      axis.text.y = element_text(size = 10, face = "italic"), 
      legend.position = "none",
      panel.grid.minor = element_blank()
    )
  
  # Print plot to screen
  print(good_plot)
  
  # 7. Safe Direct Export
  output_dir <- "results/linda/plots"
  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  safe_name <- gsub("[^a-zA-Z0-9]", "_", var_to_plot)
  ggsave(file.path(output_dir, paste0("final_linda_", safe_name, ".png")), 
         plot = good_plot, width = 11, height = 8, dpi = 300)
}