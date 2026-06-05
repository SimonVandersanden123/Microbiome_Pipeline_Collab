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

# 3. Format matrix rows into an accessible data frame matching standard plot metrics
plot_data_clean <- data.frame(
  Taxa           = rownames(raw_metrics),
  Log2FoldChange = raw_metrics$log2FoldChange,
  lfcSE          = raw_metrics$lfcSE, # Standard Error calculated by LinDA
  padj           = raw_metrics$padj,
  reject         = raw_metrics$reject
) %>%
  # 4. Filter out background noise & keep only significant threshold parameters
  filter(Taxa != "27F-1492R") %>% 
  filter(reject == TRUE) %>% 
  filter(abs(Log2FoldChange) >= log_fold_cutoff_plotting) %>%
  mutate(Taxa = reorder(Taxa, Log2FoldChange))

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
    scale_color_manual(values = c("FALSE" = colour_positive, "TRUE" = colour_negative)) + # Safe distinct ecological colors
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