# 06_differential_abundance_LinDA.R
# This script handels the plotting of the linda output results, to allow for a good customisation

# Identify the numeric position of the selected variable to plot
idx <- which(names(linda_res$output) == var_to_plot)
plots <- linda.plot(
  linda.obj = linda_res,
  variables.plot = names(linda_res$output), # This ensures all variables are plotted
  alpha = 0.05,                             # Significance threshold
  lfc.cut = 1                               # Log-Fold Change threshold for highlighting
)
# Extract the data using the index
plot_data <- plots$plot.lfc[[idx]]$data

# Clean and Filter
plot_data_clean <- plot_data %>%
  filter(bias == "Debiased") %>%
  filter(Taxa != "27F-1492R") %>% # Remove specific noise
  filter(abs(Log2FoldChange) >= 1) %>%
  mutate(Taxa = reorder(Taxa, Log2FoldChange)) # Reorder for the plot

# 6. The "Good Plot"
good_plot <- ggplot(plot_data_clean, aes(x = Taxa, y = Log2FoldChange)) +
  geom_pointrange(aes(
    ymin = Log2FoldChange - lfcSE,
    ymax = Log2FoldChange + lfcSE,
    color = Log2FoldChange < 0
  )) +
  scale_color_manual(values = c("FALSE" = "#00aeb5", "TRUE" = "indianred1")) +
  coord_flip() +
  theme_minimal() +
  labs(
    title = paste("Effect of", var_to_plot),
    subtitle = paste("Agglomerated at", target_level_linda, "level"),
    x = "", y = "Log2 Fold Change (Debiased)"
  ) +
  theme(axis.text.y = element_text(size = 11), legend.position = "none")
plot(good_plot)
# 7. Safe Saving
output_dir <- "results/linda/plots"
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

safe_name <- gsub("[^a-zA-Z0-9]", "_", var_to_plot)
ggsave(file.path(output_dir, paste0("final_linda_", safe_name, ".png")), 
       plot = good_plot, width = 8, height = 7)
