# scripts/05_beta_diversity_visualisation.R

# 1. Base Plot Creation
p_beta <- plot_ordination(ps_beta_input, ord_beta, color = color_var, shape = shape_var) + 
  geom_point(size = 4, alpha = 0.7)

# 2. CONDITIONAL: Add Grouping Ellipses
if (show_ellipses) {
  p_beta <- p_beta + 
    stat_ellipse(aes(group = .data[[group_clustering]]), linetype = 2, alpha = 0.5)
}

# 3. CONDITIONAL: Sample Labeling Logic
if (show_labels) {
  p_beta <- p_beta + 
    geom_text_repel(aes(label = sample_names(ps_beta_input)), 
                    size = 3, 
                    max.overlaps = 15, 
                    box.padding = 0.5,
                    point.padding = 0.3,
                    segment.color = 'grey50')
}

# 4. Styling and Labels (Appended to whatever layers were built above)
p_beta <- p_beta + 
  theme_bw() +
  labs(
    title = paste(toupper(beta_metric), "Ordination"),
    subtitle = paste("PERMANOVA p-value:", format.pval(beta_p_val, digits = 3)),
    color = color_var,
    caption = paste("Method:", ord_method, "| Data:", clr_variant)
  ) +
  theme(
    aspect.ratio = 1, 
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

# 5. Save Output
if(!dir.exists("results")) dir.create("results")

# Cleaned up file naming logic to match whether it is labeled or not
suffix <- if (show_labels) "_labeled" else ""
output_path <- paste0("results/beta_", beta_metric, suffix, ".png")

ggsave(output_path, p_beta, width = 12, height = 8, dpi = 300)
message("Beta diversity plot generated: ", output_path)


