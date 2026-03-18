# scripts/05_beta_diversity_visualisation.R
# 1. Base Plot Creation
p_beta <- plot_ordination(ps_beta_input, ord_beta, color = color_var, shape = shape_var) + 
  geom_point(size = 4, alpha = 0.7) + 
  
  # 2. Add Grouping Ellipses
  # Note: This uses the group_clustering variable defined in your Rmd
  stat_ellipse(aes(group = .data[[group_clustering]]), linetype = 2, alpha = 0.5) + 
  
  # 3. ADDED: Sample Labeling Logic
  # We use ggrepel to prevent labels from overlapping the points
  geom_text_repel(aes(label = sample_names(ps_beta_input)), 
                  size = 3, 
                  max.overlaps = 15, # Adjust this if too many labels disappear
                  box.padding = 0.5,
                  point.padding = 0.3,
                  segment.color = 'grey50') +
<<<<<<< HEAD
<<<<<<< HEAD

=======
  
>>>>>>> parent of 7628505 (Revert "Simon adjustments to main")
=======
  
>>>>>>> parent of 7628505 (Revert "Simon adjustments to main")
  # 4. Styling and Labels
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
ggsave(paste0("results/beta_", beta_metric, "_labeled.png"), p_beta, width = 9, height = 7)
<<<<<<< HEAD
<<<<<<< HEAD
message("Beta diversity plot with labels generated.")
ggsave(paste0("results/beta_", beta_metric, ".png"), p_beta, width = 8, height = 7)

message("Beta diversity plot generated.")
=======

message("Beta diversity plot with labels generated.")
>>>>>>> parent of 7628505 (Revert "Simon adjustments to main")
=======

message("Beta diversity plot with labels generated.")
>>>>>>> parent of 7628505 (Revert "Simon adjustments to main")
