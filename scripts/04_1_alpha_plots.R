# 03_alpha_plots.R
library(ggplot2)
library(tidyr)
library(dplyr)

# 1. Dynamic Data Preparation
# ---------------------------
plot_data <- alpha_output %>%
  # Ensure the X-axis and Color columns are treated as factors for discrete plotting
  mutate(across(all_of(c(plot_x_axis, plot_color_by)), as.factor)) %>%
  select(all_of(c(plot_x_axis, plot_color_by)), all_of(plot_metrics)) %>%
  pivot_longer(cols = all_of(plot_metrics), 
               names_to = "Metric", 
               values_to = "Value")

# 2. Build the Plot Dynamically
# -----------------------------
alpha_boxplot <- ggplot(plot_data, aes_string(x = plot_x_axis, y = "Value", fill = plot_color_by)) +
  # Use geom_boxplot for distribution
  geom_boxplot(alpha = 0.6, outlier.shape = NA, width = 0.7) +
  # Add jittered points so colleagues can see the individual replicates/N
  geom_point(position = position_jitterdodge(jitter.width = 0.1), size = 1.2, alpha = 0.4) +
  # Separate by Metric with free Y scales
  facet_wrap(~Metric, scales = "free_y", ncol = plot_cols) +
  theme_bw() +
  labs(title = "Microbial Community Trait Analysis",
       subtitle = paste("Grouped by:", plot_x_axis, "| Colored by:", plot_color_by),
       x = plot_x_axis,
       y = "Calculated Index Value",
       fill = plot_color_by) +
  theme(
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

# 3. Save Output
ggsave("../results/alpha_diversity_dynamic_plot.png", alpha_boxplot, width = 12, height = 8)